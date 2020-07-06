#include "../include/NanoporeReads.cuh"

NanoporeReads::NanoporeReads(char *fileName, int k, int n) : k(k), n(n), sketches(NULL) {
    std::ifstream infile(fileName);
    std::string line;
    numReads = 0;
    while (std::getline(infile, line)) {
//        std::cout << line << std::endl;
        size_t index = line.find(':');
        readPos.push_back(std::stol(line.substr(0, index)));
        {
            std::unique_ptr<std::string> ptr(new std::string(line.substr(index + 1)));
            editStrings.push_back(std::move(ptr));
        }
//        std::cout << readPos.back() << std::endl;
//        std::cout << editStrings.back() << std::endl;
        std::getline(infile, line);
        {
            std::unique_ptr<std::string> ptr(new std::string(line));
            readData.push_back(std::move(ptr));
        }
        numReads++;
    }
    readLen = readData[0]->length();
    readPosSorted = readPos;
    std::sort(readPosSorted.begin(), readPosSorted.end());
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "readLen " << readLen << std::endl;
}

void NanoporeReads::calculateMinHashSketches() {
    // We store all the k-mers as uint64s. This would work for all k<=32,
    // which is definitely sufficient
    const size_t numKMers = readLen - k + 1;
    kMer_t *kMers;
    // Because of memory constraints on the GPUs we cannot deal with all the reads at once.
    // So we arrange the reads into blocks of blockSize reads and only work on a single block
    // at the same time.
    const size_t blockSize = 2048;
    std::cout << "numKMers " << numKMers << std::endl;

    cudaMallocManaged(&(sketches), n * numReads * sizeof(kMer_t));

    std::random_device rd;
    std::mt19937_64 gen(rd());

    /* This is where you define the number generator for unsigned long long: */
    std::uniform_int_distribution<unsigned long long> dis;

    kMer_t *randNumbers;
    cudaMallocManaged(&randNumbers, n * sizeof(kMer_t));
    for (size_t i = 0; i < n; ++i) {
        randNumbers[i] = dis(gen);
    }

    for (size_t currentRead = 0; currentRead < numReads; currentRead += blockSize) {
        std::cout << "CurrentRead " << currentRead << std::endl;
        const long readsLeft = numReads - (long) currentRead;
        const size_t currentBlockSize = readsLeft > blockSize ? blockSize : readsLeft;

        auto generateKMers = [&]() {
            cudaMallocManaged(&kMers, currentBlockSize * numKMers * sizeof(kMer_t));

            for (size_t i = 0; i < currentBlockSize; i++) {
                size_t baseIndex = i * numKMers;
#pragma omp parallel for
                for (size_t index = baseIndex; index < numKMers + baseIndex; index++) {
                    kMers[index] =
                            kMerToInt(readData[i + currentRead]->substr(
                                    index - baseIndex, k));
                }
            }
        };

        generateKMers();

//        for (size_t i = 0; i < currentBlockSize * numKMers; ++i) {
//            std::cout << i << " " << kMers[i] << std::endl;
//        }

        // Now we generate all hashes
        // hashes is indexed by (read number, k-mer number, hash number)
        kMer_t *hashes;
        cudaMallocManaged(&hashes, n * currentBlockSize * numKMers * sizeof(kMer_t));

        const size_t blockSize = 512;
        const size_t numBlocks = 512;
        hashKMer <<< numBlocks, blockSize >>>(currentBlockSize * numKMers,
                                              n, kMers, hashes, randNumbers);
        // Finish calculating the hashes and frees unneeded memory
        cudaDeviceSynchronize();

//        for (size_t i = 0; i < currentBlockSize * numKMers * n; ++i) {
//            std::cout << i << " " << hashes[i] << std::endl;
//        }

        // Now we are going to compute the sketches which are the minimums of the hashes
        calcSketch
        <<< (currentBlockSize + blockSize - 1)
        / blockSize, blockSize >>>(currentBlockSize,
                                   currentRead, numKMers,
                                   n, hashes,
                                   sketches, kMers);
        cudaDeviceSynchronize();
        cudaFree(kMers);
        cudaFree(hashes);
    }
    cudaFree(randNumbers);
    populateHashTables();
}

kMer_t NanoporeReads::kMerToInt(const std::string &s) {
    size_t l = s.length();
    kMer_t result = 0;
    for (size_t i = 0; i < l; ++i) {
        result <<= 2;
        result |= baseToInt(s[i]);
    }
    return result;
}

char NanoporeReads::baseToInt(const char base) {
    switch (base) {
        case 'A':
            return 0;
        case 'T':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 3;
        default:
            std::cout << "Oh No!" << std::endl;
            return 0;
    }
}

__global__ void hashKMer(const size_t totalKMers, const size_t n,
                         kMer_t *kMers, kMer_t *hashes, kMer_t *randNumbers) {
    size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = blockDim.x * gridDim.x;
    for (size_t i = index; i < totalKMers; i += stride) {
        size_t hashIndex = i * n;
        kMer_t currentHash = kMers[i];
        currentHash = (currentHash * (uint64_t) HASH_C64);
        currentHash ^= randNumbers[0];
        hashes[hashIndex++] = currentHash;
        for (size_t j = 1; j < n; j++) {
            kMer_t newHash = ((currentHash >> ROTATE_BITS)
                              | (currentHash << (KMER_BITS - ROTATE_BITS)))
                             ^0xABCD32108475AC38;
            newHash = (newHash * (uint64_t) HASH_C64);
            newHash ^= randNumbers[j];
            newHash += currentHash;
            currentHash = newHash;
            hashes[hashIndex++] = currentHash;
        }
    }
}

__global__ void calcSketch(const size_t numReads, const size_t currentRead,
                           const size_t numKMers, const size_t n,
                           kMer_t *hashes, kMer_t *sketches, kMer_t *kMers) {
    size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = blockDim.x * gridDim.x;
    for (size_t i = index; i < numReads; i += stride) {
        size_t sketchIndex = (i + currentRead) * n;
        for (size_t j = 0; j < n; ++j) {
            size_t hashIndex = i * n * numKMers + j;
            kMer_t currentMin = ~(kMer_t) 0;
            size_t minIndex = 0;
            for (size_t l = 0; l < numKMers; ++l) {
                kMer_t temp = hashes[hashIndex];
                hashIndex += n;
                minIndex = currentMin < temp ? minIndex : l;
                currentMin = currentMin < temp ? currentMin : temp;
            }
            //std::cout << "thread: " << i << " hash id: " << j << std::endl;
            //std::cout << "minIndex " << minIndex << ":" << kMers[i * numKMers + minIndex];

            if (kMers)
                sketches[sketchIndex++] = kMers[i * numKMers + minIndex];
            else
                sketches[sketchIndex++] = currentMin;
        }
    }
}

NanoporeReads::~NanoporeReads() {
    cudaFree(sketches);
}

void NanoporeReads::printHashes() {
    for (size_t i = 0; i < numReads; ++i) {
        std::cout << readPos[i];
        for (size_t j = 0; j < n; ++j) {
            std::cout << ", \"" << sketches[i * n + j] << "\"";
        }
        std::cout << std::endl;
    }
}

void NanoporeReads::populateHashTables() {
    std::cout << "Starting to populate hash tables" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < n; ++i) {
        hashTables.push_back(std::map<kMer_t, std::vector<size_t >>());
    }
//#pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        std::map<kMer_t, std::vector<size_t>> &hT = hashTables[i];
        size_t currentIndex = i;
        for (size_t j = 0; j < numReads; ++j) {
            hT[sketches[currentIndex]].push_back(j);
//            std::cout << "Inserting " << j << " to " << sketches[currentIndex] << std::endl;
            currentIndex += n;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "finished populating hash tables" << std::endl;
    std::cout << duration.count() << " milliseconds passed" << std::endl;
}

filterStats NanoporeReads::getFilterStats(unsigned int overlapBaseThreshold, unsigned int overlapSketchThreshold) {
    filterStats result(overlapBaseThreshold, overlapSketchThreshold);
    // First we calculate the number of overlaps and disjoints
    for (size_t i = 0; i < numReads; ++i) {
        long curTh = readPosSorted[i] + readLen - overlapBaseThreshold;
        for (size_t j = i + 1; j < numReads; ++j) {
            if (((long) readPosSorted[j]) <= curTh)
                result.numOverlaps++;
            else
                break;
        }
    }
    result.numDisjoint = (numReads * (unsigned long long) (numReads - 1)) / 2 - result.numOverlaps;

    // Now we calculate falsePositives, falseNegatives, etc
    for (size_t i = 0; i < numReads; ++i) {
        std::multiset<size_t> matches;
        unsigned long curPos = readPos[i];
        long th = readLen - overlapBaseThreshold;
        for (size_t sketchIndex = 0; sketchIndex < n; ++sketchIndex) {
            kMer_t curHash = sketches[i * n + sketchIndex];
//            std::cout << i << " " << sketchIndex << " " << curHash << std::endl;
//            auto currentMap = hashTables[sketchIndex];
//            for (auto p : currentMap) {
//                std::cout << p.first << " ";
//            }
//            std::cout << std::endl;
            std::vector<size_t> &m = hashTables[sketchIndex].at(curHash);
            matches.insert(m.begin(), m.end());
        }
        auto end = matches.end();
        for (auto it = matches.begin(); it != end; it = matches.upper_bound(*it)) {
            if (*it <= i)
                continue;
            unsigned long pos = readPos[*it];
            if (matches.count(*it) >= overlapSketchThreshold) {
                if (abs((long) pos - (long) curPos) > th)
                    result.falsePositives++;
                result.totalPositive++;
            }
        }
//        std::cout << std::endl;
    }

    result.totalNegative = result.numOverlaps + result.numDisjoint - result.totalPositive;
    result.falseNegatives = result.numOverlaps - result.totalPositive + result.falsePositives;
    return result;
}

filterStats::filterStats(unsigned int overlapBaseThreshold, unsigned int overlapSketchThreshold) :
        overlapBaseThreshold(overlapBaseThreshold), overlapSketchThreshold(overlapSketchThreshold),
        totalPositive(0), totalNegative(0), numOverlaps(0), numDisjoint(0),
        falsePositives(0), falseNegatives(0) {
}

std::ostream &operator<<(std::ostream &out, const filterStats &o) {
    const int w = 13;
    out << std::setw(w) << "overlapBaseTh" << ","
        << std::setw(w) << "numKMerTh" << ","
        << std::setw(w) << "totalPos" << ","
        << std::setw(w) << "totalNeg" << ","
        << std::setw(w) << "numOverlaps" << ","
        << std::setw(w) << "numDisjoint" << ","
        << std::setw(w) << "falsePos" << ","
        << std::setw(w) << "falseNeg" << std::endl;
    out << std::setw(w) << o.overlapBaseThreshold << ","
        << std::setw(w) << o.overlapSketchThreshold << ","
        << std::setw(w) << o.totalPositive << ","
        << std::setw(w) << o.totalNegative << ","
        << std::setw(w) << o.numOverlaps << ","
        << std::setw(w) << o.numDisjoint << ","
        << std::setw(w) << o.falsePositives << ","
        << std::setw(w) << o.falseNegatives << std::endl;
    return out;
}