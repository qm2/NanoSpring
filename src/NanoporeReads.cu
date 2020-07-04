#include "../include/NanoporeReads.cuh"

NanoporeReads::NanoporeReads(char *fileName, int k, int n) : k(k), n(n), sketches(NULL) {
    std::ifstream infile(fileName);
    std::string line;
    numReads = 0;
    while (std::getline(infile, line)) {
//        std::cout << line << std::endl;
        size_t index = line.find(':');
        this->readPos.push_back(std::stol(line.substr(0, index)));
        {
            std::unique_ptr <std::string> ptr(new std::string(line.substr(index + 1)));
            this->editStrings.push_back(std::move(ptr));
        }
//        std::cout << this->readPos.back() << std::endl;
//        std::cout << this->editStrings.back() << std::endl;
        std::getline(infile, line);
        {
            std::unique_ptr <std::string> ptr(new std::string(line));
            this->readData.push_back(std::move(ptr));
        }
        numReads++;
    }
    this->readLen = this->readData[0]->length();
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "readLen " << readLen << std::endl;
}

void NanoporeReads::calculateMinHashSketches() {
    // We store all the k-mers as uint64s. This would work for all k<=32,
    // which is definitely sufficient
    const size_t numKMers = this->readLen - this->k + 1;
    kMer_t *kMers;
//    const size_t totalKMers = this->numReads * numKMers;
    const size_t blockSize = 2048;
    std::cout << "numKMers " << numKMers << std::endl;

    cudaMallocManaged(&(this->sketches), this->n * this->numReads * sizeof(kMer_t));

    for (size_t currentRead = 0; currentRead < this->numReads; currentRead += blockSize) {
        std::cout << "CurrentRead " << currentRead << std::endl;
        const long readsLeft = numReads - (long) currentRead;
        const size_t currentBlockSize = readsLeft > blockSize ? blockSize : readsLeft;

        auto generateKMers = [&]() {
            cudaMallocManaged(&kMers, currentBlockSize * numKMers * sizeof(kMer_t));

//#pragma omp parallel for collapse(2)
            for (size_t i = 0; i < currentBlockSize; i++) {
                size_t baseIndex = i * numKMers;
#pragma omp parallel for
                for (size_t index = baseIndex; index < numKMers + baseIndex; index++) {
                    kMers[index] =
                            kMerToInt(readData[i + currentRead]->substr(
                                    index - baseIndex, this->k));
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
        cudaMallocManaged(&hashes, this->n * currentBlockSize * numKMers * sizeof(kMer_t));

        const size_t blockSize = 512;
        const size_t numBlocks = 512;
        hashKMer <<< numBlocks, blockSize >>>(currentBlockSize * numKMers,
                                              this->n, kMers, hashes);
        // Finish calculating the hashes and frees unneeded memory
        cudaDeviceSynchronize();

        // Now we are going to compute the sketches which are the minimums of the hashes
        calcSketch
        <<< (currentBlockSize + blockSize - 1) / blockSize, blockSize >>>(currentBlockSize,
                                                                          currentRead, numKMers,
                                                                          this->n, hashes,
                                                                          this->sketches);
        cudaDeviceSynchronize();
        cudaFree(hashes);
    }
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
                         kMer_t *kMers, kMer_t *hashes) {
    size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = blockDim.x * gridDim.x;
    for (size_t i = index; i < totalKMers; i += stride) {
        size_t hashIndex = i * n;
        kMer_t currentHash = kMers[index];
        currentHash = (currentHash ^ (currentHash >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        currentHash = (currentHash ^ (currentHash >> 27)) * UINT64_C(0x94d049bb133111eb);
        currentHash = currentHash ^ (currentHash >> 31);
        hashes[hashIndex++] = currentHash;
        for (size_t j = 1; j < n; j++) {
            currentHash = ((currentHash >> ROTATE_BITS)
                          | (currentHash << (KMER_BITS - ROTATE_BITS)))
                            ^ 0xABCD32108475AC38;
            hashes[hashIndex++] = currentHash;
        }
    }
}

__global__ void calcSketch(const size_t numReads, const size_t currentRead,
                           const size_t numKMers, const size_t n,
                           kMer_t *hashes, kMer_t *sketches) {
    size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = blockDim.x * gridDim.x;
    for (size_t i = index; i < numReads; i += stride) {
        size_t sketchIndex = (index + currentRead) * n;
        for (size_t j = 0; j < n; ++j) {
            size_t hashIndex = index * n * numKMers + j;
            kMer_t currentMin = ~(kMer_t) 0;
            for (size_t l = 0; l < numKMers; ++l) {
                kMer_t temp = hashes[hashIndex];
                hashIndex += n;
                currentMin = currentMin < temp ? currentMin : temp;
            }
            sketches[sketchIndex++] = currentMin;
        }
    }
}

NanoporeReads::~NanoporeReads() {
    cudaFree(this->sketches);
}

void NanoporeReads::printHashes() {
    for (size_t i = 0; i < this->numReads; ++i) {
        std::cout << this->readPos[i];
        for (size_t j = 0; j < this->n; ++j) {
            std::cout << ", \"" << this->sketches[i * this->n + j] << "\"";
        }
        std::cout << std::endl;
    }
}