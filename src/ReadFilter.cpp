#include "ReadFilter.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>

#define ROTATE_BITS 13
#define HASH_C64 11400714819323198485ULL
#define HASH_C32 2654435769L
#define KMER_BITS 64

FilterStats::FilterStats(unsigned int overlapBaseThreshold,
                         unsigned int overlapSketchThreshold)
    : overlapBaseThreshold(overlapBaseThreshold),
      overlapSketchThreshold(overlapSketchThreshold), totalPositive(0),
      totalNegative(0), numOverlaps(0), numDisjoint(0), falsePositives(0),
      falseNegatives(0) {}

std::ostream &operator<<(std::ostream &out, const FilterStats &o) {
    const int w = 13;
    out << std::setw(w) << "overlapBaseTh"
        << "," << std::setw(w) << "numKMerTh"
        << "," << std::setw(w) << "totalPos"
        << "," << std::setw(w) << "totalNeg"
        << "," << std::setw(w) << "numOverlaps"
        << "," << std::setw(w) << "numDisjoint"
        << "," << std::setw(w) << "falsePos"
        << "," << std::setw(w) << "falseNeg" << std::endl;
    out << std::setw(w) << o.overlapBaseThreshold << "," << std::setw(w)
        << o.overlapSketchThreshold << "," << std::setw(w) << o.totalPositive
        << "," << std::setw(w) << o.totalNegative << "," << std::setw(w)
        << o.numOverlaps << "," << std::setw(w) << o.numDisjoint << ","
        << std::setw(w) << o.falsePositives << "," << std::setw(w)
        << o.falseNegatives << std::endl;
    return out;
}

ReadFilter::~ReadFilter() {}

void MinHashReadFilter::initialize(ReadData &rD) {
    this->rD = &rD;
    numReads = rD.getNumReads();
    readLen = rD.getRead(0).length();
    readPos = &rD.getReadPos();
    readPosSorted = &rD.getReadPosSorted();

    if (sketches)
        free(sketches);
    sketches = (kMer_t *)malloc(n * numReads * sizeof(kMer_t));

    generateRandomNumbers(n);

#pragma omp parallel for
    for (read_t i = 0; i < numReads; ++i)
        string2Sketch(rD.getRead(i), sketches + i * n);

    populateHashTables();
}

void MinHashReadFilter::generateRandomNumbers(size_t n) {
    if (randNumbers)
        free(randNumbers);
    randNumbers = (kMer_t *)malloc(n * sizeof(kMer_t));

    std::random_device rd;
    std::mt19937_64 gen(rd());

    /* This is where you define the number generator for unsigned long long: */
    std::uniform_int_distribution<unsigned long long> dis;

    for (size_t i = 0; i < n; ++i) {
        randNumbers[i] = dis(gen);
    }
}

void MinHashReadFilter::getFilteredReads(read_t readToFind,
                                         std::vector<read_t> &results) {
    results.clear();
    getFilteredReads(sketches + readToFind * n, results);
}

void MinHashReadFilter::getFilteredReads(kMer_t sketch[],
                                         std::vector<read_t> &results) {
    std::vector<read_t> matches;
    results.clear();
    for (size_t sketchIndex = 0; sketchIndex < n; ++sketchIndex) {
        kMer_t curHash = sketch[sketchIndex];
        // std::vector<read_t> &m = hashTables[sketchIndex].at(curHash);
        auto m = hashTables[sketchIndex].find(curHash);
        if (m == hashTables[sketchIndex].end())
            continue;
        matches.insert(matches.end(), m->second.begin(), m->second.end());
    }
    std::sort(matches.begin(), matches.end());
    auto end = matches.end();
    auto next = matches.begin();
    for (auto it = matches.begin(); it != end; it = next) {
        next = std::upper_bound(it, end, *it);

        if (next - it >= (ssize_t)overlapSketchThreshold) {
            results.push_back(*it);
        }
    }
}

void MinHashReadFilter::getFilteredReads(const std::string &s,
                                         std::vector<read_t> &results) {
    results.clear();
    kMer_t sketch[n];
    string2Sketch(s, sketch);
    getFilteredReads(sketch, results);
}

MinHashReadFilter::MinHashReadFilter() {}

FilterStats MinHashReadFilter::getFilterStats(size_t overlapBaseThreshold,
                                              size_t overlapSketchThreshold) {
    FilterStats result(overlapBaseThreshold, overlapSketchThreshold);
    // First we calculate the number of overlaps and disjoints

    size_t numOverlaps = 0;
#pragma omp parallel for reduction(+ : numOverlaps)
    for (read_t i = 0; i < numReads; ++i) {
        long curTh = (*readPosSorted)[i] + readLen - overlapBaseThreshold;
        for (read_t j = i + 1; j < numReads; ++j) {
            if (((long)(*readPosSorted)[j]) <= curTh)
                numOverlaps++;
            else
                break;
        }
    }
    result.numOverlaps = numOverlaps;
    result.numDisjoint = (numReads * (unsigned long long)(numReads - 1)) / 2 -
                         result.numOverlaps;

    // Now we calculate falsePositives, falseNegatives, etc
    size_t totalPositive = 0;
    size_t falsePositives = 0;
#pragma omp parallel for reduction(+ : totalPositive, falsePositives)
    for (read_t i = 0; i < numReads; ++i) {
        std::multiset<read_t> matches;
        unsigned long curPos = (*readPos)[i];
        long th = readLen - overlapBaseThreshold;
        for (size_t sketchIndex = 0; sketchIndex < n; ++sketchIndex) {
            kMer_t curHash = sketches[i * n + sketchIndex];
            //            std::cout << i << " " << sketchIndex << " " << curHash
            //            << std::endl; auto currentMap =
            //            hashTables[sketchIndex]; for (auto p : currentMap) {
            //                std::cout << p.first << " ";
            //            }
            //            std::cout << std::endl;
            std::vector<read_t> &m = hashTables[sketchIndex].at(curHash);
            matches.insert(m.begin(), m.end());
        }
        auto end = matches.end();
        for (auto it = matches.begin(); it != end;
             it = matches.upper_bound(*it)) {
            if (*it <= i)
                continue;
            unsigned long pos = (*readPos)[*it];
            if (matches.count(*it) >= overlapSketchThreshold) {
                if (abs((long)pos - (long)curPos) > th)
                    falsePositives++;
                totalPositive++;
            }
        }
        //        std::cout << std::endl;
    }

    result.falsePositives = falsePositives;
    result.totalPositive = totalPositive;

    result.totalNegative =
        result.numOverlaps + result.numDisjoint - result.totalPositive;
    result.falseNegatives =
        result.numOverlaps - result.totalPositive + result.falsePositives;
    return result;
}

MinHashReadFilter::kMer_t MinHashReadFilter::kMerToInt(const std::string &s) {
    size_t l = s.length();
    kMer_t result = 0;
    for (size_t i = 0; i < l; ++i) {
        result <<= 2;
        result |= baseToInt(s[i]);
    }
    return result;
}

// Using the bit operations version of this function provides a 13X improvement
// in speed
char MinHashReadFilter::baseToInt(const char base) {
    return (base & 0b10) | ((base & 0b100) >> 2);
}

template <class OutputIt>
void MinHashReadFilter::string2KMers(const std::string &s, const size_t k,
                                     OutputIt kMers) {
    ssize_t maxI = s.length() - k + 1;
    if (maxI <= 0)
        return;
    kMer_t currentKMer = kMerToInt(s.substr(0, k));
    *kMers = currentKMer;
    kMers++;
    const unsigned long long mask = (1ull << (2 * k)) - 1;
    for (size_t i = 1; i < (size_t)maxI; ++i) {
        currentKMer =
            ((currentKMer << 2) | MinHashReadFilter::baseToInt(s[i + k - 1])) &
            mask;
        *kMers = currentKMer;
        kMers++;
    }
}

template <class OutputIt>
void MinHashReadFilter::hashKMer(kMer_t kMer, OutputIt hashes) {
    kMer_t currentHash = kMer;
    currentHash = (currentHash * (uint64_t)HASH_C64);
    currentHash ^= randNumbers[0];
    *hashes = currentHash;
    ++hashes;
    for (size_t l = 1; l < n; l++) {
        kMer_t newHash = ((currentHash >> ROTATE_BITS) |
                          (currentHash << (KMER_BITS - ROTATE_BITS))) ^
                         0xABCD32108475AC38;
        newHash = (newHash * (uint64_t)HASH_C64);
        newHash ^= randNumbers[l];
        newHash += currentHash;
        currentHash = newHash;
        *hashes = currentHash;
        ++hashes;
    }
}

void MinHashReadFilter::calcSketch(const size_t numKMers, const size_t n,
                                   kMer_t *hashes, kMer_t *sketches) {
    // #pragma omp parallel for
    for (size_t hashIndex = 0; hashIndex < n; ++hashIndex) {
        kMer_t currentMin = ~(kMer_t)0;
        //#pragma omp parallel for reduction(min:currentMin)
        for (size_t l = 0; l < numKMers * n; l += n) {
            kMer_t temp = hashes[hashIndex + l];
            currentMin = std::min(currentMin, temp);
        }
        sketches[hashIndex] = currentMin;
    }
}

void MinHashReadFilter::string2Sketch(const std::string &s, kMer_t *sketch) {
    size_t numKMers = s.length() - k + 1;
    kMer_t kMers[numKMers];
    string2KMers(s, k, kMers);
    kMer_t hashes[numKMers * n];
    for (size_t i = 0; i < numKMers; ++i)
        hashKMer(kMers[i], hashes + i * n);
    calcSketch(numKMers, n, hashes, sketch);
}

MinHashReadFilter::~MinHashReadFilter() {
    if (randNumbers)
        free(randNumbers);
    if (sketches)
        free(sketches);
}

void MinHashReadFilter::populateHashTables() {
    std::cout << "Starting to populate hash tables" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < n; ++i) {
        hashTables.push_back(std::map<kMer_t, std::vector<read_t>>());
    }
#pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        std::map<kMer_t, std::vector<read_t>> &hT = hashTables[i];
        size_t currentIndex = i;
        for (read_t j = 0; j < numReads; ++j) {
            hT[sketches[currentIndex]].push_back(j);
            //            std::cout << "Inserting " << j << " to " <<
            //            sketches[currentIndex] << std::endl;
            currentIndex += n;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "finished populating hash tables" << std::endl;
    std::cout << duration.count() << " milliseconds passed" << std::endl;
}
