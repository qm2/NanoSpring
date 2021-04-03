#include "ReadFilter.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>

ReadFilter::~ReadFilter() {}

void MinHashReadFilter::initialize(ReadData &rD) {
    this->rD = &rD;
    numReads = rD.getNumReads();
    readPos = &rD.getReadPos();
    readPosSorted.clear();
    for (read_t r = 0; r < numReads; ++r) {
        readPosSorted.push_back(std::make_pair((*readPos)[r], r));
    }
    std::sort(readPosSorted.begin(), readPosSorted.end());

    std::vector<kMer_t> sketches(n * numReads);

    generateRandomNumbers(n);

    size_t maxNumkMers;
    if (rD.maxReadLen < k-1)
        maxNumkMers = 0;
    else
        maxNumkMers = rD.maxReadLen - k + 1;

#pragma omp parallel 
    {
        // We define these vectors here rather than allocate inside string2Sketch
        // to avoid thread contention during repeated allocation and deallocation.
        // Note that memory allocation typically leads to waits when multiple threads
        // do it at the same time.
        std::vector<kMer_t> kMersVec(maxNumkMers), hashesVec(n);
        std::string readStr;
#pragma omp for
        for (read_t i = 0; i < numReads; ++i){
            rD.getRead(i, readStr);
            string2Sketch(readStr, sketches.data() + i * n, kMersVec, hashesVec);
        }
    } // pragma omp parallel

    populateHashTables(sketches);
}

void MinHashReadFilter::generateRandomNumbers(size_t n) {
    if (randNumbers)
        delete[] randNumbers;
    randNumbers = new kMer_t[n];

    std::random_device rd;
    std::mt19937_64 gen(rd());

    /* This is where you define the number generator for unsigned long long: */
    std::uniform_int_distribution<unsigned long long> dis;

    for (size_t i = 0; i < n; ++i) {
        randNumbers[i] = dis(gen);
    }
}

void MinHashReadFilter::getFilteredReads(kMer_t sketch[],
                                         std::vector<read_t> &results) {
    std::vector<read_t> matches;
    results.clear();
    for (size_t sketchIndex = 0; sketchIndex < n; ++sketchIndex) {
        kMer_t curHash = sketch[sketchIndex];
        hashTables[sketchIndex].pushMatchesInVector(curHash,matches);
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
    std::vector<kMer_t> sketch(n);
    size_t numKmers;
    if (s.size() < k - 1)
        numKmers = 0;
    else
        numKmers = s.size() - k + 1;
    std::vector<kMer_t> kMersVec(numKmers), hashesVec(n); // preallocation
    string2Sketch(s, sketch.data(), kMersVec, hashesVec);
    getFilteredReads(sketch.data(), results);
}

MinHashReadFilter::MinHashReadFilter() {}

kMer_t MinHashReadFilter::kMerToInt(const std::string &s) {
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

void MinHashReadFilter::string2Sketch(const std::string &s, kMer_t *sketch, std::vector<kMer_t> &kMers, std::vector<kMer_t> &hashes) {
    ssize_t numKMers = s.length() - k + 1;
    if (numKMers < 0)
        return;
    string2KMers(s, k, kMers);
    // set sketch to -1 (max possible)
    for (size_t j = 0; j < n; j++)
        sketch[j] = (kMer_t)(-1);
    // now compute sketch
    for (size_t i = 0; i < (size_t)numKMers; ++i) {
        hashKMer(kMers[i], hashes);
        for (size_t j = 0; j < n; j++)
            sketch[j] = std::min(hashes[j],sketch[j]);
    }
}

void MinHashReadFilter::hashKMer(kMer_t kMer, std::vector<kMer_t> &hashes) {
    for (size_t l = 0; l < n; l++)
        hashes[l] = hasher(kMer^randNumbers[l]);
}

void MinHashReadFilter::string2KMers(const std::string &s, const size_t k,
                                     std::vector<kMer_t> &kMers) {
    ssize_t maxI = s.length() - k + 1;
    if (maxI <= 0)
        return;
    kMer_t currentKMer = kMerToInt(s.substr(0, k));
    kMers[0] = currentKMer;
    const unsigned long long mask = (1ull << (2 * k)) - 1;
    for (size_t i = 1; i < (size_t)maxI; ++i) {
        currentKMer =
            ((currentKMer << 2) | MinHashReadFilter::baseToInt(s[i + k - 1])) &
            mask;
        kMers[i] = currentKMer;
    }
}

MinHashReadFilter::~MinHashReadFilter() {
    if (randNumbers)
        delete[] randNumbers;
}

void MinHashReadFilter::populateHashTables(const std::vector<kMer_t> &sketches) {
    std::cout << "Starting to populate hash tables" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    hashTables.resize(n);
#pragma omp parallel for
    for (size_t i = 0; i < n; ++i)
        hashTables[i].initialize(n,numReads,sketches.data(),i,tempDir);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "finished populating hash tables" << std::endl;
    std::cout << duration.count() << " milliseconds passed" << std::endl;
}
