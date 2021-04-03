#include "ReadFilter.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>

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
    readPos = &rD.getReadPos();
    readPosSorted.clear();
    for (read_t r = 0; r < numReads; ++r) {
        readPosSorted.push_back(std::make_pair((*readPos)[r], r));
    }
    std::sort(readPosSorted.begin(), readPosSorted.end());

    if (sketches)
        delete[] sketches;
    sketches = new kMer_t[n * numReads];

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
            string2Sketch(readStr, sketches + i * n, kMersVec, hashesVec);
        }
    } // pragma omp parallel

    populateHashTables();
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

FilterStats MinHashReadFilter::getFilterStats(size_t overlapBaseThreshold,
                                              size_t overlapSketchThreshold) {
    FilterStats result(overlapBaseThreshold, overlapSketchThreshold);
    // First we calculate the number of overlaps and disjoints

    size_t numOverlaps = 0;
#pragma omp parallel for reduction(+ : numOverlaps)
    for (read_t i = 0; i < numReads; ++i) {
        read_t realI = readPosSorted[i].second;
        for (read_t j = i + 1; j < numReads; ++j) {
            read_t realJ = readPosSorted[j].second;
            std::string readStr1, readStr2;
            rD->getRead(realI, readStr1);
            rD->getRead(realJ, readStr2);
            size_t minEnd =
                std::min((*readPos)[realI] + readStr1.size(),
                         (*readPos)[realJ] + readStr2.size());
            size_t maxBegin = std::max((*readPos)[realI], (*readPos)[realJ]);
            if (minEnd > maxBegin + overlapBaseThreshold)
                numOverlaps++;
            else if ((*readPos)[realJ] > minEnd)
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
        for (size_t sketchIndex = 0; sketchIndex < n; ++sketchIndex) {
            kMer_t curHash = sketches[i * n + sketchIndex];
            //            std::cout << i << " " << sketchIndex << " " << curHash
            //            << std::endl; auto currentMap =
            //            hashTables[sketchIndex]; for (auto p : currentMap) {
            //                std::cout << p << " ";
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
            if (matches.count(*it) >= overlapSketchThreshold) {
                std::string readStr1, readStr2;
                rD->getRead(i, readStr1);
                rD->getRead(*it, readStr2);
                ssize_t minEnd =
                    std::min((*readPos)[i] + readStr1.size(),
                             (*readPos)[*it] + readStr2.size());
                ssize_t maxBegin = std::max((*readPos)[i], (*readPos)[*it]);
                if (minEnd - maxBegin < (ssize_t)overlapBaseThreshold)
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
    if (sketches)
        delete[] sketches;
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
