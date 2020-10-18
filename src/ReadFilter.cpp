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

#pragma omp parallel for
    for (read_t i = 0; i < numReads; ++i)
        string2Sketch(rD.getRead(i), sketches + i * n);

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
        read_t realI = readPosSorted[i].second;
        for (read_t j = i + 1; j < numReads; ++j) {
            read_t realJ = readPosSorted[j].second;
            size_t minEnd =
                std::min((*readPos)[realI] + rD->getRead(realI).size(),
                         (*readPos)[realJ] + rD->getRead(realJ).size());
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
                ssize_t minEnd =
                    std::min((*readPos)[i] + rD->getRead(i).size(),
                             (*readPos)[*it] + rD->getRead(*it).size());
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

template <typename InputIt, typename OutputIt>
void MinHashReadFilter::calcSketch(const size_t numKMers, const size_t n,
                                   InputIt hashes, OutputIt sketches) {
    // This code has been optimized such that
    // 1. we are accessing memory in sequence
    // 2. We only try to increment pointers
    OutputIt sketchesEnd = sketches + n;
    for (OutputIt tempSketches = sketches; tempSketches != sketchesEnd;
         ++tempSketches)
        *tempSketches = ~(kMer_t)0;
    InputIt hashesEnd = hashes + numKMers * n;
    OutputIt tempSketches = sketches;
    for (InputIt tempHashes = hashes; tempHashes != hashesEnd;
         ++tempHashes, ++tempSketches) {
        tempSketches = tempSketches == sketchesEnd ? sketches : tempSketches;
        *tempSketches = std::min(*tempSketches, *tempHashes);
    }
}

void MinHashReadFilter::string2Sketch(const std::string &s, kMer_t *sketch) {
    ssize_t numKMers = s.length() - k + 1;
    if (numKMers < 0)
        return;
    std::vector<kMer_t> kMers(numKMers);
    string2KMers(s, k, kMers.begin());
    std::vector<kMer_t> hashes(numKMers * n);
    for (size_t i = 0; i < (size_t)numKMers; ++i)
        hashKMer(kMers[i], hashes.begin() + i * n);
    calcSketch(numKMers, n, hashes.begin(), sketch);
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
