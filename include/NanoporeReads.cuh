#ifndef EXPERIMENTS_NANOPOREREADS_H
#define EXPERIMENTS_NANOPOREREADS_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <map>
#include <gperftools/profiler.h>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <set>
#include <cmath>

#define KMER_BITS 64
#define ROTATE_BITS 13
#define HASH_C64 11400714819323198485ULL
#define HASH_C32 2654435769L
#define _GPU
#undef _GPU

typedef uint64_t kMer_t;

class filterStats {
public:
    unsigned int totalPositive, totalNegative, numOverlaps,
            numDisjoint, falsePositives, falseNegatives;
    const unsigned int overlapBaseThreshold, overlapSketchThreshold;


    filterStats(unsigned int overlapBaseThreshold, unsigned int overlapSketchThreshold);

    friend std::ostream &operator<<(std::ostream &out, const filterStats &o);
};

class NanoporeReads {
public:
    std::vector<std::unique_ptr<std::string>> readData;
    std::vector<unsigned long> readPos;
    std::vector<unsigned long> readPosSorted;
    std::vector<std::unique_ptr<std::string>> editStrings;
    unsigned long numReads;
    unsigned long readLen;
    const size_t k, n;

    kMer_t *sketches;
    std::vector<std::map<kMer_t, std::vector<size_t>>> hashTables;

    /***
     * Initializes from a file of reads with arguments k and n.
     * @param fileName Name of file to read
     * @param k k-mer
     * @param n Size of sketch
     */
    NanoporeReads(const char *fileName, int k, int n);

    ~NanoporeReads();

    void calculateMinHashSketches();

    void printHashes();

    filterStats getFilterStats(unsigned int overlapBaseThreshold, unsigned int overlapSketchThreshold);

    /***
     * Turns a k-mer in string format to an int
     * @param s k-mer
     * @return int representing k-mer
     */
    static kMer_t kMerToInt(const std::string &s);

    static char baseToInt(const char base);

private:


    void populateHashTables();



//    /**
//    * Changes "ATCG" into two bools
//    * @param c Base
//    * @return First bit in representation (ATCG: 00, 01, 10, 11)
//    */
//    static bool base2Bool0(char c);
//
//    /**
//    * Changes "ATCG" into two bools
//    * @param c Base
//    * @return Second bit in representation (ATCG: 00, 01, 10, 11)
//    */
//    static bool base2Bool1(char c);
};

#ifdef _GPU

__global__
#endif

void hashKMer(const size_t totalKMers, const size_t n,
              kMer_t *kMers, kMer_t *hashes, kMer_t *randNumbers);

/***
 * Kernel to calculate the sketch. If kMers is provided then the sequence is stored;
 * if not, the hash values are stored
 * @param numReads
 * @param currentRead
 * @param numKMers
 * @param n
 * @param hashes
 * @param sketches
 * @param kMers
 */
#ifdef _GPU

__global__
#endif

void calcSketch(const size_t numReads, const size_t currentRead,
                const size_t numKMers, const size_t n,
                kMer_t *hashes, kMer_t *sketches, kMer_t *kMers);

class ReadFilter {
public:
    virtual void getFilteredReads(size_t readToFind, std::vector<size_t> &results) = 0;
};

class MinHashReadFilter : public ReadFilter {
public:
    void getFilteredReads(size_t readToFind, std::vector<size_t> &results);

    MinHashReadFilter(size_t overlapSketchThreshold, NanoporeReads &nR);

private:
    size_t overlapSketchThreshold;
    NanoporeReads &nR;
};

#endif //EXPERIMENTS_NANOPOREREADS_H
