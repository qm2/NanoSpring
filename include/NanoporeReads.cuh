#ifndef EXPERIMENTS_NANOPOREREADS_H
#define EXPERIMENTS_NANOPOREREADS_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <functional>
#include <memory>

#define KMER_BITS 64
#define ROTATE_BITS 1

typedef uint64_t kMer_t;

class NanoporeReads {
public:
    /***
     * Initializes from a file of reads with arguments k and n.
     * @param fileName Name of file to read
     * @param k k-mer
     * @param n Size of sketch
     */
    NanoporeReads(char *fileName, int k, int n);

    ~NanoporeReads();

    void calculateMinHashSketches();

    void printHashes();

private:
    const size_t k, n;
    unsigned long numReads;
    unsigned long readLen;
    kMer_t *sketches;

    std::vector <std::unique_ptr<std::string>> readData;
    std::vector<unsigned long> readPos;
    std::vector <std::unique_ptr<std::string>> editStrings;

    /***
     * Turns a k-mer in string format to an int
     * @param s k-mer
     * @return int representing k-mer
     */
    static kMer_t kMerToInt(const std::string &s);

    static char baseToInt(const char base);



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

__global__ void hashKMer(const size_t totalKMers, const size_t n,
                         kMer_t *kMers, kMer_t *hashes);

__global__ void calcSketch(const size_t numReads, const size_t currentRead,
                           const size_t numKMers, const size_t n,
                           kMer_t *hashes, kMer_t *sketches);

#endif //EXPERIMENTS_NANOPOREREADS_H
