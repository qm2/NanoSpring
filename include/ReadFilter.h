#ifndef DE8B470D_E48E_48C7_809A_B4368CEC172C
#define DE8B470D_E48E_48C7_809A_B4368CEC172C

#include "ReadData.h"
#include "Types.h"
#include "BBHashMap.h"
#include <map>
#include <string>
#include <vector>

/**
 * @brief Filters reads that are likely to overlap
 *
 */
class ReadFilter {
public:
    /**
     * @brief Rechieves reads that are likely to overlap with string s and
     * stores them in results
     *
     * @param s string to find overlapping reads with
     * @param results Vector where we store the results
     */
    virtual void getFilteredReads(const std::string &s,
                                  std::vector<read_t> &results) = 0;

    virtual void initialize(ReadData &rD) = 0;

    virtual ~ReadFilter();
};


class MinHashReadFilter : public ReadFilter {
public:

    /** [k]-mer **/
    size_t k;
    /** size of sketch **/
    size_t n;
    size_t overlapSketchThreshold;

    std::string tempDir;

    /**
     * @brief Builds the hash tables from the data in rD
     *
     * @param rD
     */
    void initialize(ReadData &rD) override;

    /**
     * @brief Generates a sequence of n kMer_t random numbers
     *
     * @param n
     */
    void generateRandomNumbers(size_t n);

    /**
     * @brief Converts a string to a sketch and stores it in sketch
     *
     * @param s
     * @param sketch
     * @param kMers preallocated for speed during multithreaded initialization 
     * (size should be at least s.size()-k+1)
     * @param hashes preallocated for speed during multithreaded initialization
     * (size should be at least n*(s.size()-k+1))
     */
    void string2Sketch(const std::string &s, kMer_t *sketch, 
                    std::vector<kMer_t> &kMers, std::vector<kMer_t> &hashes);

    /**
     * @brief Calculates n hashes of kMer and stores them in hashes
     *
     * @param kMer
     * @param hashes
     */
    void hashKMer(const kMer_t kMer, std::vector<kMer_t> &hashes);

    void getFilteredReads(const std::string &s,
                          std::vector<read_t> &results) override;


    MinHashReadFilter();

    ~MinHashReadFilter() override;

    /**
     * Turns a k-mer in string format to an int
     * @param s k-mer
     * @return int representing k-mer
     */
    static kMer_t kMerToInt(const std::string &s);

    /**
     * @brief Converts a base, one of 'A', 'T', 'C', 'G' into a two-bit int
     *
     * @param base
     * @return char
     */
    static char baseToInt(const char base);

    /**
     * @brief Converts a string to kmers and stores it in kMers
     *
     * @param s
     * @param k
     * @param kMers
     */
    static void string2KMers(const std::string &s, const size_t k,
                             std::vector<kMer_t> &kMers);

private:
    ReadData *rD;
    std::vector<size_t> *readPos;
    std::vector<std::pair<size_t, read_t>> readPosSorted;
    read_t numReads;
    std::vector<BBHashMap> hashTables; // vector of size n to store the n hash maps
    kMer_t *randNumbers = nullptr;

    std::hash<kMer_t> hasher;

    /**
     * @brief Initializes the hash tables given the data calculated in
     * sketches
     *
     * @param sketches
     */
    void populateHashTables(const std::vector<kMer_t> &sketches);

    /**
     * @brief Get the Filtered Reads likely to overlap with the read represented
     * by sketch and stores them in results
     *
     * @param sketch
     * @param results
     */
    void getFilteredReads(kMer_t sketch[], std::vector<read_t> &results);

};


#endif /* DE8B470D_E48E_48C7_809A_B4368CEC172C */
