#ifndef DE8B470D_E48E_48C7_809A_B4368CEC172C
#define DE8B470D_E48E_48C7_809A_B4368CEC172C

#include "ReadData.h"
#include <map>
#include <string>
#include <vector>

class FilterStats {
public:
    const unsigned int overlapBaseThreshold, overlapSketchThreshold;
    unsigned int totalPositive, totalNegative, numOverlaps, numDisjoint,
        falsePositives, falseNegatives;

    FilterStats(unsigned int overlapBaseThreshold,
                unsigned int overlapSketchThreshold);

    friend std::ostream &operator<<(std::ostream &out, const FilterStats &o);
};

/**
 * @brief Filters reads that are likely to overlap
 *
 */
class ReadFilter {
public:
    /**
     * @brief Rechieves reads that are likely to overlap with readToFind and
     * stores them in results
     *
     * @param readToFind id to the read we want to find overlapping reads with
     * @param results Vector where we store the results
     */
    virtual void getFilteredReads(read_t readToFind,
                                  std::vector<read_t> &results) = 0;

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

    /** For testing **/
    virtual FilterStats getFilterStats(size_t overlapBaseThreshold,
                                       size_t overlapSketchThreshold) = 0;

    virtual ~ReadFilter();
};

class MinHashReadFilter : public ReadFilter {
public:
    /** We store all the k-mers as uint64s. This would work for all k<=32,
     which is definitely sufficient **/
    typedef uint64_t kMer_t;
    static const size_t ROTATE_BITS = 13;
    static const uint64_t HASH_C64 = 11400714819323198485ULL;
    static const uint32_t HASH_C32 = 2654435769L;
    static const size_t KMER_BITS = 64;

    /** [k]-mer **/
    size_t k;
    /** size of sketch **/
    size_t n;
    size_t overlapSketchThreshold;

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
     */
    void string2Sketch(const std::string &s, kMer_t *sketch);

    /**
     * @brief Calculates n hashes of kMer and stores them in hashes
     *
     * @tparam OutputIt An output iterator to store the hashes
     * @param n
     * @param kMer
     * @param hashes
     * @param randNumbers n kMer_t random numbers
     */
    template <class OutputIt> void hashKMer(const kMer_t kMer, OutputIt hashes);

    /**
     * @brief Gets reads probably overlapping with readToFind and stores them in
     * results
     *
     * @param readToFind
     * @param results
     */
    void getFilteredReads(read_t readToFind,
                          std::vector<read_t> &results) override;

    void getFilteredReads(const std::string &s,
                          std::vector<read_t> &results) override;

    /**
     * @brief Finds out how will the filter does with the given parameters
     *
     * @param overlapBaseThreshold
     * @param overlapSketchThreshold
     * @return FilterStats
     */
    FilterStats getFilterStats(size_t overlapBaseThreshold,
                               size_t overlapSketchThreshold) override;

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
     * @tparam OutputIt the output iterator to store the kmers
     * @param s
     * @param k
     * @param kMers
     */
    template <class OutputIt>
    static void string2KMers(const std::string &s, const size_t k,
                             OutputIt kMers);

private:
    ReadData *rD;
    std::vector<size_t> *readPos;
    std::vector<std::pair<size_t, read_t>> readPosSorted;
    read_t numReads;
    kMer_t *sketches = nullptr;
    std::vector<std::map<kMer_t, std::vector<read_t>>> hashTables;
    kMer_t *randNumbers = nullptr;

    /**
     * @brief Initializes the hash tables given the data calculated in
     * sketches
     *
     */
    void populateHashTables();

    /**
     * @brief Get the Filtered Reads likely to overlap with the read represented
     * by sketch and stores them in results
     *
     * @param sketch
     * @param results
     */
    void getFilteredReads(kMer_t sketch[], std::vector<read_t> &results);

    /**
     * @brief Calculates "MinHash" sketches based on the hashes provided.
     *
     * @tparam InputIt an input random access iterator that dereferences to
     * kMer_t
     * @tparam OutputIt an output random access iterator that dereferences to
     * kMer_t
     * @param numReads The number of reads to calculate sketches for.
     * @param numKMers The number of kMers that are calculated per read
     * @param n The size of the sketch
     * @param hashes The (numReads * numKMers * n) hashes. Indexed by
     * (readId, kMerId, hashId)
     * @param sketches The location to store the sketches. Should have enough
     * space to store (numReads * n) kMer_ts. Indexed by (readId, id in sketch)
     */
    template <typename InputIt, typename OutputIt>
    static void calcSketch(const size_t numKMers, const size_t n,
                           InputIt hashes, OutputIt sketches);
};

/******************************************************************************/
/* public template methods */

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

/******************************************************************************/

#endif /* DE8B470D_E48E_48C7_809A_B4368CEC172C */