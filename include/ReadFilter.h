#ifndef DE8B470D_E48E_48C7_809A_B4368CEC172C
#define DE8B470D_E48E_48C7_809A_B4368CEC172C

#include "ReadData.h"
#include <map>
#include <string>
#include <vector>

class filterStats {
public:
    const unsigned int overlapBaseThreshold, overlapSketchThreshold;
    unsigned int totalPositive, totalNegative, numOverlaps, numDisjoint,
        falsePositives, falseNegatives;

    filterStats(unsigned int overlapBaseThreshold,
                unsigned int overlapSketchThreshold);

    friend std::ostream &operator<<(std::ostream &out, const filterStats &o);
};

/**
 * @brief Filters reads that are likely to overlap
 *
 */
class ReadFilter {
public:
    virtual void getFilteredReads(size_t readToFind,
                                  std::vector<size_t> &results) = 0;

    virtual void getFilteredReads(const std::string &s,
                                  std::vector<size_t> &results) = 0;

    virtual void initialize(ReadData &rD) = 0;

    // For testing
    virtual filterStats getFilterStats(size_t overlapBaseThreshold,
                                       size_t overlapSketchThreshold) = 0;

    virtual ~ReadFilter();
};

class MinHashReadFilter : public ReadFilter {
public:
    // We store all the k-mers as uint64s. This would work for all k<=32,
    // which is definitely sufficient
    typedef uint64_t kMer_t;

    size_t k; // k-mer
    size_t n; // size of sketch
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
    void getFilteredReads(size_t readToFind,
                          std::vector<size_t> &results) override;

    void getFilteredReads(const std::string &s,
                          std::vector<size_t> &results) override;

    /**
     * @brief Finds out how will the filter does with the given parameters
     *
     * @param overlapBaseThreshold
     * @param overlapSketchThreshold
     * @return filterStats
     */
    filterStats getFilterStats(size_t overlapBaseThreshold,
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
    size_t readLen;
    std::vector<unsigned long> *readPos;
    std::vector<unsigned long> *readPosSorted;
    size_t numReads;
    kMer_t *sketches = nullptr;
    std::vector<std::map<kMer_t, std::vector<size_t>>> hashTables;
    kMer_t *randNumbers = nullptr;

    /**
     * @brief Initializes the hash tables given the data calculated in
     * sketches
     *
     */
    void populateHashTables();

    /**
     * @brief Get the Filtered Reads likely to overlap with the read represented
     * by sketch
     *
     * @param sketch
     * @param results
     */
    void getFilteredReads(kMer_t *sketch, std::vector<size_t> &results);

    /**
     * @brief Calculates "MinHash" sketches based on the hashes provided.
     *
     * @param numReads The number of reads to calculate sketches for.
     * @param numKMers The number of kMers that are calculated per read
     * @param n The size of the sketch
     * @param hashes The (numReads * numKMers * n) hashes. Indexed by
     * (readId, kMerId, hashId)
     * @param sketches The location to store the sketches. Should have enough
     * space to store (numReads * n) kMer_ts. Indexed by (readId, id in sketch)
     */
    static void calcSketch(const size_t numKMers, const size_t n,
                           kMer_t *hashes, kMer_t *sketches);
};

#endif /* DE8B470D_E48E_48C7_809A_B4368CEC172C */