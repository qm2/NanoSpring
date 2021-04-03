#ifndef BBHASHMAP_H_
#define BBHASHMAP_H_

#include "BooPHF.h"
#include "Types.h"
#include <string>
#include <vector>

// typedefs for the boomphf library
typedef boomphf::SingleHashFunctor<kMer_t>  hasher_t;
typedef boomphf::mphf<kMer_t, hasher_t> boophf_t; 

class BBHashMap {
    boophf_t *bphf = NULL; // BBHash from BooPHF.h
    size_t numKeys; // number of keys in hash table
    std::vector<kMer_t> keys; // keys in hash table
    // Note that we need to explicitly store keys since the BBHash doesn't
    std::vector<read_t> startPosInReadIds; 
    // for each key hash, store start position in readIds vector
    // This is more efficient than using a vector of vectors
    // Length will be numKeys
    std::vector<read_t> readIds;
    // vector storing the actual read ids for each key
    // Length will be numReads
public:
    /** 
     * @brief Construct the hash table using precomputed sketches as generated 
     * in MinHashReadFilter::initialize. 
     *
     * @param n number of sketches
     * @param numReads
     * @param sketches 
     * @param j sketch number we are working on
     * @param tempDir temporary directory
     */
    void initialize(const size_t n, const read_t numReads, const kMer_t *sketches,
                    const size_t j, const std::string &tempDir);

    /**
     * @brief insert all read ids corresponding to key to end of matches
     *
     * @param key
     * @param matches
     */
    void pushMatchesInVector(const kMer_t key, std::vector<read_t> &matches);

    BBHashMap() {
        bphf = NULL;
    }
    ~BBHashMap() {
        if (bphf != NULL) delete bphf;
    }
};

#endif // BBHASHMAP_H_
