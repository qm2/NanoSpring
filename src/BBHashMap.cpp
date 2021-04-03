#include "BBHashMap.h"
#include "BooPHF.h"
#include "ReadFilter.h"
#include <omp.h>
#include <fstream>
#include <string>
#include <algorithm> // sort
#include <cstdio> // remove

void BBHashMap::initialize(const size_t n, const read_t numReads, const kMer_t *sketches,
                           const size_t j, const std::string &tempDir) {

    int tid = omp_get_thread_num(); // for naming temporary files without collision
    // create temporary array for various operations
    kMer_t *arr = new kMer_t[numReads];
    // Step 1: collect the sketches for hash table j into arr
    size_t currentIndex = j;
    for (read_t r = 0; r < numReads; r++) {
        arr[r] = sketches[currentIndex];
        currentIndex += n;
    }
    // Step 2: write arr to file (to save memory)
    std::string keysFile = tempDir + std::string("/keys.bin.") + std::to_string(tid);
    std::ofstream foutkey(keysFile, std::ios::binary);
    for (read_t r = 0; r < numReads; r++)
        foutkey.write((char*)&arr[r], sizeof(kMer_t));
    foutkey.close();

    // deduplicating arr
    std::sort(arr, arr + numReads);
    read_t cur = 0;
    for (read_t r = 1; r < numReads; r++)
      if (arr[r] != arr[cur]) arr[++cur] = arr[r];

    numKeys = cur + 1; // number of unique keys

    // Step 3: construct the minimal perfect hash
    auto data_iterator =
        boomphf::range(static_cast<const kMer_t *>(arr),
                       static_cast<const kMer_t *>(arr + numKeys));
    double gammaFactor = 5.0;  // balance between speed and memory
    bphf = new boomphf::mphf<kMer_t, hasher_t>(
        numKeys, data_iterator, 1, gammaFactor, false, false);

    delete[] arr; // not needed anymore

    // Step 4: compute hash values for all reads using the hash table
    // and write to temporary file. In addition, we store a mapping from the
    // hash value (0..numKeys) to the actual value of the key to help during lookup
    keys.resize(numKeys);
    std::string hashesFile = tempDir + std::string("/hashes.bin.") + std::to_string(tid);
    std::ofstream fouthash(hashesFile, std::ios::binary);
    std::ifstream finkey(keysFile, std::ios::binary);
    kMer_t currentkey;
    uint64_t currentHash;
    for (read_t r = 0; r < numReads; r++) {
        finkey.read((char *)&currentkey, sizeof(kMer_t));
        currentHash = bphf->lookup(currentkey);
        keys[currentHash] = currentkey;
        fouthash.write((char *)&currentHash, sizeof(uint64_t));
    }
    finkey.close();
    remove(keysFile.c_str());
    fouthash.close();

    // Step 5: Next we define the location array to store the start position in readIds 
    // array for each key 
    // fill startPosInReadIds by first storing counts for each key and then
    // doing cumulative sum
    startPosInReadIds.resize(numKeys);
    std::ifstream finhash(hashesFile, std::ios::binary);   
    for (read_t r = 0; r < numReads; r++) {
        finhash.read((char *)&currentHash, sizeof(uint64_t));
        startPosInReadIds[currentHash]++;
    }
    finhash.close();
    read_t cumSum = startPosInReadIds[0], prevCumSum;
    startPosInReadIds[0] = 0;
    for (read_t i = 1; i < numKeys; i++) {
        prevCumSum = cumSum;
        cumSum += startPosInReadIds[i];
        startPosInReadIds[i] = prevCumSum;
    }
    
    // Step 6: Finally we insert the read ids for each key in the appropriate place within 
    // the readIds array. For convenience we use a vector numFilled which denotes the 
    // number of reads already handled for each key
    std::vector<read_t> numFilled(numKeys,0);
    readIds.resize(numReads);
    finhash.open(hashesFile, std::ios::binary);
    for (read_t r = 0; r < numReads; r++) {
        finhash.read((char *)&currentHash, sizeof(uint64_t));
        readIds[startPosInReadIds[currentHash]+numFilled[currentHash]] = r;
        numFilled[currentHash]++;
    }
    finhash.close();
    remove(hashesFile.c_str());
    return;
}

void BBHashMap::pushMatchesInVector(const kMer_t key, std::vector<read_t> &matches) {
    auto hash = bphf->lookup(key);
    if (hash > numKeys)
        return; // no matches found
    if (keys[hash] != key)
        return;
    // no matches found (note bbhash can return false positives so we need to check)
    
    read_t startPos = startPosInReadIds[hash];
    read_t endPos; // endPos is either the end of vector if we are at the last key,
                   // or it is computed based on the next start pos
    if (hash == numKeys - 1)
        endPos = readIds.size();
    else
        endPos = startPosInReadIds[hash+1];
    
    // insert the matches to end of matches vector
    matches.insert(matches.end(),readIds.begin()+startPos,readIds.begin()+endPos);
    return;
}
