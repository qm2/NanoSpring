#ifndef Z_CONSENSUS_CUH
#define Z_CONSENSUS_CUH

#include "ConsensusGraph.h"
#include "Contig.h"
#include "Edits.h"
#include "OmpMutex.h"
#include "ReadData.h"
#include "StringAligner.h"
#include "omp.h"
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>


// structure for stats 
struct CountStats {
    size_t countMinHash; // number of reads passing minhash filter
    size_t countMinHashNotInGraph; // number of reads passing minhash filter that are not already in graph
    size_t countMergeSort; // number of reads passing merge sort
    size_t countAligner; // number of reads passing aligner
    CountStats() {
        countMinHash = countMinHashNotInGraph = countMergeSort = countAligner = 0;   
    }
    CountStats operator+(const CountStats &c) {
        CountStats cs;
        cs.countMinHash = countMinHash + c.countMinHash;
        cs.countMinHashNotInGraph = countMinHashNotInGraph + c.countMinHashNotInGraph;
        cs.countMergeSort = countMergeSort + c.countMergeSort;
        cs.countAligner = countAligner + c.countAligner;
        return cs;
    }
};

class Consensus {
public:
    ReadData *rD;

    ReadFilter *rF;

    /**
     * This is the read aligner used to align the reads as they are being added
     * to the Contig. This is only a rough aligner
     */
    ReadAligner *rA;

    // Directory for storing the temp files (.genome, .pos, .type, .id metaData)
    std::string tempDir = "tempRaw/";

    std::string tempFileName = "Contig";

    int numThr;

    //parameters for minimap2
    size_t  m_k, m_w, max_chain_iter, edge_threshold;

    /**
     * @brief Generates consensus, calls writeReads and writeMainPath on each 
     * of the consensus graphs, and combines their output
     */
    void generateAndWriteConsensus();

    /**
     * @brief Combine files from threads and write metadata file
     *
     * @param numReadsInContig Vector with number of reads per contig in each thread
     */
    void finishWriteConsensus(const std::vector<std::vector<read_t>>& numReadsInContig);

    /**
     * @brief Checks that the read read in cG is equal to the read in ReadData
     * rD
     *
     * Assumes read is contained inside cG
     * @param cG
     * @param read
     * @return true
     * @return false
     */
    bool checkRead(ConsensusGraph *cG, read_t read);

    Consensus();

private:
    /**
     * @brief vector for checking repetitives
     *
     */
    std::vector<uint8_t> isRepetitive;
    /**
     * @brief Whether the reads have been added to a graph
     *
     */
    std::vector<uint8_t> inGraph;

    /**
     * @brief Number of locks to use (read i protected by lock i%numLocks)
     * Use large enough number of locks to avoid lock contention
     */
    static const uint32_t numLocks = (1<<24); // 16777216 locks

    /**
     * @brief Protects inGraph
     *
     * Note this is an array of locks, with read i protected by lock i%numLocks
     */
    std::vector<OmpMutex> readStatusLock;

    read_t numReads;

    /**
     * @brief Initializes numReads, readStatusLock
     *
     */
    void initialize();

    /**
     * @brief Gets an unadded read and updates its status to added
     *
     * @param read - starting point for searching unadded reads
     * @return true
     * @return false
     */
    bool getRead(read_t &read);

    /**
     * @brief Add reads overlapping with the mainPath of cG at position curPos
     *
     * @param cG
     * @param curPos
     * @param len The length of the mainPath used to get filtered reads
     * @param cs for collecting stats
     * @param logfile for writing pass/fail info to log
     * @param contigId current contig id for logging purposes
     */
    void addRelatedReads(ConsensusGraph *cG, ssize_t curPos, int len, CountStats &cs, std::ofstream &logfile, int contigId);

    /**
     * @brief Create consensus graph by picking previously unadded read
     *
     * @param firstUnaddedRead starting point for searching unadded reads
     */
    ConsensusGraph *createGraph(read_t &firstUnaddedRead);

    /**
     * @brief check if a read string is reptitive or not
     *
     * @param readID the read id
     */
    bool checkRepetitive(read_t readID);
};

#endif // Z_CONSENSUS_CUH
