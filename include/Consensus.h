#ifndef Z_CONSENSUS_CUH
#define Z_CONSENSUS_CUH

#include "ConsensusGraph.h"
#include "Contig.h"
#include "Edits.h"
#include "ReadData.h"
#include <memory>
#include <string>
#include <vector>

class Consensus {
public:
    ReadData *rD;

    ReadFilter *rF;

    /***
     * This is the read aligner used to align the reads as they are being added
     * to the Contig. This is only a rough aligner
     */
    ReadAligner *rA;

    /**
     * @brief The detailed aligner that aligns reads to the consensus senquence.
     *
     */
    StringAligner *aligner;

    const static POS_T END_SYMBOL = ~(POS_T)0l;

    /**
     * @brief Directory for storing the temp files (.genome, .pos, .type)
     *
     */
    std::string tempDir = "tempRaw/";

    /**
     * @brief Directory for storing the compressed temp files
     *
     */
    std::string compressedTempDir = "tempCompressed/";

    void generateConsensus();

    void writeConsensus();

    ~Consensus();

private:
    /**
     * @brief Whether the reads have been added to a graph
     * !!vector<bool> is not thread safe!!
     */
    std::vector<bool> inGraph;
    /**
     * @brief Index of the first read that has not been added to any graph
     *
     */
    size_t firstUnaddedRead;
    size_t numReads;

    std::vector<std::unique_ptr<ConsensusGraph>> graphs;

    /**
     * @brief Initializes firstUnaddedRead and numReads
     *
     */
    void initialize();

    /**
     * @brief Whether there are still reads that ar not in any graph
     *
     * @return true
     * @return false
     */
    bool hasReadsLeft();

    /**
     * @brief Gets an unadded read and updates its status to added
     *
     * @param read
     * @return true
     * @return false
     */
    bool getRead(size_t &read);

    /**
     * @brief Sets a read to unadded
     *
     * @param read
     * @return true
     * @return false
     */
    void putReadBack(size_t read);

    /**
     * @brief Sets an unadded read to added
     *
     * @param read
     * @return true
     * @return false
     */
    bool addRead(size_t read);

    ConsensusGraph *createGraph();
};

#endif // Z_CONSENSUS_CUH