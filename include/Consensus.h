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

class Consensus {
public:
    ReadData *rD;

    ReadFilter *rF;

    /**
     * This is the read aligner used to align the reads as they are being added
     * to the Contig. This is only a rough aligner
     */
    ReadAligner *rA;

    /**
     * @brief The detailed aligner that aligns reads to the consensus senquence.
     *
     */
    ConsensusGraph::StringAligner_t *aligner;

    // Directory for storing the temp files (.genome, .pos, .type, .id metaData)
    std::string tempDir = "tempRaw/";

    std::string tempFileName = "Contig";

    void generateConsensus();

    /**
     * @brief Calls writeReads and writeMainPath on each of the consensus
     * graphs, and combines their otuput
     *
     */
    void writeConsensus();

    Consensus();

    ~Consensus();

private:
    /**
     * @brief Whether the reads have been added to a graph
     *
     * !!vector<bool> is not thread safe!!
     */
    std::vector<bool> inGraph;
    /**
     * @brief Index of the first read that has not been added to any graph
     *
     */
    read_t firstUnaddedRead;
    /**
     * @brief Protects inGraph and firstUnaddedRead
     *
     */
    OmpNestMutex readStatusLock;
    read_t numReads;

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
    bool getRead(read_t &read);

    /**
     * @brief Sets a read to unadded
     *
     * @param read
     * @return true
     * @return false
     */
    void putReadBack(read_t read);

    /**
     * @brief Sets an unadded read to added
     *
     * @param read
     * @return true
     * @return false
     */
    bool addRead(read_t read);

    /**
     * @brief Add reads overlapping with the mainPath of cG at position curPos
     *
     * @param cG
     * @param curPos
     * @param len The length of the mainPath used to get filtered reads
     */
    void addRelatedReads(ConsensusGraph *cG, ssize_t curPos, size_t len);

    ConsensusGraph *createGraph();
};

#endif // Z_CONSENSUS_CUH
