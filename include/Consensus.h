#ifndef Z_CONSENSUS_CUH
#define Z_CONSENSUS_CUH

#include "Contig.h"
#include "Edits.h"
#include "ReadData.h"
#include <map>
#include <set>
#include <string>
#include <vector>

typedef uint32_t POS_T;

class Node;

class Edge;

class Node {
public:
    /***
     * The base (one of ATCG) that the node corresponds to.
     */
    const char base;

    /***
     * Whether this node is on the main path
     */
    bool onMainPath = false;

    /***
     * Creates a Node with given base
     * @param base
     */
    Node(const char base);

    /***
     * Gets the edge going to Node n. Returns null if the edge is not found.
     * @param n The node to search for
     * @return The edge to Node n or null if such an edge does not exist.
     */
    Edge *getEdgeTo(Node *n);

    /***
     * Returns an edge to a node that is not on mainPath with given base.
     * Returns NULL if not found.
     * @param base
     * @return
     */
    Edge *getEdgeToSide(char base);

    /***
     * Returns the edge with the most weight going out of this node.
     * @return the edge with the most weight going out of this node. null if no
     * edge goes out of this node.
     */
    Edge *getBestEdge();

    /***
     * Adds an edge starting from this node
     * @param e The edge to add.
     */
    //    void addEdge(Edge *e);

    // The edges with this node as source
    std::map<Node *, Edge *> edgesOut;

    std::set<Edge *> edgesIn;

    size_t cumulativeWeight = 0;
    bool hasReached = false;

private:
};

/***
 * Each directed edge has a source Node and a sink Node. The field count denotes
 * the number of reads going down that edge. An Edge also stores all the reads
 * going down that path in a vector.
 */
class Edge {
public:
    Node *source;
    Node *sink;
    size_t count;
    std::set<size_t> reads;

    /***
     * Initializes the Edge from a source, a sink, and a single read.
     * @param source
     * @param sink
     * @param read
     */
    Edge(Node *source, Node *sink, size_t read);

    /***
     * Adds a read to the edge. Basically just increases count by 1 and add read
     * to the internal vector.
     * @param read
     */
    void addRead(size_t read);
};

/***
 * A Path is a sequence of Node -> Edge -> Node .... -> Edge -> Node.
 * We only need to store the edges
 */
class Path {
public:
    std::vector<Edge *> edges;
    std::string path;

    /***
     * Returns the average weight of the edges on path
     * @return the average weight of the edges on path
     */
    double getAverageWeight();

    /***
     * Clears the main path. This means setting onMainPath to false for
     * every node on this path.
     */
    void clear();
};

/***
 * We store the current reads along with their relations as a directly acyclic
 * graph. The nodes represent different bases, and each each edge means there is
 * at least one read where the base in the next node follows the current node.
 * The nodes contain the base, and the edges going in and the edges going out.
 * The edges contain a count of reads going down that path.
 */
class ConsensusGraph {
public:
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

    class Read {
    public:
        /***
         * Relative position of read in contig
         */
        long pos;
        /***
         * Starting node of read
         */
        Node *start;
        /***
         * Length of read
         */
        size_t len;

        /***
         *
         * @param pos Relative position of read in contig
         * @param start Starting node of read
         * @param len Length of read
         */
        Read(long pos, Node *start, size_t len);
    };

    ConsensusGraph(StringAligner *aligner);

    /***
     * Initializes the graph from a seeding read
     * @param seed
     * @param readId id of the seeding read
     * @param pos Relative position in contig
     */
    void initialize(const std::string &seed, size_t readId, long pos);

    /***
     * Adds a read to the consensus graph. Does not update mainPath.
     * @param s String of read to add
     * @param readId id of the read to add (id right now just means the position
     * in all the reads read)
     * @param pos Relative position of read in contig
     */
    bool addRead(const std::string &s, long pos, std::vector<Edit> &editScript,
                 ssize_t &beginOffset, ssize_t &endOffset);

    void addReads(const std::set<std::pair<long, read_t>> &reads,
                  std::vector<std::unique_ptr<std::string>> &readData);

    /***
     * Clears the old mainPath, calculates the new mainPath and adds it.
     * Uses dynamic programming.
     * @return the new mainPath
     */
    Path &calculateMainPath();

    /***
     * Clears the old mainPath, calculates the new mainPath and adds it.
     * Uses greedy algorithm.
     * @return the new mainPath
     */
    Path &calculateMainPathGreedy();

    Path &getMainPath();

    /**
     * @brief
     *
     * @param filename
     */
    void writeMainPath(const std::string &filename);

    /***
     * @param editScript outputs the editScript (consisting of insertions,
     * deletions, and unchanged)
     * @param pos outputs the relative position in mainPath
     */
    size_t read2EditScript(Read &r, size_t id, std::vector<Edit> &editScript,
                           size_t &pos);

    /**
     * @brief Write the edit strings of the reads in a single file f
     *
     * @param f
     */
    void writeReads(std::ofstream &f);

    /**
     * @brief Write the raw data into tempDir with file names .pos, .genome,
     * .type, .base, .id, .unalignedReads, .unalignedIds; compresses them into
     * tempCompressedDir;
     * The .id file is delta encoded
     *
     * @param filename
     */
    void writeReads(const std::string &filename);

    /***
     * Prints the info of the ConsensusGraph. For debugging purposes
     */
    void printStatus();

    ~ConsensusGraph();

private:
    Node *startingNode;
    Node *leftMostChangedNode;
    size_t leftMostChangedNodeOffset = 0;
    size_t numNodes = 0;
    size_t numEdges = 0;
    // Maps ID of read to (relative position of read in contig, beginning node
    // of the read)
    std::map<size_t, Read> readsInGraph;
    Path mainPath;
    // Starting and ending positions of mainPath in contig
    ssize_t startPos, endPos;
    StringAligner *aligner;
    // Maps ID of reads that cannot be successfully aligned to actual strings
    std::map<size_t, std::string> unalignedReads;

    Node *createNode(char base);

    Edge *createEdge(Node *source, Node *sink, size_t read);

    void removeReadsFromEdge(Edge *e, std::set<size_t> const &reads);

    /***
     * Removes an edge from the graph
     */
    void removeEdge(Edge *e);

    void removeNode(Node *n);

    void removeBelow(Node *n);

    void removeAbove(Node *n);

    void removeConnectedNodes(Node *n);

    void updateGraph(const std::string &s, std::vector<Edit> &editScript,
                     ssize_t beginOffset, ssize_t endOffset, size_t readId,
                     long pos);

    void clearHasReached(Node *n);

    /***
     * Remove possible cycles from mainPath. Should be called internally
     * every time the mainPath is updated. Cycles are avoided by ensuring
     * that every side path only has one edge in from the main path.
     */
    void removeCycles();

    void walkAndPrune(Edge *e);

    /***
     * Move everything below e to newPre->e (everything means all
     * the reads in reads2Split until the mainPath or a leaf node is
     * reached)
     */
    void splitPath(Node *newPre, Edge *e, std::set<size_t> const &reads2Split);

    /***
     * Write the edits trings of the reads into a single file
     * @param f
     * @param r
     * @param id
     * @return Edit distance
     */
    size_t writeRead(std::ofstream &f, Read &r, size_t id);

    /**
     * @brief
     *
     * @param posFile Here we store the number of unchanged bases before the
     * current edit. For example, the (n+1)th position stores the number of
     * uncahgned bases before the nth edit. The 0th position stores the initial
     * position of the read. If there are num edits, then we store (num + 2)
     * positions, where the last position is the number of unchanged bases after
     * the last edit. The data between different reads are separated by
     * END_SYMBOL.
     * @param editTypeFile
     * i[nsertion], d[eletion]. s[ubstitution]. Different reads are separated by
     * '\n'
     * @param editBaseFile
     * ATCG.  Different reads are separated by '\n'
     * @param r
     * @param id
     * @return size_t
     */
    size_t writeRead(std::ofstream &posFile, std::ofstream &editTypeFile,
                     std::ofstream &editBaseFile, Read &r, size_t id);

    /**
     * @brief Stores the reads in filename.unalignedReads and the ids in
     * filename.unalignedIds.
     *
     * filename.unalignedIds is delta coded.
     * filename.unalignedReads is just one line per read
     *
     * @param filename
     */
    void writeUnalignedReads(const std::string &filename);

    void clearMainPath();

    // void writeGraph(std::ofstream &f);
};

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

    std::vector<ConsensusGraph *> graphs;

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

    ConsensusGraph *createGraph();
};

#endif // Z_CONSENSUS_CUH