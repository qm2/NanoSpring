#ifndef CCCEF0A4_30C3_49E3_999D_C948EA366BEF
#define CCCEF0A4_30C3_49E3_999D_C948EA366BEF

#include "Contig.h"
#include "Edits.h"
#include "ReadData.h"
#include "StringAligner.h"
#include <deque>
#include <functional>
#include <map>
#include <set>
#include <vector>

class Node;

class Edge;

class Node {
public:
    /**
     * The base (one of ATCG) that the node corresponds to.
     */
    const char base;

    /**
     * Whether this node is on the main path
     */
    bool onMainPath = false;

    /**
     * Creates a Node with given base
     * @param base
     */
    Node(const char base);

    /**
     * Gets the edge going to Node n. Returns null if the edge is not found.
     * @param n The node to search for
     * @return The edge to Node n or null if such an edge does not exist.
     */
    Edge *getEdgeTo(Node *n);

    /**
     * Returns an edge to a node that is not on mainPath with given base.
     * Returns NULL if not found.
     * @param base
     * @return
     */
    Edge *getEdgeToSide(char base);

    /**
     * Returns the edge with the most weight going out of this node.
     * @return the edge with the most weight going out of this node. null if no
     * edge goes out of this node.
     */
    Edge *getBestEdgeOut();

    /**
     * Returns the edge with the most weight going in to this node.
     * @return the edge with the most weight going in to this node. null if no
     * edge goes in to this node.
     */
    Edge *getBestEdgeIn();

    /**
     * @brief Finds an edge going out containing the given read. Returns nullptr
     * is no such read is found.
     *
     * Assumes that the reads are already sorted.
     *
     * @param read
     * @return Edge*
     */
    Edge *getEdgeInRead(read_t read) const;

    /**
     * @brief Finds the node next in read. Returns nullptr if not found.
     *
     * @param read
     * @return Node*
     */
    Node *getNextNodeInRead(read_t read) const;

    /**
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

/**
 * Each directed edge has a source Node and a sink Node. The field count denotes
 * the number of reads going down that edge. An Edge also stores all the reads
 * going down that path in a vector.
 */
class Edge {
public:
    Node *source;
    Node *sink;
    read_t count;
    // always maintained in sorted order
    std::vector<read_t> reads;

    /**
     * Initializes the Edge from a source, a sink, and a single read.
     * @param source
     * @param sink
     * @param read
     */
    Edge(Node *source, Node *sink, read_t read);

    /**
     * Adds a read to the edge. Basically just increases count by 1 and add read
     * to the internal vector (while maintaining sorting).
     * @param read
     */
    void addRead(read_t read);
};

/**
 * A Path is a sequence of Node -> Edge -> Node .... -> Edge -> Node.
 * We only need to store the edges
 */
class Path {
public:
    std::deque<Edge *> edges;
    std::deque<char> path;

    /**
     * Returns the average weight of the edges on path
     * @return the average weight of the edges on path
     */
    double getAverageWeight();

    /**
     * Clears the main path. This means setting onMainPath to false for
     * every node on this path.
     */
    void clear();
};

/**
 * We store the current reads along with their relations as a directly acyclic
 * graph. The nodes represent different bases, and each each edge means there is
 * at least one read where the base in the next node follows the current node.
 * The nodes contain the base, and the edges going in and the edges going out.
 * The edges contain a count of reads going down that path.
 */
class ConsensusGraph {
public:
    typedef std::deque<char>::const_iterator RAItA;
    typedef const char *RAItB;
    typedef StringAligner<RAItA, RAItB> StringAligner_t;

    /** Starting and ending positions of mainPath in contig **/
    ssize_t startPos, endPos;

    /** Directory for storing the temp files (.genome, .pos, .type) **/
    std::string tempDir = "tempRaw/";

    Path mainPath;

    class Read {
    public:
        /**
         * Relative position of read in contig
         */
        long pos;
        /**
         * Starting node of read
         */
        Node *start;
        /**
         * Length of read
         */
        size_t len;

        /**
         *
         * @param pos Relative position of read in contig
         * @param start Starting node of read
         * @param len Length of read
         */
        Read(long pos, Node *start, size_t len);
    };

    /** Maps ID of read to (relative position of read in contig, beginning node
     of the read) **/
    std::map<read_t, Read> readsInGraph;

    ConsensusGraph(StringAligner_t *aligner);

    /**
     * Initializes the graph from a seeding read
     * @param seed
     * @param readId id of the seeding read
     * @param pos Relative position in contig
     */
    void initialize(const std::string &seed, read_t readId, long pos);

    /**
     * Adds a read to the consensus graph. Does not update mainPath or list of
     * reads.
     * @param s String of read to add
     * @param pos Relative position of read in contig
     */
    __attribute__((warn_unused_result)) bool
    addRead(const std::string &s, long pos, std::vector<Edit> &editScript,
            ssize_t &beginOffset, ssize_t &endOffset);

    void addReads(const std::set<std::pair<long, read_t>> &reads,
                  std::vector<std::unique_ptr<std::string>> &readData);

    /**
     * @brief Updates the graph with the new read s, and the alignment results
     * editScript, beginOffset, and endOffset
     *
     * @param s
     * @param editScript
     * @param beginOffset
     * @param endOffset
     * @param readId
     * @param pos The relative position of the read in contig
     */
    void updateGraph(const std::string &s, std::vector<Edit> &editScript,
                     ssize_t beginOffset, ssize_t endOffset, read_t readId,
                     long pos);

    /**
     * Clears the old mainPath, calculates the new mainPath and adds it.
     * Uses dynamic programming.
     * @return the new mainPath
     */
    Path &calculateMainPath();

    /**
     * Clears the old mainPath, calculates the new mainPath and adds it.
     * Uses greedy algorithm.
     * Uses the cumulative weight field to store the position of the nodes on
     * mainPath (starts at 0)
     * @return the new mainPath
     */
    Path &calculateMainPathGreedy();

    /**
     * @brief Prints the mainPath to file tempDir/filename.genome
     *
     * @param filename
     */
    void writeMainPath(const std::string &filename);

    /**
     * @param editScript outputs the editScript (consisting of insertions,
     * deletions, and unchanged)
     * For this to work properly, the cumulativeWeight field of the nodes on the
     * main path must store their positions on main path (starting with 0)
     * @param pos outputs the relative position in mainPath
     */
    size_t read2EditScript(Read &r, read_t id, std::vector<Edit> &editScript,
                           size_t &pos);

    /**
     * @brief Write the edit strings of the reads in a single file f
     *
     * @param f
     */
    void writeReads(std::ofstream &f);

    /**
     * @brief Write the raw data into tempDir with file names .pos, .genome,
     * .type, .base, .id, .unalignedReads, .unalignedIds;
     * The .id file is delta encoded
     *
     * @param filename
     */
    void writeReads(const std::string &filename);

    /**
     * Prints the info of the ConsensusGraph. For debugging purposes
     */
    void printStatus();

    read_t getNumReads();

    /**
     * @brief Gets the read identified by read and stores the bases into
     * inserter.
     *
     * @tparam Inserter And insert iterator that dereferences to char
     * @param read
     * @param inserter
     * @return true If read was found.
     * @return false If read was not found.
     */
    template <typename Inserter> bool getRead(read_t read, Inserter inserter);

    ~ConsensusGraph();

private:
    Node *startingNode;
    // These are stored such that the path only needs to be updated locally
    Node *rightMostUnchangedNode;
    size_t rightMostUnchangedNodeOffset = 0;
    Node *leftMostUnchangedNode;
    size_t leftMostUnchangedNodeOffset = 0;
    size_t numNodes = 0;
    size_t numEdges = 0;

    StringAligner_t *aligner;
    // Maps ID of reads that cannot be successfully aligned to actual strings
    std::map<read_t, std::string> unalignedReads;

    Node *createNode(char base);

    Edge *createEdge(Node *source, Node *sink, read_t read);

    // reads must be a sorted vector
    void removeReadsFromEdge(Edge *e, std::vector<read_t> const &reads);

    /**
     * Removes an edge from the graph
     */
    void removeEdge(Edge *e);

    void removeNode(Node *n);

    void removeBelow(Node *n);

    void removeAbove(Node *n);

    void removeConnectedNodes(Node *n);

    /**
     * @brief Traverses the graph starting at Node n and calls the function f on
     * each Node
     *
     * Assumes that the graph is connected; assumes that originally some portion
     * of the graph connected to Node n has hasReached == status, and the rest
     * has hasReached == !status. Will call f on all the Nodes connected to n
     * with hasReached == status, and will set their hasReached to !status when
     * the call returns.
     *
     * @tparam Functor
     * @param n
     * @param f
     * @param status
     */
    template <typename Functor>
    void traverseAndCall(Node *n, bool status, Functor f);

    /**
     * @brief
     * Remove possible cycles from mainPath. Should be called internally
     * every time the mainPath is updated. Cycles are avoided by ensuring
     * that every side path only has one edge in from the main path.
     */
    void removeCycles();

    /**
     * @brief Does a depth first search of all edges on side paths that are
     * reachable from e, and guarantees that the sink of any such edge only has
     * one edge going in.
     * @details
     * Post condition: for all Edge edge such that
     *
     * - edge is on side paths
     * - edge is reachable from e (including the case where edge == e)
     *
     * the number of edges going into edge will be <= 1.
     *
     * All edges and nodes NOT reachable from e via sidepaths are kept constant.
     *
     * In particular, e->source might receive new edges going out, and the old
     * might be deleted. But the newly created edges are guaranteed to need no
     * further pruning.
     * @param e
     */
    void walkAndPrune(Edge *e);

    /**
     * @brief
     * Move everything below e to newPre->newNode.
     * @details
     * "Everything" means all the reads in reads2Split until the mainPath or a
     * leaf node is reached. reads2Split must be sorted vector.
     *
     * - Nodes and Edges that will be deleted: \n
     *  e will be erased if it only has reads in reads2Split; e->sink will be
     * erased if it becomes disconnected from the graph. The same applies to
     * nodes and edges reachable from e via sidepaths.
     * - Nodes and Edges that will be created: \n
     *  A duplicate of "everything" below e in reads2Split will be created under
     * newPre. In particular, newPre will have a new edge going out. The newly
     * created nodes are guaranteed to have at most one Edge going in.
     */
    void splitPath(Node *newPre, Edge *e,
                   std::vector<read_t> const &reads2Split);

    /**
     * Write the edits trings of the reads into a single file
     * @param f
     * @param r
     * @param id
     * @return Edit distance
     */
    size_t writeRead(std::ofstream &f, Read &r, read_t id);

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
                     std::ofstream &editBaseFile, Read &r, read_t id);

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

/******************************************************************************/
/* Implementations of templates */
template <typename Inserter>
bool ConsensusGraph::getRead(read_t read, Inserter inserter) {
    auto readIt = readsInGraph.find(read);
    if (readIt == readsInGraph.end())
        return false;
    Node *curNode = readIt->second.start;
    while (curNode) {
        *(inserter++) = curNode->base;
        curNode = curNode->getNextNodeInRead(read);
    }
    return true;
}

#endif /* CCCEF0A4_30C3_49E3_999D_C948EA366BEF */
