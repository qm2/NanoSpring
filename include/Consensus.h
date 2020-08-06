/***
 * We store the current reads along with their relations as a directly acyclic
 * graph. The nodes represent different bases, and each each edge means there is
 * at least one read where the base in the next node follows the current node.
 * The nodes contain the base, and the edges going in and the edges going out.
 * The edges contain
 */

#ifndef Z_CONSENSUS_CUH
#define Z_CONSENSUS_CUH

#endif // Z_CONSENSUS_CUH

#include "Contig.h"
#include "Edits.h"
#include <map>
#include <set>
#include <string>
#include <vector>

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

class ConsensusGraph {
public:
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
     * @param readId id of the read to add
     * @param pos Relative position of read in contig
     */
    bool addRead(const std::string &s, size_t readId, long pos,
                 std::vector<Edit> &editScript, ssize_t &beginOffset,
                 ssize_t &endOffset);

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

    void writeMainPath(const std::string &filename);

    void writeReads(std::ofstream &f);

    void writeReads(const std::string &filename);

    /***
     * Prints the info of the ConsensusGraph. For debugging purposes
     */
    void printStatus();

    ~ConsensusGraph();

private:
    Node *startingNode;
    size_t numNodes = 0;
    size_t numEdges = 0;
    // Maps ID of read to (relative position of read in contig, beginning node
    // of the read)
    std::map<size_t, Read> readsInGraph;
    Path mainPath;
    // Starting and ending positions of mainPath in contig
    long startPos, endPos;
    StringAligner *aligner;

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
     * Move everything below oldPre->e to newPre->e (everything means all
     * the reads in reads2Split until the mainPath or a leaf node is
     * reached)
     */
    void splitPath(Node *oldPre, Node *newPre, Edge *e,
                   std::set<size_t> const &reads2Split);

    /***
     *
     * @param f
     * @param r
     * @param id
     * @return Edit distance
     */
    size_t writeRead(std::ofstream &f, Read &r, size_t id);

    size_t writeRead(std::ofstream &posFile, std::ofstream &editTypeFile,
                     std::ofstream &editBaseFile, Read &r, size_t id);

    // void writeGraph(std::ofstream &f);
};