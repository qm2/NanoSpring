/***
 * We store the current reads along with their relations as a directly acyclic graph.
 * The nodes represent different bases, and each each edge means there is at least
 * one read where the base in the next node follows the current node.
 * The nodes contain the base, and the edges going in and the edges going out.
 * The edges contain
 */

#ifndef Z_CONSENSUS_CUH
#define Z_CONSENSUS_CUH

#endif //Z_CONSENSUS_CUH

#include <vector>
#include <map>
#include <string>

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
     * Returns the edge with the most weight going out of this node.
     * @return the edge with the most weight going out of this node. null if no
     * edge goes out of this node.
     */
    Edge *getBestEdge();

    /***
     * Adds an edge starting from this node
     * @param e The edge to add.
     */
    void addEdge(Edge *e);

private:
    // The edges with this node as source
    std::map<Node *, Edge *> edgesFrom;
};

/***
 * Each directed edge has a source Node and a sink Node. The field count denotes
 * the number of reads going down that edge. An Edge also stores all the reads going down
 * that path in a vector.
 */
class Edge {
public:
    Node *source;
    Node *sink;
    unsigned int count;
    std::vector<unsigned int> reads;

    /***
     * Initializes the Edge from a source, a sink, and a single read.
     * @param source
     * @param sink
     * @param read
     */
    Edge(Node *source, Node *sink, unsigned int read);

    /***
     * Adds a read to the edge. Basically just increases count by 1 and add read
     * to the internal vector.
     * @param read
     */
    void addRead(unsigned int read);
};

/***
 * A Path is a sequence of Node -> Edge -> Node .... -> Edge -> Node.
 * We only need to store the edges
 */
class Path {
public:
    std::vector<Edge *> edges;

    /***
     * Returns the average weight of the edges on path
     * @return the average weight of the edges on path
     */
    double getAverageWeight();
};

class ConsensusGraph {
public:
    /***
     * Initializes the graph from a seeding read
     * @param seed
     */
    void initialize(const std::string &seed, size_t readId);

    Path &calculateMainPath();

    Path &getMainPath();

    /***
     * Prints the info of the ConsensusGraph. For debugging purposes
     */
    void printStatus();

    ~ConsensusGraph();

private:
    Node *startingNode;
    std::vector<Node *> nodes;
    std::vector<Edge *> edges;
    // Maps ID of read to beginning node of the read
    std::map<size_t, Node *> reads;
    Path mainPath;

    Node *createNode(char base);

    Edge *createEdge(Node *source, Node *sink, unsigned int read);
};