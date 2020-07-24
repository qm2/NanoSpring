#include <iostream>
#include "../include/Consensus.cuh"

Edge::Edge(Node *source, Node *sink, unsigned int read) : source(source), sink(sink) {
    count = 1;
    reads.push_back(read);
}

void Edge::addRead(unsigned int read) {
    count++;
    reads.push_back(read);
}

Node::Node(const char base) : base(base) {}

void Node::addEdge(Edge *e) {
    edgesFrom[e->source] = e;
}

Edge *Node::getBestEdge() {
    Edge *bestEdge = NULL;
    size_t bestCount = 0;
    auto end = edgesFrom.end();
    for (auto it = edgesFrom.begin(); it != end; ++it) {
        Edge *e = it->second;
        if (e->count > bestCount) {
            bestCount = e->count;
            bestEdge = e;
        }
    }
    return bestEdge;
}

double Path::getAverageWeight() {
    size_t totalWeight = 0;
    for (Edge *e : edges)
        totalWeight += e->count;
    return ((double) totalWeight) / edges.size();
}

void ConsensusGraph::initialize(const std::string &seed, size_t readId) {
    size_t len = seed.length();
    Node *currentNode = createNode(seed[0]);
    // We create a read that points to this node
    reads[readId] = currentNode;
    startingNode = currentNode;
    for (size_t i = 1; i < len; ++i) {
        Node *nextNode = createNode(seed[i]);
        Edge *e = createEdge(currentNode, nextNode, readId);
        currentNode->addEdge(e);
        currentNode = nextNode;
    }
}

Path &ConsensusGraph::calculateMainPath() {
    std::vector<Edge *> &edgesInPath = mainPath.edges;
    edgesInPath.clear();
    Node *currentNode = startingNode;
    Edge *edgeToAdd;
    while (edgeToAdd = currentNode->getBestEdge()) {
        edgesInPath.push_back(edgeToAdd);
        currentNode = edgeToAdd->sink;
    }
    return mainPath;
}

ConsensusGraph::~ConsensusGraph() {
    for (auto n : nodes)
        delete n;
    for (auto e : edges)
        delete e;
}

Node *ConsensusGraph::createNode(char base) {
    Node *n = new Node(base);
    nodes.push_back(n);
    return n;
}

Edge *ConsensusGraph::createEdge(Node *source, Node *sink, unsigned int read) {
    Edge *e = new Edge(source, sink, read);
    edges.push_back(e);
    return e;
}

void ConsensusGraph::printStatus() {
    std::cout << reads.size() << " reads, "
              << nodes.size() << " nodes, and "
              << edges.size() << " edges." << std::endl;
    std::cout << "mainPath len " << mainPath.edges.size() + 1
              << " average weight " << mainPath.getAverageWeight() << std::endl;
}