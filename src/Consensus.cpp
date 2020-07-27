#include <iostream>
#include "../include/Consensus.h"

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

void Path::clear() {
    size_t l = edges.size();
    if (l > 0) {
        Node *currentNode = edges[0]->source;
        currentNode->onMainPath = false;
    }
    for (size_t i = 0; i < l; i++) {
        Node *currentNode = edges[i]->sink;
        currentNode->onMainPath = false;
    }
    edges.clear();
    path.clear();
}

void ConsensusGraph::initialize(const std::string &seed, size_t readId, long pos) {
    size_t len = seed.length();
    Node *currentNode = createNode(seed[0]);
    // We create a read that points to this node
    reads.insert(std::make_pair(readId,
                                ConsensusGraph::Read(pos, currentNode, seed.length())));
    startingNode = currentNode;
    for (size_t i = 1; i < len; ++i) {
        Node *nextNode = createNode(seed[i]);
        Edge *e = createEdge(currentNode, nextNode, readId);
        currentNode->addEdge(e);
        currentNode = nextNode;
    }
}

void ConsensusGraph::addRead(const std::string &s, size_t readId, long pos) {
    std::cout << "mainPath startPos " << startPos
              << " mainPath endPos " << endPos
              << " current read startPos " << pos
              << std::endl;
    std::string &originalString = mainPath.path;
    std::vector<Edit> editScript;
    aligner->align(originalString, s, editScript);
}

Path &ConsensusGraph::calculateMainPath() {
    mainPath.clear();
    std::vector<Edge *> &edgesInPath = mainPath.edges;
    std::string &stringPath = mainPath.path;
    Node *currentNode = startingNode;
    currentNode->onMainPath = true;
    // Updating string
    stringPath.push_back(currentNode->base);
    Edge *edgeToAdd;
    while (edgeToAdd = currentNode->getBestEdge()) {
        edgesInPath.push_back(edgeToAdd);
        currentNode = edgeToAdd->sink;
        currentNode->onMainPath = true;
        stringPath.push_back(currentNode->base);
    }
    size_t startingReadId = edgesInPath.front()->reads[0];
    startPos = reads.at(startingReadId).pos;
    size_t endingReadId = edgesInPath.back()->reads[0];
    Read &endingRead = reads.at(endingReadId);
    endPos = endingRead.pos + endingRead.len;
    return mainPath;
}

ConsensusGraph::~ConsensusGraph() {
    for (auto n : nodes)
        delete n;
    for (auto e : edges)
        delete e;
    delete aligner;
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
              << edges.size() << " edges."
              << std::endl;
    std::cout << "mainPath len " << mainPath.edges.size() + 1
              << " avg weight " << mainPath.getAverageWeight()
              << " starts at " << startPos
              << " ends at " << endPos
              << " len in Contig " << endPos - startPos
              << std::endl;
}

ConsensusGraph::ConsensusGraph(StringAligner *aligner) : aligner(aligner) {}

ConsensusGraph::Read::Read(long pos, Node *start, size_t len) : pos(pos), start(start), len(len) {}
