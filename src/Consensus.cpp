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
    edgesFrom[e->sink] = e;
}

Edge *Node::getEdgeTo(Node *n) {
    try {
        return edgesFrom.at(n);
    } catch (std::out_of_range &e) {
        return NULL;
    }
}

Edge *Node::getEdgeToSide(char base) {
    auto end = edgesFrom.end();
    for (auto it = edgesFrom.begin(); it != end; it++) {
        Edge *e = it->second;
        Node *n = e->sink;
        if (!n->onMainPath && n->base == base)
            return e;
    }
    return NULL;
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
        currentNode = nextNode;
    }
}

void ConsensusGraph::addRead(const std::string &s, size_t readId, long pos) {
    std::string &originalString = mainPath.path;
    std::vector<Edit> editScript;
    const ssize_t offsetGuess = originalString.length() - endPos + pos;
    size_t editDis;
    ssize_t beginOffset;
    ssize_t endOffset;

    bool success = aligner->align(originalString, s, offsetGuess,
                                  beginOffset, endOffset, editScript, editDis);
//    std::cout << "success ? " << success << std::endl;
    if (!success) {
        std::cout << "Failed to add" << std::endl;
        std::cout << "mainPath startPos " << startPos
                  << " mainPath endPos " << endPos
                  << " current read startPos " << pos
                  << " OffsetGuess " << offsetGuess
                  << std::endl;
        return;
    }

    auto edgeInPath = mainPath.edges.begin();
    const auto edgeInPathEnd = mainPath.edges.end();
    Node *nodeInPath = (*edgeInPath)->source;
    Node *currentNode = NULL;
    Node *initialNode = NULL;
    auto advanceNodeInPath = [&]() {
        if (edgeInPath == edgeInPathEnd)
            return;

        nodeInPath = (*edgeInPath)->sink;

//        std::cout << "NumOfEdgeInPath" << edgeInPath - mainPath.edges.begin()
//                  << " OutOf " << mainPath.edges.size()
//                  << std::endl;
        edgeInPath++;
    };

    // First we deal with beginOffset
    auto initialAdvance = [&] {
//        std::cout << "initialAdvance" << std::endl;
        if (beginOffset >= 0) {
            // We need to advance in the mainPath
            for (size_t i = 0; i < beginOffset; i++) {
                advanceNodeInPath();
            }
        } else {
            // We need to insert the initial parts of the read
            // This number must be positive
            size_t numOfNodes2Insert = -beginOffset;
            size_t i = 0;
            // We create an initial node
            currentNode = createNode(s[i++]);
            initialNode = currentNode;
            for (; i < numOfNodes2Insert; ++i) {
                Node *nextNode = createNode(s[i]);
                Edge *edge = createEdge(currentNode, nextNode, readId);
                currentNode = nextNode;
            }
        }
    };

    initialAdvance();

    auto insertNode = [&](char base) {
        if (!currentNode) {
            // If there is no currentNode, we create one
            currentNode = createNode(base);
            initialNode = currentNode;
        } else {
            Edge *edge = currentNode->getEdgeToSide(base);
            if (edge) {
                edge->addRead(readId);
            } else {
                Node *n = createNode(base);
                edge = createEdge(currentNode, n, readId);
            }
            currentNode = edge->sink;
        }
    };

//    std::cout << "Edits" << std::endl;

    // Now we deal with the edits one by one
    for (const Edit &e: editScript) {
//        std::cout << e;
        if (e.editType == SAME) {
            size_t num = e.editInfo.num;

            // We deal with the first node
            if (!currentNode) {
                // If there is no currentNode, set it and initialNode
                initialNode = nodeInPath;
                currentNode = nodeInPath;
            } else {
                // Otherwise create an edge from the currentNode to nodeInPath
                Edge *edge = currentNode->getEdgeTo(nodeInPath);
                if (edge) {
                    // If there is already an edge to it
                    edge->addRead(readId);
                } else {
                    // Otherwise create an edge
                    edge = createEdge(currentNode, nodeInPath, readId);
                }
                currentNode = nodeInPath;
            }
            advanceNodeInPath();

            // We deal with the rest of the nodes
            for (size_t i = 1; i < num; ++i) {
                Edge *temp = currentNode->getEdgeTo(nodeInPath);
                currentNode->getEdgeTo(nodeInPath)->addRead(readId);
                currentNode = nodeInPath;
                advanceNodeInPath();
            }
        } else if (e.editType == DELETE) {
            // Deletion just means advancing the node in path
            advanceNodeInPath();
        } else if (e.editType == INSERT) {
            char base = e.editInfo.ins;
            insertNode(base);
        }
    }

//    std::cout << "EndOffset" << std::endl;

    // Finally we deal with endOffset
    if (endOffset > 0) {
//        std::cout << endOffset << std::endl;
        auto end = s.end();
        for (auto it = end - endOffset; it < end; it++) {
//            std::cout << *it;
            insertNode(*it);
        }
    }

    // Don't forget to add the read!
    reads.insert(std::make_pair(readId, Read(pos, initialNode, s.length())));
}

Path &ConsensusGraph::calculateMainPath() {
//    std::cout << "Calculating Main Path" << std::endl;
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
//    std::cout << "Finished Calculating Main Path" << std::endl;
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
    source->addEdge(e);
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
