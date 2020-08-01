#include <iostream>
#include <algorithm>
#include <iterator>
#include <deque>
#include <fstream>
#include "../include/Consensus.h"

Edge::Edge(Node *source, Node *sink, size_t read) : source(source), sink(sink) {
    count = 1;
    reads.insert(read);
}

void Edge::addRead(size_t read) {
    count++;
    reads.insert(read);
}

Node::Node(const char base) : base(base) {}

//void Node::addEdge(Edge *e) {
//    edgesOut[e->sink] = e;
//}

Edge *Node::getEdgeTo(Node *n) {
    const auto &it = edgesOut.find(n);
    if (it == edgesOut.end()) {
        return NULL;
    } else {
        return it->second;
    }
}

Edge *Node::getEdgeToSide(char base) {
    const auto &end = edgesOut.end();
    for (auto it = edgesOut.begin(); it != end; it++) {
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
    const auto &end = edgesOut.end();
    for (auto it = edgesOut.begin(); it != end; ++it) {
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
        std::cout << "Failed to add" << "\n";
        std::cout << "mainPath startPos " << startPos
                  << " mainPath endPos " << endPos
                  << " current read startPos " << pos
                  << " OffsetGuess " << offsetGuess
                  << std::endl;
        return;
    }

    auto edgeInPath = mainPath.edges.begin();
    const auto &edgeInPathEnd = mainPath.edges.end();
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
    mainPath.clear();

    std::deque<Node *> unfinishedNodes;
    unfinishedNodes.push_back(startingNode);
    // First we set all the hasReached fields to false
    while (!unfinishedNodes.empty()) {
        Node *currentNode = unfinishedNodes.front();
        if (!currentNode->hasReached) {
            unfinishedNodes.pop_front();
            continue;
        }

        currentNode->hasReached = false;
        unfinishedNodes.pop_front();

        for (auto edgeIt : currentNode->edgesOut) {
            Node *n = edgeIt.second->sink;
            if (n->hasReached) {
                unfinishedNodes.push_front(n);
            }
        }

        for (auto edgeIt : currentNode->edgesIn) {
            if (edgeIt->source->hasReached)
                unfinishedNodes.push_back(edgeIt->source);
        }
    }

    size_t globalMaxWeight = 0;
    Node *globalMaxWeightNode = NULL;
    unfinishedNodes.push_back(startingNode);
    while (!unfinishedNodes.empty()) {
        Node *currentNode = unfinishedNodes.front();
        if (currentNode->hasReached) {
            unfinishedNodes.pop_front();
            continue;
        }
        size_t maxWeight = 0;
        bool hasPreviousWeights = true;
        for (auto edgeIt : currentNode->edgesOut) {
            Node *n = edgeIt.second->sink;
            if (!edgeIt.second->sink->hasReached) {
                unfinishedNodes.push_front(n);
                hasPreviousWeights = false;
                break;
            }
            size_t prevWeight = edgeIt.second->sink->cumulativeWeight;
            size_t curWeight = prevWeight + edgeIt.second->count;
            maxWeight = std::max(maxWeight, curWeight);
        }
        if (!hasPreviousWeights)
            continue;
        globalMaxWeightNode = maxWeight >= globalMaxWeight ? currentNode : globalMaxWeightNode;
        globalMaxWeight = std::max(maxWeight, globalMaxWeight);
        unfinishedNodes.pop_front();
        for (auto edgeIt : currentNode->edgesIn) {
            unfinishedNodes.push_back(edgeIt->source);
        }
        currentNode->cumulativeWeight = maxWeight;
        currentNode->hasReached = true;
    }

    startingNode = globalMaxWeightNode;

    std::vector<Edge *> &edgesInPath = mainPath.edges;
    std::string &stringPath = mainPath.path;

    stringPath.push_back(globalMaxWeightNode->base);
    globalMaxWeightNode->onMainPath = true;
//    std::cout << "max weight" << globalMaxWeight << std::endl;
//    std::cout << "here" << std::endl;
    while (globalMaxWeight > 0) {
//        std::cout << "max weight" << globalMaxWeight << " "
//                  << edgesInPath.size() << std::endl;
        for (auto e : globalMaxWeightNode->edgesOut) {
            if (e.second->count + e.first->cumulativeWeight == globalMaxWeight) {
                edgesInPath.push_back(e.second);
                globalMaxWeight -= e.second->count;
                globalMaxWeightNode = e.first;
                stringPath.push_back(globalMaxWeightNode->base);
                e.first->onMainPath = true;
                break;
            };
        }
    }

    size_t startingReadId = *edgesInPath.front()->reads.begin();
    startPos = reads.at(startingReadId).pos;
    size_t endingReadId = *edgesInPath.back()->reads.begin();
    Read &endingRead = reads.at(endingReadId);
    endPos = endingRead.pos + endingRead.len;
    printStatus();
    removeCycles();
    return mainPath;
}

Path &ConsensusGraph::calculateMainPathGreedy() {
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
    size_t startingReadId = *edgesInPath.front()->reads.begin();
    startPos = reads.at(startingReadId).pos;
    size_t endingReadId = *edgesInPath.back()->reads.begin();
    Read &endingRead = reads.at(endingReadId);
    endPos = endingRead.pos + endingRead.len;
    printStatus();
    removeCycles();
    return mainPath;
}

void ConsensusGraph::removeCycles() {
    auto edgeOnPath = mainPath.edges.begin();
    auto edgeOnPathEnd = mainPath.edges.end();
    if (edgeOnPath == edgeOnPathEnd)
        return;
    Node *nodeOnPath = (*edgeOnPath)->source;
    // We first iterate over all nodes on mainPath
    while (true) {
        // Then we iterate over all edges pointing to side nodes that have other edges in
        for (auto &edgeIt : nodeOnPath->edgesOut) {
            Node *sideNode = edgeIt.first;
            if (sideNode->onMainPath || sideNode->edgesIn.size() == 1)
                continue;
            splitPath(nodeOnPath, nodeOnPath, edgeIt.second, edgeIt.second->reads);
        }
        if (edgeOnPath == edgeOnPathEnd)
            break;
        nodeOnPath = (*edgeOnPath)->sink;
        ++edgeOnPath;
    }
}

void ConsensusGraph::splitPath(Node *oldPre, Node *newPre, Edge *e,
                               std::set<size_t> const &reads2Split) {
    // The reads that need to be split going down this path
    std::set<size_t> readsInPath2Split;
    std::set_intersection(reads2Split.begin(), reads2Split.end(),
                          e->reads.begin(), e->reads.end(),
                          std::inserter(readsInPath2Split, readsInPath2Split.begin()));
    if (readsInPath2Split.empty())
        return;

    // Remove the reads from the node edge
    std::set<size_t> updatedReadsInOldEdge;
    std::set_difference(e->reads.begin(), e->reads.end(),
                        readsInPath2Split.begin(), readsInPath2Split.end(),
                        std::inserter(updatedReadsInOldEdge, updatedReadsInOldEdge.begin()));
    e->reads.swap(updatedReadsInOldEdge);
    e->count = e->reads.size();
    if (e->count == 0) {
        e->sink->edgesIn.erase(e);
        e->source->edgesOut.erase(e->sink);
    }

    Node *oldCur = e->sink;
    // just create a new edge and return if we are arriving at a node in mainPath
    if (oldCur->onMainPath) {
        Edge *newEdge = createEdge(newPre, oldCur, 0);
        newEdge->reads = readsInPath2Split;
        newEdge->count = newEdge->reads.size();
        return;
    }

    // Otherwise create a new node
    Node *newCur = createNode(oldCur->base);
    // Add an edge to this new node
    Edge *newEdge = createEdge(newPre, newCur, 0);
    newEdge->reads = readsInPath2Split;
    newEdge->count = newEdge->reads.size();
    for (auto subEdge : oldCur->edgesOut) {
        splitPath(oldCur, newCur, subEdge.second, readsInPath2Split);
    }
}

ConsensusGraph::~ConsensusGraph() {
    for (auto &n : nodes)
        delete n;
    for (auto &e : edges)
        delete e;
    delete aligner;
}

Node *ConsensusGraph::createNode(char base) {
    Node *n = new Node(base);
    nodes.push_back(n);
    return n;
}

Edge *ConsensusGraph::createEdge(Node *source, Node *sink, size_t read) {
    Edge *e = new Edge(source, sink, read);
    edges.push_back(e);
    source->edgesOut.insert(std::make_pair(sink, e));
    sink->edgesIn.insert(e);
    return e;
}

void ConsensusGraph::printStatus() {
    std::cout << reads.size() << " reads, "
              << nodes.size() << " nodes, and "
              << edges.size() << " edges."
              << "\n";
    std::cout << "mainPath len " << mainPath.edges.size() + 1
              << " avg weight " << mainPath.getAverageWeight()
              << " starts at " << startPos
              << " ends at " << endPos
              << " len in Contig " << endPos - startPos
              << "\n";
    double stat = mainPath.edges.size() * mainPath.getAverageWeight() / (reads.size() * 9999);
    std::cout << stat / (1 - 0.17) << " "
              << stat / (1 - 0.17 * 1.05) << " "
              << stat / (1 - 0.17 * 1.1) << " "
              << std::endl;
}

void ConsensusGraph::writeMainPath(std::ofstream &f) {
    f << mainPath.path << std::endl;
}

void ConsensusGraph::writeReads(std::ofstream &f) {
    // First we write the index of each character into the cumulativeWeight
    // field of the nodes on mainPath
    mainPath.edges.front()->source->cumulativeWeight = 0;
    size_t i = 0;
    for (auto e : mainPath.edges) {
        e->sink->cumulativeWeight = ++i;
    }
    size_t totalEditDis = 0;
    for (auto it : reads) {
        totalEditDis += writeRead(f, it.second, it.first);
    }
    std::cout << "AvgEditDis " << totalEditDis / (double) reads.size() << std::endl;
}

size_t ConsensusGraph::writeRead(std::ofstream &f, Read &r, size_t id) {
    // First we write the read id and initial position
    Node *curNode = r.start;
    auto advanceInRead = [id](Node *n) -> Node * {
        for (auto e : n->edgesOut) {
            if (e.second->reads.find(id) != e.second->reads.end()) {
                return e.first;
            }
        }
        return NULL;
    };
    while (!curNode->onMainPath)
        curNode = advanceInRead(curNode);

    f << std::to_string(id) << ":" << std::to_string(curNode->cumulativeWeight) << "\n";

    size_t editDis = 0;
    size_t posInMainPath = curNode->cumulativeWeight;
    curNode = r.start;
    size_t unchangedCount = 0;
    auto dealWithUnchanged = [&f, &unchangedCount]() {
        if (unchangedCount > 0) {
            f << 'u' << std::to_string(unchangedCount);
            unchangedCount = 0;
        }
    };
    do {
        if (curNode->onMainPath) {
            size_t curPos = curNode->cumulativeWeight;
            if (curPos > posInMainPath)
                dealWithUnchanged();
            for (; posInMainPath < curPos; posInMainPath++) {
                f << 'd';
                editDis++;
            }
            unchangedCount++;
            posInMainPath++;
        } else {
            dealWithUnchanged();
            // Else we have an insertion
            f << 'i' << curNode->base;
            editDis++;
        }
    } while (curNode = advanceInRead(curNode));

    dealWithUnchanged();

    f << std::endl;
    return editDis;
}

ConsensusGraph::ConsensusGraph(StringAligner *aligner) : aligner(aligner) {}

ConsensusGraph::Read::Read(long pos, Node *start, size_t len) : pos(pos), start(start), len(len) {}
