#include "ConsensusGraph.h"
#include "bsc_helper.h"
#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <stack>

Edge::Edge(Node *source, Node *sink, size_t read) : source(source), sink(sink) {
    count = 1;
    // insert into sorted vector
    reads.insert(std::lower_bound(reads.begin(),reads.end(),read), read);
}

void Edge::addRead(size_t read) {
    count++;
    // insert into sorted vector
    reads.insert(std::lower_bound(reads.begin(),reads.end(),read), read);
}

Node::Node(const char base) : base(base) {}

// void Node::addEdge(Edge *e) {
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

Edge *Node::getBestEdgeOut() {
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

Edge *Node::getBestEdgeIn() {
    Edge *bestEdge = NULL;
    size_t bestCount = 0;
    const auto &end = edgesIn.end();
    for (auto it = edgesIn.begin(); it != end; ++it) {
        Edge *e = *it;
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
    return ((double)totalWeight) / edges.size();
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

void ConsensusGraph::initialize(const std::string &seed, size_t readId,
                                long pos) {

    unalignedReads.clear();
    size_t len = seed.length();
    Node *currentNode = createNode(seed[0]);
    // We create a read that points to this node
    readsInGraph.insert(std::make_pair(
        readId, ConsensusGraph::Read(pos, currentNode, seed.length())));
    startingNode = currentNode;
    rightMostUnchangedNode = currentNode;
    rightMostUnchangedNodeOffset = 0;
    leftMostUnchangedNode = currentNode;
    leftMostUnchangedNodeOffset = 0;
    mainPath.path.push_back(currentNode->base);
    currentNode->onMainPath = true;
    currentNode->cumulativeWeight = 0;
    for (size_t i = 1; i < len; ++i) {
        Node *nextNode = createNode(seed[i]);
        createEdge(currentNode, nextNode, readId);
        currentNode = nextNode;
    }
    startPos = pos;
    endPos = pos + len;
}

bool ConsensusGraph::addRead(const std::string &s, long pos,
                             std::vector<Edit> &editScript,
                             ssize_t &beginOffset, ssize_t &endOffset) {
    auto &originalString = mainPath.path;
    const ssize_t offsetGuess = originalString.size() - endPos + pos;
    size_t editDis;

    RAItA Abegin = originalString.begin();
    RAItA Aend = Abegin + originalString.size();
    RAItB Bbegin = s.c_str();
    RAItB Bend = Bbegin + s.length();
    bool success = aligner->align(Abegin, Aend, Bbegin, Bend, offsetGuess,
                                  beginOffset, endOffset, editScript, editDis);
    //    std::cout << "success ? " << success << std::endl;
    if (!success) {
        // std::cout << "Failed to add"
        //           << "\n";
        // std::cout << "mainPath startPos " << startPos << " mainPath endPos "
        //           << endPos << " current read startPos " << pos
        //           << " OffsetGuess " << offsetGuess << std::endl;
        return false;
    }
    return true;
    // updateGraph(s, editScript, beginOffset, endOffset, readId, pos);
}

void ConsensusGraph::addReads(
    const std::set<std::pair<long, read_t>> &reads,
    std::vector<std::unique_ptr<std::string>> &readData) {
    std::set<std::pair<long, read_t>> readsInContig(reads);
    auto currentRead = readsInContig.begin();

    initialize(*readData[currentRead->second], currentRead->second,
               currentRead->first);
    size_t len = (*readData[currentRead->second]).length();
    readsInContig.erase(currentRead);
    calculateMainPathGreedy();

    while (!readsInContig.empty()) {

        // First we lengthen mainPath by len/2
        auto read2Lengthen = readsInContig.lower_bound(
            std::make_pair(endPos - len + len / 2, 0));

        if (read2Lengthen != readsInContig.begin()) {
            read2Lengthen--;
        }

        {
            std::vector<Edit> editScript;
            ssize_t beginOffset, endOffset;
            read_t readId = read2Lengthen->second;
            ssize_t pos = read2Lengthen->first;
            std::string &s = *readData[readId];
            bool success = addRead(s, pos, editScript, beginOffset, endOffset);
            if (success) {
                updateGraph(s, editScript, beginOffset, endOffset, readId, pos);
                calculateMainPathGreedy();
            } else {
                unalignedReads.insert(std::make_pair(readId, s));
            }
            readsInContig.erase(read2Lengthen);
        }

        // Then we add all reads that should overlap with mainPath
        auto endRead2Add =
            readsInContig.lower_bound(std::make_pair(endPos - len - 100, 0));
        // size_t count = 0;
        std::vector<std::pair<long, read_t>> reads2Add(readsInContig.begin(),
                                                       endRead2Add);
        size_t num = reads2Add.size();
#pragma omp parallel for
        for (size_t i = 0; i < num; ++i) {
            auto read2Add = reads2Add[i];
            std::vector<Edit> editScript;
            ssize_t beginOffset, endOffset;
            read_t readId = read2Add.second;
            ssize_t pos = read2Add.first;
            std::string &s = *readData[readId];
            bool success = addRead(s, pos, editScript, beginOffset, endOffset);
#pragma omp critical
            {
                if (success)
                    updateGraph(s, editScript, beginOffset, endOffset, readId,
                                pos);
                else
                    unalignedReads.insert(std::make_pair(readId, s));
            }

            // count++;
            // if (count % 4 == 0)
            // calculateMainPath();
        }
        calculateMainPathGreedy();
        readsInContig.erase(readsInContig.begin(), endRead2Add);

        printStatus();
    }
}

void ConsensusGraph::updateGraph(const std::string &s,
                                 std::vector<Edit> &editScript,
                                 ssize_t beginOffset, ssize_t endOffset,
                                 size_t readId, long pos) {

    auto edgeInPath = mainPath.edges.begin();
    const auto &edgeInPathEnd = mainPath.edges.end();
    Node *nodeInPath = (*edgeInPath)->source;
    Node *currentNode = NULL;
    Node *initialNode = NULL;

    // First we update leftMostUnchangedNode and rightMostUnchangedNode
    if (beginOffset >= 0 || endOffset >= 0) {
        rightMostUnchangedNodeOffset = static_cast<size_t>(std::max(
            static_cast<ssize_t>(leftMostUnchangedNodeOffset),
            std::min(static_cast<ssize_t>(rightMostUnchangedNodeOffset),
                     beginOffset)));
        if (rightMostUnchangedNodeOffset > 0)
            rightMostUnchangedNode =
                mainPath.edges[rightMostUnchangedNodeOffset - 1]->sink;
        else
            rightMostUnchangedNode = mainPath.edges[0]->source;
    } else {
        leftMostUnchangedNodeOffset =
            std::min(rightMostUnchangedNodeOffset,
                     std::max(leftMostUnchangedNodeOffset,
                              mainPath.path.size() + endOffset));
        if (leftMostUnchangedNodeOffset > 0)
            leftMostUnchangedNode =
                mainPath.edges[leftMostUnchangedNodeOffset - 1]->sink;
        else
            leftMostUnchangedNode = mainPath.edges[0]->source;
    }

    auto advanceNodeInPath = [&]() {
        if (edgeInPath == edgeInPathEnd)
            return;
        nodeInPath = (*edgeInPath)->sink;
        edgeInPath++;
    };

    // First we deal with beginOffset
    auto initialAdvance = [&] {
        //        std::cout << "initialAdvance" << std::endl;
        if (beginOffset >= 1) {
            // We need to advance in the mainPath
            // for (size_t i = 0; i < beginOffset; i++) {
            //     advanceNodeInPath();
            // }
            edgeInPath += beginOffset - 1;
            nodeInPath = (*edgeInPath)->sink;
            edgeInPath++;

        } else {
            // We need to insert the initial parts of the read
            // This number must be positive
            size_t numOfNodes2Insert = -beginOffset;
            // std::cout << numOfNodes2Insert << " ";
            size_t i = 0;
            // We create an initial node
            currentNode = createNode(s[i++]);
            initialNode = currentNode;
            for (; i < numOfNodes2Insert; ++i) {
                Node *nextNode = createNode(s[i]);
                createEdge(currentNode, nextNode, readId);
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
    for (const Edit &e : editScript) {
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
    readsInGraph.insert(
        std::make_pair(readId, Read(pos, initialNode, s.length())));
}

// TODO: Optimize this to only update the portion of the graph that has changed
Path &ConsensusGraph::calculateMainPath() {
    mainPath.clear();

    // First we set all the hasReached fields to false
    clearHasReached(startingNode);

    std::deque<Node *> unfinishedNodes;
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
            if (!n->hasReached) {
                // std::deque<Node *>::iterator temp = std::find(
                //     unfinishedNodes.begin(), unfinishedNodes.end(), n);
                // if (temp != unfinishedNodes.end()) {
                //     std::cout << std::distance(unfinishedNodes.begin(),
                //     temp)
                //               << std::endl;
                //     std::cout << unfinishedNodes.size() << std::endl;
                //     std::raise(SIGINT);
                // }
                // if (unfinishedNodes.size() > 1000000) {
                //     std::cout << unfinishedNodes.size() << ' ';
                //     std::raise(SIGINT);
                // }
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
        globalMaxWeightNode =
            maxWeight >= globalMaxWeight ? currentNode : globalMaxWeightNode;
        globalMaxWeight = std::max(maxWeight, globalMaxWeight);
        unfinishedNodes.pop_front();
        for (auto edgeIt : currentNode->edgesIn) {
            unfinishedNodes.push_back(edgeIt->source);
        }
        currentNode->cumulativeWeight = maxWeight;
        currentNode->hasReached = true;
    }

    startingNode = globalMaxWeightNode;

    auto &edgesInPath = mainPath.edges;
    auto &stringPath = mainPath.path;

    stringPath.push_back(globalMaxWeightNode->base);
    globalMaxWeightNode->onMainPath = true;
    //    std::cout << "max weight" << globalMaxWeight << std::endl;
    //    std::cout << "here" << std::endl;
    while (globalMaxWeight > 0) {
        //        std::cout << "max weight" << globalMaxWeight << " "
        //                  << edgesInPath.size() << std::endl;
        for (auto e : globalMaxWeightNode->edgesOut) {
            if (e.second->count + e.first->cumulativeWeight ==
                globalMaxWeight) {
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
    startPos = readsInGraph.at(startingReadId).pos;
    size_t endingReadId = *edgesInPath.back()->reads.begin();
    Read &endingRead = readsInGraph.at(endingReadId);
    endPos = endingRead.pos + endingRead.len;
    //    printStatus();
    removeCycles();
    return mainPath;
}

void ConsensusGraph::clearHasReached(Node *n) {
    std::deque<Node *> unfinishedNodes;
    unfinishedNodes.push_back(n);
    while (!unfinishedNodes.empty()) {
        Node *currentNode = unfinishedNodes.front();
        // We don't do anything if it has already been set to false
        if (!currentNode->hasReached) {
            unfinishedNodes.pop_front();
            continue;
        }
        // if (unfinishedNodes.size() > 1000000) {
        //     std::cout << unfinishedNodes.size() << ' ';
        //     std::raise(SIGINT);
        // }
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
}

Path &ConsensusGraph::calculateMainPathGreedy() {
    clearMainPath();

    auto &edgesInPath = mainPath.edges;
    auto &stringPath = mainPath.path;

    // Extend to the right
    {
        Node *currentNode = rightMostUnchangedNode;
        currentNode->onMainPath = true;
        Edge *edgeToAdd;
        while ((edgeToAdd = currentNode->getBestEdgeOut())) {
            edgesInPath.push_back(edgeToAdd);
            currentNode = edgeToAdd->sink;
            currentNode->onMainPath = true;
            stringPath.push_back(currentNode->base);
        }
        size_t endingReadId = *edgesInPath.back()->reads.begin();
        Read &endingRead = readsInGraph.at(endingReadId);
        endPos = endingRead.pos + endingRead.len;
    }

    // Extend to the left
    {
        Node *currentNode = leftMostUnchangedNode;
        currentNode->onMainPath = true;
        Edge *edgeToAdd;
        while ((edgeToAdd = currentNode->getBestEdgeIn())) {
            edgesInPath.insert(edgesInPath.begin(), edgeToAdd);
            currentNode = edgeToAdd->source;
            currentNode->onMainPath = true;
            stringPath.insert(stringPath.begin(), currentNode->base);
            leftMostUnchangedNodeOffset++;
            rightMostUnchangedNodeOffset++;
        }
        size_t startingReadId = *edgesInPath.front()->reads.begin();
        startPos = readsInGraph.at(startingReadId).pos;
    }

    // rightMostUnchangedNodeOffset = 0;
    // rightMostUnchangedNode = edgesInPath.front()->source;
    // printStatus();
    removeCycles();
    rightMostUnchangedNode = edgesInPath.back()->sink;
    rightMostUnchangedNodeOffset = edgesInPath.size() + 1;
    leftMostUnchangedNode = edgesInPath.front()->source;
    leftMostUnchangedNodeOffset = 0;
    return mainPath;
}

void ConsensusGraph::clearMainPath() {
    // printStatus();
    // std::cout << "clearing\n";
    size_t l = mainPath.edges.size();
    // We first erase the right tail
    for (size_t i = rightMostUnchangedNodeOffset; i < l; ++i) {
        Node *currentNode = mainPath.edges[i]->sink;
        currentNode->onMainPath = false;
    }
    if (rightMostUnchangedNodeOffset < mainPath.edges.size()) {
        mainPath.edges.erase(mainPath.edges.begin() +
                                 rightMostUnchangedNodeOffset,
                             mainPath.edges.end());
    }
    if (mainPath.path.size() > rightMostUnchangedNodeOffset + 1)
        mainPath.path.erase(mainPath.path.begin() +
                                rightMostUnchangedNodeOffset + 1,
                            mainPath.path.end());
    // Now we erase the left tail
    for (size_t i = 0; i < leftMostUnchangedNodeOffset; ++i) {
        Node *currentNode = mainPath.edges[i]->source;
        currentNode->onMainPath = false;
    }
    if (leftMostUnchangedNodeOffset > 0 &&
        leftMostUnchangedNodeOffset <= mainPath.edges.size()) {
        mainPath.edges.erase(mainPath.edges.begin(),
                             mainPath.edges.begin() +
                                 leftMostUnchangedNodeOffset);
        mainPath.path.erase(mainPath.path.begin(),
                            mainPath.path.begin() +
                                leftMostUnchangedNodeOffset);
        rightMostUnchangedNodeOffset -= leftMostUnchangedNodeOffset;
    }
    leftMostUnchangedNodeOffset = 0;
    // printStatus();
}

void ConsensusGraph::removeCycles() {
    // We first iterate over all nodes on mainPath on the right
    auto edgeOnPath = mainPath.edges.begin() + rightMostUnchangedNodeOffset;
    auto edgeOnPathEnd = mainPath.edges.end();
    Node *nodeOnPath =
        edgeOnPath < edgeOnPathEnd ? (*edgeOnPath)->source : nullptr;
    while (edgeOnPath < edgeOnPathEnd) {
        // Then we iterate over all edges pointing to side nodes that have
        // other edges in
        auto end = nodeOnPath->edgesOut.end();
        for (auto edgeIt = nodeOnPath->edgesOut.begin(); edgeIt != end;) {
            // Node *sideNode = edgeIt->first;
            Edge *e = edgeIt->second;
            edgeIt++;
            walkAndPrune(e);
        }
        if (edgeOnPath == edgeOnPathEnd)
            break;
        nodeOnPath = (*edgeOnPath)->sink;
        ++edgeOnPath;
    }
    // Now we iterate over all nodes on mainPath on the left
    for (edgeOnPath = mainPath.edges.begin() + leftMostUnchangedNodeOffset;
         edgeOnPath >= mainPath.edges.begin(); --edgeOnPath) {
        // Then we iterate over all edges pointing to side nodes that have
        // other edges in
        nodeOnPath = (*edgeOnPath)->source;
        auto end = nodeOnPath->edgesOut.end();
        for (auto edgeIt = nodeOnPath->edgesOut.begin(); edgeIt != end;) {
            // Node *sideNode = edgeIt->first;
            Edge *e = edgeIt->second;
            edgeIt++;
            walkAndPrune(e);
        }
    }
}

void ConsensusGraph::walkAndPrune(Edge *e) {
    Node *sink = e->sink;
    Node *source = e->source;
    if (sink->onMainPath)
        return;
    if (sink->edgesIn.size() > 1)
        splitPath(source, e, e->reads);
    auto end = sink->edgesOut.end();
    for (auto subEdge = sink->edgesOut.begin(); subEdge != end;) {
        Edge *edge2WorkOn = subEdge->second;
        subEdge++;
        // Here walkAndPrune will only change sink->edgesOut by 1. deleting
        // edge2WorkOn and 2. adding new edges that we don't need to prune. So
        // copying sink->edgesOut should work.
        walkAndPrune(edge2WorkOn);
    }
}

void ConsensusGraph::splitPath(Node *newPre, Edge *e,
                               std::vector<size_t> const &reads2Split) {
    // The reads that need to be split going down this path
    std::vector<size_t> readsInPath2Split;
    std::set_intersection(
        reads2Split.begin(), reads2Split.end(), e->reads.begin(),
        e->reads.end(),
        std::inserter(readsInPath2Split, readsInPath2Split.begin()));
    if (readsInPath2Split.empty())
        return;

    Node *oldCur = e->sink;

    // Remove the reads from the node edge
    removeReadsFromEdge(e, readsInPath2Split);

    // just create a new edge and return if we are arriving at a node in
    // mainPath
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

    auto end = oldCur->edgesOut.end();
    for (auto subEdge = oldCur->edgesOut.begin(); subEdge != end;) {
        Edge *edge2WorkOn = subEdge->second;
        subEdge++;
        // Here, splitPath will only affect oldCur->edgesOut by deleting
        // edge2WorkOn. So again just copying everything should work.
        splitPath(newCur, edge2WorkOn, readsInPath2Split);
    }
    if (oldCur->edgesIn.empty() && oldCur->edgesOut.empty())
        removeNode(oldCur);
}

ConsensusGraph::~ConsensusGraph() {
    if (!startingNode)
        return;
    // std::cerr << "Removing" << std::endl;

    removeConnectedNodes(startingNode);
    // std::cerr << std::to_string(numEdges) << " edges "
    //           << std::to_string(numNodes) << " nodes left\n";
}

Node *ConsensusGraph::createNode(char base) {
    Node *n = new Node(base);
    numNodes++;
    return n;
}

Edge *ConsensusGraph::createEdge(Node *source, Node *sink, size_t read) {
    Edge *e = new Edge(source, sink, read);
    source->edgesOut.insert(std::make_pair(sink, e));
    sink->edgesIn.insert(e);
    numEdges++;
    return e;
}

void ConsensusGraph::removeReadsFromEdge(Edge *e,
                                         std::vector<size_t> const &reads) {
    std::vector<size_t> updatedReadsInOldEdge;
    std::set_difference(
        e->reads.begin(), e->reads.end(), reads.begin(), reads.end(),
        std::inserter(updatedReadsInOldEdge, updatedReadsInOldEdge.begin()));
    e->reads.swap(updatedReadsInOldEdge);
    e->count = e->reads.size();
    if (e->count == 0) {
        removeEdge(e);
    }
}

void ConsensusGraph::removeEdge(Edge *e) {
    e->source->edgesOut.erase(e->sink);
    e->sink->edgesIn.erase(e);
    delete e;
    numEdges--;
}

void ConsensusGraph::removeNode(Node *n) {
    {
        std::set<Edge *>::iterator edgeIt = n->edgesIn.begin();
        std::set<Edge *>::iterator end = n->edgesIn.end();
        while (edgeIt != end)
            removeEdge(*(edgeIt++));
    }
    {
        std::map<Node *, Edge *>::iterator edgeIt = n->edgesOut.begin();
        std::map<Node *, Edge *>::iterator end = n->edgesOut.end();
        while (edgeIt != end)
            removeEdge(edgeIt++->second);
    }
    delete n;
    numNodes--;
}

void ConsensusGraph::removeBelow(Node *n) {
    std::stack<Node *> nodesToRemove;
    nodesToRemove.push(n);
    while (!nodesToRemove.empty()) {
        Node *curNode = nodesToRemove.top();
        if (curNode->edgesOut.empty()) {
            nodesToRemove.pop();
            removeNode(curNode);
        } else {
            nodesToRemove.push(curNode->edgesOut.begin()->first);
        }
    }
}

void ConsensusGraph::removeAbove(Node *n) {
    std::stack<Node *> nodesToRemove;
    nodesToRemove.push(n);
    while (!nodesToRemove.empty()) {
        Node *curNode = nodesToRemove.top();
        if (curNode->edgesIn.empty()) {
            nodesToRemove.pop();
            removeNode(curNode);
        } else {
            nodesToRemove.push((*curNode->edgesIn.begin())->source);
        }
    }
}

void ConsensusGraph::removeConnectedNodes(Node *n) {
    clearHasReached(n);
    std::stack<Node *> nodesToRemove;
    nodesToRemove.push(n);
    n->hasReached = true;
    while (!nodesToRemove.empty()) {
        Node *curNode = nodesToRemove.top();
        if (curNode->edgesOut.empty()) {
            nodesToRemove.pop();
            for (auto edgeIt : curNode->edgesIn) {
                Node *source = edgeIt->source;
                if (!source->hasReached && source->edgesOut.size() == 1) {
                    nodesToRemove.push(source);
                    source->hasReached = true;
                }
            }
            removeNode(curNode);
        } else {
            Node *node2Add = curNode->edgesOut.begin()->first;
            nodesToRemove.push(node2Add);
            node2Add->hasReached = true;
        }
    }
}

void ConsensusGraph::printStatus() {
    std::cout << unalignedReads.size() << " unaligned reads "
              << readsInGraph.size() << " reads in graph " << this << ", "
              << numNodes << " nodes, and " << numEdges << " edges."
              << "\n";
    std::cout << "mainPath len " << mainPath.edges.size() + 1 << " "
              << mainPath.path.size() << " avg weight "
              << mainPath.getAverageWeight() << " starts at " << startPos
              << " ends at " << endPos << " len in Contig " << endPos - startPos
              << "\n";
    double stat = mainPath.edges.size() * mainPath.getAverageWeight() /
                  (readsInGraph.size() * 9999);
    std::cout << stat / (1 - 0.17) << " " << stat / (1 - 0.17 * 1.05) << " "
              << stat / (1 - 0.17 * 1.1) << " " << std::endl;
}

void ConsensusGraph::writeMainPath(const std::string &filename) {
    std::ofstream f;
    std::string genomeFileName = tempDir + filename + ".genome";
    f.open(genomeFileName);
    f << std::string(mainPath.path.begin(), mainPath.path.end()) << std::endl;
    f.close();
    // std::string compressedGenomeFileName =
    //     compressedTempDir + filename + ".genomeCompressed";
    // bsc::BSC_compress(genomeFileName.c_str(),
    // compressedGenomeFileName.c_str());
}

void ConsensusGraph::writeReads(const std::string &filename) {
    // First we write the index of each character into the cumulativeWeight
    // field of the nodes on mainPath
    mainPath.edges.front()->source->cumulativeWeight = 0;
    size_t i = 0;
    for (auto e : mainPath.edges) {
        e->sink->cumulativeWeight = ++i;
    }
    size_t totalEditDis = 0;

    const std::string posFileName = tempDir + filename + ".pos";
    const std::string editTypeFileName = tempDir + filename + ".type";
    const std::string editBaseFileName = tempDir + filename + ".base";
    const std::string idFileName = tempDir + filename + ".id";
    std::ofstream posFile, editTypeFile, editBaseFile, idFile;
    posFile.open(posFileName);
    editTypeFile.open(editTypeFileName);
    editBaseFile.open(editBaseFileName);
    idFile.open(idFileName);
    size_t pasId = 0;
    for (auto it : readsInGraph) {
        {
            idFile << it.first - pasId << ':';
            pasId = it.first;
        }
        totalEditDis +=
            writeRead(posFile, editTypeFile, editBaseFile, it.second, it.first);
    }
    idFile << '\n';
    posFile.close();
    editTypeFile.close();
    editBaseFile.close();
    idFile.close();
    std::cout << "AvgEditDis " << totalEditDis / (double)readsInGraph.size()
              << std::endl;
    printStatus();
    // const std::string posFileCompressedName =
    //     compressedTempDir + filename + ".pos" + "Compressed";
    // const std::string editTypeFileCompressedName =
    //     compressedTempDir + filename + ".type" + "Compressed";
    // const std::string editBaseFileCompressedName =
    //     compressedTempDir + filename + ".base" + "Compressed";
    // const std::string idFileCompressedName =
    //     compressedTempDir + filename + ".id" + "Compressed";
    // bsc::BSC_compress(posFileName.c_str(), posFileCompressedName.c_str());
    // bsc::BSC_compress(editTypeFileName.c_str(),
    //                   editTypeFileCompressedName.c_str());
    // bsc::BSC_compress(editBaseFileName.c_str(),
    //                   editBaseFileCompressedName.c_str());
    // bsc::BSC_compress(idFileName.c_str(), idFileCompressedName.c_str());

    // writeUnalignedReads(filename);
}

void ConsensusGraph::writeUnalignedReads(const std::string &filename) {
    const std::string unalignedReadsFileName =
        tempDir + filename + ".unalignedReads";
    const std::string unalignedIdsFileName =
        tempDir + filename + ".unalignedIds";
    std::ofstream unalignedReadsFile, unalignedIdsFile;
    unalignedReadsFile.open(unalignedReadsFileName);
    unalignedIdsFile.open(unalignedIdsFileName);
    size_t pastId = 0;
    for (auto it : unalignedReads) {
        unalignedIdsFile << it.first - pastId << ":";
        pastId = it.first;
        unalignedReadsFile << it.second << '\n';
    }
    unalignedIdsFile.close();
    unalignedReadsFile.close();
    // const std::string unalignedReadsFileCompressedName =
    //     compressedTempDir + filename + ".unalignedReads" + "Compressed";
    // const std::string unalignedIdsFileCompressedName =
    //     compressedTempDir + filename + ".unalignedIds" + "Compressed";
    // bsc::BSC_compress(unalignedReadsFileName.c_str(),
    //                   unalignedReadsFileCompressedName.c_str());
    // bsc::BSC_compress(unalignedIdsFileName.c_str(),
    //                   unalignedIdsFileCompressedName.c_str());
}

size_t ConsensusGraph::getNumReads() { return readsInGraph.size(); }

size_t ConsensusGraph::read2EditScript(ConsensusGraph::Read &r, size_t id,
                                       std::vector<Edit> &editScript,
                                       size_t &pos) {
    editScript.clear();
    // First we store the initial position
    Node *curNode = r.start;
    auto advanceInRead = [id](Node *n) -> Node * {
        for (auto e : n->edgesOut) {
            if (std::binary_search(e.second->reads.begin(),e.second->reads.end(),id)) {
                return e.first;
            }
        }
        return NULL;
    };
    while (!curNode->onMainPath)
        curNode = advanceInRead(curNode);

    pos = curNode->cumulativeWeight;

    size_t editDis = 0;
    size_t posInMainPath = curNode->cumulativeWeight;
    curNode = r.start;
    size_t unchangedCount = 0;
    auto dealWithUnchanged = [&unchangedCount, &editScript]() {
        if (unchangedCount > 0) {
            editScript.push_back(Edit(SAME, unchangedCount));
            unchangedCount = 0;
        }
    };
    do {
        if (curNode->onMainPath) {
            size_t curPos = curNode->cumulativeWeight;
            if (curPos > posInMainPath)
                dealWithUnchanged();
            for (; posInMainPath < curPos; posInMainPath++) {
                editScript.push_back(Edit(DELETE, '-'));
                editDis++;
            }
            unchangedCount++;
            posInMainPath++;
        } else {
            dealWithUnchanged();
            // Else we have an insertion
            editScript.push_back(Edit(INSERT, curNode->base));
            editDis++;
        }
    } while ((curNode = advanceInRead(curNode)));

    dealWithUnchanged();

    return editDis;
}

size_t ConsensusGraph::writeRead(std::ofstream &posFile,
                                 std::ofstream &editTypeFile,
                                 std::ofstream &editBaseFile, Read &r,
                                 size_t id) {

    size_t offset;
    std::vector<Edit> editScript;
    size_t editDis = read2EditScript(r, id, editScript, offset);
    posFile << offset << ':';

    std::vector<Edit> newEditScript;
    editDis = Edit::optimizeEditScript(editScript, newEditScript);
    size_t unchangedCount = 0;
    for (Edit e : newEditScript) {
        switch (e.editType) {
        case SAME: {
            unchangedCount += e.editInfo.num;
            break;
        }
        case INSERT: {
            posFile << unchangedCount << ':';
            unchangedCount = 0;
            editTypeFile << 'i';
            editBaseFile << e.editInfo.ins;
            break;
        }
        case DELETE: {
            posFile << unchangedCount << ':';
            unchangedCount = 0;
            editTypeFile << 'd';
            break;
        }
        case SUBSTITUTION: {
            posFile << unchangedCount << ':';
            unchangedCount = 0;
            editTypeFile << 's';
            editBaseFile << e.editInfo.sub;
            break;
        }
        }
    }

    posFile << unchangedCount << ':';
    unchangedCount = 0;

    posFile << '\n';
    editTypeFile << '\n';
    editBaseFile << '\n';
    return editDis;
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
    for (auto it : readsInGraph) {
        totalEditDis += writeRead(f, it.second, it.first);
    }
    std::cout << "AvgEditDis " << totalEditDis / (double)readsInGraph.size()
              << std::endl;
}

size_t ConsensusGraph::writeRead(std::ofstream &f, Read &r, size_t id) {

    std::vector<Edit> editScript;
    size_t pos;
    size_t editDis = read2EditScript(r, id, editScript, pos);
    std::vector<Edit> newEditScript;

    editDis = Edit::optimizeEditScript(editScript, newEditScript);

    f << std::to_string(id) << ":" << std::to_string(pos) << "\n";

    for (Edit e : newEditScript) {
        switch (e.editType) {
        case SAME:
            f << 'u' << std::to_string(e.editInfo.num);
            break;
        case DELETE:
            f << 'd';
            break;
        case INSERT:
            f << 'i' << e.editInfo.ins;
            break;
        case SUBSTITUTION:
            f << 's' << e.editInfo.sub;
            break;
        }
    }

    f << std::endl;
    return editDis;
}

ConsensusGraph::ConsensusGraph(StringAligner_t *aligner) : aligner(aligner) {}

ConsensusGraph::Read::Read(long pos, Node *start, size_t len)
    : pos(pos), start(start), len(len) {}
