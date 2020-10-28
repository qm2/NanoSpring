#include "ConsensusGraph.h"
#include "bsc_helper.h"
#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <stack>

Edge::Edge(Node *source, Node *sink, read_t read) : source(source), sink(sink) {
    count = 1;
    // insert into sorted vector
    reads.insert(std::lower_bound(reads.begin(), reads.end(), read), read);
}

void Edge::addRead(read_t read) {
    count++;
    // insert into sorted vector
    reads.insert(std::lower_bound(reads.begin(), reads.end(), read), read);
}

Node::Node(const char base) : base(base) {}

// void Node::addEdge(Edge *e) {
//    edgesOut[e->sink] = e;
//}

Edge *Node::getEdgeTo(Node *n) {
    const auto &it = std::find_if(
        edgesOut.begin(), edgesOut.end(),
        [&](const std::pair<Node *, Edge *> &p) { return (p.first == n); });
    if (it == edgesOut.end()) {
        return nullptr;
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
    return nullptr;
}

Edge *Node::getBestEdgeOut() {
    Edge *bestEdge = nullptr;
    read_t bestCount = 0;
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
    Edge *bestEdge = nullptr;
    read_t bestCount = 0;
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

Edge *Node::getEdgeInRead(read_t read) const {
    for (auto e : edgesOut) {
        if (std::binary_search(e.second->reads.begin(), e.second->reads.end(),
                               read)) {
            return e.second;
        }
    }
    return nullptr;
}

Node *Node::getNextNodeInRead(read_t read) const {
    Edge *e = getEdgeInRead(read);
    return e ? e->sink : nullptr;
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

void ConsensusGraph::initialize(const std::string &seed, read_t readId,
                                long pos) {
    size_t len = seed.length();
    Node *currentNode = createNode(seed[0]);
    // We create a read that points to this node
    readsInGraph.insert(std::make_pair(
        readId, ConsensusGraph::Read(pos, currentNode, seed.length(), false)));
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
    endPos = pos + 1;
}

bool ConsensusGraph::addRead(const std::string &s, long pos,
                             std::vector<Edit> &editScript,
                             ssize_t &beginOffset, ssize_t &endOffset) {
    auto &originalString = mainPath.path;
    /// TODO: don't always use endPos for offsetGuess
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
    size_t numUnchanged = 0;
    for (auto e : editScript)
        if (e.editType == SAME)
            numUnchanged += e.editInfo.num;
    if (numUnchanged == 0)
        return false;
    return true;
    // updateGraph(s, editScript, beginOffset, endOffset, readId, pos);
}

void ConsensusGraph::updateGraph(const std::string &s,
                                 std::vector<Edit> &editScript,
                                 ssize_t beginOffset, ssize_t endOffset,
                                 read_t readId, long pos,
                                 bool reverseComplement) {

    const auto &edgeInPathEnd = mainPath.edges.end();
    auto edgeInPath = mainPath.edges.begin();
    /** nodeInPath should always be set to edgeInPath.source, unless it is the
    last node, in which case edgeInPath should be edgeInPathEnd **/
    Node *nodeInPath = (*edgeInPath)->source;
    /** The last node of this read that has been added **/
    Node *currentNode = nullptr;
    Node *initialNode = nullptr;

    size_t numUnchanged = 0;

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
                              mainPath.path.size() - 1 + endOffset));
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

        } else if (beginOffset <= -1) {
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
            numUnchanged += num;
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

    assert(numUnchanged > 0);
    // Don't forget to add the read!
    readsInGraph.insert(std::make_pair(
        readId, Read(pos, initialNode, s.length(), reverseComplement)));
}

// Deprecated
// Path &ConsensusGraph::calculateMainPath() {
// mainPath.clear();

//// First we set all the hasReached fields to false
// traverseAndCall(startingNode, true, [](Node *) {});

// std::deque<Node *> unfinishedNodes;
// size_t globalMaxWeight = 0;
// Node *globalMaxWeightNode = nullptr;
// unfinishedNodes.push_back(startingNode);
// while (!unfinishedNodes.empty()) {
// Node *currentNode = unfinishedNodes.front();
// if (currentNode->hasReached) {
// unfinishedNodes.pop_front();
// continue;
//}

// size_t maxWeight = 0;
// bool hasPreviousWeights = true;
// for (auto edgeIt : currentNode->edgesOut) {
// Node *n = edgeIt.second->sink;
// if (!n->hasReached) {
//// std::deque<Node *>::iterator temp = std::find(
////     unfinishedNodes.begin(), unfinishedNodes.end(), n);
//// if (temp != unfinishedNodes.end()) {
////     std::cout << std::distance(unfinishedNodes.begin(),
////     temp)
////               << std::endl;
////     std::cout << unfinishedNodes.size() << std::endl;
////     std::raise(SIGINT);
//// }
//// if (unfinishedNodes.size() > 1000000) {
////     std::cout << unfinishedNodes.size() << ' ';
////     std::raise(SIGINT);
//// }
// unfinishedNodes.push_front(n);
// hasPreviousWeights = false;
// break;
//}
// size_t prevWeight = edgeIt.second->sink->cumulativeWeight;
// size_t curWeight = prevWeight + edgeIt.second->count;
// maxWeight = std::max(maxWeight, curWeight);
//}
// if (!hasPreviousWeights)
// continue;
// globalMaxWeightNode =
// maxWeight >= globalMaxWeight ? currentNode : globalMaxWeightNode;
// globalMaxWeight = std::max(maxWeight, globalMaxWeight);
// unfinishedNodes.pop_front();
// for (auto edgeIt : currentNode->edgesIn) {
// unfinishedNodes.push_back(edgeIt->source);
//}
// currentNode->cumulativeWeight = maxWeight;
// currentNode->hasReached = true;
//}

// startingNode = globalMaxWeightNode;

// auto &edgesInPath = mainPath.edges;
// auto &stringPath = mainPath.path;

// stringPath.push_back(globalMaxWeightNode->base);
// globalMaxWeightNode->onMainPath = true;
////    std::cout << "max weight" << globalMaxWeight << std::endl;
////    std::cout << "here" << std::endl;
// while (globalMaxWeight > 0) {
////        std::cout << "max weight" << globalMaxWeight << " "
////                  << edgesInPath.size() << std::endl;
// for (auto e : globalMaxWeightNode->edgesOut) {
// if (e.second->count + e.first->cumulativeWeight ==
// globalMaxWeight) {
// edgesInPath.push_back(e.second);
// globalMaxWeight -= e.second->count;
// globalMaxWeightNode = e.first;
// stringPath.push_back(globalMaxWeightNode->base);
// e.first->onMainPath = true;
// break;
//};
//}
//}

// read_t startingReadId = *edgesInPath.front()->reads.begin();
// startPos = readsInGraph.at(startingReadId).pos;
// read_t endingReadId = *edgesInPath.back()->reads.begin();
// Read &endingRead = readsInGraph.at(endingReadId);
// endPos = endingRead.pos + endingRead.len;
////    printStatus();
// removeCycles();
// return mainPath;
//}

Path &ConsensusGraph::calculateMainPathGreedy() {
    clearMainPath();

    auto &edgesInPath = mainPath.edges;
    auto &stringPath = mainPath.path;

    // Extend to the right
    {
        Node *currentNode = rightMostUnchangedNode;
        assert(currentNode->onMainPath);
        Edge *edgeToAdd;
        while ((edgeToAdd = currentNode->getBestEdgeOut())) {
            edgesInPath.push_back(edgeToAdd);
            currentNode = edgeToAdd->sink;
            currentNode->onMainPath = true;
            stringPath.push_back(currentNode->base);
        }
        read_t endingReadId = *edgesInPath.back()->reads.begin();
        Read &endingRead = readsInGraph.at(endingReadId);
        endPos = endingRead.pos + endingRead.len;
    }

    // Extend to the left
    {
        Node *currentNode = leftMostUnchangedNode;
        assert(currentNode->onMainPath);
        Edge *edgeToAdd;
        while ((edgeToAdd = currentNode->getBestEdgeIn())) {
            edgesInPath.insert(edgesInPath.begin(), edgeToAdd);
            currentNode = edgeToAdd->source;
            currentNode->onMainPath = true;
            stringPath.insert(stringPath.begin(), currentNode->base);
            leftMostUnchangedNodeOffset++;
            rightMostUnchangedNodeOffset++;
        }
        read_t startingReadId = *edgesInPath.front()->reads.begin();
        startPos = readsInGraph.at(startingReadId).pos;
    }

    // rightMostUnchangedNodeOffset = 0;
    // rightMostUnchangedNode = edgesInPath.front()->source;
    // printStatus();
    removeCycles();
    rightMostUnchangedNode = edgesInPath.back()->sink;
    rightMostUnchangedNodeOffset = edgesInPath.size();
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
    if (leftMostUnchangedNodeOffset > 0) {
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
    Node *nodeOnPath = edgeOnPath < edgeOnPathEnd ? (*edgeOnPath)->source
                                                  : edgeOnPath[-1]->sink;
    while (true) {
        // Then we iterate over all edges pointing to side nodes that have
        // other edges in
        //
        // Work with copy to avoid iterator invalidation.
        auto edgesOutCopy = nodeOnPath->edgesOut;
        for (const auto &edgeIt : edgesOutCopy)
            walkAndPrune(edgeIt.second);

        if (edgeOnPath == edgeOnPathEnd)
            break;
        nodeOnPath = (*edgeOnPath)->sink;
        ++edgeOnPath;
    }
    // Now we iterate over all nodes on mainPath on the left
    edgeOnPath = mainPath.edges.begin() + leftMostUnchangedNodeOffset;
    while (true) {
        // Then we iterate over all edges pointing to side nodes that have
        // other edges in
        nodeOnPath = (*edgeOnPath)->source;
        // Work with copy to avoid iterator invalidation.
        auto edgesOutCopy = nodeOnPath->edgesOut;
        for (const auto &edgeIt : edgesOutCopy)
            walkAndPrune(edgeIt.second);

        if (edgeOnPath == mainPath.edges.begin())
            break;
        --edgeOnPath;
    }
}

void ConsensusGraph::walkAndPrune(Edge *e) {
    std::stack<Edge *> callStack; // converting recursion to iteration
    callStack.push(e);
    while (!callStack.empty()) {
        Edge *curr = callStack.top();
        callStack.pop();
        Node *sink = curr->sink;
        Node *source = curr->source;
        if (sink->onMainPath)
            continue;
        if (sink->edgesIn.size() > 1)
            splitPath(source, curr, &(curr->reads));
        // Now put sink->edgesOut into stack, use reverse order for consistency with 
        // earlier code
        for (auto it = sink->edgesOut.rbegin(); it != sink->edgesOut.rend(); ++it) {
            callStack.push(it->second);
            // OLD COMMENT:
            // Here walkAndPrune will only change sink->edgesOut by 1. deleting
            // edge2WorkOn and 2. adding new edges that we don't need to prune. So
            // copying sink->edgesOut should work.
        }
    }
}

void ConsensusGraph::splitPath(Node *newPre_, Edge *e_,
                               std::vector<read_t> *reads2Split_) {
    // converting recursion to iteration.
    std::stack<std::tuple<Node *,Edge *, std::vector<read_t> *, bool, Node *>> callStack;
    // callStack contains:
    // 1. Node *newPre
    // 2. Edge *e
    // 3. std::vector<read_t> *reads2Split
    // 4. bool value denoting if this is the last subedge in oldCur->edgesOut iteration
    // 5. Node *oldCurToRemove // only needed when bool value is true
    // When the bool value is true, two things happen at the end (cleanup for calling iteration):
    // i. the reads2Split vector is deleted
    // ii. oldCur node is removed if it has no more edges
    // If oldCur->edgesOut is empty, these happen in the current iteration itself

    callStack.push(std::make_tuple(newPre_,e_,reads2Split_,false,(Node*)NULL));
    
    while (!callStack.empty()) {
        auto newPre = std::get<0>(callStack.top());
        auto e = std::get<1>(callStack.top());
        auto reads2Split = std::get<2>(callStack.top());
        auto lastIterationFlag = std::get<3>(callStack.top());
        auto oldCurToRemove = std::get<4>(callStack.top());
        callStack.pop();
        
        // The reads that need to be split going down this path
        auto readsInPath2Split = new std::vector<read_t>;
        std::set_intersection(
            reads2Split->begin(), reads2Split->end(), e->reads.begin(),
            e->reads.end(),
            std::inserter(*readsInPath2Split, readsInPath2Split->begin()));
        if (readsInPath2Split->empty()) {
            delete readsInPath2Split;
        } else {
            Node *oldCur = e->sink;

            // Remove the reads from the node edge
            removeReadsFromEdge(e, *readsInPath2Split);

            // just create a new edge and return if we are arriving at a node in
            // mainPath
            if (oldCur->onMainPath) {
                Edge *newEdge = createEdge(newPre, oldCur, 0);
                newEdge->reads = *readsInPath2Split;
                newEdge->count = newEdge->reads.size();
                delete readsInPath2Split;
            } else { 
                // Otherwise create a new node
                Node *newCur = createNode(oldCur->base);
                // Add an edge to this new node
                Edge *newEdge = createEdge(newPre, newCur, 0);
                newEdge->reads = *readsInPath2Split;
                newEdge->count = newEdge->reads.size();

                // Now put sink->edgesOut into stack, use reverse order for consistency with 
                // earlier code
                auto numElemsToPush = oldCur->edgesOut.size();
                if (numElemsToPush > 0) {
                    for (size_t idx = numElemsToPush - 1; idx != (size_t)(-1); --idx) {
                        if (idx == numElemsToPush - 1) {
                            // this is the last one that will be called
                            // So we set lastIterationFlag to true, and include oldCur
                            callStack.push(std::make_tuple(newCur,oldCur->edgesOut[idx].second,readsInPath2Split,true,oldCur));
                        } else {
                            callStack.push(std::make_tuple(newCur,oldCur->edgesOut[idx].second,readsInPath2Split,false,(Node*)NULL));
                        }
                        // OLD COMMENT:
                        // Here, splitPath will only affect oldCur->edgesOut by deleting
                        // edge2WorkOn. So again just copying everything should work.
                    }
                } else {
                    // no edge out so perform cleanup here itself
                    delete readsInPath2Split;
                    if (oldCur->edgesIn.empty() && oldCur->edgesOut.empty())
                        removeNode(oldCur);
                }
            }
        }
        // finally, if lastIterationFlag is set, perform cleanup for calling iteration
        if (lastIterationFlag) {
            delete reads2Split;
            if (oldCurToRemove->edgesIn.empty() && oldCurToRemove->edgesOut.empty())
                removeNode(oldCurToRemove);
        }
    }
}

ConsensusGraph::~ConsensusGraph() {
    // std::cerr << "Removing" << std::endl;

    std::vector<Node *> startingNodes;
    for (auto it : readsInGraph)
        startingNodes.push_back(it.second.start);
    removeConnectedNodes(startingNodes);
    if (numEdges || numNodes)
        std::cerr << "graph " << this << " " << std::to_string(numEdges)
                  << " edges " << std::to_string(numNodes) << " nodes left"
                  << std::endl;

    assert(numEdges == 0);
    assert(numNodes == 0);
}

Node *ConsensusGraph::createNode(char base) {
    Node *n = new Node(base);
    numNodes++;
    return n;
}

Edge *ConsensusGraph::createEdge(Node *source, Node *sink, read_t read) {
    Edge *e = new Edge(source, sink, read);
    source->edgesOut.push_back(std::make_pair(sink, e));
    sink->edgesIn.push_back(e);
    numEdges++;
    return e;
}

void ConsensusGraph::removeReadsFromEdge(Edge *e,
                                         std::vector<read_t> const &reads) {
    std::vector<read_t> updatedReadsInOldEdge;
    std::set_difference(
        e->reads.begin(), e->reads.end(), reads.begin(), reads.end(),
        std::inserter(updatedReadsInOldEdge, updatedReadsInOldEdge.begin()));
    e->reads.swap(updatedReadsInOldEdge);
    e->count = e->reads.size();
    if (e->count == 0) {
        removeEdge(e);
    }
}

void ConsensusGraph::removeEdge(Edge *e,
                                bool dontRemoveFromSource /* = false */,
                                bool dontRemoveFromSink /* = false */) {
    if (!dontRemoveFromSource)
        e->source->edgesOut.erase(
            std::find_if(e->source->edgesOut.begin(), e->source->edgesOut.end(),
                         [&](const std::pair<Node *, Edge *> &p) {
                             return (p.first == e->sink);
                         }));

    if (!dontRemoveFromSink)
        e->sink->edgesIn.erase(
            std::find(e->sink->edgesIn.begin(), e->sink->edgesIn.end(), e));
    delete e;
    numEdges--;
}

void ConsensusGraph::removeNode(Node *n) {
    {
        auto edgeIt = n->edgesIn.begin();
        auto end = n->edgesIn.end();
        while (edgeIt != end)
            removeEdge(*(edgeIt++), false, true);
        // don't remove from sink node which is this node!
        // Avoids iterator invalidation and speeds things up
    }
    {
        auto edgeIt = n->edgesOut.begin();
        auto end = n->edgesOut.end();
        while (edgeIt != end)
            removeEdge(edgeIt++->second, true, false);
        // don't remove from source node which is this node!
        // Avoids iterator invalidation and speeds things up
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

void ConsensusGraph::removeConnectedNodes(std::vector<Node *> nodes) {
#ifdef CHECKS
    {
        // This code segment makes sure that the graph is connected and
        // traverseAndCall is working as desired
        size_t count = 0;
        for (Node *n : nodes)
            traverseAndCall(n, false, [&count](Node *) { ++count; });
        assert(count == numNodes);
        count = 0;
        for (Node *n : nodes)
            traverseAndCall(n, true, [&count](Node *) { ++count; });
        assert(count == numNodes);
    }
#endif

    /** Nodes that have no edges in **/
    std::stack<Node *> leafNodes;
    for (Node *n : nodes)
        traverseAndCall(n, false, [&leafNodes](Node *node) {
            if (node->edgesOut.empty())
                leafNodes.push(node);
        });
#ifdef DEBUG
    // std::cerr << "Num of leaf Nodes " << leafNodes.size() << std::endl;
#endif
    while (!leafNodes.empty()) {
        Node *curNode = leafNodes.top();
        leafNodes.pop();
        for (auto edgeIt : curNode->edgesIn) {
            Node *source = edgeIt->source;
            if (source->edgesOut.size() == 1) {
                leafNodes.push(source);
            }
        }
        removeNode(curNode);
    }
}

void ConsensusGraph::printStatus() {
    std::cout << readsInGraph.size() << " reads in graph " << this << ", "
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
    const std::string complementFileName = tempDir + filename + ".complement";
    std::ofstream posFile, editTypeFile, editBaseFile, idFile, complementFile;
    posFile.open(posFileName);
    editTypeFile.open(editTypeFileName);
    editBaseFile.open(editBaseFileName);
    idFile.open(idFileName);
    complementFile.open(complementFileName);
    read_t pasId = 0;
    for (auto it : readsInGraph) {
        {
            idFile << it.first - pasId << ':';
            complementFile << (it.second.reverseComplement ? 'c' : 'n') << ':';
            pasId = it.first;
        }
        totalEditDis +=
            writeRead(posFile, editTypeFile, editBaseFile, it.second, it.first);
    }
    idFile << '\n';
    complementFile << '\n';
    posFile.close();
    editTypeFile.close();
    editBaseFile.close();
    idFile.close();
    complementFile.close();
    std::cout << "AvgEditDis " << totalEditDis / (double)readsInGraph.size()
              << std::endl;
    printStatus();
}

read_t ConsensusGraph::getNumReads() { return readsInGraph.size(); }

size_t ConsensusGraph::read2EditScript(ConsensusGraph::Read &r, read_t id,
                                       std::vector<Edit> &editScript,
                                       size_t &pos) {
    editScript.clear();
    // First we store the initial position
    Node *curNode = r.start;
    assert(curNode);
    bool intersectWithMainPath = true;
    while (!curNode->onMainPath) {
        curNode = curNode->getNextNodeInRead(id);
        if (!curNode) {
            intersectWithMainPath = false;
            break;
        }
    }

    // When the read has no intersection with mainPath, we just store it as a
    // bunch of inserts
    if (!intersectWithMainPath) {
        pos = 0;
        size_t editDis = 0;
        Node *curNode = r.start;
        do {
            editScript.push_back(Edit(INSERT, curNode->base));
            ++editDis;
        } while ((curNode = curNode->getNextNodeInRead(id)));

        return editDis;
    }

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
    } while ((curNode = curNode->getNextNodeInRead(id)));

    dealWithUnchanged();

    return editDis;
}

size_t ConsensusGraph::writeRead(std::ofstream &posFile,
                                 std::ofstream &editTypeFile,
                                 std::ofstream &editBaseFile, Read &r,
                                 read_t id) {

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

size_t ConsensusGraph::writeRead(std::ofstream &f, Read &r, read_t id) {

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

ConsensusGraph::Read::Read(long pos, Node *start, size_t len,
                           bool reverseComplement)
    : pos(pos), start(start), len(len), reverseComplement(reverseComplement) {}

bool ConsensusGraph::checkNoCycle() {
    size_t countNodes = 0;
    size_t countEdges = 0;
    // In this function cumulativeWeight == 1 means is parent of current Node,
    // and cumulativeWeight == 0 means otherwise
    std::vector<Node *> sourceNodes;
    for (auto it : readsInGraph)
        traverseAndCall(it.second.start, false,
                        [&sourceNodes, &countNodes, &countEdges](Node *node) {
                            ++countNodes;
                            countEdges += node->edgesIn.size();
                            node->cumulativeWeight = 0;
                            if (node->edgesIn.empty())
                                sourceNodes.push_back(node);
                        });
    assert(countNodes == numNodes);
    assert(countEdges == numEdges);
    assert(!sourceNodes.empty());
    // Here all the .hasReached has been set to true
    bool status = true;
    for (Node *node : sourceNodes) {
        std::stack<Node *> nodes2Visit;
        nodes2Visit.push(node);
        node->hasReached = !status;
        while (!nodes2Visit.empty()) {
            Node *currentNode = nodes2Visit.top();
            // cumulativeWeight == 1 means parent of current Node (including
            // itself)
            currentNode->cumulativeWeight = 1;
            bool hasUnvisitedChild = false;
            for (auto it : currentNode->edgesOut) {
                // A back edge
                if (it.first->cumulativeWeight == 1)
                    return false;
                // == status means has not visited
                if (it.first->hasReached == status) {
                    nodes2Visit.push(it.first);
                    it.first->hasReached = !status;
                    hasUnvisitedChild = true;
                    continue;
                }
            }
            if (hasUnvisitedChild)
                continue;
            // All children has been visited
            currentNode->cumulativeWeight = 0;
            nodes2Visit.pop();
        }
    }
    return true;
}

template <typename Functor>
void ConsensusGraph::traverseAndCall(Node *n, bool status, Functor f) {
    std::deque<Node *> unfinishedNodes;
    unfinishedNodes.push_back(n);
    while (!unfinishedNodes.empty()) {
        Node *currentNode = unfinishedNodes.front();
        // We don't do anything if it has already been set to !status
        if (currentNode->hasReached != status) {
            unfinishedNodes.pop_front();
            continue;
        }

        currentNode->hasReached = !status;
        unfinishedNodes.pop_front();
        f(currentNode);

        for (auto edgeIt : currentNode->edgesOut) {
            Node *n = edgeIt.second->sink;
            if (n->hasReached == status) {
                unfinishedNodes.push_front(n);
            }
        }

        for (auto edgeIt : currentNode->edgesIn) {
            if (edgeIt->source->hasReached == status)
                unfinishedNodes.push_back(edgeIt->source);
        }
    }
}
