#include "ConsensusGraph.h"
#include "bsc_helper.h"
#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <stack>
#include <zlib.h>
#include "minimap.h"

Edge::Edge(Node *source, Node *sink, read_t read) : source(source), sink(sink) {
    count = 1;
    // insert into sorted vector
    reads.push_back(read);
}

Edge::Edge(Node *source, Node *sink, std::vector<read_t> &reads)
    : source(source), sink(sink), reads(reads) {
    count = reads.size();
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

ConsensusGraphWriter::ConsensusGraphWriter(const std::string &filePrefix) {
    const std::string posFileName = filePrefix + ".pos";
    const std::string editTypeFileName = filePrefix + ".type";
    const std::string editBaseFileName = filePrefix + ".base";
    const std::string idFileName = filePrefix + ".id";
    const std::string complementFileName = filePrefix + ".complement";
    const std::string genomeFileName = filePrefix + ".genome";
    posFile.open(posFileName);
    editTypeFile.open(editTypeFileName);
    editBaseFile.open(editBaseFileName);
    idFile.open(idFileName);
    complementFile.open(complementFileName);
    genomeFile.open(genomeFileName);
}

void ConsensusGraph::initialize(const std::string &seed, read_t readId,
                                long pos) {
    // pos is zero here                                
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
    startPos = pos; // = 0
    endPos = pos + 1; 
    // endPos is set to pos+1=1 here because we only inserted one base to mainPath
    // Rest will be inserted later when calculateMainPathGreedy is called
}

bool ConsensusGraph::addRead(const std::string &s, std::vector<Edit> &editScript,
            ssize_t &beginOffset, ssize_t &endOffset, size_t m_k, size_t m_w, size_t hashBits) {
    // General comments: 
    // 1. The whole idea behind startPos, endPos and Read.pos in ConsensusGraph
    //    is to provide a reference point for the main path when looking for the next 
    //    read in addRelatedReads in Consensus.cpp. This is needed because we extend the 
    //    mainPath to the left in some cases, and so curPos in Consensus::addRelatedReads
    //    doesn't make sense unless we offset it by startPos which essentially tells us
    //    where the current start position of mainPath is wrt to the start of the first read
    //    in the contig (so startPos is always <=0). endPos is the end of the mainPath wrt
    //    the start of the first read (e.g., it is equal to the length of the first read at 
    //    the very start of the contig when there is only one read). It generally increases 
    //    as we add reads to the right. endPos seems less important for us. Finally, for each 
    //    read in the graph, there is a pos variable (0 for first read) which tells us roughly 
    //    where on the mainPath the read lies (it is actually decided in Consensus::addRelatedReads
    //    itself after the sort-merge procedure). The only use for this pos variable seems 
    //    to be for computation for endPos and startPos in calculateMainPathGreedy after the 
    //    new read is added in. The pos variable is also defined wrt the start of the first read.


    auto &originalString = mainPath.path;
    /// We use either the head or tail of mainPath as a reference to obtain an
    /// offsetGuess

    // TODO: can we remove the Read.pos, startPos, endPos variables in the minimap case 
    //       to help simplify things significantly?

    size_t editDis;

    std::string originalStringCopy(originalString.begin(),originalString.end());
    const char* Abegin = originalStringCopy.c_str();
    const char* Bbegin = s.c_str();
    //we comment out the original aligner for the minimap
    // bool success = aligner->align(Abegin, Aend, Bbegin, Bend, offsetGuess,
    //                               beginOffset, endOffset, editScript, editDis);
    int hits;
    bool success = true;
    //initialize the local buffer
    mm_tbuf_t *b = mm_tbuf_init();
    //initialize the mapopt and iopt
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    //0 correpons to map-ont
    mm_set_opt(0, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR;
    mopt.flag |= MM_F_FOR_ONLY; 
    // only forward alignment, no reverse complement (which is handled elsewhere)
    //call the mm_idx_str to return the index for the reference read    
    // std::cout<<"k:"<<iopt.k<<"w:"<<iopt.w<<std::endl;
    // std::cout<<"flag:"<<iopt.flag<<"bits:"<<iopt.bucket_bits<<hits<<std::endl;      
    //the defalut parameters are: 15 10 false 14
    mm_idx_t * idx = mm_idx_str(m_k, m_w, false, hashBits, 1, &Abegin, NULL);
    mm_mapopt_update(&mopt, idx);
    //use the index to align with the current read
    //we only want forward matches: use rev in mm_reg1_t
    mm_reg1_t* reg = mm_map(idx, s.length(), s.c_str(), &hits, b, &mopt, NULL);
    editScript.clear();
    //experiment how many hits there are
    if(hits > 0) { 
        //if we have multiple hits, just stick with first hit
        mm_reg1_t *r = &reg[0];
        assert(r->p); 
        //qpos is the current position on the query read
        //rpos is the current position on the reference read
        int qpos = r->qs;
        int rpos = r->rs;
        unsigned int j, k;
        unsigned int count_same;
        int i;

        // See comments in ConsensusGraph::updateGraph for its high-level functioning

        // Based on my understanding the correct approach would be something like this:
        // - beginOffset - this represents where the alignment starts on the reference 
        //   (can be +ve or -ve). So if r->rs is +ve then beginOffset is simply r->rs. 
        //   On the other hand if r->rs = 0 (so there is part of read to left), beginOffset
        //   is -r->qs (that's how much read is shifted to left wrt reference).
        // - Now if r->rs is +ve, then we need to add the soft clipped bases as insertions in
        //   editScript. When r->rs = 0 (and hence beginOffset is -ve), the first 
        //   |beginOffset| bases in the read (which are the soft clipped bases) will be added
        //   to the graph in updateGraph directly. So we don't add them to editScript in this 
        //   case.
        // - endOffset - this represents where the alignment ends on the reference wrt end of 
        //   reference. So if r->re < originalString.size() (i.e., the alignment ends 
        //   before the last base in reference), endOffset is -ve and equal to 
        //   (r->re-originalString.size()). If r->re == originalString.size(), the read alignment 
        //   potentially extends beyond the end of reference, and endOffset is +ve and equal to 
        //   (s.length()-r->qe). 
        // - If endOffset is -ve, we need to add the soft-clipped bases at end to the editScript 
        //   as insertions. If endOffset is +ve, these bases will be automatically added in 
        //   the graph in updateGraph function, so we do not add them to editScript.


        //update the beginOffset by checking if r->rs is positive or not
        if(r->rs >0){
            beginOffset = r->rs;  
            //add the soft clipped bases as insertions
            for (i = 0; i < r->qs; ++i){
                editScript.push_back(Edit(INSERT, s[i]));             
            }       
        }
        else if(r->rs == 0){
            beginOffset = -r->qs;               
        }
        else{
            throw std::runtime_error("Encountered invalid reference start");
        }

        for (j = 0; j < r->p->n_cigar; ++j){ // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
            //std::cout<< (r->p->cigar[j]>>4) << "MIDNSH"[r->p->cigar[j]&0xf];
            switch("MIDNSH"[r->p->cigar[j]&0xf]) {
                case 'M':
                    //fix: handle the substitute and same differently
                    count_same = 0;
                    for(k=0; k<r->p->cigar[j]>>4;k++){
                        if(s[qpos] == Abegin[rpos]){
                            count_same++;
                        }
                        else{
                            if (count_same > 0)
                                editScript.push_back(Edit(SAME, count_same));
                            count_same = 0; // reset count_same
                            //add the substitute as one insert and one delete
                            editScript.push_back(Edit(DELETE,Abegin[rpos]));
                            editScript.push_back(Edit(INSERT,s[qpos]));                             
                        }
                        qpos++;
                        rpos++;
                    }
                    //check if we need to push the "same" again
                    if(count_same!=0){
                        editScript.push_back(Edit(SAME, count_same));                       
                    }
                    break; 
                case 'I':
                    for(k=0; k<r->p->cigar[j]>>4;k++){
                        editScript.push_back(Edit(INSERT, s[qpos])); 
                        qpos++; 
                    }                   
                    break; 
                case 'D':
                    for(k=0; k<r->p->cigar[j]>>4;k++){
                        editScript.push_back(Edit(DELETE, Abegin[rpos])); 
                        rpos++; 
                    }                       
                    break; 
                default: 
                    throw std::runtime_error("Encountered invalid CIGAR symbol!");
            }           
        }

        //update the endOffset by checking if r->re is positive or not
        if(r->re < originalString.size()){
            endOffset = (r->re-originalString.size());
            //add the soft clipped bases as insertions
            for (i = r->qe; i < s.length(); ++i){
                editScript.push_back(Edit(INSERT, s[i]));             
            }       
        }
        else if(r->re == originalString.size()){
            endOffset = (s.length()-r->qe);             
        }
        else{
            throw std::runtime_error("Encountered invalid reference end");
        }
        //calculate the edit distance
        editDis = r->blen - r->mlen + r->p->n_ambi;
#ifdef LOG
        std::cout<< "editDis: " << editDis<<std::endl;     
#endif
        free(r->p);
        if (hits > 1) {
            // cleanup
            for (int i = 1; i < hits; i++)
                free(reg[i].p);
            }
    }
    else{
        //return fail when there are not hits
        success = false;
    }
    //return the correct beginOffset and endOffset   
    editScript.shrink_to_fit();
    free(reg);  
    mm_tbuf_destroy(b);
    mm_idx_destroy(idx);       
    //number of hits
    //edit distance as threshold
    //length of alignment as fraction
    //DP alignment score
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
}

void ConsensusGraph::updateGraph(const std::string &s,
                                 std::vector<Edit> &editScript,
                                 ssize_t beginOffset, ssize_t endOffset,
                                 read_t readId, long pos,
                                 bool reverseComplement) {

    // How does the updateGraph function work (high level)?
    // 1. If beginOffset is -ve, create |beginOffset| new nodes in the graph 
    //    joined in a chain representing those bases in read 
    //    (the chain will be essentially joined to position 0 of mainPath in graph).
    // 2. Then begin adding nodes & edges to graph according to the editScript 
    //    at the relevant place in the graph as determined by beginOffset 
    //    (essentially starting at position 0 of mainPath when beginOffset is -ve)
    // 3. Finally, if endOffset is +ve, create |endOffset| new nodes in the graph 
    //    joined in a chain representing those bases in read (the chain begins 
    //    in the graph where step 2 ends, which will be roughly the end of mainPath 
    //    in this case). 

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
        // TODO: why take the begin (i.e., first read through this last edge in path)?
        // Is this a relic of the constant read length code?
        Read &endingRead = readsInGraph.at(endingReadId);
        endPos = endingRead.pos + endingRead.len;
    }

    // Extend to the left
    // NOTE: I have changed mainPath.path to vector to speed up alignment,
    // but that means the code below is not efficient (inserting to start
    // of vector). Fortunately, this is currently a very small contributor
    // to the total time. But we might want to fix this if issues crop up later.
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
    std::stack<Edge *> walkAndPruneCallStack; 
    // used for converting recursion to iteration.
    // Define here to avoid overhead of repeated allocations
    while (true) {
        // Then we iterate over all edges pointing to side nodes that have
        // other edges in
        //
        // Work with copy to avoid iterator invalidation.
        auto edgesOutCopy = nodeOnPath->edgesOut;
        for (const auto &edgeIt : edgesOutCopy)
            walkAndPrune(edgeIt.second, walkAndPruneCallStack);

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
            walkAndPrune(edgeIt.second, walkAndPruneCallStack);

        if (edgeOnPath == mainPath.edges.begin())
            break;
        --edgeOnPath;
    }
}

void ConsensusGraph::walkAndPrune(Edge *e, std::stack<Edge *> &callStack) {
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
        // Now put sink->edgesOut into stack
        for (auto it = sink->edgesOut.begin(); it != sink->edgesOut.end();
             ++it) {
            callStack.push(it->second);
            // OLD COMMENT:
            // Here walkAndPrune will only change sink->edgesOut by 1. deleting
            // edge2WorkOn and 2. adding new edges that we don't need to prune.
            // So copying sink->edgesOut should work.
        }
    }
}

void ConsensusGraph::splitPath(Node *newPre, Edge *e,
                               std::vector<read_t> *reads2Split) {
    /**
     * @brief The context for changing recursion to iteration.
     *
     * Every context will be visited exactly twice.
     *
     * When constructed, hasVisited is set to false, oldCur is set to
     * nullptr, and reads2Split is set to the pointer passed down by the last
     * iteration.
     *
     * After the first visit, hasVisited is set to true, reads2Split is set to
     * reads2Split created in this iteration and passed on to future iterations.
     * The new value is created by new and the destructor is responsible for
     * deleting it. And oldCur will be set to e->sink, and the destructor will
     * delete it if it becomes disconnected.
     *
     */
    class SplitPathContext {
    public:
        Node *const newPre;
        Edge *const e;
        std::vector<read_t> *reads2Split;
        bool hasVisited;
        /** we need to store this (e->sink) because e might be deleted and we
         * still need to delete e->sink **/
        Node *oldCur;

        SplitPathContext(ConsensusGraph *cG, Node *newPre, Edge *e,
                         std::vector<read_t> *reads2Split)
            : newPre(newPre), e(e), reads2Split(reads2Split), hasVisited(false),
              oldCur(nullptr), cG(cG) {}

        ~SplitPathContext() {
            delete reads2Split;
            if (oldCur && oldCur->edgesIn.empty() && oldCur->edgesOut.empty())
                cG->removeNode(oldCur);
        }

    private:
        ConsensusGraph *const cG;
    };

    std::stack<SplitPathContext> callStack;

    callStack.emplace(this, newPre, e, reads2Split);

    while (!callStack.empty()) {
        SplitPathContext &currentContext = callStack.top();

        if (currentContext.hasVisited) {
            callStack.pop();
            continue;
        }

        // The reads that need to be split going down this path
        auto readsInPath2Split = new std::vector<read_t>;
        std::set_intersection(
            currentContext.reads2Split->begin(),
            currentContext.reads2Split->end(), currentContext.e->reads.begin(),
            currentContext.e->reads.end(),
            std::inserter(*readsInPath2Split, readsInPath2Split->begin()));

        currentContext.reads2Split = readsInPath2Split;
        currentContext.hasVisited = true;

        if (readsInPath2Split->empty())
            continue;

        Node *oldCur = currentContext.e->sink;
        currentContext.oldCur = oldCur;

        // Remove the reads from the node edge
        removeReadsFromEdge(currentContext.e, *readsInPath2Split);

        // just create a new edge and return if we are arriving at a node in
        // mainPath
        if (oldCur->onMainPath) {
            createEdge(currentContext.newPre, oldCur, *readsInPath2Split);
            continue;
        }

        // Otherwise create a new node
        Node *newCur = createNode(oldCur->base);
        // Add an edge to this new node
        createEdge(currentContext.newPre, newCur, *readsInPath2Split);

        // Now put sink->edgesOut into stack
        for (auto &it : oldCur->edgesOut)
            callStack.emplace(this, newCur, it.second, readsInPath2Split);
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

Edge *ConsensusGraph::createEdge(Node *source, Node *sink,
                                 std::vector<read_t> &reads) {
    Edge *e = new Edge(source, sink, reads);
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
                e->source->edgesOut.erase(std::find_if(
                    e->source->edgesOut.begin(), e->source->edgesOut.end(),
                    [&](const std::pair<Node *, Edge *> &p) {
                        return (p.first == e->sink);
                    }));

            if (!dontRemoveFromSink)
                e->sink->edgesIn.erase(std::find(e->sink->edgesIn.begin(),
                                                 e->sink->edgesIn.end(), e));
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
                // std::cerr << "Num of leaf Nodes " << leafNodes.size() <<
                // std::endl;
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
            std::cout << readsInGraph.size() << " reads in graph " << this
                      << ", " << numNodes << " nodes, and " << numEdges
                      << " edges."
                      << "\n";
            std::cout << "mainPath len " << mainPath.edges.size() + 1 << " "
                      << mainPath.path.size() << " avg weight "
                      << mainPath.getAverageWeight() << " starts at "
                      << startPos << " ends at " << endPos << " len in Contig "
                      << endPos - startPos << "\n";
            double stat = mainPath.edges.size() * mainPath.getAverageWeight() /
                          (readsInGraph.size() * 9999);
            std::cout << stat / (1 - 0.17) << " " << stat / (1 - 0.17 * 1.05)
                      << " " << stat / (1 - 0.17 * 1.1) << " " << std::endl;
        }

        void ConsensusGraph::writeMainPath(ConsensusGraphWriter &cgw) {
            cgw.genomeFile << std::string(mainPath.path.begin(), mainPath.path.end())
              << std::endl;
            cgw.genomeFile << ".\n";
        }

        void ConsensusGraph::writeReads(ConsensusGraphWriter &cgw) {
            // First we write the index of each character into the
            // cumulativeWeight field of the nodes on mainPath
            mainPath.edges.front()->source->cumulativeWeight = 0;
            size_t i = 0;
            for (auto e : mainPath.edges) {
                e->sink->cumulativeWeight = ++i;
            }
            size_t totalEditDis = 0;

            read_t pasId = 0;
            for (auto it : readsInGraph) {
                {
                    cgw.idFile << it.first - pasId << ':';
                    cgw.complementFile << (it.second.reverseComplement ? 'c' : 'n')
                                   << ':';
                    pasId = it.first;
                }
                totalEditDis += writeRead(cgw.posFile, cgw.editTypeFile, cgw.editBaseFile,
                                          it.second, it.first);
            }
            cgw.idFile << '\n';
            cgw.complementFile << '\n';
            cgw.posFile << ".\n";
            cgw.editTypeFile << ".\n";
            cgw.editBaseFile << ".\n";
            cgw.idFile << ".\n";
            cgw.complementFile << ".\n";
#ifdef LOG
            std::cout << "AvgEditDis "
                      << totalEditDis / (double)readsInGraph.size()
                      << std::endl;
            printStatus();
#endif
        }

        read_t ConsensusGraph::getNumReads() { return readsInGraph.size(); }

        size_t ConsensusGraph::read2EditScript(ConsensusGraph::Read &r,
                                               read_t id,
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

            // When the read has no intersection with mainPath, we just store it
            // as a bunch of inserts
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
        
        /*
        // NOTE: Commenting this out since it doesn't seem to be in use.
        // Please remove if not needed. 
        // -Shubham
        void ConsensusGraph::writeReads(std::ofstream &f) {
            // First we write the index of each character into the
            // cumulativeWeight field of the nodes on mainPath
            mainPath.edges.front()->source->cumulativeWeight = 0;
            size_t i = 0;
            for (auto e : mainPath.edges) {
                e->sink->cumulativeWeight = ++i;
            }
            size_t totalEditDis = 0;
            for (auto it : readsInGraph) {
                totalEditDis += writeRead(f, it.second, it.first);
            }
            std::cout << "AvgEditDis "
                      << totalEditDis / (double)readsInGraph.size()
                      << std::endl;
        }
        */

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

        ConsensusGraph::ConsensusGraph(StringAligner_t *aligner)
            : aligner(aligner) {}

        ConsensusGraph::Read::Read(long pos, Node *start, size_t len,
                                   bool reverseComplement)
            : pos(pos), start(start), len(len),
              reverseComplement(reverseComplement) {}

        bool ConsensusGraph::checkNoCycle() {
            size_t countNodes = 0;
            size_t countEdges = 0;
            // In this function cumulativeWeight == 1 means is parent of current
            // Node, and cumulativeWeight == 0 means otherwise
            std::vector<Node *> sourceNodes;
            for (auto it : readsInGraph)
                traverseAndCall(
                    it.second.start, false,
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
                    // cumulativeWeight == 1 means parent of current Node
                    // (including itself)
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
