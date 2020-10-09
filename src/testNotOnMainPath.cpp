#define private public
#define protected public

#include "ConsensusGraph.h"
#include "LocalMyersRollBack.h"
#include "LocalMyersRollBack_impl.h"
#include <iostream>
#include <random>

/**
 * @brief Produces a random read of length len and stores in inserter
 *
 * @tparam Inserter An insert iterator to char
 * @param len
 * @param inserter
 */
template <typename Inserter> void getRandomRead(size_t len, Inserter inserter) {
    const std::string BASES = "ATCG";
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<> baseDis(0, BASES.length() - 1);
    for (size_t i = 0; i < len; ++i)
        *(inserter++) = BASES[baseDis(generator)];
}

/**
 * @brief Takes in three reads, and creates a ConsensusGraph of the three reads.
 * Assumes all reads have length >= 4
 *
 * @param cG
 * @param reads
 */
void createGraphFromReads(ConsensusGraph *cG, const std::string reads[3]) {
    // Add first read and create mainPath
    size_t readId = 0;
    size_t len = reads[readId].size();
    std::cout << len << std::endl;
    assert(len >= 4);
    Node *nodeBeforeMainPath = cG->createNode('Z');
    cG->startingNode = nodeBeforeMainPath;
    Node *lastNode = nodeBeforeMainPath;
    Node *currentNode = cG->createNode(reads[readId][0]);
    cG->createEdge(lastNode, currentNode, readId);
    currentNode->onMainPath = true;
    cG->mainPath.path.push_back(currentNode->base);
    lastNode = currentNode;
    for (size_t i = 1; i < len; ++i) {
        Node *currentNode = cG->createNode(reads[readId][i]);
        Edge *e = cG->createEdge(lastNode, currentNode, readId);
        currentNode->onMainPath = true;
        cG->mainPath.edges.push_back(e);
        cG->mainPath.path.push_back(currentNode->base);
        lastNode = currentNode;
    }
    Node *nodeAfterMainPath = cG->createNode('Z');
    cG->createEdge(lastNode, nodeAfterMainPath, readId);
    cG->readsInGraph.insert(std::make_pair(
        readId, ConsensusGraph::Read(0, nodeBeforeMainPath, len)));

    ++readId;
    len = reads[readId].size();
    assert(len > 4);
    lastNode = cG->createNode(reads[readId][0]);
    cG->readsInGraph.insert(
        std::make_pair(readId, ConsensusGraph::Read(0, lastNode, len)));
    for (size_t i = 1; i < len; ++i) {
        currentNode = cG->createNode(reads[readId][i]);
        cG->createEdge(lastNode, currentNode, readId);
        lastNode = currentNode;
    }
    cG->createEdge(lastNode, nodeBeforeMainPath, readId);

    ++readId;
    len = reads[readId].size();
    assert(len > 4);
    lastNode = cG->createNode(reads[readId][0]);
    cG->createEdge(nodeAfterMainPath, lastNode, readId);
    cG->readsInGraph.insert(std::make_pair(
        readId, ConsensusGraph::Read(0, nodeAfterMainPath, len)));
    for (size_t i = 1; i < len; ++i) {
        currentNode = cG->createNode(reads[readId][i]);
        cG->createEdge(lastNode, currentNode, readId);
        lastNode = currentNode;
    }
}

int main(int argc, char **argv) {
    LocalMyersRollBack<ConsensusGraph::RAItA, ConsensusGraph::RAItB>
        stringAligner(100, 200, 3200);
    ConsensusGraph cG(&stringAligner);

    std::string reads[3];
    size_t len = 10;
    for (auto &r : reads) {
        getRandomRead(len, std::inserter(r, r.end()));
        std::cout << r << std::endl;
    }

    createGraphFromReads(&cG, reads);

    cG.writeMainPath("Contig0");
    cG.writeReads("Contig0");

    return 0;
}