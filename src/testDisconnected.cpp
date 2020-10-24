#define private public
#define protected public

#include "ConsensusGraph.h"
#include "LocalMyersRollBack.h"
#include "LocalMyersRollBack_impl.h"
#include <iostream>
#include <random>

int main(int argc, char **argv) {
    (void)argc;
    (void)argv;
    LocalMyersRollBack<ConsensusGraph::RAItA, ConsensusGraph::RAItB>
        stringAligner(100, 200, 100, 0.21);
    {
        ConsensusGraph cG(&stringAligner);

        cG.initialize(std::string("AB"), 0, 0);
        Node *A = cG.readsInGraph.at(0).start;
        cG.calculateMainPathGreedy();

        Node *C = cG.createNode('C');
        cG.createEdge(A, C, 1);
        cG.readsInGraph.insert(
            std::make_pair(1, ConsensusGraph::Read(0, A, 2, false)));
        Node *D = cG.createNode('D');
        cG.createEdge(D, C, 2);
        cG.readsInGraph.insert(
            std::make_pair(2, ConsensusGraph::Read(0, D, 2, false)));
        cG.calculateMainPathGreedy();
        cG.checkNoCycle();
    }
    {
        ConsensusGraph cG(&stringAligner);

        cG.initialize(std::string("AB"), 0, 0);
        Node *A = cG.readsInGraph.at(0).start;
        Node *B = A->edgesOut.front().first;
        cG.calculateMainPathGreedy();

        Node *C = cG.createNode('C');
        cG.createEdge(A, C, 1);
        cG.readsInGraph.insert(
            std::make_pair(1, ConsensusGraph::Read(0, A, 2, false)));
        Node *D = cG.createNode('D');
        cG.createEdge(D, C, 2);
        cG.readsInGraph.insert(
            std::make_pair(2, ConsensusGraph::Read(0, D, 2, false)));
        cG.createEdge(B, C, 3);
        cG.readsInGraph.insert(
            std::make_pair(3, ConsensusGraph::Read(0, B, 2, false)));
        cG.checkNoCycle();
        cG.calculateMainPathGreedy();
        cG.checkNoCycle();
    }
    return 0;
}