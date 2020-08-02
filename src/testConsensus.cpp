#include "Consensus.h"
#include "Contig.h"
#include "NanoporeReads.h"
#include "ReadAligner.h"
#include "myers.h"
#include <chrono>
#include <ctime>
#include <fstream>
#include <gperftools/heap-profiler.h>
#include <gperftools/profiler.h>
#include <iostream>

int main(int argc, char **argv) {
    std::srand(unsigned(std::time(0)));
    ProfilerStart("testConsensus.prof");
    if (argc < 2) {
        std::cout << "Usage ./testConsensus filename" << std::endl;
        return 1;
    }
    {
        size_t k, n, overlapSketchThreshold;
        std::cout << "k n overlapSketchThreshold" << std::endl;
        std::cin >> k >> n >> overlapSketchThreshold;
        if (k == 0)
            return 0;
        MergeSortReadAligner rA(21, 10);
        //        MergeSortReadAligner rA(10, 1);
        NanoporeReads nR(argv[1], k, n);
        MinHashReadFilter rF(overlapSketchThreshold, nR);
        {
            auto start = std::chrono::high_resolution_clock::now();
            nR.calculateMinHashSketches();
            auto end = std::chrono::high_resolution_clock::now();
            auto duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start);
            std::cout << "Calculating MinHashes took " << duration.count()
                      << " milliseconds" << std::endl;
        }
        ContigGenerator cG(&rA, nR, &rF);
        {
            auto start = std::chrono::high_resolution_clock::now();
            cG.generateContigs();
            auto end = std::chrono::high_resolution_clock::now();
            auto duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start);
            std::cout << "Generating contigs took " << duration.count()
                      << " milliseconds" << std::endl;
        }
        std::cout << cG << std::endl;

        size_t i = 0;
        for (Contig *c : cG.contigs) {
            std::set<std::pair<long, read_t>> &readsInContig = c->reads;
            auto currentRead = readsInContig.begin();
            ConsensusGraph consensusGraph(
                new LocalMyersRollBack(100, 200, 3200));

            consensusGraph.addReads(readsInContig, cG.nR.readData);

            consensusGraph.calculateMainPathGreedy();
            consensusGraph.printStatus();
            std::ofstream f;
            f.open("testContig" + std::to_string(i) + ".genome");
            consensusGraph.writeMainPath(f);
            f.close();

            f.open("testContig" + std::to_string(i) + ".reads");
            consensusGraph.writeReads(f);
            f.close();
            ++i;
        }
    }

    ProfilerStop();
}
