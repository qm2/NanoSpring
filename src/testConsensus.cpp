#include "Contig.h"
#include "NanoporeReads.h"
#include "ReadAligner.h"
#include "Consensus.h"
#include "myers.h"
#include <gperftools/profiler.h>
#include <gperftools/heap-profiler.h>
#include <iostream>
#include <ctime>
#include <chrono>

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
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "Calculating MinHashes took " << duration.count() << " milliseconds" << std::endl;
        }
        ContigGenerator cG(&rA, nR, &rF);
        {
            auto start = std::chrono::high_resolution_clock::now();
            cG.generateContigs();
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "Generating contigs took " << duration.count() << " milliseconds" << std::endl;
        }
        std::cout << cG << std::endl;

        for (Contig *c : cG.contigs) {
            std::set<std::pair<long, read_t>> &readsInContig = c->reads;
            auto currentRead = readsInContig.begin();
            ConsensusGraph consensusGraph(new LocalMyersRollBack(100, 200, 3200));
            consensusGraph.initialize(*cG.nR.readData[currentRead->second],
                                      currentRead->second, currentRead->first);
            consensusGraph.calculateMainPath();
            auto end = readsInContig.end();
            size_t count = 0;
            for (currentRead++; currentRead != end; currentRead++) {
                consensusGraph.addRead(*cG.nR.readData[currentRead->second],
                                       currentRead->second, currentRead->first);
                consensusGraph.calculateMainPath();
                consensusGraph.printStatus();
                if (count % 100 == 0)
                    std::cout << "Added read " << count << std::endl;
                count++;
            }
            consensusGraph.calculateMainPath();
            consensusGraph.printStatus();
        }
    }


    ProfilerStop();
}
