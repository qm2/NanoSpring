#include "Contig.h"
#include "ReadAligner.h"
#include "ReadData.h"
#include "ReadFilter.h"
#include <chrono>
#include <ctime>
#include <gperftools/profiler.h>
#include <iostream>

int main(int argc, char **argv) {
    std::srand(unsigned(std::time(0)));
    ProfilerStart("testContig.prof");
    if (argc < 2) {
        std::cout << "Usage ./testContig filename" << std::endl;
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
        ReadData rD;
        rD.loadFromFile(argv[1]);

        // NanoporeReads nR(argv[1], k, n);
        MinHashReadFilter rF;
        rF.k = k;
        rF.n = n;
        rF.overlapSketchThreshold = overlapSketchThreshold;
        {
            auto start = std::chrono::high_resolution_clock::now();
            rF.initialize(rD);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start);
            std::cout << "Calculating MinHashes took " << duration.count()
                      << " milliseconds" << std::endl;
        }
        ContigGenerator cG(&rA, &rD, &rF);
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
    }
    ProfilerStop();
}