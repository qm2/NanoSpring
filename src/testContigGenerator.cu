#include "../include/testContigGenerator.cuh"
#include "../include/Contig.cuh"
#include <gperftools/profiler.h>
#include <iostream>
#include <ctime>
#include <chrono>

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
    }
    ProfilerStop();
}