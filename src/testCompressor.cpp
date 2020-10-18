#include "Compressor.h"
#include "LocalMyersRollBack.h"
#include "LocalMyersRollBack_impl.h"
#include "ReadData.h"
#include "ReadFilter.h"
#include "StringAligner.h"
#include <ctime>
#include <gperftools/profiler.h>
#include <iostream>
#include <omp.h>
// #include <gperftools/heap-profiler.h>

int main(int argc, char **argv) {
    omp_set_nested(1);
    // omp_set_num_threads(1);
    // if (!fork())
    //     return -1;
    std::srand(unsigned(std::time(0)));
    // ProfilerStart("testCompressor.prof");
    if (argc < 2) {
        std::cout << "Usage ./testCompressor filename" << std::endl;
        return 1;
    }
    {
        size_t k, n, overlapSketchThreshold;
        std::cout << "k n overlapSketchThreshold" << std::endl;
        std::cin >> k >> n >> overlapSketchThreshold;
        if (k == 0)
            return 0;
        MergeSortReadAligner rA(21, 10);
        LocalMyersRollBack<ConsensusGraph::RAItA, ConsensusGraph::RAItB>
            localMyersRollBackAligner(100, 200, 100, 0.21);
        ConsensusGraph::StringAligner_t *aligner = &localMyersRollBackAligner;
        Compressor compressor;
        compressor.k = k;
        compressor.n = n;
        compressor.overlapSketchThreshold = overlapSketchThreshold;
        compressor.rA = &rA;
        compressor.aligner = aligner;
        compressor.compress(argv[1]);
    }
    // ProfilerStop();
}
