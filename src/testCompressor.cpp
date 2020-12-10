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


    // if (!fork())
    //     return -1;
    std::srand(unsigned(std::time(0)));
    // ProfilerStart("testCompressor.prof");
    if (argc < 2) {
        std::cout << "Usage ./testCompressor filename [numThr]" << std::endl;
        return 1;
    }
    int numThr = omp_get_max_threads();
    if (argc == 3)
       numThr = atoi(argv[2]); 
    {
        size_t k, n, overlapSketchThreshold, editSlack;
        double maxErrorRate;
        std::cout << "k n overlapSketchThreshold\neditSlack maxErrorRate"
                  << std::endl;
        std::cin >> k >> n >> overlapSketchThreshold >> editSlack >>
            maxErrorRate;
        if (k == 0)
            return 0;
        std::cout << "k n overlapSketchThreshold editSlack maxErrorRate\n"
                  << k << " " << n << " " << overlapSketchThreshold << " "
                  << editSlack << " " << maxErrorRate << std::endl;
        MergeSortReadAligner rA(21, 10);
        LocalMyersRollBack<ConsensusGraph::RAItA, ConsensusGraph::RAItB>
            localMyersRollBackAligner(100, 200, editSlack, maxErrorRate);
        ConsensusGraph::StringAligner_t *aligner = &localMyersRollBackAligner;
        Compressor compressor;
        compressor.k = k;
        compressor.n = n;
        compressor.overlapSketchThreshold = overlapSketchThreshold;
        compressor.rA = &rA;
        compressor.aligner = aligner;
        const std::string filename(argv[1]);
        const std::string extension =
            filename.substr(filename.find_last_of('.') + 1);
        if (!extension.compare("reads"))
            compressor.filetype = ReadData::Filetype::READ;
        else if (!extension.compare("fastq"))
            compressor.filetype = ReadData::Filetype::FASTQ;
        compressor.compress(argv[1], numThr);
    }
    // ProfilerStop();
}
