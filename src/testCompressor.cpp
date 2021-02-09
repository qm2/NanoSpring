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
#include <chrono> 
using namespace std::chrono; 
using namespace std;
// #include <gperftools/heap-profiler.h> 


int main(int argc, char **argv) {
    auto start = high_resolution_clock::now(); 
    omp_set_nested(1);
    // omp_set_num_threads(1);

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
        //modifed the input parameters for the minimap2 version
        size_t k, n, overlapSketchThreshold, m_k, m_w, hashBits;
        //original parameter for aligner
        double maxErrorRate;
        size_t editSlack;
        std::cout << "k n overlapSketchThreshold\nminimap:k minimap:w minimap:hashBits"
                  << std::endl;
        std::cin >> k >> n >> overlapSketchThreshold >>  m_k >>
             m_w >> hashBits;
        if (k == 0)
            return 0;
        std::cout << "k n overlapSketchThreshold minimap:k minimap:w minimap:hashBits\n"
                  << k << " " << n << " " << overlapSketchThreshold << " "
                  << m_k << " " << m_w << " " << hashBits << std::endl;
        MergeSortReadAligner rA(21, 10);
        LocalMyersRollBack<ConsensusGraph::RAItA, ConsensusGraph::RAItB>
            localMyersRollBackAligner(100, 200, editSlack, maxErrorRate);
        ConsensusGraph::StringAligner_t *aligner = &localMyersRollBackAligner;
        Compressor compressor;
        compressor.k = k;
        compressor.n = n;
        compressor.overlapSketchThreshold = overlapSketchThreshold;
        compressor.m_k = m_k;
        compressor.m_w = m_w;
        compressor.hashBits = hashBits;        
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

        auto stop = high_resolution_clock::now(); 
        auto duration = duration_cast<seconds>(stop - start); 
        cout << "Time taken by function: "<< duration.count() << " seconds" << endl; 
    }
    // ProfilerStop();
}
