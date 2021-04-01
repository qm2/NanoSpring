#include "Compressor.h"
#include "ReadData.h"
#include "ReadFilter.h"
#include "StringAligner.h"
#include <ctime>
#include <iostream>
#include <omp.h>
#include <chrono> 
using namespace std::chrono; 
using namespace std;

int main(int argc, char **argv) {
    auto start = high_resolution_clock::now(); 
    omp_set_nested(1);
    std::srand(unsigned(std::time(0)));
    if (argc < 2) {
        std::cout << "Usage ./testCompressor filename [numThr]" << std::endl;
        return 1;
    }
    int numThr = omp_get_max_threads();
    if (argc == 3)
       numThr = atoi(argv[2]); 
    {
        //modifed the input parameters for the minimap2 version
        size_t k, n, overlapSketchThreshold, m_k, m_w, max_chain_iter;
        //original parameter for aligner
        double maxErrorRate;
        size_t editSlack;
        std::cout << "k n overlapSketchThreshold\nminimap:k minimap:w minimap:max_chain_iter"
                  << std::endl;
        std::cin >> k >> n >> overlapSketchThreshold >>  m_k >>
             m_w >> max_chain_iter;
        if (k == 0)
            return 0;
        std::cout << "k n overlapSketchThreshold minimap:k minimap:w minimap:max_chain_iter\n"
                  << k << " " << n << " " << overlapSketchThreshold << " "
                  << m_k << " " << m_w << " " << max_chain_iter << std::endl;
        MergeSortReadAligner rA(21, 10);
        Compressor compressor;
        compressor.k = k;
        compressor.n = n;
        compressor.overlapSketchThreshold = overlapSketchThreshold;
        compressor.m_k = m_k;
        compressor.m_w = m_w;
        compressor.max_chain_iter = max_chain_iter;        
        compressor.rA = &rA;
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
}
