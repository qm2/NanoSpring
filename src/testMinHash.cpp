#include "ReadData.h"
#include "ReadFilter.h"
#include <gperftools/profiler.h>
#include <iostream>
#include <omp.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "Usage: testMinHash fileName k n" << std::endl;
        return 1;
    }
    omp_set_nested(1);
    {
        ProfilerStart("testMinHash.prof");
        //        nanoporeReads.printHashes();
        ReadData rD;
        rD.loadFromFile(argv[1]);
        std::cout << rD.getNumReads() << " reads\n";
        MinHashReadFilter *minHashReadFilter = new MinHashReadFilter();
        minHashReadFilter->k = std::stoi(argv[2]);
        minHashReadFilter->n = std::stoi(argv[3]);
        ReadFilter *rF = minHashReadFilter;
        rF->initialize(rD);
        std::cout << rF->getFilterStats(1000, 4);
        delete rF;
        ProfilerStop();
        //        unsigned int overlapBaseThs[] = {10, 20, 50, 100, 500, 1000};
        //        unsigned int overlapBaseThs[] = {500, 1000, 2000, 5000, 8000};
        //        unsigned int overlapSketchThs[] = {1, 2, 3, 4, 5, 6, 7};
        //        for (unsigned int overlapBaseTh: overlapBaseThs) {
        //            for (unsigned int overlapSketchTh : overlapSketchThs) {
        //                std::cout <<
        //                nanoporeReads.getFilterStats(overlapBaseTh,
        //                overlapSketchTh) << std::endl;
        //
        //            }
        //        }
        std::cout << "calculating filter stats" << std::endl;
        // while (true) {
        //     unsigned int overlapBaseTh, overlapSketchTh;
        //     std::cin >> overlapBaseTh >> overlapSketchTh;
        // std::cout << nanoporeReads.getFilterStats(overlapBaseTh,
        //                                           overlapSketchTh)
        //           << std::endl;
        // }
        //                std::cout << nanoporeReads.getFilterStats(1000, 1) <<
        //                std::endl;
    }

    return 0;
}
