//
// Created by The MAC PRO on 2020/7/3.
//

#include "../include/testMinHash.cuh"


int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "Usage: testMinHash fileName k n" << std::endl;
        return 1;
    }
    {
//        ProfilerStart("testMinHash.prof");
        NanoporeReads nanoporeReads(argv[1], atol(argv[2]), atol(argv[3]));
        nanoporeReads.calculateMinHashSketches();
//        nanoporeReads.printHashes();
        std::cout << "here" << std::endl;
//        ProfilerStop();
//        unsigned int overlapBaseThs[] = {10, 20, 50, 100, 500, 1000};
//        unsigned int overlapBaseThs[] = {500, 1000, 2000, 5000, 8000};
//        unsigned int overlapSketchThs[] = {1, 2, 3, 4, 5, 6, 7};
//        for (unsigned int overlapBaseTh: overlapBaseThs) {
//            for (unsigned int overlapSketchTh : overlapSketchThs) {
//                std::cout << nanoporeReads.getFilterStats(overlapBaseTh, overlapSketchTh) << std::endl;
//
//            }
//        }
//        while (true) {
//            unsigned int overlapBaseTh, overlapSketchTh;
//            std::cin >> overlapBaseTh >> overlapSketchTh;
//            std::cout << nanoporeReads.getFilterStats(overlapBaseTh, overlapSketchTh) << std::endl;
//
//        }
//                std::cout << nanoporeReads.getFilterStats(1000, 1) << std::endl;

    }

    std::cout << "there" << std::endl;

    return 0;
}
