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
    }
    std::cout << "there" << std::endl;

    return 0;
}
