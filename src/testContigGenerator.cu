#include "../include/testContigGenerator.cuh"

int main(int argc, char **argv) {
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
//        MergeSortReadAligner rA(8, 1);
        NanoporeReads nR(argv[1], k, n);
        MinHashReadFilter rF(overlapSketchThreshold, nR);
        nR.calculateMinHashSketches();
        ContigGenerator cG(&rA, nR, &rF);
        cG.generateContigs();
        std::cout << cG << std::endl;
    }
    ProfilerStop();
}