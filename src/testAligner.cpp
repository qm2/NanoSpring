#include "testAligner.h"
#include "ReadAligner.h"
#include <ctime>
#include <gperftools/profiler.h>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv) {
    srand(time(NULL));
    ProfilerStart("testAligner.prof");
    if (argc < 2) {
        std::cout << "Usage ./testAligner filename" << std::endl;
        return 1;
    }
    TestAligner ta(argv[1]);
    size_t k, kMerNumTh, baseOverlapTh, n;
    while (true) {
        std::cout << "k kMerNumth baseOverlapTh n" << std::endl;
        std::cin >> k >> kMerNumTh >> baseOverlapTh >> n;
        if (k == 0)
            break;
        for (size_t i = 0; i < n; ++i)
            ta.test(k, kMerNumTh, baseOverlapTh);
    }
    ProfilerStop();
}

TestAligner::TestAligner(const char *fileName) {
    rD.loadFromFile(fileName);
    std::cout << "Finished Initializing TestAligner from " << fileName
              << std::endl;
}

void TestAligner::test(const size_t k, const size_t kMerNumTh,
                       const size_t baseOverlapTh) {
    ReadAligner *rA = new MergeSortReadAligner(k, kMerNumTh);
    size_t randomReadIndex = rand() % rD.getNumReads();
    const std::string &randomRead = rD.getRead(randomReadIndex);
    const long randomPos = rD.getReadPos()[randomReadIndex];
    const long th = rD.getRead(0).size() - baseOverlapTh;
    std::cout << "Pos of read is " << randomPos << std::endl;
    size_t numPositives = 0;
    size_t truePositives = 0;
    size_t falsePositives = 0;
    size_t numNegatives = 0;
    size_t falseNegatives = 0;
    size_t trueNegatives = 0;
    double posError = 0;
    size_t numReads = rD.getNumReads();
    size_t readLen = rD.getRead(0).size();
#pragma omp parallel for reduction(+ : numPositives, numNegatives, \
    falsePositives, truePositives, falseNegatives, trueNegatives, posError)
    for (size_t i = 0; i < numReads; ++i) {
        if (i == randomReadIndex)
            continue;
        //        if (i % 1000 == 0)
        //            std::cout << i << std::endl;
        ssize_t relPos;
        if (rA->align(randomRead, rD.getRead(i), relPos)) {
            // std::cout << "Real " << (long) nR.readPos[i] - (long) randomPos
            //<< " Predicted " << relPos << std::endl;
            posError +=
                abs((long)rD.getReadPos()[i] - (long)randomPos - relPos);
            numPositives++;
            if (abs(randomPos - (long)rD.getReadPos()[i]) > th) {
                falsePositives++;
                std::cout << (long)readLen -
                                 abs((long)randomPos - (long)rD.getReadPos()[i])
                          << std::endl;
            } else
                truePositives++;
        } else {
            numNegatives++;
            if (abs(randomPos - (long)rD.getReadPos()[i]) > th)
                trueNegatives++;
            else
                falseNegatives++;
        }
    }
    posError /= numPositives;
    std::cout << "Average position error is " << posError << std::endl;
    const int w = 13;
    std::cout << std::setw(w) << "k"
              << "," << std::setw(w) << "kMerNumTh"
              << "," << std::setw(w) << "baseOverlapTh"
              << "," << std::setw(w) << "totalPos"
              << "," << std::setw(w) << "totalNeg"
              << "," << std::setw(w) << "falsePos"
              << "," << std::setw(w) << "falseNeg"
              << "," << std::setw(w) << "numOverlaps" << std::endl;
    std::cout << std::setw(w) << k << "," << std::setw(w) << kMerNumTh << ","
              << std::setw(w) << baseOverlapTh << "," << std::setw(w)
              << numPositives << "," << std::setw(w) << numNegatives << ","
              << std::setw(w) << falsePositives << "," << std::setw(w)
              << falseNegatives << "," << std::setw(w)
              << truePositives + falseNegatives << std::endl;
    delete rA;
}
