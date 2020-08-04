#include "AlignerTester.h"
#include "Edits.h"
#include "myers.h"
#include <iomanip>
#include <iostream>

void testAlg(std::vector<StringAligner *> aligners);

void printString(std::string::iterator *begin, std::string::iterator *end) {
    for (; *begin != *end; (*begin)++)
        std::cout << **begin;
    std::cout << '\n';
}

int main(int argc, char **argv) {
    std::vector<StringAligner *> aligners;

    //    aligners.push_back(new LocalMyers(50));
    //    aligners.push_back(new LocalMyers(100));
    //    aligners.push_back(new LocalMyers(200));
    size_t maxEditDis = 3200;
    //    aligners.push_back(new LocalMyers(32, 64));
    //    aligners.push_back(new LocalMyersRollBack(32, 64, maxEditDis));
    // aligners.push_back(new LocalMyers(50, 100));
    // aligners.push_back(new LocalMyersRollBack(50, 100, maxEditDis));
    // aligners.push_back(new LocalMyersRollBack(50, 100, maxEditDis * 2));
    aligners.push_back(new LocalMyersRollBack(100, 200, maxEditDis));
    // aligners.push_back(new LocalMyersRollBack(100, 200, maxEditDis * 2));

    //    aligners.push_back(new LocalMyersRollBack(100, 50, maxEditDis));
    //    aligners.push_back(new LocalMyers(100, 200));
    //    aligners.push_back(new LocalMyersRollBack(100, 200, maxEditDis));
    //    aligners.push_back(new LocalMyersRollBack(200, 100, maxEditDis));
    // aligners.push_back(new MyersAligner());
    testAlg(aligners);
    for (StringAligner *aligner : aligners)
        delete aligner;
}

void testAlg(std::vector<StringAligner *> aligners) {
    AlignerTester aT;
    // ssize_t offsets2Test[] = {0,    100, -100, 200,  -200, 400,
    //                           -400, 800, -800, 1600, -1600};
    ssize_t offsets2Test[] = {0, -5, 5, -10, 10, -20, 20, -40, 40};
    //    size_t offsets2Test[] = {0, 100};

    aT.generateData(10000, 200, 1, 0.03, 0.03, 0.04);
    const size_t algW = 35;
    for (StringAligner *aligner : aligners) {
        if (aT.validate(aligner)) {
            std::cout << "Validation of " << std::setw(algW) << aligner->name
                      << " succeeded." << std::endl;
        } else {
            std::cout << "Validation of " << std::setw(algW) << aligner->name
                      << " failed." << std::endl;
        }
    }
    std::cout << std::endl;
    for (ssize_t offset : offsets2Test) {
        aT.generateData(10000, offset, 100, 0.03, 0.03, 0.04);
        //        aT.generateData(10000, offset, 8, 0.01, 0.01, 0);
        std::cout << std::setw(algW) << "Algorithm"
                  << " " << std::setw(6) << "Offset"
                  << " " << std::setw(7) << "Time"
                  << " " << std::setw(10) << "SuccessRate"
                  << " " << std::setw(10) << "AvgEditDis"
                  << " " << std::setw(14) << "AvgBeginOffset"
                  << " " << std::setw(12) << "AvgEndOffset"
                  << " " << std::endl;
        for (StringAligner *aligner : aligners) {
            double duration, successRate, avgEditDis, avgBeginOffset,
                avgEndOffset;
            aT.profile(aligner, duration, successRate, avgBeginOffset,
                       avgEndOffset, avgEditDis);
            std::cout << std::setw(algW) << aligner->name << " " << std::setw(6)
                      << offset << " " << std::setw(7) << std::fixed
                      << std::setprecision(4) << duration << "s "
                      << std::setw(10) << std::fixed << std::setprecision(3)
                      << successRate << " " << std::setw(10) << std::fixed
                      << std::setprecision(0) << avgEditDis << " "
                      << std::setw(14) << std::fixed << std::setprecision(0)
                      << avgBeginOffset << " " << std::setw(12) << std::fixed
                      << std::setprecision(0) << avgEndOffset << " "
                      << std::endl;
        }
        std::cout << std::endl;
    }
}