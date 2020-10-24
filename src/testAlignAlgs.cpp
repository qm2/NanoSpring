#include "AlignerTester.h"
#include "Edits.h"
#include "LocalMyers.h"
#include "LocalMyersRollBack.h"
#include "LocalMyersRollBack_impl.h"
#include "LocalMyers_impl.h"
#include "StringAligner.h"
//#include "../include/LocalMyersRollBackOld.h"
#include <iomanip>
#include <iostream>

void testAlg(std::vector<StringAligner<const char *> *> aligners);

void printString(std::string::iterator *begin, std::string::iterator *end) {
    for (; *begin != *end; (*begin)++)
        std::cout << **begin;
    std::cout << '\n';
}

int main(int argc, char **argv) {
    std::vector<StringAligner<const char *> *> aligners;

    //    aligners.push_back(new LocalMyers(50));
    //    aligners.push_back(new LocalMyers(100));
    //    aligners.push_back(new LocalMyers(200));
    //    aligners.push_back(new LocalMyers(32, 64));
    // aligners.push_back(new LocalMyers<const char *>(50, 100));
    aligners.push_back(
        new LocalMyersRollBack<const char *>(100, 200, 300, 0.21));
    // aligners.push_back(new MyersAligner());
    testAlg(aligners);
    for (StringAligner<const char *> *aligner : aligners)
        delete aligner;
}

void testAlg(std::vector<StringAligner<const char *> *> aligners) {
    AlignerTester aT;
    // ssize_t offsets2Test[] = {0,    100, -100, 200,  -200, 400,
    //                           -400, 800, -800, 1600, -1600};
    // ssize_t offsets2Test[] = {0, -5, 5, -10, 10, -20, 20, -40, 40};
    ssize_t offsets2Test[] = {0, 100, -100, 200, -200};

    const size_t algW = 35;
    for (ssize_t offset : offsets2Test) {
        aT.generateData(10000, offset, 1000, 0.03, 0.03, 0.04);
        for (StringAligner<const char *> *aligner : aligners) {
            if (aT.validate(aligner)) {
                std::cout << "Validation of " << std::setw(algW)
                          << aligner->name << " succeeded." << std::endl;
            } else {
                std::cout << "Validation of " << std::setw(algW)
                          << aligner->name << " failed." << std::endl;
            }
        }
    }
    std::cout << std::endl;
    for (ssize_t offset : offsets2Test) {
        aT.generateData(10000, offset, 1, 0.03, 0.03, 0.04);
        //        aT.generateData(10000, offset, 8, 0.01, 0.01, 0);
        std::cout << std::setw(algW) << "Algorithm"
                  << " " << std::setw(6) << "Offset"
                  << " " << std::setw(7) << "Time"
                  << " " << std::setw(10) << "SuccessRate"
                  << " " << std::setw(10) << "AvgEditDis"
                  << " " << std::setw(14) << "AvgBeginOffset"
                  << " " << std::setw(12) << "AvgEndOffset"
                  << " " << std::endl;
        for (StringAligner<const char *> *aligner : aligners) {
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