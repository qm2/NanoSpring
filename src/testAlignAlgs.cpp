#include <iostream>
#include <iomanip>
#include "AlignerTester.h"
#include "myers.h"
#include "Edits.h"

void testAlg(std::vector<StringAligner *> aligners);

int main(int argc, char **argv) {
    std::vector<StringAligner *> aligners;
    aligners.push_back(new PiecewiseMyers(50));
    aligners.push_back(new PiecewiseMyers(100));
    aligners.push_back(new PiecewiseMyers(200));
    aligners.push_back(new PiecewiseMyers(32, 64));
    aligners.push_back(new PiecewiseMyers(50, 100));
    aligners.push_back(new PiecewiseMyers(100, 200));
    aligners.push_back(new MyersAligner());
    testAlg(aligners);
    for (StringAligner *aligner : aligners)
        delete aligner;
}

void testAlg(std::vector<StringAligner *> aligners) {
    AlignerTester aT;
    size_t offsets2Test[] = {0, 100, 200, 400, 800};
    aT.generateData(10000, 0, 1, 0.03, 0.03, 0.04);
    const size_t algW = 35;
    for (StringAligner *aligner : aligners) {
        if (aT.validate(aligner)) {
            std::cout << "Validation of " << std::setw(algW) << aligner->name << " succeeded." << std::endl;
        } else {
            std::cout << "Validation of " << std::setw(algW) << aligner->name << " failed." << std::endl;
        }
    }

    for (size_t offset : offsets2Test) {
        aT.generateData(10000, offset, 4, 0.03, 0.03, 0.04);
//        aT.generateData(10000, offset, 8, 0.01, 0.01, 0);
        for (StringAligner *aligner : aligners) {
            double duration, successRate, avgEditDis;
            aT.profile(aligner, duration, successRate, avgEditDis);
            std::cout << "Algorithm " << std::setw(algW) << aligner->name << " "
                      << "Offset " << std::setw(6) << offset << " "
                      << "Time " << std::setw(10) << duration << "s "
                      << "Success Rate " << std::setw(6) << successRate << " "
                      << "Average Edit Distance " << std::setw(6) << avgEditDis
                      << std::endl;
        }
    }
}