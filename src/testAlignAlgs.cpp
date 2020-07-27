#include <iostream>
#include "AlignerTester.h"
#include "myers.h"
#include "Edits.h"

void testAlg(StringAligner *aligner);

int main(int argc, char **argv) {
    testAlg(new MyersAligner());
}

void testAlg(StringAligner *aligner) {
    AlignerTester aT;
    size_t offsets2Test[] = {0, 100, 200};
    aT.generateData(10000, 0, 1, 0.03, 0.03, 0.04);
    if (aT.validate(aligner)) {
        std::cout << "Validation of " << aligner->name << " successed." << std::endl;
    } else {
        std::cout << "Validation of " << aligner->name << " failed." << std::endl;
    }

    for (size_t offset : offsets2Test) {
        aT.generateData(10000, offset, 8, 0.03, 0.03, 0.04);
//        aT.generateData(10000, offset, 8, 0.01, 0.01, 0);

        double duration, successRate, avgEditDis;
        aT.profile(aligner, duration, successRate, avgEditDis);
        std::cout << "Algorithm " << aligner->name << " "
                  << "Offset " << offset << " "
                  << "Time " << duration << "ms "
                  << "Success Rate " << successRate << " "
                  << "Average Edit Distance " << avgEditDis
                  << std::endl;
    }
    delete aligner;
}