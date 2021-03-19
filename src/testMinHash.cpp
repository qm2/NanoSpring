#include "ReadData.h"
#include "ReadFilter.h"
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
        ReadData rD;
        rD.loadFromFile(argv[1]);
        std::cout << rD.getNumReads() << " reads\n";
        MinHashReadFilter *minHashReadFilter = new MinHashReadFilter();
        minHashReadFilter->k = std::stoi(argv[2]);
        minHashReadFilter->n = std::stoi(argv[3]);
        ReadFilter *rF = minHashReadFilter;
        rF->initialize(rD);
        std::cout << "calculating filter stats" << std::endl;
        while (true) {
            unsigned int overlapBaseTh, overlapSketchTh;
            std::cin >> overlapBaseTh >> overlapSketchTh;
            if (overlapBaseTh == 0)
                break;
            std::cout << rF->getFilterStats(overlapBaseTh, overlapSketchTh)
                      << std::endl;
        }
        delete rF;
    }

    return 0;
}
