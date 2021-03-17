#include "Decompressor.h"
#include <iostream>
#include <omp.h>
#include <cstdlib>

int main(int argc, char **argv) {
    omp_set_nested(1);

    if (argc < 3) {
        std::cout << "Usage ./testDecompressor inputFilename outputFilename []"
                  << std::endl;
        return 1;
    }
    int numThr = omp_get_max_threads();
    if (argc == 4)
        numThr = std::atoi(argv[3]);
    Decompressor dc;
    dc.decompress(argv[1], argv[2], numThr);
    return 0;
}
