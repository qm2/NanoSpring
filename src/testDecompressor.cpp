#include "Decompressor.h"
#include <iostream>
#include <omp.h>

int main(int argc, char **argv) {
    omp_set_nested(1);

    if (argc < 3) {
        std::cout << "Usage ./testDecompressor inputFilename outputFilename"
                  << std::endl;
        return 1;
    }
    Decompressor dc;
    dc.decompress(argv[1], argv[2]);
    return 0;
}
