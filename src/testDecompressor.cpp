#include "Decompressor.h"
#include "DirectoryUtils.h"
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

    // generate randomly named temporary directory in the working directory
    std::string temp_dir;
    while (true) {
      std::string random_str = "tmp." + DirectoryUtils::random_string(10);
      temp_dir = "./" + random_str + "/";
      if (!boost::filesystem::exists(temp_dir)) break;
    }
    if (!boost::filesystem::create_directory(temp_dir)) {
      throw std::runtime_error("Cannot create temporary directory.");
    }
    std::cout << "Temporary directory: " << temp_dir << "\n";
    Decompressor dc;
    dc.tempDir = temp_dir;
    dc.decompress(argv[1], argv[2], numThr);
    return 0;
}
