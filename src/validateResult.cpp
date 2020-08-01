#include <iostream>
#include <fstream>
#include <AlignerTester.h>

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "Usage: ./validateResult compressed.genome compressed.reads original.reads" << std::endl;
        return 1;
    }
    std::cout << "Validating\n";
    std::string genome;
    std::ifstream f;
    f.open(argv[1]);
    f >> genome;
    f.close();
    std::ifstream compressed, original;
    compressed.open(argv[2]);
    original.open(argv[3]);

    std::string read;
    size_t i = 0;
    bool success = true;
    while (std::getline(original, read)) {
        std::getline(original, read);
//        std::cout << "read " << read << "\n";
        std::string editScript;
        size_t beginPos;
        char c;
        compressed >> beginPos >> c >> beginPos;
//        std::cout << std::to_string(beginPos) << " " << c;
        compressed >> editScript;
//        std::getline(compressed, editScript);
//        std::cout << editScript << "\n";
        std::string result;
        AlignerTester::applyEditsToString(genome.substr(beginPos), editScript, result);
        if (result != read) {
            std::cout << "Read " << i << " failed\n";
            success = false;
        }
        ++i;
    }
    compressed.close();
    original.close();

    if (success) {
        std::cout << "Validation succeeded!\n";
    }
}
