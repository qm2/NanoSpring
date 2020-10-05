#include "ReadData.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>       // std::numeric_limits
#include <stdexcept>

void ReadData::loadFromFile(const char *fileName) {
    numReads = 0;
    readData.clear();
    readPos.clear();
    editStrings.clear();

    std::ifstream infile(fileName);
    std::string line;
    while (std::getline(infile, line)) {
        //        std::cout << line << std::endl;
        size_t index = line.find(':');
        readPos.push_back(std::stol(line.substr(0, index)));
        {
            std::unique_ptr<std::string> ptr(
                new std::string(line.substr(index + 1)));
            editStrings.push_back(std::move(ptr));
        }
        //        std::cout << readPos.back() << std::endl;
        //        std::cout << editStrings.back() << std::endl;
        std::getline(infile, line);
        {
            std::unique_ptr<std::string> ptr(new std::string(line));
            readData.push_back(std::move(ptr));
        }
        numReads++;
        if (numReads == std::numeric_limits<read_t>::max()) {
            throw std::runtime_error("Too many reads for read_t type to handle.");
        }
    }
    readLen = readData[0]->length();
    readPosSorted = readPos;
    std::sort(readPosSorted.begin(), readPosSorted.end());
#ifdef DEBUG
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "readLen " << readLen << std::endl;
#endif
}

read_t ReadData::getNumReads() { return numReads; }

std::string &ReadData::getRead(read_t readId) { return *readData[readId]; }

std::vector<unsigned long> &ReadData::getReadPos() { return readPos; }

std::vector<unsigned long> &ReadData::getReadPosSorted() {
    return readPosSorted;
}

std::vector<std::unique_ptr<std::string>> &ReadData::getReadData() {
    return readData;
}
