#include "ReadData.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits> // std::numeric_limits
#include <stdexcept>

void ReadData::loadFromFile(const char *fileName, enum Filetype filetype) {
    switch (filetype) {
    case READ:
        loadFromReadFile(fileName);
        break;
    case FASTQ:
        loadFromFastqFile(fileName);
        break;
    default:
        assert(false);
    }
}

void ReadData::loadFromReadFile(const char *fileName) {
    numReads = 0;
    readData.clear();
    readPos.clear();
    editStrings.clear();
    reverse.clear();

    std::ifstream infile(fileName);
    std::string line;
    size_t totalNumBases = 0;
    maxReadLen = 0;
    while (std::getline(infile, line)) {
        //        std::cout << line << std::endl;
        size_t index = line.find(':');
        assert(index != std::string::npos);
        reverse.push_back(line[index - 1] == 'c');
        size_t last_index = index;
        index = line.find(':', last_index + 1);
        assert(index != std::string::npos);
        ssize_t readLen = 0;
        readPos.push_back(std::stol(line.substr(last_index + 1, index)));
        {
            std::unique_ptr<std::string> ptr(
                new std::string(line.substr(index + 1)));
            editStrings.push_back(std::move(ptr));
        }
        //        std::cout << readPos.back() << std::endl;
        //        std::cout << editStrings.back() << std::endl;
        std::getline(infile, line);
        {
            // std::unique_ptr<std::string> ptr(new std::string(line));
            // readData.push_back(std::move(ptr));
            const char *cstr = line.c_str();
            readLen = std::strlen(cstr);
            std::unique_ptr<DnaBitset> ptr(new DnaBitset(cstr, readLen));
            readData.push_back(std::move(ptr));
            readPos.push_back(0);
        }
        // totalNumBases += readData.back()->size();
        // if (readData.back()->size() > maxReadLen)
        //     maxReadLen = readData.back()->size();
        totalNumBases += readLen;
        if (readLen > maxReadLen)
            maxReadLen = readLen;

        numReads++;
        if (numReads == std::numeric_limits<read_t>::max()) {
            throw std::runtime_error(
                "Too many reads for read_t type to handle.");
        }
    }
    assert(numReads != 0);
    avgReadLen = totalNumBases / numReads;
#ifdef DEBUG
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "avgReadLen " << avgReadLen << std::endl;
    std::cout << "maxReadLen " << maxReadLen << std::endl;
#endif
}

void ReadData::loadFromFastqFile(const char *fileName) {
    numReads = 0;
    readData.clear();
    readPos.clear();
    editStrings.clear();
    reverse.clear();

    std::ifstream infile(fileName);
    std::string line;
    size_t totalNumBases = 0;
    maxReadLen = 0;
    ssize_t readLen = 0;
    while (std::getline(infile, line)) {
        std::getline(infile, line);
        {
            // std::unique_ptr<std::string> ptr(new std::string(line));
            //convert to c style string  
            const char *cstr = line.c_str();
            // readData.push_back(std::move(ptr));
            readLen = std::strlen(cstr);
            std::unique_ptr<DnaBitset> ptr(new DnaBitset(cstr, readLen));
            readData.push_back(std::move(ptr));
            readPos.push_back(0);

        }
        {
            std::unique_ptr<std::string> ptr(new std::string());
            editStrings.push_back(std::move(ptr));
        }
        std::getline(infile, line);
        std::getline(infile, line);
        // totalNumBases += readData.back()->size();
        // if (readData.back()->size() > maxReadLen)
        //     maxReadLen = readData.back()->size();
        totalNumBases += readLen;
        if (readLen > maxReadLen)
            maxReadLen = readLen;
        numReads++;
        if (numReads == std::numeric_limits<read_t>::max()) {
            throw std::runtime_error(
                "Too many reads for read_t type to handle.");
        }
    }
    assert(numReads != 0);
    avgReadLen = totalNumBases / numReads;
#ifdef DEBUG
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "avgReadLen " << avgReadLen << std::endl;
    std::cout << "maxReadLen " << maxReadLen << std::endl;
#endif
}

read_t ReadData::getNumReads() { return numReads; }

void ReadData::getRead(read_t readId, std::string &readStr){
	readData[readId]->to_string(readStr);
}

std::vector<unsigned long> &ReadData::getReadPos() { return readPos; }

// std::vector<std::unique_ptr<std::string>> &ReadData::getReadData() {
//     return readData;
// }
std::vector<std::unique_ptr<DnaBitset>> &ReadData::getReadData() {
    return readData;
}

/// TODO: profile and maybe optimize
char ReadData::toComplement(char base) {
    switch (base) {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'C':
        return 'G';
    case 'G':
        return 'C';
    default:
        return base;
    }
}
