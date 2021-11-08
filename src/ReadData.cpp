#include "ReadData.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits> // std::numeric_limits
#include <stdexcept>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

void ReadData::loadFromFile(const char *fileName, enum Filetype filetype, bool low_mem) {
    switch (filetype) {
        case READ:
            loadFromReadFile(fileName);
            break;
        case FASTQ:
            loadFromFastqFile(fileName, false, low_mem);
            break;
        case GZIP:
            loadFromFastqFile(fileName, true, low_mem);
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
        size_t index = line.find(':');
        assert(index != std::string::npos);
        reverse.push_back(line[index - 1] == 'c');
        size_t last_index = index;
        index = line.find(':', last_index + 1);
        assert(index != std::string::npos);
        size_t readLen = 0;
        readPos.push_back(std::stol(line.substr(last_index + 1, index)));
        {
            std::unique_ptr<std::string> ptr(
                new std::string(line.substr(index + 1)));
            editStrings.push_back(std::move(ptr));
        }
        std::getline(infile, line);
        {
            readLen = line.size();
            const char *cstr = line.c_str();
            std::unique_ptr<DnaBitset> ptr(new DnaBitset(cstr, readLen));
            readData.push_back(std::move(ptr));
            readPos.push_back(0);
        }
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
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "avgReadLen " << avgReadLen << std::endl;
    std::cout << "maxReadLen " << maxReadLen << std::endl;
}

void ReadData::loadFromFastqFile(const char *fileName, bool gzip_flag, bool low_mem) {
    if (low_mem) {
        loadFromFastqFile_lowmem(fileName, gzip_flag);
    } else {
        loadFromFastqFile_highmem(fileName, gzip_flag);
    }
}

void ReadData::loadFromFastqFile_highmem(const char *fileName, bool gzip_flag) {
    reads_in_memory = true;
    numReads = 0;
    readData.clear();
    readPos.clear();
    editStrings.clear();
    reverse.clear();
    
    std::ifstream infile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> *inbuf;
    std::istream *fin = &infile;
    if (gzip_flag) {
      infile.open(fileName, std::ios_base::binary);
      inbuf =
          new boost::iostreams::filtering_streambuf<boost::iostreams::input>;
      inbuf->push(boost::iostreams::gzip_decompressor());
      inbuf->push(infile);
      fin = new std::istream(inbuf);
    } else {
      infile.open(fileName);
    }
    std::string line;
    size_t totalNumBases = 0;
    maxReadLen = 0;
    size_t numReadsPerBlock = 5000;
    std::vector<std::string> lines(numReadsPerBlock);
    size_t numReadsCurrBlock = 0;
    size_t numReadsInserted = 0;
    while (true) {
        while (std::getline(*fin, line)) {
            std::getline(*fin, lines[numReadsCurrBlock++]);
            auto readLen = lines[numReadsCurrBlock-1].size();
            totalNumBases += readLen;
            if (readLen > maxReadLen)
                maxReadLen = readLen;
            numReads++;
            if (numReads == std::numeric_limits<read_t>::max())
                throw std::runtime_error(
                        "Too many reads for read_t type to handle.");
            readPos.push_back(0);
            std::getline(*fin, line);
            std::getline(*fin, line);
            if (numReadsCurrBlock == numReadsPerBlock)
                break;
        }
        readData.resize(numReadsInserted + numReadsCurrBlock);
#pragma omp parallel for
        for (size_t i = 0; i < numReadsCurrBlock; i++) {
            std::unique_ptr<DnaBitset> ptr(new DnaBitset(
                        lines[i].c_str(), lines[i].size()));
            readData[numReadsInserted+i] = std::move(ptr);
        }
        numReadsInserted += numReadsCurrBlock;
        if (numReadsCurrBlock < numReadsPerBlock)
            break;
        numReadsCurrBlock = 0;
    }
    assert(numReads != 0);
    avgReadLen = totalNumBases / numReads;
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "avgReadLen " << avgReadLen << std::endl;
    std::cout << "maxReadLen " << maxReadLen << std::endl;
    if (gzip_flag) {
        delete fin;
        delete inbuf;
    }
    // close files
    infile.close();
}

void ReadData::loadFromFastqFile_lowmem(const char *fileName, bool gzip_flag) {
    reads_in_memory = false;
    numReads = 0;
    readData.clear();
    readPos.clear();
    editStrings.clear();
    reverse.clear();
    
    std::ifstream infile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> *inbuf;
    std::istream *fin = &infile;
    if (gzip_flag) {
      infile.open(fileName, std::ios_base::binary);
      inbuf =
          new boost::iostreams::filtering_streambuf<boost::iostreams::input>;
      inbuf->push(boost::iostreams::gzip_decompressor());
      inbuf->push(infile);
      fin = new std::istream(inbuf);
    } else {
      infile.open(fileName);
    }
    std::string line;
    size_t totalNumBases = 0;
    maxReadLen = 0;
    size_t cur_pos_in_bitset_file = 0;
    std::string bitsetFileNameFull = tempDir + "/" + bitsetFileName;
    std::ofstream bitset_fout(bitsetFileNameFull, std::ios::binary);
    while (std::getline(*fin, line)) {
        std::getline(*fin, line);
        auto readLen = line.size();
        totalNumBases += readLen;
        read_lengths.push_back(readLen);
        if (readLen > maxReadLen)
            maxReadLen = readLen;
        numReads++;
        if (numReads == std::numeric_limits<read_t>::max())
            throw std::runtime_error(
                    "Too many reads for read_t type to handle.");
        readPos.push_back(0);
        auto read_bitset = new DnaBitset(line.c_str(), readLen);
        // write to bitset file
        size_t bitset_len = read_bitset->to_file(bitset_fout);
        delete read_bitset;
        read_pos_in_file.push_back(cur_pos_in_bitset_file);
        cur_pos_in_bitset_file += bitset_len;
        std::getline(*fin, line);
        std::getline(*fin, line);
    }
    assert(numReads != 0);
    avgReadLen = totalNumBases / numReads;
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "avgReadLen " << avgReadLen << std::endl;
    std::cout << "maxReadLen " << maxReadLen << std::endl;
    if (gzip_flag) {
        delete fin;
        delete inbuf;
    }
    // close files
    infile.close();
    bitset_fout.close();

    // open bitset file
    fin_bitset = std::unique_ptr<std::ifstream>(new std::ifstream(bitsetFileNameFull, std::ios::binary));
}

read_t ReadData::getNumReads() { return numReads; }

void ReadData::getRead(read_t readId, std::string &readStr){
    if (!reads_in_memory) {
        fin_bitset_mtx.lock();
        fin_bitset->seekg(read_pos_in_file[readId]);
        auto bitset = new DnaBitset(*fin_bitset.get(), read_lengths[readId]);
        bitset->to_string(readStr);
        delete bitset;
        fin_bitset_mtx.unlock();
    } else {
    	readData[readId]->to_string(readStr);
    }
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

ReadData::~ReadData() {
    if (!reads_in_memory) {
        fin_bitset->close();
        std::string fileName = tempDir + "/" + bitsetFileName;
        remove(fileName.c_str());
    }
}
