#include "Decompressor.h"
#include "Consensus.h"
#include "bsc_helper.h"
#include <arpa/inet.h>
#include <boost/filesystem.hpp>
#include <cassert>
// #include <filesystem>
#include <fstream>

void Decompressor::decompress(const char *inputFileName,
                              const char *outputFileName) {
    prepareTempDirs();

    untar(inputFileName);

    bscDecompress();

    generateReads(outputFileName);
}

void Decompressor::prepareTempDirs() const {
    boost::system::error_code ec;
    const boost::filesystem::path tempDirPath(tempDir);
    const boost::filesystem::path compressedTempDirPath(compressedTempDir);
    boost::filesystem::remove_all(tempDirPath, ec);
    boost::filesystem::remove_all(compressedTempDirPath, ec);
    boost::filesystem::create_directory(tempDirPath, ec);
    boost::filesystem::create_directory(compressedTempDirPath, ec);
}

void Decompressor::untar(const char *inputFileName) const {
    std::cout << "Extracting tar archive ...";
    std::string infile = inputFileName;
    std::string tar_command = "tar -C " + compressedTempDir + " -xvf " + infile;
    int tar_status = std::system(tar_command.c_str());
    if (tar_status != 0)
        throw std::runtime_error(
            "Error occurred during tar archive extraction.");
    std::cout << "Tar extraction done!\n";
}

void Decompressor::bscDecompress() {
    {
        std::ifstream metaData;
        metaData.open(compressedTempDir + "metaData");
        metaData >> numContigs >> numReads;
        std::cout << "numContigs " << numContigs << " numReads " << numReads
                  << '\n';
    }
#pragma omp parallel for
    for (size_t i = 0; i < numContigs; ++i)
        bscDecompress(tempFilename + std::to_string(i));
}

void Decompressor::bscDecompress(const std::string &f) const {
    const char *const extensions[] = {
        ".genome",       ".base",          ".id", ".pos", ".type",
        ".unalignedIds", ".unalignedReads"};
    const size_t numExtensions = sizeof(extensions) / sizeof(extensions[0]);
#pragma omp parallel for
    for (size_t i = 0; i < numExtensions; ++i) {
        std::string infile =
            compressedTempDir + f + extensions[i] + "Compressed";
        std::string outfile = tempDir + f + extensions[i];
        bsc::BSC_decompress(infile.c_str(), outfile.c_str());
    }
}

void Decompressor::generateReads(const char *outputFileName) const {
    std::string reads[numReads];
    for (size_t i = 0; i < numContigs; ++i)
        generateReads(reads, i);
    std::ofstream outFile;
    outFile.open(outputFileName);
    for (size_t i = 0; i < numReads; ++i)
        outFile << reads[i] << '\n';
}

void Decompressor::generateReads(std::string *reads, size_t contigId) const {
    std::string currentFilename =
        tempDir + tempFilename + std::to_string(contigId);
    // First we read the genome
    std::string genome;
    {
        std::ifstream genomeFile;
        genomeFile.open(currentFilename + ".genome");
        genomeFile >> genome;
        genomeFile.close();
    }
    std::ifstream idFile, posFile, editTypeFile, editBaseFile;
    idFile.open(currentFilename + ".id", std::ios::binary);
    posFile.open(currentFilename + ".pos", std::ios::binary);
    editTypeFile.open(currentFilename + ".type");
    editBaseFile.open(currentFilename + ".base");
    POS_T id = 0;
    while (true) {
        POS_T idInc;
        idFile.read(reinterpret_cast<char *>(&idInc), sizeof(idInc));
        if (!idFile)
            break;
        idInc = ntohl(idInc);
        id = id + idInc;
        generateRead(genome, reads[id], posFile, editTypeFile, editBaseFile);
        // std::cout << id << " ";
    }
    idFile.close();
    posFile.close();
    editTypeFile.close();
    editBaseFile.close();
    // Now we deal with unaligned reads
    std::ifstream unalignedIdsFile, unalignedReadsFile;
    unalignedIdsFile.open(currentFilename + ".unalignedIds", std::ios::binary);
    unalignedReadsFile.open(currentFilename + ".unalignedReads");
    id = 0;
    while (true) {
        POS_T idInc;
        unalignedIdsFile.read(reinterpret_cast<char *>(&idInc), sizeof(idInc));
        if (!unalignedIdsFile)
            break;
        idInc = ntohl(idInc);
        id = id + idInc;
        unalignedReadsFile >> reads[id];
        // std::cout << id << " ";
    }

    unalignedIdsFile.close();
    unalignedReadsFile.close();
}

void Decompressor::generateRead(const std::string &genome, std::string &read,
                                std::ifstream &posFile,
                                std::ifstream &editTypeFile,
                                std::ifstream &editBaseFile) const {
    POS_T curPos;
    posFile.read(reinterpret_cast<char *>(&curPos), sizeof(POS_T));
    curPos = ntohl(curPos);
    while (true) {
        // First we handle the unchanged bases
        POS_T numUnchanged;
        posFile.read(reinterpret_cast<char *>(&numUnchanged), sizeof(POS_T));
        numUnchanged = ntohl(numUnchanged);
        for (size_t i = 0; i < numUnchanged; ++i) {
            read.push_back(genome[curPos++]);
        }

        // Now for the edit
        char editType;
        editTypeFile.get(editType);
        if (editType == '\n')
            break;
        if (editType == 'd')
            curPos++;
        else if (editType == 'i') {
            char editBase;
            editBaseFile.get(editBase);
            read.push_back(editBase);
        } else if (editType == 's') {
            curPos++;
            char editBase;
            editBaseFile.get(editBase);
            read.push_back(editBase);
        }
    }
    // We need to read extra END_SYMBOL
    POS_T temp;
    posFile.read(reinterpret_cast<char *>(&temp), sizeof(POS_T));
    // Read the line breaks in editBase
    char c;
    editBaseFile.get(c);
    assert(c == '\n');

    // std::cout << read << '\n';
}