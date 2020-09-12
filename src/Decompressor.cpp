#include "Decompressor.h"
#include "Consensus.h"
#include "bsc_helper.h"
#include <boost/filesystem.hpp>
#include <cassert>
#include <fstream>
#include <iostream>

void Decompressor::decompress(const char *inputFileName,
                              const char *outputFileName) {
    prepareTempDirs();

    // bsc::BSC_decompress(inputFileName, tarFileName.c_str());

    // Untar
    {
        std::cout << "Extracting tar archive ...";
        std::string tar_command =
            "tar -C " + tempDir + " -xvf " + inputFileName;
        int tar_status = std::system(tar_command.c_str());
        if (tar_status != 0)
            throw std::runtime_error(
                "Error occurred during tar archive extraction.");
        std::cout << "Tar extraction done!\n";
    }

    // bsc Decompress
    {
        boost::filesystem::path tempPath(tempDir);
        boost::filesystem::directory_iterator endIt;
        for (boost::filesystem::directory_iterator it(tempDir); it != endIt;
             ++it) {
            if (boost::filesystem::is_regular_file(*it)) {
                const boost::filesystem::path &fullPath = it->path();
                const std::string filename = fullPath.filename().string();
                // We only decompress files ending with "Compressed"
                const std::string ending("Compressed");
                if (filename.length() < ending.length())
                    continue;
                if (filename.compare(filename.length() - ending.length(),
                                     ending.length(), ending) != 0)
                    continue;
                const std::string &fullPathString = fullPath.string();
                std::string outPath = fullPathString.substr(
                    0, fullPathString.length() - ending.length());
                bsc::BSC_decompress(fullPath.string().c_str(), outPath.c_str());
                boost::filesystem::remove(fullPath);
            }
        }
    }

    unpack();

    generateReads(outputFileName);
}

void Decompressor::unpack() {
    std::cout << "Unpacking..." << std::endl;
    std::ifstream metaDataFile;
    metaDataFile.open(tempDir + "metaData");
    std::string line;

    while (std::getline(metaDataFile, line)) {
        size_t delimPos = line.find('=');
        if (delimPos == std::string::npos)
            continue;
        if (line.substr(0, delimPos) == "numReads") {
            numReads = std::stol(line.substr(delimPos + 1));
            std::cout << "NumReads:" << numReads << std::endl;
        } else if (line.substr(0, delimPos) == "numContigs") {
            numContigs = std::stol(line.substr(delimPos + 1));
            std::cout << "NumContigs:" << numContigs << std::endl;
        }
    }

    const char *const extensions[] = {".genome", ".base", ".id", ".pos",
                                      ".type"};
    const size_t numExtensions = sizeof(extensions) / sizeof(extensions[0]);
    for (size_t extId = 0; extId < numExtensions; extId++) {
        size_t i = 0;
        std::ifstream inFile;
        inFile.open(tempDir + tempFilename + extensions[extId]);
        std::ofstream outFile;
        std::string line;
        bool need2CreateFile = true;
        while (std::getline(inFile, line)) {
            if (need2CreateFile) {
                need2CreateFile = false;
                outFile.open(tempDir + tempFilename + std::to_string(i) +
                             extensions[extId]);
                ++i;
                outFile << line << '\n';
            } else {
                if (line == ".") {
                    outFile.close();
                    need2CreateFile = true;
                } else
                    outFile << line << '\n';
            }
        }
    }
}

void Decompressor::prepareTempDirs() const {
    boost::system::error_code ec;
    const boost::filesystem::path tempDirPath(tempDir);
    boost::filesystem::remove_all(tempDirPath, ec);
    boost::filesystem::create_directory(tempDirPath, ec);
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
    idFile.open(currentFilename + ".id");
    posFile.open(currentFilename + ".pos");
    editTypeFile.open(currentFilename + ".type");
    editBaseFile.open(currentFilename + ".base");
    size_t id = 0;
    while (true) {
        size_t idInc;
        idFile >> idInc;
        char c;
        idFile.get(c);
        if (!idFile)
            break;
        id = id + idInc;
        generateRead(genome, reads[id], posFile, editTypeFile, editBaseFile);
        // std::cout << id << " ";
    }
    idFile.close();
    posFile.close();
    editTypeFile.close();
    editBaseFile.close();

    // Now we deal with unaligned reads
    // std::ifstream unalignedIdsFile, unalignedReadsFile;
    // unalignedIdsFile.open(currentFilename + ".unalignedIds",
    // std::ios::binary); unalignedReadsFile.open(currentFilename +
    // ".unalignedReads"); id = 0; while (true) {
    //     size_t idInc;
    //     unalignedIdsFile >> idInc;
    //     char c;
    //     unalignedIdsFile.get(c);
    //     if (!unalignedIdsFile)
    //         break;
    //     id = id + idInc;
    //     unalignedReadsFile >> reads[id];
    //     // std::cout << id << " ";
    // }
    // unalignedIdsFile.close();
    // unalignedReadsFile.close();
}

void Decompressor::generateRead(const std::string &genome, std::string &read,
                                std::ifstream &posFile,
                                std::ifstream &editTypeFile,
                                std::ifstream &editBaseFile) const {
    size_t curPos;
    char c;
    posFile >> curPos;
    posFile.get(c);
    while (true) {
        // First we handle the unchanged bases
        size_t numUnchanged;
        posFile >> numUnchanged;
        posFile.get(c);
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
    // We need to read extra '\n'
    posFile.get(c);
    // Read the line breaks in editBase
    editBaseFile.get(c);
    assert(c == '\n');

    // std::cout << read << '\n';
}