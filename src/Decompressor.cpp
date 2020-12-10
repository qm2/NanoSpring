#include "Decompressor.h"
#include "Consensus.h"
#include "DirectoryUtils.h"
#include "ReadData.h"
#include "bsc_helper.h"
#include <boost/filesystem.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>

void Decompressor::decompress(const char *inputFileName,
                              const char *outputFileName) {
    DirectoryUtils::clearDir(tempDir);

    // bsc::BSC_decompress(inputFileName, tarFileName.c_str());

    // Untar
    {
        std::cout << "Extracting tar archive ...";
        std::string tar_command =
            "tar -C " + tempDir + " -xf " + inputFileName;
        int tar_status = std::system(tar_command.c_str());
        if (tar_status != 0)
            throw std::runtime_error(
                "Error occurred during tar archive extraction.");
        std::cout << "Tar extraction done!\n";
    }

    // bsc Decompress
    {
        std::cout << "BSC decompression ...";
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
        std::cout << "BSC decompression done!\n";
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
        } else if (line.substr(0, delimPos) == "numThr") {
            numEncodingThreads = std::stol(line.substr(delimPos + 1));
            std::cout << "NumEncodingThreads:" << numEncodingThreads << std::endl;
        }
    }

    // since the compressed files are on a per-thread level, we first combine the
    // files with common extension (TODO: use better solution with parallelized decompression)
    std::set<std::string> extensions;
    DirectoryUtils::getAllExtensions(
            tempDir, std::inserter(extensions, extensions.end()));
    // Now we combine the files (note delimiter already present so not added here)
    for (const std::string &ext : extensions)
        DirectoryUtils::combineFilesWithExt(tempDir + tempFilename, ext, numEncodingThreads, false);
    
    // Unpack all files in directory with extensions
    boost::filesystem::directory_iterator endIt;
    for (boost::filesystem::directory_iterator it(tempDir); it != endIt; ++it) {
        if (!boost::filesystem::is_regular_file(*it))
            continue;
        boost::filesystem::path fullPath = it->path();
        std::string filename = fullPath.filename().string();
        const auto &ext = fullPath.extension();
        if (ext.empty())
            continue;
        DirectoryUtils::unpack(fullPath.string());
    }
}

void Decompressor::generateReads(const char *outputFileName) const {
    std::string reads[numReads];
    for (size_t i = 0; i < numContigs; ++i)
        generateReads(reads, i);
    std::ofstream outFile;
    outFile.open(outputFileName);
    for (read_t i = 0; i < numReads; ++i)
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
    std::ifstream idFile, posFile, editTypeFile, editBaseFile, complementFile;
    idFile.open(currentFilename + ".id");
    posFile.open(currentFilename + ".pos");
    editTypeFile.open(currentFilename + ".type");
    editBaseFile.open(currentFilename + ".base");
    complementFile.open(currentFilename + ".complement");
    read_t id = 0;
    while (true) {
        read_t idInc;
        idFile >> idInc;
        char c;
        idFile.get(c);
        if (!idFile)
            break;
        id = id + idInc;
        generateRead(genome, reads[id], posFile, editTypeFile, editBaseFile,
                     complementFile);
        // std::cout << id << " ";
    }
    idFile.close();
    posFile.close();
    editTypeFile.close();
    editBaseFile.close();
    complementFile.close();

    // Now we deal with unaligned reads
    // std::ifstream unalignedIdsFile, unalignedReadsFile;
    // unalignedIdsFile.open(currentFilename + ".unalignedIds",
    // std::ios::binary); unalignedReadsFile.open(currentFilename +
    // ".unalignedReads"); id = 0; while (true) {
    //     read_t idInc;
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
                                std::ifstream &editBaseFile,
                                std::ifstream &complementFile) const {
    size_t curPos;
    char c;
    posFile >> curPos;
    posFile.get(c);
    complementFile.get(c);
    bool reverseComplement = (c == 'c');
    complementFile.get(c);
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
    if (reverseComplement) {
        std::string temp;
        temp.swap(read);
        ReadData::toReverseComplement(temp.begin(), temp.end(),
                                      std::inserter(read, read.end()));
    }
    // We need to read extra '\n'
    posFile.get(c);
    // Read the line breaks in editBase
    editBaseFile.get(c);
    assert(c == '\n');

    // std::cout << read << '\n';
}
