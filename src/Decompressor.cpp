#include "Decompressor.h"
#include "Consensus.h"
#include "DirectoryUtils.h"
#include "ReadData.h"
#include "bsc_helper.h"
#include "dnaToBits.h"
#include <boost/filesystem.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <omp.h>
#include <vector>
#include <string>

void Decompressor::decompress(const char *inputFileName,
                              const char *outputFileName,
                              const int numThreads) {
    std::cout << "numDecodingThreads: " << numThreads << "\n";
    omp_set_num_threads(numThreads);

    DirectoryUtils::clearDir(tempDir);

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

    // extract information from the metaDataFile
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
    
    // bsc Decompress

    std::vector<std::string> suffix = {".genome",
                                      ".id",
                                      ".pos",
                                      ".type",
                                      ".base",
                                      ".complement",
                                      ".lone"};
    std::cout << "BSC decompression ...";
    //loop through each thread
#pragma omp parallel for
    for (size_t i = 0; i < numEncodingThreads; ++i){
        for (auto &s : suffix) {
            std::string uncompressedFile = tempDir+"/"+tempFilename+".tid."+std::to_string(i)+s;
            std::string compressedFile = uncompressedFile + "Compressed";
            bsc::BSC_decompress(compressedFile.c_str(), uncompressedFile.c_str());
        }
    }
    std::cout << "BSC decompression done!\n";
    
    //create a DnaBits vector (to reduce memory consumption of decompressor)
    std::vector<DnaBitset*> reads(numReads);
    std::cout << "Generating reads...";
    //loop through each thread
#pragma omp parallel for
    for (size_t i = 0; i < numEncodingThreads; ++i){
        std::string currentFilename = tempDir + tempFilename + ".tid." + std::to_string(i);
        //open all the files for this thread
        std::ifstream genomeFile;
        genomeFile.open(currentFilename + ".genome");
        std::ifstream idFile, posFile, editTypeFile, editBaseFile, complementFile, loneFile;
        idFile.open(currentFilename + ".id", std::ios::binary);
        posFile.open(currentFilename + ".pos", std::ios::binary);
        editTypeFile.open(currentFilename + ".type");
        editBaseFile.open(currentFilename + ".base");
        complementFile.open(currentFilename + ".complement");
        loneFile.open(currentFilename + ".lone");
        // read each genome in seris
        std::string genome;
        std::string currentRead;
        while(std::getline(genomeFile, genome))
        {
            //the id for the current read
            read_t id = 0;
            while (true) {
                char c;
                complementFile.get(c);
                //break after we cover all reads in contig
                if(c == '\n'){
                    break;
                }
                //check if it is complement or not
                bool reverseComplement = (c == 'c');
                //read in the id
                read_t idInc;
                idFile.read((char*)&idInc, sizeof(read_t));
                id = id + idInc;
                generateRead(genome, currentRead, posFile, editTypeFile, editBaseFile,
                             reverseComplement);
                reads[id] = new DnaBitset(currentRead.c_str(), currentRead.size());
            }
        }
        // now do the lone reads
        std::string loneRead;
        read_t id = 0, idInc;
        while(std::getline(loneFile, loneRead)) {
            // the id for the current read
            idFile.read((char*)&idInc, sizeof(read_t));
            id = id + idInc;
            reads[id] = new DnaBitset(loneRead.c_str(),loneRead.size());
        }
        //close all files
        genomeFile.close();
        idFile.close();
        posFile.close();
        editTypeFile.close();
        editBaseFile.close();
        complementFile.close();
        loneFile.close();
    }
    std::cout << "Done!\n";

    //output the reads to a file
    std::ofstream outFile;
    outFile.open(outputFileName);
    std::cout << "Writing to file...";
    std::string currentRead;
    // TODO: can we parallelize the to_string part somehow
    for (size_t i = 0; i < numReads; i++) {
        reads[i]->to_string(currentRead);
        delete reads[i];
        outFile << currentRead << '\n';
    }
    std::cout << "Done!\n";

    boost::filesystem::remove_all(tempDir);
}

void Decompressor::generateRead(const std::string &genome, std::string &read,
                                std::ifstream &posFile,
                                std::ifstream &editTypeFile,
                                std::ifstream &editBaseFile,
                                bool reverseComplement) const {
    read.clear();
    uint32_t curPos;

    // posFile.read((char*)&curPos,sizeof(uint32_t));
    curPos = DirectoryUtils::read_var_uint32(posFile);
    while (true) {
        // First we handle the unchanged bases
        uint32_t numUnchanged;
        // posFile.read((char*)&numUnchanged,sizeof(uint32_t));
        numUnchanged = DirectoryUtils::read_var_uint32(posFile);
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
}
