#include "Compressor.h"
#include "Consensus.h"
#include "Contig.h"
#include "ReadData.h"
#include "ReadFilter.h"
#include "bsc_helper.h"
#include <boost/filesystem.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

void Compressor::compress(const char *inputFileName) const {

    ReadData rD;
    rD.loadFromFile(inputFileName);
    // (inputFileName, k, n);
    MinHashReadFilter rF;
    rF.k = k;
    rF.n = n;
    rF.overlapSketchThreshold = overlapSketchThreshold;
    {
        auto start = std::chrono::high_resolution_clock::now();
        rF.initialize(rD);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Calculating MinHashes took " << duration.count()
                  << " milliseconds" << std::endl;
    }

    // We clear the temp directories and create them if they do not
    // exist
    boost::system::error_code ec;
    const boost::filesystem::path tempDirPath(tempDir);
    boost::filesystem::remove_all(tempDirPath, ec);
    boost::filesystem::create_directory(tempDirPath, ec);

    Consensus consensus;
    {
        consensus.rD = &rD;
        consensus.rF = &rF;
        consensus.rA = rA;
        consensus.aligner = aligner;
        consensus.tempDir = tempDir;
        consensus.tempFileName = tempFileName;

        consensus.generateConsensus();
        consensus.writeConsensus();
    }

    std::cout << "Creating tar archive ..." << std::endl;
    std::string tarFileName = "originalFile";
    std::string tar_command =
        "tar -cvf " + tarFileName + " -C " + tempDir + " . ";
    int tar_status = std::system(tar_command.c_str());
    if (tar_status)
        throw std::runtime_error(
            "Error occurred during tar archive generation.");
    std::cout << "Tar archive done!\n";

    {
        std::string lsCommand = "ls -lh " + tarFileName;
        int ls_status = std::system(lsCommand.c_str());
        if (ls_status)
            throw std::runtime_error("Error occurred during ls command.");
    }

    bsc::BSC_compress(tarFileName.c_str(), outputFileName.c_str());

    {
        std::string lsCommand = "ls -lh " + outputFileName;
        int ls_status = std::system(lsCommand.c_str());
        if (ls_status)
            throw std::runtime_error("Error occurred during ls command.");
    }
}