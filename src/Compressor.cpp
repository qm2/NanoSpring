#include "Compressor.h"
#include "Consensus.h"
#include "Contig.h"
#include "DirectoryUtils.h"
#include "ReadData.h"
#include "ReadFilter.h"
#include "bsc_helper.h"
#include <boost/filesystem.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <ios>
#include <unistd.h>
#include <string>

void mem_usage(double& vm_usage, double& resident_set) {
   // from https://www.tutorialspoint.com/how-to-get-memory-usage-at-runtime-using-cplusplus
   using namespace std;
   vm_usage = 0.0;
   resident_set = 0.0;
   ifstream stat_stream("/proc/self/stat",ios_base::in); //get info from proc
   // directory
   //create some variables to get info
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;
   unsigned long vsize;
   long rss;
   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
   >> utime >> stime >> cutime >> cstime >> priority >> nice
   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care
   // about the rest
   stat_stream.close();
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured
   // to use 2MB pages
   vm_usage = vsize / 1024.0;
   resident_set = rss * page_size_kb;
   cout << "resident_set: " << resident_set << " KB\n";
}

void Compressor::compress(const char *inputFileName, const int numThr) const {
	std::cout<<"number of threads: "<<numThr<<std::endl;
	double vm_usage, resident_set;
	std::cout << "Entering compress():\n";
    mem_usage(vm_usage, resident_set);
    omp_set_num_threads(numThr);
    ReadData rD;
    rD.loadFromFile(inputFileName, filetype);
    std::cout << "After rD.loadFromFile():\n";
    mem_usage(vm_usage, resident_set);

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
    std::cout << "After minhash computation:\n";
    mem_usage(vm_usage, resident_set);
    // We clear the temp directories and create them if they do not
    // exist
    DirectoryUtils::clearDir(tempDir);

    Consensus consensus;
    {
        consensus.rD = &rD;
        consensus.rF = &rF;
        consensus.rA = rA;
        consensus.aligner = aligner;
        consensus.tempDir = tempDir;
        consensus.tempFileName = tempFileName;
        consensus.numThr = numThr;
        //parameters for minimap2
        consensus.m_k = m_k;
        consensus.m_w = m_w;
        consensus.hashBits = hashBits; 

        consensus.generateAndWriteConsensus();
    }
    
    {
// We first compress all the files in the temp directory (we compress
// only files with extensions)
#ifdef DEBUG
        std::cout << "bsc compression starts" << std::endl;
#endif
        std::set<std::string> extensions;
        DirectoryUtils::getAllExtensions(tempDir, std::inserter(extensions, extensions.end()));
        for (const std::string &ext : extensions) {
            std::vector<size_t> uncompressedSizes(numThr);
            std::vector<size_t> compressedSizes(numThr);
#pragma omp parallel num_threads(numThr)
#pragma omp for
            for (int i = 0; i < numThr; i++) {
                std::string fullPath = tempDir + tempFileName + ".tid."+ std::to_string(i) + ext;
                std::string outPath = fullPath + "Compressed";
                bsc::BSC_compress(fullPath.c_str(), outPath.c_str());
                uncompressedSizes[i] = boost::filesystem::file_size(fullPath);
                compressedSizes[i] = boost::filesystem::file_size(outPath);
                boost::filesystem::remove(fullPath);
            }
            size_t totalUncompressed = 0, totalCompressed = 0;
            for (int i = 0; i < numThr; i++) {
                totalUncompressed += uncompressedSizes[i];
                totalCompressed += compressedSizes[i];
            }
            std::cout << "Extension " << ext << ": Compressed " << totalUncompressed << " bytes to " << totalCompressed << " bytes\n";
        }
    }

    std::cout << "Creating tar archive ..." << std::endl;
    std::string tar_command =
        "tar -cf " + outputFileName + " -C " + tempDir + " . ";
    int tar_status = std::system(tar_command.c_str());
    if (tar_status)
        throw std::runtime_error(
            "Error occurred during tar archive generation.");
    std::cout << "Tar archive done!\n";

    {
        std::string lsCommand = "ls -lh " + outputFileName;
        int ls_status = std::system(lsCommand.c_str());
        if (ls_status)
            throw std::runtime_error("Error occurred during ls command.");
    }
}
