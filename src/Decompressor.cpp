#include "Decompressor.h"
#include "Consensus.h"
#include "DirectoryUtils.h"
#include "ReadData.h"
#include "bsc_helper.h"
#include "lzma2_helper.h"
#include "dnaToBits.h"
#include <boost/filesystem.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <omp.h>
#include <vector>
#include <string>
#include <chrono>

void Decompressor::decompress(const char *inputFileName,
                              const char *outputFileName,
                              const int numThreads,
                              const int decompression_memory_gb) {
    std::cout << "Decompression ...\n";
    auto overall_start = std::chrono::high_resolution_clock::now();
    std::cout << "numDecodingThreads: " << numThreads << "\n";
    std::cout << "Aiming to use close to or less than " << decompression_memory_gb << " GB peak memory.\n";
    omp_set_nested(1);
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
    std::cout << "BSC/LZMA2 decompression ...\n";
    auto bsc_start = std::chrono::high_resolution_clock::now();
    //loop through each thread
#pragma omp parallel for
    for (size_t i = 0; i < numEncodingThreads; ++i){
        for (auto &s : suffix) {
            std::string uncompressedFile = tempDir+"/"+tempFilename+".tid."+std::to_string(i)+s;
            std::string compressedFile = uncompressedFile + "Compressed";
            // use lzma2 for .base to get better compression
            if (s == std::string(".base"))
                lzma2::lzma2_decompress(compressedFile.c_str(), uncompressedFile.c_str());
            else
                bsc::BSC_decompress(compressedFile.c_str(), uncompressedFile.c_str());
        }
    }
    auto bsc_end = std::chrono::high_resolution_clock::now();
    std::cout << "BSC/LZMA2 decompression done!\n";
    auto duration = 
        std::chrono::duration_cast<std::chrono::milliseconds>(bsc_end - bsc_start);
    std::cout << "Took " << duration.count()
              << " milliseconds" << std::endl;

    
    // during decompression, we first write the reads (as bitsets) to disk
    // and then we write the reads in order to output file in multiple passes
    // to limit memory consumption.
    std::vector<uint32_t> read_lengths(numReads);
    std::vector<std::vector<uint32_t>> read_ids(numEncodingThreads);
    std::vector<std::vector<uint32_t>> read_bitset_len(numEncodingThreads);

    std::cout << "Generating reads...\n";
    auto gen_start = std::chrono::high_resolution_clock::now();
    //loop through each thread
#pragma omp parallel for
    for (size_t i = 0; i < numEncodingThreads; ++i) {
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

        // temporary file to write bitset
        std::ofstream tmpFile(currentFilename, std::ios::binary);

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
                read_lengths[id] = currentRead.size();
                auto read_bitset = new DnaBitset(currentRead.c_str(), currentRead.size());
                read_bitset_len[i].push_back(read_bitset->to_file(tmpFile));
                read_ids[i].push_back(id);
                delete read_bitset;
            }
        }
        // now do the lone reads
        std::string loneRead;
        read_t id = 0, idInc;
        while(std::getline(loneFile, loneRead)) {
            // the id for the current read
            idFile.read((char*)&idInc, sizeof(read_t));
            id = id + idInc;
            read_lengths[id] = loneRead.size();
            auto read_bitset = new DnaBitset(loneRead.c_str(), loneRead.size());
            read_bitset_len[i].push_back(read_bitset->to_file(tmpFile));
            read_ids[i].push_back(id);
            delete read_bitset;
        }
        //close all files
        genomeFile.close();
        idFile.close();
        posFile.close();
        editTypeFile.close();
        editBaseFile.close();
        complementFile.close();
        loneFile.close();
        tmpFile.close();
    }
    std::cout << "Done!\n";
    auto gen_end = std::chrono::high_resolution_clock::now();
    duration = 
        std::chrono::duration_cast<std::chrono::milliseconds>(gen_end - gen_start);
    std::cout << "Took " << duration.count()
              << " milliseconds" << std::endl;

    //output the reads to a file
    std::ofstream outFile;
    outFile.open(outputFileName);
    std::cout << "Sorting and writing to file...\n";
    auto write_start = std::chrono::high_resolution_clock::now();

    // we first calculate the start and end points for the sort passes
    size_t num_bases_per_pass = (size_t)decompression_memory_gb*4*1000*1000*1000;
    std::vector<uint32_t> start_read, end_read;
    uint32_t cur_start = 0, cur_end = 0;
    while (true) {
        size_t cumul_bases = 0;
        for (; cur_end < numReads; cur_end++) {
            cumul_bases += read_lengths[cur_end];
            if (cumul_bases > num_bases_per_pass)
                break;
        }
        start_read.push_back(cur_start);
        end_read.push_back(cur_end);
        if (cur_end == numReads)
            break;
        cur_start = cur_end;
        cur_end = cur_start;
    }

    // now we write the reads to outfile in correct order
    // by picking in small passes over the temporary files
    // TODO: can we parallelize the to_string part somehow
    std::string currentReadToWrite;
    std::vector<DnaBitset*> reads_in_pass; 
    for (uint32_t i = 0; i < start_read.size(); i++) {
        // first pick the reads in this pass and put in reads_in_pass
        // we loop over the temporary file for each thread and all reads within
        // that file, picking anything withing start and end.
        reads_in_pass.resize(end_read[i]-start_read[i]);
        for (size_t j = 0; j < numEncodingThreads; j++) {
            std::string currentFilename = tempDir + tempFilename + ".tid." + std::to_string(j);
            std::ifstream fin(currentFilename, std::ios::binary);
            size_t cur_pos = 0;
            for (uint32_t k = 0; k < read_ids[j].size(); k++) {
                auto id = read_ids[j][k];
                if (id < end_read[i] && id >= start_read[i]) {
                    fin.seekg(cur_pos);
                    reads_in_pass[id-start_read[i]] = new DnaBitset(fin, read_lengths[id]);
                }
                cur_pos += read_bitset_len[j][k];
            }
        }

        // now we write these to disk
        for (uint32_t k = 0; k < end_read[i]-start_read[i]; k++) {
            reads_in_pass[k]->to_string(currentReadToWrite);
            delete reads_in_pass[k];
            outFile << currentReadToWrite << '\n';
        }
    }

    std::cout << "Done!\n";
    auto write_end = std::chrono::high_resolution_clock::now();
    duration = 
        std::chrono::duration_cast<std::chrono::milliseconds>(write_end - write_start);
    std::cout << "Sorting and writing to file took " << duration.count()
              << " milliseconds" << std::endl;

    boost::filesystem::remove_all(tempDir);
    auto overall_end = std::chrono::high_resolution_clock::now();
    duration = 
        std::chrono::duration_cast<std::chrono::milliseconds>(overall_end - overall_start);
    std::cout << "Total time for decompressioin " << duration.count()
              << " milliseconds" << std::endl;
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
    //handle the inserted bases in beginning
    uint32_t numInsStart;
    numInsStart = DirectoryUtils::read_var_uint32(posFile);
    for (size_t i = 0; i < numInsStart; i++){
        //write the inserted bases to read
        char editBase;
        editBaseFile.get(editBase);
        read.push_back(editBase);
    }
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
        if (editType == '\n'){
            break;
        }
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
    //handle the inserted bases in end
    uint32_t numInsEnd;
    numInsEnd = DirectoryUtils::read_var_uint32(posFile);
    for (size_t i = 0; i < numInsEnd; i++){
        //write the inserted bases to read
        char editBase;
        editBaseFile.get(editBase);
        read.push_back(editBase);
    }
    if (reverseComplement) {
        std::string temp;
        temp.swap(read);
        ReadData::toReverseComplement(temp.begin(), temp.end(),
                                      std::inserter(read, read.end()));
    }
}
