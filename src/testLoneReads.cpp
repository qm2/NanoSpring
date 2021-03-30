#include "Decompressor.h"
// #include "ConsensusGraph.h"
// #include "Consensus.h"
#include "DirectoryUtils.h"
#include "ReadData.h"
#include "bsc_helper.h"
#include "dnaToBits.h"
#include "minimap.h"
#include <cstdlib>
#include <boost/filesystem.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <omp.h>
#include <vector>
#include <string>

int main(int argc, char **argv) {
    omp_set_nested(1);

    if (argc < 2) {
        std::cout << "Usage ./testLoneReads inputFilename []"
                  << std::endl;
        return 1;
    }
    int numThreads = omp_get_max_threads();
    if (argc == 3)
        numThreads = std::atoi(argv[2]);
    std::cout << "numDecodingThreads: " << numThreads << "\n";
    omp_set_num_threads(numThreads);
    
    std::string tempDir = "tempRaw/";
    std::string tempFilename = "Stream";
    size_t numEncodingThreads;
    size_t numContigs;
    read_t numReads;

    DirectoryUtils::clearDir(tempDir);

    // Untar
    {
        std::cout << "Extracting tar archive ...";
        std::string tar_command =
            "tar -C " + tempDir + " -xf " + argv[1];
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
#pragma omp parallel for
    for (size_t i = 0; i < numEncodingThreads; ++i){
        for (auto &s : suffix) {
            std::string uncompressedFile = tempDir+"/"+tempFilename+".tid."+std::to_string(i)+s;
            std::string compressedFile = uncompressedFile + "Compressed";
            bsc::BSC_decompress(compressedFile.c_str(), uncompressedFile.c_str());
        }
    }
    std::cout << "BSC decompression done!\n";

    // create a vector to temporarily hold the strings
    //combine all reference genomes
    std::vector<std::string> genome_vector;
    for (size_t i = 0; i < numEncodingThreads; ++i){
        std::string currentFilename = tempDir + tempFilename + ".tid." + std::to_string(i);
        //open all the files for this thread
        std::ifstream genomeFile;
        genomeFile.open(currentFilename + ".genome");
        std::string genome;
        while(std::getline(genomeFile, genome)){
            genome_vector.push_back(genome);
        }
    }
    //convert the string vector to c-style char**
    const char **genome_cstrings = new const char*[genome_vector.size()];
    for (int i = 0; i < genome_vector.size(); ++i){
        genome_cstrings[i] = genome_vector[i].c_str();
    }



    //create minimap index for all reference genomes
    //initialize the local buffer
    mm_tbuf_t *b = mm_tbuf_init();
    //initialize the mapopt and iopt
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    //0 correpons to map-ont
    mm_set_opt(0, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR;
    // mopt.flag |= MM_F_FOR_ONLY; 
    //call the mm_idx_str to return the index for the reference read    
    mm_idx_t * idx = mm_idx_str(10, 15, false, 14, 1, genome_cstrings, NULL);
    mm_mapopt_update(&mopt, idx);
    //use the index to align with the current read
    //we only want forward matches: use rev in mm_reg1_t
    // create a vector to temporarily hold the strings
    //combine all reference genomes
    int total_lone_reads = 0;
    int recovered_lone_reads = 0;
    size_t total_aligned_length = 0;
    size_t total_edit_distance = 0;
    for (size_t i = 0; i < numEncodingThreads; ++i){
        std::string currentFilename = tempDir + tempFilename + ".tid." + std::to_string(i);
        //open all the files for this thread
        std::ifstream loneFile;
        loneFile.open(currentFilename + ".lone");
        std::string lone;
        int hits;
        while(std::getline(loneFile, lone)){
            total_lone_reads++;
            mm_reg1_t* reg = mm_map(idx, lone.length(), lone.c_str(), &hits, b, &mopt, NULL); 
            if(hits > 0){ 
                recovered_lone_reads++;
                //if we have multiple hits, just stick with first hit
                mm_reg1_t *r = &reg[0];
                assert(r->p); 
                //qpos is the current position on the query read
                //rpos is the current position on the reference read
                int alignedLen;
                int editDis; 
                //calculate the edit distance
                editDis = r->blen - r->mlen + r->p->n_ambi;
                //calculate the aligned length; notice that I use the aligned length for the query read here
                alignedLen = r->qe - r->qs; 
                total_aligned_length += alignedLen; 
                total_edit_distance += editDis; 
                // std::cout<<"aligned length: "<<alignedLen<<" edit distance: "<<editDis<<std::endl;
                free(r->p);
                if (hits > 1) {
                    // cleanup
                    for (int i = 1; i < hits; i++)
                        free(reg[i].p);
                }
                free(reg);  
            }
            else{
                    free(reg);  
            }
        }
    }
    std::cout<<"total lone reads: "<<total_lone_reads<<" recovered lone reads: "<<recovered_lone_reads<<std::endl; 
    std::cout<<"total aligned length: "<<total_aligned_length<<std::endl;
    std::cout<<"average aligned length: "<<(double)total_aligned_length/recovered_lone_reads<<std::endl;     
    std::cout<<"average edit distance: "<<(double)total_edit_distance/recovered_lone_reads<<std::endl;   
    mm_tbuf_destroy(b);
    mm_idx_destroy(idx);
    delete [] genome_cstrings;
    boost::filesystem::remove_all(tempDir);
    return 0;
}