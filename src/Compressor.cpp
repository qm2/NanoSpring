#include "Compressor.h"
#include "Consensus.h"
#include "Contig.h"
#include <boost/filesystem.hpp>
#include <chrono>
// #include <filesystem>
#include <fstream>

void Compressor::compress(const char *inputFileName) const {

    //        MergeSortReadAligner rA(10, 1);
    NanoporeReads nR(inputFileName, k, n);
    MinHashReadFilter rF(overlapSketchThreshold, nR);
    {
        auto start = std::chrono::high_resolution_clock::now();
        nR.calculateMinHashSketches();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Calculating MinHashes took " << duration.count()
                  << " milliseconds" << std::endl;
    }
    ContigGenerator cG(rA, nR, &rF);
    {
        auto start = std::chrono::high_resolution_clock::now();
        cG.generateContigs();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Generating contigs took " << duration.count()
                  << " milliseconds" << std::endl;
    }
    std::cout << cG << std::endl;

    // We clear the temp directories and create them if they do not
    // exist
    boost::system::error_code ec;
    const boost::filesystem::path tempDirPath(tempDir);
    const boost::filesystem::path compressedTempDirPath(compressedTempDir);
    boost::filesystem::remove_all(tempDirPath, ec);
    boost::filesystem::remove_all(compressedTempDirPath, ec);
    boost::filesystem::create_directory(tempDirPath, ec);
    boost::filesystem::create_directory(compressedTempDirPath, ec);

    size_t numContigs = cG.contigs.size();
    std::vector<Contig *> contigs(cG.contigs.begin(), cG.contigs.end());
#pragma omp parallel for
    for (size_t i = 0; i < numContigs; ++i) {
        Contig *c = contigs[i];
        std::set<std::pair<long, read_t>> &readsInContig = c->reads;
        ConsensusGraph consensusGraph(aligner);
        consensusGraph.tempDir = tempDir;
        consensusGraph.compressedTempDir = compressedTempDir;
        consensusGraph.addReads(readsInContig, cG.nR.readData);

        consensusGraph.calculateMainPathGreedy();
        consensusGraph.printStatus();
        std::string filename = tempFilename + std::to_string(i);
        consensusGraph.writeMainPath(filename);
        consensusGraph.writeReads(filename);
    }

    std::ofstream metaData;
    metaData.open(compressedTempDir + "metaData");
    metaData << cG.contigs.size() << '\n';
    metaData << cG.nR.readData.size() << '\n';
    metaData.close();

    std::cout << "Creating tar archive ..." << std::endl;
    std::string outfile = "compressedFile";
    std::string tar_command =
        "tar -cvf " + outfile + " -C " + compressedTempDir + " . ";
    int tar_status = std::system(tar_command.c_str());
    if (tar_status)
        throw std::runtime_error(
            "Error occurred during tar archive generation.");
    std::cout << "Tar archive done!\n";
    std::string lsCommand = "ls -lh " + outfile;
    int ls_status = std::system(lsCommand.c_str());
    if (ls_status)
        throw std::runtime_error("Error occurred during ls command.");
}