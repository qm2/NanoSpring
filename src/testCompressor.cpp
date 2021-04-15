#include "Compressor.h"
#include "ReadData.h"
#include "ReadFilter.h"
#include "StringAligner.h"
#include <ctime>
#include <iostream>
#include <omp.h>
#include <chrono> 
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
using namespace std::chrono; 
using namespace std;

int main(int argc, char **argv) {
    auto start = high_resolution_clock::now(); 
    omp_set_nested(1);
    std::srand(unsigned(std::time(0)));
    //program options
    namespace po = boost::program_options;
    bool help_flag = false, gzip_flag = false;
    std::string infile;
    int num_thr;
    //modifed the input parameters for the minimap2 version
    size_t k, n, overlapSketchThreshold, m_k, m_w, max_chain_iter;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", po::bool_switch(&help_flag), "produce help message")(
        "input-file,i",
        po::value<std::string>(&infile),
        "input file name")(
        "num-threads,t", po::value<int>(&num_thr)->default_value(16),
        "number of threads (default 16)")(
        "kmer,k", po::value<size_t>(&k)->default_value(23),
        "kmer size for the minhash (default 23)")(
        "num-hash,n", po::value<size_t>(&n)->default_value(64),
        "number of hash functions for minhash (default 64)")(
        "overlap-sketch-thr", po::value<size_t>(&overlapSketchThreshold)->default_value(8),
        "the overlap sketch threhold for minhash (default 8)")(
        "minimap-k", po::value<size_t>(&m_k)->default_value(20),
        "kmer size for the minimap2 (default 20)")(
        "minimap-w, w", po::value<size_t>(&m_w)->default_value(50),
        "window size for the minimap2 (default 50)")(
        "max-chain-iter", po::value<size_t>(&max_chain_iter)->default_value(400),
        "the max number of partial chains during chaining for minimap2 (default 400)")(
        "gzipped-fastq,g", po::bool_switch(&gzip_flag),
        "enable if compression input is gzipped fastq");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (help_flag) {
        std::cout << desc << "\n";
        return 0;
    }
    {
        // //original parameter for aligner
        // double maxErrorRate;
        // size_t editSlack;
        std::cout << "k n overlapSketchThreshold\nminimap:k minimap:w minimap:max_chain_iter"
                  << std::endl;
        if (k == 0)
            return 0;
        std::cout << "k n overlapSketchThreshold minimap:k minimap:w minimap:max_chain_iter\n"
                  << k << " " << n << " " << overlapSketchThreshold << " "
                  << m_k << " " << m_w << " " << max_chain_iter << std::endl;
        MergeSortReadAligner rA(21, 10);
        Compressor compressor;
        compressor.k = k;
        compressor.n = n;
        compressor.overlapSketchThreshold = overlapSketchThreshold;
        compressor.m_k = m_k;
        compressor.m_w = m_w;
        compressor.max_chain_iter = max_chain_iter;        
        compressor.rA = &rA;
        const std::string filename(infile);
        const std::string extension =
            filename.substr(filename.find_last_of('.') + 1);
        if (gzip_flag)
            compressor.filetype = ReadData::Filetype::GZIP;
        else if (!extension.compare("fastq"))
            compressor.filetype = ReadData::Filetype::FASTQ;
        compressor.compress(infile.c_str(), num_thr);

        auto stop = high_resolution_clock::now(); 
        auto duration = duration_cast<seconds>(stop - start); 
        cout << "Time taken by function: "<< duration.count() << " seconds" << endl; 
    }
}
