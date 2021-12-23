#include "Compressor.h"
#include "Decompressor.h"
#include "ReadData.h"
#include "ReadFilter.h"
#include "StringAligner.h"
#include "DirectoryUtils.h"
#include <csignal>
#include <ctime>
#include <iostream>
#include <omp.h>
#include <chrono> 
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
using namespace std::chrono; 
using namespace std;

std::string temp_dir_global;  // for interrupt handling
bool temp_dir_flag_global = false;

void signalHandler(int signum) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";
  std::cout << "Program terminated unexpectedly\n";
  if (temp_dir_flag_global) {
    std::cout << "Deleting temporary directory:" << temp_dir_global << "\n";
    boost::filesystem::remove_all(temp_dir_global);
  }
  exit(signum);
}


int main(int argc, char **argv) {
    // register signal SIGINT and signal handler
    signal(SIGINT, signalHandler);
    auto start = high_resolution_clock::now(); 
    omp_set_nested(1);
    std::srand(unsigned(std::time(0)));
    //program options
    namespace po = boost::program_options;
    bool help_flag = false, compress_flag = false, decompress_flag = false;
    bool low_mem = true; // fix to low mem mode that uses temporary file to store read bitsets
    std::string infile, outfile;
    int num_thr, decompression_memory_gb;
	std::string working_dir;
    //modifed the input parameters for the minimap2 version
    size_t k, n, overlapSketchThreshold, m_k, m_w, max_chain_iter, edge_threshold;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", po::bool_switch(&help_flag), "produce help message")(
        "compress,c", po::bool_switch(&compress_flag), "compress")(
        "decompress,d", po::bool_switch(&decompress_flag), "decompress")(
        "input-file,i",
        po::value<std::string>(&infile),
        "input file name")(
        "output-file,o",po::value<std::string>(&outfile),
        "output file name")(
        "num-threads,t", po::value<int>(&num_thr)->default_value(20),
        "number of threads (default 20)")(
        "kmer,k", po::value<size_t>(&k)->default_value(23),
        "kmer size for the minhash (default 23)")(
        "num-hash,n", po::value<size_t>(&n)->default_value(60),
        "number of hash functions for minhash (default 60)")(
        "overlap-sketch-thr", po::value<size_t>(&overlapSketchThreshold)->default_value(6),
        "the overlap sketch threshold for minhash (default 6)")(
        "minimap-k", po::value<size_t>(&m_k)->default_value(20),
        "kmer size for the minimap2 (default 20)")(
        "minimap-w, w", po::value<size_t>(&m_w)->default_value(50),
        "window size for the minimap2 (default 50)")(
        "max-chain-iter", po::value<size_t>(&max_chain_iter)->default_value(400),
        "the max number of partial chains during chaining for minimap2 (default 400)")(
        "edge-thr", po::value<size_t>(&edge_threshold)->default_value(4000000),
        "the max number of edges allowed in a consensus graph (default 4000000)")(
        "working-dir,w", po::value<std::string>(&working_dir)->default_value("."),
        "directory to create temporary files (default current directory)")(
        "decompression-memory", po::value<int>(&decompression_memory_gb)->default_value(5),
        "attempt to set peak memory usage for decompression in GB (default 5 GB) by "
        "using disk-based sort for writing reads in the correct order. This is only "
        "approximate and might have no effect at very low settings or with large "
        "number of threads when another decompressor stage is the biggest memory "
        "contributor. Very low values might lead to slight reduction in speed.");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (help_flag) {
        std::cout << desc << "\n";
        return 0;
    }
    if(infile.empty()){
    	std::cout<<"No input file specified"<<"\n";
    	std::cout << desc << "\n";
    	return 0;
    } 
    if(outfile.empty()){
    	std::cout<<"No output file specified"<<"\n";
		std::cout << desc << "\n";
		return 0;
    }

    if ((!compress_flag && !decompress_flag) ||
      (compress_flag && decompress_flag)) {
    std::cout
        << "Exactly one of compress or decompress needs to be specified \n";
    std::cout << desc << "\n";
    return 1;
    }
    // generate randomly named temporary directory in the working directory
    std::string temp_dir;
    while (true) {
      std::string random_str = "tmp." + DirectoryUtils::random_string(10);
      temp_dir = working_dir + "/" + random_str + "/";
      if (!boost::filesystem::exists(temp_dir)) break;
    }
    if (!boost::filesystem::create_directory(temp_dir)) {
      throw std::runtime_error("Cannot create temporary directory.");
    }
    std::cout << "Temporary directory: " << temp_dir << "\n";

    temp_dir_global = temp_dir;
    temp_dir_flag_global = true;
    try{
        if (compress_flag){
            if (k == 0)
                throw std::runtime_error("Invalid k");
            std::cout << "k n overlapSketchThreshold minimap:k minimap:w minimap:max_chain_iter edge_threshold\n"
                      << k << " " << n << " " << overlapSketchThreshold << " "
                      << m_k << " " << m_w << " " << max_chain_iter <<" " << edge_threshold<< std::endl;
            MergeSortReadAligner rA(21, 10);
            Compressor compressor;
            compressor.k = k;
            compressor.n = n;
            compressor.overlapSketchThreshold = overlapSketchThreshold;
            compressor.m_k = m_k;
            compressor.m_w = m_w;
            compressor.max_chain_iter = max_chain_iter;  
            compressor.edge_threshold = edge_threshold;       
            compressor.rA = &rA;
            compressor.tempDir = temp_dir;
            compressor.outputFileName = outfile;
            compressor.low_mem = low_mem;
            const std::string filename(infile);
            const std::string extension =
                filename.substr(filename.find_last_of('.') + 1);
            if (!extension.compare("gz"))
                compressor.filetype = ReadData::Filetype::GZIP;
            else if (!extension.compare("fastq"))
                compressor.filetype = ReadData::Filetype::FASTQ;
            compressor.compress(infile.c_str(), num_thr);

            auto stop = high_resolution_clock::now(); 
            auto duration = duration_cast<seconds>(stop - start); 
            cout << "Time taken by function: "<< duration.count() << " seconds" << endl; 
        }
        else{
            Decompressor dc;
            dc.tempDir = temp_dir;
            if (decompression_memory_gb < 1) {
                throw std::runtime_error("Invalid decompression-memory parameter: must be >= 1.");
            }
            dc.decompress(infile.c_str(), outfile.c_str(), num_thr, decompression_memory_gb);
        }
    }
    // Error handling
    catch (std::runtime_error& e) {
        std::cout << "Program terminated unexpectedly with error: " << e.what()
                  << "\n";
        std::cout << "Deleting temporary directory...\n";
        boost::filesystem::remove_all(temp_dir);
        temp_dir_flag_global = false;
        std::cout << desc << "\n";
        return 1;
    } catch (...) {
        std::cout << "Program terminated unexpectedly\n";
        std::cout << "Deleting temporary directory...\n";
        boost::filesystem::remove_all(temp_dir);
        temp_dir_flag_global = false;
        std::cout << desc << "\n";
        return 1;
    }
    boost::filesystem::remove_all(temp_dir);
    temp_dir_flag_global = false;
    return 0;
}
