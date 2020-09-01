#ifndef D5822833_5BFC_4986_B7F9_20EEA68B3568
#define D5822833_5BFC_4986_B7F9_20EEA68B3568

#include "Edits.h"
#include "ReadAligner.h"
#include <string>

class Compressor {
public:
    // Parameters for filtering
    size_t k, n, overlapSketchThreshold;
    ReadAligner *rA;
    // Parameters for consensus
    StringAligner *aligner;
    // The temp directories
    std::string tempDir = "tempRaw/";
    std::string compressedTempDir = "tempCompressed/";
    // The output filenames to use in temp directories
    std::string tempFileName = "Contig";

    void compress(const char *inputFileName) const;
};

#endif /* D5822833_5BFC_4986_B7F9_20EEA68B3568 */