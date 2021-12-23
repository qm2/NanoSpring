#ifndef D5822833_5BFC_4986_B7F9_20EEA68B3568
#define D5822833_5BFC_4986_B7F9_20EEA68B3568

#include "Consensus.h"
#include "Edits.h"
#include "ReadAligner.h"
#include "StringAligner.h"
#include <string>

class Compressor {
public:
    /** Parameters for filtering **/
    size_t k, n, overlapSketchThreshold, m_k, m_w, max_chain_iter, edge_threshold;
    ReadData::Filetype filetype = ReadData::Filetype::READ;
    ReadAligner *rA;
    std::string tempDir;
    /** The temp directories **/
    // std::string tempDir = "tempRaw/";
    /** The output filenames to use in temp directories **/
    std::string tempFileName = "Stream";

    std::string outputFileName = "compressedFile";

    bool low_mem;
    // std::string tarFileName = "originalFile";

    /**
     * @brief Compresses the read data file inputFileName and stores the result
     * in this->outputFileName
     *
     * @param inputFileName
     */
    void compress(const char *inputFileName, const int numThr) const;

};

#endif /* D5822833_5BFC_4986_B7F9_20EEA68B3568 */
