#ifndef C7EF1708_1404_4E39_A140_F2F9914B7F6A
#define C7EF1708_1404_4E39_A140_F2F9914B7F6A

#include "Edits.h"
#include "ReadAligner.h"
#include <string>

class Decompressor {
public:
    // The temp directories. Must have trailing "/"
    std::string tempDir = "tempRaw/";
    std::string compressedTempDir = "tempCompressed/";

    // The output filenames to use in temp directories
    std::string tempFilename = "Contig";

    void decompress(const char *inputFileName, const char *outputFileName);

private:
    // Number of contigs and number of reads; initialized by bscDecompress
    size_t numContigs;
    size_t numReads;

    /**
     * @brief Extract tar archive into compressedTempDir
     *
     * @param inputFileName
     */
    void untar(const char *inputFileName) const;

    /**
     * @brief  We clear the temp directories and create them if they do not
     * exist
     *
     */
    void prepareTempDirs() const;

    /**
     * @brief Decompressed the *Compressed files generated using bsc and
     * initializes numContigs and numReads
     *
     */
    void bscDecompress();

    /**
     * @brief Decompressed the files associated with a specific contig
     *
     * @param f the name of that contig
     */
    void bscDecompress(const std::string &f) const;

    /**
     * @brief This function generates the outputfile from the files in tempDir
     *
     * @param outputFileName
     */
    void generateReads(const char *outputFileName) const;

    /**
     * @brief Stores the reads in contig contigId into reads
     *
     * @param reads
     * @param contigId
     */
    void generateReads(std::string *reads, size_t contigId) const;

    void generateRead(const std::string &genome, std::string &read,
                      std::ifstream &posFile, std::ifstream &editTypeFile,
                      std::ifstream &editBaseFile) const;
};

#endif /* C7EF1708_1404_4E39_A140_F2F9914B7F6A */