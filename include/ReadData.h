#ifndef C3C2DD2A_9114_4E1C_B89E_18F31D007DCB
#define C3C2DD2A_9114_4E1C_B89E_18F31D007DCB

#include "Types.h"
#include "dnaToBits.h"
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <mutex>

/**
 * @brief This class is responsible for reading the reads from a file, storing
 * them, and providing them to other users
 *
 */
class ReadData {
public:
    /**
     * @brief The filetype to read
     *
     */
    enum Filetype { FASTQ, READ, GZIP};

    /** The average length of the reads **/
    size_t avgReadLen;

    /** The length of the longest read **/
    size_t maxReadLen;
    
    std::string tempDir;



    /**
     * @brief Loads the read data from a file (.reads format)
     *
     * low_mem decides whether we use low memory mode where reads are kept in a temporary file
     * @param fileName
     */
    void loadFromFile(const char *fileName, enum Filetype filetype = FASTQ, bool low_mem = false);

    /**
     * @brief Get the number of reads
     *
     * @return read_t
     */
    read_t getNumReads();

    /**
     * @brief Returns the read with readId as a string
     *
     * @param readId, readStr
     * @return none
     */
    void getRead(read_t readId, std::string &readStr);

    /// TODO: don't expose these three funcitons. change the interface to hide
    /// the raw data.
    std::vector<unsigned long> &getReadPos();

    // std::vector<std::unique_ptr<std::string>> &getReadData();
    std::vector<std::unique_ptr<DnaBitset>> &getReadData();
    /**
     * @brief Turns a base to its complement base
     *
     * A <-> T
     * C <-> G
     *
     * @param base
     * @return char
     */
    static char toComplement(char base);

    /**
     * @brief Turns a DNA strand to its reverse complement.
     *
     * @tparam Iterator Bidirectional random access iterator that dereferences
     * to char; in particular, dereferences to one of 'A', 'T', 'C', and 'G'
     * @tparam Inserter Insert iterator used to store the reverse complement
     * DNA. Should dereference to char.
     * @param originalBegin Beginning of original string
     * @param originalEnd End of original string (one past end)
     * @param reverseComplement Insert iterator used to store the reverse
     * complement string
     */
    template <typename Iterator, typename Inserter>
    static void toReverseComplement(Iterator originalBegin,
                                    Iterator originalEnd,
                                    Inserter reverseComplement);

    /* Destructor to delete bitset file when applicable */
    ~ReadData();

private:
    read_t numReads;
    // std::vector<std::unique_ptr<std::string>> readData;
    /* read data (when reads_in_memory is true) */
    std::vector<std::unique_ptr<DnaBitset>> readData;

    // The following are for testing
    std::vector<unsigned long> readPos;
    std::vector<std::unique_ptr<std::string>> editStrings;
    /** Whether the read is a reverse complement **/
    std::vector<bool> reverse;
    // End testing

    /* whether reads are in memory or in temporary file on disk */
    bool reads_in_memory;
    
    /* if reads_in_memory false, bitset used for loading data from temporary file */
    DnaBitset read_bitset;
    
    /* if reads_in_memory false, store the file pointer to the temporary file */
    std::unique_ptr<std::ifstream> fin_bitset;

    /* if reads_in_memory false, mutex to enable exclusive access to fin and read_bitset */
    std::mutex fin_bitset_mtx;

    /* if reads_in_memory false, store the position in temporary file for each read */
    std::vector<size_t> read_pos_in_file;

    /* if reads_in_memory false, store the read length for each read */
    std::vector<size_t> read_lengths;

    /* if reads_in_memory false, file where we put temporary bitsets */
    std::string bitsetFileName = "readBitset";
    /**
     * @brief Loads the read data from a file (.reads format)
     *
     * Assumes the file has format
     * [nc]:readPos:editString
     * read
     *
     * [nc] stands for [n]ormal and reverse [c]omplement
     *
     * @param fileName
     */
    void loadFromReadFile(const char *fileName);

    /**
     * @brief Loads the read data from a file (.reads format)
     *
     * For every four lines, the second line is the read. The rest we disgard
     *
     * If the gzip_flag is set to True, then the input file is in gzip version 
     *
     * low_mem flag sets whether we want to store reads in RAM or in temporary file
     * and is used to call either loadFromFastqFile_highmem or loadFromFastqFile_lowmem
     *
     * @param fileName
     */
    void loadFromFastqFile(const char *fileName, bool gzip_flag, bool low_mem);
    
    void loadFromFastqFile_lowmem(const char *fileName, bool gzip_flag);
    void loadFromFastqFile_highmem(const char *fileName, bool gzip_flag);
};

/******************************************************************************/
/* Implementation of public template functions */
template <typename Iterator, typename Inserter>

void ReadData::toReverseComplement(Iterator originalBegin, Iterator originalEnd,
                              Inserter reverseComplement) {
    if (originalBegin == originalEnd)
        return;
    Iterator it = originalEnd;
    do {
        --it;
        *(reverseComplement++) = toComplement(*it);
    } while (it != originalBegin);
}

#endif /* C3C2DD2A_9114_4E1C_B89E_18F31D007DCB */
