#ifndef C3C2DD2A_9114_4E1C_B89E_18F31D007DCB
#define C3C2DD2A_9114_4E1C_B89E_18F31D007DCB

#include "Types.h"

#include <memory>
#include <string>
#include <vector>

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
    enum Filetype { FASTQ, READ };

    /** The average length of the reads **/
    size_t avgReadLen;

    /**
     * @brief Loads the read data from a file (.reads format)
     *
     * @param fileName
     */
    void loadFromFile(const char *fileName, enum Filetype filetype = READ);

    /**
     * @brief Get the number of reads
     *
     * @return read_t
     */
    read_t getNumReads();

    /**
     * @brief Returns the read with readId as a string
     *
     * @param readId
     * @return std::string
     */
    std::string &getRead(read_t readId);

    /// TODO: don't expose these three funcitons. change the interface to hide
    /// the raw data.
    std::vector<unsigned long> &getReadPos();

    std::vector<std::unique_ptr<std::string>> &getReadData();

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

private:
    read_t numReads;
    std::vector<std::unique_ptr<std::string>> readData;

    // The following are for testing
    std::vector<unsigned long> readPos;
    std::vector<std::unique_ptr<std::string>> editStrings;
    /** Whether the read is a reverse complement **/
    std::vector<bool> reverse;
    // End testing

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
     * @param fileName
     */
    void loadFromFastqFile(const char *fileName);
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
