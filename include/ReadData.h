#ifndef C3C2DD2A_9114_4E1C_B89E_18F31D007DCB
#define C3C2DD2A_9114_4E1C_B89E_18F31D007DCB

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
     * @brief Loads the read data from a file
     *
     * @param fileName
     */
    void loadFromFile(const char *fileName);

    /**
     * @brief Get the number of reads
     *
     * @return size_t
     */
    size_t getNumReads();

    /**
     * @brief Returns the read with readId as a string
     *
     * @param readId
     * @return std::string
     */
    std::string &getRead(size_t readId);

    /// TODO: don't expose these three funcitons. change the interface to hide
    /// the raw data.
    std::vector<unsigned long> &getReadPos();

    std::vector<unsigned long> &getReadPosSorted();

    std::vector<std::unique_ptr<std::string>> &getReadData();

private:
    size_t numReads;
    std::vector<std::unique_ptr<std::string>> readData;

    // The following are for testing
    std::vector<unsigned long> readPos;
    std::vector<unsigned long> readPosSorted;
    std::vector<std::unique_ptr<std::string>> editStrings;
    size_t readLen;
    // End testing
};

#endif /* C3C2DD2A_9114_4E1C_B89E_18F31D007DCB */