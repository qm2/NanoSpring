#ifndef __dnatobits_h__
#define __dnatobits_h__

#include <string>
#include <fstream>

class DnaBitset
{
public:
    /**
     * @brief Constructs a compressed representation of a DNA sequence.
     * @param dna_str A string holding a DNA sequence (e.g. "ATGCACG").
     * @param dna_len The length of the DNA sequence.
     */
    DnaBitset(const char* dna_str, const size_t dna_len);

    /**
     * @brief Constructs a DnaBitset from a file with the bitset representation.
     * @param fin ifstream from where to read the binary m_data.
     * @param dna_len The length of the DNA sequence.
     */
    DnaBitset(std::ifstream &fin, const size_t dna_len);

    /**
     * @brief Destructor.
     */
    ~DnaBitset();

    /**
     * @brief Returns the stored DNA sequence as an ASCII string.
     */
    void to_string(std::string &readStr);

    /**
     * @brief Writes the stored DNA bitset to a binary file.
     * @return Number of bytes written
     */    
    size_t to_file(std::ofstream &fout);

private:
    uint8_t* m_data;
    size_t m_len;
    size_t dna_bytes;
};

#endif /* __dntobits_h__ */
