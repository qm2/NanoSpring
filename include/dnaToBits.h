#ifndef __dnatobits_h__
#define __dnatobits_h__

#include <string>

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
     * @brief Destructor.
     */
    ~DnaBitset();

    /**
     * @brief Returns the stored DNA sequence as an ASCII string.
     */
    void to_string(std::string &readStr);

private:
    uint8_t* m_data;
    size_t m_len;
};

#endif /* __dntobits_h__ */
