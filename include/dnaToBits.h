#ifndef __dnatobits_h__
#define __dnatobits_h__

#include <string>

#define BASE_MASK 0x3 /* binary: 11 */

/* useful constants */
enum
{
    BASE_A = 0x0, /* binary: 00 */
    BASE_C = 0x1, /* binary: 01 */
    BASE_G = 0x2, /* binary: 10 */
    BASE_T = 0x3, /* binary: 11 */
};

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
