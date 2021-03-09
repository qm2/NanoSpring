#ifndef __dna_h__
#define __dna_h__

#include <cstring>   /* std::memset */
#include <stdexcept> /* std::invalid_argument */

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
    DnaBitset(const char* dna_str, const size_t dna_len)
    {
        m_len = dna_len;

        /* number of bytes necessary to store dna_str as a bitset */
        size_t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);

        m_data = new uint8_t[dna_bytes];

        std::memset(m_data, 0, dna_bytes);

        /* for each base of the DNA sequence */
        for (size_t i = 0; i < dna_len; ++i)
        {
            uint8_t shift = 6 - 2 * (i % 4);

            switch (dna_str[i])
            {
                case 'A':
                    m_data[i / 4] |= BASE_A << shift;
                    break;
                case 'C':
                    m_data[i / 4] |= BASE_C << shift;
                    break;
                case 'G':
                    m_data[i / 4] |= BASE_G << shift;
                    break;
                case 'T':
                    m_data[i / 4] |= BASE_T << shift;
                    break;
                default:
                    throw std::invalid_argument("invalid DNA base");
            }

            shift = (shift == 0) ? 6 : shift - 2;
        }
    }

    /**
     * @brief Destructor.
     */
    ~DnaBitset()
    {
        delete[] m_data;
    }

    /**
     * @brief Returns the stored DNA sequence as an ASCII string.
     */
    void to_string(std::string &readStr)
    {
        readStr.resize(m_len);
        /* for each base of the DNA sequence */
        for (size_t i = 0; i < m_len; ++i)
        {
            uint8_t shift = 6 - 2 * (i % 4);
            uint8_t mask = BASE_MASK << shift;

            /* get the i-th DNA base */
            uint8_t base = (m_data[i / 4] & mask) >> shift;

            switch (base)
            {
                case BASE_A:
                    readStr[i] = 'A';
                    break;
                case BASE_C:
                    readStr[i] = 'C';
                    break;
                case BASE_G:
                    readStr[i] = 'G';
                    break;
                case BASE_T:
                    readStr[i] = 'T';
                    break;
                default:
                    throw std::runtime_error("invalid DNA base");
            }
        }
    }

private:
    uint8_t* m_data;
    size_t m_len;
};

#endif /* __dna_h__ */