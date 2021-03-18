#include <cstring> /* std::memset */
#include <stdexcept> /* std::invalid_argument */
#include "dnaToBits.h"

uint8_t baseToInt(const char base) {
    // A->0, T->1, C->2, G->3
    return (base & 0b10) | ((base & 0b100) >> 2);
}

DnaBitset::DnaBitset(const char* dna_str, const size_t dna_len) {
	m_len = dna_len;

	/* number of bytes necessary to store dna_str as a bitset */
	size_t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);

	m_data = new uint8_t[dna_bytes];

	std::memset(m_data, 0, dna_bytes);

	/* for each base of the DNA sequence */
	for (size_t i = 0; i < m_len / 4; i++) {
        for (uint8_t j = 0; j < 4; j++) {
            m_data[i] <<= 2;
            m_data[i] |= baseToInt(dna_str[4*i+j]);
        }
    }
    if (m_len % 4 != 0) {
        int i = m_len / 4;
        for (uint8_t j = 0; j < m_len % 4; j++) {
            m_data[i] <<= 2;
            m_data[i] |= baseToInt(dna_str[4*i+j]);
        }
        for (uint8_t j = m_len % 4; j < 4; j++)
            m_data[i] <<= 2;
    }
}

DnaBitset::~DnaBitset() {
    delete[] m_data;
}

void DnaBitset::to_string(std::string &readStr) {
    const char int2dna[4] = {'A','T','C','G'}; // inverse of baseToInt
    const uint8_t masks[4] = {0x3<<6,0x3<<4,0x3<<2,0x3}; // 11 shifted
	readStr.resize(m_len);
	/* for each base of the DNA sequence */
    size_t pos = 0;
	for (size_t i = 0; i < m_len / 4; i++) {
        for (uint8_t j = 0; j < 4; j++) {
            readStr[pos++] = int2dna[(masks[j] & m_data[i])>>(6-2*j)];
        }
    }
    if (m_len % 4 != 0) {
        size_t i = m_len / 4;
        for (uint8_t j = 0; j < m_len % 4; j++) {
            readStr[pos++] = int2dna[(masks[j] & m_data[i])>>(6-2*j)];
        }
    }
}
