#ifndef LZMA2_HELPER_H_
#define LZMA2_HELPER_H_

namespace lzma2 {

const int PRESET = 6; // preset between 1-10

void lzma2_compress(const char *infile, const char *outfile,
                  const int preset = PRESET);

void lzma2_decompress(const char *infile, const char *outfile);

} // namespace lzma2

#endif /* LZMA2_HELPER_H_ */
