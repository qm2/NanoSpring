#ifndef C9E093FF_F3FD_4F67_9206_519C995B0ED1
#define C9E093FF_F3FD_4F67_9206_519C995B0ED1

namespace bsc {

const int BSC_BLOCK_SIZE = 48; // 48 MB

void BSC_compress(const char *infile, const char *outfile,
                  const int bsize = BSC_BLOCK_SIZE);

void BSC_decompress(const char *infile, const char *outfile);

} // namespace bsc

#endif /* C9E093FF_F3FD_4F67_9206_519C995B0ED1 */
