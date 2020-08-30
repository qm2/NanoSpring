#ifndef EXPERIMENTS_TESTALIGNER_CUH
#define EXPERIMENTS_TESTALIGNER_CUH

#include "ReadData.h"

/**
 * @brief Tests aligner algorithm used for contig generation.
 *
 */
class TestAligner {
public:
    TestAligner(const char *fileName);

    void test(const size_t k, const size_t kMerNumTh,
              const size_t baseOverlapTh);

private:
    ReadData rD;
};

#endif // EXPERIMENTS_TESTALIGNER_CUH
