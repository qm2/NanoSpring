#ifndef EXPERIMENTS_TESTALIGNER_CUH
#define EXPERIMENTS_TESTALIGNER_CUH

#include "NanoporeReads.cuh"

class TestAligner {
public:
    TestAligner(const char *fileName);

    void test(const size_t k, const size_t kMerNumTh, const size_t baseOverlapTh);

private:
    NanoporeReads nR;
};

#endif //EXPERIMENTS_TESTALIGNER_CUH
