#ifndef EXPERIMENTS_READALIGNER_CUH
#define EXPERIMENTS_READALIGNER_CUH

#include "ReadFilter.h"
#include <string>
#include <vector>

class ReadAligner {
public:
    /***
     * Interface for roughly aligning two reads.
     * @param r1 First read to align
     * @param r2 Second read to align
     * @param relPos The rough position of the alignment. (Positive if r2 goes
     * after r1 and negative otherwise)
     * @return Whether alignment succeeded
     */
    virtual bool align(const std::string &r1, const std::string &r2,
                       ssize_t &relPos) = 0;
};

class MergeSortReadAligner : public ReadAligner {
public:
    const size_t k;
    const size_t kMerNumTh;

    /***
     * Returns true if the number of matching k-mers is greater than kMerNumTh.
     *
     * On simulate data, for overlaps of size 2000, error rate = 10%, k = 21 and
     * kMerNumTh = 10 works really well (No false negatives, false positives
     * capture overlaps of size >= 400)
     * @param k
     * @param kMerNumTh
     */
    MergeSortReadAligner(size_t k, size_t kMerNumTh);

    bool align(const std::string &r1, const std::string &r2, ssize_t &relPos);

private:
    void stringToSortedKMers(
        const std::string &s,
        std::vector<std::pair<MinHashReadFilter::kMer_t, size_t>> &v);
};

#endif // EXPERIMENTS_READALIGNER_CUH
