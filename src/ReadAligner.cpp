#include "ReadAligner.h"
#include <algorithm>

MergeSortReadAligner::MergeSortReadAligner(size_t k, size_t kMerNumTh)
    : k(k), kMerNumTh(kMerNumTh) {}

bool MergeSortReadAligner::align(const std::string &r1, const std::string &r2,
                                 ssize_t &relPos) {
    std::vector<std::pair<kMer_t, size_t>> v1, v2;
    stringToSortedKMers(r1, v1);
    stringToSortedKMers(r2, v2);
    auto begin1 = v1.begin();
    auto end1 = v1.end();
    auto k2 = v2.begin();
    auto end2 = v2.end();
    size_t numMatching = 0;
    double sum = 0;
    for (auto k1 = begin1; k1 < end1;) {
        if (k2 >= end2)
            break;
        kMer_t kMer1 = k1->first;
        kMer_t kMer2 = k2->first;
        if (kMer1 < kMer2)
            k1 = std::upper_bound(k1, end1, *k1);
        else if (kMer1 > kMer2)
            k2 = std::upper_bound(k2, end2, *k2);
        else {
            numMatching++;
            sum += (long)k1->second - (long)k2->second;
            k1 = std::upper_bound(k1, end1, *k1);
            k2 = std::upper_bound(k2, end2, *k2);
        }
    }
    if (numMatching > 0)
        relPos = sum / numMatching;
    return numMatching >= kMerNumTh;
}

void MergeSortReadAligner::stringToSortedKMers(
    const std::string &s,
    std::vector<std::pair<kMer_t, size_t>> &v) {
    ssize_t maxI = s.length() - k + 1;
    if (maxI <= 0)
        return;
    v.resize(maxI);
    std::vector<kMer_t> kmers(maxI);
    MinHashReadFilter::string2KMers(s, k, kmers);
    for (size_t i = 0; i < (size_t)maxI; ++i) {
        v[i] = (std::make_pair(kmers[i], i));
    }
    std::sort(v.begin(), v.end());
}
