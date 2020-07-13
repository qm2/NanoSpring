//
// Created by yifan on 7/10/20.
//

#ifndef EXPERIMENTS_CONTIG_CUH
#define EXPERIMENTS_CONTIG_CUH

#include <set>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <omp.h>
#include "ReadAligner.cuh"

// We use the indices of the reads to represent them
typedef size_t read_t;

/***
 * A contig is just an ordered set of reads along with their positions
 */
class Contig {
public:
    // The first element means position, the second means read.
    std::set<std::pair<long, read_t>> reads;
};

class ContigGenerator {
public:
    ContigGenerator(ReadAligner *rA, NanoporeReads &nR, ReadFilter *rF);

    void generateContigs();

    friend std::ostream &operator<<(std::ostream &out, const ContigGenerator &o);

private:
    /***
     * This is the read aligner used to align the reads as they are being added to the Contig
     */
    ReadAligner *rA;
    NanoporeReads &nR;
    ReadFilter *rF;
    std::vector<Contig> contigs;
};

#endif //EXPERIMENTS_CONTIG_CUH
