#ifndef EXPERIMENTS_CONTIG_CUH
#define EXPERIMENTS_CONTIG_CUH

#include "ReadAligner.h"
#include "ReadData.h"
#include "ReadFilter.h"
#include "Types.h"
#include <set>
#include <unordered_map>
#include <unordered_set>

/**
 * The content of this file is no longer used
 */
/// TODO: Remove this file from the project

/**
 * A contig is just an ordered set of reads along with their positions
 */
class Contig {
public:
    // The first element means position, the second means read.
    std::set<std::pair<long, read_t>> reads;
};

class ContigGenerator {
public:
    std::set<Contig *> contigs;

    ContigGenerator(ReadAligner *rA, ReadData *rD, ReadFilter *rF);

    ~ContigGenerator();

    void generateContigs();

    friend std::ostream &operator<<(std::ostream &out,
                                    const ContigGenerator &o);

private:
    ReadData *rD;

    /**
     * This is the read aligner used to align the reads as they are being added
     * to the Contig
     */
    ReadAligner *rA;

    ReadFilter *rF;

    bool hasActiveContig = false;
    Contig *activeContig = nullptr;

    // We maintain a list of all reads that have not been treated and a map from
    // reads to contigs
    std::unordered_map<read_t, std::pair<Contig *, long>> readsInContig;
    std::unordered_set<read_t> reads2Contig;

    bool addRelatedReads(const std::pair<long, read_t> r);

    /**
     * Merges the second contig into the first one
     * @param c1
     * @param c2
     * @param pos
     */
    void mergeContigs(Contig *c1, Contig *c2, long pos);

    void initialize();
};

#endif // EXPERIMENTS_CONTIG_CUH
