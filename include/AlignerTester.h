#ifndef Z_ALIGNERTESTER_H
#define Z_ALIGNERTESTER_H

#include "Edits.h"
#include "StringAligner.h"
#include <cstddef>
#include <string>
#include <vector>

class AlignerTester {
public:
    /**
     * Create random reads.
     * @param readLen read length
     * @param offset offset between the start of the two reads
     * @param num number of pairs of reads to generate
     * @param pIns insertion error rate
     * @param pDel deletion error rate
     * @param pS substitution error rate
     */
    void generateData(size_t readLen, ssize_t offset, size_t num,
                      double pIns = 0.03, double pDel = 0.03, double pS = 0.04);

    /**
     * Runs the algorithm contained in aligner on the data generatad
     * @param aligner
     * @param duration Will store the average time taken in seconds
     * @param successRate Will store probability that the aligner succeeds
     * @param avgEditDis the average edit distance returned by the aligner
     */
    void profile(StringAligner<const char *> *aligner, double &duration,
                 double &successRate, double &avgEditDis);

    /**
     * Runs the algorithm contained in aligner on the data generatad
     * @param aligner
     * @param duration Will store the average time taken in seconds
     * @param successRate Will store probability that the aligner succeeds
     * @param avgBeginOffset the average beginOffset
     * @param avgEndOffset the average endOffset
     * @param avgEditDis the average edit distance returned by the aligner
     */
    void profile(StringAligner<const char *> *aligner, double &duration,
                 double &successRate, double &avgBeginOffset,
                 double &avgEndOffset, double &avgEditDis);

    /**
     * Validates on the generated data that aligner->align returns true
     * and that the editScript is correct
     * @param aligner
     * @return
     */
    bool validate(StringAligner<const char *> *aligner);

private:
    std::vector<std::string> readsA, readsB;
};

#endif // Z_ALIGNERTESTER_H