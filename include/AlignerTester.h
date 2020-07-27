#ifndef Z_ALIGNERTESTER_H
#define Z_ALIGNERTESTER_H

#include <cstddef>
#include <vector>
#include <string>
#include "Edits.h"

class AlignerTester {
public:
    /***
     * Create random reads.
     * @param readLen read length
     * @param offset offset between the start of the two reads
     * @param num number of pairs of reads to generate
     * @param pIns insertion error rate
     * @param pDel deletion error rate
     * @param pS substitution error rate
     */
    void generateData(size_t readLen, size_t offset, size_t num, double pIns = 0.03,
                      double pDel = 0.03,
                      double pS = 0.04
    );

    /***
     * Runs the algorithm contained in aligner on the data generatad
     * @param aligner
     * @param duration Will store the average time taken in seconds
     * @param successRate Will store probability that the aligner succeeds
     * @param avgEditDis the average edit distance returned by the aligner
     */
    void profile(StringAligner *aligner, double &duration, double &successRate,
                 double &avgEditDis);

    /***
     * Validates on the generated data that aligner->align returns true
     * and that the editScript is correct
     * @param aligner
     * @return
     */
    bool validate(StringAligner *aligner);

private:
    std::vector<std::string> readsA, readsB;

    /***
     * Applies editScript to origString and stores in result
     * @param origString The original string to transform
     * @param editScript The edit script to apply
     * @param result stores the resulting string
     * @return
     */
    static void applyEditsToString(const std::string &origString,
                                   const std::vector<Edit> &editScript, std::string &result);
};

#endif //Z_ALIGNERTESTER_H