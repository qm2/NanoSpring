#ifndef F1ADC39F_6A7F_4671_B4B5_77EF545E2B5C
#define F1ADC39F_6A7F_4671_B4B5_77EF545E2B5C

#include "Edits.h"
#include "LocalMyers.h"
#include "MyersAligner.h"

/***
 * This class basically runs the local myers algorithm from left the right to
 * obtain a good alignment at the right most end, and then uses that alignment
 * to run the local myers algorithm from right to left. In internal
 * calculations, the cost of substituion is hardcoded as 3, while the costs of
 * insertion and deletion are 2. Or equivalently, the cost of insertion and
 * deletion is 1, and the cost of substitution is 1.5. However, in the actual
 * returned editDis, the Levenshtein distance is used, i.e., insertions,
 * deletions, and substitutions all have edit distance 1.
 */

class LocalMyersRollBack : public LocalMyers {
public:
    /***
     * The maximum edit distance to search for among the overlapping portions
     */
    const size_t maxEditDis;

    const double errorRate = 0.2;

    virtual bool align(const std::string &s1, const std::string &s2,
                       std::vector<Edit> &editScript, size_t &editDis) override;

    virtual bool align(const std::string &s1, const std::string &s2,
                       const ssize_t offsetGuess, ssize_t &beginOffset,
                       ssize_t &endOffset, std::vector<Edit> &editScript,
                       size_t &editDis) override;

    LocalMyersRollBack(const size_t lenA, const size_t lenB,
                       const size_t maxEditDis);

    /**
     * @brief
     * Performs alignment using an modification of myers algorithm where
     * substitution cost is 1.5. Alignment stops when one of the strings reaches
     * its end.
     * @param Abegin
     * @param Aend
     * @param Bbegin
     * @param Bend
     * @param max
     * @param editScript
     * @param editDis This will store the Levenshtein distance is used, i.e.,
     * insertions, deletions, and substitutions all have edit distance 1.
     * @param forward True if we are aligning from left to right, and false if
     * we are aligning from right to left.
     * @return true
     * @return false
     */
    static bool localAlign(const char *&Abegin, const char *const Aend,
                           const char *&Bbegin, const char *const Bend,
                           const size_t max, std::vector<Edit> &editScript,
                           size_t &editDis, bool forward);

private:
    bool alignReverse(const std::string &s1, const std::string &s2,
                      const ssize_t offsetGuess, ssize_t &beginOffset,
                      ssize_t &endOffset, std::vector<Edit> &editScript,
                      size_t &editDis);
};
#endif /* F1ADC39F_6A7F_4671_B4B5_77EF545E2B5C */
