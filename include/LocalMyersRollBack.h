#ifndef F1ADC39F_6A7F_4671_B4B5_77EF545E2B5C
#define F1ADC39F_6A7F_4671_B4B5_77EF545E2B5C

#include "Edits.h"
#include "LocalMyers.h"
#include "StringAligner.h"

/**
 * This class basically runs the local myers algorithm from left the right to
 * obtain a good alignment at the right most end, and then uses that alignment
 * to run the local myers algorithm from right to left. In internal
 * calculations, the cost of substituion is hardcoded as 3, while the costs of
 * insertion and deletion are 2. Or equivalently, the cost of insertion and
 * deletion is 1, and the cost of substitution is 1.5. However, in the actual
 * returned editDis, the Levenshtein distance is used, i.e., insertions,
 * deletions, and substitutions all have edit distance 1.
 */

/**
 * @brief
 *
 * @tparam RandomAccessIt Random Access Read Iterator that dereferences to
 * char.
 */
template <typename RandomAccessItA, typename RandomAccessItB = RandomAccessItA>
class LocalMyersRollBack : public LocalMyers<RandomAccessItA, RandomAccessItB> {
public:
    /**
     * The maximum edit distance to search for among the overlapping portions
     */
    const size_t maxEditDis;

    const double errorRate = 0.2;

    /**
     * Calculates a good edit script from string s1 to string s2 and stores it
     * in a vector of edits. Here the edit script only transforms the
     * overlapping portion of s1 to the overlapping portion of s2. The tail and
     * head that are not aligned are stored as beginOffset and endOffset.
     * Returns whether the alignment succeeded.
     * @param Abegin Start of original string
     * @param Aend End of original string
     * @param Bbegin Start of target string
     * @param Bend End of target string
     * @param offsetGuess An initial guess of beginOffset
     * @param beginOffset The alignment of the left most ends of the two reads.
     * Positive means the second read goes after the first, and negative means
     * the second read goes before the first.
     * @param endOffset The alignment of the right most ends of the two reads.
     * Positive means the second read goes after the first, and negative means
     * the second read goes before the first.
     * @param editScript stores a vector of edits that will transform the
     * overlapping portion of the first string into the second.
     * @param editDis the edit distance obtained by this algorithm
     * (this edit distance is only of the overlapping portions)s
     * @return whether alignment succeeded
     */
    virtual bool align(RandomAccessItA Abegin, RandomAccessItA Aend,
                       RandomAccessItB Bbegin, RandomAccessItB Bend,
                       const ssize_t offsetGuess, ssize_t &beginOffset,
                       ssize_t &endOffset, std::vector<Edit> &editScript,
                       size_t &editDis) override;

    LocalMyersRollBack(const size_t lenA, const size_t lenB,
                       const size_t maxEditDis);

protected:
    /**
     * @brief
     * Performs alignment using an modification of myers algorithm where
     * substitution cost is 1.5. Alignment stops when one of the strings reaches
     * its end.
     * @tparam RIt Random Access Read Iterator that dereferences to
     * char. This may be the reverse iterator of RandomAccessIt
     * @param Abegin
     * @param Aend
     * @param Bbegin
     * @param Bend
     * @param max
     * @param editScript
     * @param editDis This will store the Levenshtein distance used, i.e.,
     * insertions, deletions, and substitutions all have edit distance 1.
     * @return true
     * @return false
     */
    template <typename RItA, typename RItB>
    static bool localAlign(RItA &Abegin, RItA Aend, RItB &Bbegin, RItB Bend,
                           const size_t max, std::vector<Edit> &editScript,
                           size_t &editDis);

    /**
     * @brief Performs global alignment algorithm once
     *
     * @tparam RIt Random Access Read Iterator that dereferences to
     * char. This may be the reverse iterator of RandomAccessIt
     * @param Abegin
     * @param Aend
     * @param Bbegin
     * @param Bend
     * @param offsetGuess
     * @param beginOffset
     * @param endOffset
     * @param editScript
     * @param editDis
     * @return true
     * @return false
     */
    template <typename RItA, typename RItB>
    bool alignOnce(RItA Abegin, RItA Aend, RItB Bbegin, RItB Bend,
                   const ssize_t offsetGuess, ssize_t &beginOffset,
                   ssize_t &endOffset, std::vector<Edit> &editScript,
                   size_t &editDis);

    /**
     * @brief Performs local alignment, advances Abegin and Bbegin, and report
     * whether we should continue
     *
     * @tparam RItA
     * @tparam RItB
     * @param Abegin
     * @param Aend
     * @param Bbegin
     * @param Bend
     * @param lenAString
     * @param lenBString
     * @param dirSuccess Whether we should continue based on maxEditDis and the
     * current editDis
     * @param editDis
     * @param max Maximum edit distance in a single alignment
     */
    template <typename RItA, typename RItB>
    void advance(RItA &Abegin, RItA Aend, RItB &Bbegin, RItB Bend,
                 const size_t lenAString, const size_t lenBString,
                 bool &dirSuccess, size_t &editDis, const size_t max);
};
#endif /* F1ADC39F_6A7F_4671_B4B5_77EF545E2B5C */
