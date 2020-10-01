#ifndef CA5F280E_564A_4065_BECC_348C1184C9DC
#define CA5F280E_564A_4065_BECC_348C1184C9DC

#include "Edits.h"

/**
 * @brief
 *
 * @tparam RandomAccessIt Random Access Read Iterator that dereferences to
 * char.
 */
template <typename RandomAccessItA, typename RandomAccessItB = RandomAccessItA>
class StringAligner {
public:
    /***
     * Name of the alignment algorithm
     */
    const std::string name;

    /***
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
    __attribute__((warn_unused_result)) virtual bool
    align(RandomAccessItA Abegin, RandomAccessItA Aend, RandomAccessItB Bbegin,
          RandomAccessItB Bend, const ssize_t offsetGuess, ssize_t &beginOffset,
          ssize_t &endOffset, std::vector<Edit> &editScript,
          size_t &editDis) = 0;

    StringAligner(const std::string &name) : name(name){};

    virtual ~StringAligner(){};
};

#endif /* CA5F280E_564A_4065_BECC_348C1184C9DC */