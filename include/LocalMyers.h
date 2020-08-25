#ifndef E8B69C0C_68E6_42E2_B825_13682A52A2DD
#define E8B69C0C_68E6_42E2_B825_13682A52A2DD

#include "Edits.h"
#include "MyersAligner.h"

/***
 * Given strings A and B, let a0 and b0 be initial substrings with length len.
 * Run myers on a0 and b0 until at least one of them is fully aligned.
 * Advance a0 and b0 to new substrings of length len and repeat.
 */
class LocalMyers : public StringAligner {
public:
    /***
     * The length of substrings to run myers algorithm on
     */
    const size_t lenA;
    const size_t lenB;

    LocalMyers(const size_t len);

    LocalMyers(const size_t lenA, const size_t lenB);

    virtual bool align(const std::string &s1, const std::string &s2,
                       std::vector<Edit> &editScript, size_t &editDis) override;

protected:
    /***
     * Performs alignment via myers algorithm on two strings and returns when
     * one of them is fully aligned
     * @param Abegin The start of the first string. Will be updated to reflect
     * alignment
     * @param Aend The end of the first string.
     * @param Bbegin The start of the second string. Will be updated to reflect
     * alignment
     * @param Bend The end of the second string.
     * @param max Maximum number of edits to search
     * @param editScript The vector to store the edits. The edits will be stored
     * backwards
     * @param editDis The edit distance
     * @return Whether alignment succeeded within max steps
     */
    static bool localAlign(const char *&Abegin, const char *const Aend,
                           const char *&Bbegin, const char *const Bend,
                           const size_t max, std::vector<Edit> &editScript,
                           size_t &editDis);

    LocalMyers(const std::string &name, const size_t lenA, const size_t lenB);
};

#endif /* E8B69C0C_68E6_42E2_B825_13682A52A2DD */