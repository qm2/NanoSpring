#ifndef Z_MYERS_H
#define Z_MYERS_H

#include <vector>
#include <string>
#include "Edits.h"


class MyersAligner : public StringAligner {
public:

    /***
    * Uses the myers algorithm to calculate the optimal edit script from the first
    * string to the second string
    * @param s1 The first string (the original string)
    * @param s2 The second string (the target string)
    * @param editScript stores a vector of edits that will transform the first
    * string into the second.
    * @param editDis the edit distance obtained by this algorithm
    * @return whether alignment succeeded
    */
    bool align(const std::string &s1, const std::string &s2,
               std::vector<Edit> &editScript, size_t &editDis);

    MyersAligner();

    /***
     * Runs the myers algorithm to compare the two strings and returns an edit script
     * @param A
     * First string
     * @param B
     * Second string
     * @param editScript
     * @param editDis the edit distance obtained by this algorithm
     * a vector of edits that is supposed to be read from back to front.
     */
    static void myers(const std::string &s1, const std::string &s2,
                      std::vector<Edit> &editScript, size_t &editDis);

};

/***
 * Given strings A and B, let a0 and b0 be initial substrings with length len.
 * Run myers on a0 and b0 until at least one of them is fully aligned.
 * Advance a0 and b0 to new substrings of length len and repeat.
 */
class PiecewiseMyers : public StringAligner {
public:
    /***
     * The length of substrings to run myers algorithm on
     */
    const size_t lenA;
    const size_t lenB;

    PiecewiseMyers(const size_t len);

    PiecewiseMyers(const size_t lenA, const size_t lenB);

    bool align(const std::string &s1, const std::string &s2,
               std::vector<Edit> &editScript, size_t &editDis);

private:

    /***
     * Performs alignment via myers algorithm on two strings and returns when
     * one of them is fully aligned
     * @param Abegin The start of the first string. Will be updated to reflect alignment
     * @param Aend The end of the first string.
     * @param Bbegin The start of the second string. Will be updated to reflect alignment
     * @param Bend The end of the second string.
     * @param max Maximum number of edits to search
     * @param editScript The vector to store the edits
     * @param editDis The edit distance
     * @return Whether alignment succeeded within max steps
     */
    static bool localAlign(const char *&Abegin, const char *const Aend,
                           const char *&Bbegin, const char *const Bend,
                           const size_t max,
                           std::vector<Edit> &editScript, size_t &editDis);
};

#endif //Z_MYERS_H
