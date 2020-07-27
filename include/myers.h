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
     * @return whether alignment succeeded
     */
    bool align(const std::string &s1, const std::string &s2,
               std::vector<Edit> &editScript);

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
    const size_t len;

    PiecewiseMyers(const size_t len);
};

#endif //Z_MYERS_H
