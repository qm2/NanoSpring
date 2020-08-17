#ifndef Z_MYERS_H
#define Z_MYERS_H

#include "Edits.h"
#include <string>
#include <vector>

/***
 * An edit path in the myers algorithm consists of a starting point (from
 * another furthest reaching path), and a horizontal or vertical move to the mid
 * point, and finally a snake (all diagonal moves) to the end point.
 */
class EditPath {
public:
    // start point
    int xStart;
    int yStart;

    // mid point
    int xMid;
    int yMid;

    // snake length
    size_t snakeLen;

    EditPath(int xStart, int yStart, int xMid, int yMid, size_t snakeLen);
};

class MyersAligner : public StringAligner {
public:
    /***
     * Uses the myers algorithm to calculate the optimal edit script from the
     * first string to the second string
     * @param s1 The first string (the original string)
     * @param s2 The second string (the target string)
     * @param editScript stores a vector of edits that will transform the first
     * string into the second.
     * @param editDis the edit distance obtained by this algorithm
     * @return whether alignment succeeded
     */
    virtual bool align(const std::string &s1, const std::string &s2,
                       std::vector<Edit> &editScript, size_t &editDis) override;

    MyersAligner();

    /***
     * Runs the myers algorithm to compare the two strings and returns an edit
     * script
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

#endif // Z_MYERS_H