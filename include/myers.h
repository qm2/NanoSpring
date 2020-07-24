#ifndef Z_MYERS_H
#define Z_MYERS_H

#include <vector>
#include <string>
#include "Edits.h"

/***
 * Runs the myers algorithm to compare the two strings and returns an edit script
 * @param A
 * First string
 * @param B
 * Second string
 * @return
 * Returns a vector of edits that is supposed to be read from back to front. The vector
 * is dynamically allocated and the user is responsible for freeing it.
 */
std::vector<Edit> *myers(const std::string &s1, const std::string &s2);

class MyersAligner : public StringAligner {
public:
    /***
     * Uses the myers algorithm to calculate the optimal edit script from the first
     * string to the second string
     * @param s1 The first string (the original string)
     * @param s2 The second string (the target string)
     * @return a vector of edits that will transform the first string into the second.
     * The user is responsible for freeing the vector.
     */
    std::vector<Edit> *align(const std::string &s1, const std::string &s2);
};

#endif //Z_MYERS_H
