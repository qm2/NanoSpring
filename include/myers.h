#ifndef Z_MYERS_H
#define Z_MYERS_H

#include <vector>
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

#endif //Z_MYERS_H
