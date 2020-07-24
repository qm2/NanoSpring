#ifndef EXPERIMENTS_EDITS_H
#define EXPERIMENTS_EDITS_H

#include <ostream>
#include <vector>
#include <string>

typedef enum {
    SAME, INSERT, DELETE
} EDIT_TYPE;

/***
 * Class denoting a single edit. editType is either SAME (no change), or INSERT, or DELETE.
 * If editType is SAME, we use the num field in editInfo to denote the number of
 * consecutive characters left unchanged.
 * If editType is INSERT, we use the ins field to denote the inserted character
 * If editType is DELETE, we use the del field to donote the deleted character
 */
class Edit {
public:
    EDIT_TYPE editType;
    union {
        unsigned int num;
        char del;
        char ins;
    } editInfo;

    Edit(EDIT_TYPE editType, char c);

    Edit(EDIT_TYPE editType, unsigned int num);

    friend std::ostream &operator<<(std::ostream &out, const Edit &o);
};


class StringAligner {
public:
    /***
     * Calculates a good edit script from string s1 to string s2 and returns it in
     * a vector of edits
     * @param s1 The first string (the original string)
     * @param s2 The second string (the target string)
     * @return a vector of edits that will transform the first string into the second.
     * The user is responsible for freeing the vector
     */
    virtual std::vector<Edit> *align(const std::string &s1, const std::string &s2) = 0;
};

#endif //EXPERIMENTS_EDITS_H