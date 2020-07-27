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
     * Name of the alignment algorithm
     */
    const std::string name;

    /***
     * Calculates a good edit script from string s1 to string s2 and stores it in
     * a vector of edits. Returns whether the alignment succeeded
     * @param s1 The first string (the original string)
     * @param s2 The second string (the target string)
     * @param editScript stores a vector of edits that will transform the first
     * string into the second.
     * @return whether alignment succeeded
     */
    bool align(const std::string &s1, const std::string &s2,
               std::vector<Edit> &editScript);

    /***
     * Calculates a good edit script from string s1 to string s2 and stores it in
     * a vector of edits. Returns whether the alignment succeeded
     * @param s1 The first string (the original string)
     * @param s2 The second string (the target string)
     * @param editScript stores a vector of edits that will transform the first
     * string into the second.
     * @param editDis the edit distance obtained by this algorithm
     * @return whether alignment succeeded
     */
    virtual bool align(const std::string &s1, const std::string &s2,
                       std::vector<Edit> &editScript, size_t &editDis) = 0;

    StringAligner(const std::string &name);
};

#endif //EXPERIMENTS_EDITS_H