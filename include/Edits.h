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
    virtual bool align(const std::string &s1, const std::string &s2,
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

    /***
     * Calculates a good edit script from string s1 to string s2 and stores it in
     * a vector of edits. Here the edit script only transforms the overlapping
     * portion of s1 to the overlapping portion of s2. The tail and head that are
     * not aligned are stored as beginOffset and endOffset.
     * Returns whether the alignment succeeded.
     * @param s1 The first string (the original string)
     * @param s2 The second string (the target string)
     * @param offsetGuess An initial guess of beginOffset
     * @param beginOffset The alignment of the left most ends of the two reads.
     * Positive means the second read goes after the first, and negative means the
     * second read goes before the first.
     * @param endOffset The alignment of the right most ends of the two reads.
     * Positive means the second read goes after the first, and negative means the
     * second read goes before the first.
     * @param editScript stores a vector of edits that will transform the
     * overlapping portion of the first string into the second.
     * @param editDis the edit distance obtained by this algorithm
     * (this edit distance is only of the overlapping portions)s
     * @return whether alignment succeeded
     */
    virtual bool align(const std::string &s1, const std::string &s2,
                       const ssize_t offsetGuess,
                       ssize_t &beginOffset, ssize_t &endOffset,
                       std::vector<Edit> &editScript, size_t &editDis);

    StringAligner(const std::string &name);

    virtual ~StringAligner();
};

#endif //EXPERIMENTS_EDITS_H