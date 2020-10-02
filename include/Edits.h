#ifndef EXPERIMENTS_EDITS_H
#define EXPERIMENTS_EDITS_H

#include <ostream>
#include <string>
#include <vector>

typedef enum { SAME, INSERT, DELETE, SUBSTITUTION } EDIT_TYPE;

/**
 * Class denoting a single edit. editType is either SAME (no change), or INSERT,
 * or DELETE. If editType is SAME, we use the num field in editInfo to denote
 * the number of consecutive characters left unchanged. If editType is INSERT,
 * we use the ins field to denote the inserted character If editType is DELETE,
 * we use the del field to donote the deleted character
 */
class Edit {
public:
    EDIT_TYPE editType;
    union {
        size_t num;
        char del;
        char ins;
        char sub;
    } editInfo;

    Edit(EDIT_TYPE editType, size_t num);

    friend std::ostream &operator<<(std::ostream &out, const Edit &o);

    /**
     * Changes adjacent insertions and deletions to substitutions. The
     * oldEditScript can only have insertions and deletions. It can't have
     * substitutions
     * @return The new edit distance
     */
    static size_t optimizeEditScript(std::vector<Edit> &oldEditScript,
                                     std::vector<Edit> &newEditScript);
};

/**
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

#endif // EXPERIMENTS_EDITS_H