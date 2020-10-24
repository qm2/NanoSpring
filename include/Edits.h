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
    EditPath() {}; // default constructor needed for initialization of matrix
};

namespace Edits {
/**
 * @brief
 *
 * @tparam Iterator Random access read iterator that dereferences to char
 * @tparam Inserter Write inserter that dereferences to char
 * @param original
 * @param editScript
 * @param result
 */
template <typename Iterator, typename Inserter>
void applyEdits(Iterator original, const std::vector<Edit> &editScript,
                Inserter result) {
    for (const Edit &e : editScript) {
        switch (e.editType) {
        case SAME:
            for (size_t i = 0; i < e.editInfo.num; ++i)
                *(result++) = *(original++);
            break;
        case INSERT:
            *(result++) = e.editInfo.ins;
            break;
        case DELETE:
            ++original;
            break;
        case SUBSTITUTION:
            *(result++) = e.editInfo.sub;
            ++original;
            break;
        }
    }
}
} // namespace Edits

#endif // EXPERIMENTS_EDITS_H
