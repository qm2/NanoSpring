#ifndef EXPERIMENTS_EDITS_H
#define EXPERIMENTS_EDITS_H

#include <ostream>

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

#endif //EXPERIMENTS_EDITS_H