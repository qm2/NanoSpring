#ifndef EEF78146_670F_4E04_ABCE_4FD978240D5F
#define EEF78146_670F_4E04_ABCE_4FD978240D5F

#include "Edits.h"

Edit::Edit(EDIT_TYPE editType, size_t num) : editType(editType) {
    switch (editType) {
    case INSERT:
        this->editInfo.ins = num;
        break;
    case DELETE:
        this->editInfo.del = num;
        break;
    case SUBSTITUTION:
        this->editInfo.sub = num;
        break;
    case SAME:
        this->editInfo.num = num;
        break;
    }
}

size_t Edit::optimizeEditScript(std::vector<Edit> &oldEditScript,
                                std::vector<Edit> &newEditScript) {
    size_t editDis = 0;
    newEditScript.clear();
    std::vector<Edit>::iterator editIt = oldEditScript.begin();
    std::vector<Edit>::iterator end = oldEditScript.end();
    while (editIt < end) {
        while (editIt < end && editIt->editType == SAME) {
            newEditScript.push_back(*editIt);
            editIt++;
        }
        std::string insertChars;
        size_t delNum = 0;
        while (editIt < end && editIt->editType != SAME) {
            if (editIt->editType == INSERT)
                insertChars.push_back(editIt->editInfo.ins);
            else
                delNum++;
            editIt++;
        }
        size_t insNum = insertChars.size();
        editDis += std::max(delNum, insNum);
        size_t numSubs = std::min(insNum, delNum);
        size_t i;
        for (i = 0; i < numSubs; ++i) {
            newEditScript.push_back(Edit(SUBSTITUTION, insertChars[i]));
        }
        if (insNum > delNum) {
            for (; i < insNum; ++i)
                newEditScript.push_back(Edit(INSERT, insertChars[i]));
        } else {
            // std::cout << "more d\n";
            for (; i < delNum; ++i)
                newEditScript.push_back(Edit(DELETE, '-'));
        }
    }
    return editDis;
}

std::ostream &operator<<(std::ostream &out, const Edit &o) {
    switch (o.editType) {
    case SAME:
        out << 'u' << o.editInfo.num;
        break;
    case INSERT:
        out << 'i' << o.editInfo.ins;
        break;
    case DELETE:
        out << 'd' << o.editInfo.del;
        break;
    case SUBSTITUTION:
        out << 's' << o.editInfo.sub;
        break;
    }
    return out;
}

EditPath::EditPath(int xStart, int yStart, int xMid, int yMid, size_t snakeLen)
    : xStart(xStart), yStart(yStart), xMid(xMid), yMid(yMid),
      snakeLen(snakeLen) {}

#endif /* EEF78146_670F_4E04_ABCE_4FD978240D5F */