#include "Edits.h"

Edit::Edit(EDIT_TYPE editType, unsigned int num) : editType(editType) {
    this->editInfo.num = num;
}

Edit::Edit(EDIT_TYPE editType, char c) : editType(editType) {
    if (editType == INSERT) {
        this->editInfo.ins = c;
    } else {
        this->editInfo.del = c;
    }
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
    }
    return out;
}


StringAligner::StringAligner(const std::string &name) : name(name) {}

bool StringAligner::align(const std::string &s1, const std::string &s2, std::vector<Edit> &editScript) {
    size_t editDis;
    return align(s1, s2, editScript, editDis);
}

bool StringAligner::align(const std::string &s1, const std::string &s2,
                          const ssize_t offsetGuess,
                          ssize_t &beginOffset, ssize_t &endOffset,
                          std::vector<Edit> &editScript, size_t &editDis) {
    beginOffset = 0;
    endOffset = 0;
    return align(s1, s2, editScript, editDis);
}