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
