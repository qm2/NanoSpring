#include "myers.h"
#include <algorithm>
#include <iostream>
#include <map>

/***
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
    unsigned int snakeLen;

    EditPath(int xStart, int yStart, int xMid, int yMid, unsigned int snakeLen);
};

EditPath::EditPath(int xStart, int yStart, int xMid, int yMid,
                   unsigned int snakeLen)
    : xStart(xStart), yStart(yStart), xMid(xMid), yMid(yMid),
      snakeLen(snakeLen) {}

void MyersAligner::myers(const std::string &s1, const std::string &s2,
                         std::vector<Edit> &editScript, size_t &editDis) {
    unsigned const int len1 = s1.length();
    unsigned const int len2 = s2.length();
    unsigned const int max = len1 + len2;
    bool foundEdit = false;

    int VTemp[2 * max + 1];

    int *V = VTemp + max;

    V[1] = 0;

    editScript.clear();

    // Each member maps end point to a portion of the path leading to that end
    // point
    std::map<std::pair<int, int>, EditPath> editInfo;

    for (int d = 0; d <= max; d++) {
        for (int k = -d; k <= d; k += 2) {
            // we either choose to go down or go right
            // We choose to go down if 1. we have no choice or 2. going down is
            // the better choice
            bool down = (k == -d || (k < d && V[k - 1] < V[k + 1]));

            int kPrev = down ? k + 1 : k - 1;

            // starting point
            int xStart = V[kPrev];
            int yStart = xStart - kPrev;
            // In the first iteration the starting point calculated would be 0,
            // -1. We want to make it 0, 0
            if (yStart < 0)
                yStart = 0;

            // middle point
            int xMid = down ? xStart : xStart + 1;
            int yMid = xMid - k;

            // end point
            int xEnd = xMid;
            int yEnd = yMid;

            // we go along the diagonal until we fail
            unsigned int snakeLen = 0;
            while (xEnd < len1 && yEnd < len2 && s1[xEnd] == s2[yEnd]) {
                ++xEnd;
                ++yEnd;
                ++snakeLen;
            }

            // we save this furthest reaching end point
            V[k] = xEnd;

            //            std::cout << xStart << " " << yStart << std::endl;
            editInfo.insert(
                std::make_pair(std::make_pair(xEnd, yEnd),
                               EditPath(xStart, yStart, xMid, yMid, snakeLen)));

            // whether we have found an edit
            if (xEnd >= len1 && yEnd >= len2) {
                //                std::cout << "Minimum edits required to
                //                convert string A into B is: " << d <<
                //                std::endl;
                editDis = d;
                foundEdit = true;
                break;
            }
        }

        if (foundEdit)
            break;
    }

    int currentX = len1;
    int currentY = len2;
    editScript.reserve(max);
    while (currentX > 0 || currentY > 0) {
        //        std::cout << currentX << " " << currentY << std::endl;
        EditPath &e = editInfo.at(std::make_pair(currentX, currentY));
        if (e.snakeLen > 0) {
            // If there is a snake (a series of diagonals)
            editScript.push_back(Edit(SAME, e.snakeLen));
        }
        if (e.xMid > e.xStart) {
            // If we moved right, then this is a deletion
            editScript.push_back(Edit(DELETE, s1[e.xStart]));
        } else if (e.yMid > e.yStart) {
            // If we moved down and there is an insertion
            editScript.push_back(Edit(INSERT, s2[e.yStart]));
        }
        currentX = e.xStart;
        currentY = e.yStart;
    }

    editScript.shrink_to_fit();
}

bool MyersAligner::align(const std::string &s1, const std::string &s2,
                         std::vector<Edit> &editScript, size_t &editDis) {
    myers(s1, s2, editScript, editDis);
    std::reverse(editScript.begin(), editScript.end());
    return true;
}

MyersAligner::MyersAligner() : StringAligner("Myers Algorithm") {}

LocalMyers::LocalMyers(const size_t len)
    : StringAligner("LocalMyers Len " + std::to_string(len)), lenA(len),
      lenB(len) {}

LocalMyers::LocalMyers(const size_t lenA, const size_t lenB)
    : StringAligner("LocalMyers Lens " + std::to_string(lenA) + " and " +
                    std::to_string(lenB)),
      lenA(lenA), lenB(lenB) {}

bool LocalMyers::align(const std::string &s1, const std::string &s2,
                       std::vector<Edit> &editScript, size_t &editDis) {
    const char *Abegin = s1.c_str();
    const char *const Aend = Abegin + s1.length();
    const char *Bbegin = s2.c_str();
    const char *const Bend = Bbegin + s2.length();
    editDis = 0;
    editScript.clear();
    while (Abegin < Aend && Bbegin < Bend) {
        const size_t max = lenA < lenB ? lenA : lenB;
        const char *const ALocalEnd =
            Abegin + lenA < Aend ? Abegin + lenA : Aend;
        const char *const BLocalEnd =
            Bbegin + lenB < Bend ? Bbegin + lenB : Bend;
        std::vector<Edit> localEditScript;
        size_t localEditDis;
        bool success = localAlign(Abegin, ALocalEnd, Bbegin, BLocalEnd, max,
                                  localEditScript, localEditDis);
        if (!success)
            return false;
        //        std::cout << Aend - Abegin << " "
        //                  << Bend - Bbegin << std::endl;
        std::reverse(localEditScript.begin(), localEditScript.end());
        editScript.insert(editScript.end(), localEditScript.begin(),
                          localEditScript.end());
        editDis += localEditDis;
    }
    if (Abegin != Aend) {
        // We still need to delete the rest
        while (Abegin < Aend) {
            editDis++;
            editScript.push_back(Edit(DELETE, *Abegin++));
        }
    } else if (Bbegin != Bend) {
        // We still need to insert the rest
        while (Bbegin < Bend) {
            editDis++;
            editScript.push_back(Edit(INSERT, *Bbegin++));
        }
    }
    return true;
}

bool LocalMyers::localAlign(const char *&Abegin, const char *const Aend,
                            const char *&Bbegin, const char *const Bend,
                            const size_t max, std::vector<Edit> &editScript,
                            size_t &editDis) {
    unsigned const int lenA = Aend - Abegin;
    unsigned const int lenB = Bend - Bbegin;
    bool foundEdit = false;

    int VTemp[2 * max + 1];

    int *V = VTemp + max;

    V[1] = 0;

    editScript.clear();

    // Each member maps end point to a portion of the path leading to that end
    // point
    std::map<std::pair<int, int>, EditPath> editInfo;

    int Xreached = lenA;
    int Yreached = lenB;
    for (int d = 0; d <= max; d++) {
        for (int k = -d; k <= d; k += 2) {
            // we either choose to go down or go right
            // We choose to go down if 1. we have no choice or 2. going down is
            // the better choice
            bool down = (k == -d || (k < d && V[k - 1] < V[k + 1]));

            int kPrev = down ? k + 1 : k - 1;

            // starting point
            int xStart = V[kPrev];
            int yStart = xStart - kPrev;
            // In the first iteration the starting point calculated would be 0,
            // -1. We want to make it 0, 0
            if (yStart < 0)
                yStart = 0;

            // middle point
            int xMid = down ? xStart : xStart + 1;
            int yMid = xMid - k;

            // end point
            int xEnd = xMid;
            int yEnd = yMid;

            // we go along the diagonal until we fail
            unsigned int snakeLen = 0;
            while (xEnd < lenA && yEnd < lenB && Abegin[xEnd] == Bbegin[yEnd]) {
                ++xEnd;
                ++yEnd;
                ++snakeLen;
            }

            // we save this furthest reaching end point
            V[k] = xEnd;

            //            std::cout << xStart << " " << yStart << std::endl;
            editInfo.insert(
                std::make_pair(std::make_pair(xEnd, yEnd),
                               EditPath(xStart, yStart, xMid, yMid, snakeLen)));

            // whether we have found an edit
            if (xEnd >= lenA || yEnd >= lenB) {
                //                std::cout << "Minimum edits required to
                //                convert string A into B is: " << d <<
                //                std::endl;
                Xreached = xEnd;
                Yreached = yEnd;
                editDis = d;
                foundEdit = true;
                break;
            }
        }

        if (foundEdit)
            break;
    }

    if (!foundEdit)
        return false;

    int currentX = Xreached;
    int currentY = Yreached;
    editScript.reserve(max);
    while (currentX > 0 || currentY > 0) {
        //        std::cout << currentX << " " << currentY << std::endl;
        EditPath &e = editInfo.at(std::make_pair(currentX, currentY));
        if (e.snakeLen > 0) {
            // If there is a snake (a series of diagonals)
            editScript.push_back(Edit(SAME, e.snakeLen));
        }
        if (e.xMid > e.xStart) {
            // If we moved right, then this is a deletion
            editScript.push_back(Edit(DELETE, Abegin[e.xStart]));
        } else if (e.yMid > e.yStart) {
            // If we moved down and there is an insertion
            editScript.push_back(Edit(INSERT, Bbegin[e.yStart]));
        }
        currentX = e.xStart;
        currentY = e.yStart;
    }

    editScript.shrink_to_fit();
    // Update Abegin and Bbegin
    Abegin += Xreached;
    Bbegin += Yreached;
    return true;
}

LocalMyers::LocalMyers(const std::string &name, const size_t lenA,
                       const size_t lenB)
    : StringAligner(name), lenA(lenA), lenB(lenB) {}

bool LocalMyersRollBack::localAlign(const char *&Abegin, const char *const Aend,
                                    const char *&Bbegin, const char *const Bend,
                                    const size_t max,
                                    std::vector<Edit> &editScript,
                                    size_t &editDis, bool forward) {
    unsigned const int lenA = forward ? Aend - Abegin : Abegin - Aend;
    unsigned const int lenB = forward ? Bend - Bbegin : Bbegin - Bend;
    bool foundEdit = false;

    int VTemp[2 * max + 1];

    int *V = VTemp + max;

    V[1] = 0;

    editScript.clear();

    // Each member maps end point to a portion of the path leading to that end
    // point
    std::map<std::pair<int, int>, EditPath> editInfo;

    int Xreached = lenA;
    int Yreached = lenB;
    for (int d = 0; d <= max; d++) {
        for (int k = -d; k <= d; k += 2) {
            // we either choose to go down or go right
            // We choose to go down if 1. we have no choice or 2. going down is
            // the better choice
            bool down = (k == -d || (k < d && V[k - 1] < V[k + 1]));

            int kPrev = down ? k + 1 : k - 1;

            // starting point
            int xStart = V[kPrev];
            int yStart = xStart - kPrev;
            // In the first iteration the starting point calculated would be 0,
            // -1. We want to make it 0, 0
            if (yStart < 0)
                yStart = 0;

            // middle point
            int xMid = down ? xStart : xStart + 1;
            int yMid = xMid - k;

            // end point
            int xEnd = xMid;
            int yEnd = yMid;

            // we go along the diagonal until we fail
            unsigned int snakeLen = 0;
            while (xEnd < lenA && yEnd < lenB &&
                   ((forward && Abegin[xEnd] == Bbegin[yEnd]) ||
                    (!forward && Abegin[-xEnd] == Bbegin[-yEnd]))) {
                ++xEnd;
                ++yEnd;
                ++snakeLen;
            }

            // we save this furthest reaching end point
            V[k] = xEnd;

            //            std::cout << xStart << " " << yStart << std::endl;
            editInfo.insert(
                std::make_pair(std::make_pair(xEnd, yEnd),
                               EditPath(xStart, yStart, xMid, yMid, snakeLen)));

            // whether we have found an edit
            if (xEnd >= lenA || yEnd >= lenB) {
                //                std::cout << "Minimum edits required to
                //                convert string A into B is: " << d <<
                //                std::endl;
                Xreached = xEnd;
                Yreached = yEnd;
                editDis = d;
                foundEdit = true;
                break;
            }
        }

        if (foundEdit)
            break;
    }

    if (!foundEdit)
        return false;

    int currentX = Xreached;
    int currentY = Yreached;
    editScript.reserve(max);
    while (currentX > 0 || currentY > 0) {
        //        std::cout << currentX << " " << currentY << std::endl;
        EditPath &e = editInfo.at(std::make_pair(currentX, currentY));
        if (e.snakeLen > 0) {
            // If there is a snake (a series of diagonals)
            editScript.push_back(Edit(SAME, e.snakeLen));
        }
        if (e.xMid > e.xStart) {
            // If we moved right, then this is a deletion
            editScript.push_back(
                Edit(DELETE, forward ? Abegin[e.xStart] : Abegin[-e.xStart]));
        } else if (e.yMid > e.yStart) {
            // If we moved down and there is an insertion
            editScript.push_back(
                Edit(INSERT, forward ? Bbegin[e.yStart] : Bbegin[-e.yStart]));
        }
        currentX = e.xStart;
        currentY = e.yStart;
    }

    editScript.shrink_to_fit();
    // Update Abegin and Bbegin
    Abegin += forward ? Xreached : -Xreached;
    Bbegin += forward ? Yreached : -Yreached;
    return true;
}

LocalMyersRollBack::LocalMyersRollBack(const size_t lenA, const size_t lenB,
                                       const size_t maxEditDis)
    : LocalMyers("LMRB " + std::to_string(lenA) + " " + std::to_string(lenB) +
                     " MaxEditDis " + std::to_string(maxEditDis),
                 lenA, lenB),
      maxEditDis(maxEditDis) {}

bool LocalMyersRollBack::alignReverse(const std::string &s1,
                                      const std::string &s2,
                                      const ssize_t offsetGuess,
                                      ssize_t &beginOffset, ssize_t &endOffset,
                                      std::vector<Edit> &editScript,
                                      size_t &editDis) {
    const char *const Aend = s1.c_str() - 1;
    const char *Abegin = Aend + s1.length();
    const char *const Bend = s2.c_str() - 1;
    const char *Bbegin = Bend + s2.length();
    beginOffset = offsetGuess;
    endOffset = 0;
    editDis = 0;
    editScript.clear();
    if (offsetGuess > 0) {
        Abegin -= offsetGuess;
        if (Abegin <= Aend)
            return false;
    } else {
        Bbegin -= (-offsetGuess);
        if (Bbegin <= Bend)
            return false;
    }
    const size_t max = std::min(lenA, lenB);

    while (Abegin != Aend && Bbegin != Bend) {
        const char *const ALocalEnd = std::max(Abegin - lenA, Aend);
        const char *const BLocalEnd = std::max(Bbegin - lenB, Bend);
        std::vector<Edit> localEditScript;
        size_t localEditDis;
        bool success = localAlign(Abegin, ALocalEnd, Bbegin, BLocalEnd, max,
                                  localEditScript, localEditDis, false);
        if (!success)
            return false;
        editDis += localEditDis;
        if (editDis + errorRate * std::min(Abegin - Aend, Bbegin - Bend) >
            maxEditDis)
            return false;
        std::reverse(localEditScript.begin(), localEditScript.end());
        editScript.insert(editScript.end(), localEditScript.begin(),
                          localEditScript.end());
    }
    if (Abegin != Aend) {
        // We still need to delete the rest
        endOffset = Aend - Abegin;
    } else if (Bbegin != Bend) {
        // We still need to insert the rest
        endOffset = Bbegin - Bend;
    }
    return true;
}

bool LocalMyersRollBack::align(const std::string &s1, const std::string &s2,
                               const ssize_t offsetGuess, ssize_t &beginOffset,
                               ssize_t &endOffset,
                               std::vector<Edit> &editScript, size_t &editDis) {
    const char *Abegin1 = s1.c_str();
    const char *const Aend = Abegin1 + s1.length();
    const char *Bbegin1 = s2.c_str();
    const char *const Bend = Bbegin1 + s2.length();
    ssize_t beginOffset1 = offsetGuess;
    ssize_t beginOffset2 = -offsetGuess;
    ssize_t endOffset1 = 0;
    ssize_t endOffset2 = 0;
    size_t editDis1 = 0;
    size_t editDis2 = 0;

    if (offsetGuess > 0) {
        Abegin1 += offsetGuess;
        if (Abegin1 >= Aend)
            return false;
    } else {
        Bbegin1 += (-offsetGuess);
        if (Bbegin1 > Bend)
            return false;
    }
    const char *Abegin2 = Abegin1;
    const char *Bbegin2 = Bbegin1;
    bool dir1Success = true;
    bool dir2Success = true;
    const size_t max = std::min(lenA, lenB);
    auto const advance = [max,
                          this](const char *&Abegin, const char *const Aend,
                                const char *&Bbegin, const char *const Bend,
                                size_t lenA, size_t lenB, bool &dirSuccess,
                                size_t &editDis) {
        const char *const ALocalEnd = std::min(Abegin + lenA, Aend);
        const char *const BLocalEnd = std::min(Bbegin + lenB, Bend);
        std::vector<Edit> localEditScript;
        size_t localEditDis;
        bool success =
            LocalMyers::localAlign(Abegin, ALocalEnd, Bbegin, BLocalEnd, max,
                                   localEditScript, localEditDis);
        if (success) {
            editDis += localEditDis;
            if (editDis + errorRate * std::min(Aend - Abegin, Bend - Bbegin) >
                maxEditDis)
                dirSuccess = false;
        } else {
            dirSuccess = false;
        }
    };
    while (Abegin1 < Aend && Bbegin1 < Bend && Abegin2 < Aend &&
           Bbegin2 < Bend && (dir1Success || dir2Success)) {
        double expectedEditDis1 =
            editDis1 + errorRate * std::min(Aend - Abegin1, Bend - Bbegin1);
        double expectedEditDis2 =
            editDis2 + errorRate * std::min(Aend - Abegin2, Bend - Bbegin2);

        if ((!dir2Success) ||
            (dir1Success && expectedEditDis1 < expectedEditDis2)) {

            advance(Abegin1, Aend, Bbegin1, Bend, lenA, lenB, dir1Success,
                    editDis1);
        } else {
            advance(Bbegin2, Bend, Abegin2, Aend, lenA, lenB, dir2Success,
                    editDis2);
        }
    }
    if (!dir1Success && !dir2Success)
        return false;
    ssize_t forwardPassEndOffset = 0;
    if (Abegin1 == Aend || Bbegin1 == Bend) {
        // std::cout << "Left\n";
        if (Abegin1 != Aend) {
            // We still need to delete the rest
            forwardPassEndOffset = Abegin1 - Aend;
        } else if (Bbegin1 != Bend) {
            // We still need to insert the rest
            forwardPassEndOffset = Bend - Bbegin1;
        }
    } else {
        // std::cout << "Right\n";
        if (Abegin2 != Aend) {
            // We still need to delete the rest
            forwardPassEndOffset = Abegin2 - Aend;
        } else if (Bbegin2 != Bend) {
            // We still need to insert the rest
            forwardPassEndOffset = Bend - Bbegin2;
        }
    }
    ssize_t backwardPassBeginOffset, backwardPassEndOffset;
    bool secondPassSuccess =
        alignReverse(s1, s2, -forwardPassEndOffset, backwardPassBeginOffset,
                     backwardPassEndOffset, editScript, editDis);
    if (secondPassSuccess) {
        beginOffset = -backwardPassEndOffset;
        endOffset = -backwardPassBeginOffset;
        std::reverse(editScript.begin(), editScript.end());
        return true;
    }
    return false;
    // firstPassSuccess = align(s2, s1, -offsetGuess,
    // forwardPassBeginOffset,
    //                          forwardPassEndOffset, false, editScript,
    //                          editDis);
    // if (!firstPassSuccess)
    //     return false;
    // bool secondPassSuccess = align(
    //     s1Reverse, s2Reverse, forwardPassEndOffset,
    //     backwardPassBeginOffset, backwardPassEndOffset, true, editScript,
    //     editDis);
    // if (!secondPassSuccess)
    //     return false;
    // beginOffset = -backwardPassEndOffset;
    // endOffset = -backwardPassBeginOffset;
    // std::reverse(editScript.begin(), editScript.end());

    // if (Abegin != Aend) {
    //     // We still need to delete the rest
    //     endOffset = Abegin - Aend;
    // } else if (Bbegin != Bend) {
    //     // We still need to insert the rest
    //     endOffset = Bend - Bbegin;
    // }
    // std::string s1Reverse = s1;
    // std::string s2Reverse = s2;
    // std::reverse(s1Reverse.begin(), s1Reverse.end());
    // std::reverse(s2Reverse.begin(), s2Reverse.end());
    // ssize_t forwardPassBeginOffset, forwardPassEndOffset;
    // ssize_t backwardPassBeginOffset, backwardPassEndOffset;

    // bool firstPassSuccess =
    //     align(s1, s2, offsetGuess, forwardPassBeginOffset,
    //     forwardPassEndOffset,
    //           false, editScript, editDis);

    // if (firstPassSuccess) {
    //     bool secondPassSuccess =
    //         align(s1Reverse, s2Reverse, -forwardPassEndOffset,
    //               backwardPassBeginOffset, backwardPassEndOffset, true,
    //               editScript, editDis);
    //     if (secondPassSuccess) {
    //         beginOffset = -backwardPassEndOffset;
    //         endOffset = -backwardPassBeginOffset;
    //         std::reverse(editScript.begin(), editScript.end());
    //         return true;
    //     }
    // }

    return true;
}

bool LocalMyersRollBack::align(const std::string &s1, const std::string &s2,
                               std::vector<Edit> &editScript, size_t &editDis) {
    const ssize_t offsetGuess = 0;
    ssize_t beginOffset;
    ssize_t endOffset;
    return align(s1, s2, offsetGuess, beginOffset, endOffset, editScript,
                 editDis);
}
