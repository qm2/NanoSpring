#include "LocalMyersRollBack.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>

template <typename RandomAccessIt>
bool LocalMyersRollBack::localAlign(RandomAccessIt &Abegin, RandomAccessIt Aend,
                                    RandomAccessIt &Bbegin, RandomAccessIt Bend,
                                    const size_t max,
                                    std::vector<Edit> &editScript,
                                    size_t &editDis) {
    unsigned const int lenAString = Aend - Abegin;
    unsigned const int lenBString = Bend - Bbegin;
    bool foundEdit = false;

    editScript.clear();

    // Each member maps end point to a portion of the path leading to that end
    // point
    std::map<std::pair<int, int>, EditPath> editInfo;

    int Xreached = lenAString;
    int Yreached = lenBString;
    auto findEdit = [lenAString, lenBString, Abegin, Bbegin, max, &Xreached,
                     &Yreached, &editInfo, &editDis, &foundEdit]() {
        const size_t VSize = max + 1;
        int VTemp[VSize * 4];
        int *V[4];
        {
            int *cur = VTemp + (int)(max / 2);
            for (int i = 0; i < 4; cur += VSize, ++i)
                V[i] = cur;
        }
        auto findSnakeAndUpdate = [lenAString, lenBString, V, Abegin, Bbegin,
                                   &Xreached, &Yreached, &editInfo, &foundEdit](
                                      const int xStart, const int yStart,
                                      const int xMid, const int yMid,
                                      const int d, const int k) {
            // end point
            int xEnd = xMid;
            int yEnd = yMid;
            // we go along the diagonal until we fail
            size_t snakeLen = 0;
            while (xEnd < (int)lenAString && yEnd < (int)lenBString &&
                   Abegin[xEnd] == Bbegin[yEnd]) {
                ++xEnd;
                ++yEnd;
                ++snakeLen;
            }

            // we save this furthest reaching end point
            V[d % 4][k] = xEnd;

            editInfo.insert(
                std::make_pair(std::make_pair(xEnd, yEnd),
                               EditPath(xStart, yStart, xMid, yMid, snakeLen)));

            // std::cout << d << ',' << k << ':' << xStart << ' ' << yStart << '
            // '
            //<< xMid << ' ' << yMid << ' ' << xEnd << ' ' << yEnd
            //<< '\n';

            // whether we have found an edit
            if (xEnd >= (int)lenAString || yEnd >= (int)lenBString) {
                // std::cout
                //     << "Minimum edits required to convert string A into B
                //     is : "
                //     << d << std::endl;
                Xreached = xEnd;
                Yreached = yEnd;
                // std::cout << Xreached << " " << Yreached << '\n';
                foundEdit = true;
            }
        };

        const int MIN_VALUE = -2;

        // First we deal with the case where d = 0
        int d = 0;
        {
            int k = 0;

            int xStart = 0;
            int yStart = 0;

            // middle point
            int xMid = 0;
            int yMid = 0;

            findSnakeAndUpdate(xStart, yStart, xMid, yMid, d, k);
        }

        // Now the case where d = 1;
        d++;
        if (foundEdit || d > (int)max)
            return;

        V[d][0] = MIN_VALUE;

        // Now the case where d = 2;
        d++;
        if (foundEdit || d > (int)max)
            return;
        {
            // k ranges from -1 to 1
            // First the case where k = -1
            {
                int k = -1;
                // we need to go down
                int VDown = V[(d + 2) % 4][k + 1];
                int kPrev = k + 1;

                int xStart = VDown;
                int yStart = xStart - kPrev;

                int xMid = xStart;
                int yMid = xMid - k;

                findSnakeAndUpdate(xStart, yStart, xMid, yMid, d, k);
                if (foundEdit)
                    return;
            }

            // Now the case where k = 0
            {
                int k = 0;
                V[d][k] = MIN_VALUE;
            }

            // Finally the case where k = 1
            {
                int k = 1;
                // we need to go right
                int VRight = V[(d + 2) % 4][k - 1];
                int kPrev = k - 1;

                int xStart = VRight;
                int yStart = xStart - kPrev;

                int xMid = xStart + 1;
                int yMid = xMid - k;

                findSnakeAndUpdate(xStart, yStart, xMid, yMid, d, k);
                if (foundEdit)
                    return;
            }
        }

        // Now for d>=3
        ++d;
        for (; d <= (int)max; d++) {
            int halfD = d / 2;
            for (int k = -halfD; k <= halfD; ++k) {
                int VDown = 2 * (k + 1) >= -(d - 2) && 2 * (k + 1) <= d - 2
                                ? V[(d + 2) % 4][k + 1]
                                : MIN_VALUE;
                int VRight = 2 * (k - 1) >= -(d - 2) && 2 * (k - 1) <= d - 2
                                 ? V[(d + 2) % 4][k - 1]
                                 : MIN_VALUE;
                int VS = 2 * k >= -(d - 3) && 2 * k <= d - 3 ? V[(d + 1) % 4][k]
                                                             : MIN_VALUE;
                // 1 means VDown; -1 means VRight; 0 means VS
                int bestOldPath;
                int xStart;
                if (VDown > VRight) {
                    if (VDown > VS) {
                        bestOldPath = 1;
                        xStart = VDown;
                    } else {
                        bestOldPath = 0;
                        xStart = VS;
                    }
                } else {
                    if (VRight > VS) {
                        bestOldPath = -1;
                        xStart = VRight;
                    } else {
                        bestOldPath = 0;
                        xStart = VS;
                    }
                }
                if (xStart < 0) {
                    V[d % 4][k] = MIN_VALUE;
                    continue;
                }
                int kPrev = k + bestOldPath;

                // starting point
                int yStart = xStart - kPrev;

                // middle point
                int xMid = bestOldPath == 1 ? xStart : xStart + 1;
                int yMid = xMid - k;

                findSnakeAndUpdate(xStart, yStart, xMid, yMid, d, k);
                if (foundEdit)
                    return;
            }
        }
    };

    findEdit();

    if (!foundEdit)
        return false;

    editDis = 0;
    int currentX = Xreached;
    int currentY = Yreached;
    editScript.reserve(max);
    while (currentX > 0 || currentY > 0) {
        // std::cout << currentX << " " << currentY << std::endl;
        EditPath &e = editInfo.at(std::make_pair(currentX, currentY));
        if (e.snakeLen > 0) {
            // If there is a snake (a series of diagonals)
            editScript.push_back(Edit(SAME, e.snakeLen));
        }
        if (e.xMid > e.xStart) {
            editDis++;
            if (e.yMid > e.yStart) {
                // There is a substitution
                // editScript.push_back(
                // Edit(SUBSTITUTION,
                //      forward ? Bbegin[e.yStart] : Bbegin[-e.yStart]));
                editScript.push_back(Edit(DELETE, Abegin[e.xStart]));
                editScript.push_back(Edit(INSERT, Bbegin[e.yStart]));
            } else {
                // If we moved right, then this is a deletion
                editScript.push_back(Edit(DELETE, Abegin[e.xStart]));
            }
        } else if (e.yMid > e.yStart) {
            editDis++;
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

    // {
    //     std::cout << editDis << std::endl;
    //     for (auto e = editScript.begin(); e < editScript.end(); e++) {
    //         std::cout << *e;
    //     }
    //     std::cout << std::endl;
    // }
    return true;
}

LocalMyersRollBack::LocalMyersRollBack(const size_t lenA, const size_t lenB,
                                       const size_t maxEditDis)
    : LocalMyers("LMRB " + std::to_string(lenA) + " " + std::to_string(lenB) +
                     " MaxEditDis " + std::to_string(maxEditDis),
                 lenA, lenB),
      maxEditDis(maxEditDis) {}

template <typename RandomAccessIt>
bool LocalMyersRollBack::alignOnce(RandomAccessIt Abegin, RandomAccessIt Aend,
                                   RandomAccessIt Bbegin, RandomAccessIt Bend,
                                   const ssize_t offsetGuess,
                                   ssize_t &beginOffset, ssize_t &endOffset,
                                   std::vector<Edit> &editScript,
                                   size_t &editDis) {
    beginOffset = offsetGuess;
    endOffset = 0;
    editDis = 0;
    editScript.clear();

    if (offsetGuess > 0) {
        if (Aend - Abegin <= offsetGuess)
            return false;
        Abegin += offsetGuess;
    } else {
        if (Bend - Bbegin <= -offsetGuess)
            return false;
        Bbegin += (-offsetGuess);
    }
    const size_t max = std::min(lenA, lenB) * 2;

    while (Abegin != Aend && Bbegin != Bend) {
        RandomAccessIt ALocalEnd =
            Aend - Abegin > (ssize_t)lenA ? Abegin + lenA : Aend;
        RandomAccessIt BLocalEnd =
            Bend - Bbegin > (ssize_t)lenB ? Bbegin + lenB : Bend;
        std::vector<Edit> localEditScript;
        size_t localEditDis;
        bool success = localAlign(Abegin, ALocalEnd, Bbegin, BLocalEnd, max,
                                  localEditScript, localEditDis);
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
        endOffset = Abegin - Aend;
    } else if (Bbegin != Bend) {
        // We still need to insert the rest
        endOffset = Bend - Bbegin;
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
    const size_t max = std::min(lenA, lenB) * 2;
    auto const advance = [max,
                          this](const char *&Abegin, const char *const Aend,
                                const char *&Bbegin, const char *const Bend,
                                size_t lenAString, size_t lenBString,
                                bool &dirSuccess, size_t &editDis) {
        const char *const ALocalEnd = std::min(Abegin + lenAString, Aend);
        const char *const BLocalEnd = std::min(Bbegin + lenBString, Bend);
        std::vector<Edit> localEditScript;
        size_t localEditDis;
        bool success = localAlign(Abegin, ALocalEnd, Bbegin, BLocalEnd, max,
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
    if (!dir1Success && !dir2Success) {
#ifdef DEBUG
        std::cerr << "initial alignment failed\n";
#endif
        return false;
    }

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

    // std::cout << "forwardPassEndOffset " << forwardPassEndOffset << '\n';

    ssize_t backwardPassBeginOffset, backwardPassEndOffset;

    bool secondPassSuccess = alignOnce(
        std::reverse_iterator<const char *>(Aend),
        std::reverse_iterator<const char *>(s1.c_str()),
        std::reverse_iterator<const char *>(Bend),
        std::reverse_iterator<const char *>(s2.c_str()), -forwardPassEndOffset,
        backwardPassBeginOffset, backwardPassEndOffset, editScript, editDis);
    if (secondPassSuccess) {
        beginOffset = -backwardPassEndOffset;
        endOffset = -backwardPassBeginOffset;
        std::reverse(editScript.begin(), editScript.end());
        // std::vector<Edit> newEditScript;
        // editDis = Edit::optimizeEditScript(editScript, newEditScript);
        // editScript.swap(newEditScript);
        // {
        //     std::cout << "BeginOffset " << beginOffset << " EndOffset "
        //               << endOffset << '\n';
        //     std::cout << editDis << std::endl;
        //     for (auto e = editScript.begin(); e < editScript.end(); e++) {
        //         std::cout << *e;
        //     }
        //     std::cout << std::endl;
        // }
        return true;
    }
    // std::cerr << "Backward alignment failed\n";
    return false;
}

bool LocalMyersRollBack::align(const std::string &s1, const std::string &s2,
                               std::vector<Edit> &editScript, size_t &editDis) {
    const ssize_t offsetGuess = 0;
    ssize_t beginOffset;
    ssize_t endOffset;
    return align(s1, s2, offsetGuess, beginOffset, endOffset, editScript,
                 editDis);
}
