#ifndef DF93CB52_1F2B_4A32_9D3B_FEA65A4B8F96
#define DF93CB52_1F2B_4A32_9D3B_FEA65A4B8F96
#include "LocalMyers.h"
#include <algorithm>
#include <iostream>
#include <map>

template <typename RandomAccessItA, typename RandomAccessItB>
LocalMyers<RandomAccessItA, RandomAccessItB>::LocalMyers(const size_t len)
    : StringAligner<RandomAccessItA, RandomAccessItB>("LocalMyers Len " +
                                                      std::to_string(len)),
      lenA(len), lenB(len) {}

template <typename RandomAccessItA, typename RandomAccessItB>
LocalMyers<RandomAccessItA, RandomAccessItB>::LocalMyers(const size_t lenA,
                                                         const size_t lenB)
    : StringAligner<RandomAccessItA, RandomAccessItB>(
          "LocalMyers Lens " + std::to_string(lenA) + " and " +
          std::to_string(lenB)),
      lenA(lenA), lenB(lenB) {}

template <typename RandomAccessItA, typename RandomAccessItB>
bool LocalMyers<RandomAccessItA, RandomAccessItB>::align(
    RandomAccessItA s1Begin, RandomAccessItA s1End, RandomAccessItB s2Begin,
    RandomAccessItB s2End, const ssize_t offsetGuess, ssize_t &beginOffset,
    ssize_t &endOffset, std::vector<Edit> &editScript, size_t &editDis) {
    (void)offsetGuess;
    beginOffset = 0;
    endOffset = 0;
    RandomAccessItA Abegin = s1Begin;
    RandomAccessItA Aend = s1End;
    RandomAccessItB Bbegin = s2Begin;
    RandomAccessItB Bend = s2End;
    editDis = 0;
    editScript.clear();
    while (Abegin < Aend && Bbegin < Bend) {
        const size_t max = std::min(lenA, lenB) * 2;
        RandomAccessItA ALocalEnd =
            Aend - Abegin > (ssize_t)lenA ? Abegin + lenA : Aend;
        RandomAccessItB BLocalEnd =
            Bend - Bbegin > (ssize_t)lenB ? Bbegin + lenB : Bend;
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

template <typename RandomAccessItA, typename RandomAccessItB>
bool LocalMyers<RandomAccessItA, RandomAccessItB>::localAlign(
    RandomAccessItA &Abegin, RandomAccessItA Aend, RandomAccessItB &Bbegin,
    RandomAccessItB Bend, const size_t max, std::vector<Edit> &editScript,
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
    for (int d = 0; d <= (int)max; d++) {
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
            size_t snakeLen = 0;
            while (xEnd < (int)lenA && yEnd < (int)lenB &&
                   Abegin[xEnd] == Bbegin[yEnd]) {
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
            if (xEnd >= (int)lenA || yEnd >= (int)lenB) {
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

template <typename RandomAccessItA, typename RandomAccessItB>
LocalMyers<RandomAccessItA, RandomAccessItB>::LocalMyers(
    const std::string &name, const size_t lenA, const size_t lenB)
    : StringAligner<RandomAccessItA, RandomAccessItB>(name), lenA(lenA),
      lenB(lenB) {}
#endif /* DF93CB52_1F2B_4A32_9D3B_FEA65A4B8F96 */
