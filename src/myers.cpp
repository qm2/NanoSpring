#include "myers.h"
#include <map>
#include <iostream>
#include <algorithm>

/***
 * An edit path in the myers algorithm consists of a starting point (from another furthest reaching path),
 * and a horizontal or vertical move to the mid point, and finally a snake (all diagonal moves) to the
 * end point.
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

EditPath::EditPath(int xStart, int yStart, int xMid, int yMid, unsigned int snakeLen)
        : xStart(xStart), yStart(yStart), xMid(xMid), yMid(yMid), snakeLen(snakeLen) {}

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

    // Each member maps end point to a portion of the path leading to that end point
    std::map<std::pair<int, int>, EditPath> editInfo;

    for (int d = 0; d <= max; d++) {
        for (int k = -d; k <= d; k += 2) {
            // we either choose to go down or go right
            // We choose to go down if 1. we have no choice or 2. going down is the better choice
            bool down = (k == -d || (k < d && V[k - 1] < V[k + 1]));

            int kPrev = down ? k + 1 : k - 1;

            // starting point
            int xStart = V[kPrev];
            int yStart = xStart - kPrev;
            // In the first iteration the starting point calculated would be 0, -1.
            // We want to make it 0, 0
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
            editInfo.insert(std::make_pair(std::make_pair(xEnd, yEnd),
                                           EditPath(xStart, yStart, xMid, yMid, snakeLen)));

            // whether we have found an edit
            if (xEnd >= len1 && yEnd >= len2) {
//                std::cout << "Minimum edits required to convert string A into B is: " << d << std::endl;
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

PiecewiseMyers::PiecewiseMyers(const size_t len) :
        StringAligner("Piecewise Myers Len " + std::to_string(len)), lenA(len), lenB(len) {}

PiecewiseMyers::PiecewiseMyers(const size_t lenA, const size_t lenB) :
        StringAligner("Piecewise Myers Lens "
                      + std::to_string(lenA)
                      + " and "
                      + std::to_string(lenB)),
        lenA(lenA), lenB(lenB) {}

bool PiecewiseMyers::align(const std::string &s1, const std::string &s2,
                           std::vector<Edit> &editScript, size_t &editDis) {
    const char *Abegin = s1.c_str();
    const char *const Aend = Abegin + s1.length();
    const char *Bbegin = s2.c_str();
    const char *const Bend = Bbegin + s2.length();
    editDis = 0;
    editScript.clear();
    while (Abegin != Aend && Bbegin != Bend) {
        const size_t max = lenA < lenB ? lenA : lenB;
        const char *const ALocalEnd = Abegin + lenA < Aend ? Abegin + lenA : Aend;
        const char *const BLocalEnd = Bbegin + lenB < Bend ? Bbegin + lenB : Bend;
        std::vector<Edit> localEditScript;
        size_t localEditDis;
        bool success = localAlign(Abegin, ALocalEnd, Bbegin, BLocalEnd,
                                  max, localEditScript, localEditDis);
        if (!success)
            return false;
//        std::cout << Aend - Abegin << " "
//                  << Bend - Bbegin << std::endl;
        editScript.insert(editScript.end(), localEditScript.begin(), localEditScript.end());
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

bool PiecewiseMyers::localAlign(const char *&Abegin, const char *const Aend,
                                const char *&Bbegin, const char *const Bend,
                                const size_t max,
                                std::vector<Edit> &editScript, size_t &editDis) {
    unsigned const int lenA = Aend - Abegin;
    unsigned const int lenB = Bend - Bbegin;
    bool foundEdit = false;

    int VTemp[2 * max + 1];

    int *V = VTemp + max;

    V[1] = 0;

    editScript.clear();

    // Each member maps end point to a portion of the path leading to that end point
    std::map<std::pair<int, int>, EditPath> editInfo;

    int Xreached = lenA;
    int Yreached = lenB;
    for (int d = 0; d <= max; d++) {
        for (int k = -d; k <= d; k += 2) {
            // we either choose to go down or go right
            // We choose to go down if 1. we have no choice or 2. going down is the better choice
            bool down = (k == -d || (k < d && V[k - 1] < V[k + 1]));

            int kPrev = down ? k + 1 : k - 1;

            // starting point
            int xStart = V[kPrev];
            int yStart = xStart - kPrev;
            // In the first iteration the starting point calculated would be 0, -1.
            // We want to make it 0, 0
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
            editInfo.insert(std::make_pair(std::make_pair(xEnd, yEnd),
                                           EditPath(xStart, yStart, xMid, yMid, snakeLen)));

            // whether we have found an edit
            if (xEnd >= lenA || yEnd >= lenB) {
//                std::cout << "Minimum edits required to convert string A into B is: " << d << std::endl;
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
    std::reverse(editScript.begin(), editScript.end());
    return true;
}