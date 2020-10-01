#include "Edits.h"
#include "LocalMyers.h"
#include "LocalMyersRollBack.h"
#include "LocalMyersRollBack_impl.h"
#include <iostream>
#include <map>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

int main() {

    std::string s1, s2;

    // Input strings
    std::cin >> s1 >> s2;

    // {
    //     std::vector<Edit> editScript;
    //     size_t editDis;
    //     MyersAligner::myers(s1, s2, editScript, editDis);
    //     std::stringstream result;
    //     for (auto e = editScript.end() - 1; e >= editScript.begin(); e--) {
    //         result << *e;
    //     }
    //     std::cout << result.str() << std::endl;
    // }

    // {
    //     std::vector<Edit> editScript;
    //     size_t editDis;
    //     LocalMyers pM(2);
    //     pM.align(s1, s2, editScript, editDis);
    //     std::stringstream result;
    //     for (auto e = editScript.begin(); e < editScript.end(); e++) {
    //         result << *e;
    //     }
    //     std::cout << editDis << std::endl;
    //     std::cout << result.str() << std::endl;
    // }

    {
        const char *const s1Start = s1.c_str();
        const char *const s1End = s1.c_str() + s1.length();
        const char *const s2Start = s2.c_str();
        const char *const s2End = s2.c_str() + s2.length();
        std::vector<Edit> editScript;
        size_t editDis;
        LocalMyersRollBack<const char *> pM(100, 100, 1000);
        ssize_t beginOffset, endOffset;
        pM.align(s1Start, s1End, s2Start, s2End, 0, beginOffset, endOffset,
                 editScript, editDis);
        std::stringstream result;
        for (auto e = editScript.begin(); e < editScript.end(); e++) {
            result << *e;
        }
        std::cout << editDis << " beginOffset " << beginOffset << " endOffset "
                  << endOffset << std::endl;
        std::cout << result.str() << std::endl;
    }

    return 0;
}
