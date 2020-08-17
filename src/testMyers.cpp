#include "Edits.h"
#include "LocalMyers.h"
#include "LocalMyersRollBack.h"
#include "LocalMyersRollBackOld.h"
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
        std::vector<Edit> editScript;
        size_t editDis;
        LocalMyersRollBack pM(100, 100, 1000);
        pM.align(s1, s2, editScript, editDis);
        std::stringstream result;
        for (auto e = editScript.begin(); e < editScript.end(); e++) {
            result << *e;
        }
        std::cout << editDis << std::endl;
        std::cout << result.str() << std::endl;
    }

    return 0;
}
