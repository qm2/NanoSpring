#include "Edits.h"
#include <unordered_map>
#include <map>
#include <iostream>
#include <string>
#include <stack>
#include <sstream>
#include <vector>
#include <myers.h>


int main() {

    std::string s1, s2;

    // Input strings
    std::cin >> s1 >> s2;

    std::vector<Edit> *editScript = myers(s1, s2);
    std::stringstream result;
    for (auto e = editScript->end() - 1; e >= editScript->begin(); e--) {
        result << *e;
    }
    std::cout << result.str() << std::endl;
    delete editScript;
    return 0;

}

