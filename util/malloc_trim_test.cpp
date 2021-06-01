/*
Short code for testing if system supports malloc_trim(0);
If not, then we disable it to allow successful compilation.
*/

#include <malloc.h>

int main() {
    malloc_trim(0);
    return 0;
}
