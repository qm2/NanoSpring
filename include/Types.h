#ifndef TYPES_H_
#define TYPES_H_

#include <cstdint>
// We use the indices of the reads to represent them
typedef uint32_t read_t;

/** We store all the k-mers as uint64s. This would work for all k<=32,
 which is definitely sufficient **/
typedef uint64_t kMer_t;

#endif // TYPES_H_
