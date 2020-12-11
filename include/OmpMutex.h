#ifndef AC67E3E1_065C_493C_A86B_661928304612
#define AC67E3E1_065C_493C_A86B_661928304612

#include <omp.h>

/**
 * @brief Mutex-wrapper for omp_lock_t. Provides the .lock() and .unlock()
 * methods
 *
 */
class OmpMutex {
public:
    void lock();
    void unlock();
    bool try_lock();
    OmpMutex();
    ~OmpMutex();

private:
    omp_lock_t internalLock;
};

/**
 * @brief Mutex-wrapper for omp_nest_lock_t. Provides the .lock() and .unlock()
 * methods
 *
 */
class OmpNestMutex {
public:
    void lock();
    void unlock();
    bool try_lock();
    OmpNestMutex();
    ~OmpNestMutex();

private:
    omp_nest_lock_t internalLock;
};
#endif /* AC67E3E1_065C_493C_A86B_661928304612 */
