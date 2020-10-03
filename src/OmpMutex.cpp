#include "OmpMutex.h"

OmpMutex::OmpMutex() { omp_init_lock(&internalLock); }

OmpMutex::~OmpMutex() { omp_destroy_lock(&internalLock); }

void OmpMutex::lock() { omp_set_lock(&internalLock); }

void OmpMutex::test() { omp_test_lock(&internalLock); }

void OmpMutex::unlock() { omp_unset_lock(&internalLock); }

OmpNestMutex::OmpNestMutex() { omp_init_nest_lock(&internalLock); }

OmpNestMutex::~OmpNestMutex() { omp_destroy_nest_lock(&internalLock); }

void OmpNestMutex::lock() { omp_set_nest_lock(&internalLock); }

void OmpNestMutex::test() { omp_test_nest_lock(&internalLock); }

void OmpNestMutex::unlock() { omp_unset_nest_lock(&internalLock); }