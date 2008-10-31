#include <pthread.h>

#ifndef BARRIER_T
#define BARRIER_T

// Threads things. Threads is a way of distributing workload across
// different workers to solve a problem faster. Nevertheless, sometimes
// we need synchronization between threads to avoid undesired race
// conditions and other problems. Here we are an implementation of a barrier
// that allows putting all threads to wait at a given point untill all of them
// have reached such point and can continue working. Barriers are usually
// available through pthreads system library. Nonetheless, sometimes it is not
// so we have to implement it here.

typedef struct mybarrier_t {
    int needed;
    int called;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
}barrier_t; 

int barrier_init(barrier_t *barrier,int needed);
int barrier_destroy(barrier_t *barrier);
int barrier_wait(barrier_t *barrier);

#endif
