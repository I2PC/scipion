#include <pthread.h>

#ifndef BARRIER_T
#define BARRIER_T

/**@defgroup Threads Threads
   @ingroup DataLibrary

   Threads is a way of distributing workload across
   different workers to solve a problem faster. Nevertheless, sometimes
   we need synchronization between threads to avoid undesired race
   conditions and other problems. Here we are an implementation of a barrier
   that allows putting all threads to wait at a given point until all of them
   have reached such point and can continue working. Barriers are usually
   available through pthreads system library. Nonetheless, sometimes it is not
   so we have to implement it here.
*/
//@{

/** Barrier structure */
typedef struct mybarrier_t {
	/// How many threads should be awaited
    int needed;
    /// How many threads already arrived
    int called;
    /// Mutex to update this structure
    pthread_mutex_t mutex;
    /// Condition on which the threads are waiting
    pthread_cond_t cond;
} barrier_t;

/** Barrier initialization */
int barrier_init(barrier_t *barrier, int needed);
/** Barrier destruction */
int barrier_destroy(barrier_t *barrier);
/** Wait at the barrier */
int barrier_wait(barrier_t *barrier);
//@}
#endif
