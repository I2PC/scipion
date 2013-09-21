/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef THREADS_T
#define THREADS_T

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

class ThreadManager;
class ThreadArgument;
class Condition;

/* Prototype of functions for threads works. */
typedef void (*ThreadFunction) (ThreadArgument &arg);

//TODO (MARIANA) Please give more documentation and in a good structure e.g. @name (see args.h as example)

/** @defgroup Threads Threads
 *  @ingroup ParallelLibrary
 //@{
 */

/** Class wrapping around the pthreads mutex.
 * This class will provide a more object oriented implementation
 * of a mutex, to ensure the unique access to critical regions of
 * code and other synchronization problems.
 */
class Mutex
{
private:
    pthread_mutex_t mutex; //our pthread mutex

public:
    /** Default constructor.
     * This constructor just initialize the pthread_mutex_t structure
     * with its defaults values, just like static initialization with PTHREAD_MUTEX_INITIALIZER
     */
    Mutex();

    /** Destructor. */
    virtual ~Mutex();

    /** Function to get the access to the mutex.
     * If the some thread has the mutex and other
     * ask to lock will be waiting until the first one
     * release the mutex
     */
    virtual void lock();

    /** Function to release the mutex.
     * This allow the access to the mutex to other
     * threads that are waiting for it.
     */
    virtual void unlock();

    friend class Condition;
}
;//end of class Mutex

/** Class wrapping around the pthreads condition.
 * This class will provide a more object oriented implementation
 * of a condition variable to achieve syncronization between threads.
 */
class Condition
{
private:
    Mutex *mutex; //our pthread mutex
    pthread_cond_t cond; //our pthread cond

public:
    /** Default constructor.
     * This constructor just initialize the pthread_mutex_t structure
     * with its defaults values, just like static initialization with PTHREAD_MUTEX_INITIALIZER
     */
    Condition();

    /** Destructor. */
     ~Condition();

    /** Function to get the access to the mutex.
     */
     void lock();

    /** Function to release the mutex.
     */
     void unlock();

    /** Function to be call from a thread to wait on
     * the condition. This function should be called after
     * acquiring the Condition lock...after that, will block
     * the threads until the condition be signaled.
     */
     void wait();

    /** Function to notify the condition was met.
     * Thread that can be waiting in the condition
     * will be awaked.
     */
     void signal();

    /** Send the signal to all waiting threads. */
     void broadcast();

}
;//end of class Condition

/** Class to synchronize several threads in some point of execution.
 * Threads is a way of distributing workload across
 * different workers to solve a problem faster. Nevertheless, sometimes
 * we need synchronization between threads to avoid undesired race
 * conditions and other problems. Here we are an implementation of a barrier
 * that allows putting all threads to wait at a given point until all of them
 * have reached such point and can continue working. Barriers are usually
 * available through pthreads system library. Nonetheless, sometimes it is not
 * so we have to implement it here.
 * @code
 * Mutex mutex;
 *
 * //Then in each thread to access the critical section:
 *  mutex.lock();
 *  //...Do critical section stuff
 *  mutex.unlock();
 *
   @endcode
 */
class Barrier
{
private:
    int needed; ///< How many threads should arraive to meet point
    int called; ///< How many threads already arrived
    Condition * condition; ///< Condition on which the threads are waiting

public:
    /** Constructor of the barrier to initialize the object.
     * You should pass the number of working threads that
     * you want to wait on the barrier. The internal counter
     * of the barrier will be initialized with numberOfThreads + 1
     * taking into account the main thread, so it need to wait
     * also in the barrier with the worker threads to all
     * can move on.
     * @code
     *  //For synchronize 10 threads created by a main thread
     *  //you can create the barrier from the main thread
     *  Barrier * barrier = new Barrier(10);
     *  //...
     *  //In the synchronization point
     *  barrier->wait();
     * @endcode
     * */
    Barrier(int numberOfThreads);

    /** Destructor to free all memory used */
    ~Barrier();

    /** Request to wait in this meet point.
     * For each thread calling this function the execution will
     * be paused untill all threads arrive this point.
     */
    void wait();

}
;//end of class Barrier

void * _threadMain(void * data);

/** Class for manage a group of threads performing one or several tasks.
 * This class is very useful when we have some function that can be executed
 * in parrallel by threads. The threads are created in the contructor of the object
 * and released in destructor. This way threads can execute different
 * functions at diffent moments and exit at the end of manager life. Also, the
 * wait() function allow in the main thread to wait until all threads have
 * finish working on a task and maybe then execute another one.
 * This class is supposed to be used only in the main thread.
 */
class ThreadManager
{
public:
    int threads; ///< number of working threads.
private:
    pthread_t * ids; ///< pthreads identifiers
    ThreadArgument * arguments; ///< Arguments passed to threads
    Barrier * barrier; ///< barrier for synchronized between tasks.
    /// Pointer to the function to work on,
    /// if null threads should exit
    ThreadFunction workFunction;
    bool started;
    void * workClass;

    /** Function to create threads structure and each thread
     * will be waiting to start working. Will be called on the first use.
     */
    void createThreads();

public:
    /** Set data for working threads.
     * If nThread = -1 then data is set for all threads.
     */
    void setData(void * data, int nThread = -1);

    /** Constructor, number of working threads should be supplied */
    ThreadManager(int numberOfThreads, void * workClass = NULL);

    /** Destructor, free memory and exit threads */
    ~ThreadManager();

    /** Function to start working in a task.
     * The function that you want to execute in parallel
     * by the working threads should be passed as argument.
     * If data is passed, then it is set to all threads.
     * Functions that can be executed by thread should by of the
     * type ThreadFunction, i.e., return void * and only
     * one argument of type ThreadArgument.
     * The call of this function will block the main thread
     * until all workers finish their job, if you dont want to block
     * use runAsync instead, and later can call wait for waiting
     * until threads are done.
     * @code
     *
     *  //Global variables, so it are visible in 'processSeveralImages()'
     *  ParallelTaskDistributor * td;
     *  //function to perform some operation
     *  //to N images executed in parellel
     *  void * processImages(ThreadArgument & data)
     *  {
     *      int thread_id = arg.thread_id;
     *
     *      size_t firstImage, lastImage;
     *      while (td->getTasks(firstImage, lastImage))
     *          for (int image = firstImage; image <= lastImage; ++image)
     *          {
     *              //...
     *              processOneImage(image);
     *              //...
     *          }
     *  }
     *
     *  int main()
     *  {
     *  //...
     *  //Distribute 1000 tasks in blocks of 100.
     *  td = new ThreadTaskDistributor(1000, 100);
     *  //Start 2 threads to work on it
     *  ThreadManager * tm = new ThreadManager(2);
     *  tm.run(processImages);
     *  //...
     *  //Same threads can work in other function
     *  tm.run(processVolumes);
     *  }
     * @endcode
     */
    void run(ThreadFunction function, void * data = NULL);

    /** Same as run but without blocking. */
    void runAsync(ThreadFunction function, void * data = NULL);

    /** Function that should be called to wait until all threads finished work */
    void wait();

    /** function to start running the threads.
     * Should be external and declared as friend */
    friend void * _threadMain(void * data);

    /** Get number of threads */
    int getNumberOfThreads() {return threads;}
}
;//end of class ThreadManager

/** Class to pass arguments to threads functions.
 * The argument passed can be obtained casting
 * the void * data received in the function.
 * @see ThreadManager
 */
class ThreadArgument
{
private:
    ThreadManager * manager;
public:
    int thread_id; ///< The thread id
    int threads;   ///< Number of threads
    void * workClass; ///< The class in wich threads will be working
    void * data; // Pointer to void *

    ThreadArgument();
    ThreadArgument(int id, ThreadManager * manager = NULL, void * data = NULL);

    friend class ThreadManager;
    friend void * _threadMain(void * data);
    int getNumberOfThreads() {return manager->getNumberOfThreads();}
};

/** This class distributes dynamically N tasks between parallel workers.
 * @ingroup ParallelLibrary
 * This class is a generalization of a common task in a parallel
 * environment of dynamically distribute N tasks between workers(threads or mpi proccess).
 * Each worker will ask for a group of tasks, proccess it and ask for more tasks
 * until there is not more task to process.
 *
 * This class is abstract and only serves as base for
 * concrete implementations, which will provides the
 * specific lock mechanisms and the way of distribution.
 */
class ParallelTaskDistributor
{
protected:
    //How many tasks give in each request
    size_t blockSize;
    //The number of tasks that have been assigned
    size_t assignedTasks;

public:
    //The total number of tasks to be distributed
    size_t numberOfTasks;
    /** Constructor.
     * The number of jobs and block size should be provided.
     */
    /** Constructor for Master node.
     */
    ParallelTaskDistributor(size_t nTasks, size_t bSize);

    /** Destructor.
     */
    virtual ~ParallelTaskDistributor()
    {}
    ;

    /** Restart the number of assigned tasks and distribution again.
     * This method should only be called in the main thread
     * before start distributing the tasks between the workers
     * threads.
     */
    void clear();

    /** Set the number of tasks assigned in each request */
    void setBlockSize(size_t bSize);

    /** Return the number of tasks assigned in each request */
    int getBlockSize() const;

    /** Gets parallel tasks.
     *  @ingroup ParallelJobHandler
     *  This function will be called by workers for asking tasks
     *  until there are not more tasks to process.
     *  Example:
     *  @code
     *  //...
     *  ParallelTaskDistributor * td = new ThreadTaskDistributor(1000, 100);
     *  //...
     *  //function to perform some operation
     *  //to N images executed in parellel
     *  void processSeveralImages()
     *  {
     *      size_t firstImage, lastImage;
     *      while (td->getTasks(firstImage, lastImage))
     *          for (size_t image = firstImage; image <= lastImage; ++image)
     *          {
     *              //...
     *              processOneImage(image);
     *              //...
     *          }
     *  }
     *  @endcode
     */
    bool getTasks(size_t &first, size_t &last); // False = no more jobs, true = more jobs
    /* This function set the number of completed tasks.
     * Usually this not need to be called. Its more useful
     * for restarting work, when usually the master detect
     * the number of tasks already done.
     */
    bool setAssignedTasks(size_t tasks);


protected:
    //Virtual functions that should be implemented in
    //subclasses, providing a mechanism of lock and
    //the specific way of distribute tasks.
    virtual void lock() = 0;
    virtual void unlock() = 0;
    virtual bool distribute(size_t &first, size_t &last) = 0;

}
;//class ParallelTaskDistributor

/** This class is a concrete implementation of ParallelTaskDistributor for POSIX threads.
 * It use mutex as the locking mechanism
 * and distributes tasks from 0 to numberOfTasks.
 */
class ThreadTaskDistributor: public ParallelTaskDistributor
{
public:
    ThreadTaskDistributor(size_t nTasks, size_t bSize):ParallelTaskDistributor(nTasks, bSize)
    {}
    virtual ~ThreadTaskDistributor()
    {}
    ;
protected:
    Mutex mutex; ///< Mutex to synchronize access to critical region
    virtual void lock();
    virtual void unlock();
    virtual bool distribute(size_t &first, size_t &last);
public:
    virtual void reset() { setAssignedTasks(0); };
}
;//end of class ThreadTaskDistributor

/** @name Old parallel stuff. */
/** Barrier structure */
//@{
typedef struct mybarrier_t
{
    /// How many threads should be awaited
    int needed;
    /// How many threads already arrived
    int called;
    /// Mutex to update this structure
    pthread_mutex_t mutex;
    /// Condition on which the threads are waiting
    pthread_cond_t cond;
}
barrier_t;

/** Barrier initialization */
int barrier_init(barrier_t *barrier, int needed);
/** Barrier destruction */
int barrier_destroy(barrier_t *barrier);
/** Wait at the barrier */
int barrier_wait(barrier_t *barrier);
//@}
//@}
#endif
