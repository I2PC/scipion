#include "xmipp_threads.h"
#include "xmipp_error.h"
#include <stdio.h>
#include <iostream>


// ================= MUTEX ==========================

Mutex::Mutex()
{
    pthread_mutex_init(&mutex, NULL);
}

Mutex::~Mutex()
{
    pthread_mutex_destroy(&mutex);
}

void Mutex::lock()
{
    pthread_mutex_lock(&mutex);
}

void Mutex::unlock()
{
    pthread_mutex_unlock(&mutex);
}

// ================= CONDITION ==========================

Condition::Condition()
{
    mutex = new Mutex();
    pthread_cond_init(&cond, NULL);
}

Condition::~Condition()
{
    delete mutex;
    pthread_cond_destroy(&cond);
}

void Condition::lock()
{
    mutex->lock();
}

void Condition::unlock()
{
    mutex->unlock();
}

void Condition::wait()
{
    pthread_cond_wait(&cond, &(mutex->mutex));
}

void Condition::signal()
{
    pthread_cond_signal(&cond);
}

void Condition::broadcast()
{
    pthread_cond_broadcast(&cond);
}

// ================= BARRIER ==========================

Barrier::Barrier(int numberOfThreads)
{
    needed = numberOfThreads;
    called = 0;
    condition = new Condition();
}

Barrier::~Barrier()
{
    delete condition;
}

void Barrier::wait()
{
    condition->lock();
    ++called;
    if (called == needed)
    {
        called = 0;
        condition->broadcast();
    }
    else
        condition->wait();
    condition->unlock();
}


// ================= THREAD =======================

Thread::Thread()
{
}

Thread::~Thread()
{
    pthread_join(thId, NULL);
}

void Thread::start()
{
    int result = pthread_create(&thId, NULL, _singleThreadMain, (void*)this);

    if (result != 0)
    {
        std::cerr << "Thread: can't start thread." << std::endl;
        exit(1);
    }
}

void * _singleThreadMain(void * data){
  Thread * thread = (Thread*) data;
  thread->run();
  return NULL;
}

// ================= THREAD MANAGER =======================

ThreadArgument::ThreadArgument()
{
    thread_id = -1;
    threads = -1;
    manager = NULL;
    data = NULL;
    workClass = NULL;
}

ThreadArgument::ThreadArgument(int id, ThreadManager * manager, void * data)
{
    this->thread_id = id;
    this->threads = manager->threads;
    this->manager = manager;
    this->data = data;
    this->workClass = NULL;
}


void * _threadMain(void * data)
{
    ThreadArgument * thArg = (ThreadArgument*) data;
    ThreadManager * thMgr = thArg->manager;

    while (true)
    {
        //Wait for start working or leave
        thMgr->wait();
        //After awaked check what to do
        if (thMgr->workFunction != NULL)
        {
            try
            {
                thMgr->workFunction(*thArg);
                thMgr->wait(); //wait for finish together
            }
            catch (XmippError &xe)
            {
                std::cerr << xe << std::endl
                << "In thread " << thArg->thread_id << std::endl;
                pthread_exit(NULL);
            }
        }
        else //exit thread
        {
            pthread_exit(NULL);
        }
    }
}

ThreadManager::ThreadManager(int numberOfThreads, void * workClass)
{
    threads = numberOfThreads;
    barrier = new Barrier(threads + 1);
    workFunction = NULL;
    ids = new pthread_t[threads];
    arguments = new ThreadArgument[threads];
    started = false;
    this->workClass = workClass;
}

void ThreadManager::setData(void * data, int idxThread)
{
    if (idxThread == -1)
        for (int i = 0; i < threads; ++i)
            arguments[i].data = data;
    else
        arguments[idxThread].data = data;
}

void ThreadManager::createThreads()
{
    //Create threads
    int result;

    for (int i = 0; i < threads; ++i)
    {
        arguments[i].thread_id = i;
        arguments[i].threads = threads;
        arguments[i].manager = this;
        arguments[i].workClass = workClass;

        result = pthread_create(ids + i, NULL, _threadMain, (void*) (arguments + i));

        if (result != 0)
        {
            std::cerr << "ThreadManager: can't create threads." << std::endl;
            exit(1);
        }
    }
    started = true;
}

ThreadManager::~ThreadManager()
{
    //Destroy the threads
    workFunction = NULL;
    if (started)
    {
        wait();
        for (int i = 0; i < threads; ++i)
            pthread_join(ids[i], NULL);
    }

    delete barrier;
    delete[] ids;
    delete[] arguments;
}

void ThreadManager::run(ThreadFunction function, void * data)
{
    runAsync(function, data);
    //Wait on barrier to wait for threads finish
    wait();
}

void ThreadManager::runAsync(ThreadFunction function, void * data)
{
    if (data != NULL)
        setData(data);
    workFunction = function;
    if (!started)
        createThreads();
    //Wait on barrier to threads starts working
    wait();
}

void ThreadManager::wait()
{
    barrier->wait();
}

// =================== TASK_DISTRIBUTOR ============================

ParallelTaskDistributor::ParallelTaskDistributor(size_t nTasks, size_t bSize)
{
    if (!(nTasks && bSize && bSize <= nTasks))
        REPORT_ERROR(ERR_ARG_INCORRECT, "nTasks and bSize should be > 0, also bSize <= nTasks");

    numberOfTasks = nTasks;
    blockSize = bSize;
    assignedTasks = 0;
}

void ParallelTaskDistributor::clear()
{
    lock();
    assignedTasks = 0;
    unlock();
}

void ParallelTaskDistributor::setBlockSize(size_t bSize)
{
    lock();
    blockSize = bSize;
    unlock();
}

int ParallelTaskDistributor::getBlockSize() const
{
    return blockSize;
}

bool ParallelTaskDistributor::getTasks(size_t &first, size_t &last)
{
    lock();
    bool result = distribute(first, last);
    unlock();
    return result;
}

bool ParallelTaskDistributor::setAssignedTasks(size_t tasks)
{
    if (tasks < 0 || tasks >= numberOfTasks)
        return false;
    lock();
    assignedTasks = tasks;
    unlock();
    return true;
}

void ThreadTaskDistributor::lock()
{
    mutex.lock();
}

void ThreadTaskDistributor::unlock()
{
    mutex.unlock();
}

bool ThreadTaskDistributor::distribute(size_t &first, size_t &last)
{
    bool result = true;
    first = last = 0;
    if (assignedTasks >= numberOfTasks)
    {
        result = false;
    }
    else
    {
        first = assignedTasks;
        assignedTasks
        = (assignedTasks + blockSize < numberOfTasks) ? (assignedTasks
                + blockSize) : numberOfTasks;
        last = assignedTasks - 1;
    }
    return result;
}

// =================== OLD THREADS IMPLEMENTATION ============================
int barrier_init(barrier_t *barrier,int needed)
{
    barrier->needed = needed;
    barrier->called = 0;
    pthread_mutex_init(&barrier->mutex, NULL);
    pthread_cond_init(&barrier->cond, NULL);
    return 0;
}

int barrier_destroy(barrier_t *barrier)
{
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->cond);
    return 0;
}

int barrier_wait(barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->mutex);
    barrier->called++;
    if (barrier->called == barrier->needed)
    {
        barrier->called = 0;
        pthread_cond_broadcast(&barrier->cond);
    }
    else
    {
        pthread_cond_wait(&barrier->cond,&barrier->mutex);
    }
    pthread_mutex_unlock(&barrier->mutex);
    return 0;
}





