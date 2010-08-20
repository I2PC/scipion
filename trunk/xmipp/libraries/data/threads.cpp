#include "threads.h"
#include <stdio.h>
#include <iostream>

// ================= MUTEX ==========================
Mutex::Mutex()
{
    pthread_mutex_init(&mutex, NULL);
}

void Mutex::lock()
{
    pthread_mutex_lock(&mutex);
}

Mutex::~Mutex()
{
    pthread_mutex_destroy(&mutex);
}

void Mutex::unlock()
{
    pthread_mutex_unlock(&mutex);
}

// ================= BARRIER ==========================
Barrier::Barrier(int numberOfThreads)
{
    needed = numberOfThreads + 1;
    called = 0;
    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond, NULL);
}

Barrier::~Barrier()
{
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&cond);
}

void Barrier::wait()
{
    pthread_mutex_lock(&mutex);
    ++called;
    if (called == needed)
    {
        called = 0;
        pthread_cond_broadcast(&cond);
    }
    else
    {
        pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
}

// ================= THREAD MANAGER =======================

ThreadArgument::ThreadArgument()
{
    thread_id = -1;
    manager = NULL;
    data = NULL;
}

ThreadArgument::ThreadArgument(int id, ThreadManager * manager, void * data)
{
    this->thread_id = id;
    this->manager = manager;
    this->data = data;
}

void * _threadMain(void * data)
{
    ThreadArgument * thArg = (ThreadArgument*) data;
    ThreadManager * thMgr = thArg->manager;

    while (true)
    {
        //Wait for start working or leave
        thMgr->barrier->wait();
        //After awaked check what to do
        if (thMgr->workFunction != NULL)
        {
            thMgr->workFunction(*thArg);
            thMgr->barrier->wait(); //wait for finish together
        }
        else //exit thread
        {
            pthread_exit(NULL);
            return NULL;
        }
    }
}

ThreadManager::ThreadManager(int numberOfThreads, void * workClass)
{
    threads = numberOfThreads;
    barrier = new Barrier(threads);
    workFunction = NULL;
    ids = new pthread_t[threads];
    arguments = new ThreadArgument[threads];
    //Create threads
    int result;
    for (int i = 0; i < numberOfThreads; ++i)
    {
        arguments[i].thread_id = i;
        arguments[i].manager = this;
        arguments[i].data = NULL;
        arguments[i].workClass = workClass;

        result = pthread_create(ids + i, NULL, _threadMain, (void*) (arguments
                + i));

        if (result != 0)
        {
            std::cerr << "ThreadManager Constructor: can't create threads.";
            exit(1);
        }
    }

}

ThreadManager::~ThreadManager()
{
    //Destroy the threads
    workFunction = NULL;
    barrier->wait();

    delete barrier;
    delete[] ids;
    delete[] arguments;
}

void ThreadManager::run(ThreadFunction function)
{
    runAsync(function);
    //Wait on barrier to wait for threads finish
    wait();
}

void ThreadManager::runAsync(ThreadFunction function)
{
    workFunction = function;
    //Wait on barrier to threads starts working
    barrier->wait();
}

void ThreadManager::wait()
{
    barrier->wait();
}

// =================== TASK_DISTRIBUTOR ============================

ParallelTaskDistributor::ParallelTaskDistributor(longint nTasks, longint bSize)
{
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

void ParallelTaskDistributor::setBlockSize(longint bSize)
{
    lock();
    blockSize = bSize;
    unlock();
}

int ParallelTaskDistributor::getBlockSize() const
{
    return blockSize;
}

bool ParallelTaskDistributor::getTasks(longint &first, longint &last)
{
    lock();
    bool result = distribute(first, last);
    unlock();
    return result;
}

void ThreadTaskDistributor::lock()
{
    mutex.lock();
}

void ThreadTaskDistributor::unlock()
{
    mutex.unlock();
}

bool ThreadTaskDistributor::distribute(longint &first, longint &last)
{
    bool result = true;
    first = last = -1;
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





