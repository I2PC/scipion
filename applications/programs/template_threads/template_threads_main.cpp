/***************************************************************************
 *
 * Authors: J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
 *
 ****************************************************************************/


/****************************************************************************
* This program is a simple example of using the class ParallelJobHandler
* with MPI for calculating PI by sampling discrete
*  points of a quarter of a circle with radius R and aproximating
*  the area by the number of points inside.
***************************************************************************/


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <iomanip> //for print pi with more decimals
#include <data/xmipp_threads.h>
#include <data/xmipp_image.h>

using namespace std;

//Some useful macros
double PI25DT = 3.14159265358979323842643;
#define X2 (x*x) //((x-R)*(x-R))
#define Y2 (y*y) //((y-R)*(y-R))
#define R2 (R*R)

typedef char STRING[256];

//For time tests
double  T = 0;
double T2 = 0;

//Some global variables
Mutex mutex, taskMutex;
Barrier * barrier;
ParallelTaskDistributor * td;

int numberOfThreads = 1;

size_t R = 100000;
size_t blockSize = 1000;
size_t totalCounter = 0;

/*Here the number of jobs is the radius R that will be sampled */

void nothingFunctionNew(ThreadArgument &thArg)
{
    std::cerr << "Hello from thread " <<  thArg.thread_id;
}

/*Function to count the number of points that lies
 * inside a circle of radius R, just for calculating PI.
 * This is parallelized distributing some points for each threads
 */
void computePoints(int thread_id)
{
    size_t first = -1, last = -1;
    size_t insideCounter = 0;
   //ThreadManager tm(2);

    while (td->getTasks(first, last))
    {
        std::cerr << "th" << thread_id << ": working from " << first << " to " << last <<std::endl;
        for (size_t x = 0; x <= R; x++)//just for more work to do
            for (size_t y = first; y <= last; y++)
            {
                if (X2 + Y2 <= R2)
                    insideCounter++;

            }
        std::cerr << "th" << thread_id << ": asking more jobs" << std::endl;

    }
    std::cerr << "Updating counter thread " << thread_id << std::endl;
    //Lock for update the total counter
    mutex.lock();
    totalCounter += insideCounter;
    mutex.unlock();
        std::cerr << "Leaving counter thread " << thread_id << std::endl;
}

void threadFunctionNew(ThreadArgument &thArg)
{
    computePoints(thArg.thread_id);
}

void * threadFunctionOld(void * data)
{
    long id = (long)data;
    computePoints((int)id);
}

int main(int argc, char **argv)
{

    if (argc > 1) //read number of threads as first argument
        numberOfThreads = atoi(argv[1]);

    if (argc > 2) //read R as second argument
        R = atoll(argv[2]);

    if (argc > 3) //read blockSize
        blockSize =  atoll(argv[3]);

    bool oldMethod = false;
    if (argc > 4) //old or new method
        oldMethod = true;

    ThreadManager * thMgr = NULL;//new  ThreadManager(numberOfThreads);

    for (int i = 0; i < 1; ++i)
    {
        totalCounter = 0;
        //Create the job handler to distribute jobs
        //jobHandler = new ParallelJobHandler(R, blockSize);
        td = new ThreadTaskDistributor(R, blockSize);

        if (!oldMethod)
        {
            //Create threads to start working
            thMgr = new ThreadManager(numberOfThreads);
            std::cerr << "Creating ThreadManager...." << std::endl;
            thMgr->run(threadFunctionNew);
        }
        else
        {
            /////////////////////////////////////////////////////////////////////////
            /* Following is the old way to do the same thing with posix functions*/
            pthread_t threads[numberOfThreads];
            int rc;
            long t;
            for(t = 0; t < numberOfThreads; t++)
            {
                printf("In main: creating thread %ld\n", t);
                rc = pthread_create(&threads[t], NULL, threadFunctionOld, (void *)t);
                if (rc)
                {
                    printf("ERROR; return code from pthread_create() is %d\n", rc);
                    exit(-1);
                }
            }

            for (t = 0; t < numberOfThreads; ++t)
            {
                rc = pthread_join(threads[t], NULL);
            }

        }
        ////////////////////////////////////////////////////////////////////////
        //Calculate PI based on the number of points inside circle
        double myPI = (double)(totalCounter * 4) / R2;
        std::cout.precision(20);
        std::cout << "PI: " << std::fixed << myPI <<std::endl;

        //Terminate threads and free memory
        delete td;
        delete thMgr;
    }
    //delete thMgr;
    return 0;
}
