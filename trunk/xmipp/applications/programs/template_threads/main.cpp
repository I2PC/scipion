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
#include <iomanip> //for print pi with more decimals
#include "data/threads.h"

using namespace std;

//Some useful macros
double PI25DT = 3.14159265358979323842643;
#define X2 (x*x) //((x-R)*(x-R))
#define Y2 (y*y) //((y-R)*(y-R))
#define R2 (R*R)
#define TH_EXIT 0
#define TH_CHECK_POINTS 1

typedef char STRING[256];

//For time tests
double  T = 0;
double T2 = 0;

//Some global variables
Mutex mutex, taskMutex;
Barrier * barrier;
ParallelTaskDistributor * td;

int numberOfThreads = 1;
pthread_t * th_structs;
int * th_ids;
int threadTask = -1;
long long int R = 100000;
long long int blockSize = 1000;
long long int totalCounter = 0;

long long int assignedJobs = 0;
long long int &numberOfJobs = R;

void computePoints(ThreadArgument &thArg)
{
    //ThreadArgument * thArg = (ThreadArgument*)data;
    int thread_id = thArg.thread_id;

    long long int first = -1, last = -1;
    long long int insideCounter = 0;

    while (td->getTasks(first, last))
    {
        std::cerr << "th" << thread_id << ": working from " << first << " to " << last <<std::endl;
        for (long long int x = 0; x <= R; x++)//just for more work to do
            for (long long int y = first; y <= last; y++)
            {
                if (X2 + Y2 <= R2)
                    insideCounter++;
            }

    }
    //Lock for update the total counter
    mutex.lock();
    totalCounter += insideCounter;
    mutex.unlock();
}

int main(int argc, char **argv)
{

    if (argc > 1) //read number of threads as first argument
        numberOfThreads = atoi(argv[1]);

    if (argc > 2) //read R as second argument
        R = atoll(argv[2]);

    if (argc > 3) //read blockSize
        blockSize =  atoll(argv[3]);

    //Create the job handler to distribute jobs
    //jobHandler = new ParallelJobHandler(R, blockSize);
    td = new ThreadTaskDistributor(R, blockSize);

    //Create threads to start working
    //createThreads();
    ThreadManager * thMgr = new ThreadManager(numberOfThreads);
    thMgr->run(computePoints);
    //Terminate threads and free memory
    delete td;
    delete thMgr;

    //Calculate PI based on the number of points inside circle
    double myPI = (double)(totalCounter * 4) / R2;
    std::cout.precision(20);
    std::cout << "PI: " << std::fixed << myPI <<std::endl;

    return 0;
}
