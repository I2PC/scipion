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

#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip> //for print pi with more decimals
#include "data/parallel_job_handler.h"

//Some useful macros
double PI25DT = 3.14159265358979323842643;
#define X2 (x*x) //((x-R)*(x-R))
#define Y2 (y*y) //((y-R)*(y-R))
#define R2 (R*R)
#define IS_MASTER (node == 0)

//mpi macros
#define TAG_WORK   0
#define TAG_STOP   1
#define TAG_WAIT   2

//For time tests
elapsedTime lockTime;
elapsedTime processingTime;
double  T = 0;
double T2 = 0;

int main(int argc, char **argv)
{
    ParallelJobHandler *jobHandler;
    char fName[L_tmpnam];

    long long int R = 10000;
    long long int blockSize = 1000;
    int node = 0, size = 0;
    long long int workBuffer[2];
    /** status after am MPI call */
    MPI_Status status;

    //MPI Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc > 1)
        R = atoll(argv[1]);

    if (argc > 2)
        blockSize = atoll(argv[2]);

    long long int first = -1, last = -1;
    long long int totalLocks = 1;
    long long int insideCounter = 0;
    bool moreJobs = true;
    float T = 0.;
    bool checkON = true;

    //Create the job handler in the master
    if (IS_MASTER)
    {
        fName[0] = '\0';
        jobHandler = new ParallelJobHandler(R, blockSize, fName);


        int finalizedWorkers = 0;

        while (finalizedWorkers < size - 1)
        {
            //wait for request form workers
            MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            jobHandler->getJobs(workBuffer[0], workBuffer[1]);
            if (workBuffer[0] == -1) //no more jobs
                finalizedWorkers++;
            //send work
            MPI_Send(workBuffer, 2, MPI_LONG_LONG_INT, status.MPI_SOURCE, TAG_WORK, MPI_COMM_WORLD);
        }

    }
    else
    {

        while (moreJobs)
        {
            lockTime.setStartTime();

            //any message from the master, is tag is TAG_STOP then stop
            MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(workBuffer, 2, MPI_LONG_LONG_INT, 0, TAG_WORK, MPI_COMM_WORLD, &status);

            lockTime.setEndTime(); //this is only for test the locking time
            T += lockTime.getElapsedTime();
            totalLocks++;

            first = workBuffer[0];
            last = workBuffer[1];
            moreJobs = (first != -1);

            if (moreJobs)
            {
                processingTime.setStartTime(); //start counting processing time

                for (long long int x = 0; x <= R; x++)//just for more work to do
                    for (long long int y = first; y <= last; y++)
                    {
                        if (X2 + Y2 <= R2)
                            insideCounter++;
                    }

                processingTime.setEndTime();
                T2 += processingTime.getElapsedTime();

                lockTime.setStartTime();
            }
        }
    }


    // Print out results subsequently
    for (int n = 0; n < size; n++)
    {
        if (n == node)
        {
            if(checkON && IS_MASTER)
            {
                if (!processingTime.saneInterval())
                    std::cout << "WARNING: increase job size (mpi_job_size)" <<std::endl
                    << "at present each block takes about " << processingTime.getElapsedTime()
                    << " seconds. At least it should take a few seconds, minutes will be even better"
                    << std::endl;
            }
            std::cout << "Node" << node
            << ": locks: " << totalLocks
            << " total locktime " << T
            << " total processingTime " << T2
            << " avg locktime " << (T/totalLocks)
            << " avg processingTime " << (T2/totalLocks)
            << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    //Send all counters to the master
    long long int totalInsideCounter;
    MPI_Reduce(&insideCounter, &totalInsideCounter, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    //Calculate PI on the master node
    if (IS_MASTER)
    {
        double myPI = (double)(totalInsideCounter*4)/R2;
        std::cout.precision(20);
        std::cout << "PI: " << std::fixed << myPI <<std::endl;
        delete jobHandler;
    }



    MPI_Finalize();

    return 0;
}
