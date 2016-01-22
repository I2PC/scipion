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
 * This program is a simple example of using the MPI parallel library
 * in which some needed MPI calls are encapsulated in the MpiNode class.
 * In the constructor the MPI initialization are done and the finalize
 * in destructor.
 *
 * In this example PI is calculated by sampling discrete
 * points of a quarter of a circle with radius R and aproximating
 * the area by the number of points inside the circle.
 *
 * Also the points to be process(parallel tasks) are dynamically distrubuted
 * between the MPI nodes. Two distributors are available, one base on
 * filesystem locking and another in which one of the MPI nodes will be
 * distributing the work.
 ***************************************************************************/

#include "parallel/xmipp_mpi.h"

//Some useful macros
double PI25DT = 3.14159265358979323842643;
#define X2 (x*x) //((x-R)*(x-R))
#define Y2 (y*y) //((y-R)*(y-R))
#define R2 (R*R)
#define IS_MASTER (node->rank == 0)

//For time tests
double T = 0;
double T2 = 0;

ParallelTaskDistributor * distributor;
MpiNode * node;

int main(int argc, char **argv)
{

    long long int R = 10000;
    long long int blockSize = 1000;
    /** status after am MPI call */

    //MPI Initialization
    node = new MpiNode(argc, argv);
    //MpiNode node(argc, argv);

    if (argc > 1)
        R = atoll(argv[1]);

    if (argc > 2)
        blockSize = atoll(argv[2]);

    longint first = -1, last = -1;
    longint totalLocks = 1;
    longint insideCounter = 0;
    bool moreJobs = true;
    float T = 0.;

    //MPI distributor using a node to distribute
    //distributor = new MpiTaskDistributor(R, blockSize, &node);
    //distributor = new MpiTaskDistributor(R, blockSize, node);
    //File distributor using a file
    distributor = new FileTaskDistributor(R, blockSize, node);

    std::cerr << "hello from node(" << node->rank << ")!!!" << std::endl;

    //All nodes working, there is not master!!!
    //Liberté, égalité, fraternité!!!
    while ((moreJobs = distributor->getTasks(first, last)))
    {
        std::cerr << "node:" << node->rank << " working from " << first
                << " to " << last << std::endl;
        for (long long int x = 0; x <= R; x++)//just for more work to do
            for (long long int y = first; y <= last; y++)
            {
                if (X2 + Y2 <= R2)
                    insideCounter++;
            }
    }

    delete distributor;
    MPI_Barrier(MPI_COMM_WORLD);

    // Print out results subsequently
    for (int n = 0; n < node->size; n++)
    {
        if (n == node->rank)
        {
            std::cout << "Node" << node->rank << ": locks: " << totalLocks
                    << " total locktime " << T << " total processingTime "
                    << T2 << " avg locktime " << (T / totalLocks)
                    << " avg processingTime " << (T2 / totalLocks) << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    //Send all counters to the master
    long long int totalInsideCounter;
    MPI_Reduce(&insideCounter, &totalInsideCounter, 1, MPI_LONG_LONG_INT,
            MPI_SUM, 0, MPI_COMM_WORLD);

    //Calculate PI on the master node
    if (node->isMaster())
    {
        double myPI = (double) (totalInsideCounter * 4) / R2;
        std::cout.precision(20);
        std::cout << "PI: " << std::fixed << myPI << std::endl;
    }

    delete node;
    return 0;
}
