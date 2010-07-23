//
//------------------------------------------------*/

#include <mpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include "data/parallel_job_handler.h"
//#include "funcs.h"

double PI25DT = 3.14159265358979323842643;

//For time tests
struct timeval start_time, end_time,end_time2;
float T = 0;
float T2 = 0;

float elapsed(struct timeval start, struct timeval end)
{
    return  (float) (end.tv_sec  - start.tv_sec ) +
            ((float) (end.tv_usec - start.tv_usec)/1000000);
}




int main(int argc, char **argv)
{
    ParallelJobHandler *jobHandler;
    char fName[L_tmpnam];

    long long int number_of_intervals = 10000000;
    long long int blockSize = 10000;

    int node = 0;

    //MPI Initialization
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    //Random number generator initialization
    srand ( time(NULL)*node );

    //Read arguments if Master
    if (node == 0)
    {
        if (argc > 1)
            number_of_intervals = atoll(argv[1]);
        if (argc > 2)
            blockSize = atoll(argv[2]);

        fName[0] = '\0';
        jobHandler = new ParallelJobHandler(number_of_intervals, blockSize, fName);
    }

    //Send the temporary filename to the slaves nodes
    MPI_Bcast(fName, L_tmpnam, MPI_CHAR, 0, MPI_COMM_WORLD);

    //MPI_Barrier(MPI_COMM_WORLD);

    if (node != 0)
        jobHandler = new ParallelJobHandler(fName);

    long long int first = -1, last = -1;
    long long int totalLocks = 0;
    long long int insideCounter = 0;
    bool moreJobs = true;
    float T = 0.;

    while (moreJobs)
    {
        gettimeofday(&start_time, 0);

        moreJobs = jobHandler->getJobs(first, last);

        gettimeofday(&end_time, 0);

#ifdef DEBUG_P
        if (moreJobs)
            std::cerr << "Node" << node <<" from " << first << " to " << last <<std::endl;
        else
            std::cerr << "Node" << node <<" no more jobs "<<std::endl;
#endif
        T += elapsed(start_time, end_time);
        totalLocks++;

        for (long long int rr = 0; rr < 1000; rr++)//just for more work to do
            for (long long int r = first; r <= last; r++)
            {
                double x = (double)rand() / RAND_MAX;
                double y = (double)rand() / RAND_MAX;
                double distance2 = x*x + y*y;
                if (distance2 <= 1)
                    insideCounter++;
            }
        //int r = rand();
        //usleep((r % 5000000) + (2 - 1)*4000000);
        gettimeofday(&end_time2, 0);
        T2 += elapsed(start_time, end_time2);

    }

    long long int totalInsideCounter;
    std::cout << "Node" << node
              << ": locks: " << totalLocks
              << " total time " << T
              << " total time2 " << T2
              << " avg time " << (T/totalLocks)
              << " avg time2 " << (T2/totalLocks)
              << std::endl;
    MPI_Reduce(&insideCounter, &totalInsideCounter, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (node == 0)
    {
        double myPI = (double)(totalInsideCounter*4)/(1000*number_of_intervals);
        std::cout << "PI: " << myPI <<std::endl;
    }
    //std::cerr << "th=  delta= "<< node << " " << T <<std::endl;

    delete jobHandler;

    MPI_Finalize();

    return 0;
}
