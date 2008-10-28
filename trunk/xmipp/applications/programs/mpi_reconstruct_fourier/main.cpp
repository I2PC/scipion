/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

//#include "mpi_run.h"

#include <data/args.h>
#include <reconstruction/reconstruct_fourier.h>
#include <data/header.h>
#include <reconstruction/projection.h>
#include <cstring>
#include <cstdlib>
#include <data/funcs.h>
#include <data/matrix2d.h>
#include <sys/time.h>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>  

#define TAG_WORKFORWORKER   0
#define TAG_STOP   1
#define TAG_TRANSFER 2
#define TAG_FREEWORKER   3

class Prog_mpi_RecFourier_prm:Prog_RecFourier_prm
{
    public:

	FileName          fn_img;            

	/** Fourier transform size for volumes **/
	long int sizeout;
		
        /** Number of Procesors **/
        int nProcs;
        
        /** Dvide the job in this number block with this number of images */
        int mpi_job_size;
		
        /** Number of independent MPI jobs **/
        int numberOfJobs;
        
        /** computing node number. Master=0 */
        int rank;
		
        /** status after am MPI call */
        MPI_Status status;
		
        /** verbose mode on/off.  */
        bool verbose;
		
	/*  constructor ------------------------------------------------------- */
	Prog_mpi_RecFourier_prm()
	{
		//parent class constructor will be called by deault without parameters
		MPI_Comm_size(MPI_COMM_WORLD, &(nProcs));
		MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
		//if (nProcs < 2)
		//	error_exit("This program cannot be executed in a single working node");
		//Blocks until all process have reached this routine.
		//very likelly this is
		MPI_Barrier(MPI_COMM_WORLD);
	}


	/* Read parameters --------------------------------------------------------- */
	void read(int argc, char **argv)
	{
		Prog_RecFourier_prm::read(argc,argv);
		mpi_job_size=textToInteger(getParameter(argc,argv,"-mpi_job_size","5"));
	}

	/* Usage ------------------------------------------------------------------- */
	void usage()
	{
		Prog_RecFourier_prm::usage();
		std::cerr << " [ -mpi_job_size default=5]    : Number of images sent to a cpu in a single job \n";
		std::cerr << "                                  10 may be a good value";
		std::cerr << "                                  if  -1 the computer will put the maximum";
		std::cerr << "                                  posible value that may not be the best option";
	}


	/* Show -------------------------------------------------------------------- */
	void show()
	{
		Prog_RecFourier_prm::show();
		std::cerr << " Size of mpi jobs " << mpi_job_size <<std::endl;
	}

	/* Pre Run PreRun for all nodes but not for all works */
	void preRun()
	{
		if (rank == 0) 
		{
			show();
		}
		
		produce_Side_info();

		SF.go_beginning();
		
		//leer sel file / dividir por mpi_job_size 
		numberOfJobs=ceil(SF.ImgNo()/mpi_job_size);

		//only one node will write in the console
		if (rank == 0 )
		{
//#define DEBUG
#ifdef DEBUG
			std::cerr << "numberOfJobs: " << numberOfJobs << std::endl <<std::endl;
#endif
#undef DEBUG
		}
		else
			verbose=false;

	}
	/* Run --------------------------------------------------------------------- */
	void run()
	{   
        	struct timeval start_time, end_time;
  		long int total_usecs;
		double total_time;

		double * fourierVolume = (double *)VoutFourier.data;
		double * fourierWeights = FourierWeights.data;

		sizeout = MULTIDIM_SIZE(FourierWeights);

		if (rank == 0)
		{	
		 	gettimeofday(&start_time,NULL);

			if( verb )
				init_progress_bar(numberOfJobs);

			int stopTagsSent =0;
			for (int i=0;i<numberOfJobs;)
			{
				MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
							 MPI_COMM_WORLD, &status);
				
				if( status.MPI_TAG != TAG_FREEWORKER )
				{
					cout << "Unexpected TAG, please contact developers " << endl;
					exit(-1);
				}
				
				MPI_Send(&i,
					 1,
					 MPI_INT,
					 status.MPI_SOURCE,
					 TAG_WORKFORWORKER,
					 MPI_COMM_WORLD);
				
				i++; //increase job number 
				if( verb)
					progress_bar(i);      
       			}

  			gettimeofday(&end_time,NULL);
  
  			total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
  			total_time=(double)total_usecs/(double)1000000;		

			std::cout << std::flush;
			std::cout << "\nProcessing time: " << total_time << " secs." << std::endl;

			int currentSource;

			gettimeofday(&start_time,NULL);

			// Start collecting results
			while(1)
			{
				for( int works=0; works < (nProcs-1); works++)
				{ 
					MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
						 MPI_COMM_WORLD, &status);
					currentSource = status.MPI_SOURCE;

					// Signal the worker to start sending back the results
					MPI_Send(0, 0, MPI_INT, currentSource, TAG_TRANSFER, MPI_COMM_WORLD);
				}
						
				MPI_Reduce( MPI_IN_PLACE, fourierVolume, 2 * sizeout, MPI_DOUBLE, MPI_SUM,
				0, MPI_COMM_WORLD);
				
				MPI_Reduce( MPI_IN_PLACE, fourierWeights, sizeout, MPI_DOUBLE, MPI_SUM,
				0, MPI_COMM_WORLD);
				
				gettimeofday(&end_time,NULL);

  				total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
  				total_time=(double)total_usecs/(double)1000000;		

				std::cout << "Transfers time: " << total_time << " secs." << std::endl;
				
				// Normalize global volume and store data
				gettimeofday(&start_time,NULL);
				finishComputations();
    				gettimeofday(&end_time,NULL);

  				total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
  				total_time=(double)total_usecs/(double)1000000;		

				std::cout << "Weighting time: " << total_time << " secs." << std::endl;
				
				std::cout << "Execution completed successfully\n" << std::endl;
				break;
			}   
		}
		else
		{
			// Select only relevant part of selfile for this rank
			// job number
			// job size
			// aux variable
			while (1)
			{
				int jobNumber;
				//I am free
				MPI_Send(0, 0, MPI_INT, 0, TAG_FREEWORKER, MPI_COMM_WORLD);
				
				MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				if (status.MPI_TAG == TAG_TRANSFER)
            			{
					//If I  do not read this tag
					//master will no further process
					//a posibility is a non-blocking send
					MPI_Recv(0, 0, MPI_INT, 0, TAG_TRANSFER,
							 MPI_COMM_WORLD, &status);
#ifdef DEBUG
					std::cerr << "Wr" << rank << " " << "TAG_STOP" << std::endl;
#endif
					MPI_Reduce( fourierVolume, NULL, 2*sizeout, MPI_DOUBLE, MPI_SUM,
					0, MPI_COMM_WORLD);
			
					MPI_Reduce( fourierWeights, NULL, sizeout, MPI_DOUBLE, MPI_SUM,
					0, MPI_COMM_WORLD);
			
					break;
            			}
				else if (status.MPI_TAG == TAG_WORKFORWORKER)
					//there is still some work to be done    
            			{		
					//get the jobs number
					MPI_Recv(&jobNumber, 1, MPI_INT, 0, TAG_WORKFORWORKER, MPI_COMM_WORLD, &status);
					int min_i, max_i;

					min_i = jobNumber*mpi_job_size;
					max_i = min_i + mpi_job_size;

					// Process all images
					for( int i = min_i ; i < max_i ; i ++ )
					{
						// Check whether all projections have already 
						// been processed. If so, break loop
						if( i > SF.ImgNo() )
							break;

						fn_img = SF.get_file_number(i);

						processImage(fn_img );
					}
            			}
				else
				{
					std::cerr << "3) Received unknown TAG I quit" << std::endl;
					exit(0);
				}           
			}
		}
		MPI_Finalize();
	}

	/* a short function to print a message and exit */
	void error_exit(char * msg)
	{
		fprintf(stderr, "%s", msg);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

};

int main(int argc, char *argv[])
{
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        fprintf(stderr, "MPI initialization error\n");
        exit(EXIT_FAILURE);
    }

    Prog_mpi_RecFourier_prm prm;
    try
    {
        prm.read(argc, argv);
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE;
        prm.usage();
        exit(1);
    }
		
    try
    {
        prm.preRun();
        prm.run();
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE;
        exit(1);
    }
    
    exit(0);
}


