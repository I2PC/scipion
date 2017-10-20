/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao (jrbcast@ace.ual.es)
 *         Roberto Marabini (roberto@cnb.csic.es)
 *         Vahid Abrishami (vabrishamoi@cnb.csic.es)
 *         David Strelak (davidstrelak@gmail.com)
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
#include "mpi_reconstruct_fourier_accel.h"

/*  constructor ------------------------------------------------------- */
ProgMPIRecFourierAccel::ProgMPIRecFourierAccel(int argc, char *argv[])
{
    this->read(argc, argv);
}

/* constructor providing an MpiNode
 * this is useful for using this programs from others
 */
ProgMPIRecFourierAccel::ProgMPIRecFourierAccel(MpiNode * node)
{
    this->setNode(node);
}

/* Special way of reading to sync all nodes */
void ProgMPIRecFourierAccel::read(int argc, char** argv)
{
    XmippMpiProgram::read(argc, argv);
    ProgRecFourierAccel::read(argc, (const char **)argv);
}

/* Usage ------------------------------------------------------------------- */
void ProgMPIRecFourierAccel::defineParams()
{
	ProgRecFourierAccel::defineParams();
    addParamsLine("  [--mpi_job_size <size=100>]    : Number of images sent to a cpu in a single job ");
    addParamsLine("                                : 100 may be a good value");
    addParamsLine("                                : if  -1 the computer will put the maximum");
    addParamsLine("                                : posible value that may not be the best option");
}

/* Read parameters --------------------------------------------------------- */
void ProgMPIRecFourierAccel::readParams()
{
	ProgRecFourierAccel::readParams();
    mpi_job_size=getIntParam("--mpi_job_size");
}

/* Pre Run PreRun for all nodes but not for all works */
void ProgMPIRecFourierAccel::preRun()
{
    if (nProcs < 2)
        REPORT_ERROR(ERR_ARG_INCORRECT,"This program cannot be executed in a single working node");

    if (node->isMaster())
    {
        show();
        SF.read(fn_in);
        //Send verbose level to node 1
        MPI_Send(&verbose, 1, MPI_INT, 1, TAG_SETVERBOSE, MPI_COMM_WORLD);
    }
    else
    {
        produceSideinfo();
        SF.firstObject();
    }

    //read projection file / divide by mpi_job_size
    numberOfJobs=(size_t)ceil((double)SF.size()/mpi_job_size);

    //only one node will write in the console
    if (node->rank == 1 )
    {
        // Get verbose status
        MPI_Recv(&verbose, 1, MPI_INT, 0, TAG_SETVERBOSE, MPI_COMM_WORLD, &status);
        //#define DEBUG
#ifdef DEBUG

        std::cerr << "SF.ImgNo() mpi_job_size "
        << SF.ImgNo() << " "
        << mpi_job_size
        << std::endl;
        std::cerr << "numberOfJobs: " << numberOfJobs << std::endl <<std::endl;
#endif
#undef DEBUGDOUBLE

    }
}
/* Run --------------------------------------------------------------------- */
void ProgMPIRecFourierAccel::run()
{
    preRun();
    struct timeval start_time, end_time;
    MPI_Group  orig_group, new_group;
    MPI_Comm   new_comm;
    long int total_usecs;
    double total_time_processing=0., total_time_communicating=0., total_time;
    int * ranks;

    // Real workers, rank=0 is the master, does not work
    nProcs = nProcs - 1;

    if( numberOfJobs < nProcs )
    {
        if( node->isMaster() )
        {
            std::cerr << "\nReducing the number of MPI workers from " <<
            nProcs << " to " <<
            numberOfJobs << std::endl;
        }

        nProcs = numberOfJobs;

        // Unused nodes are removed from the MPI communicator
        node->active = (node->rank <= numberOfJobs);
    }

    // Generate a new group to do all reduce without the master
    ranks = new int [nProcs];
    for (int i=0;i<(int)nProcs;i++)
        ranks[i]=i+1;
    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    MPI_Group_incl(orig_group, nProcs, ranks, &new_group);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);

	if (node->isMaster())
	{
		gettimeofday(&start_time,NULL);
		std::cerr<<"Computing volume"<<std::endl;
		if ( verbose )
			init_progress_bar(numberOfJobs);

		for (size_t i=0;i<numberOfJobs;i++)
		{
//          #define DEBUG
#ifdef DEBUG
		std::cerr << "master-recv  i=" << i << std::endl;
		std::cerr << "numberOfJobs: " << numberOfJobs << std::endl <<std::endl;
#endif
#undef DEBUG

			MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
					 MPI_COMM_WORLD, &status);

			if ( status.MPI_TAG != TAG_FREEWORKER )
				REPORT_ERROR(ERR_ARG_INCORRECT,"Unexpected TAG, please contact developers");

		//	#define DEBUG
#ifdef DEBUG
			std::cerr << "master-send i=" << i << std::endl;
#endif
#undef DEBUG
			MPI_Send(&i,
					 1,
					 MPI_INT,
					 status.MPI_SOURCE,
					 TAG_WORKFORWORKER,
					 MPI_COMM_WORLD);

			if (verbose)
				progress_bar(i);
		}

		// Wait for all processes to finish processing current jobs
		// so time statistics are correct
		for ( size_t i = 1 ; i <= nProcs ; i ++ )
		{
			MPI_Recv(0,
					 0,
					 MPI_INT,
					 MPI_ANY_SOURCE,
					 TAG_FREEWORKER,
					 MPI_COMM_WORLD,
					 &status);
		}

		// update progress
		gettimeofday(&end_time,NULL);
		total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
		total_time_processing += ((double)total_usecs/(double)1000000);

		gettimeofday(&start_time,NULL);
		// Start collecting results
		for ( size_t i = 1 ; i <= nProcs ; i ++ )
		{
			MPI_Send(0,
					 0,
					 MPI_INT,
					 i,
					 TAG_TRANSFER,
					 MPI_COMM_WORLD );
		}

		gettimeofday(&end_time,NULL);
		total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
		total_time_communicating += ((double)total_usecs/(double)1000000);

		if (verbose > 0)
		{
			std::cout << "\n\nProcessing time: " << total_time_processing << " secs." << std::endl;
			std::cout << "Transfers time: " << total_time_communicating << " secs." << std::endl;
			std::cout << "Execution completed successfully"<< std::endl;
		}
	}
	else if( node->active ) {
		createLoadingThread();
		while (1)
		{
			int jobNumber;
			//#define DEBUG
#ifdef DEBUG
			std::cerr << "slave-send TAG_FREEWORKER rank=" << node->rank << std::endl;
#endif
     #undef DEBUG
			//I am free
			MPI_Send(0, 0, MPI_INT, 0, TAG_FREEWORKER, MPI_COMM_WORLD);
			MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if (status.MPI_TAG == TAG_TRANSFER)
			{
				// reduce the data
				mirrorAndCropTempSpaces();

				//If I  do not read this tag
				//master will no further process
				MPI_Recv(0, 0, MPI_INT, 0, TAG_TRANSFER, MPI_COMM_WORLD, &status);
#ifdef DEBUG
				std::cerr << "Wr" << node->rank << " " << "TAG_STOP" << std::endl;
#endif
				for(int z = 0; z <= maxVolumeIndexYZ; z++) {
					for(int y = 0; y <= maxVolumeIndexYZ; y++) {
						if (node->rank == 1) {
							MPI_Reduce(MPI_IN_PLACE,&(tempVolume[z][y][0]),
									(maxVolumeIndexX+1)*2, MPI_FLOAT,
															MPI_SUM, 0, new_comm);
							MPI_Reduce(MPI_IN_PLACE,&(tempWeights[z][y][0]),
									(maxVolumeIndexX+1), MPI_FLOAT,
															MPI_SUM, 0, new_comm);
						} else {
							MPI_Reduce(&(tempVolume[z][y][0]),&(tempVolume[z][y][0]),
									(maxVolumeIndexX+1)*2, MPI_FLOAT,
															MPI_SUM, 0, new_comm);
							MPI_Reduce(&(tempWeights[z][y][0]),&(tempWeights[z][y][0]),
									(maxVolumeIndexX+1), MPI_FLOAT,
															MPI_SUM, 0, new_comm);
						}
					}
				}

				if ( node->rank == 1 )
				{
					gettimeofday(&end_time,NULL);
					finishComputations(fn_out);
					break;
				}
				else
				{
					// release data and finish the thread
					releaseTempSpaces();
					break;
				}
			}
			else if (status.MPI_TAG == TAG_WORKFORWORKER)
			{
				//get the job number
				MPI_Recv(&jobNumber, 1, MPI_INT, 0, TAG_WORKFORWORKER, MPI_COMM_WORLD, &status);

				size_t min_i, max_i;
				min_i = jobNumber*mpi_job_size;
				max_i = min_i + mpi_job_size - 1;

				if ( max_i >= SF.size())
					max_i  = SF.size()-1;
				// do the job
				processImages( min_i, max_i);
			} else {
				std::cerr << "3) Received unknown TAG I quit" << std::endl;
				exit(0);
			}
		}
	}

	// Kill threads used on workers
	if ( node->active && !node->isMaster() )
	{
		cleanLoadingThread();
	}
	delete[] ranks;
}
