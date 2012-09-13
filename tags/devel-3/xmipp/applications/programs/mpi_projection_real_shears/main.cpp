/***************************************************************************
 *
 * Authors:     Slavica JONIC (slavica.jonic@impmc.jussieu.fr, slavica.jonic@a3.epfl.ch)
 *			 Jean-Noel PIOCHE (jnp95@hotmail.com) 
 *		
 * Biomedical Imaging Group, EPFL (Lausanne, Suisse).
 * Structures des Assemblages Macromoleculaires, IMPMC UMR 7590 (Paris, France).
 * IUT de Reims-Chlons-Charleville (Reims, France).
 *
 * Last modifications by JNP the 27/05/2009 15:53:04 
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

#include <mpi.h>
#include <reconstruction/projection_real_shears.h>

#define WORKTAG	1
#define DIETAG		2

class MPI_projection_real_shears : public Projection_real_shears
{
	/*---------------- Fields ----------------*/	
	public :
		int nTasks;
		double work_toSlave[7];

	/*---------------- Functions ----------------*/
	public :
		int master()
		{		
			MPI_Status status;

			if(start_to_process() == ERROR)
			{
				del_VolumeStruct(Data);
				return (ERROR);
			}

			long sizeVolume = Data.nx_Volume * Data.ny_Volume * Data.nz_Volume;
			long sizeResult = (Data.nx_Volume * Data.ny_Volume) +1 ; 
			int numFile_fromSlave;

			long n_Volume[3];
			n_Volume[0] = Data.nx_Volume;
			n_Volume[1] = Data.ny_Volume;
			n_Volume[2] = Data.nz_Volume;
			
			//Data.Output = (double*) malloc((size_t) Data.nx_Volume * Data.ny_Volume * sizeof(double)); MODIF

			//Sends all fixed parameters of Data (VolumeStruct)
			for(int rank = 1; rank < nTasks; ++rank)
			{
				//Dimensions of Volume (it will be used to fill n[x-y-z]_Volume of Data)
				MPI_Send(n_Volume,            // message buffer 
					    3,                   // buffer size 
					    MPI_LONG,            // data item is an integer 
					    rank,                // destination process rank 
					    WORKTAG,             // user chosen message tag 
					    MPI_COMM_WORLD);     // default communicator 

				//Data.Volume 
				MPI_Send(Data.Volume,        // message buffer 
					    sizeVolume,         // buffer size 
					    MPI_DOUBLE,         // data item is an integer 
					    rank,               // destination process rank 
					    WORKTAG,            // user chosen message tag 
					    MPI_COMM_WORLD);    // default communicator  
			}

			// Seed the slaves; send one unit of work to each slave. 
			for (int rank = 1; rank < nTasks; ++rank) 
			{
				if(next_work() == ERROR)
				{
					del_VolumeStruct(Data);
					return (ERROR);
				}			

				// Send it to each rank 
				MPI_Send(work_toSlave,     // message buffer 
					    7,                // buffer size
					    MPI_DOUBLE,       // data item is an integer 
					    rank,             // destination process rank 
					    WORKTAG,          // user chosen message tag 
					    MPI_COMM_WORLD);  // default communicator 
 
				num_file++;
				DF.next_data_line(); 
			}

			// Loop over getting new work requests until there is no more work to be done  
			while (!DF.eof()) 
			{
				// Receive results from a slave 
				MPI_Recv(Data.Output,       // message buffer 
					    sizeResult,        // buffer size 
					    MPI_DOUBLE,        // of type double real 
					    MPI_ANY_SOURCE,    // receive from any sender 
					    MPI_ANY_TAG,       // any type of message 
					    MPI_COMM_WORLD,    // default communicator 
					    &status);          // info about the received message 
					
				if(write_projection_file((int)Data.Output[sizeResult-1]) == ERROR)
				{
					del_VolumeStruct(Data);
					return (ERROR);
				}

				// Get the next unit of work to be done
				if(next_work() == ERROR)
				{
					del_VolumeStruct(Data);
					return (ERROR);
				}	

		    		// Send the slave a new work unit 
		    		MPI_Send(work_toSlave,       // message buffer 
				   		7,                 // buffer size 
				   		MPI_DOUBLE,        // data item is an integer 
				   		status.MPI_SOURCE, // to who we just received from 
				   		WORKTAG,           // user chosen message tag 
				   		MPI_COMM_WORLD);   // default communicator 

		    		 
				num_file++;
				DF.next_data_line(); 
			}

		  	// There's no more work to be done, so receive all the outstanding results from the slaves. 
		  	for (int rank = 1; rank < nTasks; ++rank) 
			{
		    		// Receive results from a slave 
				MPI_Recv(Data.Output,       // message buffer 
					    sizeResult,        // buffer size 
					    MPI_DOUBLE,        // of type double real 
					    MPI_ANY_SOURCE,    // receive from any sender 
					    MPI_ANY_TAG,       // any type of message 
					    MPI_COMM_WORLD,    // default communicator 
					    &status);          // info about the received message

				if(write_projection_file((int)Data.Output[sizeResult-1]) == ERROR)
				{
					del_VolumeStruct(Data);
					return (ERROR);
				}
		  	}

			SF = SF.sort_by_filenames();

			if(finish_to_process() == ERROR)
				return (ERROR);

			//free(Data.Output);  MODIF
			EndTasks();

			return (!ERROR);
		}

		int slave(int myRank)
		{		
			MPI_Status status;
			VolumeStruct Data2;
			int sizeProj;
			double work[7];

			//------------------------- Data initialization ------------------------- 
			long n_Volume[3];

			// Receive a message from the master 
		    	MPI_Recv(n_Volume, 3, MPI_LONG, 0, MPI_ANY_TAG,
				   	MPI_COMM_WORLD, &status);

			// Check the tag of the received message. 
		    	if (status.MPI_TAG == DIETAG) 
			{
				del_VolumeStruct(Data2);
				return(!ERROR);
			}

			Data2.nx_Volume = n_Volume[0];
			Data2.ny_Volume = n_Volume[1];
			Data2.nz_Volume = n_Volume[2];
			sizeProj = Data2.nx_Volume * Data2.ny_Volume;

			allocAndInit_VolumeStruct(Data2);

			// Receive a message from the master 
		    	MPI_Recv(Data2.Volume, sizeProj*Data2.nz_Volume, MPI_DOUBLE, 0, MPI_ANY_TAG,
				   	MPI_COMM_WORLD, &status);

			// Check the tag of the received message. 
		    	if (status.MPI_TAG == DIETAG) 
			{
				del_VolumeStruct(Data2);
				return(!ERROR);
			}
			//-------------------------------------------------------------------

		  	while (true) 
			{
		    		// Receive a message from the master 
		    		MPI_Recv(work, 7, MPI_DOUBLE, 0, MPI_ANY_TAG,
				   		MPI_COMM_WORLD, &status);

		    		// Check the tag of the received message. 
		    		if (status.MPI_TAG == DIETAG) 
				{
					del_VolumeStruct(Data2);
					return(!ERROR);
				}

				//result is a pack constitued by Data.Output and the num_file
				double result[sizeProj +1];

				Data2.InitPsiThetaPhi[0] = work[0];
				Data2.InitPsiThetaPhi[1] = work[1];
				Data2.InitPsiThetaPhi[2] = work[2];
				Data2.InitDelta123[0] = work[3];
				Data2.InitDelta123[1] = work[4];
				Data2.InitDelta123[2] = work[5];
				//Add the num_file at the index ((Data2.nx_Volume * Data2.ny_Volume) +1 -1)
				result[sizeProj] = work[6];

		    		// Do the work 
				if(do_oneProjection(Data2) == ERROR)
				{
					del_VolumeStruct(Data2);
					return (ERROR);
				}

				//Transfer from Data2.Output to result
				for(int i=0; i<sizeProj; i++)
					result[i] = Data2.Output[i];

		    		// Send the result back 
		    		MPI_Send(result, sizeProj+1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

				//free(Data2.Output);  MODIF
		  	}
		}

		void EndTasks()
		{
			// Tell all the slaves to exit by sending an empty message with the DIETAG. 
		  	for (int rank = 1; rank < nTasks; ++rank) 
		    		MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
		}

		int next_work()
		{
			if(read_a_DocLine() == ERROR)
				return (ERROR);

			work_toSlave[0] = Data.InitPsiThetaPhi[0];
			work_toSlave[1] = Data.InitPsiThetaPhi[1];
			work_toSlave[2] = Data.InitPsiThetaPhi[2];
			work_toSlave[3] = Data.InitDelta123[0];
			work_toSlave[4] = Data.InitDelta123[1];
			work_toSlave[5] = Data.InitDelta123[2];
			work_toSlave[6] = (double)num_file;

			return (!ERROR);	
		}
};

int main(int argc, char *argv[])
{
	int myRank;
	MPI_projection_real_shears mpi_proj;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    	{
        	fprintf(stderr, "MPI initialization error\n");
        	exit(EXIT_FAILURE);
    	}

	// Find out my identity in the default communicator 
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	if (myRank == 0)
	{
		// Find out how many processes there are in the default communicator 
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_proj.nTasks);

		//This MPI program doesn't allow to work with only one task
		if (mpi_proj.nTasks<2)
    		{
			std::cout<<"\n\tERROR : Impossible to launch this MPI program with only one task !"<<std::endl<<std::endl;
			mpi_proj.EndTasks();
			exit(EXIT_FAILURE);
    		}

		// Check the command line
    		try
    		{
        		mpi_proj.read(argc, argv);
    		}
    		catch (Xmipp_error &XE)
    		{
        		std::cout << XE <<std::endl;
        		mpi_proj.usage();
			mpi_proj.EndTasks();
        		MPI_Finalize();
			return 1;
    		}

        	// Really project
		try
		{
			mpi_proj.master();
		}
		catch(Xmipp_error &XE)
		{
			std::cout << XE <<std::endl;
			mpi_proj.EndTasks();
        		MPI_Finalize();
			return 1;
    		}
	}
	else
	{
		try
		{
			mpi_proj.slave(myRank);
		}
		catch(Xmipp_error &XE)
		{
			std::cout << XE <<std::endl;
			mpi_proj.EndTasks();
        		MPI_Finalize();
			return 1;
    		}
	}

	// Shut down MPI
	MPI_Finalize();
	
	return 0;	
}

