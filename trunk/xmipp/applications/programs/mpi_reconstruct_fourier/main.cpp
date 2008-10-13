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
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>  

#define TAG_WORKFORWORKER   0
#define TAG_STOP   1
#define TAG_TRANSFER 2
#define TAG_FREEWORKER   3

// Really big things to transfer are 
// divided into smaller chunks of size ...
#define BUFFSIZE 10000000

class Prog_mpi_RecFourier_prm:Prog_RecFourier_prm
	{
    public:
	//int rank, size, num_img_tot;

	/** Contains the resulting volume in Fourier space **/ 
	double * FourierVol;
	int paddim_proj;
	int paddim_vol;
	/** Contains normalization (weight) values for 'FourierVol' **/
	double * FourierVolWeight;

	Projection        proj;

	/** Symmetrization variables **/
	Matrix2D<double>  L, R, A;
	Mask_Params       mask_prm;
	FileName          fn_img;            

	ImageXmipp paddedImg;

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
		L.resize(4,4);
		R.resize(4,4);
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
		int               c, nn, imgno;
		double            rot, tilt, psi, newrot, newtilt, newpsi, xoff, yoff, flip, weight;

		//vol().resize(dim, dim, dim);
		//vol().initZeros();
		//alloc memory for volume in Fourier space
		//the allocation is a bit wier but can be reused for fftw
		//I do not want to create a fourier object yet
		//because I only need it at the end of the program
		//and requires to alloc the doble amount of memory
		paddim_vol=ROUND(dim*padding_factor_vol);
		if(paddim_vol%2!=0)
			REPORT_ERROR(1, "\n\nvolume dim * padding_factor must be even.\n\n");
		padding_factor_vol=(double)paddim_vol/(double)dim;
		{
			int ndim    = 3;
			int * myfN;
			try
			{
				myfN = new int [ndim];
			}
			catch (std::bad_alloc&)
			{
				std::cout << "Error allocating memory." << std::endl;
				exit(1);
			}    
			myfN[0]=paddim_vol;
			myfN[1]=paddim_vol;
			myfN[2]=paddim_vol;
			// Set padding dimension
			long int fTotalSize=(long int)myfN[0]*(long int)myfN[1]*(long int)myfN[2];
			sizeout = long(double(fTotalSize)*(long int)(myfN[ndim-1]/2+1)/(long int)myfN[ndim-1]);
			try
			{
				FourierVol = new double[2*sizeout];
			}
			catch (std::bad_alloc&)
			{
				std::cout << "Error allocating memory." << std::endl;
				exit(1);
			}
			for (int i=0; i< 2*sizeout; i++) FourierVol[i]=0.;

			//we will need a volume for weight with same size in x and y and half in z
			try
			{                                 //the 2 factor is not needed
				FourierVolWeight = new double[sizeout];
			}
			catch (std::bad_alloc&)
			{
				std::cout << "Error allocating memory." << std::endl;
				exit(1);
			} 
		}
		//
		// init volumes with zeros
		std::complex<double> * FOURIERVOL;
		FOURIERVOL = (std::complex<double> *)FourierVol;
		for (int i=0; i<sizeout;i++)
		{
			FOURIERVOL[i]=(std::complex<double>)0.0;
			FourierVolWeight[i]=0.;
		}

		SF.go_beginning();
		imgno = 0;
		//create Fourier object of padded size
		//I want to reuse it so I will create it by hand
		paddim_proj=ROUND(dim*padding_factor_proj);
		if(paddim_proj%2!=0)
			REPORT_ERROR(1, "\n\nprojection dim * padding_factor must be even.\n\n");
		padding_factor_proj=(double)paddim_proj/(double)dim;
		paddedImg = ImageXmipp(paddim_proj,paddim_proj);
		paddedImg().initZeros();

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
                xmippFftw fftPaddedImg(paddedImg(), false);
		
		std::complex<double> * FOURIERPROJ;
		FOURIERPROJ = (std::complex<double> *)(fftPaddedImg.fOut);

		std::complex<double> * FOURIERVOL;
		FOURIERVOL = (std::complex<double> *)FourierVol;
		
		fftPaddedImg.Init("ES",FFTW_FORWARD,false);

		if (rank == 0)
		{			
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


			int transfersCompleted = 0;

			double * recBuffer = (double *) malloc( sizeof(double) * BUFFSIZE );
			int receivedSize;
			int currentSource;
			double * pointer;

			// Start collecting results
			while(1)
			{
				MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
						 MPI_COMM_WORLD, &status);

				currentSource = status.MPI_SOURCE;

				// Signal the worker to start sending back the results
				MPI_Send(0, 0, MPI_INT, currentSource, TAG_TRANSFER, MPI_COMM_WORLD);
				
				// Just for convenience, use another "sliding" pointer
				pointer = FourierVol;
				while(1)
				{
					MPI_Probe( currentSource, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					
					// worker is free
					if (status.MPI_TAG == TAG_FREEWORKER)
					{	
						MPI_Recv(0, 0, MPI_INT, currentSource, TAG_FREEWORKER,
							 MPI_COMM_WORLD, &status);
						
						break;
					}

					MPI_Recv( recBuffer, 
							 BUFFSIZE, 
							 MPI_DOUBLE, 
							 currentSource, 
							 MPI_ANY_TAG, 
							 MPI_COMM_WORLD, 
							 &status );
					
					MPI_Get_count(&status,MPI_DOUBLE,&receivedSize);
						
					for( int i = 0; i < receivedSize ; i ++ )
					{
						pointer[i] += recBuffer[i];
					}

					pointer += receivedSize; 
				}
				
				pointer = FourierVolWeight;

				while(1)
				{
					MPI_Probe( currentSource, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

					// worker is free
					if (status.MPI_TAG == TAG_FREEWORKER)
					{
						// Clear buffer
						MPI_Recv(0, 0, MPI_INT, currentSource, TAG_FREEWORKER,
							 MPI_COMM_WORLD, &status);
						
						break;
					}
					
					MPI_Recv( recBuffer, 
							 BUFFSIZE, 
							 MPI_DOUBLE, 
							 currentSource, 
							 MPI_ANY_TAG, 
							 MPI_COMM_WORLD, 
							 &status );

					MPI_Get_count(&status,MPI_DOUBLE,&receivedSize);                    

					for( int i = 0; i < receivedSize ; i ++ )
					{
						pointer[i] += recBuffer[i];
					}

					pointer += receivedSize;
				}         
		
				transfersCompleted++;

				if( transfersCompleted == (nProcs-1) )
				{
					// All threads have finished their job, calculating
					// final volume and storing it.
					
					int Xdim,Xsize;
    					int Ydim=1,Ysize=1;
    					int Zdim=1,Zsize=1;
    					Zdim=Zsize=paddim_vol;
    					Ydim=Ysize=paddim_vol;
    					Xsize=paddim_vol;
    					Xdim = (int)  paddim_vol/2  +1;

					double dZdim = (double) Zdim;
    					for ( int z=0, i=0; z<Zdim; z++ ) 
						for ( int y=0; y<Ydim; y++ ) 
							for ( int x=0; x<Xdim; x++, i++ ) 
            						{
               							if( FourierVolWeight[i] > 0)//MINIMUMWEIGHT)
									FOURIERVOL[i] /= FourierVolWeight[i] * dZdim ;
            						}
					
    					// release memory for weights and create xmippfourier object for volume
    					delete [] FourierVolWeight;
    					VolumeXmipp vol;
					{
        					int ndim    = 3;
        					int * myfN;
        					try
        					{	
            						myfN = new int [ndim];
        					}
        					catch (std::bad_alloc&)
        					{
          						std::cout << "Error allocating memory." << std::endl;
          						exit(1);
        					}    
        					myfN[0]=paddim_vol;
        					myfN[1]=paddim_vol;
        					myfN[2]=paddim_vol;
        					// Set padding dimension
        					int fTotalSize=myfN[0]*myfN[1]*myfN[2];
        					bool inplace = false;
        					xmippFftw Volfft(ndim, myfN, inplace,FourierVol);
        					Volfft.Init("ES",FFTW_BACKWARD,false);
        					Volfft.CenterRealImageInFourierSpace(false) ;//Change phases
        					Volfft.Transform();
        					Volfft.CenterFourierTransformInRealSpace(true);//ok MOVE TRANSFORM TO THE CENTER

        					/* create table with 3D blob fourier transform */
        					//precompute blob fourier transform values close to origin
        					fourier_blob_table = new double [BLOB_TABLE_SIZE];
        					/*it should be divided by 2 not by 4
          					but 4 seems to work. I do not know why */
        					double ww = 1/(4*(double)(BLOB_TABLE_SIZE-1));
        
        					//acces to the blob table should take into account 
        					//pad factor but not sampling
        					double w;
        					double w0= blob_Fourier_val(0., blob);
        					for (int i=0; i<BLOB_TABLE_SIZE; i++)
        					{
          						w = ww*(double)i;
          						fourier_blob_table[i] =  blob_Fourier_val(w, blob)/w0;
          						//#define DEBUG
          						#ifdef DEBUG
          						std::cout.setf(std::ios::scientific); 
          						std::cout.width(12); 
          						std::cout /*<< i << " " 
                    						*/<< w*paddim_vol <<  "\t" 
                    						<< "\t"<< fourier_blob_table[i] << " "
                    						<<  blob.radius
                    						<< std::endl;
          						#endif
          						#undef DEBUG
        					}
        					//copy volume to original volume
        					int xdim=dim;
        					int ydim=dim;
        					int zdim=dim;
        					int center_shiftx = (paddim_vol-xdim)/2+xdim%2;// add one if odd
        					int center_shifty = (paddim_vol-ydim)/2+ydim%2;// add one if odd
        					int center_shiftz = (paddim_vol-zdim)/2+zdim%2;// add one if odd
        					
						vol().resize(dim,dim,dim);
        					for(int i=0;i<zdim;i++)
           						for(int j=0;j<ydim;j++)
            					   		for(int k=0;k<xdim;k++)
                						{           
                							vol(i,j,k)=Volfft.fOut[(center_shiftz + k) + 
                                       					(center_shifty + j) * paddim_vol +
                                       					(center_shiftx + i) * paddim_vol * paddim_vol ];
 	                                              		}
        					vol().setXmippOrigin();

        					for (int k = STARTINGZ(vol()); k <= FINISHINGZ(vol()); k++)
        					{
            						for (int i = STARTINGY(vol()); i <= FINISHINGY(vol()); i++)
            						{
                						for (int j = STARTINGX(vol()); j <= FINISHINGX(vol()); j++)
                						{
                    							double r      = sqrt(k*k+i*i+j*j);
                    							if(r>paddim_vol/2)
                        							vol(i,j,k)=0.;
                    							else
                    							{
                        							double factor = fourier_blob_table[(int)(r*BLOB_TABLE_SIZE/(paddim_vol/2.))];
                        							if (factor > 0.001)
                        							{
                            								vol(i,j,k)   /=  factor;
                        							}
                   							}
                						}
            						}
        					}
        					/*        
        					CENTER XMIPP VOLUME   
        					APPLY FILTER
        					*/     
	
						#ifdef NEVERDEFINED
        					vol().resize(paddim_vol,paddim_vol,paddim_vol);
        					for(int i=0;i<paddim_vol;i++)
           						for(int j=0;j<paddim_vol;j++)
               							for(int k=0;k<paddim_vol;k++)
                						{           
                							vol(i,j,k)=Volfft.fOut[(  k) + 
                                       						(  j) * paddim_vol +
                                       						(  i) * paddim_vol * paddim_vol ];
                						}
						#endif
    					}    
    
    					//symmetrize in reaL SPACE IF Needed 
    					if (fn_sym_vol == "")
    					{
        					vol.write(fn_out);
    					}
    					else 
    					{
        					if (verb > 0) 
            						std::cout << std::endl << "Symmetrizing volume (using Bsplines)" << std::endl;
        					VolumeXmipp     V_out;
        					symmetrize_Bspline(SL_vol, vol, V_out, 3, false, true);
        					V_out.write(fn_out);
    					}
    
    					cout << "Execution completed successfully\n" << endl;
					break;
				}
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
					double * pointer = FourierVol;
					int totalSize = 2 * sizeout;	
					int numChunks = ceil( (double)totalSize / (double)BUFFSIZE );
					int packetSize;

					// Send Fourier Transform
					for( int i = 0 ; i < numChunks ; i ++ )
					{	
						if( i == (numChunks-1))
						{
							packetSize = totalSize-i*BUFFSIZE;
						}
						else
						{
							packetSize = BUFFSIZE;
						}
						
						MPI_Send( pointer, 
								 packetSize,
								 MPI_DOUBLE,
								 0,
								 0,
								 MPI_COMM_WORLD
								 );
						pointer += BUFFSIZE;
					}
					
					MPI_Send(0, 0, MPI_INT, 0, TAG_FREEWORKER, MPI_COMM_WORLD);
					
					pointer = FourierVolWeight;
					totalSize = sizeout;
					numChunks = ceil( (double)totalSize / (double)BUFFSIZE);
					// Send Weights

					for( int i = 0 ; i < numChunks ; i ++ )
					{	
						if( i == (numChunks-1))
						{
							packetSize = totalSize-i*BUFFSIZE;
						}
						else
						{
							packetSize = BUFFSIZE;
						}

						MPI_Send( pointer, 
								 packetSize,
								 MPI_DOUBLE,
								 0,
								 0,
								 MPI_COMM_WORLD
								 );

						pointer += BUFFSIZE;
					}
					
					MPI_Send(0, 0, MPI_INT, 0, TAG_FREEWORKER, MPI_COMM_WORLD);

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

						ProcessOneImage(fn_img,
								fftPaddedImg,
								paddim_proj,
								paddim_vol,
								proj,
								FOURIERVOL,
								FourierVolWeight,
								FOURIERPROJ);
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
    std::cerr << "0\n";
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        fprintf(stderr, "MPI initialization error\n");
        exit(EXIT_FAILURE);
    }
    //size of the mpi block, number of images
    //mpi_job_size=!checkParameter(argc,argv,"-mpi_job_size","-1");

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


