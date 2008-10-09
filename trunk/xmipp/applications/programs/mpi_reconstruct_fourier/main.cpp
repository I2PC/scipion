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

#include <cstring>
#include <cstdlib>
#include <data/funcs.h>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>  

#define TAG_WORKFORWORKER   0
#define TAG_STOP   1
#define TAG_FREEWORKER   3

class Prog_mpi_RecFourier_prm:Prog_RecFourier_prm
{
    public:
    //int rank, size, num_img_tot;

	/** **/ 
	double * FourierVol;
	double * FourierVolWeight;
	Projection        proj;
	Matrix2D<double>  L(4, 4), R(4, 4), A;
	Mask_Params       mask_prm;
	FileName          fn_img;            
	long int sizeout;

        /** Number of Procesors **/
        int nProcs;
        
        /** Dvide the job in this number block with this number of images */
        int mpi_job_size;

        /** Number of jobs **/
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
        if (nProcs < 2)
            error_exit("This program cannot be executed in a single working node");
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
	int               c, nn, imgno;
	double            rot, tilt, psi, newrot, newtilt, newpsi, xoff, yoff, flip, weight;
	
	//vol().resize(dim, dim, dim);
	//vol().initZeros();
	//alloc memory for volume in Fourier space
	//the allocation is a bit wier but can be reused for fftw
	//I do not want to create a fourier object yet
	//because I only need it at the end of the program
	//and requires to alloc the doble amount of memory
	int paddim_vol=ROUND(dim*padding_factor_vol);
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
	int paddim_proj=ROUND(dim*padding_factor_proj);
	if(paddim_proj%2!=0)
            REPORT_ERROR(1, "\n\nprojection dim * padding_factor must be even.\n\n");
	padding_factor_proj=(double)paddim_proj/(double)dim;
	ImageXmipp paddedImg(paddim_proj,paddim_proj);
	paddedImg().initZeros();
	xmippFftw fftPaddedImg(paddedImg(), false);
	fftPaddedImg.Init("ES",FFTW_FORWARD,false);

	// pointer for output matrix
	std::complex<double> * FOURIERPROJ;
	FOURIERPROJ = (std::complex<double> *)fftPaddedImg.fOut;
	int i_counter=0;
	int screen_size=XMIPP_MIN(60,SF.ImgNo());
	int step_counter=(int)(SF.ImgNo()/screen_size);
	if (verb)
	{
            if (step_counter<=0) step_counter=1;
            init_progress_bar(screen_size);
	}

	//Bcast must be seem by all processors
        if (rank != 0) 
        {

        }        
	//leer sel file / dividir por mpi_job_size 
        numberOfJobs=ceil(SF.ImgNo()/mpi_job_size);
           
        //only one node will write in the console
        verbose=false;
        if (rank == 1)
        {
            if(quiet)
                verbose=false;
            else    
                verbose = true;
            #define DEBUG
            #ifdef DEBUG
            std::cerr << "numberOfJobs: " << numberOfJobs << std::endl
                 << "number of projections to be created: " <<  mysampling.no_redundant_sampling_points_angles.size()
                 <<std::endl;
            #endif
            #undef DEBUG
        }
        
    }
    /* Run --------------------------------------------------------------------- */
    void run()
    {   
	std::complex<double> * FOURIERPROJ;
	FOURIERPROJ = (std::complex<double> *)fftPaddedImg.fOut;
	
	std::complex<double> * FOURIERVOL;
	FOURIERVOL = (std::complex<double> *)FourierVol.fOut;
	
        double * FOURIERVOL_real = FourierVol.fOut;
	double * FOURIERVOLWEIGHT_real = FourierVolWeight.fOut;
	
	if (rank == 0)
        {
            int stopTagsSent =0;
            for (int i=0;i<numberOfJobs;)
            {
                //collect data if available
                //be aware that mpi_Probe will block the program untill a message is received
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                // worker is free
                if (status.MPI_TAG == TAG_FREEWORKER)
                   {
                   MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                         MPI_COMM_WORLD, &status);
//#define DEBUG
#ifdef DEBUG
std::cerr << "Mr_f received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE << std::endl;
#endif
                   //send work
                   MPI_Send(&i,
                            1,
                            MPI_INT,
                            status.MPI_SOURCE,
                            TAG_WORKFORWORKER,
                            MPI_COMM_WORLD);
                    i++; //increase job number       
//#define DEBUG
#ifdef DEBUG
std::cerr << "Ms_f sent TAG_WORKFORWORKER to worker " <<  status.MPI_SOURCE << std::endl;
std::cerr << "Sent jobNo " <<  i << std::endl;
#endif
//#undef DEBUG
                    }
                 else
                    {
                    std::cerr << "M_f Received unknown TAG" << std::endl;
                    exit(0);
                    }           
            }
            

        //send TAG_STOP
        while (stopTagsSent < (nProcs-1))
            {
            MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                  MPI_COMM_WORLD, &status);
    #ifdef DEBUG
    std::cerr << "Mr received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE << std::endl;
    std::cerr << "Ms sent TAG_STOP to worker" << status.MPI_SOURCE << std::endl;
    #endif
    #undef DEBUG
            MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
            stopTagsSent++;
            }         
	ALLREDUCE() DOS VOLUMENES volumen y pesos
        
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
//#define DEBUG
#ifdef DEBUG
std::cerr << "W" << rank << " " << "sent TAG_FREEWORKER to master " << std::endl;
#endif
#undef DEBUG
                //get yor next task
                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#ifdef DEBUG
std::cerr << "W" << rank << " " << "probe MPI_ANY_TAG " << std::endl;
#endif
                if (status.MPI_TAG == TAG_STOP)//no more jobs exit
                    {
                   //If I  do not read this tag
                   //master will no further process
                   //a posibility is a non-blocking send
                   MPI_Recv(0, 0, MPI_INT, 0, TAG_STOP,
                         MPI_COMM_WORLD, &status);
#ifdef DEBUG
std::cerr << "Wr" << rank << " " << "TAG_STOP" << std::endl;
#endif
 	           MPI_Allreduce( (void *) FOURIERVOL, (void *) FOURIERVOL)
                   break;
                    }
                if (status.MPI_TAG == TAG_WORKFORWORKER)
                //there is still some work to be done    
                    {
                    //get the jobs number
                    MPI_Recv(&jobNumber, 1, MPI_INT, 0, TAG_WORKFORWORKER, MPI_COMM_WORLD, &status);
#ifdef DEBUG  
std::cerr << "Wr" << rank << " " << "TAG_WORKFORWORKER" << std::endl;
    std::cerr <<    "jobNumber "  << jobNumber << std::endl;  
#endif
#undef DEBUG
                    // Process all images
		    bucle con las imagenes que son (i*mpi_job_size)
		    comprobar es si el numero i*mpi_job_size + j > SF.numIMG
                    ProcessOneImage
                    //get yor next task
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


