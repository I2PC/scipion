/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao (jrbcast@ace.ual.es)
 *				Roberto Marabini (roberto@cnb.uam.es)
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
#include <iomanip>  

#define TAG_WORKFORWORKER   0
#define TAG_STOP   1
#define TAG_TRANSFER 2
#define TAG_FREEWORKER   3
#define TAG_COLLECT_FOR_FSC 4

#define BUFFSIZE 10000000
/*
Discuss existing mpi.h files differences ¶

    Types               HP                     IBM      Microsoft       MPICH2	OpenMPI 	                ABI
    MPI_Datatype 	struct hpmp_dtype_s*	int 	int 	        int 	struct ompi_datatype_t* 	TBD
    MPI_Op 	        struct hpmp_op_s* 	int 	int 	        int 	struct ompi_op_t* 	        TBD
    MPI_Comm 	        struct hpmp_comm_s* 	int 	int 	        int 	struct ompi_communicator_t* 	TBD
    MPI_Errhandler	struct hpmp_err_s* 	int 	int 	        int 	struct ompi_errhandler_t* 	TBD
source: https://svn.mpi-forum.org/trac/mpi-forum-web/wiki/Compare_mpi_h
*/    
#ifdef HP_MPI
 #define my_MPI_COMM  hpmp_comm_s* 
#elif OPENMPI
 #define my_MPI_COMM  ompi_mpi_comm_t*   
#else
 #define my_MPI_COMM  int
#endif

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
        if ( rank > 0)
        {
            produce_Side_info();

            SF.go_beginning();
        }
        else
        {
            show();

            SF.read(fn_sel);
        }

        //leer sel file / dividir por mpi_job_size 
        numberOfJobs=ceil((double)SF.ImgNo()/mpi_job_size);

        //only one node will write in the console
        if (rank == 1 )
        {
            //#define DEBUG
#ifdef DEBUG
            std::cerr << "SF.ImgNo() mpi_job_size " 
            << SF.ImgNo() << " "
            << mpi_job_size
            << std::endl;
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

        if (rank == 0)
        {
            gettimeofday(&start_time,NULL);

            if ( verb )
                init_progress_bar(numberOfJobs);

            int FSC=numberOfJobs/2;

            for (int i=0;i<numberOfJobs;)
            {

//#define DEBUG
#ifdef DEBUG
                std::cerr << "master-recv  i=" << i << std::endl;
                std::cerr << "numberOfJobs: " << numberOfJobs << std::endl <<std::endl;
#endif
#undef DEBUG
                MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                         MPI_COMM_WORLD, &status);

                if ( status.MPI_TAG != TAG_FREEWORKER )
                {
                    std::cout << "Unexpected TAG, please contact developers " << std::endl;
                    exit(-1);
                }

//#define DEBUG
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

                
                i++; //increase job number 
                if( i == FSC && fn_fsc != "" )
                {                
                    // sending every worker COLLECT_FOR_FSC
                    for ( int worker = 1 ; worker < nProcs ; worker ++ )
                    {
                        MPI_Send(0,
                                  0,
                                  MPI_INT,
                                  worker,
                                  TAG_COLLECT_FOR_FSC,
                                  MPI_COMM_WORLD); 
                    }   
                }

                if ( verb)
                    progress_bar(i);
            }

            gettimeofday(&end_time,NULL);

            total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
            total_time=(double)total_usecs/(double)1000000;     

            std::cout << std::flush << std::endl;
            std::cout << "Processing time: " << total_time << " secs." << std::endl;

            int currentSource;

            // Start collecting results
            for ( int i = 1 ; i < nProcs ; i ++ )
            {
                MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                         MPI_COMM_WORLD, &status);

                currentSource = status.MPI_SOURCE;

                // Signal the worker to start sending back the results to
                // Worker whose rank==1
                MPI_Send(0, 0, MPI_INT, currentSource, TAG_TRANSFER, MPI_COMM_WORLD);
            }   
        }
        else
        {
            // Select only relevant part of selfile for this rank
            // job number
            // job size
            // aux variable
            double * fourierVolume = (double *)VoutFourier.data;
            double * fourierWeights = FourierWeights.data;

            sizeout = MULTIDIM_SIZE(FourierWeights);

            barrier_init( &barrier, numThreads+1);
            pthread_mutex_init( &workLoadMutex, NULL );
            statusArray = NULL;
            th_ids = (pthread_t *)malloc(numThreads * sizeof(pthread_t));
            th_args = (ImageThreadParams *)malloc(numThreads * sizeof(ImageThreadParams));

            for ( int nt = 0 ; nt < numThreads ; nt++ )
            {
                th_args[nt].parent=this;
                th_args[nt].myThreadID = nt;
                th_args[nt].docFile = new DocFile( DF );
                pthread_create((th_ids+nt),NULL,processImageThread,(void*)(th_args+nt));
            }

            while (1)
            {
                int jobNumber;
                
//#define DEBUG
#ifdef DEBUG
                std::cerr << "slave-send TAG_FREEWORKER rank=" << rank << std::endl;
#endif
#undef DEBUG
                //I am free
                MPI_Send(0, 0, MPI_INT, 0, TAG_FREEWORKER, MPI_COMM_WORLD);

                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                if (status.MPI_TAG == TAG_COLLECT_FOR_FSC)
                {
                    //If I  do not read this tag
                    //master will no further process
                    //a posibility is a non-blocking send
                    MPI_Recv(0, 0, MPI_INT, 0, TAG_COLLECT_FOR_FSC, MPI_COMM_WORLD, &status);
                    
                    if( rank == 1 )
                    {
                        // Reserve memory for the receive buffer
                        double * recBuffer = (double *) malloc (sizeof(double)*BUFFSIZE);
                        int receivedSize;
                        double * pointer;
                        pointer = fourierVolume;
                        int currentSource;

                        if ( nProcs > 2 )
                        {
                            // Receive from other workers
                            for ( int i = 0 ; i < (nProcs-2) ; i++)
                            {
                                MPI_Recv(0,0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                                         MPI_COMM_WORLD, &status);

                                currentSource = status.MPI_SOURCE;

                                pointer = fourierVolume;

                                while (1)
                                {
                                    MPI_Probe( currentSource, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

                                    if ( status.MPI_TAG == TAG_FREEWORKER )
                                    {
                                        MPI_Recv(0,0, MPI_INT, currentSource, TAG_FREEWORKER, MPI_COMM_WORLD, &status );

                                        break;
                                    }

                                    MPI_Recv( recBuffer,
                                              BUFFSIZE,
                                              MPI_DOUBLE,
                                              currentSource,
                                              MPI_ANY_TAG,
                                              MPI_COMM_WORLD,
                                              &status );

                                    MPI_Get_count( &status, MPI_DOUBLE, &receivedSize );

                                    for ( int i = 0 ; i < receivedSize ; i ++ )
                                    {
                                        pointer[i] += recBuffer[i];
                                    }

                                    pointer += receivedSize;    
                                }

                                pointer = fourierWeights;

                                while (1)
                                {
                                    MPI_Probe( currentSource, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

                                    if ( status.MPI_TAG == TAG_FREEWORKER )
                                    {
                                        MPI_Recv( 0,0,MPI_INT,currentSource,TAG_FREEWORKER, MPI_COMM_WORLD,&status );

                                        break;
                                    }

                                    MPI_Recv( recBuffer,
                                              BUFFSIZE,
                                              MPI_DOUBLE,
                                              currentSource,
                                              MPI_ANY_TAG,
                                              MPI_COMM_WORLD,
                                              &status );

                                    MPI_Get_count( &status, MPI_DOUBLE, &receivedSize );

                                    for ( int i = 0 ; i < receivedSize ; i ++ )
                                    {
                                        pointer[i] += recBuffer[i];
                                    }

                                    pointer += receivedSize;
                                }
                            }
                        }
                        free( recBuffer );

                        VolumeXmipp save;
                        save().alias( FourierWeights );
                        save.write((std::string)fn_fsc + "_1_Weights.vol");
                        
                        FourierVolume save2;
                        save2().alias( VoutFourier );
                        save2.write((std::string) fn_fsc + "_1_Fourier.vol");
                        
                        // Normalize global volume and store data
                        finishComputations(FileName((std::string) fn_fsc + "_1_recons.vol"));
                        
                        Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                        transformerVol.setReal(Vout());
                        Vout().clear();
                        transformerVol.getFourierAlias(VoutFourier);
                        FourierWeights.initZeros(VoutFourier);
                        VoutFourier.initZeros();
                    }
                    else
                    {
                        MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, MPI_COMM_WORLD );
          
                        sendDataInChunks( fourierVolume, 1, 2*sizeout, BUFFSIZE, MPI_COMM_WORLD );

                        MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, MPI_COMM_WORLD );

                        sendDataInChunks( fourierWeights, 1, sizeout, BUFFSIZE, MPI_COMM_WORLD);

                        MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER,MPI_COMM_WORLD);
                            
                        Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                        transformerVol.setReal(Vout());
                        Vout().clear();
                        transformerVol.getFourierAlias(VoutFourier);
                        FourierWeights.initZeros(VoutFourier);
                        VoutFourier.initZeros();
                    }
                }
                else if (status.MPI_TAG == TAG_TRANSFER)
                {
                    //If I  do not read this tag
                    //master will no further process
                    //a posibility is a non-blocking send
                    MPI_Recv(0, 0, MPI_INT, 0, TAG_TRANSFER, MPI_COMM_WORLD, &status);
#ifdef DEBUG
                    std::cerr << "Wr" << rank << " " << "TAG_STOP" << std::endl;
#endif
                    if ( rank == 1 )
                    {
                        // Reserve memory for the receive buffer
                        double * recBuffer = (double *) malloc (sizeof(double)*BUFFSIZE);
                        int receivedSize;
                        double * pointer;
                        pointer = fourierVolume;
                        int currentSource;

                        gettimeofday(&start_time,NULL);

                        if ( nProcs > 2 )
                        {
                            // Receive from other workers
                            for ( int i = 0 ; i < (nProcs-2) ; i++)
                            {
                                MPI_Recv(0,0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                                         MPI_COMM_WORLD, &status);

                                currentSource = status.MPI_SOURCE;

                                pointer = fourierVolume;

                                while (1)
                                {
                                    MPI_Probe( currentSource, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

                                    if ( status.MPI_TAG == TAG_FREEWORKER )
                                    {
                                        MPI_Recv(0,0, MPI_INT, currentSource, TAG_FREEWORKER, MPI_COMM_WORLD, &status );

                                        break;
                                    }

                                    MPI_Recv( recBuffer,
                                              BUFFSIZE,
                                              MPI_DOUBLE,
                                              currentSource,
                                              MPI_ANY_TAG,
                                              MPI_COMM_WORLD,
                                              &status );

                                    MPI_Get_count( &status, MPI_DOUBLE, &receivedSize );

                                    for ( int i = 0 ; i < receivedSize ; i ++ )
                                    {
                                        pointer[i] += recBuffer[i];
                                    }

                                    pointer += receivedSize;    
                                }

                                pointer = fourierWeights;

                                while (1)
                                {
                                    MPI_Probe( currentSource, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

                                    if ( status.MPI_TAG == TAG_FREEWORKER )
                                    {
                                        MPI_Recv( 0,0,MPI_INT,currentSource,TAG_FREEWORKER, MPI_COMM_WORLD,&status );

                                        break;
                                    }

                                    MPI_Recv( recBuffer,
                                              BUFFSIZE,
                                              MPI_DOUBLE,
                                              currentSource,
                                              MPI_ANY_TAG,
                                              MPI_COMM_WORLD,
                                              &status );

                                    MPI_Get_count( &status, MPI_DOUBLE, &receivedSize );

                                    for ( int i = 0 ; i < receivedSize ; i ++ )
                                    {
                                        pointer[i] += recBuffer[i];
                                    }

                                    pointer += receivedSize;
                                }
                            }
                        }

                        free( recBuffer );
                        gettimeofday(&end_time,NULL);

                        total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                        total_time=(double)total_usecs/(double)1000000;     

                        std::cout << "Transfers time: " << total_time << " secs." << std::endl;

                        if( fn_fsc != "" )
                        {
                            VolumeXmipp auxVolume;
                            auxVolume().alias( FourierWeights );
                            auxVolume.write((std::string)fn_fsc + "_2_Weights.vol");
                        
                            FourierVolume auxFourierVolume;
                            auxFourierVolume().alias( VoutFourier );
                            auxFourierVolume.write((std::string) fn_fsc + "_2_Fourier.vol");
                        
                            // Normalize global volume and store data
                            finishComputations(FileName((std::string) fn_fsc + "_2_recons.vol"));
                        
                            Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                            transformerVol.setReal(Vout());
                            Vout().clear();
                            transformerVol.getFourierAlias(VoutFourier);
                            FourierWeights.initZeros(VoutFourier);
                            VoutFourier.initZeros();   
                            
                            int x,y,z;
                            
                            auxVolume.sumWithFile(FileName((std::string) fn_fsc + "_1_Weights.vol"));
                            auxVolume.sumWithFile(FileName((std::string) fn_fsc + "_2_Weights.vol"));
                            
                            VoutFourier.getDimension(y,x,z);
                            
                            auxFourierVolume.sumWithFile((std::string) fn_fsc + "_1_Fourier.vol",z,y,x,false, VCOMPLEX);
                            auxFourierVolume.sumWithFile((std::string) fn_fsc + "_2_Fourier.vol",z,y,x,false, VCOMPLEX);
                        }

                        // Normalize global volume and store data
                        gettimeofday(&start_time,NULL);
                        finishComputations(fn_out);
                        gettimeofday(&end_time,NULL);

                        total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                        total_time=(double)total_usecs/(double)1000000;     

                        std::cout << "Weighting time: " << total_time << " secs." << std::endl;
                        std::cout << "Execution completed successfully\n" << std::endl;
                        break;  
                    }
                    else
                    {
                        MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, MPI_COMM_WORLD );
                        
                        sendDataInChunks( fourierVolume, 1, 2 * sizeout, BUFFSIZE, MPI_COMM_WORLD);

                        MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, MPI_COMM_WORLD );

                        sendDataInChunks( fourierWeights, 1, sizeout, BUFFSIZE, MPI_COMM_WORLD);

                        MPI_Send(0,0,MPI_INT,1,TAG_FREEWORKER,MPI_COMM_WORLD);

                        break;
                    }
                }
                else if (status.MPI_TAG == TAG_WORKFORWORKER)
                {
                    threadOpCode=PROCESS_IMAGE; 

                    //get the jobs number
                    MPI_Recv(&jobNumber, 1, MPI_INT, 0, TAG_WORKFORWORKER, MPI_COMM_WORLD, &status);
                    int min_i, max_i;

                    min_i = jobNumber*mpi_job_size;
                    max_i = min_i + mpi_job_size - 1;

                    if ( max_i >= SF.ImgNo())
                        max_i  = SF.ImgNo()-1;

                    processImages( min_i, max_i );
                }
                else
                {
                    std::cerr << "3) Received unknown TAG I quit" << std::endl;
                    exit(0);
                }           
            }
        }

        // Synchronize all workers
        MPI_Barrier( MPI_COMM_WORLD );

        // Kill threads used on workers
        if ( rank > 0 )
        {
            threadOpCode=EXIT_THREAD;
            barrier_wait( &barrier );

            for ( int nt=0; nt<numThreads; nt++)
            {
                pthread_join(*(th_ids+nt),NULL);
            }

            barrier_destroy( &barrier );
        }

        MPI_Finalize();
    }

    /* a short function to print a message and exit */
    void error_exit(char * msg)
    {
        fprintf(stderr, "%s", msg);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int  sendDataInChunks( double * pointer, int dest, int totalSize, int buffSize, my_MPI_COMM comm )
    {
        double * localPointer = pointer;

        int numChunks = ceil((double)totalSize/(double)buffSize);
        int packetSize;
        int err=0;

        for ( int i = 0 ; i < numChunks ; i ++ )
        {
            if ( i == (numChunks-1))
            {
                packetSize = totalSize-i*buffSize;
            }
            else
            {
                packetSize = buffSize;
            }

            if ( err = MPI_Send( localPointer,
                                 packetSize,
                                 MPI_DOUBLE,
                                 dest,
                                 0,
                                 comm
                               ) != MPI_SUCCESS )
            {
                break;
            }

            localPointer += packetSize;
        }

        return err;
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


