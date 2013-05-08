/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao (jrbcast@ace.ual.es)
 *       Roberto Marabini (roberto@cnb.csic.es)
 *       Vahid Abrishami (vabrishamoi@cnb.csic.es)
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
#include "mpi_reconstruct_fourier.h"

/** Empty constructor */
//ProgMPIRecFourier::ProgMPIRecFourier()
//{
//    node = NULL;
//}

/*  constructor ------------------------------------------------------- */
ProgMPIRecFourier::ProgMPIRecFourier(int argc, char *argv[])
{
    this->read(argc, argv);
}

/* constructor providing an MpiNode
 * this is usefull for using this programs from others
 */
ProgMPIRecFourier::ProgMPIRecFourier(MpiNode * node)
{
    this->setNode(node);
}

/* Special way of reading to sync all nodes */
void ProgMPIRecFourier::read(int argc, char** argv)
{
    XmippMpiProgram::read(argc, argv);
    ProgRecFourier::read(argc, (const char **)argv);
}

/* Usage ------------------------------------------------------------------- */
void ProgMPIRecFourier::defineParams()
{
    ProgRecFourier::defineParams();
    addParamsLine("  [--mpi_job_size <size=10>]    : Number of images sent to a cpu in a single job ");
    addParamsLine("                                : 10 may be a good value");
    addParamsLine("                                : if  -1 the computer will put the maximum");
    addParamsLine("                                : posible value that may not be the best option");
}

/* Read parameters --------------------------------------------------------- */
void ProgMPIRecFourier::readParams()
{
    ProgRecFourier::readParams();
    mpi_job_size=getIntParam("--mpi_job_size");
}

/* Pre Run PreRun for all nodes but not for all works */
void ProgMPIRecFourier::preRun()
{
    if (nProcs < 2)
        REPORT_ERROR(ERR_ARG_INCORRECT,"This program cannot be executed in a single working node");

    if (node->isMaster())
    {
        show();
        SF.read(fn_sel);

        //Send verbose level to node 1
        MPI_Send(&verbose, 1, MPI_INT, 1, TAG_SETVERBOSE, MPI_COMM_WORLD);
    }
    else
    {
        produceSideinfo();
        SF.firstObject();
    }

    //leer sel file / dividir por mpi_job_size
    numberOfJobs=(size_t)ceil((double)SF.size()/mpi_job_size);

    //only one node will write in the console
    if (node->rank == 1 )
    {
        // Get verbose status
        MPI_Recv(&verbose, 1, MPI_INT, 0, TAG_SETVERBOSE, MPI_COMM_WORLD, &status);

        //use threads for volume inverse fourier transform, plan is created in setReal()
        //only rank=1 makes inverse Fourier trnasform
        transformerVol.setThreadsNumber(numThreads);

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
void ProgMPIRecFourier::run()
{
    preRun();
    struct timeval start_time, end_time;
    MPI_Group  orig_group, new_group;
    MPI_Comm   new_comm;
    long int total_usecs;
    double total_time_processing, total_time_weightening, total_time_communicating, total_time;
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

            std::cerr << std::flush;
        }

        nProcs = numberOfJobs;

        // Unused nodes are removed from the MPI communicator
        node->active = (node->rank <= numberOfJobs);
        ////////////////////ROB node->updateComm();
    }

    // Generate a new group to do all reduce without the master
    ranks = new int [nProcs];
    for (int i=0;i<nProcs;i++)
        ranks[i]=i+1;
    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    MPI_Group_incl(orig_group, nProcs, ranks, &new_group);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);

    for (int iter=0;iter<=NiterWeight;iter++)
    {
        if (node->isMaster())
        {
            gettimeofday(&start_time,NULL);

            std::cerr<<std::endl;
            if (iter != NiterWeight)
                std::cerr<<"Computing weights "<<iter+1<<"/"<<NiterWeight<<std::endl;
            else
                std::cerr<<"Computing volume"<<std::endl;

            if ( verbose )
                init_progress_bar(numberOfJobs);

            size_t FSC=numberOfJobs/2;
            for (size_t i=0;i<numberOfJobs;i++)
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
                    REPORT_ERROR(ERR_ARG_INCORRECT,"Unexpected TAG, please contact developers");

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

                if (iter == NiterWeight)
                {
                    if( i == FSC && fn_fsc != "")
                    {
                        // sending every worker COLLECT_FOR_FSC
                        for ( size_t worker = 1 ; worker <= nProcs ; worker ++ )
                        {
                            MPI_Recv(0,
                                     0,
                                     MPI_INT,
                                     MPI_ANY_SOURCE,
                                     TAG_FREEWORKER,
                                     MPI_COMM_WORLD,
                                     &status);

                            MPI_Send( 0,
                                      0,
                                      MPI_INT,
                                      status.MPI_SOURCE,
                                      TAG_COLLECT_FOR_FSC,
                                      MPI_COMM_WORLD);
                        }
                    }
                }

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

            if (iter != NiterWeight)
            {
                gettimeofday(&end_time,NULL);

                total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                total_time_weightening += ((double)total_usecs/(double)1000000);
            }
            else
            {
                gettimeofday(&end_time,NULL);

                total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                total_time_processing += ((double)total_usecs/(double)1000000);
            }

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

            if (iter == NiterWeight)
                if (verbose > 0)
                {
                    std::cout << "\n\nProcessing time: " << total_time_processing << " secs." << std::endl;
                    std::cout << "Transfers time: " << total_time_communicating << " secs." << std::endl;
                    std::cout << "Weighting time: " << total_time_weightening << " secs." << std::endl;
                    std::cout << "Execution completed successfully"<< std::endl;
                }
        }
        else if( node->active )
        {
            // Select only relevant part of selfile for this rank
            // job number
            // job size
            // aux variable
            double * fourierVolume = (double *)VoutFourier.data;
            double * fourierWeights = FourierWeights.data;

            sizeout = MULTIDIM_SIZE(FourierWeights);

            //First
            barrier_init( &barrier, numThreads+1);
            pthread_mutex_init( &workLoadMutex, NULL );
            statusArray = NULL;
            th_ids = (pthread_t *)malloc(numThreads * sizeof(pthread_t));
            th_args = (ImageThreadParams *)malloc(numThreads * sizeof(ImageThreadParams));

            for ( int nt = 0 ; nt < numThreads ; nt++ )
            {
                th_args[nt].parent=this;
                th_args[nt].myThreadID = nt;
                th_args[nt].selFile = new MetaData(SF);
                pthread_create((th_ids+nt),NULL,processImageThread,(void*)(th_args+nt));
            }

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

                if (status.MPI_TAG == TAG_COLLECT_FOR_FSC)
                {
                    //If I  do not read this tag
                    //master will no further process
                    //a posibility is a non-blocking send
                    MPI_Recv(0, 0, MPI_INT, 0, TAG_COLLECT_FOR_FSC, MPI_COMM_WORLD, &status);

                    if( node->rank == 1 )
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
                            for ( size_t i = 2 ; i <= nProcs ; i++)
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
                            }
                        }
                        free( recBuffer );

                        Image<double> auxVolume1;
                        auxVolume1().alias( FourierWeights );
                        auxVolume1.write((std::string)fn_fsc + "_1_Weights.vol");

                        Image< std::complex<double> > auxFourierVolume1;
                        auxFourierVolume1().alias( VoutFourier );
                        auxFourierVolume1.write((std::string) fn_fsc + "_1_Fourier.vol");


                        // Normalize global volume and store data
                        finishComputations(FileName((std::string) fn_fsc + "_split_1.vol"));

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
                    MPI_Recv(0, 0, MPI_INT, 0, TAG_TRANSFER, MPI_COMM_WORLD, &status);
#ifdef DEBUG

                    std::cerr << "Wr" << node->rank << " " << "TAG_STOP" << std::endl;
#endif

                    if (iter != NiterWeight)
                    {
                        int err;
                        err = MPI_Allreduce(MPI_IN_PLACE, fourierWeights,
                                            sizeout, MPI_DOUBLE,
                                            MPI_SUM, new_comm);
                        forceWeightSymmetry(FourierWeights);
                        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(VoutFourier)
                        {
                            double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, k,i,j));
                            ptrOut[0] /= A3D_ELEM(FourierWeights,k,i,j);
                        }
                        if (iter == NiterWeight-1)
                        {
                            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(VoutFourier)
                            {
                                double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, k,i,j));
                                A3D_ELEM(FourierWeights,k,i,j) = ptrOut[0];
                                ptrOut[0] = 0;
                            }
                        }
                        break;
                    }

                    else if ( node->rank == 1 )
                    {
                        // Reserve memory for the receive buffer
                        double * recBuffer = (double *) malloc (sizeof(double)*BUFFSIZE);
                        int receivedSize;
                        double * pointer;
                        pointer = fourierVolume;
                        int currentSource;

                        gettimeofday(&start_time,NULL);

                        if ( nProcs > 1 )
                        {
                            // Receive from other workers

                            for (size_t i = 0 ; i <= (nProcs-2) ; i++)
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
                            }
                        }

                        free( recBuffer );
                        gettimeofday(&end_time,NULL);

                        if( fn_fsc != "")
                        {

                            Image<double> auxVolume2;
                            auxVolume2().alias( FourierWeights );
                            auxVolume2.write((std::string)fn_fsc + "_2_Weights.vol");

                            Image< std::complex<double> > auxFourierVolume2;
                            auxFourierVolume2().alias( VoutFourier );
                            auxFourierVolume2.write((std::string) fn_fsc + "_2_Fourier.vol");


                            // Normalize global volume and store data
                            finishComputations(FileName((std::string) fn_fsc + "_split_2.vol"));

                            Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                            transformerVol.setReal(Vout());
                            Vout().clear();
                            transformerVol.getFourierAlias(VoutFourier);
                            FourierWeights.initZeros(VoutFourier);
                            VoutFourier.initZeros();

                            //int x,y,z;

                            //FourierWeights.getDimension(y,x,z);
                            gettimeofday(&start_time,NULL);

                            auxVolume2.sumWithFile((std::string) fn_fsc + "_1_Weights.vol");

                            gettimeofday(&end_time,NULL);
                            total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                            total_time=(double)total_usecs/(double)1000000;
                            if (verbose > 0)
                                std::cout << "SumFile1: " << total_time << " secs." << std::endl;

                            auxVolume2.sumWithFile((std::string) fn_fsc + "_2_Weights.vol");

                            //VoutFourier.getDimension(y,x,z);
                            auxFourierVolume2.sumWithFile((std::string) fn_fsc + "_1_Fourier.vol");
                            auxFourierVolume2.sumWithFile((std::string) fn_fsc + "_2_Fourier.vol");

                            //remove temporary files
                            remove(((std::string) fn_fsc + "_1_Weights.vol").c_str());
                            remove(((std::string) fn_fsc + "_2_Weights.vol").c_str());
                            remove(((std::string) fn_fsc + "_1_Fourier.vol").c_str());
                            remove(((std::string) fn_fsc + "_2_Fourier.vol").c_str());
                            gettimeofday(&end_time,NULL);
                            total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                            total_time=(double)total_usecs/(double)1000000;
                            if (verbose > 0)
                                std::cout << "SumFile: " << total_time << " secs." << std::endl;

                            /*Save SUM
                                                        //this is an image but not an xmipp image
                                                        auxFourierVolume.write((std::string)fn_fsc + "_all_Fourier.vol",
                                                                false,VDOUBLE);
                                                        auxVolume.write((std::string)fn_fsc + "_all_Weights.vol",
                                                                false,VDOUBLE);
                            */
                        }

                        // Normalize global volume and store data
                        finishComputations(fn_out);
                        break;
                    }
                    else
                    {
                        MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, MPI_COMM_WORLD );

                        sendDataInChunks( fourierVolume, 1, 2 * sizeout, BUFFSIZE, MPI_COMM_WORLD);

                        MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, MPI_COMM_WORLD );

                        break;
                    }
                }
                else if (status.MPI_TAG == TAG_WORKFORWORKER)
                {
                    //get the job number
                    MPI_Recv(&jobNumber, 1, MPI_INT, 0, TAG_WORKFORWORKER, MPI_COMM_WORLD, &status);
                    //LABEL
                    //(if jobNumber == -1) break;
                    threadOpCode=PROCESS_IMAGE;

                    size_t min_i, max_i;

                    min_i = jobNumber*mpi_job_size;
                    max_i = min_i + mpi_job_size - 1;

                    if ( max_i >= SF.size())
                        max_i  = SF.size()-1;
                    if (iter == NiterWeight)
                        processImages( min_i, max_i, false, false);
                    else
                        processImages( min_i, max_i, false, true);
                }
                else
                {
                    std::cerr << "3) Received unknown TAG I quit" << std::endl;
                    exit(0);
                }
            }
        }

        // Kill threads used on workers
        if ( node->active && !node->isMaster() )
        {
            threadOpCode = EXIT_THREAD;
            barrier_wait( &barrier );

            for ( int nt=0; nt<numThreads; nt++)
            {
                pthread_join(*(th_ids+nt),NULL);
            }
            barrier_destroy( &barrier );
        }
    }
}

int  ProgMPIRecFourier::sendDataInChunks( double * pointer, int dest, int totalSize, int buffSize, MPI_Comm comm )
{
    double * localPointer = pointer;

    int numChunks =(int)ceil((double)totalSize/(double)buffSize);
    int packetSize;
    int err=0;

    for ( int i = 0 ; i < numChunks ; i ++ )
    {
        if ( i == (numChunks-1))
            packetSize = totalSize-i*buffSize;
        else
            packetSize = buffSize;

        if ( (err = MPI_Send( localPointer, packetSize,
                              MPI_DOUBLE, dest, 0, MPI_COMM_WORLD ))
             != MPI_SUCCESS )
        {
            break;
        }

        localPointer += packetSize;
    }

    return err;
}





