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
#include <reconstruction/angular_projection_matching_project.h>
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
#define TAG_WAIT   2
#define TAG_FREEWORKER   3

class Prog_mpi_projection_matching_project_prm:Prog_Angular_Projection_Matching_Project_Parameters
{
    public:
    //int rank, size, num_img_tot;


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
    Prog_mpi_projection_matching_project_prm()
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
        Prog_Angular_Projection_Matching_Project_Parameters::read(argc,argv);
        mpi_job_size=textToInteger(getParameter(argc,argv,"-mpi_job_size","-1"));
        
    }

    /* Usage ------------------------------------------------------------------- */
    void usage()
    {
        Prog_Angular_Projection_Matching_Project_Parameters::usage();
        cerr << " [ -mpi_job_size default=-1]    : Number of images sent to a cpu in a single job \n";
        cerr << "                                 if  -1 the computer will fill the value for you";
    }


    /* Show -------------------------------------------------------------------- */
    void show()
    {
        Prog_Angular_Projection_Matching_Project_Parameters::show();
	cerr << " Size of mpi jobs " << mpi_job_size <<endl;
    }

    /* Pre Run --------------------------------------------------------------------- */
    void preRun()
    {
        if (rank == 0) 
        {
            show();
        }
        //all ranks
        mysampling.SetSampling(sampling);
        //mysampling.SetNeighborhoodRadius(0.);//irelevant
        //true -> half_sphere
        mysampling.Compute_sampling_points(false,max_tilt_angle,min_tilt_angle);
        mysampling.remove_redundant_points(symmetry, sym_order);
        remove_points_not_close_to_experimental_points();
        if (rank == 0) 
        {
            mysampling.create_sym_file(symmetry, sym_order);
            mysampling.create_asym_unit_file(output_file_root);
        }
        
        if (rank != 0) 
        {
        inputVol.read(input_volume);
        inputVol().setXmippOrigin();
        Xdim = XSIZE(inputVol());
        Ydim = YSIZE(inputVol());
        }        
        if (mpi_job_size != -1)
        {   
            numberOfJobs = ceil((double)(close_points_angles.size())/mpi_job_size);
        }
        else
        {   
            numberOfJobs=nProcs-1;//one node is the master
            mpi_job_size=ceil((double)close_points_angles.size()/numberOfJobs);
        } 
           
        //only one node will write in the console
        verbose=false;
        if (rank == 1)
        {
            verbose = true;
            #define DEBUG
            #ifdef DEBUG
            cerr << "numberOfJobs " << numberOfJobs << endl
                 << "mpi_job_size " << mpi_job_size << endl
                 << "close_points_angles.size()" <<  close_points_angles.size()
                 <<endl;
            #endif
            #undef DEBUG
        }
        
    }
    /* Run --------------------------------------------------------------------- */
    void run()
    {   
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
cerr << "Mr_f received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE << endl;
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
cerr << "Ms_f sent TAG_WORKFORWORKER to worker " <<  status.MPI_SOURCE << endl;
cerr << "Sent jobNo " <<  i << endl;
#endif
//#undef DEBUG
                    }
                 else
                    {
                    cerr << "M_f Recived unknown TAG" << endl;
                    exit(0);
                    }           
            }
            

        //send TAG_STOP
        while (stopTagsSent < (nProcs-1))
            {
            MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                  MPI_COMM_WORLD, &status);
    #ifdef DEBUG
    cerr << "Mr received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE << endl;
    cerr << "Ms sent TAG_STOP to worker" << status.MPI_SOURCE << endl;
    #endif
    #undef DEBUG
            MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
            stopTagsSent++;
            }         
        //only rank 0 create sel file
        if(rank==0)
        {
            SelFile  mySF;
            FileName fn_temp;
            int myCounter=0;

            for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
               for (int i=0;i<=close_points_angles.size()-1;i++)
               { 
                fn_temp.compose(output_file_root, myCounter++,"xmp");
                mySF.insert(fn_temp);
               }
            fn_temp=output_file_root+".sel";   
            mySF.write(fn_temp);         
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
//#define DEBUG
#ifdef DEBUG
cerr << "W" << rank << " " << "sent TAG_FREEWORKER to master " << endl;
#endif
#undef DEBUG
                //get yor next task
                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#ifdef DEBUG
cerr << "W" << rank << " " << "probe MPI_ANY_TAG " << endl;
#endif
                if (status.MPI_TAG == TAG_STOP)//no more jobs exit
                    {
                   //If I  do not read this tag
                   //master will no further process
                   //a posibility is a non-blocking send
                   MPI_Recv(0, 0, MPI_INT, 0, TAG_STOP,
                         MPI_COMM_WORLD, &status);
#ifdef DEBUG
cerr << "Wr" << rank << " " << "TAG_STOP" << endl;
#endif
                    break;
                    }
                if (status.MPI_TAG == TAG_WORKFORWORKER)
                //there is still some work to be done    
                    {
                    //get the jobs number
                    MPI_Recv(&jobNumber, 1, MPI_INT, 0, TAG_WORKFORWORKER, MPI_COMM_WORLD, &status);
#ifdef DEBUG  
cerr << "Wr" << rank << " " << "TAG_WORKFORWORKER" << endl;
    cerr <<    "jobNumber "  << jobNumber << endl;  
#endif
#undef DEBUG
                    // Process all images
                     project_angle_vector(jobNumber*mpi_job_size,
                     XMIPP_MIN((jobNumber+1)* mpi_job_size -1 , 
                                close_points_angles.size()-1)
                                 );
                    //get yor next task
                    }
                else
                   {
                   cerr << "3) Recived unknown TAG I quit" << endl;
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
    //size of the mpi block, number of images
    //mpi_job_size=!checkParameter(argc,argv,"-mpi_job_size","-1");

    Prog_mpi_projection_matching_project_prm prm;
    try
    {
        prm.read(argc, argv);
    }

    catch (Xmipp_error XE)
    {
        cerr << XE;
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
        cerr << XE;
        exit(1);
    }
    
    exit(0);
}


