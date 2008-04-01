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
#include <reconstruction/new_projmatch.h>
#include <data/header.h>
#include <reconstruction/sampling.h>
#include <reconstruction/symmetries.h>

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
#define TAG_WORKFROMWORKER   4

int     * input_images ;
int       input_images_size;
double  * output_values;
int       output_values_size;
//consider
//class Prog_mpi_projection_matching_prm:public Prog_projection_matching_prm
//to access parent variables
class Prog_mpi_new_projection_matching_prm:Prog_new_projection_matching_prm
{
    public:
    //int rank, size, num_img_tot;

        /** Number of Procesors **/
        int nProcs;
        
        /** Dvide the job in this number block with this number of images */
        int mpi_job_size;

        /** classify the experimental data making voronoi regions
            with an hexagonal grid mapped onto a sphere surface */
        int chuck_angular_distance;
               
        /** computing node number. Master=0 */
        int rank;

        /** status after am MPI call */
        MPI_Status status;
                
        /** total number of images */
        int num_img_tot;

        /** symmetry file */
        FileName        fn_sym;
        
        /** sampling object */
        XmippSampling mysampling;

        /** Symmetry. One of the 17 possible symmetries in
            single particle electron microscopy.
             */
        int symmetry;

        /** For infinite groups symmetry order*/
        int sym_order;

    /*  constructor ------------------------------------------------------- */
    Prog_mpi_new_projection_matching_prm()
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
        Prog_new_projection_matching_prm::read(argc,argv);
        mpi_job_size=textToInteger(getParameter(argc,argv,"-mpi_job_size","10"));
        chuck_angular_distance = checkParameter(argc, argv,"-chuck_angular_distance");
        chuck_angular_distance = textToFloat(getParameter(argc, argv,"-chuck_angular_distance"));
        fn_sym = getParameter(argc, argv, "-sym","c1");
    }

    /* Usage ------------------------------------------------------------------- */
    void usage()
    {
        Prog_new_projection_matching_prm::usage();
        std::cerr << " [ -mpi_job_size default=-1]    : Number of images sent to a cpu in a single job \n";
        std::cerr << "                                  10 may be a good value\n";
        std::cerr << "                                 if  -1 the computer will fill the value for you\n";
        std::cerr << " [ -chuck_angular_distance N]    : sample the projection sphere with this \n";
        std::cerr << "                                   sampling rate and create subsets of experimental\n";
        std::cerr << "                                   using the voronoi regions\n";
        std::cerr << "  [-sym cn]   :One of the 17 possible symmetries in\n"
                  << "                                single particle electronmicroscopy\n"
                  << "                                i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv,\n "
                  << "                                dnh, t, td, th, o, oh, i1 (default MDB), i2, i3, i4, ih\n"
                  << "                                i1h (default MDB), i2h, i3h, i4h\n"
                  << "                               : where n may change from 1 to 99\n"
                  ;
     }


    /* Show -------------------------------------------------------------------- */
    void show()
    {
        Prog_new_projection_matching_prm::show();
	std::cerr << " Size of mpi jobs " << mpi_job_size <<std::endl
              << " Sampling rate(chuck_angular_distance): " << chuck_angular_distance    << std::endl
              << " Symmetry group:            " << fn_sym << std::endl
              ;
    }

    /* Pre Run --------------------------------------------------------------------- */
    void preRun()
    {
        int max_number_of_images_in_around_a_sampling_point=0;
        if (rank == 0) 
        {
            show();
            //read experimental doc file
            DFexp.read(fn_exp);
            //first set sampling rate
            mysampling.SetSampling(chuck_angular_distance);
            //create sampling points in the whole sphere
            mysampling.Compute_sampling_points(false,91.,-91.);
            //process the symmetry file
            if (!mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
                 REPORT_ERROR(3005, (std::string)"mpi_angular_proj_match::prerun Invalid symmetry: " +  fn_sym);
            mysampling.SL.read_sym_file(fn_sym);
            //store symmetry matrices, this is faster than computing them 
            //each time we need them
            mysampling.fill_L_R_repository();
            //precompute product between symmetry matrices 
            //and experimental data
            mysampling.fill_exp_data_projection_direction_by_L_R(fn_exp);
            //remove redundant sampling points: symmetry
            mysampling.remove_redundant_points(symmetry, sym_order);
           //remove sampling points too far away from experimental data
            mysampling.remove_points_far_away_from_experimental_data(fn_exp);
            //for each sampling point find the experimental images
            //closer to that point than to any other
	        mysampling.find_closest_experimental_point(fn_exp);
            //print number of points per node
            std::cerr << "voronoi region, number of elements" << std::endl;
            for (int j = 0;
                  j < mysampling.my_exp_img_per_sampling_point.size();
                  j++)   
                    std::cerr << j 
                              << " " 
                              << mysampling.my_exp_img_per_sampling_point[j].size()
                              << std::endl;
            for (int j = 0;
                  j < mysampling.my_exp_img_per_sampling_point.size();
                  j++)   
                    {
                     if (max_number_of_images_in_around_a_sampling_point  
                         < mysampling.my_exp_img_per_sampling_point[j].size())
                         max_number_of_images_in_around_a_sampling_point
                         = mysampling.my_exp_img_per_sampling_point[j].size();       
                    }
            std::cerr << "biggest subset " 
                      << max_number_of_images_in_around_a_sampling_point 
                      << std::endl;
            //alloc memory for buffer          
           if (mpi_job_size == -1)
            {   
                int numberOfJobs=nProcs-1;//one node is the master
                mpi_job_size=ceil((double)DFexp.dataLineNo()/numberOfJobs);
            }    
        } 
        MPI_Bcast(&max_number_of_images_in_around_a_sampling_point, 
                  1, MPI_INT, 0, MPI_COMM_WORLD);
        input_images_size = max_number_of_images_in_around_a_sampling_point+1 +1;
        input_images  = (int *)    malloc(input_images_size*sizeof(int));
        output_values_size=MY_OUPUT_SIZE*max_number_of_images_in_around_a_sampling_point+1;
        output_values = (double *) malloc(output_values_size*sizeof(double));
           
        //only one node will write in the console
        if (rank != 1)
            verb = 0;
        else
            verb = 1;
        //initialize each node, this shoud be out of run
        //because is made once per node but not one per packet
        
        //many sequential programs free object alloc in side_info
        //becareful with that
        if (rank != 0)
        {
            produceSideInfo();
        }
    }

    /* Run --------------------------------------------------------------------- */
    void run()
    {   
        if (rank == 0)
        {
            int N = mysampling.my_exp_img_per_sampling_point.size();
            int killed_jobs=0;
            int index=0;
            int index_counter=0;
            int tip=-1;
            int number_of_processed_images=0;
            int stopTagsSent =0;
            int total_number_of_images=DFexp.dataLineNo();
            Matrix1D<double>                 dataline(8);
            //DocFile DFo;
            FileName                         fn_tmp;
            
            DFexp.go_beginning();
            DFexp.remove_current();
            DFexp.go_beginning();
            DFexp.insert_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Refno (6), Flip (7), maxCC (8)");
            DFexp.set(1,7,0);
            while(1)
            {
                //Wait until any message arrives
                //be aware that mpi_Probe will block the program untill a message is received
                //#define DEBUG
                #ifdef DEBUG
                std::cerr << "Mp1 waiting for any  message " << std::endl;
                #endif
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                #ifdef DEBUG
                std::cerr << "Mp2 received tag from worker " <<  status.MPI_SOURCE << std::endl;
                #endif
                // worker sends work
                if (status.MPI_TAG == TAG_WORKFROMWORKER)
                {
                    MPI_Recv(output_values, 
                             output_values_size, 
                             MPI_DOUBLE, 
                             MPI_ANY_SOURCE, 
                             TAG_WORKFROMWORKER,
                             MPI_COMM_WORLD, 
                             &status);
                    int number= round(output_values[0]/MY_OUPUT_SIZE);        
                    //create doc file
                    for (int i = 0; i < number; i++)
	                {
                        int lineNumber=round(output_values[i*MY_OUPUT_SIZE+1]+1);
	                    DFexp.locate(lineNumber);
	                    DFexp.set(0,output_values[i*MY_OUPUT_SIZE+2]);
	                    DFexp.set(1,output_values[i*MY_OUPUT_SIZE+3]);
	                    DFexp.set(2,output_values[i*MY_OUPUT_SIZE+4]);
	                    DFexp.set(3,output_values[i*MY_OUPUT_SIZE+5]);
	                    DFexp.set(4,output_values[i*MY_OUPUT_SIZE+6]);
	                    DFexp.set(5,output_values[i*MY_OUPUT_SIZE+7]+1);
	                    DFexp.set(6,output_values[i*MY_OUPUT_SIZE+8]);
	                    DFexp.set(7,output_values[i*MY_OUPUT_SIZE+9]);
                    }         
                }
                // worker is free
                else if (status.MPI_TAG == TAG_FREEWORKER)
                {
                    MPI_Recv(&tip, 1, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                         MPI_COMM_WORLD, &status);
                    //#define DEBUG     
                    #ifdef DEBUG
                    std::cerr << "Mr3 received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE 
                              << "with tip= " << tip
                              << "\nnumber_of_processed_images "
                              << number_of_processed_images
                              << "\ntotal_number_of_images "
                              << total_number_of_images
                              << std::endl;
                    #endif
                    #undef DEBUG
                    if(number_of_processed_images>=total_number_of_images)
                        {
                        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
                        stopTagsSent++;
                        #ifdef DEBUG
                        std::cerr << "Ms4 sent stop tag to worker " <<  status.MPI_SOURCE << std::endl
                                 << std::endl;
                        #endif
                        break;
                        }
                    if(tip==-1)
                    {
                        index=index_counter;
                    }
                    else
                        index=tip;
                    while( mysampling.my_exp_img_per_sampling_point[index].size()<=0)
                    {
                        index_counter = index_counter++;
                        index_counter = index_counter%N;
                        index = index_counter;
                    }
                    int number_of_images_to_transfer=XMIPP_MIN(
                                                mysampling.my_exp_img_per_sampling_point[index].size(),
                                                mpi_job_size);

                    input_images[0]=index;
                    input_images[1]=number_of_images_to_transfer;


                    for(int k=2;k<number_of_images_to_transfer+2;k++)
                       {
                       number_of_processed_images++;
                       input_images[k]=mysampling.my_exp_img_per_sampling_point[index].back();
                       mysampling.my_exp_img_per_sampling_point[index].pop_back();
                       }

                    MPI_Send(input_images,
                             number_of_images_to_transfer+2,
                             MPI_INT, 
                             status.MPI_SOURCE,
                             TAG_WORKFORWORKER,
                             MPI_COMM_WORLD);
                    #ifdef DEBUG
                    std::cerr << "Ms_s send work for worker " 
                              <<  status.MPI_SOURCE 
                              << "with index " << index%N
                              << std::endl;
                    #endif
                    }//TAG_FREEWORKER
            }//while       

            while (stopTagsSent < (nProcs-1))
            {
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == TAG_WORKFROMWORKER)
                {
                    MPI_Recv(output_values, 
                             output_values_size, 
                             MPI_DOUBLE, 
                             MPI_ANY_SOURCE, 
                             TAG_WORKFROMWORKER,
                             MPI_COMM_WORLD, 
                             &status);
                    int number= round(output_values[0]/MY_OUPUT_SIZE);        
                    //create doc file
                    for (int i = 0; i < number; i++)
	                {
                        int lineNumber=round(output_values[i*MY_OUPUT_SIZE+1]+1);
	                    DFexp.locate(lineNumber);
	                    DFexp.set(0,output_values[i*MY_OUPUT_SIZE+2]);
	                    DFexp.set(1,output_values[i*MY_OUPUT_SIZE+3]);
	                    DFexp.set(2,output_values[i*MY_OUPUT_SIZE+4]);
	                    DFexp.set(3,output_values[i*MY_OUPUT_SIZE+5]);
	                    DFexp.set(4,output_values[i*MY_OUPUT_SIZE+6]);
	                    DFexp.set(5,output_values[i*MY_OUPUT_SIZE+7]+1);
	                    DFexp.set(6,output_values[i*MY_OUPUT_SIZE+8]);
	                    DFexp.set(7,output_values[i*MY_OUPUT_SIZE+9]);
                    }         
                }
                else if (status.MPI_TAG == TAG_FREEWORKER)
                {
                    MPI_Recv(&tip, 1, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                          MPI_COMM_WORLD, &status);
                    #ifdef DEBUG
                    std::cerr << "Mr received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE << std::endl;
                    std::cerr << "Ms sent TAG_STOP to worker" << status.MPI_SOURCE << std::endl;
                    #endif
                    MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
                    stopTagsSent++;
                }
                else
                {
                    error_exit("Received unknown TAG I quit (master)");
                }           
                  
            }
            //close temperoal file with results               
	        DFexp.write(fn_root + ".doc");


        }
        else //rank !=0
        {
        // Select only relevant part of selfile for this rank
        // job number
        // job size
        // aux variable
        int worker_tip=-1;;
            while (1)
            {
                int jobNumber=0;
                MPI_Send(&worker_tip, 1, MPI_INT, 0, TAG_FREEWORKER, MPI_COMM_WORLD);
                #define DEBUG
                #ifdef DEBUG
                std::cerr << "W" << rank << " " << "sent TAG_FREEWORKER to master " << std::endl;
                #endif
                //get your next task
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
                    std::cerr << "Wr" << rank 
                              << " " << "TAG_STOP" << std::endl;
                    #endif
                    break;
                    }
                if (status.MPI_TAG == TAG_WORKFORWORKER)
                //there is still some work to be done    
                    {
                    //get the jobs number
                    MPI_Recv(input_images, 
                             input_images_size, 
                             MPI_INT, 
                             0, 
                             TAG_WORKFORWORKER, 
                             MPI_COMM_WORLD, 
                             &status);
                    worker_tip =  input_images[0];
                    int number_of_transfered_images =  input_images[1];
                    #ifdef DEBUG
                    std::cerr << "Wr " << rank << " " << "TAG_WORKFROMWORKER" << std::endl;
                    std::cerr << "\n rank, tip, input_images_size " 
                              << rank 
                              << " " 
                              << worker_tip 
                              << " "
                              << input_images[1]
                              << std::endl;
                    for (int i=2; i < number_of_transfered_images+2 ; i++)
                        std::cerr << input_images[i] << " " ;   
                    std::cerr << std::endl;  
                    #endif
                    /////////////
                    processSomeImages(&(input_images[1]),output_values);
                    #ifdef DEBUG
                    std::cerr << "Ws " << rank << " " << "TAG_WORKFROMWORKER" << std::endl;
                    std::cerr << "\n rank, size " 
                              << rank 
                              << " " 
                              << output_values[0]
                              << std::endl;
                    std::cerr << std::endl;    
                    for (int i=0; i < output_values[0]; i++)
                        std::cerr << output_values[i] << " " ;   
                    std::cerr << std::endl;  
                    #endif
                    MPI_Send(output_values,
                             output_values_size,
                             MPI_DOUBLE, 
                             0,
                             TAG_WORKFROMWORKER,
                             MPI_COMM_WORLD);
                    ///////////////
                    }
                else
                    {
                    error_exit("Received unknown TAG I quit (worker)");
                    }           
             }//while(1)
        }//worker    
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

    Prog_mpi_new_projection_matching_prm prm;
    //if (prm.rank == 0)
    {    
        try
        {
            prm.read(argc, argv);
        }

        catch (Xmipp_error XE)
        {
            std::cerr << XE;
            if (prm.rank == 0)
                prm.usage();
            MPI_Finalize();
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (prm.rank != 0)
    {    
        try
        {
            prm.read(argc, argv);
        }

        catch (Xmipp_error XE)
        {
            std::cerr << XE;
            prm.usage();
            MPI_Finalize();
        }
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
