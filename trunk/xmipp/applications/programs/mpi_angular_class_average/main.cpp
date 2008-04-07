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


#include <data/args.h>
#include <data/image.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <reconstruction/angular_class_average.h>

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

class Prog_mpi_angular_class_average:Prog_angular_class_average_prm
{
    public:
    //int rank, size, num_img_tot;

        /** Number of Procesors **/
        int nProcs;
        
        /** Dvide the job in this number block with this number of images */
        //int mpi_job_size;
               
        /** computing node number. Master=0 */
        int rank;

        /** status after am MPI call */
        MPI_Status status;
                
    /*  constructor ------------------------------------------------------- */
    Prog_mpi_angular_class_average()
    {
        //parent class constructor will be called by deault without parameters
        MPI_Comm_size(MPI_COMM_WORLD, &(nProcs));
        MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
        if (nProcs < 2)
            error_exit("This program cannot be executed in a single working node");
        //Blocks until all process have reached this routine.
        //very likelly this is useless
        MPI_Barrier(MPI_COMM_WORLD);
    }


    /* Read parameters --------------------------------------------------------- */
    void read(int argc, char **argv)
    {
        Prog_angular_class_average_prm::read(argc,argv);
        //mpi_job_size=textToInteger(getParameter(argc,argv,"-mpi_job_size","10"));
    }

    /* Usage ------------------------------------------------------------------- */
    void usage()
    {
        Prog_angular_class_average_prm::usage();
        //std::cerr << " [ -mpi_job_size default=-1]    : Number of images sent to a cpu in a single job \n";
        //std::cerr << "                                  10 may be a good value\n";
        //std::cerr << "                                 if  -1 the computer will fill the value for you\n";
     }


    /* Show -------------------------------------------------------------------- */
    void show()
    {
        Prog_angular_class_average_prm::show();
	    //std::cerr << " Size of mpi jobs " << mpi_job_size <<std::endl
        //      ;
    }

    /* Pre Run --------------------------------------------------------------------- */
    void preRun()
    {
        produceSideInfo();
//        MPI_Bcast(&max_number_of_images_in_around_a_sampling_point, 
//                  1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    /* Run --------------------------------------------------------------------- */
    void run()
    {   
    
        double myw[4];
        int number_of_references_image=1;
        if (rank == 0)
        {
	        int nr_ref;
            SelFile          SFclasses, SFclasses1, SFclasses2;
            FileName         fn_tmp;

	        SFclasses.clear(); 
	        SFclasses1.clear(); 
	        SFclasses2.clear(); 
            nr_ref = DFlib.dataLineNo();
	        init_progress_bar(nr_ref);  



            int stopTagsSent =0;
	        init_progress_bar(nr_ref);// master
            while(1)
            {
                //Wait until any message arrives
                //be aware that mpi_Probe will block the program untill a message is received
                #define DEBUG
                #ifdef DEBUG
                std::cerr << "Mp1 waiting for any  message " << std::endl;
                #endif
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                #ifdef DEBUG
                std::cerr << "Mp2 received tag from worker " <<  status.MPI_SOURCE << std::endl;
                #endif
                // worker is free
                // worker sends work
                if (status.MPI_TAG == TAG_WORKFROMWORKER)
                {
                    MPI_Recv(myw, 
                             4, 
                             MPI_DOUBLE, 
                             MPI_ANY_SOURCE, 
                             TAG_WORKFROMWORKER,
                             MPI_COMM_WORLD, 
                             &status);
                    double w, w1,w2;
                    int myref_number;         
                    w  = myw[0];
                    w1 = myw[1];
                    w2 = myw[2]; 
                    myref_number = round(myw[3]);        
                    #ifdef DEBUG
                    std::cerr << "Mr2.5 received work from worker " <<  status.MPI_SOURCE << std::endl;
                    std::cerr << " w w1 w2 myref_number" 
                              << w  << " "
                              << w1 << " "
                              << w2 << " "
                              << myref_number <<  std::endl;  
                    #endif
	                if (w > 0.)
	                {
		                fn_tmp.compose(fn_out,myref_number,"xmp");
		                SFclasses.insert(fn_tmp);
	                }
	                if (do_split)
	                {
		                if (w1 > 0.)
		                {
		                    fn_tmp.compose(fn_out1,number_of_references_image,"xmp");
		                    SFclasses1.insert(fn_tmp);
		                }
		                if (w2 > 0.)
		                {
		                    fn_tmp.compose(fn_out2,number_of_references_image,"xmp");
		                    SFclasses2.insert(fn_tmp);
		                }
	                }
                }//TAG_WORKFROMWORKER
                // worker is free
                if (status.MPI_TAG == TAG_FREEWORKER)
                {
                    MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                         MPI_COMM_WORLD, &status);
                    #ifdef DEBUG
                    std::cerr << "Mr3 received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE 
                              << std::endl;
                    #endif
                    if(number_of_references_image>nr_ref)
                    {
                        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
                        stopTagsSent++;
                        #ifdef DEBUG
                        std::cerr << "Ms4 sent stop tag to worker " <<  status.MPI_SOURCE << std::endl
                                 << std::endl;
                        #endif
                        break;
                    }
                    else
                    {
	                    //DFlib.adjust_to_data_line();
                        //dataforworker=DFlib(ABS(col_ref) - 1);
                        MPI_Send(&number_of_references_image, 
                                 1, 
                                 MPI_INT, 
                                 status.MPI_SOURCE, 
                                 TAG_WORKFORWORKER, 
                                 MPI_COMM_WORLD);
	                    //prm.DFlib.next();
                     
                    }   
            #ifdef DEBUG
            std::cerr << "Ms5 sent TAG_WORKFORWORKER for " <<  status.MPI_SOURCE << std::endl
                      << std::endl;
            #endif
            number_of_references_image++;//////////////////////
            progress_bar(number_of_references_image);
                }//TAG_FREEWORKER
            }//while       
            progress_bar(nr_ref);


            while (stopTagsSent < (nProcs-1))
            {
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == TAG_WORKFROMWORKER)
                {
                    MPI_Recv(myw, 
                             4, 
                             MPI_DOUBLE, 
                             MPI_ANY_SOURCE, 
                             TAG_WORKFROMWORKER,
                             MPI_COMM_WORLD, 
                             &status);
                    double w, w1,w2;
                    int myref_number;         
                    w  = myw[0];
                    w1 = myw[1];
                    w2 = myw[2];         
                    myref_number = round(myw[3]);        
                    #ifdef DEBUG
                    std::cerr << "Mr2.5 received work from worker " <<  status.MPI_SOURCE << std::endl;
                    std::cerr << " w w1 w2 myref_number" 
                              << w  << " "
                              << w1 << " "
                              << w2 << " "
                              << myref_number <<  std::endl;  
                    #endif
	                if (w > 0.)
	                {
		                fn_tmp.compose(fn_out,myref_number,"xmp");
		                SFclasses.insert(fn_tmp);
	                }
	                if (do_split)
	                {
		                if (w1 > 0.)
		                {
		                    fn_tmp.compose(fn_out1,number_of_references_image,"xmp");
		                    SFclasses1.insert(fn_tmp);
		                }
		                if (w2 > 0.)
		                {
		                    fn_tmp.compose(fn_out2,number_of_references_image,"xmp");
		                    SFclasses2.insert(fn_tmp);
		                }
	                }
                }
                else if (status.MPI_TAG == TAG_FREEWORKER)
                {
                    MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                          MPI_COMM_WORLD, &status);
                    #ifdef DEBUG
                    std::cerr << "Mr received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE << std::endl;
                    #endif
                    MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
                    stopTagsSent++;
                }    
            }         

	        SelFile          auxSF;
            fn_tmp=fn_out+"es.sel";
            auxSF=SFclasses.sort_by_filenames();
	        auxSF.write(fn_tmp);
	        if (do_split)
	        {
	            fn_tmp=fn_out1+"es.sel";
                auxSF=SFclasses1.sort_by_filenames();
	            auxSF.write(fn_tmp);
	            fn_tmp=fn_out2+"es.sel";
                auxSF=SFclasses2.sort_by_filenames();
	            auxSF.write(fn_tmp);
	        }

        }
        else //rank !=0
        {
        // Select only relevant part of selfile for this rank
        // job number
        // job size
        // aux variable
            while (1)
            {
                MPI_Send(0, 0, MPI_INT, 0, TAG_FREEWORKER, MPI_COMM_WORLD);
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
                    break;
                    }
                if (status.MPI_TAG == TAG_WORKFORWORKER)
                //there is still some work to be done    
                    {
                    //get the jobs number
                    MPI_Recv(&number_of_references_image, 
                             1, 
                             MPI_INT, 
                             0, 
                             TAG_WORKFORWORKER, 
                             MPI_COMM_WORLD, 
                             &status);
                    #ifdef DEBUG
                    std::cerr << "Wr" << rank << " " << "TAG_WORKFORWORKER" << std::endl;
                    #endif
                    processOneClass(number_of_references_image, 
                                         myw[0], 
                                         myw[1], 
                                         myw[2]);//slaves
                    myw[3]= (double)number_of_references_image;                    
                    MPI_Send(myw,
                             4,
                             MPI_DOUBLE, 
                             0,
                             TAG_WORKFROMWORKER,
                             MPI_COMM_WORLD);
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

    Prog_mpi_angular_class_average prm;
    if (prm.rank == 0)
    {    
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
            exit(1);
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
