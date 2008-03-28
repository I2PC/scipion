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

#define TAG_FIRSTWORKFORWORKER   -1
#define TAG_WORKFORWORKER   0
#define TAG_STOP   1
#define TAG_WAIT   2
#define TAG_FREEWORKER   3

#define RESULTS_SIZE      500000
char results[RESULTS_SIZE];

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
        
        /** Number of jobs **/
        int numberOfJobs;
        
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
        mpi_job_size=textToInteger(getParameter(argc,argv,"-mpi_job_size","-1"));
        chuck_angular_distance = checkParameter(argc, argv,"-chuck_angular_distance");
        fn_sym = getParameter(argc, argv, "-sym","c1");
    }

    /* Usage ------------------------------------------------------------------- */
    void usage()
    {
        Prog_new_projection_matching_prm::usage();
        std::cerr << " [ -mpi_job_size default=-1]    : Number of images sent to a cpu in a single job \n";
        std::cerr << "                                  10 may be a good value";
        std::cerr << "                                 if  -1 the computer will fill the value for you";
        std::cerr << " [ -chuck_angular_distance N]    : sample the projection sphere with this \n";
        std::cerr << "                                   sampling rate and create subsets of experimental";
        std::cerr << "                                   using the voronoi regions";
        std::cerr << "  [-sym cn]   :One of the 17 possible symmetries in\n"
                  << "                                single particle electronmicroscopy\n"
                  << "                                i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv, "
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
        if (rank == 0) 
        {
            show();
            //read experimental doc file
            //DFexp.read(fn_exp);
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
#ifdef NEVERDEFINED
            if (mpi_job_size != -1)
            {   
                numberOfJobs = ceil((double)(SF.ImgNo())/mpi_job_size);
            }
            else
            {   
                numberOfJobs=nProcs-1;//one node is the master
                mpi_job_size=ceil((double)SF.ImgNo()/numberOfJobs);
            }    
            //only one node will write in the console
#endif
        }
#ifdef NEVERDEFINED
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
        //sent to the master the number of projection directions  alias nr_dir
        if (rank == 1)
        {
            MPI_Send(&nr_dir, 1, MPI_INT, 0, TAG_NUMBEROFPROJECTIONDIRECTIONS, MPI_COMM_WORLD);
        }
        if (rank == 0)
        {
            MPI_Status status;
            MPI_Recv(&nr_dir, 1, MPI_INT, 1, TAG_NUMBEROFPROJECTIONDIRECTIONS,
                     MPI_COMM_WORLD, &status);
//#define DEBUG
#ifdef DEBUG
std::cerr << "nr_dir " <<  nr_dir << std::endl;
#endif
        }
#endif        
    }

    /* Run --------------------------------------------------------------------- */
    void run()
    {   

exit(0);
#ifdef NEVERDEFINED
        if (rank == 0)
        {
            int stopTagsSent =0;
            for (int i=0;i<numberOfJobs;)
            {
                //collect data if available
                //be aware that mpi_Probe will block the program untill a message is received
#ifdef DEBUG
std::cerr << "Mp1 waiting for any  message " << std::endl;
#endif
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#ifdef DEBUG
std::cerr << "Mp2 received tag from worker " <<  status.MPI_SOURCE << std::endl;
#endif
               // croscorrelation coheficient
                if (status.MPI_TAG == TAG_SUMCC)
                   {
                   MPI_Recv(&auxSumCC, 1, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_SUMCC, MPI_COMM_WORLD, &status);
                   sumCcCounter++;
//#define DEBUG
#ifdef DEBUG
std::cerr << "Mr received TAG_SUMCC from worker " <<  status.MPI_SOURCE << std::endl;
std::cerr << "sumCcCounter " <<  sumCcCounter << std::endl;
#endif
#undef DEBUG
                   sumCC += auxSumCC;
//#define DEBUG
#ifdef DEBUG
std::cerr << "sumCC auxSumCC" <<  sumCC << " " <<  auxSumCC << std::endl;
#endif
#undef DEBUG
                   }
                //doc file with angles and shifts

                else if (status.MPI_TAG == TAG_DOCFILE)
                   {
                   int iNumber;
                   MPI_Recv(results, RESULTS_SIZE, MPI_CHAR, 
                                                MPI_ANY_SOURCE, 
                                                TAG_DOCFILE, 
                                                MPI_COMM_WORLD, &status); 
                   MPI_Get_count(&status,MPI_CHAR,&iNumber);
                   results[iNumber]='\0';
#ifdef DEBUG
std::cerr << "Mr received TAG_DOCFILE from worker " <<  status.MPI_SOURCE << std::endl;
#endif
                   docCounter++;
                   myDocFile<<results ;
                   }

                //sel file with images asigned to classes
                //we  have a std::vector because there is a sel file for each class

                else if (status.MPI_TAG == TAG_SELFILE)
                   {
                   int iNumber;
                   MPI_Recv(results, RESULTS_SIZE, MPI_CHAR, 
                                                MPI_ANY_SOURCE, 
                                                TAG_SELFILE, 
                                                MPI_COMM_WORLD, &status); 
                   //extract reference_library number;
                   MPI_Get_count(&status,MPI_CHAR,&iNumber);
                   results[iNumber]='\0';
//#define DEBUG
#ifdef DEBUG
std::cerr << "Mr1 received TAG_SELFILE from worker " <<  status.MPI_SOURCE << std::endl;
//std::cerr << "results value: " <<  results << std::endl;
#endif
#undef DEBUG
                   char auxChar[11], *charPointer;
                   strncpy(auxChar, results ,10);
                   auxChar[10]='\0';
                   charPointer = &(results[10]);
                   int auxInt=textToInteger(auxChar);
                   if (selData.size()< (1+auxInt))
                       selData.resize(1+auxInt,"\0");
                   selData[auxInt].append(charPointer);
                   //we still need to save this sel files   
                   } 
                //allsel files for a particular job have been sent,
                //increase counter
                else if (status.MPI_TAG == TAG_SELEND)
                   {
                   //remove mesage
                   MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_SELEND,
                         MPI_COMM_WORLD, &status);
                   selCounter++;
#ifdef DEBUG
std::cerr << "Mr_f received TAG_SELEND from worker " <<  status.MPI_SOURCE << std::endl;
#endif
                    }
                // worker is free
                else if (status.MPI_TAG == TAG_FREEWORKER)
                   {
                   MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                         MPI_COMM_WORLD, &status);
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
#undef DEBUG
                    }
                 else
                    {
                    std::cerr << "M_f Recived unknown TAG" << std::endl;
                    exit(0);
                    }           
            }
            
            num_img_tot=SF.ImgNo();
            while (sumCcCounter< numberOfJobs 
                  || docCounter < numberOfJobs
                  || selCounter < numberOfJobs
                  ) //add other data here with ||
                {
#ifdef DEBUG
std::cerr << "_sumCcCounter " <<  sumCcCounter << std::endl;
std::cerr << "_docCounter "   <<  docCounter << std::endl;
std::cerr << "_numberOfJobs " <<  numberOfJobs << std::endl;
#endif

#ifdef DEBUG
std::cerr << "Mp2 waiting for any  message " << std::endl;
#endif
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#ifdef DEBUG
std::cerr << "Mp received tag from worker " <<  status.MPI_SOURCE << std::endl;
std::cerr << "Mp tag is " <<  status.MPI_TAG << std::endl;
#endif
                // croscorrelation coheficient
                if (status.MPI_TAG == TAG_SUMCC)
                   {
                   MPI_Recv(&auxSumCC, 1, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_SUMCC, MPI_COMM_WORLD, &status);
//#define DEBUG
#ifdef DEBUG
std::cerr << "Mr received TAG_SUMCC from worker " <<  status.MPI_SOURCE << std::endl;
#endif
#undef DEBUG
                   sumCcCounter++;
                   sumCC += auxSumCC;
//#define DEBUG
#ifdef DEBUG
std::cerr << "sumCC auxSumCC" <<  sumCC << " " <<  auxSumCC << std::endl;
#endif
#undef DEBUG
                   }

                else if (status.MPI_TAG == TAG_DOCFILE)
                   {
                   int iNumber;
                   MPI_Recv(results, RESULTS_SIZE, MPI_CHAR, 
                                                MPI_ANY_SOURCE, 
                                                TAG_DOCFILE, 
                                                MPI_COMM_WORLD, &status);
                   MPI_Get_count(&status,MPI_CHAR,&iNumber);
                   results[iNumber]='\0';
#ifdef DEBUG
std::cerr << "Mr received TAG_DOCFILE from worker " <<  status.MPI_SOURCE << std::endl;
#endif
                   docCounter++;
                   myDocFile<<results ;
                   }

                else if (status.MPI_TAG == TAG_SELFILE)
                   {
                   int iNumber;
                   MPI_Recv(results, RESULTS_SIZE, MPI_CHAR, 
                                                MPI_ANY_SOURCE, 
                                                TAG_SELFILE, 
                                                MPI_COMM_WORLD, &status); 
                   //extract reference_library number;
                   MPI_Get_count(&status,MPI_CHAR,&iNumber);
                   results[iNumber]='\0';
//#define DEBUG
#ifdef DEBUG
std::cerr << "Mr2 received TAG_SELFILE from worker " <<  status.MPI_SOURCE << std::endl;
std::cerr << "results value: " <<  results << std::endl;
#endif
#undef DEBUG
                   char auxChar[11], *charPointer;
                   strncpy(auxChar, results ,10);
                   auxChar[10]='\0';
                   charPointer = &(results[10]);
                   int auxInt=textToInteger(auxChar);
                   if (selData.size()< (1+auxInt))
                       selData.resize(1+auxInt,"\0");
                   selData[auxInt].append(charPointer);   
                   //we still need to save this sel files   
                   } 
                //allsel files for a particular job have been sent,
                //increase counter
                else if (status.MPI_TAG == TAG_SELEND)
                   {
                   //remove mesage
                   MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_SELEND,
                         MPI_COMM_WORLD, &status);
                   selCounter++;
#ifdef DEBUG
std::cerr << "Mr_w received TAG_SELEND from worker " <<  status.MPI_SOURCE << std::endl;
#endif
                    }

                   
                else if (status.MPI_TAG == TAG_FREEWORKER)
                   {
                   MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                         MPI_COMM_WORLD, &status);
#ifdef DEBUG
std::cerr << "Mr_w received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE << std::endl;
std::cerr << "Ms_w sent TAG_STOP to worker" << status.MPI_SOURCE << std::endl;
#endif
                   MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
                   stopTagsSent++;
                   }
                 else
                    {
                    std::cerr << "M_w Recived unknown TAG" << std::endl;
                    exit(0);
                    }           

                }
        //some workers did not got their TAG_STOP
        while (stopTagsSent < (nProcs-1))
        {
            MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                  MPI_COMM_WORLD, &status);
    #ifdef DEBUG
    std::cerr << "Mr received TAG_FREEWORKER from worker " <<  status.MPI_SOURCE << std::endl;
    std::cerr << "Ms sent TAG_STOP to worker" << status.MPI_SOURCE << std::endl;
    #endif
            MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
            stopTagsSent++;
        }         
        sumCC /=  (double) num_img_tot;
        //save doc_file and renumber it
        myDocFile.close();
        DFo.read(fn_tmp);
        DFo.renum();
        DFo.write(fn_tmp);
        {
        //save sel files
        FileName fn_base,fn_sel,fn_xmp;
        fn_base=fn_root+"_class";
        SelFile SFxmp,SFClass;
        SFxmp.clear();
        SFClass.clear();
        std::ofstream mySelFile;
        for (int dirno = 0; dirno < selData.size(); dirno++)
             {
             //compute sel files
             fn_xmp.compose(fn_base,dirno+1,"xmp");
             SFClass.insert(fn_xmp);
             fn_sel.compose(fn_base,dirno+1,"sel");
             mySelFile.open (fn_sel.c_str());
             mySelFile << selData[dirno];
             mySelFile.close();
             //compute averages
             }
        //create empty extra sel files     
        for (int dirno = selData.size();  dirno < nr_dir; dirno++)
             {
             fn_xmp.compose(fn_base,dirno+1,"xmp");
             SFClass.insert(fn_xmp);
             fn_sel.compose(fn_base,dirno+1,"sel");
             mySelFile.open (fn_sel.c_str());
             mySelFile.close();
             }
        fn_base+="es.sel";
        SFClass.write(fn_base);
        }
        // compute averages
        //create empty image
        {
        ImageXmipp empty,img;
        FileName fn_base,fn_sel,fn_doc,fn_average;

        //read projection to get right size
        FileName fn_tmp, fn_refs,fn_img;
        fn_refs=fn_root+"_lib";
        fn_tmp.compose(fn_refs,0+1,"proj");
        empty.read(fn_tmp);
        //creat sel file name
        fn_base=fn_root+"_class";
        DocFile         DF;
        SelFile         SF_tmp;
        fn_doc = fn_root + ".doc";
        DF.read(fn_doc);
        DF.go_beginning();
        //#define DEBUG
        #ifdef DEBUG
        std::cerr << "doc dile name" << fn_doc << std::endl;
        std::cerr << "nmax" << DF.dataLineNo() << std::endl;
        #endif
        #undef DEBUG
        int nmax = DF.dataLineNo();
        for (int dirno = 0;  dirno < nr_dir; dirno++)
            {
            fn_sel.compose(fn_base,dirno+1,"sel");//sel file
            fn_average.compose(fn_base,dirno+1,"xmp");//average file
            SF_tmp.read(fn_sel);
            SF_tmp.go_beginning();
            empty().init_constant(0.);
	    empty.weight() = 0;
            empty.set_originOffsets(0.f,0.f);
            empty.set_eulerAngles1(0.,0.,0.); 
            if(SF_tmp.ImgNo()==0)
                {
                ;
                }
            else
                {
                    while (!SF_tmp.eof())
                    {
                        fn_img = SF_tmp.NextImg();
                        img.read(fn_img);
                        img().selfApplyGeometryBSpline(img.get_transformation_matrix(),3,IS_INV,WRAP);
                        empty() += img();
                        empty.weight() += 1;
                    }
                }
            empty.rot()=img.rot();
	    empty.tilt()=img.tilt();
            empty.write(fn_average);
            }        
        }        
        //for (int dirno = 1; dirno < nr_dir; dirno++)
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
                    break;
                    }
                if (status.MPI_TAG == TAG_WORKFORWORKER)
                //there is still some work to be done    
                    {
                    //get the jobs number
                    MPI_Recv(&jobNumber, 1, MPI_INT, 0, TAG_WORKFORWORKER, MPI_COMM_WORLD, &status);
                    //#define DEBUG1 
                    #ifdef DEBUG1  
                        std::cerr <<    "SF.ImgNo() " << SF.ImgNo() << std::endl;     
                        std::cerr <<    "jobNumber "  << jobNumber << std::endl;  
                        std::cerr <<    "rank "  << rank << std::endl;  
                    #endif
                    auxSF=SF;
                    //create local sel file
                    auxSF.mpi_select_part2(jobNumber, numberOfJobs, num_img_tot,mpi_job_size);
                    //#define DEBUG1
                    #ifdef DEBUG1
                        std::cerr << "Rank " << rank << "chunck " <<  jobNumber << std::endl; 
                        std::cerr <<    "auxSF.ImgNo() " << auxSF.ImgNo() << std::endl;     
                        std::cerr <<    "num_img_tot "   << num_img_tot << std::endl;     
                        char aux_str[5];
                        sprintf(aux_str, "%05d", jobNumber);
                        std::string sAux;
                        sAux=aux_str;
                        FileName   fnAux;
                        fnAux="tmp_"+sAux+".sel";
                        auxSF.write(fnAux);
                    #endif
                    #undef DEBUG1
                    DFo.clear();
                    // Process all images
                    MPI_Request request;
                    //process all the images in local sel
                    PM_loop_over_all_images(auxSF, DFo, sumCC);
                    //work done now send back the results to the master
                    //this is going to be tricky
                    //since sending complex structures is difficult
                    //we wil use the ioerator << to write them in a
                    //string and the send the char array back
                    //the master will reasamble cafully all the data
                    
                    //DOCFILE
                    // redirect standard output to a string
                    std::ostringstream doc;

                    //print docfile to string s
                    doc << DFo;
                    // name and values
                    int s_size=  doc.str().size();
                    if (s_size >= RESULTS_SIZE)
                    {
                        std::cerr << "docfile to long, reduce mpi_job_size" << std::endl;
                        exit(0);
                    }
                    strncpy(results,doc.str().c_str(),s_size);
                    results[s_size]='\0';

                    /* MPI_Isend(results,s_size, MPI_CHAR, 0,TAG_DOCFILE, MPI_COMM_WORLD,&request);*/
                    MPI_Send(results, s_size, MPI_CHAR, 0, TAG_DOCFILE, MPI_COMM_WORLD);

                    //SELFILE
                    // redirect standard output to a string
                    //note that there is a sel file per projectin direction
                    
                    for (int dirno = 0; dirno < nr_dir; dirno++)
                    {   //do not bother about classes without assigned images
                        if(class_selfiles[dirno].ImgNo()==0)
                             continue;
                        std::ostringstream sel;
                        //#define DEBUG2
                        #ifdef DEBUG2
                        FileName fn_img;
                        fn_img.compose("kk",rank*1000+dirno+1,"sel");
                        class_selfiles[dirno].write(fn_img);
                        //std::cerr << fn_img << class_selfiles[dirno] <<std::endl;
                        #endif
                        #undef DEBUG2
                       //first 10 characters are the projection library number
                        sel << std::setw(10) << dirno ;
                        sel << class_selfiles[dirno];
                        // name and values
                        int s_size=  sel.str().size();
                        if (s_size >= RESULTS_SIZE)
                        {
                            std::cerr << "selfile to long, reduce mpi_job_size" << std::endl;
                            std::cerr << " (or increase RESULTS_SIZE and recompile) " ;
                            exit(0);
                        }
                        //#define  DEBUG_TAG_SELFILE
                        #ifdef   DEBUG_TAG_SELFILE 
                        std::cerr << "Ws-" << rank << " sel" << sel.str().c_str() << std::endl;
                        #endif                
                        #undef   DEBUG_TAG_SELFILE                      
                        strncpy(results,sel.str().c_str(),s_size);
                        //this is not needed since we only pass 
                        //s_size but I like it
                        results[s_size]='\0';
                        //#define  DEBUG_TAG_SELFILE
                        #ifdef   DEBUG_TAG_SELFILE 
                        std::cerr << "Ws-" << rank << " results" << results << std::endl;
                        #endif                
                        #undef   DEBUG_TAG_SELFILE                      
                        MPI_Send(results, s_size, MPI_CHAR, 0, TAG_SELFILE, MPI_COMM_WORLD);
                    }
                    //All sel files created for this job have been sent                
                    MPI_Send(0, 0, MPI_INT, 0, TAG_SELEND, MPI_COMM_WORLD);
                    //CC COHEFICIENT
                    // may be send no bloking sumCC
                    //that will be faster but the program logic became more complex
                    /*MPI_Isend(&sumCC, 1, MPI_DOUBLE, 0, TAG_SUMCC, MPI_COMM_WORLD,&request);*/
                    MPI_Send(&sumCC, 1, MPI_DOUBLE, 0, TAG_SUMCC, MPI_COMM_WORLD);
                    //#define DEBUGCC
                    #ifdef DEBUGCC
                    std::cerr << "Ws_" << rank << " sumCC" << sumCC << std::endl; 
                    #endif
                    #undef DEBUGCC
                    //get yor next task
                }
                else
                   {
                   std::cerr << "3) Recived unknown TAG I quit" << std::endl;
                   exit(0);
                   }           
            }
        }
        MPI_Finalize();
#endif
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
            prm.usage();
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
