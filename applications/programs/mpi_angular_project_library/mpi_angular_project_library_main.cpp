/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

//#include "mpi_run.h"

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <reconstruction/angular_project_library.h>

#define TAG_WORKFORWORKER   0
#define TAG_STOP   1
#define TAG_WAIT   2
#define TAG_FREEWORKER   3

#define DEBUG
class ProgMpiAngularProjectLibrary: public ProgAngularProjectLibrary
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
    ProgMpiAngularProjectLibrary()
    {
        //parent class constructor will be called by deault without parameters
        MPI_Comm_size(MPI_COMM_WORLD, &(nProcs));
        MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
        //Blocks until all process have reached this routine.
        //very likelly this is
        MPI_Barrier(MPI_COMM_WORLD);
    }


    /* Read parameters --------------------------------------------------------- */
    void readParams()
    {
        if (nProcs < 2)
            error_exit("This program cannot be executed in a single working node");
        ProgAngularProjectLibrary::readParams();
        mpi_job_size = getIntParam("--mpi_job_size");
    }

    /* Usage ------------------------------------------------------------------- */
    void defineParams()
    {
        ProgAngularProjectLibrary::defineParams();
        addParamsLine("  [--mpi_job_size <size=10>]: Number of images sent to a cpu in a single job");
        addParamsLine("                            : 10 may be a good value");
        addParamsLine("                            : if  -1 the computer will put the maximum");
        addParamsLine("                            : posible value that may not be the best option");
    }


    /* Show -------------------------------------------------------------------- */
    void show()
    {
        ProgAngularProjectLibrary::show();
        std::cerr << " Size of mpi jobs " << mpi_job_size <<std::endl;
    }

    /* Pre Run PreRun for all nodes but not for all works */
    void preRun()
    {
        //do not use rank 0 here, since rank 1 will manipulate this file later
        if (rank == 1)
        {
            unlink(output_file.c_str());
        }


        //#define DEBUGTIME
#ifdef  DEBUGTIME
#include <ctime>

        time_t start,end;
        double time_dif;
        time (&start);
        time (&end);
        time_dif = difftime (end,start);
        start=end;
        std::cerr<<" starting prerun rank= "<<rank<<std::endl;
#endif

        int my_seed;
        if (rank == 0)
        {
            show();
            //randon numbers must be the same in all nodes
            srand ( time(NULL) );
            if(perturb_projection_vector!=0)
            {
                my_seed=rand();
            }
        }
#ifdef  DEBUGTIME
        time (&end);
        time_dif = difftime (end,start);
        start=end;
        std::cerr<<" set rand seed rank= "<<rank<<std::endl;
#endif
        //Bcast must be seem by all processors
        if(perturb_projection_vector!=0)
        {
            MPI_Bcast (&my_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
            mysampling.setNoise(perturb_projection_vector,my_seed);
#ifdef  DEBUGTIME

            time (&end);
            time_dif = difftime (end,start);
            start=end;
            std::cerr<<" after perturb rank= "<<rank<<std::endl;
#endif

        }
        //all ranks
        mysampling.setSampling(sampling);
        //symmetry for sampling may be different from neighbourhs
        if (!mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
            REPORT_ERROR(ERR_NUMERICAL, (std::string)"angular_project_library::run Invalid symmetry" +  fn_sym);//set sampling must go before set noise
        if(angular_distance_bool!=0)
            mysampling.setNeighborhoodRadius(angular_distance);//irelevant

#ifdef  DEBUGTIME

        time (&end);
        time_dif = difftime (end,start);
        start=end;
        std::cerr<<" setsampling rank= "<<rank<<std::endl;
#endif
        //true -> half_sphere
        mysampling.computeSamplingPoints(false,max_tilt_angle,min_tilt_angle);
        mysampling.SL.readSymmetryFile(fn_sym);
        //store symmetry matrices, this is faster than computing them each time
        mysampling.fillLRRepository();

#ifdef  DEBUGTIME

        time (&end);
        time_dif = difftime (end,start);
        start=end;
        std::cerr<<" compute sampling points rank= "<<rank<<std::endl;
#endif

        // We first sample The  whole sphere
        // Then we remove point redundant due to sampling symmetry
        // use old symmetry, this is geometric does not use L_R
        mysampling.removeRedundantPoints(symmetry, sym_order);

        //=========================
        //======================
        //recompute symmetry with neigh symmetry
        if (!mysampling.SL.isSymmetryGroup(fn_sym_neigh, symmetry, sym_order))
            REPORT_ERROR(ERR_NUMERICAL, (std::string)"angular_project_library::run Invalid neig symmetry" +  fn_sym_neigh);
        mysampling.SL.readSymmetryFile(fn_sym_neigh);
        mysampling.fillLRRepository();
        //precompute product between symmetry matrices and experimental data
        if (FnexperimentalImages.size() > 0)
            mysampling.fillExpDataProjectionDirectionByLR(FnexperimentalImages);
#ifdef  DEBUGTIME

        time (&end);
        time_dif = difftime (end,start);
        start=end;
        std::cerr<<" remove redundant rank= "<<rank<<std::endl;
#endif

        if (FnexperimentalImages.size() > 0 &&
            remove_points_far_away_from_experimental_data_bool)
        {
            mysampling.removePointsFarAwayFromExperimentalData();
#ifdef  DEBUGTIME

            time (&end);
            time_dif = difftime (end,start);
            start=end;
            std::cerr<<" remove points far away rank= "<<rank<<std::endl;
#endif

        }

        /* save files */
        if (rank == 0)
        {
            if(compute_closer_sampling_point_bool)
            {
                //find sampling point closer to experimental point (only 0) and bool
                //and save docfile with this information
                mysampling.findClosestSamplingPoint(FnexperimentalImages,output_file_root);
            }
            //mysampling.createSymFile(symmetry, sym_order);
            mysampling.createAsymUnitFile(output_file_root);
        }

        if (rank != 0)
        {
            try
            {
                inputVol.read(input_volume);
            }
            catch (XmippError &XE)
            {
                std::cout << XE;
                error_exit("Error reading reference volume\n\n");
            }
            inputVol().setXmippOrigin();
            Xdim = XSIZE(inputVol());
            Ydim = YSIZE(inputVol());
        }
        if (rank == 0)
        {
            if (compute_neighbors_bool)
            {
                mysampling.computeNeighbors(only_winner);
                mysampling.saveSamplingFile(output_file_root,false);
            }
        }
        //release some memory
        mysampling.exp_data_projection_direction_by_L_R.clear();

        if (mpi_job_size != -1)
        {
            numberOfJobs = (int)ceil((double)(mysampling.no_redundant_sampling_points_angles.size())/mpi_job_size);
        }
        else
        {
            numberOfJobs=nProcs-1;//one node is the master
            mpi_job_size=(int)ceil((double)mysampling.no_redundant_sampling_points_angles.size()/numberOfJobs);
        }

        verbose=false;
        #define DEBUG
#ifdef DEBUG

        if (rank == 1)
        {
            std::cerr << "numberOfJobs: " << numberOfJobs << std::endl
            << "number of projections to be created: " <<  mysampling.no_redundant_sampling_points_angles.size()
            <<std::endl;
        }
#endif
        #undef DEBUG

        //create blanck outfile so main header does not need to be re-writen
        //rank 1 has already xdim and ydim.
        if (rank == 1)
        {
            int  numberProjections=0;
            for (double mypsi=0;mypsi<360;mypsi += psi_sampling)
            {
                    ++numberProjections;
            }
            numberProjections *= mysampling.no_redundant_sampling_points_angles.size();
            std::cerr   << "creating Blank file: "
            << output_file << " for "
            << numberProjections << "  projections."
            << std::endl;
            createEmptyFile(output_file,Xdim,Ydim,1,numberProjections,true,WRITE_REPLACE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    /* Run --------------------------------------------------------------------- */
    void process()
    {
        if (rank == 0)
        {
            int stopTagsSent =0;
            int c = XMIPP_MAX(1, numberOfJobs / 60);
            init_progress_bar(numberOfJobs);
            for (int i=0;i<numberOfJobs;)
            {
                //collect data if available
                //be aware that mpi_Probe will block the program until a message is received
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

                if (i % c == 0)
                    progress_bar(i);
            }
            progress_bar(numberOfJobs);

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
            //only rank 0 create sel file
            if(rank==0)
            {

                MetaData  mySFin;
                mySFin.read(output_file_root + "_angles.doc");
#define ANGLESDOC
#ifdef ANGLESDOC
//       std::cerr << "DEBUG_ROB, output_file_root + angle.doc:"
//    		     << output_file_root + "_angles.doc" << std::endl;
#endif
                MetaData  mySF;
                FileName fn_temp;

                size_t myCounter = 0;
                size_t id;
                int ref;
                for (double mypsi=0;mypsi<360;mypsi += psi_sampling)
                    //for (int i=0;i<=mysampling.no_redundant_sampling_points_angles.size()-1;i++)
                    FOR_ALL_OBJECTS_IN_METADATA(mySFin)
                {

                    double x,y,z, rot, tilt, psi;
                    mySFin.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
                    mySFin.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
                    mySFin.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
                    mySFin.getValue(MDL_X,x,__iter.objId);
                    mySFin.getValue(MDL_Y,y,__iter.objId);
                    mySFin.getValue(MDL_Z,z,__iter.objId);
                    mySFin.getValue(MDL_REF,ref,__iter.objId);

                    //FIXME, do I have order?
                    fn_temp.compose( ++myCounter,output_file);
                    id = mySF.addObject();
                    mySF.setValue(MDL_IMAGE,fn_temp, id);
                    mySF.setValue(MDL_ENABLED,1, id);

                    mySF.setValue(MDL_ANGLE_ROT,rot, id);
                    mySF.setValue(MDL_ANGLE_TILT,tilt, id);
                    mySF.setValue(MDL_ANGLE_PSI,psi+mypsi, id);
                    mySF.setValue(MDL_X,x, id);
                    mySF.setValue(MDL_Y,y, id);
                    mySF.setValue(MDL_Z,z, id);
                    mySF.setValue(MDL_SCALE,1.0,id);
                    mySF.setValue(MDL_REF,ref,id);
                }
                fn_temp=output_file_root+".doc";
                mySF.setComment("x,y,z refer to the coordinates of the unitary vector at direction given by the euler angles");
                mySF.write(fn_temp);
#ifndef ANGLESDOC
                unlink((output_file_root+"_angles.doc").c_str());
#endif
#undef ANGLESDOC
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
#ifdef DEBUG

                    std::cerr << "Wr" << rank << " " << "TAG_WORKFORWORKER" << std::endl;
                    std::cerr <<    "jobNumber "  << jobNumber << std::endl;
#endif
#undef DEBUG
                    // Process all images
                    project_angle_vector(jobNumber*mpi_job_size,
                                         std::min((size_t)((jobNumber+1)* mpi_job_size -1),
                                                   mysampling.no_redundant_sampling_points_angles.size()-1),
                                         verbose);
                    //get yor next task
                }
                else
                {
                    std::cerr << "3) Received unknown TAG I quit" << std::endl;
                    exit(0);
                }
            }
        }
    }

    void createGroupSamplingFiles(void)
    {
        if (rank==0 && fn_groups!="")
            ProgAngularProjectLibrary::createGroupSamplingFiles();
    }

    /* a short function to print a message and exit */
    void error_exit(const char * msg)
    {
        fprintf(stderr, "%s", msg);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    void run()
    {
      try
      {
          preRun();
          process();
          createGroupSamplingFiles();
      }
      catch (XmippError &XE)
      {
          std::cerr << "Error!" <<std::endl;
          std::cerr << XE;
          MPI_Finalize();
          exit(1);
      }

      if (rank == 0)
          std::cout << "Done!" <<std::endl;
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

    ProgMpiAngularProjectLibrary program;
    if (program.rank == 0)
    {
        try
        {
            program.read(argc, argv);
        }
        catch (XmippError &XE)
        {
            std::cerr << XE;
            MPI_Finalize();
            exit(1);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (program.rank != 0)
    {
        try
        {
            program.read(argc, argv);
        }
        catch (XmippError &XE)
        {
            std::cerr << XE;
            MPI_Finalize();
            exit(1);
        }
    }

    program.tryRun();
    MPI_Finalize();
    exit(0);

}


