/***************************************************************************
 *
 * Authors:  Slavica Jonic slavica.jonic@impmc.jussieu.fr  
 *           Carlos Oscar Sanchez Sorzano coss.eps@ceu.es
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 * Lab. de Bioingenieria, Univ. San Pablo CEU
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
//#include <mpi.h>
#include <reconstruction/nma_alignment.h>
#include <parallel/mpi.h>

/** Class to perfom the NMA Alignment with  MPI parallelization */
class MpiProgNMA: public ProgNmaAlignment
{
private:
    MpiNode *node;
    FileTaskDistributor *distributor;
    std::vector<long int> imgsId;

public:
    /** Destructor */
    ~MpiProgNMA()
    {
        delete node;
    }

    /** Redefine read to initialize MPI environment */
    void read(int argc, char **argv)
    {
        node = new MpiNode(argc, argv);
        ProgNmaAlignment::read(argc, argv);
    }
    /** main body */
    void createWorkFiles()
    {
      //Master node should prepare some stuff before start working
      if (node->isMaster())
      {
        ProgNmaAlignment::createWorkFiles();
        mdIn.write("nmaTodo.xmd");
      }
      node->barrierWait();//Sync all before start working
      mdIn.read("nmaTodo.xmd");
      mdIn.findObjects(imgsId);//get objects ids
      rangen = node->rank;
      std::cerr << "Creating file task distributor......."<< std::endl;
      distributor = new FileTaskDistributor(mdIn.size(), 1, node);
    }
    //Only master do starting progress bar stuff
    void startProcessing()
    {
      if (node->isMaster())
        ProgNmaAlignment::startProcessing();
    }
    //Only master do finishing progress bar stuff
    void finishProcessing()
    {
      if (node->isMaster())
        ProgNmaAlignment::finishProcessing();
    }
    //Only master show progress
    void showProgress()
    {
      if (node->isMaster())
        ProgNmaAlignment::showProgress();
    }
    //Now use the distributor to grasp images
    long int getImageToProcess()
    {
      longint first, last;
      bool moreTasks = distributor->getTasks(first, last);
      if (moreTasks)
      {
        time_bar_done = first + 1;
        return imgsId[first];
      }
      time_bar_done = mdIn.size();
      return -1;
    }

    void postProcess()
    {
      //All nodes wait for each other
      node->barrierWait();
    }
}
;//end of class MpiProgNMA

int main(int argc, char **argv)
{
  MpiProgNMA program;
  std::cerr << "reading " << std::endl;
  program.read(argc, argv);
  std::cerr << "running " << std::endl;
  program.tryRun();

    // Initialize MPI
  /*
    int rank, NProcessors;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &NProcessors);

    ProgNmaAlignment prm;
    try
    {
        // Read input parameters
        prm.MPIversion=true;
        prm.read(argc, argv);
    }
    catch (XmippError XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    try
    {
        // Prepare side information
        prm.preProcess();
        MetaData SF_in(prm.fn_in);

        // Divide the selfile in chunks
        int imgNbr = prm.fn_in.size();
        int Nchunk = (int)((float)imgNbr / (float)(NProcessors - 1));
        int myFirst = (rank - 1) * Nchunk;
        int myLast = rank * Nchunk - 1;
        if (rank == NProcessors - 1)
            myLast = imgNbr - 1;

        // Make the alignment, rank=0 receives all the assignments
        // The rest of the ranks compute the angular parameters for their
        // assigned images

        int numberofparam = 7 + prm.modeList.size();

        double v[numberofparam];
        if (rank == 0)
        {
            int i=0;
            FOR_ALL_OBJECTS_IN_METADATA(SF_in)
            {
                FileName fnImg;
                SF_in.getValue(MDL_IMAGE,fnImg);
                prm.img_names.push_back(fnImg);
                Matrix1D<double> dummy;
                prm.assignments.push_back(dummy);
                i++;
            }
            int toGo = imgNbr;
            MPI_Status status;
            //std::cerr << "Assigning modes and angles ...\n";
            init_progress_bar(imgNbr);
            while (toGo > 0)
            {
                MPI_Recv(v, numberofparam, MPI_DOUBLE, MPI_ANY_SOURCE,
                         MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                int i = (int)v[0];
                Matrix1D<double> aux(numberofparam-1);

                for (int j=1; j<numberofparam; j++)
                    aux(j-1)=v[j];
                prm.assignments[i]=aux;
                toGo--;
            }
            progress_bar(imgNbr);
        }
        else
        {
            for (int i = myFirst; i <= myLast; i++)
            {
                // Read image and estimate angular parameters
                Image<double> I;
                FileName tempname;
                SF_in.getValue(MDL_IMAGE,tempname,i+1);

                //Old (no metadata etc) : I.read(tempname, false, false, false);
                //I.read(tempname);

                I().setXmippOrigin();
                prm.processImage(tempname, prm.fnOut, SF_in.getActiveObject());

                // Send the alignment parameters to the master
                v[0] = i;
                for (int j=1; j<numberofparam; j++)
                {
                    v[j]=prm.parameters(j-1);
                }
                MPI_Send(v, numberofparam, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }

        if (rank == 0)
            prm.postProcess();
        MPI_Finalize();
        return 0 ;
    }

    catch (XmippError XE)
    {
        std::cout << XE << std::endl;
        MPI_Finalize();
        return 1 ;
    } */
}
