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

#include <parallel/xmipp_mpi.h>
#include <data/xmipp_funcs.h>
#include <data/xmipp_program.h>

//Tags already defined in xmipp
//#define TAG_WORK                     0
//#define TAG_STOP                     1
#define TAG_I_FINISH_WRITTING        12
#define TAG_MAY_I_WRITE              13
#define TAG_YES_YOU_MAY_WRITE        14
#define TAG_DO_NOT_DARE_TO_WRITE     15
#define TAG_I_AM_FREE                16

class ProgMpiAngularClassAverage: public XmippProgram
{
public:
    /**Number of job */
    int nJobs;

    /** status after am MPI call */
    MPI_Status status;

    // Mpi node
    MpiNode *node;

    ProgMpiAngularClassAverage(int argc, char **argv)
    {
        node = new MpiNode(argc, argv);
        if (!node->isMaster())
            verbose = 0;
    }

    void readParams()
    {
        nJobs = getIntParam("--nJobs");
        if (node->size < 2)
            REPORT_ERROR(ERR_ARG_INCORRECT,
                         "This program cannot be executed in a single working node");
    }

    void defineParams()
    {
        addUsageLine("I do not know");
        addParamsLine(
            "--nJobs <nJobs=1000>    : File with commands in different lines");
    }

    void show()
    {
        if (!verbose)
            return;
        std::cout << "number of jobs...: " << nJobs << std::endl;
    }

    /* Run --------------------------------------------------------------------- */
#define MAX_LINE 2048
    char szline[MAX_LINE + 1];
    void run()
    {
        //number of jobs
        //Lock structure
        bool* lockArray=new bool[nJobs+1];
        int lockIndex;
        for (lockIndex = 0; lockIndex < nJobs; ++lockIndex)
            lockArray[lockIndex]=false;
#define ArraySize 4

        int * Def_3Dref_2Dref_JobNo = new int[ArraySize];

        if (node->rank == 0)
        {
            //for (int iCounter = 0; iCounter < nJobs; )//increase counter after I am free
            int iCounter = 0;
            int finishedNodes = 1;
            bool whileLoop = true;

            while (whileLoop)
            {

                //wait until a worker is available
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                switch (status.MPI_TAG)
                {
                case TAG_I_AM_FREE:
                    //Some test values for defocus, 3D reference and projection direction
                    Def_3Dref_2Dref_JobNo[0] = rnd_unif( 0, 9);
                    Def_3Dref_2Dref_JobNo[1] = rnd_unif(10,19);
                    Def_3Dref_2Dref_JobNo[2] = rnd_unif(20,29);
                    Def_3Dref_2Dref_JobNo[3] = iCounter;

                    //increase counter after sending work
                    ++iCounter;

                    //read worker call just to remove it from the queue
                    MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_I_AM_FREE, MPI_COMM_WORLD,
                             &status);
                    std::cerr << "Master: recv TAG_I_AM_FREE from node "
                    << status.MPI_SOURCE<<std::endl;
                    if(iCounter < nJobs)
                    {
                        //send work, first int defocus, second 3D reference, 3rd projection
                        // direction and job number
                        MPI_Send(Def_3Dref_2Dref_JobNo, ArraySize, MPI_INT, status.MPI_SOURCE,
                                 TAG_WORK, MPI_COMM_WORLD);
                        std::cerr << "Master_a: send TAG_WORK to node "
                        << status.MPI_SOURCE<<std::endl;
                    }
                    else
                    {
                        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
                        finishedNodes ++;
                        std::cerr << "finishednodes node->size"
                        		<< finishedNodes << " " << node->size <<std::endl;
                        if (finishedNodes >= node->size)
                        	whileLoop=false;
                    }
                    break;
                case TAG_MAY_I_WRITE:
                    //where do you want to write?
                    std::cerr << "Master: recv TAG_MAY_I_WRITE"
                    << " from node " << status.MPI_SOURCE <<std::endl;
                    MPI_Recv(&lockIndex, 1, MPI_INT, MPI_ANY_SOURCE, TAG_MAY_I_WRITE,
                             MPI_COMM_WORLD, &status);
                    if (lockArray[lockIndex])
                    {//Locked
                        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE,
                                 TAG_DO_NOT_DARE_TO_WRITE, MPI_COMM_WORLD);
                    }
                    else
                    {//Unlocked
                        lockArray[lockIndex]=true;
                        std::cerr << "Master: send TAG_YES_YOU_MAY_WRITE"
                        << " to node " << status.MPI_SOURCE <<std::endl;
                        MPI_Send(&lockIndex, 1, MPI_INT, status.MPI_SOURCE,
                                 TAG_YES_YOU_MAY_WRITE, MPI_COMM_WORLD);
                    }
                    break;
                case TAG_I_FINISH_WRITTING:
                    //release which lock?
                    std::cerr << "MASTER recv TAG_I_FINISH_WRITTING" <<std::endl;
                    MPI_Recv(&lockIndex, 1, MPI_INT, MPI_ANY_SOURCE, TAG_I_FINISH_WRITTING,
                             MPI_COMM_WORLD, &status);
                    lockArray[lockIndex]=false;
                    break;
                default:
                    std::cerr << "WRONG TAG RECEIVED" << std::endl;
                    break;
                }
            }
            std::cerr << "out while" <<std::endl;
        }
        else
        {
            bool whileLoop = true;
            while (whileLoop)
            {
                //I am free
                MPI_Send(0, 0, MPI_INT, 0, TAG_I_AM_FREE, MPI_COMM_WORLD);
                //wait for message
                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                switch (status.MPI_TAG)
                {
                case TAG_STOP://I am free
                    std::cerr << "Bye at " << node->rank << std::endl;
                    MPI_Recv(0, 0, MPI_INT, 0, TAG_STOP,
                             MPI_COMM_WORLD, &status);
                    whileLoop=false;
                    break;
                case TAG_WORK://work to do
                    MPI_Recv(Def_3Dref_2Dref_JobNo, ArraySize, MPI_INT, 0, TAG_WORK,
                             MPI_COMM_WORLD, &status);
                    process(Def_3Dref_2Dref_JobNo);
                    break;
                default:
                    std::cerr << "WRONG TAG RECEIVED at " << node->rank
                    << " tag is " << status.MPI_TAG << std::endl;
                    break;
                }
            }
        }
        std::cerr << "finalize:"  << node->rank <<std::endl;
        MPI_Finalize();
    }
    void process(int * Def_3Dref_2Dref_JobNo)
    {
        std::cerr
        << " Sprocessing def "  << Def_3Dref_2Dref_JobNo[0]
        << " S3Dref "           << Def_3Dref_2Dref_JobNo[1]
        << " S2D_ref "          << Def_3Dref_2Dref_JobNo[2]
        << " SJobNo "           << Def_3Dref_2Dref_JobNo[3]
        << " Sat node: "        << node->rank
        << std::endl;

        usleep(1000);
        //may I write?
        int lockIndex =rnd_unif(100,110);
        MPI_Send(&lockIndex, 1, MPI_INT, 0, TAG_MAY_I_WRITE, MPI_COMM_WORLD);
        bool whileLoop = true;
        std::cerr << "Worker " << node->rank << " TAG_MAY_I_WRITE sent "
        << std::endl;
        while (whileLoop)
        {
            MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            std::cerr << "process Worker " << node->rank << " Received tag  " << status.MPI_TAG
            << std::endl;

            switch (status.MPI_TAG)
            {
            case TAG_DO_NOT_DARE_TO_WRITE://I am free
                MPI_Recv(0, 0, MPI_INT, 0, MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
                std::cerr << "Sleeping 0.1 seg at " << node->rank << std::endl;
                usleep(100);
                break;
            case TAG_YES_YOU_MAY_WRITE://I am free
                MPI_Recv(&lockIndex, 1, MPI_INT, 0, MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
                std::cerr << "writing at node: "  << node->rank << std::endl;
                whileLoop=false;
                break;
            default:
                std::cerr << "process WRONG TAG RECEIVED at " << node->rank << std::endl;
                break;
            }
        }
        std::cerr << "process send TAG_I_FINISH_WRITTING at " << node->rank << std::endl;
        MPI_Send(&lockIndex, 1, MPI_INT, 0, TAG_I_FINISH_WRITTING, MPI_COMM_WORLD);
    }
}
;

int main(int argc, char *argv[])
{
    ProgMpiAngularClassAverage prm(argc, argv);
    prm.read(argc, argv);
    return prm.tryRun();
}
