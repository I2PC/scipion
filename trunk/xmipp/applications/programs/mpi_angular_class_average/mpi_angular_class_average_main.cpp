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
//mpirun -np 5 xmipp_mpi_angular_class_average --nJobs 70
#include <parallel/xmipp_mpi.h>
#include <data/xmipp_funcs.h>
#include <data/xmipp_program.h>
#include <reconstruction/angular_class_average.h>
#include <data/metadata.h>


//Tags already defined in xmipp
//#define TAG_WORK                     0
//#define TAG_STOP                     1
#define TAG_I_FINISH_WRITTING        12
#define TAG_MAY_I_WRITE              13
#define TAG_YES_YOU_MAY_WRITE        14
#define TAG_DO_NOT_DARE_TO_WRITE     15
#define TAG_I_AM_FREE                16

class ProgMpiAngularClassAverage: public ProgAngularClassAverage
{
public:
    /**Number of job */
    int nJobs;

    /** status after am MPI call */
    MPI_Status status;

    // Mpi node
    MpiNode *node;

    // Metadata with the list of jobs
    MetaData mdJobList;

    ProgMpiAngularClassAverage(int argc, char **argv)
    {
        node = new MpiNode(argc, argv);
        if (!node->isMaster())
            verbose = 0;
    }

    void readParams()
    {
        ProgAngularClassAverage::readParams();

        if (node->size < 2)
            REPORT_ERROR(ERR_ARG_INCORRECT,
                         "This program cannot be executed in a single working node");
    }

    void defineParams()
    {
        ProgAngularClassAverage::defineParams();
    }

    void show()
    {
        if (!verbose)
            return;
    }

    /* Run --------------------------------------------------------------------- */
#define ArraySize 5
    void run()
    {
        //number of jobs
        //Lock structure
        bool* lockArray=new bool[nJobs+1];
        int lockIndex;
        for (lockIndex = 0; lockIndex < nJobs; ++lockIndex)
            lockArray[lockIndex]=false;

        int * Def_3Dref_2Dref_JobNo = new int[ArraySize];

        if (node->rank == 0)
        {

            rootInit();

            //for (int iCounter = 0; iCounter < nJobs; )//increase counter after I am free
            int iCounter = 0;
            int finishedNodes = 1;
            bool whileLoop = true;

            size_t id, order, count;
            int ctfGroup, ref3d;
            id = mdJobList.firstObject();
            MDIterator __iterJobs(mdJobList);
            while (whileLoop)
            {

                //wait until a worker is available
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                switch (status.MPI_TAG)
                {
                case TAG_I_AM_FREE:
                    //Some test values for defocus, 3D reference and projection direction
                    mdJobList.getValue(MDL_DEFGROUP, ctfGroup, __iterJobs.objId);
                    Def_3Dref_2Dref_JobNo[0] = ctfGroup;
                    mdJobList.getValue(MDL_REF3D, ref3d,  __iterJobs.objId);
                    Def_3Dref_2Dref_JobNo[1] = ref3d;
                    mdJobList.getValue(MDL_ORDER, order,  __iterJobs.objId);
                    Def_3Dref_2Dref_JobNo[2] = order;
                    mdJobList.getValue(MDL_COUNT, count,  __iterJobs.objId);
                    Def_3Dref_2Dref_JobNo[3] = count;
                    Def_3Dref_2Dref_JobNo[4] = iCounter;

                    //increase counter after sending work
                    ++iCounter;

                    //read worker call just to remove it from the queue
                    MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_I_AM_FREE, MPI_COMM_WORLD,
                             &status);
                    if(iCounter < nJobs)
                    {
                        //send work, first int defocus, second 3D reference, 3rd projection
                        // direction and job number
                        MPI_Send(Def_3Dref_2Dref_JobNo, ArraySize, MPI_INT, status.MPI_SOURCE,
                                 TAG_WORK, MPI_COMM_WORLD);
                        __iterJobs.moveNext();
                    }
                    else
                    {
                        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
                        finishedNodes ++;
                        if (finishedNodes >= node->size)
                            whileLoop=false;
                    }
                    break;
                case TAG_MAY_I_WRITE:
                    //where do you want to write?
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
                        MPI_Send(&lockIndex, 1, MPI_INT, status.MPI_SOURCE,
                                 TAG_YES_YOU_MAY_WRITE, MPI_COMM_WORLD);
                    }
                    break;
                case TAG_I_FINISH_WRITTING:
                    //release which lock?
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
                    break;
                }
            }
        }
        MPI_Finalize();
    }
    void process(int * Def_3Dref_2Dref_JobNo)
    {
        std::cerr
        << " DefGroup: "  << Def_3Dref_2Dref_JobNo[0]
        << " 3DRef:    "  << Def_3Dref_2Dref_JobNo[1]
        << " Order:    "  << Def_3Dref_2Dref_JobNo[2]
        << " Count:    "  << Def_3Dref_2Dref_JobNo[3]
        << " iCounter: "  << Def_3Dref_2Dref_JobNo[4]
        << " Sat node: "  << node->rank
        << std::endl;

        //may I write?
        int lockIndex =rnd_unif(1,5);
        std::cerr << "lockIndex" << lockIndex <<std::endl;
        bool whileLoop = true;
        while (whileLoop)
        {
            MPI_Send(&lockIndex, 1, MPI_INT, 0, TAG_MAY_I_WRITE, MPI_COMM_WORLD);
            MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            switch (status.MPI_TAG)
            {
            case TAG_DO_NOT_DARE_TO_WRITE://I am free
                MPI_Recv(0, 0, MPI_INT, 0, MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
                std::cerr << "Sleeping 0.1 seg at " << node->rank << std::endl;
                usleep(100000);//microsecond
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
        //REMOVE THIS DELAY!!!!!!!!!
        usleep(1000);
        MPI_Send(&lockIndex, 1, MPI_INT, 0, TAG_I_FINISH_WRITTING, MPI_COMM_WORLD);
    }

    void produceSideInfo()
    {

        ProgAngularClassAverage::produceSideInfo();
        if (node->rank == 0)
        {
            // check Wiener filter image has correct size and store the filter, we will use it later
            //This program is called once for each CTF group so there is a single wienner filter involved
            if (fn_wien != "")
            {
                // Get padding dimensions
                paddim = ROUND(pad * dim);
                Image<double> auxImg;
                //auxImg.read(fn_wien, HEADER);
                auxImg.read(fn_wien);
                Mwien = auxImg();
                if (XSIZE(Mwien) != paddim)
                {
                    std::cerr << "image size= " << dim << " padding factor= " << pad
                    << " padded image size= " << paddim
                    << " Wiener filter size= " << XSIZE(Mwien) << std::endl;
                    REPORT_ERROR(ERR_VALUE_INCORRECT,
                                 "Incompatible padding factor for this Wiener filter");
                }
            }

            // Set ring defaults
            if (Ri < 1)
                Ri = 1;
            if (Ro < 0)
                Ro = (dim / 2) - 1;

            // Set limitR
            if (do_limitR0 || do_limitRF)
            {
                MetaData tmpMT(DF);
                std::vector<double> vals;
                MDLabel codifyLabel = MDL::str2Label(col_select);
                FOR_ALL_OBJECTS_IN_METADATA(DF)
                {
                    double auxval;
                    DF.getValue(codifyLabel, auxval, __iter.objId);
                    vals.push_back(auxval);
                }
                int nn = vals.size();
                std::sort(vals.begin(), vals.end());
                if (do_limitR0)
                {
                    double val = vals[ROUND((limitR/100.) * vals.size())];
                    if (do_limit0)
                        limit0 = XMIPP_MAX(limit0, val);
                    else
                    {
                        limit0 = val;
                        do_limit0 = true;
                    }
                }
                else if (do_limitRF)
                {
                    double val = vals[ROUND(((100. - limitR)/100.) * vals.size())];
                    if (do_limitF)
                        limitF = XMIPP_MIN(limitF, val);
                    else
                    {
                        limitF = val;
                        do_limitF = true;
                    }
                }
            }


        }
    }
    void rootInit()
    {

        preprocess();
        createJobList();
    }

    void createJobList()
    {
        MetaData md;

        md.read((String)"ctfGroup[0-9][0-9][0-9][0-9][0-9][0-9]$@" + inFile);

        md.write("createJobList.xmd");
        const MDLabel myGroupByLabels[] =
            {
                MDL_DEFGROUP, MDL_REF3D, MDL_ORDER
            };
        std::vector<MDLabel> groupbyLabels(myGroupByLabels,myGroupByLabels+3);
        mdJobList.aggregateGroupBy(md, AGGR_COUNT, groupbyLabels, MDL_ORDER, MDL_COUNT);

        nJobs = mdJobList.size();
    }

}
;

int main(int argc, char *argv[])
{
    ProgMpiAngularClassAverage prm(argc, argv);
    prm.read(argc, argv);
    return prm.tryRun();
}
