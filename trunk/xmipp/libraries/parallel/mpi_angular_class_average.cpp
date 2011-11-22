/***************************************************************************
 * Authors:     AUTHOR_NAME (aerey@cnb.csic.es)
 *
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

#include "mpi_angular_class_average.h"

ProgMpiAngularClassAverage::ProgMpiAngularClassAverage(int argc, char **argv)
{
    node = new MpiNode(argc, argv);
    if (!node->isMaster())
        verbose = 0;
}

// Read arguments ==========================================================
void ProgMpiAngularClassAverage::readParams()
{
    // Read command line
    inFile = getParam("-i");
    DF.read(inFile);
    DFlib.read(getParam("--lib"));
    fn_out = getParam("-o");
    col_select = getParam("--select");
    if (checkParam("--limitR"))
    {
        limitR = getDoubleParam("--limitR");
        if (limitR < -100. || limitR > 100.)
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         "limitR should be a percentage: provide values between -100 and 100.");
        if (limitR > 0.)
            do_limitR0 = true;
        else if (limitR < 0.)
        {
            limitR *= -1.;
            do_limitRF = true;
        }
    }
    do_limit0 = checkParam("--limit0");
    if (do_limit0)
    {
        limit0 = getDoubleParam("--limit0");
    }
    do_limitF = checkParam("--limitF");
    if (do_limitF)
    {
        limitF = getDoubleParam("--limitF");
    }

    // Perform splitting of the data?
    do_split = checkParameter(argc, argv, "--split");

    // Perform Wiener filtering of average?
    fn_wien = getParameter(argc, argv, "--wien", "");
    pad = XMIPP_MAX(1.,getDoubleParam("--pad"));

    // Write selfiles
    write_selfiles = checkParam("--write_selfiles");

    // Internal re-alignment of the class averages
    Ri = getIntParam("--Ri");
    Ro = getIntParam("--Ro");
    nr_iter = getIntParam("--iter");
    max_shift = getDoubleParam("--max_shift");
    max_shift_change = getDoubleParam("--max_shift_change");
    max_psi_change = getDoubleParam("--max_psi_change");
    do_mirrors = true;

    ctfNum = getIntParam("--ctfNum");
    ref3dNum = getIntParam("--ref3dNum");

}

// Define parameters ==========================================================
void ProgMpiAngularClassAverage::defineParams()
{
    addUsageLine(
        "Make class average images and corresponding selfiles from angular_projection_matching docfiles.");

    addSeeAlsoLine("angular_project_library, angular_projection_matching");

    addParamsLine(
        "    -i <doc_file>          : Docfile with assigned angles for all experimental particles");
    addParamsLine(
        "    --lib <doc_file>       : Docfile with angles used to generate the projection matching library");
    addParamsLine(
        "    -o <root_name>         : Output rootname for class averages and selfiles");
    addParamsLine(
        "   [--split ]              : Also output averages of random halves of the data");
    addParamsLine(
        "   [--wien <img=\"\"> ]    : Apply this Wiener filter to the averages");
    addParamsLine(
        "   [--pad <factor=1.> ]    : Padding factor for Wiener correction");
    addParamsLine("   [--write_selfiles]      : Write class selfiles to disc");
    addParamsLine(
        "   [--number_3dreferences <n>] : Number of 3D references (only used with flag --postprocess_metadata).");
    addParamsLine("   [--ctfNum <n=1>]        : Ctf group number");
    addParamsLine("   [--ref3dNum <n=1>]      : 3D reference number");

    addParamsLine(
        "==+ IMAGE SELECTION BASED ON INPUT DOCFILE (select one between: limit 0, F and R ==");
    addParamsLine(
        "   [--select <col=\"maxCC\">]     : Column to use for image selection (limit0, limitF or limitR)");
    addParamsLine("   [--limit0 <l0>]         : Discard images below <l0>");
    addParamsLine("   [--limitF <lF>]         : Discard images above <lF>");
    addParamsLine(
        "   [--limitR <lR>]         : if (lR>0 && lR< 100): discard lowest  <lR> % in each class");
    addParamsLine(
        "                           : if (lR<0 && lR>-100): discard highest <lR> % in each class");

    addParamsLine("==+ REALIGNMENT OF CLASSES ==");
    addParamsLine(
        "   [--iter <nr_iter=0>]      : Number of iterations for re-alignment");
    addParamsLine(
        "   [--Ri <ri=1>]             : Inner radius to limit rotational search");
    addParamsLine(
        "   [--Ro <r0=-1>]            : Outer radius to limit rotational search");
    addParamsLine("                           : ro = -1 -> dim/2-1");
    addParamsLine(
        "   [--max_shift <ms=999.>]        : Maximum shift (larger shifts will be set to 0)");
    addParamsLine(
        "   [--max_shift_change <msc=999.>] : Discard images that change shift more in the last iteration ");
    addParamsLine(
        "   [--max_psi_change <mps=360.>]   : Discard images that change psi more in the last iteration ");

    addExampleLine(
        "Sample at default values and calculating output averages of random halves of the data",
        false);
    addExampleLine(
        "xmipp_angular_class_average -i proj_match.doc --lib ref_angles.doc -o out_dir --split");

    addKeywords("class average images");
}


/* Run --------------------------------------------------------------------- */
void ProgMpiAngularClassAverage::run()
{
    mpi_preprocess();

    //number of jobs
    //Lock structure
    bool* lockArray=new bool[nJobs+1];
    int lockIndex;
    for (lockIndex = 0; lockIndex < nJobs; ++lockIndex)
        lockArray[lockIndex]=false;

    double * Def_3Dref_2Dref_JobNo = new double[ArraySize];
    std::cerr << "*bp01*" <<std::endl;

    if (node->rank == 0)
    {


        //for (int iCounter = 0; iCounter < nJobs; )//increase counter after I am free
        int iCounter = 0;
        int finishedNodes = 1;
        bool whileLoop = true;
        double rot, tilt;

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
                Def_3Dref_2Dref_JobNo[0] = (double) ctfGroup;
                mdJobList.getValue(MDL_REF3D, ref3d,  __iterJobs.objId);
                Def_3Dref_2Dref_JobNo[1] = (double) ref3d;
                mdJobList.getValue(MDL_ORDER, order,  __iterJobs.objId);
                Def_3Dref_2Dref_JobNo[2] = (double) order;
                mdJobList.getValue(MDL_COUNT, count,  __iterJobs.objId);
                Def_3Dref_2Dref_JobNo[3] = (double) count;
                Def_3Dref_2Dref_JobNo[4] = (double) iCounter;
                mdJobList.getValue(MDL_ANGLEROT, rot,  __iterJobs.objId);
                Def_3Dref_2Dref_JobNo[5] = rot;
                mdJobList.getValue(MDL_ANGLETILT, tilt,  __iterJobs.objId);
                Def_3Dref_2Dref_JobNo[6] = tilt;

                //increase counter after sending work
                ++iCounter;

                //read worker call just to remove it from the queue
                MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_I_AM_FREE, MPI_COMM_WORLD,
                         &status);
                if(iCounter < nJobs)
                {
                    //send work, first int defocus, second 3D reference, 3rd projection
                    // direction and job number
                    MPI_Send(Def_3Dref_2Dref_JobNo, ArraySize, MPI_DOUBLE, status.MPI_SOURCE,
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
                MPI_Recv(Def_3Dref_2Dref_JobNo, ArraySize, MPI_DOUBLE, 0, TAG_WORK,
                         MPI_COMM_WORLD, &status);
                mpi_process(Def_3Dref_2Dref_JobNo);
                break;
            default:
                break;
            }
        }
    }
    mpi_postprocess();

    MPI_Finalize();
}

void ProgMpiAngularClassAverage::mpi_process(double * Def_3Dref_2Dref_JobNo)
{
    std::cerr
    << " DefGroup: "  << ROUND(Def_3Dref_2Dref_JobNo[0])
    << " 3DRef:    "  << ROUND(Def_3Dref_2Dref_JobNo[1])
    << " Order:    "  << ROUND(Def_3Dref_2Dref_JobNo[2])
    << " Count:    "  << ROUND(Def_3Dref_2Dref_JobNo[3])
    << " iCounter: "  << ROUND(Def_3Dref_2Dref_JobNo[4])
    << " rot:      "  << Def_3Dref_2Dref_JobNo[5]
    << " tilt:     "  << Def_3Dref_2Dref_JobNo[6]
    << " Sat node: "  << node->rank
    << std::endl;

    //processOneClass(dirno, output_values);

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
    //usleep(1000);
    MPI_Send(&lockIndex, 1, MPI_INT, 0, TAG_I_FINISH_WRITTING, MPI_COMM_WORLD);
}

void ProgMpiAngularClassAverage::mpi_produceSideInfo()
{
    // Set up output rootnames
    if (do_split)
    {
        fn_out1 = fn_out + "_split_1";
        fn_out2 = fn_out + "_split_2";
    }

    //init with 0 by default through memset
    Iempty().resizeNoCopy(Xdim, Ydim);
    Iempty().setXmippOrigin();

    // Randomization
    if (do_split)
        randomize_random_generator();

    if (Ri<1)
        Ri=1;
    if (Ro<0)
        Ro=(Xdim/2)-1;
    // Set up FFTW transformers
    //Is this needed if  no alignment is required?
    MultidimArray<double> Maux;
    Polar<double> P;
    Polar<std::complex<double> > fP;

    produceSplineCoefficients(BSPLINE3, Maux, Iempty());
    P.getPolarFromCartesianBSpline(Maux, Ri, Ro);
    P.calculateFftwPlans(global_plans);
    fourierTransformRings(P, fP, global_plans, false);
    corr.resize(P.getSampleNoOuterRing());
    global_transformer.setReal(corr);
    global_transformer.FourierTransform();


    if (node->rank == 0)
    {
        // check Wiener filter image has correct size and store the filter, we will use it later
        //This program is called once for each CTF group so there is a single wienner filter involved
        if (fn_wien != "")
        {
            // Get padding dimensions
            paddim = ROUND(pad * Xdim);
            Image<double> auxImg;
            //auxImg.read(fn_wien, HEADER);
            auxImg.read(fn_wien);
            Mwien = auxImg();
            if (XSIZE(Mwien) != paddim)
            {
                std::cerr << "image size= " << Xdim << " padding factor= " << pad
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
            Ro = (Xdim / 2) - 1;

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
void ProgMpiAngularClassAverage::mpi_preprocess()
{
    if (node->rank==0)
    {
        createJobList();
        deleteOutputFiles();
    }
    MPI_Bcast(&Xdim,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Ydim,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Zdim,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Ndim,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nJobs,1,MPI_INT,0,MPI_COMM_WORLD);
    mpi_produceSideInfo();
}

void ProgMpiAngularClassAverage::deleteOutputFiles()
{
    //alloc space for output files
    FileName fn_tmp;

    getImageSize(DF, Xdim, Ydim, Zdim, Ndim);

    // Output stack size (number of valid projection directions)
    MDObject mdValueOut2(MDL_ORDER);
    mdJobList.aggregateSingleSizeT(mdValueOut2, AGGR_MAX ,MDL_ORDER);
    mdValueOut2.getValue(Ndim);


    std::cerr << "X|Y|Z|N: " << Xdim << " " << Ydim << " " << Zdim << " " << Ndim << std::endl;
    for (int i = 1; i <= number_3dref; i++)
    {
        formatStringFast(fn_tmp, "_refGroup%06lu", i);

        unlink((fn_out + fn_tmp + ".xmd").c_str());
        unlink((fn_out + fn_tmp + ".stk").c_str());
        createEmptyFile(fn_out + fn_tmp + ".stk", Xdim, Ydim, Zdim, Ndim, true,
                        WRITE_OVERWRITE);
        if (do_split)
        {
            unlink((fn_out1 + fn_tmp + ".xmd").c_str());
            unlink((fn_out1 + fn_tmp + ".stk").c_str());
            createEmptyFile(fn_out1 + fn_tmp + ".stk", Xdim, Ydim, Zdim, Ndim,
                            true, WRITE_OVERWRITE);
            unlink((fn_out2 + fn_tmp + ".xmd").c_str());
            unlink((fn_out2 + fn_tmp + ".stk").c_str());
            createEmptyFile(fn_out2 + fn_tmp + ".stk", Xdim, Ydim, Zdim, Ndim,
                            true, WRITE_OVERWRITE);
        }
        unlink((fn_out + fn_tmp + "_discarded.xmd").c_str());

    }

}
void ProgMpiAngularClassAverage::mpi_postprocess()
{}

void ProgMpiAngularClassAverage::createJobList()
{
    MetaData md;

    md.read((String)"ctfGroup[0-9][0-9][0-9][0-9][0-9][0-9]$@" + inFile);

    const MDLabel myGroupByLabels[] =
        {
            MDL_DEFGROUP, MDL_REF3D, MDL_ORDER, MDL_ANGLEROT, MDL_ANGLETILT
        };
    std::vector<MDLabel> groupbyLabels(myGroupByLabels,myGroupByLabels+3);
    mdJobList.aggregateGroupBy(md, AGGR_COUNT, groupbyLabels, MDL_ORDER, MDL_COUNT);

    //mdJobList.write("createJobList.xmd");

    nJobs = mdJobList.size();
}

