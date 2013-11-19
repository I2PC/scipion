/***************************************************************************
 * Authors:     AUTHOR_NAME (aerey@cnb.csic.es)
 *        (jvargas@cnb.csic.es)
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


MpiProgAngularClassAverage::MpiProgAngularClassAverage()
{}

MpiProgAngularClassAverage::MpiProgAngularClassAverage(int argc, char **argv)
{
    this->read(argc, argv);
}

// Read arguments ==========================================================
void MpiProgAngularClassAverage::readParams()
{
    do_limitR0per = do_limitR0class = do_limitRFclass = do_limit0 = do_limitF = false;

    // Read command line
    fn_out = getParam("-o");
    fn_ref  = getParam("--lib");

    col_select = getParam("--select");

    do_limit0 = checkParam("--limit0");
    if (do_limit0)
        limit0 = getDoubleParam("--limit0");
    do_limitF = checkParam("--limitF");
    if (do_limitF)
        limitF = getDoubleParam("--limitF");

    inFile = getParam("-i");

    do_pcaSorting = false;
    do_pcaSorting= checkParam("--pcaSorting");

    do_limitRFclass = false;
    do_limitR0class = false;
    do_limitR0per = false;
    do_limitRFper = false;
    if (checkParam("--limitRclass"))
    {
        limitRclass = getDoubleParam("--limitRclass")/100.;
        if (limitRclass < -1. || limitRclass > 1.)
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         "limitRclass should be a percentage: provide values between -100 and 100.");
        if (limitRclass > 0.)
            do_limitR0class = true;
        else if (limitRclass < 0.)
        {
            limitRclass *= -1.;
            do_limitRFclass = true;
        }
    }


    if (checkParam("--limitRper"))
    {
        limitRper = getDoubleParam("--limitRper");
        if (limitRper < -100. || limitRper > 100.)
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         "limitRper should be a percentage: provide values between -100 and 100.");
        if (limitRper > 0.)
            do_limitR0per = true;
        else if (limitRper < 0.)
        {
            limitRper *= -1.;
            do_limitRFper = true;
        }
    }

    if ((do_limitR0per && (do_limitR0class ||  do_limitRFclass)) ||
        (do_limitR0per && (do_limit0 || do_limitF )) ||
        ((do_limitR0class ||  do_limitRFclass) && (do_limit0 || do_limitF )))
        REPORT_ERROR(ERR_VALUE_INCORRECT, "You can not use different kind of limits at the same time.");

    // Perform splitting of the data?
    do_split = checkParam("--split");

    // Perform Wiener filtering of average?
    fn_wien = getParam("--wien");
    pad = XMIPP_MAX(1.,getDoubleParam("--pad"));

    // Internal re-alignment of the class averages
    Ri = getIntParam("--Ri");
    Ro = getIntParam("--Ro");
    nr_iter = getIntParam("--iter");
    do_mirrors = true;

    do_save_images_assigned_to_classes = checkParam("--save_images_assigned_to_classes");
    mpi_job_size = getIntParam("--mpi_job_size");
}

// Define parameters ==========================================================
void MpiProgAngularClassAverage::defineParams()
{
    addUsageLine("Make class average images and corresponding selfiles from angular_projection_matching docfiles.");
    addSeeAlsoLine("angular_project_library, angular_projection_matching");
    addKeywords("class average images");

    addParamsLine("    -i <doc_file>          : Docfile with assigned angles for all experimental particles");
    addParamsLine("    --lib <doc_file>       : Docfile with angles used to generate the projection matching library");
    addParamsLine("    -o <root_name>         : Output rootname for class averages and selfiles");
    addParamsLine("   [--split ]              : Also output averages of random halves of the data");
    addParamsLine("   [--wien <img=\"\"> ]    : Apply this Wiener filter to the averages");
    addParamsLine("   [--pad <factor=1.> ]    : Padding factor for Wiener correction");
    addParamsLine("   [--save_images_assigned_to_classes]    : Save images assigned te each class in output metadatas");
    addParamsLine("alias --siatc;");
    addParamsLine("==+ IMAGE SELECTION BASED ON INPUT DOCFILE (select one between: limit 0, F and R ==");
    addParamsLine("   [--select <col=\"maxCC\">]     : Column to use for image selection (limit0, limitF or limitR)");
    addParamsLine("   [--limit0 <l0>]         : Discard images below <l0>");
    addParamsLine("   [--limitF <lF>]         : Discard images above <lF>");
    addParamsLine("   [--limitRclass <lRc>]         : if (lRc>0 && lRc< 100): discard lowest  <lRc> % in each class");
    addParamsLine("                           : if (lRc<0 && lR>-100): discard highest <lRc> % in each class");
    addParamsLine("   [--limitRper <lRp>]         : if (lRp>0 && lRp< 100): discard lowest  <lRa> %");
    addParamsLine("                           : if (lRp<0 && lRp>-100): discard highest <lRa> %");
    addParamsLine("   [--pcaSorting ]         : Perform PCA sorting to obtain the average classes");

    addParamsLine("==+ REALIGNMENT OF CLASSES ==");
    addParamsLine("   [--iter <nr_iter=0>]      : Number of iterations for re-alignment");
    addParamsLine("   [--Ri <ri=1>]             : Inner radius to limit rotational search");
    addParamsLine("   [--Ro <r0=-1>]            : Outer radius to limit rotational search");
    addParamsLine("                           : ro = -1 -> dim/2-1");
    addParamsLine("  [--mpi_job_size <size=10>]   : Number of images sent to a cpu in a single job ");
    addParamsLine("                                : 10 may be a good value");

    addExampleLine("Sample at default values and calculating output averages of random halves of the data",false);
    addExampleLine("xmipp_angular_class_average -i proj_match.doc --lib ref_angles.doc -o out_dir --split");

}


/* Run --------------------------------------------------------------------- */
void MpiProgAngularClassAverage::run()
{
    mpi_preprocess();

    int lockIndex;
    size_t order_index;
    int ref3d_index;

    double lockWeightIndexes[lockWeightIndexesSize];

    double * jobListRows = new double[ArraySize * mpi_job_size + 1];

    if (node->rank == 0)
    {
        //for (int iCounter = 0; iCounter < nJobs; )//increase counter after I am free
        size_t jobId = 0, size;
        size_t finishedNodes = 1;
        bool whileLoop = true;

        size_t order, count;
        int ctfGroup, ref3d, ref2d;
        MDIterator __iterJobs(mdJobList);
        while (whileLoop)
        {
            //wait until a worker is available
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            switch (status.MPI_TAG)
            {
            case TAG_I_AM_FREE:

                size = 0;
                for (size_t i=0;i<mpi_job_size && jobId < numberOfJobs;i++,size++,jobId++)
                {
                    //Some test values for defocus, 3D reference and projection direction
                    mdJobList.getValue(MDL_REF3D, ref3d,  __iterJobs.objId);
                    jobListRows[index_3DRef + ArraySize * i + 1] = (double) ref3d;
                    mdJobList.getValue(MDL_DEFGROUP, ctfGroup, __iterJobs.objId);
                    jobListRows[index_DefGroup + ArraySize * i + 1] = (double) ctfGroup;
                    mdJobList.getValue(MDL_ORDER, order,  __iterJobs.objId);
                    jobListRows[index_Order + ArraySize * i + 1] = (double) order;
                    mdJobList.getValue(MDL_COUNT, count,  __iterJobs.objId);
                    jobListRows[index_Count + ArraySize * i + 1] = (double) count;
                    mdJobList.getValue(MDL_REF, ref2d,  __iterJobs.objId);
                    jobListRows[index_2DRef + ArraySize * i + 1] = (double) ref2d;
                    jobListRows[index_jobId + ArraySize * i + 1] = (double) __iterJobs.objId;
                    mdJobList.getValue(MDL_ANGLE_ROT, jobListRows[index_Rot + ArraySize * i + 1],  __iterJobs.objId);
                    mdJobList.getValue(MDL_ANGLE_TILT, jobListRows[index_Tilt + ArraySize * i + 1],  __iterJobs.objId);
                    __iterJobs.moveNext();
                }
                jobListRows[0]=size;

                //read worker call just to remove it from the queue
                MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_I_AM_FREE, MPI_COMM_WORLD,
                         &status);
                if(size > 0)
                {
                    //send work, first int defocus, second 3D reference, 3rd projection
                    // direction and job number
                    //#define DEBUG_MPI
#ifdef DEBUG_MPI
                    usleep(10000);
                    std::cerr << "Sending job to worker " << status.MPI_SOURCE <<std::endl;
#endif

                    MPI_Send(jobListRows, ArraySize * mpi_job_size + 1, MPI_DOUBLE, status.MPI_SOURCE,
                             TAG_WORK, MPI_COMM_WORLD);
                    //__iterJobs.moveNext();
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
                MPI_Recv(&lockIndex, 1, MPI_INT, MPI_ANY_SOURCE, TAG_MAY_I_WRITE, MPI_COMM_WORLD, &status);

                mdJobList.getValue(MDL_ORDER, order_index, lockIndex);
                mdJobList.getValue(MDL_REF3D, ref3d_index, lockIndex);

                //#define DEBUG_MPI
#ifdef DEBUG_MPI

                std::cerr << "Blocking. lockIndex: " << lockIndex << " | status.MPI_SOURCE: " << status.MPI_SOURCE
                << " | lockArray[" << order_index<< "," <<ref3d_index<< "]: " << dAij(lockArray,order_index,ref3d_index) << std::endl;
#endif

                if (dAij(lockArray,order_index,ref3d_index))
                {//Locked
                    MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE,
                             TAG_DO_NOT_DARE_TO_WRITE, MPI_COMM_WORLD);
                }
                else
                {//Unlocked
                    //JV : mirar esto imprimiendo
                    dAij(lockArray,order_index,ref3d_index)=true;
                    // Send the old weight
                    lockWeightIndexes[index_lockIndex] = lockIndex;
                    lockWeightIndexes[index_weight] = dAkij(weightArray,0,order_index,ref3d_index);
                    lockWeightIndexes[index_weights1] = dAkij(weightArrays1,0,order_index,ref3d_index);
                    lockWeightIndexes[index_weights2] = dAkij(weightArrays2,0,order_index,ref3d_index);
                    lockWeightIndexes[index_ref3d] = ref3d_index;
                    //lockWeightIndexes[index_weight] = dAij(weightArray,order_index,ref3d_index);
                    //lockWeightIndexes[index_weights1] = dAij(weightArrays1,order_index,ref3d_index);
                    //lockWeightIndexes[index_weights2] = dAij(weightArrays2,order_index,ref3d_index);
                    //lockWeightIndexes[index_ref3d] = ref3d_index;
#ifdef DEBUG_MPI

                    std::cerr << "[" << node->rank << "] TAG_YES_YOU_MAY_WRITE.lockIndex: " << lockIndex <<std::endl;
                    std::cerr <<  "lockWeightIndexes[index_weight]" << lockWeightIndexes[index_weight]<<std::endl;
#endif
#undef  DEBUG_MPI

                    MPI_Send(lockWeightIndexes, lockWeightIndexesSize, MPI_DOUBLE, status.MPI_SOURCE,
                             TAG_YES_YOU_MAY_WRITE, MPI_COMM_WORLD);
                }
                break;
            case TAG_I_FINISH_WRITTING:
                //release which lock?
                MPI_Recv(lockWeightIndexes, lockWeightIndexesSize, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_I_FINISH_WRITTING,
                         MPI_COMM_WORLD, &status);

                lockIndex = (int)lockWeightIndexes[index_lockIndex];
                //#define DEBUG_MPI
#ifdef DEBUG_MPI

                std::cerr << "Unblocking. lockIndex: " << lockIndex << " | status.MPI_SOURCE: " << status.MPI_SOURCE << std::endl;
#endif

                mdJobList.getValue(MDL_ORDER, order_index, lockIndex);
                mdJobList.getValue(MDL_REF3D, ref3d_index, lockIndex);
                dAij(lockArray,order_index,ref3d_index)=false;
                //dAij(weightArray,order_index,ref3d_index) += lockWeightIndexes[index_weight];
                //dAij(weightArrays1,order_index,ref3d_index) += lockWeightIndexes[index_weights1];
                //dAij(weightArrays2,order_index,ref3d_index) += lockWeightIndexes[index_weights2];

                dAkij(weightArray,0,order_index,ref3d_index) += lockWeightIndexes[index_weight];
                dAkij(weightArrays1,0,order_index,ref3d_index) += lockWeightIndexes[index_weights1];
                dAkij(weightArrays2,0,order_index,ref3d_index) += lockWeightIndexes[index_weights2];
                break;
            default:
                std::cerr << "WRONG TAG RECEIVED" << std::endl;
                break;
            }
        }
    }
    else
    {
        bool whileLoop = true;
        while (whileLoop)
        {
            //#define DEBUG_MPI
#ifdef DEBUG_MPI
            std::cerr << "[" << node->rank << "] Asking for a job " <<std::endl;
#endif
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
                MPI_Recv(jobListRows, ArraySize * mpi_job_size + 1, MPI_DOUBLE, 0, TAG_WORK,
                         MPI_COMM_WORLD, &status);
                mpi_process_loop(jobListRows);
                break;
            default:
                break;
            }
        }
    }
    FileName gatherFile;
    formatStringFast( gatherFile, "%s_GatherMetadata.xmd", fn_out.c_str());
    node->gatherMetadatas(DFscore, gatherFile);

    if (node->rank == 0)
    {
        FileName fn = inFile.withoutExtension()+"_weight."+inFile.getExtension();
        fn = "all_exp_images@"+fn.removeBlockName();
        DFscore.write(fn);
        mpi_postprocess();
    }

    //    MPI_Finalize();
}


void MpiProgAngularClassAverage::mpi_process_loop(double * Def_3Dref_2Dref_JobNo)
{
    int size = ROUND(Def_3Dref_2Dref_JobNo[0]);

    for(int i=0;i<size;i++)
    {
        mpi_process(Def_3Dref_2Dref_JobNo+i*ArraySize+1);
    }

}

void MpiProgAngularClassAverage::mpi_process(double * Def_3Dref_2Dref_JobNo)
{
#define DEBUG
#ifdef DEBUG
    std::cerr<<"["<<node->rank<<"]"
    << " 3DRef:    "  << ROUND(Def_3Dref_2Dref_JobNo[index_3DRef])
    << " DefGroup: "  << ROUND(Def_3Dref_2Dref_JobNo[index_DefGroup])
    << " 2DRef:    "  << ROUND(Def_3Dref_2Dref_JobNo[index_2DRef])
    << " Order:    "  << ROUND(Def_3Dref_2Dref_JobNo[index_Order])
    << " Count:    "  << ROUND(Def_3Dref_2Dref_JobNo[index_Count])
    << " jobId:    "  << ROUND(Def_3Dref_2Dref_JobNo[index_jobId])
    << " rot:      "  << Def_3Dref_2Dref_JobNo[index_Rot]
    << " tilt:     "  << Def_3Dref_2Dref_JobNo[index_Tilt]
    << " Sat node: "  << node->rank
    << std::endl;
#endif
#undef DEBUG

    Image<double> img, img_ref, avg, avg1, avg2;
    FileName fn_img, fn_tmp;
    MetaData SFclass, SFclass1, SFclass2;
    MetaData SFclassDiscarded;
    double psi, xshift, yshift, w, w1, w2, scale;
    bool mirror;
    int ref_number, this_image, ref3d, defGroup;
    static int defGroup_last = 0;
    int isplit, lockIndex;
    MetaData _DF;
    size_t id;
    size_t order_number;

    w = 0.;
    w1 = 0.;
    w2 = 0.;
    this_image = 0;

    order_number = ROUND(Def_3Dref_2Dref_JobNo[index_Order]);
    ref_number   = ROUND(Def_3Dref_2Dref_JobNo[index_2DRef]);
    defGroup     = ROUND(Def_3Dref_2Dref_JobNo[index_DefGroup]);
    ref3d        = ROUND(Def_3Dref_2Dref_JobNo[index_3DRef]);
    std::cerr << "DEBUG_ROB: order_number: " << order_number << std::endl;
    lockIndex    = ROUND(Def_3Dref_2Dref_JobNo[index_jobId]);

    if (fn_wien != "" && defGroup_last != defGroup)
    {
        // Read wiener filter
        FileName fn_wfilter;
        Image<double> auxImg;

        fn_wfilter.compose(defGroup,fn_wien);
        auxImg.read(fn_wfilter);
        Mwien = auxImg();
        defGroup_last = defGroup;
    }

    //std::cerr << "DEBUG_JM: BEFORE MDValueEQ" <<std::endl;
    MDValueEQ eq1(MDL_REF3D, ref3d);
    MDValueEQ eq2(MDL_DEFGROUP, defGroup);
    MDValueEQ eq3(MDL_REF, ref_number);
    MDMultiQuery multi;
    multi.addAndQuery(eq1);
    multi.addAndQuery(eq2);
    multi.addAndQuery(eq3);

    _DF.importObjects(DF, multi);

    //std::cerr << "DEBUG_JM: AFTER MDValueEQ" <<std::endl;

    if (_DF.size() == 0)
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,
                     "Program should never execute this line, something went wrong");


    MetaData _DF_ref((fn_ref.removeAllExtensions()+".doc"));
    MetaData _DF_temp;
    MDValueEQ eq4(MDL_REF, ref_number);
    _DF_temp.importObjects(_DF_ref, eq4);
    int noRef = 0;
    if (_DF_temp.size() > 1)
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,
                     "Program should never execute this line, something went wrong");
    if (_DF_temp.size() == 0)
        noRef =1;

    if (noRef!=1)
    {
        FOR_ALL_OBJECTS_IN_METADATA(_DF_temp)
        {
            _DF_temp.getValue(MDL_IMAGE, fn_img,__iter.objId);
            img_ref.read(fn_img);
        }
    }

    Matrix2D<double> A(3, 3);
    std::vector<int> exp_number, exp_split;
    std::vector<Image<double> > exp_imgs;

    Iempty.setEulerAngles(Def_3Dref_2Dref_JobNo[index_Rot], Def_3Dref_2Dref_JobNo[index_Tilt], 0.);
    Iempty.setShifts(0., 0.);
    Iempty.setFlip(0.);
    avg = Iempty;
    avg1 = Iempty;
    avg2 = Iempty;

    PCAMahalanobisAnalyzer pcaAnalyzerSplit1,pcaAnalyzerSplit2,pcaAnalyzer;
    pcaAnalyzer.clear();
    pcaAnalyzerSplit1.clear();
    pcaAnalyzerSplit2.clear();

    // Loop over all images in the input docfile
    FOR_ALL_OBJECTS_IN_METADATA(_DF)
    {
        _DF.getValue(MDL_IMAGE, fn_img, __iter.objId);
        this_image++;
        _DF.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
        _DF.getValue(MDL_SHIFT_X, xshift, __iter.objId);
        _DF.getValue(MDL_SHIFT_Y, yshift, __iter.objId);
        if (do_mirrors)
            _DF.getValue(MDL_FLIP, mirror, __iter.objId);
        _DF.getValue(MDL_SCALE, scale, __iter.objId);

        img.read(fn_img);
        img().setXmippOrigin();
        img.setEulerAngles(0., 0., psi);
        img.setShifts(-xshift, -yshift);

        if (do_mirrors)
            img.setFlip(mirror);
        img.setScale(scale);

        if (do_split)
            isplit = ROUND(rnd_unif());
        else
            isplit = 0;
        // For re-alignment of class: store all images in memory
        if (nr_iter > 0)
        {
            exp_imgs.push_back(img);
            exp_number.push_back(this_image);
        }

        if ( (nr_iter > 0) || do_pcaSorting)
        {
            exp_split.push_back(isplit);
        }

        // Apply in-plane transformation
        img.getTransformationMatrix(A);
        if (!A.isIdentity())
            selfApplyGeometry(BSPLINE3, img(), A, IS_INV, DONT_WRAP);

        MultidimArray<float> auxImg(img().nzyxdim);
        const MultidimArray<double> &mImg=img();
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mImg)
        DIRECT_MULTIDIM_ELEM(auxImg,n)=(float)DIRECT_MULTIDIM_ELEM(mImg,n);

        if (do_pcaSorting)
        {
            pcaAnalyzer.addVector(auxImg);
        }

        // Add to average
        if (isplit == 0)
        {
            avg1() += img();
            w1 += 1.;
            id = SFclass1.addObject();
            SFclass1.setValue(MDL_IMAGE, fn_img, id);
            SFclass1.setValue(MDL_ANGLE_ROT, Def_3Dref_2Dref_JobNo[index_Rot], id);
            SFclass1.setValue(MDL_ANGLE_TILT, Def_3Dref_2Dref_JobNo[index_Tilt], id);
            SFclass1.setValue(MDL_REF, ref_number, id);
            SFclass1.setValue(MDL_REF3D, ref3d, id);
            SFclass1.setValue(MDL_DEFGROUP, defGroup, id);
            SFclass1.setValue(MDL_ORDER, order_number, id);
            if (do_pcaSorting)
            {
                pcaAnalyzerSplit1.addVector(auxImg);
            }
        }
        else
        {
            avg2() += img();
            w2 += 1.;
            id = SFclass2.addObject();
            SFclass2.setValue(MDL_IMAGE, fn_img, id);
            SFclass2.setValue(MDL_ANGLE_ROT, Def_3Dref_2Dref_JobNo[index_Rot], id);
            SFclass2.setValue(MDL_ANGLE_TILT, Def_3Dref_2Dref_JobNo[index_Tilt], id);
            SFclass2.setValue(MDL_REF, ref_number, id);
            SFclass2.setValue(MDL_REF3D, ref3d, id);
            SFclass2.setValue(MDL_DEFGROUP, defGroup, id);
            SFclass2.setValue(MDL_ORDER, order_number, id);
            if (do_pcaSorting)
            {
                pcaAnalyzerSplit2.addVector(auxImg);
            }
        }

        //#define DEBUG
#ifdef DEBUG
        //WRITE IMAGES TO AVERAGE
        FileName fn_tmp1;
        int static static_i=0;
        static_i++;
        formatStringFast(fn_tmp1, "test_img_%06d.spi", static_i);
        img.write(fn_tmp1);
        formatStringFast(fn_tmp1, "test_avg_%06d.spi", static_i);
        avg1.write(fn_tmp1);
        std::cout << fn_img
        << psi
        << A <<std::endl;
        if (static_i> 25)
        {
            std::cerr << "static_i:" << static_i << std::endl;
            _DF.write("DF.xmd");
            exit(1);
        }
#endif
#undef DEBUG

    }

    //this_image = 0;
    if (do_pcaSorting)
    {
        avg().initZeros();
        avg1().initZeros();
        avg2().initZeros();
        double max_cc = 0;

        w1 = 0;
        w2 = 0;

        if (do_split == 0)
        {
            size_t index = 0;
            if (_DF.size() <= 5 )
            {
                //pcaAnalyzer.evaluateZScore(1,20);
                FOR_ALL_OBJECTS_IN_METADATA(_DF)
                {
                    _DF.getValue(MDL_MAXCC, max_cc, __iter.objId);
                    _DF.setValue(MDL_WEIGHT,max_cc, __iter.objId);
                    w1 += max_cc;
                    avg1() += img()*(max_cc);
                    index++;
                }
            }
            else
            {
                pcaAnalyzer.evaluateZScore(1,20);
                //We must tour the metadata again:
                FOR_ALL_OBJECTS_IN_METADATA(_DF)
                {
                    _DF.getValue(MDL_IMAGE, fn_img, __iter.objId);
                    _DF.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
                    _DF.getValue(MDL_SHIFT_X, xshift, __iter.objId);
                    _DF.getValue(MDL_SHIFT_Y, yshift, __iter.objId);
                    if (do_mirrors)
                        _DF.getValue(MDL_FLIP, mirror, __iter.objId);
                    _DF.getValue(MDL_SCALE, scale, __iter.objId);

                    _DF.getValue(MDL_MAXCC, max_cc, __iter.objId);

                    img.read(fn_img);
                    img().setXmippOrigin();
                    img.setEulerAngles(0., 0., psi);
                    img.setShifts(-xshift, -yshift);

                    if (do_mirrors)
                        img.setFlip(mirror);
                    img.setScale(scale);

                    // Apply in-plane transformation
                    img.getTransformationMatrix(A);
                    if (!A.isIdentity())
                        selfApplyGeometry(BSPLINE3, img(), A, IS_INV, DONT_WRAP);

                    double score = pcaAnalyzer.getZscore(index);

                    if (score < 5)
                    {
                        avg1() += img()*(5-score)*max_cc;
                        w1 +=(5-score)*max_cc;
                        _DF.setValue(MDL_WEIGHT,(5-score)*max_cc, __iter.objId);
                    }
                    else
                        _DF.setValue(MDL_WEIGHT,0, __iter.objId);

                    index++;
                }
            }
        }
        else
        {
            pcaAnalyzerSplit1.evaluateZScore(1,20);
            pcaAnalyzerSplit2.evaluateZScore(1,20);
            size_t index = 0;
            size_t index1 = 0;
            size_t index2 = 0;

            FOR_ALL_OBJECTS_IN_METADATA(_DF)
            {
                isplit = exp_split[index];
                _DF.getValue(MDL_IMAGE, fn_img, __iter.objId);
                //this_image++;
                _DF.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
                _DF.getValue(MDL_SHIFT_X, xshift, __iter.objId);
                _DF.getValue(MDL_SHIFT_Y, yshift, __iter.objId);
                if (do_mirrors)
                    _DF.getValue(MDL_FLIP, mirror, __iter.objId);
                _DF.getValue(MDL_SCALE, scale, __iter.objId);

                _DF.getValue(MDL_MAXCC, max_cc, __iter.objId);

                img.read(fn_img);
                img().setXmippOrigin();
                img.setEulerAngles(0., 0., psi);
                img.setShifts(-xshift, -yshift);

                if (do_mirrors)
                    img.setFlip(mirror);
                img.setScale(scale);

                // Apply in-plane transformation
                img.getTransformationMatrix(A);
                if (!A.isIdentity())
                    selfApplyGeometry(BSPLINE3, img(), A, IS_INV, DONT_WRAP);

                if (isplit == 0)
                {
                    if (pcaAnalyzerSplit1.v.size() > 6 )
                    {
                        double score = pcaAnalyzerSplit1.getZscore(index1);
                        if (score < 5)
                        {
                            avg1() += img()*(5-score)*max_cc;
                            w1 +=(5-score)*max_cc;
                            index1++;
                            _DF.setValue(MDL_WEIGHT,(5-score)*max_cc, __iter.objId);
                        }
                        else
                        {
                            index1++;
                            _DF.setValue(MDL_WEIGHT,0, __iter.objId);
                        }
                    }
                    else
                    {
                        _DF.setValue(MDL_WEIGHT,max_cc, __iter.objId);
                        w1 += max_cc;
                        avg1() += img()*(w1);
                        index++;
                    }

                }
                else
                {
                    if (pcaAnalyzerSplit2.v.size() > 6 )
                    {
                        double score = pcaAnalyzerSplit2.getZscore(index2);
                        if (score < 5)
                        {
                            avg2() += img()*(5-score)*max_cc;
                            w2 +=(5-score)*max_cc;
                            index2++;
                            _DF.setValue(MDL_WEIGHT,(5-score)*max_cc, __iter.objId);
                        }
                        else
                        {
                            index2++;
                            _DF.setValue(MDL_WEIGHT,0, __iter.objId);
                        }
                    }
                    else
                    {
                        _DF.setValue(MDL_WEIGHT,max_cc, __iter.objId);
                        w2 += max_cc;
                        avg2() += img()*(w2);
                        index++;
                    }
                }
                index++;
            }
        }
    }

    // Re-alignment of the class
    if (nr_iter > 0)
    {
        SFclass = SFclass1;
        SFclass.unionAll(SFclass2);
        avg() = avg1() + avg2();
        w = w1 + w2;
        avg.setWeight(w);

        // Reserve memory for output from class realignment
        int reserve = DF.size();
        double *my_output=new double[AVG_OUPUT_SIZE * reserve + 1];

        reAlignClass(avg1, avg2, SFclass1, SFclass2, exp_imgs, exp_split,
                     exp_number, order_number, my_output);
        w1 = avg1.weight();
        w2 = avg2.weight();
        delete []my_output;
    }

    // Apply Wiener filters
    if (fn_wien != "")
    {
        if (w1 > 0)
            applyWienerFilter(avg1());
        if (w2 > 0)
            applyWienerFilter(avg2());
    }

    // Output total and split averages and selfiles to disc
    SFclass = SFclass1;
    SFclass.unionAll(SFclass2);

    avg() = avg1() + avg2();
    w = w1+w2;

    double t1;
    double t2;
    double t;
    if (noRef!=1)
    {
        img_ref().setXmippOrigin();
        if (w1 != 0)
        {
            t1 = correlationIndex(img_ref(),avg1()/w1);
            avg1()*=(t1/w1);
            w1 = t1;
        }
        else
        {
            t1 = 0.;
            w1 = 0.;
            avg1().initZeros();
        }

        if (w2 != 0)
        {
            t2 = correlationIndex(img_ref(),avg2()/w2);
            avg2()*=(t2/w2);
            w2 = t2;
        }
        else
        {
            t2 = 0.;
            w2 = 0.;
            avg2().initZeros();

        }

        if (w != 0)
        {
            t = correlationIndex(img_ref(),avg()/w);
            avg()*=(t/w);
            w=t;
        }
        else
        {
            t = 0.;
            w = 0.;
            avg().initZeros();
        }
    }

    avg.setWeight(w);
    avg1.setWeight(w1);
    avg2.setWeight(w2);

    DFscore.unionAll(_DF);
    mpi_writeController(order_number, avg, avg1, avg2, SFclass, SFclass1, SFclass2,
                        SFclassDiscarded,_DF, w1, w2, w, lockIndex);

}



void MpiProgAngularClassAverage::mpi_write(
    size_t dirno,
    int ref3dIndex,
    Image<double> avg,
    Image<double> avg1,
    Image<double> avg2,
    MetaData SFclass,
    MetaData SFclass1,
    MetaData SFclass2,
    MetaData SFclassDiscarded,
    double w1,
    double w2,
    double old_w,
    double old_w1,
    double old_w2)
{
    FileName fileNameXmd, fileNameStk;

    formatStringFast(fileNameStk, "%s_Ref3D_%03lu.stk", fn_out.c_str(), ref3dIndex);
    mpi_writeFile(avg, dirno, fileNameStk, old_w);

    if (do_split)
    {
        if (w1 > 0)
        {
            formatStringFast(fileNameStk, "%s_Ref3D_%03lu.stk",
                             fn_out1.c_str(), ref3dIndex);
            mpi_writeFile(avg1, dirno, fileNameStk, old_w1);
        }
        if (w2 > 0)
        {
            formatStringFast(fileNameStk, "%s_Ref3D_%03lu.stk",
                             fn_out2.c_str(), ref3dIndex);
            mpi_writeFile(avg2, dirno, fileNameStk, old_w2);
        }
    }
}


void MpiProgAngularClassAverage::mpi_writeController(
    size_t dirno,
    Image<double> avg,
    Image<double> avg1,
    Image<double> avg2,
    MetaData SFclass,
    MetaData SFclass1,
    MetaData SFclass2,
    MetaData SFclassDiscarded,
    MetaData _DF,
    double w1,
    double w2,
    double w,
    int lockIndex)
{

    double weight_old, weights1_old, weights2_old;
    double lockWeightIndexes[lockWeightIndexesSize];
    int ref3dIndex;
    bool whileLoop = true;
    while (whileLoop)
    {
#ifdef DEBUG_MPI
        std::cerr << "[" << node->rank << "] May I write. lockIndex: " << lockIndex <<std::endl;
#endif

        MPI_Send(&lockIndex, 1, MPI_INT, 0, TAG_MAY_I_WRITE, MPI_COMM_WORLD);
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        switch (status.MPI_TAG)
        {
        case TAG_DO_NOT_DARE_TO_WRITE://I am free
            MPI_Recv(0, 0, MPI_INT, 0, MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);
            //usleep(100000);//microsecond
            break;
        case TAG_YES_YOU_MAY_WRITE://I am free

            MPI_Recv(lockWeightIndexes, lockWeightIndexesSize, MPI_DOUBLE, 0, MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);

            lockIndex = ROUND(lockWeightIndexes[index_lockIndex]);

#ifdef DEBUG_MPI

            std::cerr << "[" << node->rank << "] TAG_YES_YOU_MAY_WRITE.lockIndex: " << lockIndex <<std::endl;
            std::cerr <<  "lockWeightIndexes[index_weight]" << lockWeightIndexes[index_weight]<<std::endl;
#endif

            weight_old   = lockWeightIndexes[index_weight];
            weights1_old = lockWeightIndexes[index_weights1];
            weights2_old = lockWeightIndexes[index_weights2];
            ref3dIndex = ROUND(lockWeightIndexes[index_ref3d]);

            mpi_write(dirno, ref3dIndex, avg, avg1, avg2, SFclass, SFclass1, SFclass2,
                      SFclassDiscarded, w1, w2, weight_old, weights1_old, weights2_old);
            whileLoop=false;
            break;
        default:
            std::cerr << "process WRONG TAG RECEIVED at " << node->rank << std::endl;
            break;
        }
    }
#ifdef DEBUG_MPI
    std::cerr << "[" << node->rank << "] Finish writer " <<std::endl;
#endif

    lockWeightIndexes[index_lockIndex] = lockIndex;
    //lockWeightIndexes[index_weight]= w1 + w2;
    lockWeightIndexes[index_weight]=    w;
    lockWeightIndexes[index_weights1] = w1;
    lockWeightIndexes[index_weights2] = w2;
    lockWeightIndexes[index_ref3d] = ref3dIndex;

    MPI_Send(lockWeightIndexes, lockWeightIndexesSize, MPI_DOUBLE, 0, TAG_I_FINISH_WRITTING, MPI_COMM_WORLD);

}


void MpiProgAngularClassAverage::mpi_writeFile(
    Image<double> avg,
    size_t dirno,
    FileName fileNameStk,
    double w_old)
{
    FileName fn_tmp;
    double w = avg.weight();
    Image<double> old;

    if (w > 0.)
    {
        if (fileNameStk.exists())
        {

            std::cerr << "DEBUG_JM: Composing"     << std::endl;
            std::cerr << "DEBUG_JM: dirno: "       << dirno << std::endl;
            std::cerr << "DEBUG_JM: fileNameStk: " << fileNameStk << std::endl;
            fn_tmp.compose(dirno, fileNameStk);
            std::cerr << "DEBUG_JM: fn_tmp: "      << fn_tmp << std::endl;
            old.read(fn_tmp);

            //w_old = old.weight();
            //std::cerr << "DEBUG_JM: SECOND Imagew dirno old  new master branch: " <<  fn_tmp << " " << dirno << " "
            //                                      <<  w_old << " "
            //                                      << w << " "
            //                                      << std::endl;
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(old())
            {
                dAij(old(),i,j) = (w_old * dAij(old(),i,j) + dAij(avg(),i,j)) / (w_old + w);
            }
            old.setWeight(w_old + w);
            old.write(fileNameStk, dirno, true, WRITE_REPLACE);
        }
        else
        {
            //std::cerr << "DEBUG_JM: FIRST Imagew dirno old  new master branch: " <<  fn_tmp << " " << dirno << " "
            //                                      <<  w_old << " "
            //                                      << w << " "
            //                                      << std::endl;
            if (w != 1.)
                avg() /= w;
            avg.write(fileNameStk, dirno, true, WRITE_REPLACE);
        }
    }
}


void MpiProgAngularClassAverage::mpi_produceSideInfo()
{
    node->barrierWait();

    //init with 0 by default through memset

    Iempty().resizeNoCopy(Xdim, Ydim);
    Iempty().setXmippOrigin();
    node->barrierWait();

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
    corr.resizeNoCopy(P.getSampleNoOuterRing());
    rotAux.local_transformer.setReal(corr);
    rotAux.local_transformer.FourierTransform();

    // Set ring defaults
    if (Ri < 1)
        Ri = 1;
    if (Ro < 0)
        Ro = (Xdim / 2) - 1;
}

void MpiProgAngularClassAverage::mpi_preprocess()
{
    initFileNames();

    if (node->rank==0)
    {
        master_seed = randomize_random_generator();
    }
    MPI_Bcast(&master_seed,1,MPI_UNSIGNED ,0,MPI_COMM_WORLD);
    init_random_generator(master_seed);

    filterInputMetadata();

    if (node->rank==0)
    {
        std::cerr << "DEBUG_JM: saveDiscardedImages" <<std::endl;
        saveDiscardedImages();
        std::cerr << "DEBUG_JM: createJobList" <<std::endl;
        createJobList();
        std::cerr << "DEBUG_JM: initDimentions" <<std::endl;
        initDimentions();
        std::cerr << "DEBUG_JM: initWeights" <<std::endl;
        initWeights();
        std::cerr << "DEBUG_JM: initOutputFiles" <<std::endl;
        initOutputFiles();
    }

    MPI_Bcast(&Xdim,1,XMIPP_MPI_SIZE_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&Ydim,1,XMIPP_MPI_SIZE_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&Zdim,1,XMIPP_MPI_SIZE_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&Ndim,1,XMIPP_MPI_SIZE_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&numberOfJobs,1,XMIPP_MPI_SIZE_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&ref3dNum,1,XMIPP_MPI_SIZE_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&ctfNum,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&paddim,1,XMIPP_MPI_SIZE_T,0,MPI_COMM_WORLD);

    mpi_produceSideInfo();

    //    if (node->rank == 0)
    //    {
    //      std::cerr << "DEBUG_JM: end of mpi_preprocess" <<std::endl;
    //    }

    node->barrierWait();
}


void MpiProgAngularClassAverage::initFileNames()
{

    // Set up output rootnames
    if (do_split)
    {
        fn_out1 = fn_out + "_split_1";
        fn_out2 = fn_out + "_split_2";
    }

}

void MpiProgAngularClassAverage::filterInputMetadata()
{
    MetaData auxDF,auxF1;

    //auxDF.read((String)"ctfGroup[0-9][0-9][0-9][0-9][0-9][0-9]$@" + inFile);
    //std::cerr << "DEBUG_JM: inFile: " << inFile << std::endl;
    auxDF.read(inFile);
    if (!auxDF.containsLabel(MDL_REF3D))
        auxDF.fillConstant(MDL_REF3D, "1");
    if (!auxDF.containsLabel(MDL_DEFGROUP))
        auxDF.fillConstant(MDL_DEFGROUP, "1");
    //    if (!auxDF.containsLabel(MDL_ORDER))
    //    {
    //
    //        String cmd = formatString("%s=%s+%d", MDL::label2Str(MDL_ORDER).c_str(),
    //                                  MDL::label2Str(MDL_REF).c_str(), FIRST_IMAGE);
    //        auxDF.addLabel(MDL_ORDER);
    //        auxDF.operate(cmd);
    //    }
    if (!auxDF.containsLabel(MDL_ORDER))
        auxDF.addLabel(MDL_ORDER);
    {

        String cmd = formatString("%s=%s+%d", MDL::label2Str(MDL_ORDER).c_str(),
                                  MDL::label2Str(MDL_REF).c_str(), FIRST_IMAGE);
        auxDF.operate(cmd);
    }
    std::cerr << "DEBUG_JM: Read inFile" <<std::endl;

    MDMultiQuery multi;
    MDValueGE eq1(MDL::str2Label(col_select), limit0);
    MDValueLE eq2(MDL::str2Label(col_select), limitF);

    // remove percent of images
    if (do_limitR0per || do_limitRFper)
    {
        bool asc= !do_limitR0per;
        MDLabel codifyLabel = MDL::str2Label(col_select);
        int size = auxDF.size();
        int limit = size - ROUND((limitRper/100.) * size);
        auxF1.sort(auxDF,codifyLabel,asc,limit,0);
    }
    //remove images bellow (above) these limits
    else if(do_limit0 || do_limitF)
    {
        if (do_limit0)
            multi.addAndQuery(eq1);
        if (do_limitF)
            multi.addAndQuery(eq2);
        auxF1.importObjects(auxDF, multi);
    }
    // remove percentage of images from each class
    else if (do_limitR0class || do_limitRFclass)
    {
        //Delete a percentage of images in each class
        MetaData auxMdJobList;

        //!a take 1: using a copy of mdJobList
        const MDLabel myGroupByLabels[] =
            {
                MDL_REF3D, MDL_DEFGROUP, MDL_ORDER, MDL_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT
            };
        std::vector<MDLabel> groupbyLabels(myGroupByLabels,myGroupByLabels+6);
        auxMdJobList.aggregateGroupBy(auxDF, AGGR_COUNT, groupbyLabels, MDL_ORDER, MDL_COUNT);

        // Output stack size (number of valid projection directions)
        MDObject mdValueOut1(MDL_ORDER);
        auxMdJobList.aggregateSingleSizeT(mdValueOut1, AGGR_MAX ,MDL_ORDER);
        mdValueOut1.getValue(Ndim);

        MDObject mdValueOut2(MDL_DEFGROUP);
        auxMdJobList.aggregateSingleInt(mdValueOut2, AGGR_MAX ,MDL_DEFGROUP);
        mdValueOut2.getValue(ctfNum);

        MDObject mdValueOut3(MDL_REF3D);
        auxMdJobList.aggregateSingleInt(mdValueOut3, AGGR_MAX ,MDL_REF3D);
        mdValueOut3.getValue(ref3dNum);
        MultidimArray <int> multiCounter(Ndim+1, ctfNum+1, ref3dNum+1);
        multiCounter.initZeros();

        int ref3d, defgroup;
        size_t order, jobCount, jobCount2;

        FOR_ALL_OBJECTS_IN_METADATA(auxMdJobList)
        {
            auxMdJobList.getValue(MDL_REF3D, ref3d, __iter.objId);
            auxMdJobList.getValue(MDL_DEFGROUP, defgroup, __iter.objId);
            auxMdJobList.getValue(MDL_ORDER, order, __iter.objId);

            auxMdJobList.getValue(MDL_COUNT,jobCount, __iter.objId);

            jobCount2 = ROUND(limitRclass * jobCount);

            if (jobCount2 == 0)
                if(rnd_unif(0,1)<limitRclass)
                    jobCount2 = 1;

            dAkij(multiCounter, order, defgroup, ref3d) = jobCount2;
        }

        // sort
        if(do_limitR0class)
            auxF1.sort(auxDF, MDL::str2Label(col_select), true);
        else
            auxF1.sort(auxDF, MDL::str2Label(col_select), false);

        auxF1.removeDisabled();
        auxF1.addLabel(MDL_ENABLED);
        auxF1.setValueCol(MDL_ENABLED, 1);

        FOR_ALL_OBJECTS_IN_METADATA(auxF1)
        {
            auxF1.getValue(MDL_REF3D, ref3d, __iter.objId);
            auxF1.getValue(MDL_DEFGROUP, defgroup, __iter.objId);
            auxF1.getValue(MDL_ORDER, order, __iter.objId);

            if (dAkij(multiCounter, order, defgroup, ref3d) > 0)
            {
                auxF1.setValue(MDL_ENABLED,-1,__iter.objId);
                dAkij(multiCounter, order, defgroup, ref3d)--;
            }
        }
        auxF1.removeDisabled();
    }
    else
    {
        auxF1 = auxDF;
    }
    DF.write("/tmp/inputfileAfterREmove.xmd");
    DF.sort(auxF1, MDL_IMAGE);
}

void MpiProgAngularClassAverage::saveDiscardedImages()
{

    MetaData auxDF, auxDFsort;
    FileName fileNameXmd;
    std::stringstream comment;

    auxDF.read(inFile);

    comment << "Discarded images";
    if (do_limit0)
        comment << ". Min value = " << limit0;
    if (do_limitF)
        comment << ". Max value = " << limitF;
    if (do_limitR0per)
        comment << ". Drop " << limitRper*100 << "% of images with lower " << col_select;
    if (do_limitRFper)
        comment << ". Drop " << limitRper*100 << "% of images with higher " << col_select;
    if (do_limitR0class)
        comment << ". Drop " << limitRclass*100 << "% of images (per class) with lower " << col_select
        << ". If the ROUND(num_images_per_class * limitRclass / 100) == 0 then images are randomly dropped"
        << " so the percentage is satisfied";
    if (do_limitRFclass)
        comment << ". Drop " << limitRclass*100 << "% of images (per class) with higher " << col_select
        << ". If the ROUND(num_images_per_class * limitRclass / 100) == 0 then images are randomly dropped"
        << " so the percentage is satisfied";

    auxDF.subtraction(DF,MDL_IMAGE);
    auxDF.setComment(comment.str());
    formatStringFast(fileNameXmd, "discarded@%s_discarded.xmd", fn_out.c_str());

    auxDFsort.sort(auxDF,MDL::str2Label(col_select));
    auxDFsort.write(fileNameXmd);

}

void MpiProgAngularClassAverage::initDimentions()
{
    getImageSize(DF, Xdim, Ydim, Zdim, Ndim);

    // Output stack size (number of valid projection directions)
    MDObject mdValueOut1(MDL_ORDER);
    mdJobList.aggregateSingleSizeT(mdValueOut1, AGGR_MAX ,MDL_ORDER);
    mdValueOut1.getValue(Ndim);

    MDObject mdValueOut2(MDL_DEFGROUP);
    mdJobList.aggregateSingleInt(mdValueOut2, AGGR_MAX ,MDL_DEFGROUP);
    mdValueOut2.getValue(ctfNum);

    MDObject mdValueOut3(MDL_REF3D);
    mdJobList.aggregateSingleInt(mdValueOut3, AGGR_MAX ,MDL_REF3D);
    mdValueOut3.getValue(ref3dNum);

    //Check Wiener filter image has correct size
    if (fn_wien != "")
    {
        size_t x,y,z, n;
        getImageSize(fn_wien,x,y,z,n);

        // Get and check padding dimensions
        paddim = ROUND(pad * Xdim);
        if (x != paddim)
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         "Incompatible padding factor for this Wiener filter");
    }

}

void MpiProgAngularClassAverage::initWeights()
{
    lockArray.initZeros(Ndim+1,ref3dNum+1);//this should be 3D or 2D is enought ROB
    weightArray.initZeros(ctfNum+1,Ndim+1,ref3dNum+1);
    weightArrays1.initZeros(ctfNum+1,Ndim+1,ref3dNum+1);
    weightArrays2.initZeros(ctfNum+1,Ndim+1,ref3dNum+1);

}

void MpiProgAngularClassAverage::initOutputFiles()
{
    //alloc space for output files
    FileName fn_tmp;

    for (int i = 1; i <= ref3dNum; i++)
    {
        formatStringFast(fn_tmp, "_Ref3D_%03lu", i);

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


void MpiProgAngularClassAverage::mpi_postprocess()
{
    // Write class selfile to disc (even if its empty)
    //std::cerr << "DEBUG_JM: mpi_postprocess: " << std::endl;
    FileName imageName, fileNameXmd, blockNameXmd;
    FileName imageNames1, fileNameXmds1, imageNames2, fileNameXmds2;
    MetaData auxMd,auxMd2,auxMd3;
    size_t order_number;
    int ref3d;
    double weights1, weights2, weights;

    MDValueEQ eq1(MDL_REF3D, 0), eq2(MDL_ORDER, (size_t) 0);
    MDMultiQuery multi;

    String
    comment =
        (String) "This file contains a list of class averages with direction projections and weights.";
    auxMd2.setComment(comment);

    for (int i=1; i<=ref3dNum; i++)
    {

        auxMd.clear();
        auxMd2.clear();

        formatStringFast(fileNameXmd,
                         "Ref3D_%03lu@%s_Ref3D_%03lu.xmd", i, fn_out.c_str(), i);



        auxMd.importObjects(mdJobList, MDValueEQ(MDL_REF3D, i));

        const MDLabel myGroupByLabels[] =
            {
                MDL_REF3D, MDL_ORDER, MDL_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT
            };
        std::vector<MDLabel> groupbyLabels(myGroupByLabels,myGroupByLabels+5);
        auxMd2.aggregateGroupBy(auxMd, AGGR_SUM, groupbyLabels, MDL_COUNT, MDL_WEIGHT);

        auxMd2.addLabel(MDL_IMAGE, 1);
        auxMd2.addLabel(MDL_ANGLE_PSI);
        auxMd2.setValueCol(MDL_ANGLE_PSI,0.);

        FOR_ALL_OBJECTS_IN_METADATA(auxMd2)
        {
            auxMd2.getValue(MDL_ORDER, order_number,__iter.objId);
            formatStringFast(imageName, "%06lu@%s_Ref3D_%03lu.stk", order_number,  fn_out.c_str(), i);
            auxMd2.setValue(MDL_IMAGE, imageName,__iter.objId);
            weights = 0.;
            weights += dAkij(weightArray,0,order_number, i);
            int enable;
            if(weights == 0)
                enable= -1;
            else
                enable=  1;
            auxMd2.setValue(MDL_ENABLED,enable,__iter.objId);
            auxMd2.setValue(MDL_WEIGHT, weights,__iter.objId);

        }
        auxMd2.write(fileNameXmd);

        if (do_save_images_assigned_to_classes)
        {
            for (size_t j= 1; j <= Ndim; j++)
            {
                eq1.setValue(i);
                eq2.setValue(j);
                multi.clear();

                multi.addAndQuery(eq1);
                multi.addAndQuery(eq2);

                auxMd3.importObjects(DF, multi);

                if(auxMd3.size() != 0)
                {
                    formatStringFast(blockNameXmd, "orderGroup%06lu_%s", j, fileNameXmd.c_str());
                    auxMd3.write(blockNameXmd, MD_APPEND);
                }
            }
        }
        auxMd2.addLabel(MDL_ENABLED);
        auxMd2.setValueCol(MDL_ENABLED, 1);

        if(do_split)
        {
            MetaData auxMds1(auxMd2);
            MetaData auxMds2(auxMd2);

            formatStringFast(fileNameXmds1,
                             "Ref3D_%03lu@%s_Ref3D_%03lu.xmd", i, fn_out1.c_str(), i);

            formatStringFast(fileNameXmds2,
                             "Ref3D_%03lu@%s_Ref3D_%03lu.xmd", i, fn_out2.c_str(), i);

            FOR_ALL_OBJECTS_IN_METADATA2(auxMds1, auxMds2)
            {
                auxMds1.getValue(MDL_ORDER, order_number,__iter.objId);
                auxMds1.getValue(MDL_REF3D, ref3d,__iter.objId);

                //weights1 = dAij(weightArrays1,order_number, ref3d);
                weights1 = 0.;
                weights2 = 0.;
                //for (size_t var = 0; var < weightArrays1.zdim; ++var)
                {
                    weights1 += dAkij(weightArrays1,0,order_number, i);
                    weights2 += dAkij(weightArrays2,0,order_number, i);
                }


                if(weights1 == 0)
                    auxMds1.setValue(MDL_ENABLED,-1,__iter.objId);
                else
                {
                    auxMds1.setValue(MDL_WEIGHT, weights1,__iter.objId);
                    formatStringFast(imageName, "%06lu@%s_Ref3D_%03lu.stk", order_number,  fn_out1.c_str(), i);
                    auxMds1.setValue(MDL_IMAGE, imageName,__iter.objId);
                }

                //weights2 = dAij(weightArrays2,order_number, ref3d);
                if(weights2 == 0)
                    auxMds2.setValue(MDL_ENABLED,-1,__iter2.objId);
                else
                {
                    auxMds2.setValue(MDL_WEIGHT, weights2,__iter2.objId);
                    formatStringFast(imageName, "%06lu@%s_Ref3D_%03lu.stk", order_number,  fn_out2.c_str(), i);
                    auxMds2.setValue(MDL_IMAGE, imageName,__iter2.objId);
                }


            }
            auxMds1.removeDisabled();
            auxMds1.write(fileNameXmds1);
            auxMds2.removeDisabled();
            auxMds2.write(fileNameXmds2);
        }
    }
}

void MpiProgAngularClassAverage::createJobList()
{
    const MDLabel myGroupByLabels[] =
        {
            MDL_REF3D, MDL_DEFGROUP, MDL_ORDER, MDL_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT
        };
    std::vector<MDLabel> groupbyLabels(myGroupByLabels,myGroupByLabels+6);
    mdJobList.aggregateGroupBy(DF, AGGR_COUNT, groupbyLabels, MDL_ORDER, MDL_COUNT);
#define DEBUG
#ifdef DEBUG
    DF.write("/tmp/kk_DF.xmd");
    mdJobList.write("/tmp/kk_mdJobList.xmd");
#endif
#undef DEBUG

    numberOfJobs = mdJobList.size();
}


void MpiProgAngularClassAverage::getPolar(MultidimArray<double> &img,
        Polar<std::complex<double> > &fP, bool conjugated, float xoff,
        float yoff)
{
    MultidimArray<double> Maux;
    Polar<double> P;

    // Calculate FTs of polar rings and its stddev
    produceSplineCoefficients(BSPLINE3, Maux, img);
    P.getPolarFromCartesianBSpline(Maux, Ri, Ro, 3, xoff, yoff);
    fourierTransformRings(P, fP, global_plans, conjugated);
}

void MpiProgAngularClassAverage::reAlignClass(Image<double> &avg1,
        Image<double> &avg2, MetaData &SFclass1, MetaData &SFclass2,
        std::vector<Image<double> > imgs, std::vector<int> splits,
        std::vector<int> numbers, size_t dirno, double * my_output)
{
    Polar<std::complex<double> > fPref, fPrefm, fPimg;
    std::vector<double> ccfs(splits.size());
    MultidimArray<double> ang;
    MultidimArray<double> Mimg, Mref, Maux;
    double maxcorr, new_xoff=0., new_yoff=0.;
    double w1, w2, opt_flip = 0., opt_psi = 0.;
    bool do_discard;

    SFclass1.clear();
    SFclass2.clear();
    Mref = avg1() + avg2();
    //#define DEBUG
#ifdef DEBUG

    Image<double> auxImg;
    auxImg() = Mref;
    auxImg.write("ref.xmp");
#endif

    for (int iter = 0; iter < nr_iter; iter++)
    {
        // Initialize iteration
        getPolar(Mref, fPref, true);
        getPolar(Mref, fPrefm, false);
        avg1().initZeros();
        avg2().initZeros();
        w1 = w2 = 0.;

#ifdef DEBUG

        std::cerr<<" entering iter "<<iter<<std::endl;
#endif

        for (size_t imgno = 0; imgno < imgs.size(); imgno++)
        {
            do_discard = false;
            maxcorr = -99.e99;
            // Rotationally align
            getPolar(imgs[imgno](), fPimg, false, (float) -imgs[imgno].Xoff(),
                     (float) -imgs[imgno].Yoff());
            // A. Check straight image
            rotationalCorrelation(fPimg, fPref, ang, rotAux);
            for (size_t k = 0; k < XSIZE(corr); k++)
            {
                if (corr(k) > maxcorr)
                {
                    maxcorr = corr(k);
                    opt_psi = ang(k);
                    opt_flip = 0.;
                }
            }

            // B. Check mirrored image
            rotationalCorrelation(fPimg, fPrefm, ang, rotAux);
            for (size_t k = 0; k < XSIZE(corr); k++)
            {
                if (corr(k) > maxcorr)
                {
                    maxcorr = corr(k);
                    opt_psi = realWRAP(360. - ang(k), -180., 180.);
                    opt_flip = 1.;
                }
            }

            // Translationally align
            if (!do_discard)
            {
                if (opt_flip == 1.)
                {
                    // Flip experimental image
                    Matrix2D<double> A(3, 3);
                    A.initIdentity();
                    A(0, 0) *= -1.;
                    A(0, 1) *= -1.;
                    applyGeometry(LINEAR, Mimg, imgs[imgno](), A, IS_INV,
                                  DONT_WRAP);
                    selfRotate(BSPLINE3, Mimg, opt_psi, DONT_WRAP);
                }
                else
                    rotate(BSPLINE3, Mimg, imgs[imgno](), opt_psi, DONT_WRAP);
            }

            if (!do_discard)
            {
                ccfs[imgno] = correlationIndex(Mref, Mimg);
                imgs[imgno].setPsi(opt_psi);
                imgs[imgno].setFlip(opt_flip);
                imgs[imgno].setShifts(new_xoff, new_yoff);
                if (opt_flip == 1.)
                    imgs[imgno].setShifts(-new_xoff, new_yoff);

                // Add to averages
                if (splits[imgno] == 0)
                {
                    w1 += 1.;
                    avg1() += Mimg;
                }
                else if (splits[imgno] == 1)
                {
                    w2 += 1.;
                    avg2() += Mimg;
                }
            }
            else
            {
                splits[imgno] = -1;
                ccfs[imgno] = 0.;
            }
        }
        Mref = avg1() + avg2();

    }

    avg1.setWeight(w1);
    avg2.setWeight(w2);

    // Report the new angles, offsets and selfiles
    my_output[4] = imgs.size() * AVG_OUPUT_SIZE;
    for (size_t imgno = 0; imgno < imgs.size(); imgno++)
    {
        if (splits[imgno] < 0)
            my_output[imgno * AVG_OUPUT_SIZE + 5] = -numbers[imgno];
        else
            my_output[imgno * AVG_OUPUT_SIZE + 5] = numbers[imgno];
        my_output[imgno * AVG_OUPUT_SIZE + 6] = avg1.rot();
        my_output[imgno * AVG_OUPUT_SIZE + 7] = avg1.tilt();
        my_output[imgno * AVG_OUPUT_SIZE + 8] = imgs[imgno].psi();
        my_output[imgno * AVG_OUPUT_SIZE + 9] = imgs[imgno].Xoff();
        my_output[imgno * AVG_OUPUT_SIZE + 10] = imgs[imgno].Yoff();
        my_output[imgno * AVG_OUPUT_SIZE + 11] = (double) dirno;
        my_output[imgno * AVG_OUPUT_SIZE + 12] = imgs[imgno].flip();
        my_output[imgno * AVG_OUPUT_SIZE + 13] = ccfs[imgno];

        if (splits[imgno] == 0)
        {
            SFclass1.setValue(MDL_IMAGE, imgs[imgno].name(),
                              SFclass1.addObject());
        }
        else if (splits[imgno] == 1)
        {
            SFclass2.setValue(MDL_IMAGE, imgs[imgno].name(),
                              SFclass2.addObject());
        }
    }
}


void MpiProgAngularClassAverage::applyWienerFilter(MultidimArray<double> &img)
{
    MultidimArray<std::complex<double> > Faux;

    if (paddim > Xdim)
    {
        // pad real-space image
        int x0 = FIRST_XMIPP_INDEX(paddim);
        int xF = LAST_XMIPP_INDEX(paddim);
        img.selfWindow(x0, x0, xF, xF);
    }
    FourierTransform(img, Faux);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mwien)
    {
        dAij(Faux,i,j) *= dAij(Mwien,i,j);
    }
    InverseFourierTransform(Faux, img);
    if (paddim > Xdim)
    {
        // de-pad real-space image
        int x0 = FIRST_XMIPP_INDEX(Xdim);
        int xF = LAST_XMIPP_INDEX(Xdim);
        img.selfWindow(x0, x0, xF, xF);
    }
}

