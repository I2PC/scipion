/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.csic.es (2008)
 *
 * Unidad de Bioinformatica del Centro Nacional de Biotecnologia , CSIC
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

#include "angular_class_average.h"

// Read arguments ==========================================================
void ProgAngularClassAverage::readParams()
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

    if (checkParam("--preprocess"))
    {
        do_preprocess = true;
        number_3dref = getIntParam("--number_3dreferences");
    }
    else
        do_preprocess = false;

    if (checkParam("--postprocess"))
    {
        do_postprocess = true;
        number_3dref = getIntParam("--number_3dreferences");
    }
    else
        do_postprocess = false;

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
void ProgAngularClassAverage::defineParams()
{
    addUsageLine("Make class average images and corresponding selfiles from angular_projection_matching docfiles.");
    addSeeAlsoLine("angular_project_library, angular_projection_matching");

    addParamsLine("    -i <doc_file>          : Docfile with assigned angles for all experimental particles");
    addParamsLine("    --lib <doc_file>       : Docfile with angles used to generate the projection matching library");
    addParamsLine("    -o <root_name>         : Output rootname for class averages and selfiles");
    addParamsLine("   [--split ]              : Also output averages of random halves of the data");
    addParamsLine("   [--wien <img=\"\"> ]    : Apply this Wiener filter to the averages");
    addParamsLine("   [--pad <factor=1.> ]    : Padding factor for Wiener correction");
    addParamsLine("   [--write_selfiles]      : Write class selfiles to disc");
    addParamsLine("   [--number_3dreferences <n>] : Number of 3D references (only used with flag --postprocess_metadata).");
    addParamsLine("   [--postprocess]: Create block with average images filenames.");
    addParamsLine("   requires --number_3dreferences;");
    addParamsLine("   [--preprocess] : Delete auxiliary files from previous execution. Alloc disk space for output stacks");
    addParamsLine("   requires --number_3dreferences;");
    addParamsLine("   [--ctfNum <n=1>]        : Ctf group number");
    addParamsLine("   [--ref3dNum <n=1>]      : 3D reference number");

    addParamsLine("==+ IMAGE SELECTION BASED ON INPUT DOCFILE (select one between: limit 0, F and R ==");
    addParamsLine("   [--select <col=\"maxCC\">]     : Column to use for image selection (limit0, limitF or limitR)");
    addParamsLine("   [--limit0 <l0>]         : Discard images below <l0>");
    addParamsLine("   [--limitF <lF>]         : Discard images above <lF>");
    addParamsLine("   [--limitR <lR>]         : if (lR>0 && lR< 100): discard lowest  <lR> % in each class");
    addParamsLine("                           : if (lR<0 && lR>-100): discard highest <lR> % in each class");

    addParamsLine("==+ REALIGNMENT OF CLASSES ==");
    addParamsLine("   [--iter <nr_iter=0>]      : Number of iterations for re-alignment");
    addParamsLine("   [--Ri <ri=1>]             : Inner radius to limit rotational search");
    addParamsLine("   [--Ro <r0=-1>]            : Outer radius to limit rotational search");
    addParamsLine("                           : ro = -1 -> dim/2-1");
    addParamsLine("   [--max_shift <ms=999.>]        : Maximum shift (larger shifts will be set to 0)");
    addParamsLine("   [--max_shift_change <msc=999.>] : Discard images that change shift more in the last iteration ");
    addParamsLine("   [--max_psi_change <mps=360.>]   : Discard images that change psi more in the last iteration ");

    addExampleLine("Sample at default values and calculating output averages of random halves of the data",false);
    addExampleLine("xmipp_angular_class_average -i proj_match.doc --lib ref_angles.doc -o out_dir --split");

    addKeywords("class average images");
}

// Run ====================================================================
void ProgAngularClassAverage::run()
{
    int i, nmax, nr_ref, nr_images, reserve;
    double rot, tilt, psi, xshift, yshift;
    bool mirror;
    double w, w1, w2;
    Matrix1D<double> dataline(8);

    produceSideInfo();

    if (do_preprocess)
    {
        preprocess();
        return;
    }

    if (do_postprocess)
    {
        postprocess();
        return;
    }

    // Only for do_add: append input docfile to add_to docfile
    //do_Add is enabled if we are using different ctf groups
    //    if (do_add)
    //        //if (false)
    //    {
    //        FileName fn_tmp=fn_out+".xmd";
    //        if (fn_tmp.exists())
    //        {
    //            MetaData DFaux = DF;
    //            MetaData DFaux2(fn_tmp);
    //            // Don't do any fancy merging or sorting because those functions are really slow...
    //            DFaux.write("dfaux_run.xmd");
    //            DFaux2.write("dfaux2_run.xmd");
    //            DFaux.unionAll(DFaux2);
    //            DFaux2.clear();
    //            DFaux2.removeDuplicates(DFaux);
    //            DFaux2.write(fn_tmp);
    //        }
    //        else
    //        {
    //            DF.write(fn_tmp);
    //        }
    //    }


    // Making class averages

    // Reserve memory for output from class realignment
    if (nr_iter > 0)
        reserve = DF.size();
    else
        reserve = 0;
    //output_values looks horrible but it is useful for the MPI version
    double output_values[AVG_OUPUT_SIZE * reserve + 1];

    nr_ref = DFlib.size();
    init_progress_bar(nr_ref);

    // Loop over all classes
    //FIXME
    //What for, better check which classes has experimental images
    //most of the classes are empty!!!!!ROB
    //for (int dirno = 1; dirno <= nr_ref; dirno++)
    size_t dirno;
    int ref_number;

    //char ch;
    FOR_ALL_OBJECTS_IN_METADATA(DFclassesExp)
    {

        // Do the actual work
        DFclassesExp.getValue(MDL_ORDER, dirno, __iter.objId);
        processOneClass(dirno, output_values);

        // Output classes sel and doc files
        w = output_values[1];
        w1 = output_values[2];
        w2 = output_values[3];
        addClassAverage(dirno, w, w1, w2);

        // Fill new docfile (with params after realignment)
        size_t id;
        bool auxBool;
        if (nr_iter > 0)
        {
            nr_images = ROUND(output_values[4] / AVG_OUPUT_SIZE);
            for (int i = 0; i < nr_images; i++)
            {
                int this_image = ROUND(output_values[i*AVG_OUPUT_SIZE+5]);
                if (!(this_image < 0))
                {
                    //FIXME: The next line has no sense since the MDL_IMAGE is string
                    // and 'this_image' is of type int...
                    REPORT_ERROR(ERR_UNCLASSIFIED,
                                 "The next line has no sense since the MDL_IMAGE is string \
                                 and 'this_image' is of type int...");
                    id = DF.firstObject(MDValueEQ(MDL_IMAGE, this_image));

                    DF.setValue(MDL_ANGLE_ROT,
                                output_values[i * AVG_OUPUT_SIZE + 6], id);
                    DF.setValue(MDL_ANGLE_TILT,
                                output_values[i * AVG_OUPUT_SIZE + 7], id);
                    DF.setValue(MDL_ANGLE_PSI,
                                output_values[i * AVG_OUPUT_SIZE + 8], id);
                    DF.setValue(MDL_SHITF_X,
                                output_values[i * AVG_OUPUT_SIZE + 9], id);
                    DF.setValue(MDL_SHITF_Y,
                                output_values[i * AVG_OUPUT_SIZE + 10], id);
                    DF.setValue(MDL_REF,
                                output_values[i * AVG_OUPUT_SIZE + 11], id);
                    if (output_values[i * AVG_OUPUT_SIZE + 12] == 0)
                        auxBool = false;
                    else
                        auxBool = true;
                    DF.setValue(MDL_FLIP, auxBool, id);
                    DF.setValue(MDL_MAXCC,
                                output_values[i * AVG_OUPUT_SIZE + 13], id);
                }
            }
        }

        progress_bar((long) dirno);

    }
    progress_bar(nr_ref);

    // Write selfiles and docfiles with all class averages
    //if (write_selfiles)
    //    finalWriteToDisc();

}

// Side info stuff ===================================================================
void ProgAngularClassAverage::produceSideInfo()
{

    // Set up output rootnames
    if (do_split)
    {
        fn_out1 = fn_out + "_split_1";
        fn_out2 = fn_out + "_split_2";
    }

    int dummyI;
    size_t dummyT;
    getImageSize(DF, dim, dummyI, dummyI, dummyT);
    //init with 0 by default through memset
    Iempty().resizeNoCopy(dim, dim);
    Iempty().setXmippOrigin();

    // Randomization
    if (do_split)
        randomize_random_generator();

    // Set up FFTW transformers
    //FIXME
    //Is this needed if  no alignment is required?
    MultidimArray<double> Maux;
    Polar<double> P;
    Polar<std::complex<double> > fP;

    produceSplineCoefficients(BSPLINE3, Maux, Iempty());
    if (Ro < 0)
    	Ro = XSIZE(Maux)/2 - 1;
    P.getPolarFromCartesianBSpline(Maux, Ri, Ro);
    P.calculateFftwPlans(global_plans);
    fourierTransformRings(P, fP, global_plans, false);
    corr.resize(P.getSampleNoOuterRing());
    rotAux.local_transformer.setReal(corr);
    rotAux.local_transformer.FourierTransform();

}

void ProgAngularClassAverage::preprocess()
{
    //alloc space for output files
    int Xdim, Ydim, Zdim;
    size_t Ndim;
    FileName fn_tmp;

    getImageSize(DF, Xdim, Ydim, Zdim, Ndim);

    Ndim = DFclassesExp.size();

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
/*
 * We need two metadata files. Number 1 with each image assigned to each 3D reference
 * 2 with the winner corrected and averaged images plus weiths for 3D reconstruction (one per t3D reference)
 */
void ProgAngularClassAverage::postprocess()
{

    MetaData MD, MDsplit1, MDsplit2, MDdiscarted;
    MetaData auxMD;
    FileName fn_tmp, fn_tmp1, fn_tmp2, fn, fn1, fn2;
    size_t order = -1, last_order = -1, id;
    double rot, tilt;
    //for a particular 3D reference
    for (int iref = 1; iref <= number_3dref; iref++)
    {
        MD.clear();
        MDsplit1.clear();
        MDsplit2.clear();
        MDdiscarted.clear();

        //for a particular ctfgroup
        for (int ictf = 1; ictf <= ctfNum; ictf++)
        {
            formatStringFast(fn_tmp, "%s_ctfGroup%06lu_refGroup%06lu.xmd",
                             fn_out.c_str(), ictf, iref);
            //open file
            if (fn_tmp.exists())
            {
                StringVector blockList;
                getBlocksInMetaDataFile(fn_tmp, blockList);
                //different projection directions
                for (StringVector::iterator it = blockList.begin(); it
                     != blockList.end(); ++it)
                {
                    if ((*it).find("_") != std::string::npos)
                        continue;
                    auxMD.read(*it + '@' + fn_tmp);
                    MD.unionAll(auxMD);
                }
                //unlink(fn_tmp.c_str());
            }

            if (do_split)
            {
                formatStringFast(fn_tmp, "%s_ctfGroup%06lu_refGroup%06lu.xmd",
                                 fn_out1.c_str(), ictf, iref);

                if (fn_tmp.exists())
                {
                    StringVector blockList;
                    getBlocksInMetaDataFile(fn_tmp, blockList);
                    for (StringVector::iterator it = blockList.begin(); it
                         != blockList.end(); ++it)
                    {
                        if ((*it).find("_") != std::string::npos)
                            continue;
                        auxMD.read(*it + '@' + fn_tmp);
                        MDsplit1.unionAll(auxMD);
                    }
                    //unlink(fn_tmp.c_str());

                }
                formatStringFast(fn_tmp, "%s_ctfGroup%06lu_refGroup%06lu.xmd",
                                 fn_out2.c_str(), ictf, iref);

                if (fn_tmp.exists())
                {
                    StringVector blockList;
                    getBlocksInMetaDataFile(fn_tmp, blockList);
                    for (StringVector::iterator it = blockList.begin(); it
                         != blockList.end(); ++it)
                    {
                        if ((*it).find("_") != std::string::npos)
                            continue;
                        auxMD.read(*it + '@' + fn_tmp);
                        MDsplit2.unionAll(auxMD);
                    }
                    //unlink(fn_tmp.c_str());
                }
            }

            formatStringFast(fn_tmp,
                             "%s_ctfGroup%06lu_refGroup%06lu_discarded.xmd",
                             fn_out.c_str(), ictf, iref);

            if (fn_tmp.exists())
            {
                StringVector blockList;
                getBlocksInMetaDataFile(fn_tmp, blockList);
                for (StringVector::iterator it = blockList.begin(); it
                     != blockList.end(); ++it)
                {
                    if ((*it).find("_") != std::string::npos)
                        continue;
                    auxMD.read(*it + '@' + fn_tmp);
                    MDdiscarted.unionAll(auxMD);
                }
                //unlink(fn_tmp.c_str());
            }

        }

        formatStringFast(fn_tmp, "%s_refGroup%06lu.xmd", fn_out.c_str(), iref);
        MD.write(fn_tmp);
        if (do_split)
        {
            formatStringFast(fn_tmp, "%s_refGroup%06lu.xmd", fn_out1.c_str(),
                             iref);
            MDsplit1.write(fn_tmp);
            formatStringFast(fn_tmp, "%s_refGroup%06lu.xmd", fn_out2.c_str(),
                             iref);
            MDsplit2.write(fn_tmp);
        }
        formatStringFast(fn_tmp, "%s_refGroup%06lu_discarded.xmd",
                         fn_out.c_str(), iref);
        MDdiscarted.write(fn_tmp);
    }
}

void ProgAngularClassAverage::getPolar(MultidimArray<double> &img,
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

void ProgAngularClassAverage::reAlignClass(Image<double> &avg1,
        Image<double> &avg2, MetaData &SFclass1, MetaData &SFclass2,
        std::vector<Image<double> > imgs, std::vector<int> splits,
        std::vector<int> numbers, size_t dirno, double * my_output)
{
    Polar<std::complex<double> > fPref, fPrefm, fPimg;
    std::vector<double> ccfs(splits.size());
    MultidimArray<double> ang;
    MultidimArray<double> Mimg, Mref, Maux;
    double maxcorr, diff_psi, diff_shift, new_xoff, new_yoff;
    double w1, w2, opt_flip = 0., opt_psi = 0., opt_xoff = 0., opt_yoff = 0.;
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

    CorrelationAux aux;
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

        for (int imgno = 0; imgno < imgs.size(); imgno++)
        {
            do_discard = false;
            maxcorr = -99.e99;
            // Rotationally align
            getPolar(imgs[imgno](), fPimg, false, (float) -imgs[imgno].Xoff(),
                     (float) -imgs[imgno].Yoff());
            // A. Check straight image
            rotationalCorrelation(fPimg, fPref, ang, rotAux);
            for (int k = 0; k < XSIZE(corr); k++)
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
            for (int k = 0; k < XSIZE(corr); k++)
            {
                if (corr(k) > maxcorr)
                {
                    maxcorr = corr(k);
                    opt_psi = realWRAP(360. - ang(k), -180., 180.);
                    opt_flip = 1.;
                }
            }

            // Check max_psi_change in last iteration
            if (iter == nr_iter - 1)
            {
                diff_psi
                = ABS(realWRAP(imgs[imgno].psi() - opt_psi, -180., 180.));
                if (diff_psi > max_psi_change)
                {
                    do_discard = true;
#ifdef DEBUG
                    //std::cerr<<"discard psi: "<<diff_psi<<opt_psi<<" "<<opt_flip<<" "<<imgs[imgno].psi()<<std::endl;
#endif

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

                if (max_shift > 0)
                {
                    bestShift(Mref, Mimg, opt_xoff, opt_yoff, aux);
                    if (opt_xoff * opt_xoff + opt_yoff * opt_yoff > max_shift
                        * max_shift)
                    {
                        new_xoff = new_yoff = 0.;
                    }
                    else
                    {
                        selfTranslate(BSPLINE3, Mimg,
                                      vectorR2(opt_xoff, opt_yoff), true);
                        new_xoff = opt_xoff * COSD(opt_psi) + opt_yoff
                                   *SIND(opt_psi);
                        new_yoff = -opt_xoff * SIND(opt_psi) + opt_yoff
                                   *COSD(opt_psi);
                    }

                    // Check max_shift_change in last iteration
                    if (iter == nr_iter - 1)
                    {
                        opt_yoff = imgs[imgno].Yoff();
                        opt_xoff = imgs[imgno].Xoff();
                        if (imgs[imgno].flip() == 1.)
                            opt_xoff *= -1.;
                        diff_shift = (new_xoff - opt_xoff) * (new_xoff
                                                              - opt_xoff) + (new_yoff - opt_yoff) * (new_yoff
                                                                                                     - opt_yoff);
                        if (diff_shift > max_shift_change * max_shift_change)
                        {
                            do_discard = true;
#ifdef DEBUG

                            std::cerr <<"discard shift: "<<diff_shift<<" "<<new_xoff<<" "<<opt_xoff<<" "<<imgs[imgno].Xoff()<<" "<<new_yoff<<" "<<opt_yoff<<" "<<imgs[imgno].Yoff()<<std::endl;
#endif

                        }
                    }
                }
            }

            if (!do_discard)
            {
                ccfs[imgno] = correlationIndex(Mref, Mimg);
                imgs[imgno].setPsi(opt_psi);
                imgs[imgno].setFlip(opt_flip);
                imgs[imgno].setShifts(new_xoff, new_yoff);
                if (opt_flip == 1.)
                    imgs[imgno].setShifts(-new_xoff, new_yoff);

                // Check max_shift_change in last iteration
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
    for (int imgno = 0; imgno < imgs.size(); imgno++)
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

void ProgAngularClassAverage::applyWienerFilter(MultidimArray<double> &img)
{
    MultidimArray<std::complex<double> > Faux;
    if (paddim > dim)
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
    if (paddim > dim)
    {
        // de-pad real-space image
        int x0 = FIRST_XMIPP_INDEX(dim);
        int xF = LAST_XMIPP_INDEX(dim);
        img.selfWindow(x0, x0, xF, xF);
    }
}

void ProgAngularClassAverage::processOneClass(size_t &dirno, double * my_output)
{
    Image<double> img, avg, avg1, avg2;
    FileName fn_img, fn_tmp;
    MetaData SFclass, SFclass1, SFclass2;
    MetaData SFclassDiscarded;
    double rot, tilt, psi, xshift, yshift, val, w, w1, w2, my_limitR, scale;
    bool mirror;
    int ref_number, this_image;
    int isplit;
    MetaData _DF;
    size_t id;
    double current_rot, current_tilt;
    size_t order_number;

    w = 0.;
    w1 = 0.;
    w2 = 0.;
    this_image = 0;

    //_DF.importObjects(DF,MDValueEQ(MDL_REF,(int)dirno));
    _DF.importObjects(DF, MDValueEQ(MDL_ORDER, dirno));

    if (_DF.size() == 0)
    {//no images assigned to this class
        //this is possible because the input metadata contains
        //several blocks
        my_output[0] = (double) dirno;
        my_output[1] = w;
        my_output[2] = w1;
        my_output[3] = w2;
        return;
    }

    Matrix2D<double> A(3, 3);
    std::vector<int> exp_number, exp_split;
    std::vector<Image<double> > exp_imgs;
    //CHECK ANGLE REFS
    // Get reference angles and preset to averages
    DFlib.getValue(MDL_ANGLE_ROT, rot, dirno);
    DFlib.getValue(MDL_ANGLE_TILT, tilt, dirno);
    Iempty.setEulerAngles(rot, tilt, 0.);
    Iempty.setShifts(0., 0.);
    Iempty.setFlip(0.);
    avg = Iempty;
    avg1 = Iempty;
    avg2 = Iempty;

    // Loop over all images in the input docfile
    FOR_ALL_OBJECTS_IN_METADATA(_DF)
    {
        _DF.getValue(MDL_IMAGE, fn_img, __iter.objId);
        this_image++;

        _DF.getValue(MDL_REF, ref_number, __iter.objId);
        _DF.getValue(MDL_ANGLE_ROT, current_rot, __iter.objId);
        _DF.getValue(MDL_ANGLE_TILT, current_tilt, __iter.objId);
        _DF.getValue(MDL_ORDER, order_number, __iter.objId);
        //if (ref_number == dirno)
        {
            bool is_selected = true;
            double auxval;
            _DF.getValue(MDL::str2Label(col_select), auxval, __iter.objId);
            if (do_limit0)
            {
                if (auxval < limit0)
                    is_selected = false;
            }
            if (do_limitF)
            {
                if (auxval > limitF)
                    is_selected = false;
            }
            if (is_selected)
            {
                _DF.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
                _DF.getValue(MDL_SHITF_X, xshift, __iter.objId);
                _DF.getValue(MDL_SHITF_Y, yshift, __iter.objId);
                if (do_mirrors)
                    _DF.getValue(MDL_FLIP, mirror, __iter.objId);
                _DF.getValue(MDL_SCALE, scale, __iter.objId);

                //TODO: Check this????
                img.read(fn_img);
                img().setXmippOrigin();
                img.setEulerAngles(0., 0., psi);
                img.setShifts(xshift, yshift);
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
                    exp_split.push_back(isplit);
                }

                // Apply in-plane transformation
                img.getTransformationMatrix(A);
                if (!A.isIdentity())
                    selfApplyGeometry(BSPLINE3, img(), A, IS_INV, WRAP);

                // Add to average
//#define DEBUG
#ifdef DEBUG
                //WRITE IMAGES TO AVERAGE
                FileName fn_tmp1;
                int static static_i=0;
                static_i++;
                formatStringFast(fn_tmp1, "test_%06d", static_i);
                img.write(fn_tmp1);
                if (static_i> 25)
                	exit(1);
#endif
#undef DEBUG
                if (isplit == 0)
                {
                    avg1() += img();
                    w1 += 1.;
                    id = SFclass1.addObject();
                    SFclass1.setValue(MDL_IMAGE, fn_img, id);
                    SFclass1.setValue(MDL_ANGLE_ROT, current_rot, id);
                    SFclass1.setValue(MDL_ANGLE_TILT, current_tilt, id);
                    SFclass1.setValue(MDL_REF, ref_number, id);
                    SFclass1.setValue(MDL_ORDER, dirno, id);
                }
                else
                {
                    avg2() += img();
                    w2 += 1.;
                    id = SFclass2.addObject();
                    SFclass2.setValue(MDL_IMAGE, fn_img, id);
                    SFclass2.setValue(MDL_ANGLE_ROT, current_rot, id);
                    SFclass2.setValue(MDL_ANGLE_TILT, current_tilt, id);
                    SFclass2.setValue(MDL_REF, ref_number, id);
                    SFclass2.setValue(MDL_ORDER, dirno, id);
                }
            }
            else
            {
                id = SFclassDiscarded.addObject();
                SFclassDiscarded.setValue(MDL_IMAGE, fn_img, id);
                SFclassDiscarded.setValue(MDL_ANGLE_ROT, current_rot, id);
                SFclassDiscarded.setValue(MDL_ANGLE_TILT, current_tilt, id);
                SFclassDiscarded.setValue(MDL_REF, ref_number, id);
                SFclassDiscarded.setValue(MDL_ORDER, dirno, id);
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
        writeToDisc(avg, dirno, SFclass, fn_out, false, "ref.xmp");
        reAlignClass(avg1, avg2, SFclass1, SFclass2, exp_imgs, exp_split,
                     exp_number, dirno, my_output);
        w1 = avg1.weight();
        w2 = avg2.weight();
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
    w = w1 + w2;
    avg.setWeight(w);
    avg1.setWeight(w1);
    avg2.setWeight(w2);

    //ROB WRITE DISK
    // blocks
    FileName fileNameXmd, fileNameStk;

    formatStringFast(fileNameXmd,
                     "classGroup%06lu_refGroup%06lu@%s_ctfGroup%06lu_refGroup%06lu.xmd",
                     dirno, ref3dNum, fn_out.c_str(), ctfNum, ref3dNum);
    formatStringFast(fileNameStk, "%s_refGroup%06lu.stk", fn_out.c_str(),
                     ref3dNum);
    writeToDisc(avg, dirno, SFclass, fileNameXmd, fileNameStk, write_selfiles);

    if (do_split)
    {
        if (w1 > 0)
        {
            formatStringFast(
                fileNameXmd,
                "classGroup%06lu_refGroup%06lu@%s_ctfGroup%06lu_refGroup%06lu.xmd",
                dirno, ref3dNum, fn_out1.c_str(), ctfNum, ref3dNum);
            formatStringFast(fileNameStk, "%s_refGroup%06lu.stk",
                             fn_out1.c_str(), ref3dNum);
            writeToDisc(avg1, dirno, SFclass1, fileNameXmd, fileNameStk,
                        write_selfiles);
        }
        if (w2 > 0)
        {
            formatStringFast(
                fileNameXmd,
                "classGroup%06lu_refGroup%06lu@%s_ctfGroup%06lu_refGroup%06lu.xmd",
                dirno, ref3dNum, fn_out2.c_str(), ctfNum, ref3dNum);
            formatStringFast(fileNameStk, "%s_refGroup%06lu.stk",
                             fn_out2.c_str(), ref3dNum);
            writeToDisc(avg2, dirno, SFclass2, fileNameXmd, fileNameStk,
                        write_selfiles);
        }
    }

    formatStringFast(
        fileNameXmd,
        "classGroup%06lu_refGroup%06lu@%s_ctfGroup%06lu_refGroup%06lu_discarded.xmd",
        dirno, ref3dNum, fn_out.c_str(), ctfNum, ref3dNum);
    SFclassDiscarded.write(fileNameXmd, MD_APPEND);

    my_output[0] = (double) dirno;
    my_output[1] = w;
    my_output[2] = w1;
    my_output[3] = w2;
}

void ProgAngularClassAverage::writeToDisc(Image<double> avg, size_t dirno,
        MetaData SF, FileName fileNameXmd, FileName fileNameStk,
        bool write_selfile, FileName oext)
{
    FileName fn_tmp;
    double w = avg.weight(), w_old = 0;
    Image<double> old;
    MetaData SFold;

    if (w > 0.)
    {
        if (w != 1.)
            avg() /= w;
        //How are we going to handle weights, are they in the header, I do not think so....
        //in spider should be OK in general...!!
        //A more independent approach would be nice
        if (fileNameStk.exists())
        {
            fn_tmp.compose(dirno, fileNameStk);
            old.read(fn_tmp);
            w_old = old.weight();
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(old())
            {
                dAij(old(),i,j) = (w_old * dAij(old(),i,j) + w
                                   * dAij(avg(),i,j)) / (w_old + w);
            }
            old.setWeight(w_old + w);
            old.write(fileNameStk, dirno, true, WRITE_REPLACE);
        }
        else
        {
            avg.write(fileNameStk, dirno, true, WRITE_REPLACE);
        }
    }

    // Write class selfile to disc (even if its empty)
    if (write_selfile)
    {

        MetaData auxMd;
        String
        comment =
            (String) "This file contains one block per each pair (Projection direction, "
            + "Reference volume), groupClassXXXXXX refers to the projection direction, "
            + "refGroupXXXXXX refers to the reference volume. Additionally the blocks named "
            + "refGroupXXXXXX contains the images needed to reconstruct the reference XXXXXX.";
        SF.setComment(comment);
        //        }
        SF.write(fileNameXmd, MD_APPEND);
    }

}

void ProgAngularClassAverage::addClassAverage(size_t dirno, double w,
        double w1, double w2)
{
    double rot, tilt;
    FileName fn_tmp;

    DFlib.getValue(MDL_ANGLE_ROT, rot, dirno);
    DFlib.getValue(MDL_ANGLE_TILT, tilt, dirno);
    double d = 0.;
    bool f = false;
    size_t id;
    int ref2d;

    DFclassesExp.getValue(MDL_REF, ref2d, dirno);

    if (w > 0.)
    {
        fn_tmp.compose(dirno, fn_out + "_class", "stk");
        id = SFclasses.addObject();
        SFclasses.setValue(MDL_IMAGE, fn_tmp, id);
        SFclasses.setValue(MDL_ANGLE_ROT, rot, id);
        SFclasses.setValue(MDL_ANGLE_TILT, tilt, id);
        //this may help to keep compatibility with the next program
        SFclasses.setValue(MDL_ANGLE_PSI, d, id);
        SFclasses.setValue(MDL_SHITF_X, d, id);
        SFclasses.setValue(MDL_SHITF_Y, d, id);

        SFclasses.setValue(MDL_WEIGHT, w, id);
        //this may help to keep compatibility with the next program
        SFclasses.setValue(MDL_FLIP, f, id);

        SFclasses.setValue(MDL_REF, ref2d, id);

    }
    if (do_split)
    {
        if (w1 > 0.)
        {
            fn_tmp.compose(dirno, fn_out1 + "_class", "stk");
            id = SFclasses1.addObject();
            SFclasses1.setValue(MDL_IMAGE, fn_tmp, id);
            SFclasses1.setValue(MDL_ANGLE_ROT, rot, id);
            SFclasses1.setValue(MDL_ANGLE_TILT, tilt, id);
            //this may help to keep compatibility with the next program
            SFclasses1.setValue(MDL_ANGLE_PSI, d, id);
            SFclasses1.setValue(MDL_SHITF_X, d, id);
            SFclasses1.setValue(MDL_SHITF_Y, d, id);

            SFclasses1.setValue(MDL_WEIGHT, w1, id);

            //this may help to keep compatibility with the next program
            SFclasses1.setValue(MDL_FLIP, f, id);

            SFclasses1.setValue(MDL_REF, ref2d, id);

        }
        if (w2 > 0.)
        {
            fn_tmp.compose(dirno, fn_out2 + "_class", "stk");
            id = SFclasses2.addObject();
            SFclasses2.setValue(MDL_IMAGE, fn_tmp, id);
            SFclasses2.setValue(MDL_ANGLE_ROT, rot, id);
            SFclasses2.setValue(MDL_ANGLE_TILT, tilt, id);
            //this may help to keep compatibility with the next program
            SFclasses2.setValue(MDL_ANGLE_PSI, d, id);
            SFclasses2.setValue(MDL_SHITF_X, d, id);
            SFclasses2.setValue(MDL_SHITF_Y, d, id);

            SFclasses2.setValue(MDL_WEIGHT, w2, id);

            //this may help to keep compatibility with the next program
            SFclasses2.setValue(MDL_FLIP, f, id);

            SFclasses2.setValue(MDL_REF, ref2d, id);
        }
    }
}

//void ProgAngularClassAverage::finalWriteToDisc()
//{
//    FileName fn_tmp;
//    MetaData  auxSF, auxDF;
//
//    // Write docfiles with angles and weights of all classes
//    fn_tmp=fn_out+"_classes.xmd";
//    if (do_add && fn_tmp.exists())
//    {
//        MetaData MDaux;
//        MDaux.read(fn_tmp);
//        MDaux.unionAll(SFclasses);
//        SFclasses.aggregate(MDaux, AGGR_SUM, MDL_IMAGE, MDL_WEIGHT, MDL_SUM);
//    }
//    SFclasses.write(fn_tmp);
//    if (do_split)
//    {
//        fn_tmp=fn_out1+"_classes.xmd";
//        if (do_add && fn_tmp.exists())
//        {
//            MetaData MDaux;
//            MDaux.read(fn_tmp);
//            MDaux.unionAll(SFclasses1);
//            SFclasses1.aggregate(MDaux, AGGR_SUM,MDL_IMAGE,MDL_WEIGHT,MDL_SUM);
//        }
//        SFclasses1.write(fn_tmp);
//        fn_tmp=fn_out2+"_classes.xmd";
//        if (do_add && fn_tmp.exists())
//        {
//            MetaData MDaux;
//            MDaux.read(fn_tmp);
//            MDaux.unionAll(SFclasses2);
//            SFclasses2.aggregate(MDaux, AGGR_SUM,MDL_IMAGE,MDL_WEIGHT,MDL_SUM);
//        }
//        SFclasses2.write(fn_tmp);
//    }
//
//    // Write docfile with data for all realigned individual images
//    if (nr_iter > 0)
//    {
//        fn_tmp=fn_out+"_realigned.doc";
//        if (do_add && fn_tmp.exists())
//        {
//            DF.merge(fn_tmp);
//        }
//        DF.write(fn_tmp);
//    }
//}
