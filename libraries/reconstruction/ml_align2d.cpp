/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2007)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
#include "ml_align2d.h"
//#define DEBUG_JM

//Mutex for each thread update sums
pthread_mutex_t update_mutex =
    PTHREAD_MUTEX_INITIALIZER;
//Mutex for each thread get next refno
pthread_mutex_t refno_mutex =
    PTHREAD_MUTEX_INITIALIZER;

// Constructor
ProgML2D::ProgML2D()
{
    do_ML3D = false;
    refs_per_class = 1;
}


void ProgML2D::defineParams()
{
    addUsageLine("Perform (multi-reference) 2D-alignment using a maximum-likelihood (ML) target function.");
    addUsageLine("+Our recommended way of performing ML alignment is to introduce as little bias in the intial reference(s) as possible.");
    addUsageLine("+This can be done by calculting average images of random subsets of the (unaligned!) input experimental images, using the --nref option.");
    addUsageLine("+Note that the estimates for the standard deviation in the noise and in the origin offsets are re-estimated every iteration,");
    addUsageLine("+so that the initial values should not matter too much, as long as they are \"reasonable\". For Xmipp-normalized images,");
    addUsageLine("+the standard deviation in the noise can be assumed to be 1. For reasonably centered particles the default value of 3 for");
    addUsageLine("+the offsets should do the job.");
    addUsageLine("+");
    addUsageLine("+The output of the program consists of the refined reference images (weighted averages over all experimental images).");
    addUsageLine("+The experimental images are not altered at all. In terms of the ML approach, optimal transformations and references for");
    addUsageLine("+each image do not play the same role as in the conventional cross-correlation (or least-sqaures) approach. This program");
    addUsageLine("+can also be used for reference-free 2D-alignment using only a single reference: just supply =--nref 1= .");
    addUsageLine("+Although the calculations can be rather time-consuming (especially for many, large experimental images and a large number of references),");
    addUsageLine("+we strongly recommend to let the calculations converge. In our experience this takes in the order of 10-100 iterations, depending on the");
    addUsageLine("+number images, the amount of noise, etc. The default stopping criterium has yielded satisfactory results in our experience. A parallel ");
    addUsageLine("+version of this program has been implemented.");
    addSeeAlsoLine("mpi_ml_align2d");

    defineBasicParams(this);

    defineAdditionalParams(this, "==+ Additional options ==");
    defineHiddenParams(this);

    addExampleLine("A typical use of this program is:", false);
    addExampleLine("xmipp_ml_align2d -i input/images_some.stk --ref input/seeds2.stk --oroot output/ml2d --fast --mirror");
}

// Read arguments ==========================================================
void ProgML2D::readParams()
{
    // Generate new command line for restart procedure

    cline = "";
    int argc2 = 0;
    char ** argv2 = NULL;

    double restart_noise=0, restart_offset=0;
    FileName restart_imgmd, restart_refmd;
    int restart_iter=0, restart_seed=0;

    if (!do_ML3D  && checkParam("--restart"))
    {
        //TODO-------- Think later -----------
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Not implemented restart in ml2d for now");
        //        do_restart = true;
        //        MetaData MDrestart;
        //        char *copy  = NULL;
        //
        //        MDrestart.read(getParam("-restart"));
        //        cline = MDrestart.getComment();
        //        MDrestart.getValue(MDL_SIGMANOISE, restart_noise);
        //        MDrestart.getValue(MDL_SIGMAOFFSET, restart_offset);
        //        MDrestart.getValue(MDL_IMGMD, restart_imgmd);
        //        MDrestart.getValue(MDL_REFMD, restart_refmd);
        //        MDrestart.getValue(MDL_ITER, restart_iter);
        //        MDrestart.getValue(MDL_RANDOMSEED, restart_seed);
        //        generateCommandLine(cline, argc2, argv2, copy);
        //        //Take argument from the restarting command line
        //        progDef->read(argc2, argv2);
    }
    else
    {
        // no restart, just copy argc to argc2 and argv to argv2
        do_restart = false;
        for (int i = 1; i < argc; i++)
            cline = cline + (String) argv[i] + " ";
    }

    factor_nref = getIntParam("--nref");
    fn_ref = getParam("--ref");
    fn_img = getParam("-i");
    fn_root = getParam("--oroot");
    psi_step = getDoubleParam("--psi_step");
    Niter = getIntParam("--iter");
    //Fixed values for the 3D case
    if (do_ML3D)
    {
        istart = 1;
        fast_mode = save_mem2 = do_mirror = true;
    }
    else
    {
        //istart = getIntParam("--istart");
        istart = 1;
        do_mirror = checkParam("--mirror");
        save_mem2 = checkParam("--save_memB");
        fast_mode = checkParam("--fast");
    }
    model.sigma_noise = getDoubleParam("--noise");
    model.sigma_offset = getDoubleParam("--offset");

    eps = getDoubleParam("--eps");
    fn_frac = getParam("--frac");
    fix_fractions = checkParam("--fix_fractions");
    fix_sigma_offset = checkParam("--fix_sigma_offset");
    fix_sigma_noise = checkParam("--fix_sigma_noise");

    C_fast = getDoubleParam("-C");
    save_mem1 = checkParam("--save_memA");

    zero_offsets = checkParam("--zero_offsets");
    model.do_student = checkParam("--student");
    df = getDoubleParam("--student");
    model.do_norm = checkParam("--norm");

    // Number of threads
    threads = getIntParam("--thr");
    //testing the thread load in refno
    refno_load_param = getIntParam("--load");
    // Hidden arguments
    fn_scratch = getParameter(argc2, (const char **)argv2, "--scratch", "");
    debug = getIntParam("--debug");
    model.do_student_sigma_trick = !checkParam("--no_sigma_trick");
    trymindiff_factor = getDoubleParam("--trymindiff_factor");
    // Random seed to use for image randomization and creation
    // of initial references and blocks order
    // could be passed for restart or for debugging
    seed = getIntParam("--random_seed");

    if (seed == -1)
        seed = time(NULL);
    else if (!do_restart && verbose)
        std::cerr << "WARNING: *** Using a non random seed and not in restarting ***" <<std::endl;

    // Only for interaction with refine3d:
    search_rot = getDoubleParam("--search_rot");
    //IEM stuff
    blocks = getIntParam("--iem");

    // Now reset some stuff for restart
    if (do_restart)
    {
        fn_img = restart_imgmd;
        fn_ref = restart_refmd;
        model.n_ref = 0; // Just to be sure (not strictly necessary)
        model.sigma_noise = restart_noise;
        model.sigma_offset = restart_offset;
        seed = restart_seed;
        istart = restart_iter + 1;
        factor_nref = 1;
    }

    no_iem = checkParam("--no_iem");

    //std::cerr << "DEBUG_JM: exiting after readParams..." <<std::endl;
    //exit(1);
}

// Show ====================================================================
void ProgML2D::show()
{
    if (verbose)
    {
        // To screen
        if (!do_ML3D)
        {
            std::cout
            << " -----------------------------------------------------------------" << std::endl
            << " | Read more about this program in the following publications:   |" << std::endl
            << " |  Scheres ea. (2005) J.Mol.Biol. 348(1), 139-49                |" << std::endl
            << " |  Scheres ea. (2005) Bioinform. 21(suppl.2), ii243-4   (-fast) |" << std::endl
            << " |                                                               |"<<  std::endl
            << " |  *** Please cite them if this program is of use to you! ***   |"<<  std::endl
            << " -----------------------------------------------------------------"<<  std::endl;
        }

        std::cout    << "--> Maximum-likelihood multi-reference refinement " << std::endl
        << formatString("  Input images            : %s (%lu)\n", fn_img.c_str(), nr_images_global);

        if (fn_ref != "")
            std::cout << formatString("  Reference image(s)      : %s (%d)\n", fn_ref.c_str(), model.n_ref);
        if (factor_nref > 1)
            std::cout << "  Reference expanding factor   : " << factor_nref << std::endl;
        std::cout << "  Number of references:   : " << model.n_ref * factor_nref << std::endl;

        std::cout
        << "  Output rootname         : " << fn_root << std::endl
        << "  Stopping criterium      : " << eps << std::endl
        << "  Initial sigma noise     : " << model.sigma_noise << std::endl
        << "  Initial sigma offset    : " << model.sigma_offset << std::endl
        << "  Psi sampling interval   : " << psi_step << std::endl
        << "  Check mirrors           : " << (do_mirror ? "true" : "false") << std::endl;

        if (!fn_frac.empty())
            std::cout << "  Initial model fractions : " << fn_frac << std::endl;

        if (fast_mode)
        {
            std::cout << "  -> Use fast, reduced search-space approach with C = " << C_fast << std::endl;
            if (zero_offsets)
                std::cout << "    + Start from all-zero translations" << std::endl;
        }

        if (search_rot < 180.)
            std::cout << formatString("    + Limit orientational search to +/- %f degrees", search_rot) << std::endl;

        if (save_mem1)
            std::cout << "  -> Save_memory A: recalculate real-space rotations in -fast" << std::endl;

        if (save_mem2)
            std::cout << "  -> Save_memory B: limit translations to 3 sigma_offset " << std::endl;

        if (fix_fractions)
            std::cout << "  -> Do not update estimates of model fractions." << std::endl;

        if (fix_sigma_offset)
            std::cout << "  -> Do not update sigma-estimate of origin offsets." << std::endl;

        if (fix_sigma_noise)
            std::cout << "  -> Do not update sigma-estimate of noise." << std::endl;

        if (model.do_student)
        {
            std::cout << "  -> Use t-student distribution with df = " << df << std::endl;

            if (model.do_student_sigma_trick)
                std::cout << "  -> Use sigma-trick for t-student distributions" << std::endl;
        }

        if (model.do_norm)
            std::cout << "  -> Refine normalization for each experimental image" << std::endl;

        if (threads > 1)
            std::cout << "  -> Using " << threads << " parallel threads" << std::endl;

        if (blocks > 1)
            std::cout << "  -> Doing IEM with " << blocks << " blocks" <<std::endl;

        std::cout << " -----------------------------------------------------------------" << std::endl;

    }

}

void ProgML2D::printModel(const String &msg, const ModelML2D & model)
{
    std::cerr << "================> " << msg << std::endl;
    model.print();
}


//#define LOG(str)
// Trying to merge produceSideInfo 1 y 2
void ProgML2D::produceSideInfo()
{
    LOG("      ProgML2D::produceSideInfo: start");
    // Read selfile with experimental images
    // and set some global variables
    LOG("      ProgML2D::produceSideInfo: reading MDimg");
    MDimg.read(fn_img);
    // Remove disabled images
    MDimg.removeDisabled();
    nr_images_global = MDimg.size();
    // By default set myFirst and myLast equal to 0 and N
    // respectively, this should be changed when using MPI
    // by calling setWorkingImages before produceSideInfo2
    myFirstImg = 0;
    myLastImg = nr_images_global - 1;

    //Initialize blocks
    current_block = 0;

    // Create a vector of objectIDs, which may be randomized later on
    MDimg.findObjects(img_id);
    // Get original image size
    size_t idum, idumLong;
    LOG("      ProgML2D::produceSideInfo: setting dimensions");
    getImageSize(MDimg, dim, idum, idum, idumLong);
    model.dim = dim;
    hdim = dim / 2;
    dim2 = dim * dim;
    ddim2 = (double) dim2;
    double sigma_noise2 = model.sigma_noise * model.sigma_noise;

    if (model.do_student)
    {
        df2 = -(df + ddim2) / 2.;
        dfsigma2 = df * sigma_noise2;
    }

    if (fn_ref.empty())
    {
        //generate an initial reference just by averaging the experimental images
        LOG("      ProgML2D::produceSideInfo: generate an initial reference just by averaging the experimental images");
        FileName fn_tmp;
        Image<double> img, avg(dim, dim);
        avg().initZeros();
        avg().setXmippOrigin();

        FOR_ALL_OBJECTS_IN_METADATA(MDimg)
        {
            MDimg.getValue(MDL_IMAGE, fn_tmp, __iter.objId);
            img.read(fn_tmp);
            img().setXmippOrigin();
            avg() += img();
        }

        avg() /= nr_images_global;
        model.setNRef(1);
        model.Iref[0] = avg;
        fn_ref = fn_root + "_images_average.xmp";
        avg.write(fn_ref);
    }

    // Print some output to screen
    LOG("      ProgML2D::produceSideInfo show");
    show();
    LOG("      ProgML2D::produceSideInfo end");
}

void ProgML2D::produceSideInfo2()
{
    Image<double> img;
    FileName fn_tmp;

    // Read in all reference images in memory
    MDref.read(fn_ref);
    MDref.removeDisabled();

    model.setNRef(MDref.size());
    int refno = 0;
    double fraction = (double) 1./ model.n_ref;
    double weight = fraction * (double) nr_images_global;
    double rot = 0., tilt = 0.;

    FOR_ALL_OBJECTS_IN_METADATA(MDref)
    {
        MDref.getValue(MDL_IMAGE, fn_tmp, __iter.objId);
        MDref.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
        MDref.getValue(MDL_ANGLE_TILT, rot, __iter.objId);
        img.read(fn_tmp);
        img.setEulerAngles(rot, tilt, 0.);
        img().setXmippOrigin();
        img.setWeight(weight);
        model.Iref[refno] = img;
        //model.WsumMref[refno] = img;
        // Default start is all equal model fractions
        model.alpha_k[refno] = fraction;
        //model.WsumMref[refno]() *= weight;
        // Default start is half-half mirrored images
        model.mirror_fraction[refno] = (do_mirror ? 0.5 : 0.);
        ++refno;
    }

    setNumberOfLocalImages();
    // prepare masks for rotated references
    mask.resize(dim, dim);
    mask.setXmippOrigin();
    BinaryCircularMask(mask, hdim, INNER_MASK);
    omask.resize(dim, dim);
    omask.setXmippOrigin();
    BinaryCircularMask(omask, hdim, OUTSIDE_MASK);

    // Construct matrices for 0, 90, 180 & 270 degree flipping and mirrors
    initSamplingStuff();

    // Set sigdim, i.e. the number of pixels that will be considered in the translations
    sigdim = 2 * (size_t)ceil(model.sigma_offset * (save_mem2 ? 3 : 6));
    ++sigdim; // (to get uneven number)
    sigdim = XMIPP_MIN(dim, sigdim);

    //Some vectors and matrixes initialization
    int num_output_refs = model.n_ref * factor_nref;
    //std::cerr << "DEBUG_JM: num_output_refs: " << num_output_refs << std::endl;
    refw.resize(num_output_refs);
    refw2.resize(num_output_refs);
    refwsc2.resize(num_output_refs);
    refw_mirror.resize(num_output_refs);
    sumw_refpsi.resize(num_output_refs * nr_psi);
    A2.resize(num_output_refs);
    fref.resize(num_output_refs * nr_psi);
    mref.resize(num_output_refs * nr_psi);
    wsum_Mref.resize(num_output_refs);
    mysumimgs.resize(num_output_refs * nr_psi);
    Iold.resize(num_output_refs);

    if (fast_mode)
    {
        int mysize = num_output_refs * (do_mirror ? 2 : 1);
        //if (do_mirror)
        //    mysize *= 2;
        ioptx_ref.resize(mysize);
        iopty_ref.resize(mysize);
        ioptflip_ref.resize(mysize);
    }

    randomizeImagesOrder();

    // Initialize trymindiff for all images
    imgs_optrefno.clear();
    imgs_trymindiff.clear();
    imgs_scale.clear();
    imgs_bgmean.clear();
    imgs_offsets.clear();
    imgs_oldphi.clear();
    imgs_oldtheta.clear();
    model.scale.clear();

    if (model.do_norm)
    {
        average_scale = 1.;
        for (int refno = 0; refno < model.n_ref; refno++)
        {
            model.scale.push_back(1.);
        }
    }

    // Initialize imgs_offsets vectors
    std::vector<double> Vdum;
    double offx = (zero_offsets ? 0. : -999);
    int idum = (do_mirror ? 4 : 2) * factor_nref * model.n_ref;

    FOR_ALL_LOCAL_IMAGES()
    {
        imgs_optrefno.push_back(0);
        imgs_trymindiff.push_back(-1.);
        imgs_offsets.push_back(Vdum);
        for (int refno = 0; refno < idum; refno++)
        {
            imgs_offsets[IMG_LOCAL_INDEX].push_back(offx);
        }
        if (model.do_norm)
        {
            // Initialize scale and bgmean for all images
            // (for now initialize to 1 and 0, below also include doc)
            imgs_bgmean.push_back(0.);
            imgs_scale.push_back(1.);
        }
        if (limit_rot)
        {
            // For limited orientational search: initialize imgs_oldphi & imgs_oldtheta to -999.
            imgs_oldphi.push_back(-999.);
            imgs_oldtheta.push_back(-999.);
        }
    }

    //TODO: THINK RESTART LATER, NOW NOT WORKING
    if (do_restart)
    {
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Restart option not implement yet in ml2d");

        // Read optimal image-parameters
        //        FOR_ALL_LOCAL_IMAGES()
        //        {
        //            if (limit_rot)
        //            {
        //                MDimg.getValue(MDL_ANGLE_ROT, imgs_oldphi[IMG_LOCAL_INDEX]);
        //                MDimg.getValue(MDL_ANGLE_TILT, imgs_oldtheta[IMG_LOCAL_INDEX]);
        //            }
        //            if (zero_offsets)
        //            {
        //                idum = (do_mirror ? 2 : 1) * model.n_ref;
        //                double xx, yy;
        //                MDimg.getValue(MDL_SHIFT_X, xx);
        //                MDimg.getValue(MDL_SHIFT_X, yy);
        //                for (int refno = 0; refno < idum; refno++)
        //                {
        //                    imgs_offsets[IMG_LOCAL_INDEX][2 * refno] = xx;
        //                    imgs_offsets[IMG_LOCAL_INDEX][2 * refno + 1] = yy;
        //                }
        //            }
        //            if (model.do_norm)
        //            {
        //                MDimg.getValue(MDL_BGMEAN, imgs_bgmean[IMG_LOCAL_INDEX]);
        //                MDimg.getValue(MDL_INTSCALE, imgs_scale[IMG_LOCAL_INDEX]);
        //            }
        //        }
        //        // read Model parameters
        //        int refno = 0;
        //        FOR_ALL_OBJECTS_IN_METADATA(MDref)
        //        {
        //            MDref.getValue(MDL_MODELFRAC, model.alpha_k[refno]);
        //            if (do_mirror)
        //                MDref.getValue(MDL_MIRRORFRAC,
        //                               model.mirror_fraction[refno]);
        //            if (model.do_norm)
        //                MDref.getValue(MDL_INTSCALE, model.scale[refno]);
        //            refno++;
        //        }

    }
    // Call the first iteration 0 if generating initial references from random subsets
    else if (factor_nref > 1)
    {
        istart = SPECIAL_ITER;

        //Expand MDref because it will be re-used in writeOutputFiles
        MetaData MDaux(MDref);

        for (int group = 1; group < factor_nref; ++group)
        {
            FOR_ALL_OBJECTS_IN_METADATA(MDaux)
            {
                MDaux.setValue(MDL_REF3D, group + 1, __iter.objId);
            }
            MDref.unionAll(MDaux);
        }

        //std::cerr << "DEBUG_JM: MDref: " << std::endl;
        //MDref.write(std::cerr);
    }

    //--------Setup for Docfile -----------
    docfiledata.resize(nr_images_local, DATALINELENGTH);

}//close function produceSideInfo

// Calculate probability density function of all in-plane transformations phi
void ProgML2D::calculatePdfInplane()
{
#ifdef DEBUG
    std::cerr << "Entering calculatePdfInplane" <<std::endl;
#endif

    double x, y, r2, pdfpix, sum;
    P_phi.resize(dim, dim);
    P_phi.setXmippOrigin();
    Mr2.resize(dim, dim);
    Mr2.setXmippOrigin();

    sum = 0.;
    double _2sigma_offset2 = 2 * model.sigma_offset * model.sigma_offset;
    double _sigma_offset_denominator =  2 * PI * model.sigma_offset * model.sigma_offset * nr_psi * nr_nomirror_flips;

    FOR_ALL_ELEMENTS_IN_ARRAY2D(P_phi)
    {
        x = (double) j;
        y = (double) i;
        r2 = x * x + y * y;

        if (model.sigma_offset > 0.)
        {
            pdfpix = exp(-r2 / _2sigma_offset2);
            pdfpix /= _sigma_offset_denominator;
        }
        else
        {
            if (j == 0 && i == 0)
                pdfpix = 1.;
            else
                pdfpix = 0.;
        }

        A2D_ELEM(P_phi, i, j) = pdfpix;
        A2D_ELEM(Mr2, i, j) = (double) r2;
        sum += pdfpix;
    }

    // Normalization
    P_phi /= sum;
#ifdef DEBUG

    std::cerr << "Leaving calculatePdfInplane" <<std::endl;
#endif

}

// Rotate reference for all models and rotations and fill Fref vectors =============
void ProgML2D::rotateReference()
{
#ifdef DEBUG
    std::cerr<<"entering rotateReference"<<std::endl;
#endif

    awakeThreads(TH_RR_REFNO, 0, refno_load_param);

#ifdef DEBUG

    std::cerr<<"leaving rotateReference"<<std::endl;
#endif
}

// Collect all rotations and sum to update Iref() for all models ==========
void ProgML2D::reverseRotateReference()
{

#ifdef DEBUG
    std::cerr<<"entering reverseRotateReference"<<std::endl;
#endif

    awakeThreads(TH_RRR_REFNO, 0, refno_load_param);

#ifdef DEBUG

    std::cerr<<"leaving reverseRotateReference"<<std::endl;
#endif

}

void ProgML2D::preselectLimitedDirections(double &phi, double &theta)
{

    double phi_ref, theta_ref, angle, angle2;
    Matrix1D<double> u, v;

    pdf_directions.clear();
    pdf_directions.resize(model.n_ref);

    for (int refno = 0; refno < model.n_ref; refno++)
    {
        if (!limit_rot || (phi == -999. && theta == -999.))
            pdf_directions[refno] = 1.;
        else
        {
            phi_ref = model.Iref[refno].rot();
            theta_ref = model.Iref[refno].tilt();
            Euler_direction(phi, theta, 0., u);
            Euler_direction(phi_ref, theta_ref, 0., v);
            u.selfNormalize();
            v.selfNormalize();
            angle = RAD2DEG(acos(dotProduct(u, v)));
            angle = fabs(realWRAP(angle, -180, 180));
            // also check mirror
            angle2 = 180. + angle;
            angle2 = fabs(realWRAP(angle2, -180, 180));
            angle = XMIPP_MIN(angle, angle2);

            if (fabs(angle) > search_rot)
                pdf_directions[refno] = 0.;
            else
                pdf_directions[refno] = 1.;
        }
    }

}

// Pre-selection of significant refno and ipsi, based on current optimal translation =======


void ProgML2D::preselectFastSignificant()
{

#ifdef DEBUG
    std::cerr<<"entering preselectFastSignificant"<<std::endl;
#endif

    // Initialize Msignificant to all zeros
    // TODO: check whether this is strictly necessary? Probably not...
    Msignificant.initZeros();
    pfs_mindiff = 99.e99;
    pfs_maxweight.resizeNoCopy((do_mirror ? 2 : 1), model.n_ref);
    pfs_maxweight.initConstant(-99.e99);
    pfs_weight.resizeNoCopy(model.n_ref, nr_psi * nr_flip);
    pfs_weight.initZeros();
    pfs_count = threads;
    awakeThreads(TH_PFS_REFNO, 0, refno_load_param);

#ifdef DEBUG

    std::cerr<<"leaving preselectFastSignificant"<<std::endl;
#endif
}

// Maximum Likelihood calculation for one image ============================================
// Integration over all translation, given  model and in-plane rotation
void ProgML2D::expectationSingleImage(Matrix1D<double> &opt_offsets)
{
#ifdef TIMING
    timer.tic(ESI_E1);
#endif

    MultidimArray<double> Maux, Mweight;
    MultidimArray<std::complex<double> > Faux;
    double my_mindiff;
    bool is_ok_trymindiff = false;
    double sigma_noise2 = model.sigma_noise * model.sigma_noise;
    FourierTransformer local_transformer;
    ioptx = iopty = 0;

    // Update sigdim, i.e. the number of pixels that will be considered in the translations
    sigdim = 2 * CEIL(XMIPP_MAX(1,model.sigma_offset) * (save_mem2 ? 3 : 6));
    sigdim++; // (to get uneven number)
    sigdim = XMIPP_MIN(dim, sigdim);

    // Setup matrices
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Mweight.initZeros(sigdim, sigdim);
    Mweight.setXmippOrigin();

    if (!model.do_norm)
        opt_scale = 1.;

    // precalculate all flipped versions of the image
    Fimg_flip.clear();

    for (size_t iflip = 0; iflip < nr_flip; iflip++)
    {
        Maux.setXmippOrigin();
        applyGeometry(LINEAR, Maux, Mimg, F[iflip], IS_INV, WRAP);
        local_transformer.FourierTransform(Maux, Faux, false);

        if (model.do_norm)
            dAij(Faux,0,0) -= bgmean;

        Fimg_flip.push_back(Faux);

    }

    // The real stuff: loop over all references, rotations and translations
    int redo_counter = 0;

#ifdef TIMING

    timer.toc(ESI_E1);

    //timer.tic("WHILE_");
#endif

    while (!is_ok_trymindiff)
    {
        // Initialize mindiff, weighted sums and maxweights
        mindiff = 99.e99;
        wsum_corr = wsum_offset = wsum_sc = wsum_sc2 = 0.;
        maxweight = maxweight2 = sum_refw = sum_refw2 = 0.;

        awakeThreads(TH_ESI_REFNO, opt_refno, refno_load_param);

        // Now check whether our trymindiff was OK.
        // The limit of the exp-function lies around
        // exp(700)=1.01423e+304, exp(800)=inf; exp(-700) = 9.85968e-305; exp(-88) = 0
        // Use 500 to be on the save side?

        if (ABS((mindiff - trymindiff) / sigma_noise2) > 500.)
            //force always redo to use real mindiff for check about LL problem
            //if (redo_counter==0)
        {
            // Re-do whole calculation now with the real mindiff
            trymindiff = mindiff;
            redo_counter++;
            // On iteration 0 images that will go to references other than first
            // will store optimus references but not yet expanded number of references
            if (iter == 0)
                opt_refno = (opt_refno % model.n_ref);

            // Never re-do more than once!
            if (redo_counter > 1)
                REPORT_ERROR(ERR_VALUE_INCORRECT, "ml_align2d BUG% redo_counter > 1");
        }
        else
        {
            is_ok_trymindiff = true;
            my_mindiff = trymindiff;
            trymindiff = mindiff;
        }

    }//close while

    fracweight = maxweight / sum_refw;

    wsum_sc /= sum_refw;

    wsum_sc2 /= sum_refw;

    // Calculate optimal transformation parameters
    opt_psi = -psi_step * (iopt_flip * nr_psi + iopt_psi) - SMALLANGLE;

    opt_offsets(0) = -(double) ioptx * MAT_ELEM(F[iopt_flip], 0, 0)
                     - (double) iopty * MAT_ELEM(F[iopt_flip], 0, 1);

    opt_offsets(1) = -(double) ioptx * MAT_ELEM(F[iopt_flip], 1, 0)
                     - (double) iopty * MAT_ELEM(F[iopt_flip], 1, 1);

    // Update normalization parameters
    if (model.do_norm)
    {
        // 1. Calculate optimal setting of Mimg
        MultidimArray<double> Maux2 = Mimg;
        selfTranslate(LINEAR, Maux2, opt_offsets, true);
        selfApplyGeometry(LINEAR, Maux2, F[iopt_flip], IS_INV, WRAP);
        // 2. Calculate optimal setting of Mref
        int refnoipsi = (opt_refno % model.n_ref) * nr_psi + iopt_psi;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Faux)
        {
            dAij(Faux,i,j) = conj(dAij(fref[refnoipsi],i,j));
            dAij(Faux,i,j) *= opt_scale;
        }

        // Still take input from Faux and leave output in Maux
        local_transformer.inverseFourierTransform();
        Maux2 = Maux2 - Maux;

        if (debug == 12)
        {
            std::cout << std::endl;
            std::cout << "scale= " << opt_scale << " changes to " << wsum_sc
            / wsum_sc2 << std::endl;
            std::cout << "bgmean= " << bgmean << " changes to "
            << Maux2.computeAvg() << std::endl;
        }

        // non-ML update of bgmean (this is much cheaper than true-ML update...)
        old_bgmean = bgmean;

        bgmean = Maux2.computeAvg();

        // ML-update of opt_scale
        opt_scale = wsum_sc / wsum_sc2;
    }

#ifdef TIMING
    timer.toc(ESI_E5);

    timer.tic(ESI_E6TH);

#endif
    // Update all global weighted sums after division by sum_refw
    wsum_sigma_noise += (2 * wsum_corr / sum_refw);

    wsum_sigma_offset += (wsum_offset / sum_refw);

    sumfracweight += fracweight;


    //    std::cerr << "-------------------- wsum_corr: " << wsum_corr << std::endl;
    //    std::cerr << "-------------------- sum_refw: " << sum_refw << std::endl;
    //    std::cerr << std::endl << "-------------------- wsum_sigma_offset: " << wsum_sigma_offset << std::endl;
    //    std::cerr << "-------------------- wsum_sigma_noise: " << wsum_sigma_noise << std::endl << std::endl;
    //    if (wsum_sigma_noise < 0)
    //    {
    //      std::cerr << "Negative wsum_sigma_noise....exiting"  <<std::endl;
    //      exit(1);
    //    }

    awakeThreads(TH_ESI_UPDATE_REFNO, 0, refno_load_param);

    if (!model.do_student)
        // 1st term: log(refw_i)
        // 2nd term: for subtracting mindiff
        // 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
        dLL = log(sum_refw) - my_mindiff / sigma_noise2 - ddim2 * log(sqrt(2.
                * PI * sigma_noise2));
    else
        // 1st term: log(refw_i)
        // 2nd term: for dividing by (1 + 2. * mindiff/dfsigma2)^df2
        // 3rd term: for sigma-dependent normalization term in t-student distribution
        // 4th&5th terms: gamma functions in t-distribution
        dLL = log(sum_refw) + df2 * log(1. + (2. * my_mindiff / dfsigma2))
              - ddim2 * log(sqrt(PI * df * sigma_noise2)) + gammln(-df2)
              - gammln(df / 2.);

    //#define DEBUG_JM1
#ifdef DEBUG_JM1

    //if (iter>1)
    {
        //static std::ios::openmode mode = std::ios::out|std::ios::trunc;
        //std::ofstream log;
        //FileName fn = fn_root + (current_image % 2 == 0 ? "_images.even" : "_images.odd");
        //log.open(fn.c_str(), mode);
        //std::cerr << "    IMAGE " << current_image << "----------------------------->>>" << std::endl;
        std::cerr << "----------------------------->>>" << std::endl;
        std::cerr << "                             dLL: " << dLL << std::endl;
        std::cerr << "                        sum_refw: " << sum_refw << std::endl;
        std::cerr << "                      my_mindiff: " << my_mindiff << std::endl;
        std::cerr << "                    sigma_noise2: " << sigma_noise2 << std::endl;
        std::cerr << "                           ddim2: " << ddim2 << std::endl;
        //std::cerr << "                        dfsigma2: " << dfsigma2 << std::endl;
        std::cerr << "                       wsum_corr: " << wsum_corr << std::endl;
        std::cerr << "                     wsum_offset: " << wsum_offset << std::endl;
        //        std::cerr << "                            refw: ";
        //        for (int refno = 0; refno < model.n_ref; ++refno)
        //            std::cerr << std::setw(15) << refw[refno];
        //        std::cerr<< std::endl;
        //        std::cerr << "                     refw_mirror: ";
        //        for (int refno = 0; refno < model.n_ref; ++refno)
        //            std::cerr << std::setw(15) << refw_mirror[refno];
        //        std::cerr<< std::endl;
        //        log.close();
        //        mode = std::ios::out|std::ios::app;
    }
#endif
#undef DEBUG_JM1

    LL += dLL;

#ifdef TIMING

    timer.toc(ESI_E6TH);

#endif
}//close function expectationSingleImage

/** Function to create threads that will work later */
void ProgML2D::createThreads()
{

    //Initialize some variables for using for threads

    th_ids = new pthread_t[threads];
    threads_d = new structThreadTasks[threads];

    barrier_init(&barrier, threads + 1);
    barrier_init(&barrier2, threads + 1);
    barrier_init(&barrier3, threads);//Main thread will not wait here

    for (int i = 0; i < threads; i++)
    {
        threads_d[i].thread_id = i;
        threads_d[i].prm = this;

        int result = pthread_create((th_ids + i), NULL, doThreadsTasks,
                                    (void *) (threads_d + i));

        if (result != 0)
        {
            REPORT_ERROR(ERR_THREADS_NOTINIT, "");
        }
    }

}//close function createThreads

/** Free threads memory and exit */
void ProgML2D::destroyThreads()
{
    threadTask = TH_EXIT;
    barrier_wait(&barrier);
    delete[] th_ids;
    delete[] threads_d;
}

/// Function for threads do different tasks
void * doThreadsTasks(void * data)
{
    structThreadTasks * thread_data = (structThreadTasks *) data;

    ProgML2D * prm = thread_data->prm;

    barrier_t & barrier = prm->barrier;
    barrier_t & barrier2 = prm->barrier2;

    //Loop until the threadTask become TH_EXIT

    do
    {
        //Wait until main threads order to start
        //debug_print("Waiting on barrier, th ", thread_id);
        barrier_wait(&barrier);

        //Check task to do

        switch (prm->threadTask)
        {

        case TH_PFS_REFNO:
            prm->doThreadPreselectFastSignificantRefno();
            break;

        case TH_ESI_REFNO:
            prm->doThreadExpectationSingleImageRefno();
            break;

        case TH_ESI_UPDATE_REFNO:
            prm->doThreadESIUpdateRefno();
            break;

        case TH_RR_REFNO:
            prm->doThreadRotateReferenceRefno();
            break;

        case TH_RRR_REFNO:
            prm->doThreadReverseRotateReferenceRefno();
            break;

        case TH_EXIT:
            pthread_exit(NULL);
            break;

        }

        barrier_wait(&barrier2);

    }
    while (1);

}//close function doThreadsTasks


/// Function to assign refno jobs to threads
/// the starting refno is passed through the out refno parameter
/// and is returned the number of refno's to do, 0 if no more refno's.
int ProgML2D::getThreadRefnoJob(int &refno)
{
    int load = 0;

    pthread_mutex_lock(&refno_mutex);

    if (refno_count < model.n_ref)
    {
        load = XMIPP_MIN(refno_load, model.n_ref - refno_count);
        refno = refno_index;
        refno_index = (refno_index + load) % model.n_ref;
        refno_count += load;
    }

    pthread_mutex_unlock(&refno_mutex);

    return load;
}//close function getThreadRefnoJob

///Function for awake threads for different tasks
void ProgML2D::awakeThreads(ThreadTask task, int start_refno, int load)
{
    threadTask = task;
    refno_index = start_refno;
    refno_count = 0;
    refno_load = load;
    barrier_wait(&barrier);
    //Wait until done
    barrier_wait(&barrier2);
}//close function awakeThreads


void ProgML2D::doThreadRotateReferenceRefno()
{
#ifdef DEBUG
    std::cerr << "entering doThreadRotateReference " << std::endl;
#endif

    double AA, stdAA=0., psi, dum, avg;
    MultidimArray<double> Maux(dim, dim);
    MultidimArray<std::complex<double> > Faux;
    FourierTransformer local_transformer;
    int refnoipsi;

    Maux.setXmippOrigin();

    FOR_ALL_THREAD_REFNO()
    {
        computeStats_within_binary_mask(omask, model.Iref[refno](), dum,
                                        dum, avg, dum);
        for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
        {
            refnoipsi = refno * nr_psi + ipsi;
            // Add arbitrary number (small_angle) to avoid 0-degree rotation (lacking interpolation)
            psi = (double) (ipsi * psi_max / nr_psi) + SMALLANGLE;
            rotate(BSPLINE3, Maux, model.Iref[refno](), -psi, 'Z', WRAP);
            apply_binary_mask(mask, Maux, Maux, avg);
            // Normalize the magnitude of the rotated references to 1st rot of that ref
            // This is necessary because interpolation due to rotation can lead to lower overall Fref
            // This would result in lower probabilities for those rotations
            AA = Maux.sum2();
            if (ipsi == 0)
            {
                stdAA = AA;
                A2[refno] = AA;
            }

            if (AA > 0)
                Maux *= sqrt(stdAA / AA);

            if (fast_mode)
                mref[refnoipsi] = Maux;

            // Do the forward FFT
            local_transformer.FourierTransform(Maux, Faux, false);

            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Faux)
            {
                dAij(Faux, i, j) = conj(dAij(Faux, i, j));
            }

            fref[refnoipsi] = Faux;
        }

        // If we dont use save_mem1 Iref[refno] is useless from here on
        //FIXME: Segmentation fault with blocks
        //if (!save_mem1)
        //    model.Iref[refno]().resize(0, 0);

    }//close while

}//close function doThreadRotateReferenceRefno

void ProgML2D::doThreadReverseRotateReferenceRefno()
{
    double psi, dum, avg;
    MultidimArray<double> Maux(dim, dim), Maux2(dim, dim), Maux3(dim, dim);
    MultidimArray<std::complex<double> > Faux;
    FourierTransformer local_transformer;

    Maux.setXmippOrigin();
    Maux2.setXmippOrigin();
    Maux3.setXmippOrigin();

    FOR_ALL_THREAD_REFNO()
    {
        Maux.initZeros();
        wsum_Mref[refno] = Maux;
        for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
        {
            // Add arbitrary number to avoid 0-degree rotation without interpolation effects
            psi = (double) (ipsi * psi_max / nr_psi) + SMALLANGLE;
            int refnoipsi = refno * nr_psi + ipsi;
            // Do the backward FFT
            // The construction with Faux and Maux3 should perhaps not be necessary,
            // but I am having irreproducible segmentation faults for the iFFT
            Faux = wsumimgs[refnoipsi];

            local_transformer.inverseFourierTransform(Faux, Maux);
            Maux3 = Maux;
            //CenterFFT(Maux3, true);
            centerFFT2(Maux3);
            computeStats_within_binary_mask(omask, Maux3, dum, dum, avg, dum);
            rotate(BSPLINE3, Maux2, Maux3, psi, 'Z', WRAP);
            apply_binary_mask(mask, Maux2, Maux2, avg);
            wsum_Mref[refno] += Maux2;
        }

    }//close while refno

}//close function doThreadReverseRotateReference

///Some macro definitions for the following function
#define IIFLIP (imirror * nr_nomirror_flips + iflip)
#define IROT (IIFLIP * nr_psi + ipsi)
#define IREFMIR ()
#define WEIGHT (dAij(pfs_weight, refno, IROT))
#define MAX_WEIGHT (dAij(pfs_maxweight, imirror, refno))
#define MSIGNIFICANT (dAij(Msignificant, refno, IROT))

void ProgML2D::doThreadPreselectFastSignificantRefno()
{
    MultidimArray<double> Mtrans, Mflip;
    double ropt, aux, diff, pdf, fracpdf;
    double A2_plus_Xi2;
    int irefmir;
    Matrix1D<double> trans(2);
    double local_mindiff;
    double sigma_noise2 = model.sigma_noise * model.sigma_noise;
    int nr_mirror = (do_mirror) ? 2 : 1;

    Mtrans.resize(dim, dim);
    Mtrans.setXmippOrigin();
    Mflip.resize(dim, dim);
    Mflip.setXmippOrigin();

    local_mindiff = 99.e99;

    // A. Translate image and calculate probabilities for every rotation
    FOR_ALL_THREAD_REFNO()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            A2_plus_Xi2 = 0.5 * (A2[refno] + Xi2);
            for (int imirror = 0; imirror < nr_mirror; imirror++)
            {
                irefmir = imirror * model.n_ref + refno;
                // Get optimal offsets
                trans(0) = allref_offsets[2 * irefmir];
                trans(1) = allref_offsets[2 * irefmir + 1];
                ropt = sqrt(trans(0) * trans(0) + trans(1) * trans(1));
                // Do not trust optimal offsets if they are larger than 3*sigma_offset:
                if (ropt > 3 * model.sigma_offset)
                {
                    for (size_t iflip = 0; iflip < nr_nomirror_flips; iflip++)
                        for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
                            MSIGNIFICANT = 1;
                }
                else
                {
                    translate(LINEAR, Mtrans, Mimg, trans, true);
                    for (size_t iflip = 0; iflip < nr_nomirror_flips; iflip++)
                    {
                        applyGeometry(LINEAR, Mflip, Mtrans, F[IIFLIP], IS_INV, WRAP);
                        for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
                        {
                            diff = A2_plus_Xi2;
                            MultidimArray<double> &mref_ref = mref[refno*nr_psi + ipsi];
                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mflip)
                            {
                                diff -= DIRECT_MULTIDIM_ELEM(Mflip, n) * DIRECT_MULTIDIM_ELEM(mref_ref, n);
                            }
                            WEIGHT = diff;
                            if (diff < local_mindiff)
                                local_mindiff = diff;
                        }
                    }
                }//close else if ropt > ...
            }//close for imirror
        }//close if !limit_rot ...
    }//close for_all refno

    ///Update the real mindiff
    pthread_mutex_lock(&update_mutex);
    pfs_mindiff = XMIPP_MIN(pfs_mindiff, local_mindiff);
    refno_index = refno_count = 0;
    pthread_mutex_unlock(&update_mutex);

    ///Wait for all threads update mindiff
    barrier_wait(&barrier3);
    local_mindiff = pfs_mindiff;

    // B. Now that we have local_mindiff, calculate the weights
    FOR_ALL_THREAD_REFNO_NODECL()
    {

        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            for (int imirror = 0; imirror < nr_mirror; imirror++)
            {
                irefmir = imirror * model.n_ref + refno;
                // Get optimal offsets
                trans(0) = allref_offsets[2 * irefmir];
                trans(1) = allref_offsets[2 * irefmir + 1];
                ///Calculate max_weight for this refno-mirror combination
                for (size_t iflip = 0; iflip < nr_nomirror_flips; iflip++)
                    for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
                    {
                        if (!MSIGNIFICANT)
                        {
                            fracpdf = model.alpha_k[refno] *
                                      (imirror ? model.mirror_fraction[refno] : (1.- model.mirror_fraction[refno]));
                            pdf = fracpdf * A2D_ELEM(P_phi, (int)trans(1), (int)trans(0));
                            if (!model.do_student)
                            {
                                // normal distribution
                                aux = (WEIGHT - local_mindiff) / sigma_noise2;
                                // next line because of numerical precision of exp-function
                                WEIGHT = (aux > 1000. ? 0. : exp(-aux) * pdf);
                            }
                            else
                            {
                                // t-student distribution
                                aux = (dfsigma2 + 2. * WEIGHT) / (dfsigma2 + 2. * local_mindiff);
                                WEIGHT = pow(aux, df2) * pdf;
                            }
                            if (WEIGHT > MAX_WEIGHT)
                                MAX_WEIGHT = WEIGHT;
                        }
                    } // close ipsi
                ///Now we have max_weight, set Msignificant values
                for (size_t iflip = 0; iflip < nr_nomirror_flips; iflip++)
                    for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
                        if (!MSIGNIFICANT)
                            MSIGNIFICANT = (WEIGHT >= C_fast * MAX_WEIGHT) ? 1 : 0;
            }//close for imirror
        } //endif limit_rot and pdf_directions
    } //end for_all refno

}//close function doThreadPreselectFastSignificantRefno

void ProgML2D::doThreadExpectationSingleImageRefno()
{
    double diff;
    double aux, pdf, fracpdf, A2_plus_Xi2;
    double weight, stored_weight, weight2, my_maxweight;
    double my_sumweight, my_sumstoredweight, ref_scale = 1.;
    int irot, output_irefmir, refnoipsi, output_refnoipsi;
    //Some local variables to store partial sums of global sums variables
    double local_mindiff, local_wsum_corr, local_wsum_offset, maxw_ref;
    double local_wsum_sc, local_wsum_sc2, local_maxweight, local_maxweight2=0.;
    double sigma_noise2 = model.sigma_noise * model.sigma_noise;
    int local_iopty=0, local_ioptx=0, local_iopt_psi=0, local_iopt_flip=0,
    local_opt_refno=0;

    MultidimArray<double> Maux, Mweight;
    MultidimArray<std::complex<double> > Faux, Fzero(dim, hdim + 1);
    FourierTransformer local_transformer;

    // Setup matrices
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Mweight.resize(sigdim, sigdim);
    Mweight.setXmippOrigin();
    Fzero.initZeros();
    local_transformer.setReal(Maux);
    local_transformer.getFourierAlias(Faux);

    // Start the loop over all refno at old_optrefno (=opt_refno from the previous iteration).
    // This will speed-up things because we will find Pmax probably right away,
    // and this will make the if-statement that checks SIGNIFICANT_WEIGHT_LOW
    // effective right from the start
    //std::cerr << "DEBUG_JM: doThreadExpectationSingleImageRefno: " << std::endl;
    FOR_ALL_THREAD_REFNO()
    {

        int output_refno = mygroup * model.n_ref + refno;
        //        std::cerr << "DEBUG_JM:           refno: " << refno << std::endl;
        //        std::cerr << "DEBUG_JM:     model.n_ref: " <<    model.n_ref << std::endl;
        //        std::cerr << "DEBUG_JM:    output_refno: " <<    output_refno << std::endl;
        refw[output_refno] = refw2[output_refno] = refw_mirror[output_refno] = 0.;
        local_maxweight = -99.e99;
        local_mindiff = 99.e99;
        local_wsum_sc = local_wsum_sc2 = local_wsum_corr = local_wsum_offset = 0;
        // Initialize my weighted sums
        for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
        {
            output_refnoipsi = output_refno * nr_psi + ipsi;
            mysumimgs[output_refnoipsi] = Fzero;
            sumw_refpsi[output_refnoipsi] = 0.;
        }

        // This if is for limited rotation options
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            if (model.do_norm)
                ref_scale = opt_scale / model.scale[refno];

            A2_plus_Xi2 = 0.5 * (ref_scale * ref_scale * A2[refno] + Xi2);

            maxw_ref = -99.e99;
            for (size_t iflip = 0; iflip < nr_flip; iflip++)
            {
                if (iflip == nr_nomirror_flips)
                    maxw_ref = -99.e99;
                for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
                {
                    refnoipsi = refno * nr_psi + ipsi;
                    output_refnoipsi = output_refno * nr_psi + ipsi;
                    irot = iflip * nr_psi + ipsi;
                    output_irefmir  = (int)floor(iflip / nr_nomirror_flips)
                                      * factor_nref * model.n_ref + refno;

                    //#define DEBUG_JM2
#ifdef DEBUG_JM2

                    if (iter > 1)
                    {
                        std::cerr << iter << " iflip, ipsi, refno, irot: " << iflip << " " << ipsi << " " << refno << " " << irot << std::endl;
                        std::cerr << iter << " A2_plus_Xi2: " << A2_plus_Xi2 << std::endl;
                        std::cerr << "dAij(Msignificant): " << dAij(Msignificant, refno, irot) << std::endl;
                    }
#endif
                    // This if is the speed-up caused by the -fast options
                    if (dAij(Msignificant, refno, irot))
                    {
                        if (iflip < nr_nomirror_flips)
                            fracpdf = model.alpha_k[refno] * (1. - model.mirror_fraction[refno]);
                        else
                            fracpdf = model.alpha_k[refno] * model.mirror_fraction[refno];

                        // A. Backward FFT to calculate weights in real-space
                        //Set this references to avoid indexing inside the heavy loop
                        //MultidimArray<std::complex<double> > & Fimg_flip_aux = Fimg_flip[iflip];
                        Faux = Fimg_flip[iflip];
                        MultidimArray<std::complex<double> > & fref_aux = fref[refnoipsi];
                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                        {
                            DIRECT_MULTIDIM_ELEM(Faux, n) *= DIRECT_MULTIDIM_ELEM(fref_aux, n);
                        }
                        // Takes the input from Faux, and leaves the output in Maux
                        local_transformer.inverseFourierTransform();
                        //CenterFFT(Maux, true);
                        centerFFT2(Maux);


#ifdef DEBUG_JM2

                        if (iter > 1)
                            std::cerr
                            << "Maux: " <<  std::endl << Maux
                            << "Fimg_flip[iflip]" <<  std::endl<< Fimg_flip[iflip]
                            << "fref[refnoipsi]" <<  std::endl<< fref[refnoipsi] << std::endl;
#endif

                        // B. Calculate weights for each pixel within sigdim (Mweight)
                        my_sumweight = my_sumstoredweight = my_maxweight = 0.;

                        FOR_ALL_ELEMENTS_IN_ARRAY2D(Mweight)
                        {
                            diff = A2_plus_Xi2 - ref_scale * A2D_ELEM(Maux, i, j) * ddim2;
                            pdf = fracpdf * A2D_ELEM(P_phi, i, j);

#ifdef DEBUG_JM2

                            if (iter >= 2 && current_image == myFirstImg)
                                std::cerr << "---------------------------------------" << std::endl
                                << "   pdf " << pdf << std::endl
                                << "   A2_plus_Xi2 " << A2_plus_Xi2 << std::endl
                                << "   A2D_ELEM(Maux, i, j) " << A2D_ELEM(Maux, i, j) << std::endl
                                << "   ref_scale " << ref_scale << std::endl
                                << "   ddim2 " << ddim2 << std::endl
                                << "diff " << diff << std::endl
                                << "trymindiff " << trymindiff << std::endl
                                << "sigma_noise2 " << sigma_noise2 << std::endl;
#endif

                            if (!model.do_student)
                            {
                                // Normal distribution
                                aux = (diff - trymindiff) / sigma_noise2;
                                // next line because of numerical precision of exp-function
                                weight = (aux > 1000.) ? 0. : exp(-aux) * pdf;
                                //#define DEBUG_JM2
#ifdef DEBUG_JM2

                                if (iter >=2 && current_image == myFirstImg && pdf > 0)
                                    std::cerr << "aux = (diff - trymindiff) / sigma_noise2: " << aux << std::endl
                                    << "weight: " << weight << std::endl;
#endif
#undef DEBUG_JM2
                                // store weight
                                A2D_ELEM(Mweight, i, j) = stored_weight = weight;
                                // calculate weighted sum of (X-A)^2 for sigma_noise update
                                local_wsum_corr += weight * diff;
                            }
                            else
                            {
                                // t-student distribution
                                // pdf = (1 + diff2/sigma2*df)^df2
                                // Correcting for mindiff:
                                // pdfc = (1 + diff2/sigma2*df)^df2 / (1 + mindiff/sigma2*df)^df2
                                //      = ( (1 + diff2/sigma2*df)/(1 + mindiff/sigma2*df) )^df2
                                //      = ( (sigma2*df + diff2) / (sigma2*df + mindiff) )^df2
                                // Extra factor two because we saved 0.5*diff2!!
                                aux = (dfsigma2 + 2. * diff)
                                      / (dfsigma2 + 2. * trymindiff);
                                weight = pow(aux, df2) * pdf;
                                // Calculate extra weight acc. to Eq (10) Wang et al.
                                // Patt. Recognition Lett. 25, 701-710 (2004)
                                weight2 = (df + ddim2) / (df + (2. * diff / sigma_noise2));
                                // Store probability weights
                                stored_weight = weight * weight2;
                                A2D_ELEM(Mweight, i, j) = stored_weight;
                                // calculate weighted sum of (X-A)^2 for sigma_noise update
                                local_wsum_corr += stored_weight * diff;
                                refw2[output_refno] += stored_weight;
                            }

                            local_mindiff = XMIPP_MIN(local_mindiff, diff);

                            // Accumulate sum weights for this (my) matrix
                            my_sumweight += weight;
                            my_sumstoredweight += stored_weight;
                            // calculated weighted sum of offsets as well
                            local_wsum_offset += weight * A2D_ELEM(Mr2, i, j);

                            if (model.do_norm)
                            {
                                // weighted sum of Sum_j ( X_ij*A_kj )
                                local_wsum_sc += stored_weight * (A2_plus_Xi2 - diff) / ref_scale;
                                // weighted sum of Sum_j ( A_kj*A_kj )
                                local_wsum_sc2 += stored_weight * A2[refno];
                            }

                            // keep track of optimal parameters
                            my_maxweight = XMIPP_MAX(my_maxweight, weight);

                            if (weight > local_maxweight)
                            {
                                if (model.do_student)
                                    local_maxweight2 = weight2;
                                local_maxweight = weight;
                                local_iopty = i;
                                local_ioptx = j;
                                local_iopt_psi = ipsi;
                                local_iopt_flip = iflip;
                                local_opt_refno = output_refno;
                            }

                            if (fast_mode && weight > maxw_ref)
                            {
                                maxw_ref = weight;
                                iopty_ref[output_irefmir] = i;
                                ioptx_ref[output_irefmir] = j;
                                ioptflip_ref[output_irefmir] = iflip;
                            }

                        } // close for over all elements in Mweight

#ifdef DEBUG_JM3
                        if (iter == 3 && refno==0)
                        {
                            std::cerr << Mweight << std::endl;
                            std::cerr << "maxweight: " << my_maxweight << " " << my_sumstoredweight << " " << my_sumweight << std::endl;
                            exit(1);
                        }
#endif
                        // C. only for significant settings, store weighted sums
                        if (my_maxweight > SIGNIFICANT_WEIGHT_LOW * maxweight)
                        {
                            sumw_refpsi[output_refno * nr_psi + ipsi] += my_sumstoredweight;

                            if (iflip < nr_nomirror_flips)
                                refw[output_refno] += my_sumweight;
                            else
                                refw_mirror[output_refno] += my_sumweight;

                            // Back from smaller Mweight to original size of Maux
                            Maux.initZeros();

                            FOR_ALL_ELEMENTS_IN_ARRAY2D(Mweight)
                            {
                                A2D_ELEM(Maux, i, j) = A2D_ELEM(Mweight, i, j);
                            }

                            // Use forward FFT in convolution theorem again
                            // Takes the input from Maux and leaves it in Faux
                            local_transformer.FourierTransform();

                            MultidimArray< std::complex<double> > &mysumimgs_ref = mysumimgs[output_refnoipsi];
                            MultidimArray< std::complex<double> > &Fimg_flip_ref = Fimg_flip[iflip];

                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                            {
                                DIRECT_MULTIDIM_ELEM(mysumimgs_ref, n) +=
                                    conj(DIRECT_MULTIDIM_ELEM(Faux,n)) * DIRECT_MULTIDIM_ELEM(Fimg_flip_ref,n);
                            }
                        }
                    } // close if Msignificant
                } // close for ipsi
            } // close for iflip
        } // close if pdf_directions

        pthread_mutex_lock(&update_mutex);
        //Update maxweight
        if (local_maxweight > maxweight)
        {
            maxweight = local_maxweight;

            if (model.do_student)
                maxweight2 = local_maxweight2;
            iopty = local_iopty;
            ioptx = local_ioptx;
            iopt_psi = local_iopt_psi;
            iopt_flip = local_iopt_flip;
            opt_refno = local_opt_refno;
        }

        //Update sums
        sum_refw += refw[output_refno] + refw_mirror[output_refno];
        wsum_offset += local_wsum_offset;
        wsum_corr += local_wsum_corr;



        mindiff = XMIPP_MIN(mindiff, local_mindiff);

        if (model.do_norm)
        {
            wsum_sc += local_wsum_sc;
            wsum_sc2 += local_wsum_sc2;
        }
        pthread_mutex_unlock(&update_mutex);

        //Ask for next job
    } // close while refno

#define DP(x) << std::setw(15) << x
    //#define DEBUG_JM
#ifdef DEBUG_JM
    if (iter > 1)
    {
        std::cerr << "DEBUG_JM: ====== iter: " << iter << "======= block: " << current_block << std::endl;

        std::cerr << formatString("====> img: %lu\n", current_image);
        std::cerr << "Xi2: " << Xi2 << std::endl;
        std::cerr << "sum_refw: " << sum_refw << std::endl;

        std::cerr DP("A2") DP("refw") DP("refw_mirror") << std::endl;
        for (int refno = 0; refno < model.n_ref; ++refno)
            std::cerr DP(A2[refno]) DP(refw[refno]) DP(refw_mirror[refno]) << std::endl;
        //                        std::cerr << "local_wsum_corr: " << local_wsum_corr << std::endl;
        //                        std::cerr << "wsum_corr: " << wsum_corr << std::endl;
        if (iter > 2)
        {
            std::cerr << "DEBUG_JM: ..............EXITING....................: "  << std::endl;

            exit(1);
        }
    }
#endif
#undef DEBUG_JM



}//close function doThreadExpectationSingleImage

void ProgML2D::doThreadESIUpdateRefno()
{

    double scale_dim2_sumw = (opt_scale * ddim2) / sum_refw;
    int num_refs = model.n_ref * factor_nref;

    FOR_ALL_THREAD_REFNO()
    {
        int output_refno = mygroup * model.n_ref + refno;

        if (fast_mode)
        {
            for (int group = 0; group < factor_nref; group++)
            {
                int group_refno = group * model.n_ref + refno;

                // Update optimal offsets for refno (and its mirror)
                allref_offsets[2 * group_refno] = -(double) ioptx_ref[output_refno]
                                                  * MAT_ELEM(F[ioptflip_ref[output_refno]], 0, 0)
                                                  - (double) iopty_ref[output_refno]
                                                  * MAT_ELEM(F[ioptflip_ref[output_refno]], 0, 1);
                allref_offsets[2 * group_refno + 1] = -(double) ioptx_ref[output_refno]
                                                      * MAT_ELEM(F[ioptflip_ref[output_refno]], 1, 0)
                                                      - (double) iopty_ref[output_refno]
                                                      * MAT_ELEM(F[ioptflip_ref[output_refno]], 1, 1);
                if (do_mirror)
                {
                    allref_offsets[2 * (num_refs + group_refno)]
                    = -(double) ioptx_ref[num_refs + output_refno]
                      * MAT_ELEM(F[ioptflip_ref[num_refs + output_refno]], 0, 0)
                      - (double) iopty_ref[num_refs + output_refno]
                      * MAT_ELEM(F[ioptflip_ref[num_refs + output_refno]], 0, 1);
                    allref_offsets[2 * (num_refs + group_refno) + 1]
                    = -(double) ioptx_ref[num_refs + output_refno]
                      * MAT_ELEM(F[ioptflip_ref[num_refs + output_refno]], 1, 0)
                      - (double) iopty_ref[num_refs + output_refno]
                      * MAT_ELEM(F[ioptflip_ref[num_refs + output_refno]], 1, 1);
                }
            }
        }

        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            sumw[output_refno] += (refw[output_refno] + refw_mirror[output_refno]) / sum_refw;
            sumw2[output_refno] += refw2[output_refno] / sum_refw;
            sumw_mirror[output_refno] += refw_mirror[output_refno] / sum_refw;

            if (model.do_student)
            {
                sumwsc[output_refno] += refw2[output_refno] * (opt_scale) / sum_refw;
                sumwsc2[output_refno] += refw2[output_refno] * (opt_scale * opt_scale)
                                         / sum_refw;
            }
            else
            {
                sumwsc[output_refno] += (refw[output_refno] + refw_mirror[output_refno])
                                        * (opt_scale) / sum_refw;
                sumwsc2[output_refno] += (refw[output_refno] + refw_mirror[output_refno])
                                         * (opt_scale * opt_scale) / sum_refw;
            }

            std::complex<double> cscale_dim2_sumw=scale_dim2_sumw;
            for (size_t ipsi = 0; ipsi < nr_psi; ipsi++)
            {
                int refnoipsi = output_refno * nr_psi + ipsi;
                // Correct weighted sum of images for new bgmean (only first element=origin in Fimg)
                if (model.do_norm)
                    dAij(mysumimgs[refnoipsi],0,0) -= sumw_refpsi[refnoipsi] * (bgmean - old_bgmean) / ddim2;
                // Sum mysumimgs to the global weighted sum
                wsumimgs[refnoipsi] += (cscale_dim2_sumw * mysumimgs[refnoipsi]);
            }
        }

    }//close while refno

}//close function doThreadESIUpdateRefno

void ProgML2D::iteration()
{
    for (current_block = 0; current_block < blocks; current_block++)
    {
        // Integrate over all images
        expectation();
        //Maximize the model
        maximization();
    }//close for blocks
}


void ProgML2D::expectation()
{
    MultidimArray<std::complex<double> > Fdzero(dim, hdim + 1);
    int num_output_refs = factor_nref * model.n_ref;

    LOG("      ProgML2D::expectation BEGIN");
#ifdef DEBUG

    std::cerr<<"entering expectation"<<std::endl;
#endif

#ifdef TIMING

    timer.tic(E_RR);
#endif

    rotateReference();
#ifdef TIMING

    timer.toc(E_RR);
    timer.tic(E_PRE);
#endif
    // Pre-calculate pdf of all in-plane transformations
    calculatePdfInplane();

    // Initialize weighted sums
    LL = 0.;
    wsumimgs.clear();

    sumw.assign(num_output_refs, 0.);
    sumw2.assign(num_output_refs, 0.);
    sumwsc.assign(num_output_refs, 0.);
    sumwsc2.assign(num_output_refs, 0.);
    sumw_mirror.assign(num_output_refs, 0.);

    wsum_sigma_noise = 0.;
    wsum_sigma_offset = 0.;
    sumfracweight = 0.;

    Fdzero.initZeros();
    wsumimgs.assign(num_output_refs * nr_psi, Fdzero);

    // Local variables of old threadExpectationSingleImage
    Image<double> img;
    FileName fn_img, fn_tkrans;
    Matrix1D<double> opt_offsets(2);
    double old_phi = -999., old_theta = -999.;
    double opt_flip;

    //Some initializations
    opt_scale = 1., bgmean = 0.;

    Msignificant.resizeNoCopy(model.n_ref, nr_psi * nr_flip);

    static size_t img_done;
    if (current_block == 0) //when not iem current block is always 0
    {
        initProgress(nr_images_local);
        img_done = 0;
    }

    String _msg = formatString("Images: %lu, first: %lu, last: %lu", nr_images_local, myFirstImg, myLastImg);
    LOG(_msg.c_str());
    //std::cerr << "-----xmipp_current: Expectation, iter " << iter << "-------" << std::endl;
    //for (int imgno = 0, img_done = 0; imgno < nn; imgno++)
    // Loop over all images
    FOR_ALL_LOCAL_IMAGES()
    {

        if (IMG_BLOCK(imgno) == current_block)
        {
            current_image = imgno;
            //std::cerr << "\n ======>>> imgno: " << imgno << std::endl;
            mygroup = (factor_nref > 1) ? divide_equally_group(nr_images_global, factor_nref, imgno) : 0;

            MDimg.getValue(MDL_IMAGE, fn_img, img_id[imgno]);
            img.read(fn_img);
            img().setXmippOrigin();
            Xi2 = img().sum2();
            Mimg = img();


            //#define DEBUG_JM1
#ifdef DEBUG_JM1

            //if (iter >= 2 && current_image == myFirstImg)
            printf("   ====================>>> Iter: %02d Image: %06lu: \n", iter, imgno);
            printf("                                     fn_img: %s", fn_img.c_str());
            //printf("                                    mygroup: %d", mygroup);
#endif
#undef DEBUG_JM1

            // These two parameters speed up expectationSingleImage
            opt_refno = imgs_optrefno[IMG_LOCAL_INDEX];
            trymindiff = imgs_trymindiff[IMG_LOCAL_INDEX];

            if (trymindiff < 0.)
                // 90% of Xi2 may be a good idea (factor half because 0.5*diff is calculated)
                trymindiff = trymindiff_factor * 0.5 * Xi2;

            if (model.do_norm)
            {
                bgmean = imgs_bgmean[IMG_LOCAL_INDEX];
                opt_scale = imgs_scale[IMG_LOCAL_INDEX];
            }

            // Get optimal offsets for all references
            if (fast_mode)
            {
                allref_offsets = imgs_offsets[IMG_LOCAL_INDEX];
            }

            // Read optimal orientations from memory
            if (limit_rot)
            {
                old_phi = imgs_oldphi[IMG_LOCAL_INDEX];
                old_theta = imgs_oldtheta[IMG_LOCAL_INDEX];
            }
            // For limited orientational search: preselect relevant directions
            preselectLimitedDirections(old_phi, old_theta);

            // Use a maximum-likelihood target function in real space
            // with complete or reduced-space translational searches (-fast)
            if (fast_mode)
                preselectFastSignificant();
            else
                Msignificant.initConstant(1);

            expectationSingleImage(opt_offsets);

            // Write optimal offsets for all references to disc
            if (fast_mode)
            {
                imgs_offsets[IMG_LOCAL_INDEX] = allref_offsets;
            }

            // Store mindiff for next iteration
            imgs_trymindiff[IMG_LOCAL_INDEX] = trymindiff;

            // Store opt_refno for next iteration
            imgs_optrefno[IMG_LOCAL_INDEX] = opt_refno;

            // Store optimal phi and theta in memory
            if (limit_rot)
            {
                imgs_oldphi[IMG_LOCAL_INDEX] = model.Iref[opt_refno % model.n_ref].rot();
                imgs_oldtheta[IMG_LOCAL_INDEX] = model.Iref[opt_refno % model.n_ref].tilt();
            }

            // Store optimal normalization parameters in memory
            if (model.do_norm)
            {
                imgs_scale[IMG_LOCAL_INDEX] = opt_scale;
                imgs_bgmean[IMG_LOCAL_INDEX] = bgmean;
            }

            // Output docfile
            opt_flip = 0.;
            if (-opt_psi > 360.)
            {
                opt_psi += 360.;
                opt_flip = 1.;
            }

            dAij(docfiledata,IMG_LOCAL_INDEX,0)
            = model.Iref[opt_refno % model.n_ref].rot(); // rot
            dAij(docfiledata,IMG_LOCAL_INDEX,1)
            = model.Iref[opt_refno % model.n_ref].tilt(); // tilt
            dAij(docfiledata,IMG_LOCAL_INDEX,2) = opt_psi + 360.; // psi
            dAij(docfiledata,IMG_LOCAL_INDEX,3) = opt_offsets(0); // Xoff
            dAij(docfiledata,IMG_LOCAL_INDEX,4) = opt_offsets(1); // Yoff
            dAij(docfiledata,IMG_LOCAL_INDEX,5) = (double) (opt_refno + 1); // Ref
            dAij(docfiledata,IMG_LOCAL_INDEX,6) = opt_flip; // Mirror
            dAij(docfiledata,IMG_LOCAL_INDEX,7) = fracweight; // P_max/P_tot
            dAij(docfiledata,IMG_LOCAL_INDEX,8) = dLL; // log-likelihood
            if (model.do_norm)
            {
                dAij(docfiledata,IMG_LOCAL_INDEX,9) = bgmean; // background mean
                dAij(docfiledata,IMG_LOCAL_INDEX,10) = opt_scale; // image scale
            }
            if (model.do_student)
            {
                dAij(docfiledata,IMG_LOCAL_INDEX,11) = maxweight2; // Robustness weight
            }

            //Report progress and increment the images done
            setProgress(++img_done);

            //#define DEBUG_JM1
#ifdef DEBUG_JM1
            //            {
            //              //std::cerr << "---------------------- DEBUG_JM: current_image: " << current_image << std::endl;
            std::cerr << "                              LL: " << LL << std::endl;
            std::cerr << "                wsum_sigma_noise: " << wsum_sigma_noise << std::endl;
            std::cerr << "               wsum_sigma_offset: " << wsum_sigma_offset << std::endl;
            std::cerr << "                   sumfracweight: " << sumfracweight << std::endl;
            //            }
#endif
#undef DEBUG_JM1

        }//close if current_block, also close of for all images
    }//close for all images

    if (current_block == (blocks - 1))
        endProgress();

    //Changes temporally the model n_ref for the
    //refno loop, but not yet n_ref because in iem
    //isn't yet the end of iteration
    model.n_ref *= factor_nref;
    // Rotate back and calculate weighted sums
    reverseRotateReference();
    //Restore back the model.n_ref
    model.n_ref /= factor_nref;
    LOG("      ProgML2D::expectation END");

}//function expectation


// Update all model parameters
void ProgML2D::maximizeModel(ModelML2D &local_model)
{

#ifdef DEBUG
    std::cerr<<"entering maximization"<<std::endl;
#endif

    Matrix1D<double> rmean_sigma2, rmean_signal2;
    Matrix1D<int> center(2), radial_count;
    MultidimArray<std::complex<double> > Faux, Faux2;
    MultidimArray<double> Maux;
    FileName fn_tmp;

    // After iteration 0, factor_nref will ALWAYS be one
    if (factor_nref > 1)
    {
        int old_refno = local_model.n_ref;
        local_model.setNRef(local_model.n_ref * factor_nref);
        // Now also expand the Iref vector to contains factor_nref times more images
        // Make sure that headers are the same by using = assignments
        for (int group = 1; group < factor_nref; group++)
        {
            for (int refno = 0; refno < old_refno; refno++)
            {

                local_model.Iref[group * old_refno + refno] = local_model.Iref[refno];
            }
        }
    }

    // Update the reference images
    local_model.sumw_allrefs = 0.;
    local_model.sumw_allrefs2 = 0.;


    //#define ASSIGN(var) if (var > SIGNIFICANT_WEIGHT_LOW) local_model.var = var
#define ASSIGN(var) local_model.var = var

    ASSIGN(dim);
    ASSIGN(sumfracweight);
    local_model.LL = LL;
    // Update sigma of the origin offsets
    if (!fix_sigma_offset)
        ASSIGN(wsum_sigma_offset);
    // Update the noise parameters
    if (!fix_sigma_noise)
        ASSIGN(wsum_sigma_noise);


    for (int refno = 0; refno < local_model.n_ref; refno++)
    {
        double weight = sumw[refno];
        if (weight > 0.)
        {
            if (local_model.do_student)
            {
                weight = sumw2[refno];
                local_model.sumw_allrefs2 += sumw2[refno];
            }

            double avg, stddev, min, max;
            wsum_Mref[refno].computeStats(avg, stddev, min, max);
            //std::cerr << formatString("DEBUG_JM: avg %f, stddev %f, min %f, max %f", avg, stddev, min, max) <<std::endl;
            local_model.WsumMref[refno]() = wsum_Mref[refno];
            //local_model.Iref[refno]() *= 1/weight;//fixme: sumwsc2[refno];
            local_model.sumw_allrefs += sumw[refno];
            local_model.WsumMref[refno].setWeight(weight);
        }
        else
        {
            local_model.WsumMref[refno].setWeight(0.);
            local_model.WsumMref[refno]().initZeros(dim, dim);
            local_model.WsumMref[refno]().setXmippOrigin();
        }

        // Adjust average scale (nr_classes will be smaller than model.n_ref for the 3D case!)
        if (local_model.do_norm)
            ASSIGN(sumwsc[refno]);

        if (!fix_fractions)
            ASSIGN(sumw_mirror[refno]);
    }

#ifdef DEBUG

    std::cerr<<"leaving maximization"<<std::endl;

#endif
}//close function maximization

void ProgML2D::maximization()
{
    LOG("   ProgML2D::maximization BEGIN");
    if (blocks == 1) //ie not IEM, normal maximization
    {
        maximizeModel(model);
        model.update();
        // After iteration 0, factor_nref will ALWAYS be one
        factor_nref = 1;
    }
    else //do IEM
    {
        bool special_first = (!do_restart && iter == istart);
        ModelML2D block_model(model.n_ref);

        LOG("      ProgML2D::maximization:   readModel and substractModel");
        if (!special_first)
        {
            readModel(block_model, current_block);
            model.substractModel(block_model);
        }
        LOG("      ProgML2D::maximization:   maximizeModel");
        maximizeModel(block_model);
        LOG("      ProgML2D::maximization:   maximizeModel");
        writeOutputFiles(block_model, OUT_BLOCK);
        LOG("      ProgML2D::maximization:   addModels");
        if (!special_first)
        {
            model.addModel(block_model);
            model.update();
        }
        else if (current_block == blocks - 1) //last block
        {
            for (current_block = 0; current_block < blocks; current_block++)
            {
                readModel(block_model, current_block);

                if (current_block == 0)
                    model = block_model;
                else
                    model.addModel(block_model);
            }
            model.update();
            //printModel("GLOBAL model: after addition and update: ", model);
            // After iteration 0, factor_nref will ALWAYS be one
            factor_nref = 1;
            //restore the value of current block
            --current_block;
        }
    }

    if (model.do_norm)
        correctScaleAverage();

    LOG("      ProgML2D::maximization END");

    //    static int times = 0;
    //    if (times > 1)
    //      exit(1);
    //    ++times;
}//close function maximizationBlocks

void ProgML2D::correctScaleAverage()
{
    //FIXME: This function needs re-implementation for the do_norm parameter
    //now the ref3d comes in metadata, refs_per_class can be avoided
    //    int iclass, nr_classes = ROUND(model.n_ref / refs_per_class);
    //    std::vector<double> wsum_scale(nr_classes), sumw_scale(nr_classes);
    //    ldiv_t temp;
    //    average_scale = 0.;
    //
    //    for (int refno = 0; refno < model.n_ref; refno++)
    //    {
    //        average_scale += model.getSumwsc(refno);
    //        temp = ldiv(refno, refs_per_class);
    //        iclass = ROUND(temp.quot);
    //        wsum_scale[iclass] += model.getSumwsc(refno);
    //        sumw_scale[iclass] += model.getSumw(refno);
    //    }
    //    for (int refno = 0; refno < model.n_ref; refno++)
    //    {
    //        temp = ldiv(refno, refs_per_class);
    //        iclass = ROUND(temp.quot);
    //        if (sumw_scale[iclass] > 0.)
    //        {
    //            model.scale[refno] = wsum_scale[iclass] / sumw_scale[iclass];
    //            model.WsumMref[refno]() *= model.scale[refno];
    //        }
    //        else
    //        {
    //            model.scale[refno] = 1.;
    //        }
    //    }
    //    average_scale /= model.sumw_allrefs;
}//close function correctScaleAverage

/// Add docfiledata to docfile
void ProgML2D::addPartialDocfileData(const MultidimArray<double> &data,
                                     size_t first, size_t last)
{
#ifdef DEBUG
    std::cerr << "Entering addPartialDocfileData" <<std::endl;
#endif

    for (size_t imgno = first; imgno <= last; imgno++)
    {
        size_t index = imgno - first;
        size_t id = img_id[imgno];
        if (do_ML3D)
        {
            MDimg.setValue(MDL_ANGLE_ROT, dAij(data, index, 0), id);
            MDimg.setValue(MDL_ANGLE_TILT, dAij(data, index, 1), id);
        }
        //Here we change the sign of the angle becase in the code
        //is represent the rotation of the reference and we want
        //to store the rotation of the image
        double psi = -dAij(data, index, 2);
        MDimg.setValue(MDL_ANGLE_PSI, psi, id);
        MDimg.setValue(MDL_SHIFT_X, dAij(data, index, 3), id);
        MDimg.setValue(MDL_SHIFT_Y, dAij(data, index, 4), id);
        MDimg.setValue(MDL_REF, (int)round(dAij(data, index, 5)), id);
        if (do_mirror)
        {
            MDimg.setValue(MDL_FLIP, dAij(data, index, 6) != 0., id);
        }
        MDimg.setValue(MDL_PMAX, dAij(data, index, 7), id);
        MDimg.setValue(MDL_LL, dAij(data, index, 8), id);
        if (model.do_norm)
        {
            MDimg.setValue(MDL_BGMEAN, dAij(data, index, 9), id);
            MDimg.setValue(MDL_INTSCALE, dAij(data, index, 10), id);
        }
        if (model.do_student)
        {
            MDimg.setValue(MDL_WROBUST, dAij(data, index, 11), id);
        }
    }

#ifdef DEBUG
    std::cerr << "Leaving addPartialDocfileData" <<std::endl;
#endif
}//close function addDocfileData

//Some macros
#define ITER_PREFIX "iter"//formatString("iter%06d", iter)
#define FINAL_PREFIX "final"
#define IS_FINAL outputType == OUT_FINAL

void ProgML2D::writeOutputFiles(const ModelML2D &model, OutputType outputType)
{
    FileName fn_tmp, fn_prefix, fn_base;
    Image<double> Itmp;
    MetaData MDo;
    bool write_img_xmd = true, write_refs_log = true, write_conv = !do_ML3D;
    bool write_norm = model.do_norm;

    if (iter == 0)
        write_conv = false;

    //By default write down the parameters
    std::vector<Image<double> > const * ptrImages = &(model.Iref);
    std::vector<double> const * ptrMirror = &(model.mirror_fraction);
    double avePmax = model.avePmax;
    double sigma_noise = model.sigma_noise;
    double sigma_offset = model.sigma_offset;

    switch (outputType)
    {
    case OUT_BLOCK:
        write_img_xmd = false;
        write_conv = false;
        // Sjors 18may2010: why not write scales in blocks? We need this now, don't we?
        //write_norm = false;
        fn_prefix = formatString("block%03d", current_block + 1);
        // For blocks we will write out the summations instead of parameters
        ptrImages = &(model.WsumMref);
        ptrMirror = &(model.sumw_mirror);
        avePmax = model.sumfracweight;
        sigma_noise = model.wsum_sigma_noise;
        sigma_offset = model.wsum_sigma_offset;
        break;
    case OUT_ITER:
        //fn_base = getBaseName("_it", iter);
        fn_prefix = ITER_PREFIX;
        break;
    case OUT_FINAL:
        //fn_base = fn_root;
        fn_prefix = FINAL_PREFIX;
        break;
    case OUT_IMGS:
        //std::cerr << "OUT_IMGS" <<std::endl;
        write_refs_log = false;
        //fn_base = getBaseName("_it", iter);
        fn_prefix = ITER_PREFIX;
        break;
    case OUT_REFS:
        //std::cerr << "OUT_REFS" <<std::endl;
        LOG("OUT_REFS");
        write_img_xmd = false;
        write_conv = false;
        //fn_base = getBaseName("_it", iter);
        fn_prefix = ITER_PREFIX;
        break;
    }

    fn_base = (fn_prefix == ITER_PREFIX) ? //All intermediate iteration files should go to "extra" folder
              getIterExtraPath(fn_root, iter) : fn_root;

    if (write_img_xmd)
    {
        //static WriteModeMetaData mode = MD_OVERWRITE;
        //Write image metadata, for each iteration a new block will be written
        fn_tmp = FN_IMAGES_MD(fn_base);
        MDimg.write(fn_tmp);
        //        if (fn_prefix == ITER_PREFIX)
        //            fn_tmp = formatString("iter%06d@%s_iter_images.xmd", iter, rootStr);
        //        else
        //            fn_tmp = formatString("%s_%s_images.xmd", rootStr, prefixStr);
        //MDimg.write(fn_tmp, mode);
        //mode = MD_APPEND;
    }


    if (write_refs_log)
    {
        // Write out current reference images and fill sel & log-file
        // Re-use the MDref metadata that was read in produceSideInfo2
        // This way. MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_REF etc are treated ok for do_ML3D
        int refno = 0;
        size_t objId;
        int select_img;
        double weight;
        size_t count;

        fn_tmp = FN_CLASSES_STK(fn_base);
        //
        //        if (fn_prefix == ITER_PREFIX)
        //            fn_tmp = formatString("%s_iter%06d_refs.stk", rootStr, iter);
        //        else
        //            fn_tmp = formatString("%s_%s_refs.stk", rootStr, prefixStr);

        FOR_ALL_OBJECTS_IN_METADATA(MDref)
        {
            objId = __iter.objId;
            select_img = refno + 1;
            //Itmp = model.Iref[refno];
            //std::cerr << "DEBUG_JM: refno: " << refno << std::endl;
            Itmp = (*ptrImages)[refno];
            //std::cerr << "DEBUG_JM: fn_tmp: " << fn_tmp << std::endl;
            //std::cerr << "DEBUG_JM: select_img: " << select_img << std::endl;
            Itmp.write(fn_tmp, select_img, true, WRITE_REPLACE);
            //MDref.setValue(MDL_ITER, iter, objId);//write out iteration number
            MDref.setValue(MDL_REF, select_img, objId); //Also write reference number
            MDref.setValue(MDL_IMAGE, formatString("%06d@%s", select_img, fn_tmp.c_str()), objId);
            MDref.setValue(MDL_ENABLED, 1, objId);
            weight = model.WsumMref[refno].weight();
            MDref.setValue(MDL_WEIGHT, weight, objId);
            count = ROUND(weight);
            MDref.setValue(MDL_CLASS_COUNT, count, objId);
            if (do_mirror)
                MDref.setValue(MDL_MIRRORFRAC, (*ptrMirror)[refno], objId);
            if (write_conv)
                MDref.setValue(MDL_SIGNALCHANGE, conv[refno]*1000, objId);
            if (write_norm)
                MDref.setValue(MDL_INTSCALE, model.scale[refno], objId);
            ++refno;
        }
        //fn_tmp.copyFile(formatString("%s_output_block%d.stk", fn_tmp.c_str(), current_block));

        //        fn_tmp = formatString("%s_%s", rootStr, prefixStr);
        //        FileName fn_ref = fn_tmp + "_refs.xmd";
        //
        //        if (!fn_prefix.contains("block"))
        //            fn_ref = formatString("iter%06d@%s", iter, fn_ref.c_str());
        //        MDref.write(fn_ref, mode);
        fn_tmp = FN_CLASSES_MD(fn_base);
        MDref.write(fn_tmp);

        if (outputType == OUT_REFS)
          outRefsMd = fn_tmp;

        // Write out log-file
        MetaData mdLog;
        objId = mdLog.addObject();
        mdLog.setValue(MDL_ITER, iter, objId);
        mdLog.setValue(MDL_LL, model.LL, objId);
        mdLog.setValue(MDL_PMAX, avePmax, objId);
        mdLog.setValue(MDL_SIGMANOISE, sigma_noise, objId);
        mdLog.setValue(MDL_SIGMAOFFSET, sigma_offset, objId);
        mdLog.setValue(MDL_RANDOMSEED, seed, objId);
        if (write_norm)
            mdLog.setValue(MDL_INTSCALE, average_scale, objId);

        mdLog.write(FN_LOGS_MD(fn_base), MD_APPEND);

        if (write_img_xmd)
        {
            MetaData mdImgs;
            size_t n = MDref.size();
            for (size_t ref = 1; ref <= n; ++ref)
            {
                mdImgs.importObjects(MDimg, MDValueEQ(MDL_REF, (int)ref));
                mdImgs.write(FN_CLASS_IMAGES_MD(fn_base, ref), MD_APPEND);
            }
        }

        //        if (fn_prefix.contains("block") || mode == MD_OVERWRITE)
        //            mdLog.write(fn_tmp);
        //        else
        //            mdLog.append(fn_tmp);
        //        mode = MD_APPEND;
    }

}//close function writeModel

void ProgML2D::readModel(ModelML2D &model, int block)
{

    // First read general model parameters from _log.xmd
    FileName fn_base = formatString("%s_block%03d", fn_root.c_str(), (block + 1));
    MetaData MDi;
    MDi.read(fn_base + "_logs.xmd");
    model.dim = dim;
    size_t id = MDi.firstObject();
    MDi.getValue(MDL_LL, model.LL, id);
    MDi.getValue(MDL_PMAX, model.sumfracweight, id);
    MDi.getValue(MDL_SIGMANOISE, model.wsum_sigma_noise, id);
    MDi.getValue(MDL_SIGMAOFFSET, model.wsum_sigma_offset, id);

    //Then read reference model parameters from _ref.xmd
    FileName fn_img;
    Image<double> img;
    MDi.clear();
    MDi.read(fn_base + "_refs.xmd");
    int refno = 0;
    double weight = 0;

    model.sumw_allrefs = 0.;
    FOR_ALL_OBJECTS_IN_METADATA(MDi)
    {
        id = __iter.objId;
        MDi.getValue(MDL_IMAGE, fn_img, id);
        img.read(fn_img);
        img().setXmippOrigin();
        model.WsumMref[refno] = img;
        MDi.getValue(MDL_WEIGHT, weight, id);
        model.WsumMref[refno].setWeight(weight);
        model.sumw_allrefs += weight;
        model.sumw_mirror[refno] = 0.;
        if (do_mirror)
            MDi.getValue(MDL_MIRRORFRAC, model.sumw_mirror[refno], id);
        refno++;
    }
}//close function readModel

FileName ProgML2D::getBaseName(String suffix, int number)
{
    FileName fn_base = fn_root + suffix;
    if (number >= 0)
        fn_base.compose(fn_base, number, "");
    return fn_base;
}

