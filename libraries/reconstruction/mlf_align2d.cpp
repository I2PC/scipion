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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#include "mlf_align2d.h"

// Constructor ===============================================
ProgMLF2D::ProgMLF2D(int nr_vols, int rank, int size)
{
    defaultRoot = "mlf2d";
    if (nr_vols == 0)
    {
        do_ML3D = false;
        this->nr_vols = 1;
        refs_per_class = 1;
    }
    else
        do_ML3D = true;

    this->rank = rank;
    this->size = size;
}

void  ProgMLF2D::defineBasicParams(XmippProgram * prog)
{
    ML2DBaseProgram::defineBasicParams(prog);
    prog->addParamsLine("[--no_ctf]          : do not use any CTF correction, you should provide sampling rate");
    prog->addParamsLine("     requires --sampling_rate;");
    //prog->addParamsLine("                    : by defaut the CTF info is read from input images metadata");
    prog->addParamsLine("   [--sampling_rate <Ts>]    : If provided, this sampling rate will overwrites the one in the ctf file");
}

void ProgMLF2D::defineAdditionalParams(XmippProgram * prog, const char * sectionLine)
{
    ML2DBaseProgram::defineAdditionalParams(prog, sectionLine);
    //even more additional params
    prog->addParamsLine(" [ --search_shift <int=3>]      : Limited translational searches (in pixels) ");
    prog->addParamsLine(" [ --reduce_snr <factor=1> ]    : Use a value smaller than one to decrease the estimated SSNRs ");
    prog->addParamsLine(" [ --not_phase_flipped ]        : Use this if the experimental images have not been phase flipped ");
    prog->addParamsLine(" [ --ctf_affected_refs ]        : Use this if the references (-ref) are not CTF-deconvoluted ");
    prog->addParamsLine(" [ --limit_resolution <first_high=0> <high=0> <low=999>]: Exclude frequencies from P-calculations (in Ang)");
    prog->addParamsLine("                               : First value is highest frequency during first iteration.");
    prog->addParamsLine("                               : Second is the highest in following iterations and third is lowest");
    prog->addParamsLine(" [ --fix_high <float=-1>]       : ");
    prog->addParamsLine(" [ --include_allfreqs ] ");
}

// Fourier mode usage ==============================================================
void ProgMLF2D::defineParams()
{
    //add usage

    //params
    defaultRoot = "mlf2d";
    allowRestart = true;
    allowIEM = false;
    defineBasicParams(this);
    defineAdditionalParams(this, "==+ Additional options ==");
    //hidden params
    defineHiddenParams(this);
    addParamsLine(" [--var_psi]");
    addParamsLine(" [--var_trans]");
    addParamsLine(" [--kstest]");
    addParamsLine(" [--iter_histogram <int=-1>]");
}

// Read arguments ==========================================================
void ProgMLF2D::readParams()
{

    // Generate new command line for restart procedure
    cline = "";
    int argc2 = 0;
    char ** argv2 = NULL;

    double restart_offset;
    FileName restart_imgmd, restart_refmd;
    int restart_iter, restart_seed;

    if (checkParam("--restart"))
    {
        do_restart = true;
        MetaData MDrestart;
        char *copy  = NULL;

        MDrestart.read(getParameter(argc, argv, "--restart"));
        cline = MDrestart.getComment();
        size_t id = MDrestart.firstObject();
        MDrestart.getValue(MDL_SIGMAOFFSET, restart_offset, id);
        MDrestart.getValue(MDL_IMGMD, restart_imgmd, id);
        MDrestart.getValue(MDL_REFMD, restart_refmd, id);
        MDrestart.getValue(MDL_ITER, restart_iter, id);
        MDrestart.getValue(MDL_RANDOMSEED, restart_seed, id);
        //MDrestart.getValue(MDL_SIGMANOISE, restart_noise);
        generateCommandLine(cline, argc2, argv2, copy);
    }
    else
    {
        // no restart, just copy argc to argc2 and argv to argv2
        do_restart = false;
        argc2 = argc;
        argv2 = (char **)argv;
        for (int i = 1; i < argc2; i++)
        {
            cline = cline + (String)argv2[i] + " ";
        }
    }
    // Number of threads
    threads = getIntParam("--thr");

    // Main parameters
    model.n_ref = getIntParam("--nref");
    fn_ref = getParam("--ref");
    fn_img = getParam("-i");
    bool a = checkParam("--no_ctf");
    do_ctf_correction = !a;//checkParam("--no_ctf");

    sampling = checkParam("--sampling_rate") ? getDoubleParam("--sampling_rate") : -1;
    fn_root = getParam("--oroot");
    search_shift = getIntParam("--search_shift");
    psi_step = getDoubleParam("--psi_step");
    do_mirror = checkParam("--mirror");
    ini_highres_limit = getIntParam("--limit_resolution", 0);
    highres_limit = getIntParam("--limit_resolution", 1);
    lowres_limit = getIntParam("--limit_resolution", 2);
    phase_flipped = !checkParam("--not_phase_flipped");
    reduce_snr = getDoubleParam("--reduce_snr");
    first_iter_noctf = checkParam("--ctf_affected_refs");

    // Less common stuff
    Niter = getIntParam("--iter");
    istart = do_ML3D ? 1 : getIntParam("--restart");
    sigma_offset = getDoubleParam("--offset");
    eps = getDoubleParam("--eps");
    fn_frac = getParam("--frac");
    fix_fractions = checkParam("--fix_fractions");
    fix_sigma_offset = checkParam("--fix_sigma_offset");
    fix_sigma_noise = checkParam("--fix_sigma_noise");
    C_fast = getDoubleParam("-C");
    //fn_doc = getParam("--doc");
    do_include_allfreqs = checkParam("--include_allfreqs");
    fix_high = getDoubleParam("--fix_high");

    search_rot = getDoubleParam("--search_rot");

    // Hidden arguments
    debug = getIntParam("--debug");
    do_variable_psi = checkParam("--var_psi");
    do_variable_trans = checkParam("--var_trans");
    do_norm = checkParam("--norm");
    do_student = checkParam("--student");
    df = getDoubleParam("--student");
    do_student_sigma_trick = !checkParam("--no_sigma_trick");
    do_kstest = checkParam("--kstest");
    iter_write_histograms = getIntParam("--iter_histogram");
    seed = getIntParam("--random_seed");

    // Now reset some stuff for restart
    if (do_restart)
    {
        fn_img = restart_imgmd;
        fn_ref = restart_refmd;
        model.n_ref = 0; // Just to be sure (not strictly necessary)
        sigma_offset = restart_offset;
        //sigma_noise = restart_noise;
        seed = restart_seed;
        istart = restart_iter + 1;
    }

    if (seed == -1)
        seed = time(NULL);

}

// Show ====================================================================
void ProgMLF2D::show(bool ML3D)
{

    if (verbose)
    {
        // To screen
        if (!do_ML3D)
        {
            std::cout
            << " -----------------------------------------------------------------" << std::endl
            << " | Read more about this program in the following publication:    |" << std::endl
            << " |  Scheres ea. (2007) Structure, 15, 1167-1177                  |" << std::endl
            << " |                                                               |" << std::endl
            << " |   *** Please cite it if this program is of use to you! ***    |" << std::endl
            << " -----------------------------------------------------------------" << std::endl;
        }
        std::cout
        << "--> Multi-reference refinement " << std::endl
        << "--> using a maximum-likelihood in Fourier-space (MLF) target " <<std::endl
        << (do_ctf_correction ? "--> with CTF correction " : "--> ignoring CTF effects ")<<std::endl
        << "  Input images            : " << fn_img << " (" << nr_images_global << ")" << std::endl;

        if (!fn_ref.empty())
            std::cout << "  Reference image(s)      : " << fn_ref << std::endl;
        else
            std::cout << "  Number of references:   : " << model.n_ref << std::endl;

        std::cout
        << "  Output rootname         : " << fn_root << std::endl
        << "  Stopping criterium      : " << eps << std::endl
        << "  initial sigma offset    : " << sigma_offset << std::endl
        << "  Psi sampling interval   : " << psi_step << " degrees" << std::endl
        << "  Translational searches  : " << search_shift << " pixels" << std::endl
        << "  Low resolution limit    : " << lowres_limit << " Ang" << std::endl
        << "  High resolution limit   : " << highres_limit << " Ang" << std::endl;

        if (reduce_snr != 1.)
            std::cout << "  Multiply estimated SNR  : " << reduce_snr << std::endl;
        if (reduce_snr > 1.)
            std::cerr << "  --> WARNING!! With reduce_snr>1 you may likely overfit the noise!" << std::endl;
        std::cout << "  Check mirrors           : " << (do_mirror ? "true" : "false") << std::endl;
        if (!fn_frac.empty())
            std::cout << "  Initial model fractions : " << fn_frac << std::endl;
        if (do_ctf_correction)
        {
            std::cout << "    + Assuming images have " << (phase_flipped ? "" : "not") << "been phase flipped " << std::endl;

            FOR_ALL_DEFOCUS_GROUPS()
            {
                std::cout << formatString("    + CTF group %d contains %d images", ifocus + 1, count_defocus[ifocus]) << std::endl;
            }
        }
        if (ini_highres_limit > 0.)
            std::cout << "    + High resolution limit for 1st iteration set to " << ini_highres_limit << "Ang"<<std::endl;
        if (search_rot < 180.)
            std::cout << "    + Limit orientational search to +/- " << search_rot << " degrees" << std::endl;
        if (do_variable_psi)
            std::cout << "    + Vary in-plane rotational sampling with resolution " << std::endl;
        if (do_variable_trans)
            std::cout << "    + Vary in-plane translational sampling with resolution " << std::endl;

        // Hidden stuff
        if (fix_fractions)
        {
            std::cout << "    + Do not update estimates of model fractions." << std::endl;
        }
        if (fix_sigma_offset)
        {
            std::cout << "    + Do not update sigma-estimate of origin offsets." << std::endl;
        }
        if (fix_sigma_noise)
        {
            std::cout << "    + Do not update estimated noise spectra." << std::endl;
        }
        if (do_student)
        {
            std::cout << "  -> Use t-student distribution with df = " <<df<< std::endl;
            if (do_student_sigma_trick)
            {
                std::cout << "  -> Use sigma-trick for t-student distributions" << std::endl;
            }
        }
        if (do_norm)
        {
            std::cout << "  -> Developmental: refine normalization internally "<<std::endl;
        }
        if (do_kstest)
        {
            std::cout << "  -> Developmental: perform KS-test on noise distributions "<<std::endl;
            if (iter_write_histograms >0.)
                std::cout << "  -> Developmental: write noise histograms at iteration "<<iter_write_histograms<<std::endl;
        }
        std::cout << " -----------------------------------------------------------------" << std::endl;

    }

}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produceSideInfo2)
void ProgMLF2D::produceSideInfo()
{
    //LOG_LEVEL(produceSideInfo);
    LOG_FUNCTION();

    FileName                    fn_tmp, fn_base, fn_tmp2;
    Image<double>              img;
    CTFDescription                    ctf;
    MultidimArray<double>           dum, rmean_ctf;
    Matrix2D<double>           A(3, 3);
    Matrix1D<double>           offsets(2);
    MultidimArray<double>      Maux, Maux2;
    MultidimArray<std::complex<double> >  Faux, ctfmask; //2D
    MultidimArray<int>        radial_count; //1D
    Matrix1D<int>               center(2);
    std::vector<int>            tmppointp, tmppointp_nolow, tmppointi, tmppointj;
    double                      Q0=0.;

    // Read selfile with experimental images
    MDimg.read(fn_img);
    MDimg.removeDisabled(); // Discard disabled images
    nr_images_global = MDimg.size();

    // Create a vector of objectIDs, which may be randomized later on
    MDimg.findObjects(img_id);

    // Get image sizes and total number of images
    size_t idum, ndum;
    getImageSize(MDimg, dim, idum, idum, ndum);
    hdim = dim / 2;
    dim2 = dim * dim;

    // Set sampling stuff: flipping matrices, psi_step etc.
    initSamplingStuff();
    max_nr_psi = nr_psi;

    // Initialization of resolhist
    if (do_kstest)
    {
        resolhist.resize(hdim);
        for (size_t ires = 0; ires < hdim; ires++)
        {
            resolhist[ires].init(HISTMIN, HISTMAX, HISTSTEPS);
        }
    }

    // Fill limited translation shift-vectors
    nr_trans = 0;
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    int ss = search_shift * search_shift;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(Maux)
    {
        int r2 = i * i + j * j;
        if (r2 <= ss)
        {
            XX(offsets) = (double)j;
            YY(offsets) = (double)i;
            Vtrans.push_back(offsets);
            if (i == 0 && j == 0)
                zero_trans = nr_trans;
            nr_trans++;
        }
    }

    FileName fnt_img, fnt;
    std::vector<FileName> all_fn_ctfs;

    count_defocus.clear();
    Vctf.clear();
    Vdec.clear();
    Vsig.clear();
    if (!do_ctf_correction)
    {
        nr_focus = 1;
        count_defocus.push_back(nr_images_global);
        Vctf.push_back(dum);
        Vctf[0].resize(hdim);
        Vctf[0].initConstant(1.); //fill with 1.
        Vdec.push_back(Vctf[0]); //just copy
        Vsig.push_back(Vctf[0]);
    }
    else
    {
        // Read ctfdat and determine the number of CTF groups
        nr_focus = 0;
        all_fn_ctfs.clear();
        FileName fn_ctf;

        //CTF info now comes on input metadatas
        MetaData mdCTF;
        groupCTFMetaData(MDimg, mdCTF);

        if (sampling > 0)
          mdCTF.setValueCol(MDL_CTF_SAMPLING_RATE, sampling);

        nr_focus = mdCTF.size();

        //todo: set defocus group for each image

        // Check the number of images in each CTF group
        // and read CTF-parameters from disc
        //FOR_ALL_DEFOCUS_GROUPS()
        int ifocus = 0;
        count_defocus.resize(nr_focus);
        size_t id;
        MultidimArray<double> dum(hdim);

        FOR_ALL_OBJECTS_IN_METADATA(mdCTF)
        {
            id = __iter.objId;
            //Set defocus group
            mdCTF.setValue(MDL_DEFGROUP, ifocus, id);
            //Read number of images in group
            size_t c;
            mdCTF.getValue(MDL_COUNT, c, id);
            count_defocus[ifocus] = (int)c;
            if (count_defocus[ifocus] < 50 && verbose)
                std::cerr << "WARNING%% CTF group " << (ifocus + 1) << " contains less than 50 images!" << std::endl;

            //Read ctf from disk
            ctf.readFromMetadataRow(mdCTF, id);

            double astigmCTFFactor = fabs( (ctf.DeltafV - ctf.DeltafU) / (std::max(ctf.DeltafV, ctf.DeltafU)) );
            // we discard the CTF with a normalized diference between deltaU and deltaV of 10%
            if (astigmCTFFactor > 0.1)
            {
                FileName ctfname;
                //REPORT_ERROR(ERR_NUMERICAL, "Prog_MLFalign2D-ERROR%% Only non-astigmatic CTFs are allowed!");
                std::cerr << "CTF is too astigmatic. We ill ignore it." << std::endl;
                std::cerr << "   defocus U: " << ctf.DeltafU << "   defocus V: " << ctf.DeltafV << std::endl;
                std::cerr << "   astigmCTFFactor " << astigmCTFFactor <<  std::endl;
                //continue;
            }
            ctf.DeltafV =  (ctf.DeltafU+ctf.DeltafV)*0.5;
            ctf.DeltafU =  ctf.DeltafV;

            ctf.K = 1.;
            ctf.enable_CTF = true;
            ctf.produceSideInfo();
            ctf.generateCTF(dim, dim, ctfmask);
            if (ifocus == 0)
            {
                sampling = ctf.Tm;
                Q0 = ctf.Q0;
            }
            else
            {
                if (sampling != ctf.Tm)
                    REPORT_ERROR(ERR_NUMERICAL, "Prog_MLFalign2D-ERROR%% Different sampling rates in CTF parameter files!");
                if (Q0 != ctf.Q0 )
                    REPORT_ERROR(ERR_NUMERICAL, "Prog_MLFalign2D-ERROR%% Avoid different Q0 values in the CTF parameter files!");
            }
            Maux.resize(dim, dim);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Maux)
            {
                if (phase_flipped)
                    dAij(Maux, i, j) = fabs(dAij(ctfmask, i, j).real());
                else
                    dAij(Maux, i, j) = dAij(ctfmask, i, j).real();
            }
            CenterFFT(Maux, true);
            center.initZeros();
            rmean_ctf.initZeros();
            Maux.setXmippOrigin();
            radialAverage(Maux, center, rmean_ctf, radial_count, true);
            Vctf.push_back(dum);
            Vdec.push_back(dum);
            Vsig.push_back(dum);
            Vctf[ifocus].resize(hdim);
            FOR_ALL_DIGITAL_FREQS()
            {
                VCTF_ITEM = dAi(rmean_ctf, irr);
            }
            Vdec[ifocus].resize(Vctf[ifocus]);
            Vsig[ifocus].resize(Vctf[ifocus]);
            ++ifocus;
        }

        MetaData md(MDimg);
        MDimg.joinNatural(md, mdCTF);
        if (IS_MASTER)
        {
        mdCTF.write("tmp_ctf.xmd");
        md.write("tmp_md.xmd");
        MDimg.write("tmp_images.xmd");
        exit(1);
        }

    }

    // Get a resolution pointer in Fourier-space
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(Maux)
    {
        A2D_ELEM(Maux, i, j) = XMIPP_MIN((double) (hdim - 1), sqrt((double)(i * i + j * j)));
    }
    CenterFFT(Maux, false);
    Mresol_int.resize(dim, dim);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mresol_int)
    {
        dAij(Mresol_int, i, j) = ROUND(dAij(Maux, i, j));
    }

    // Check whether to generate new references
    do_generate_refs = false;
    if (fn_ref.empty())
    {
        if (model.n_ref > 0)
        {
            fn_ref = FN_CLASSES_MD(getIterExtraPath(fn_root, 0));
            do_generate_refs = true;

            if (rank == 0)
                generateInitialReferences();
        }
        else
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "Please provide -ref or -nref larger than zero");
    }

    show();
    estimateInitialNoiseSpectra();
}

// Read reference images to memory and initialize offset vectors
// This side info is NOT general, i.e. in parallel mode it is NOT the
// same for all processors! (in contrast to produce_Side_info)
void ProgMLF2D::produceSideInfo2()
{
    LOG_LEVEL(produceSideInfo2);
    int                       c, idum, refno = 0;
    double                    aux = 0.;
    FileName                  fn_tmp;
    Image<double>                img;
    std::vector<double>       Vdum, sumw_defocus;

    // Read in all reference images in memory
    MDref.read(fn_ref);

    model.n_ref = MDref.size();
    refno = 0;
    double ref_fraction = (double)1. / model.n_ref;
    double ref_weight = (double)nr_images_global / model.n_ref;

    FOR_ALL_OBJECTS_IN_METADATA(MDref)
    {
        MDref.getValue(MDL_IMAGE, fn_tmp, __iter.objId);
        img.read(fn_tmp);
        img().setXmippOrigin();
        model.Iref.push_back(img);
        Iold.push_back(img);
        Ictf.push_back(img);
        // Default start is all equal model fractions
        alpha_k.push_back(ref_fraction);
        model.Iref[refno].setWeight(ref_weight);
        // Default start is half-half mirrored images
        mirror_fraction.push_back((do_mirror ? 0.5 : 0.));
        refno++;
    }

    //This will differ from nr_images_global if MPI
    nr_images_local = divide_equally(nr_images_global, size, rank, myFirstImg, myLastImg);
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "nr_images_local= "<< nr_images_local<<std::endl;
    std::cerr << "myFirstImg= "<< myFirstImg<<std::endl;
    std::cerr << "myLastImg= "<< myLastImg<<std::endl;
    std::cerr <<"size="<<size<<"rank="<<rank<<std::endl;

#endif
    //--------Setup for Docfile -----------
    docfiledata.resize(nr_images_local, DATALINELENGTH);

    if (do_norm)
    {
        average_scale = 1.;
        refs_avgscale.assign(model.n_ref, 1.);
        imgs_scale.assign(nr_images_local, 1.);
    }

    // Fill imgs_offsets vectors with zeros
    idum = (do_mirror ? 4 : 2) * model.n_ref;

    FOR_ALL_LOCAL_IMAGES()
    {
        imgs_offsets.push_back(Vdum);
        for (int refno = 0; refno < idum; refno++)
        {
            imgs_offsets[IMG_LOCAL_INDEX].push_back(0.);
        }
    }

    // For limited orientational search: fill imgs_oldphi & imgs_oldtheta
    // (either read from fn_doc or initialize to -999.)
    if (limit_rot)
    {
        imgs_oldphi.assign(nr_images_local, -999.);
        imgs_oldtheta.assign(nr_images_local, -999.);
    }

    if (do_restart)
    {
        // Read optimal image-parameters
        //        FOR_ALL_LOCAL_IMAGES()
        //        {
        //            if (limit_rot || do_norm)
        //            {
        //                MDimg.getValue(MDL_ANGLE_ROT, imgs_oldphi[IMG_LOCAL_INDEX]);
        //                MDimg.getValue(MDL_ANGLE_TILT, imgs_oldtheta[IMG_LOCAL_INDEX]);
        //            }
        //
        //            idum = (do_mirror ? 2 : 1) * model.n_ref;
        //            double xoff, yoff;
        //            MDimg.getValue(MDL_SHIFT_X, xoff);
        //            MDimg.getValue(MDL_SHIFT_Y, yoff);
        //            for (int refno = 0; refno < idum; refno++)
        //            {
        //                imgs_offsets[IMG_LOCAL_INDEX][2 * refno] = xoff;
        //                imgs_offsets[IMG_LOCAL_INDEX][2 * refno + 1] = yoff;
        //            }
        //
        //            if (do_norm)
        //            {
        //                MDimg.getValue(MDL_INTSCALE, imgs_scale[IMG_LOCAL_INDEX]);
        //            }
        //        }

        // read Model parameters
        refno = 0;
        double sumw = 0.;
        FOR_ALL_OBJECTS_IN_METADATA(MDref)
        {
            MDref.getValue(MDL_WEIGHT, alpha_k[refno], __iter.objId);
            sumw += alpha_k[refno];
            if (do_mirror)
                MDref.getValue(MDL_MIRRORFRAC,
                               mirror_fraction[refno], __iter.objId);
            if (do_norm)
                MDref.getValue(MDL_INTSCALE, refs_avgscale[refno], __iter.objId);
            refno++;
        }
        FOR_ALL_MODELS()
        {
            alpha_k[refno] /= sumw;
        }
    }

    MultidimArray<std::complex<double> >   Faux;
    MultidimArray<double>       Maux;
    MultidimArray<double>       dum, rmean_signal2, spectral_signal;
    Matrix1D<int>                center(2);
    MultidimArray<int>           radial_count;
    std::ifstream                fh;



    center.initZeros();
    // Calculate average spectral signal
    c = 0;
    FOR_ALL_MODELS()
    {
        if (alpha_k[refno] > 0.)
        {
            FourierTransform(model.Iref[refno](), Faux);
            FFT_magnitude(Faux, Maux);
            CenterFFT(Maux, true);
            Maux *= Maux;
            Maux *= alpha_k[refno] * (double)nr_images_global;
            Maux.setXmippOrigin();
            rmean_signal2.initZeros();
            radialAverage(Maux, center, rmean_signal2, radial_count, true);
            if (c == 0)
                spectral_signal = rmean_signal2;
            else
                spectral_signal += rmean_signal2;
            c++;
        }
    }
    if (do_ML3D)
    {
        // Divide by the number of reference volumes
        // But I don't know alpha (from 3DSSNR) yet:
        // Introduce a fudge-factor of 2 to prevent over-estimation ...
        // I think this is irrelevant, as the spectral_signal is
        // re-calculated in ml_refine3d.cpp
        double inv_2_nr_vol = 1. / (double)(nr_vols * 2);
        spectral_signal *= inv_2_nr_vol;
    }
    else
    {
        // Divide by the number of reference images
        double inv_n_ref = 1. / (double) model.n_ref;
        spectral_signal *= inv_n_ref;
    }

    // Read in Vsig-vectors with fixed file names
    FileName fn_base = FN_ITER_BASE(istart - 1);
    MetaData md;

    FOR_ALL_DEFOCUS_GROUPS()
    {
        fn_tmp = FN_VSIG(fn_base, ifocus, "noise.xmd");
        md.read(fn_tmp);

        FOR_ALL_OBJECTS_IN_METADATA(md)
        {
            md.getValue(MDL_MLF_NOISE, aux, __iter.objId);
            dAi(Vsig[ifocus], __iter.objIndex) = aux;
        }

#ifdef OLD_FILE
        fh.open((fn_tmp).c_str(), std::ios::in);
        if (!fh)
            REPORT_ERROR(ERR_IO_NOTOPEN, (String)"Prog_MLFalign2D_prm: Cannot read file: " + fn_tmp);
        else
        {
            FOR_ALL_DIGITAL_FREQS()
            {
                fh >> aux;
                if (ABS(aux - ((double)irr/(sampling*dim)) ) > 0.01 )
                {
                    std::cerr<<"aux= "<<aux<<" resol= "<<(double)irr/(sampling*dim)<<std::endl;
                    REPORT_ERROR(ERR_NUMERICAL, (String)"Prog_MLFalign2D_prm: Wrong format: " + fn_tmp);
                }
                fh >> aux;
                VSIG_ITEM = aux;
            }
        }
        fh.close();
#endif
        // Initially set sumw_defocus equal to count_defocus
        sumw_defocus.push_back((double)count_defocus[ifocus]);
    }
    // Calculate all Wiener filters
    updateWienerFilters(spectral_signal, sumw_defocus, istart - 1);

}

// For initial noise variances
void ProgMLF2D::estimateInitialNoiseSpectra()
{
  LOG_FUNCTION();
    // For first iteration only: calculate sigma2 (& spectral noise) from power
    // spectra of all images and subtract power spectrum of average image
    // (to take away low-res frequencies where the signal dominates!)
    if (istart == 1)
    {

        MultidimArray<double>            Maux, Mallave;
        MultidimArray<std::complex<double> >  Fimg, Faux, Fave;
        MultidimArray<double>            rmean_noise, rmean_signal, rmean_avesignal;
        std::vector<MultidimArray<double> >   Msigma2, Mave;
        Matrix1D<int>               center(2);
        MultidimArray<int>           radial_count;
        Image<double>                       img;
        FileName                    fn_tmp;
        std::ofstream                    fh;

        center.initZeros();

        int focus = 0;

        if (verbose > 0)
            std::cout << "--> Estimating initial noise models from average power spectra ..." << std::endl;

        initProgress(MDimg.size());

        Msigma2.clear();
        Msigma2.resize(nr_focus);
        Mave.clear();
        Mave.resize(nr_focus);

        FOR_ALL_DEFOCUS_GROUPS()
        {
            Msigma2[ifocus].initZeros(dim, dim);
            Msigma2[ifocus].setXmippOrigin();
            Mave[ifocus].initZeros(dim, dim);
            Mave[ifocus].setXmippOrigin();
        }

        Maux.initZeros(dim, dim);
        int imgno = 0;

        FOR_ALL_OBJECTS_IN_METADATA(MDimg)
        {
            focus = 0;
            if (do_ctf_correction)
                MDimg.getValue(MDL_DEFGROUP, focus, __iter.objId);
            MDimg.getValue(MDL_IMAGE, fn_tmp, __iter.objId);
            //img.read(fn_tmp, false, false, false, false);
            //TODO: Check this????
            img.read(fn_tmp);
            img().setXmippOrigin();
            FourierTransform(img(), Fimg);
            FFT_magnitude(Fimg, Maux);
            Maux.setXmippOrigin();
            Maux *= Maux;
            Msigma2[focus] += Maux;
            Mave[focus] += img();
            setProgress(++imgno);
        }

        endProgress();

        FileName fn_base = FN_ITER_BASE(istart - 1);

        // Calculate Vsig vectors and write them to disc

        FOR_ALL_DEFOCUS_GROUPS()
        {
            Mave[ifocus] /= (double)count_defocus[ifocus];
            FourierTransform(Mave[ifocus], Fave);
            FFT_magnitude(Fave, Maux);
            Maux *= Maux;
            CenterFFT(Maux, true);
            Maux.setXmippOrigin();
            rmean_signal.initZeros();
            radialAverage(Maux, center, rmean_signal, radial_count, true);
            CenterFFT(Msigma2[ifocus], true);
            Msigma2[ifocus].setXmippOrigin();
            rmean_noise.initZeros();
            radialAverage(Msigma2[ifocus], center, rmean_noise, radial_count, true);
            rmean_noise /= (double)count_defocus[ifocus];
            // Subtract signal terms
            // Divide by factor 2 because of the 2D-Gaussian distribution!
            FOR_ALL_DIGITAL_FREQS()
            VSIG_ITEM = (dAi(rmean_noise, irr) - dAi(rmean_signal, irr)) / 2.;


            // write Vsig vector to disc (only master)
            if (IS_MASTER)
                writeNoiseFile(fn_base, ifocus);

            fn_tmp = FN_VSIG(fn_base, ifocus, "noise.xmd");

#ifdef OLD_FILE

            fh.open((fn_tmp).c_str(), std::ios::out);
            if (!fh)
                REPORT_ERROR(ERR_IO_NOTOPEN, (String)"Prog_MLFalign2D_prm: Cannot write file: " + fn_tmp);
            FOR_ALL_DIGITAL_FREQS()
            {
                fh << (double)irr/(sampling*dim) << " " << VSIG_ITEM << "\n";
            }
            fh.close();
#endif

        }
        Msigma2.clear();
        Mave.clear();
    }

}

void ProgMLF2D::updateWienerFilters(const MultidimArray<double> &spectral_signal,
                                    std::vector<double> &sumw_defocus,
                                    int iter)
{
    LOG_FUNCTION();
    // Use formula 2.32b on p60 from Frank's book 2nd ed.,
    // Assume that Vctf, Vdec and Vsig exist (with their right sizes)
    // and that Vctf, Vsig are already filled with the right values

    std::vector<MultidimArray<double> >  Vsnr;
    MultidimArray<double>           Vzero, Vdenom, Vavgctf2;
    MultidimArray<double>           Maux;
    MultidimArray<std::complex<double> > Faux;
    std::ofstream                   fh;
    size_t                     maxres = 0;
    double                     noise, sum_sumw_defocus = 0.;
    size_t                     int_lowres_limit, int_highres_limit, int_ini_highres_limit;
    size_t                     current_probres_limit;
    FileName                   fn_base, fn_tmp;

    // integer resolution limits (in shells)
    int_lowres_limit  = (size_t)(sampling * dim / lowres_limit);
    int_highres_limit = (highres_limit > 0.) ? ROUND(sampling * dim / highres_limit) : hdim;
    int_ini_highres_limit =  (ini_highres_limit > 0.) ? ROUND(sampling * dim / ini_highres_limit) : hdim;

    // Pre-calculate average CTF^2 and initialize Vsnr
    Vavgctf2.initZeros(hdim);
    Vzero.initZeros(hdim);
    Vsnr.resize(nr_focus, Vzero);

    FOR_ALL_DEFOCUS_GROUPS()
    {
        double & defocus_weight = sumw_defocus[ifocus];
        sum_sumw_defocus += defocus_weight;
        FOR_ALL_DIGITAL_FREQS()
        dAi(Vavgctf2, irr) += VCTF_ITEM * VCTF_ITEM * defocus_weight;
    }
    Vavgctf2 /= sum_sumw_defocus;

    // std::cerr << "DEBUG_JM, updateWienerFilters: spectral_signal: " << spectral_signal << std::endl;

    // Calculate SSNR for all CTF groups
    // For each group the spectral noise is estimated via (2*Vsig)/(sumw_defocus-1)
    // The spectral signal is not split into CTF groups
    // Therefore, affect the average spectral_signal for each defocu
    // group with its CTF^2 and divide by the average CTF^2
    FOR_ALL_DEFOCUS_GROUPS()
    {
        FOR_ALL_DIGITAL_FREQS()
        {
            //if (verbose) std::cerr << "DEBUG_JM: ifocus: " << ifocus << std::endl;
            //if (verbose) std::cerr << "DEBUG_JM: irr: " << irr << std::endl;

            noise = 2. * VSIG_ITEM * sumw_defocus[ifocus];
            noise /= sumw_defocus[ifocus] - 1;

            //if (verbose) std::cerr << "DEBUG_JM: noise: " << noise << std::endl;
            //if (verbose) std::cerr << "DEBUG_JM: VSNR_ITEM: " << VSNR_ITEM << std::endl;
            VSNR_ITEM = 0.;

            if (noise > 1e-20 && dAi(Vavgctf2, irr) > 0.)
            {
                VSNR_ITEM = reduce_snr * VCTF_ITEM * VCTF_ITEM *
                            dAi(spectral_signal, irr) / (dAi(Vavgctf2, irr) * noise);
                //if (verbose) std::cerr << "DEBUG_JM: inside cond: (noise > 1e-20 && dAi(Vavgctf2, irr)" <<std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: VCTF_ITEM: " << VCTF_ITEM << std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: dAi(spectral_signal, irr): " << dAi(spectral_signal, irr) << std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: dAi(Vavgctf2, irr): " << dAi(Vavgctf2, irr) << std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: noise: " << noise << std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: VSNR_ITEM: " << VSNR_ITEM << std::endl;
            }
            // For start from already CTF-deconvoluted references:
            if ((iter == istart - 1) && !first_iter_noctf)
            {
                VSNR_ITEM *= dAi(Vavgctf2, irr);
                //if (verbose) std::cerr << "DEBUG_JM: inside cond: (noise > 1e-20 && dAi(Vavgctf2, irr)" <<std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: VSNR_ITEM: " << VSNR_ITEM << std::endl;
            }
            // Take ini_highres_limit into account (only for first iteration)
            if ( iter == 0 && ini_highres_limit > 0. && irr > int_ini_highres_limit )
            {
                VSNR_ITEM = 0.;
                //if (verbose) std::cerr << "DEBUG_JM: inside cond: iter == 0 && ini_highres_limit > 0. && irr > int_ini_highres_limit" <<std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: VSNR_ITEM: " << VSNR_ITEM << std::endl;
            }
            // Subtract 1 according Unser et al.
            VSNR_ITEM = XMIPP_MAX(0., VSNR_ITEM - 1.);
            // Prevent spurious high-frequency significant SNRs from random averages
            if (iter == 0 && do_generate_refs)
            {
                VSNR_ITEM = XMIPP_MAX(0., VSNR_ITEM - 2.);
                //if (verbose) std::cerr << "DEBUG_JM: inside cond: iter == 0 && do_generate_refs" <<std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: VSNR_ITEM: " << VSNR_ITEM << std::endl;
            }
            if (VSNR_ITEM > 0. && irr > maxres)
            {
                maxres = irr;
                //if (verbose) std::cerr << "DEBUG_JM: inside cond: VSNR_ITEM > 0. && irr > maxres" <<std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: VSNR_ITEM: " << VSNR_ITEM << std::endl;
                //if (verbose) std::cerr << "DEBUG_JM: maxres: " << maxres << std::endl;
            }
        }
    }

    // Check that at least some frequencies have non-zero SSNR...
    if (maxres == 0)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "ProgMLF2D: All frequencies have zero spectral SNRs... (increase -reduce_snr) ");

    if (do_ctf_correction)
    {
        // Pre-calculate denominator of eq 2.32b of Frank's book (2nd ed.)
        Vdenom.initZeros(hdim);
        FOR_ALL_DIGITAL_FREQS()
        {
            FOR_ALL_DEFOCUS_GROUPS()
            {
                dAi(Vdenom, irr) += sumw_defocus[ifocus] * VSNR_ITEM * VCTF_ITEM * VCTF_ITEM;
            }
            dAi(Vdenom, irr) += 1.;
            dAi(Vdenom, irr) /= sum_sumw_defocus;
        }

        // Calculate Wiener filters
        FOR_ALL_DEFOCUS_GROUPS()
        {
            FOR_ALL_DIGITAL_FREQS()
            {
                if (VSNR_ITEM > 0.)
                {
                    VDEC_ITEM = VSNR_ITEM * VCTF_ITEM / dAi(Vdenom, irr);
                    // Prevent too strong Wiener filter artefacts
                    VDEC_ITEM = XMIPP_MIN(10., VDEC_ITEM);
                }
                else
                {
                    VDEC_ITEM = 0.;
                }
            }
        }
    }

    // Write Wiener filters and spectral SNR to text files
    double resol_step = 1. / (sampling * dim);
    double resol_freq = 0;
    double resol_real = 999.;
    MDRow row;

    if (verbose > 0)
    {
        fn_base = FN_ITER_BASE(iter);
        // CTF group-specific Wiener filter files
        FOR_ALL_DEFOCUS_GROUPS()
        {
            resol_freq = 0;
            MetaData md;
            FOR_ALL_DIGITAL_FREQS()
            {
                noise = 2. * VSIG_ITEM * sumw_defocus[ifocus];
                noise /= (sumw_defocus[ifocus] - 1);
                if (irr > 0)
                    resol_real = (sampling*dim)/(double)irr;
                row.setValue(MDL_RESOLUTION_FREQ, resol_freq);
                row.setValue(MDL_RESOLUTION_SSNR, XMIPP_MIN(25., VSNR_ITEM));
                row.setValue(MDL_MLF_CTF, VCTF_ITEM);
                row.setValue(MDL_MLF_WIENER, VDEC_ITEM);
                row.setValue(MDL_MLF_SIGNAL, dAi(spectral_signal, irr));
                row.setValue(MDL_MLF_NOISE, noise);
                row.setValue(MDL_RESOLUTION_FREQREAL, resol_real);
                md.addRow(row);
                resol_freq += resol_step;
            }
            fn_tmp = FN_VSIG(fn_base, ifocus, "ssnr.xmd");
            md.write(fn_tmp, (ifocus > 0 ? MD_APPEND : MD_OVERWRITE));

        }
    }

    // Set the current resolution limits
    if (do_include_allfreqs)
    {
        current_probres_limit = hdim-1;
        current_highres_limit = hdim-1;

    }
    else if (fix_high > 0.)
    {
        current_probres_limit = ROUND((sampling*dim)/fix_high);
        current_highres_limit = ROUND((sampling*dim)/fix_high);

    }
    else
    {
        current_probres_limit = maxres;
        current_highres_limit = maxres + 5; // hard-code increase_highres_limit to 5

    }

    // Set overall high resolution limit
    current_probres_limit = XMIPP_MIN(current_probres_limit, int_highres_limit);
    current_highres_limit = XMIPP_MIN(current_highres_limit, int_highres_limit);

    current_probres_limit = XMIPP_MIN(current_probres_limit, hdim);
    current_highres_limit = XMIPP_MIN(current_highres_limit, hdim);
    //std::cerr << "DEBUG_JM(6): current_probres_limit: " << current_probres_limit << std::endl;

    if (debug>0)
    {
        std::cerr
        << "current_probres_limit dependencies: " << std::endl
        << " hdim: " << hdim << " dim: " << dim << std::endl
        << " sampling: " << sampling << " fix_high: " << fix_high << std::endl
        << " maxres: " << maxres << " int_highres_limit: " << int_highres_limit << std::endl;

        std::cerr<<" Current resolution limits: "<<std::endl;
        std::cerr<<" + low   res= "<<lowres_limit<<" Ang ("<<int_lowres_limit<<" shell)"<<std::endl;
        std::cerr<<" + prob. res= "<<sampling*dim/current_probres_limit<<" Ang ("<<current_probres_limit<<" shell)"<<std::endl;
        std::cerr<<" + extra res= "<<sampling*dim/current_highres_limit<<" Ang ("<<current_highres_limit<<" shell)"<<std::endl;
        std::cerr<<" + high  res= "<<highres_limit<<" Ang ("<<int_highres_limit<<" shell)"<<std::endl;
    }

    // Get the new pointers to all pixels in FourierTransformHalf
    Maux.initZeros(dim, dim);
    Maux.setXmippOrigin();
    FourierTransformHalf(Maux, Faux);
    pointer_2d.clear();
    pointer_i.clear();
    pointer_j.clear();

    // First, get the pixels to use in the probability calculations:
    // These are within [lowres_limit,current_probres_limit]
    nr_points_prob = 0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Faux)
    {
    	size_t ires = dAij(Mresol_int, i, j);
        if (ires > int_lowres_limit &&
            ires <= current_probres_limit &&
            !(i == 0 && j > hdim) ) // exclude first half row in FourierTransformHalf
        {
            pointer_2d.push_back(i*XSIZE(Maux) + j);
            pointer_i.push_back(i);
            pointer_j.push_back(j);
            ++nr_points_prob;
        }
    }
    // Second, get the rest of the currently relevant pixels:
    // These are within [0,lowres_limit> and between <current_probres_limit,current_highres_limit]
    nr_points_2d = nr_points_prob;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Faux)
    {
        size_t ires = dAij(Mresol_int, i, j);
        if ( (ires <= int_lowres_limit || ires > current_probres_limit) &&
             ires <= current_highres_limit &&
             !(i == 0 && j > hdim) ) // exclude first half row in FourierTransformHalf
        {
            pointer_2d.push_back(i*XSIZE(Maux) + j);
            pointer_i.push_back(i);
            pointer_j.push_back(j);
            nr_points_2d++;
        }
    }
    dnr_points_2d = 2 * nr_points_2d;

    if (debug>0)
    {
        std::cerr<<"nr_points_2d= "<<nr_points_2d<<" nr_points_prob= "<<nr_points_prob<<std::endl;
    }

    // For variable in-plane sampling rates
    setCurrentSamplingRates(sampling*dim/current_probres_limit);

}


// Vary psi_step and trans_step with resolution =============================================
void ProgMLF2D::setCurrentSamplingRates(double current_probres_limit)
{

    int trans_step = 1;

    if (do_variable_psi)
    {
        // Sample the in-plane rotations 3x the current resolution
        // Note that SIND(0.5*psi_step) = Delta / dim;
        psi_step = (2.* ASIND(current_probres_limit/(dim*sampling))) / 3.;
        nr_psi = CEIL(psi_max / psi_step);
        // Use user-provided psi_step as minimum sampling
        nr_psi = XMIPP_MIN(nr_psi, max_nr_psi);
        psi_step = psi_max / nr_psi;
    }

    if (do_variable_trans)
    {
        Matrix1D<double>  offsets(2);
        nr_trans = 0;
        Vtrans.clear();
        // Sample the in-plane translations 3x the current resolution
        trans_step = ROUND(current_probres_limit/(3.*sampling));
        // Use trans_step=1 as minimum sampling
        trans_step = XMIPP_MAX(1,trans_step);
        for (int ix = -search_shift*trans_step; ix <= search_shift*trans_step; ix+=trans_step)
        {
            for (int iy = -search_shift*trans_step; iy <= search_shift*trans_step; iy+=trans_step)
            {
                double rr = sqrt((double)(ix * ix + iy * iy));
                if (rr <= (double)trans_step*search_shift)
                {
                    XX(offsets) = ix;
                    YY(offsets) = iy;
                    Vtrans.push_back(offsets);
                    nr_trans++;
                    if (ix == 0 && iy == 0)
                    {
                        zero_trans = nr_trans;
                        // For coarser samplings, always add (-1,0) (1,0) (0,-1) and (0,1)
                        if (trans_step > 1)
                        {
                            XX(offsets) = 1;
                            YY(offsets) = 0;
                            Vtrans.push_back(offsets);
                            nr_trans++;
                            XX(offsets) = -1;
                            YY(offsets) = 0;
                            Vtrans.push_back(offsets);
                            nr_trans++;
                            XX(offsets) = 0;
                            YY(offsets) = 1;
                            Vtrans.push_back(offsets);
                            nr_trans++;
                            XX(offsets) = 0;
                            YY(offsets) = -1;
                            Vtrans.push_back(offsets);
                            nr_trans++;
                        }
                    }
                }
            }
        }
    }
    if (verbose > 0 && (do_variable_psi || do_variable_trans))
    {
        std::cout<<" Current resolution= "<<current_probres_limit<<" Ang; current psi_step = "<<psi_step<<" current trans_step = "<<trans_step<<std::endl;

    }
}

// Generate initial references =============================================
void ProgMLF2D::generateInitialReferences()
{

    if (verbose > 0)
    {
        std::cout << "  Generating initial references by averaging over random subsets" << std::endl;
        init_progress_bar(model.n_ref);
    }

    Image<double> IRef, ITemp;
    FileName fn_tmp;
    FileName fn_base = getIterExtraPath(fn_root, 0);

    //randomizeImagesOrder();

    MDref.clear();
    size_t nsub, first, last;
    size_t id;

    for (int refno = 0; refno < model.n_ref; refno++)
    {
        nsub = divide_equally(nr_images_global, model.n_ref, refno, first, last);
        //Clear images
        IRef().initZeros(dim, dim);
        IRef().setXmippOrigin();

        for (size_t imgno = first; imgno <= last; imgno++)
        {
            MDimg.getValue(MDL_IMAGE, fn_tmp, img_id[imgno]);
            std::cerr << "DEBUG_JM: fn_tmp: " << fn_tmp << std::endl;
            ITemp.read(fn_tmp);
            ITemp().setXmippOrigin();
            IRef() += ITemp();
        }

        fn_tmp = FN_REF(fn_base, refno + 1);

        IRef() /= nsub;
        IRef.write(fn_tmp);
        id = MDref.addObject();
        MDref.setValue(MDL_IMAGE, fn_tmp, id);
        MDref.setValue(MDL_ENABLED, 1, id);

        if (verbose > 0)
            progress_bar(refno);
    }//close for refno

    if (verbose > 0)
        progress_bar(model.n_ref);

    fn_ref = FN_CLASSES_MD(fn_base);
    MDref.write(fn_ref);
}


// Calculate probability density function of all in-plane transformations phi
void ProgMLF2D::calculatePdfInplane()
{

    double r2, pdfpix, sum;
    P_phi.resize(dim, dim);
    P_phi.setXmippOrigin();
    Mr2.resize(dim, dim);
    Mr2.setXmippOrigin();

    sum=0.;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(P_phi)
    {
        r2 = (double)(j * j + i * i);
        if (sigma_offset > 0.)
        {
            pdfpix = exp(-r2 / (2 * sigma_offset * sigma_offset));
            pdfpix /= 2 * PI * sigma_offset * sigma_offset * nr_psi * nr_nomirror_flips;
        }
        else
        {
            if (j == 0 && i == 0)
                pdfpix = 1.;
            else
                pdfpix = 0.;
        }
        A2D_ELEM(P_phi, i, j) = pdfpix;
        A2D_ELEM(Mr2, i, j) = (float)r2;
        sum+=pdfpix;
    }
    // Normalization
    P_phi/=sum;

}

void ProgMLF2D::appendFTtoVector(const MultidimArray<std::complex<double> > &in,
                                 std::vector<double> &out)
{

    // First, store the points used in the probability calculations
    std::complex<double> * tmp_in = MULTIDIM_ARRAY(in);
    for (size_t ipoint = 0; ipoint < nr_points_2d; ipoint++)
    {
        int ii = pointer_2d[ipoint];
        out.push_back(tmp_in[ii].real());
        out.push_back(tmp_in[ii].imag());
    }
}

void ProgMLF2D::getFTfromVector(const std::vector<double> &in,
                                const int start_point,
                                MultidimArray<std::complex<double> > &out,
                                bool only_real)
{

    out.resize(hdim + 1, dim);
    out.initZeros();
    std::complex<double> * tmp_out = MULTIDIM_ARRAY(out);
    if (only_real)
    {
        for (size_t ipoint = 0; ipoint < nr_points_2d; ipoint++)
        {
            int ii = pointer_2d[ipoint];
            std::complex<double> aux(in[start_point+ipoint], 0.);
            tmp_out[ii] = aux;
        }
    }
    else
    {
        for (size_t ipoint = 0; ipoint < nr_points_2d; ipoint++)
        {
            int ii = pointer_2d[ipoint];
            std::complex<double> aux(in[start_point+2*ipoint],in[start_point+2*ipoint+1]);
            tmp_out[ii] = aux;
        }
    }


}

// Rotate reference for all models and rotations and fill Fref vectors =============
void ProgMLF2D::rotateReference(std::vector<double> &out)
{
    out.clear();

    double AA, stdAA, psi;
    MultidimArray<double> Maux;
    MultidimArray<std::complex<double> > Faux;

    Maux.initZeros(dim, dim);
    Maux.setXmippOrigin();

    FOR_ALL_MODELS()
    {
        FOR_ALL_ROTATIONS()
        {
            // Add arbitrary number (small_angle) to avoid 0-degree rotation (lacking interpolation)
            psi = (double)(ipsi * psi_max / nr_psi) + SMALLANGLE;
            //model.Iref[refno]().rotateBSpline(3, psi, Maux, WRAP);
            rotate(BSPLINE3, Maux, model.Iref[refno](), -psi, 'Z', WRAP);
            FourierTransformHalf(Maux, Faux);

            // Normalize the magnitude of the rotated references to 1st rot of that ref
            // This is necessary because interpolation due to rotation can lead to lower overall Fref
            // This would result in lower probabilities for those rotations
            AA = Maux.sum2();
            if (ipsi == 0)
            {
                stdAA = AA;
            }
            if (AA > 0)
            {
                double sqrtVal = sqrt(stdAA/AA);
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                dAi(Faux, n) *= sqrtVal;
            }
            // Add all points as doubles to the vector
            appendFTtoVector(Faux, out);
        }
        // Free memory
        model.Iref[refno]().resize(0,0);
    }
}


// Collect all rotations and sum to update Iref() for all models ==========
void ProgMLF2D::reverseRotateReference(const std::vector<double> &in,
                                       std::vector<MultidimArray<double > > &out)
{

    double psi;
    MultidimArray<double> Maux, Maux2;
    MultidimArray<std::complex<double> > Faux;
    Maux.resize(dim, dim);
    Maux2.resize(dim, dim);
    Maux.setXmippOrigin();
    Maux2.setXmippOrigin();

    out.clear();
    FOR_ALL_MODELS()
    {
        Maux.initZeros();
        out.push_back(Maux);
        FOR_ALL_ROTATIONS()
        {
            // Add arbitrary number to avoid 0-degree rotation without interpolation effects
            psi = (double)(ipsi * psi_max / nr_psi) + SMALLANGLE;
            getFTfromVector(in, refno*nr_psi*dnr_points_2d + ipsi*dnr_points_2d, Faux);
            InverseFourierTransformHalf(Faux, Maux, dim);
            //Maux.rotateBSpline(3, -psi, Maux2, WRAP);
            rotate(BSPLINE3, Maux2, Maux, psi, 'Z', WRAP);
            out[refno] += Maux2;
        }
    }

}

void ProgMLF2D::preselectDirections(float &phi, float &theta,
                                    std::vector<double> &pdf_directions)
{

    float phi_ref, theta_ref, angle, angle2;
    Matrix1D<double> u, v;

    pdf_directions.clear();
    pdf_directions.resize(model.n_ref);
    FOR_ALL_MODELS()
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

void ProgMLF2D::fourierTranslate2D(const std::vector<double> &in,
                                   MultidimArray<double> &trans,
                                   std::vector<double> &out,
                                   int point_start)
{
    double xx, yy, xxshift, yyshift, dotp;
    double a, b, c, d, ac, bd, ab_cd;
    xxshift = -trans(0) / (double)dim;
    yyshift = -trans(1) / (double)dim;
    //Not very clean, but very fast
    const double * ptrIn = &(in[point_start]);
    int * ptrJ = &(pointer_j[0]);
    int * ptrI = &(pointer_i[0]);

    for (size_t i = 0; i < nr_points_2d; i++)
    {
        xx = (double)ptrJ[i];
        yy = (double)ptrI[i];
        dotp = 2 * PI * (xx * xxshift + yyshift * yy);
        a = cos(dotp);
        b = sin(dotp);
        c = *ptrIn++;//in[point_start + 2*i];
        d = *ptrIn++;//in[point_start + 2*i+1];
        ac = a * c;
        bd = b * d;
        ab_cd = (a + b) * (c + d);
        out.push_back(ac - bd); // real
        out.push_back(ab_cd - ac - bd); // imag
        //*ptrOut = ac - bd; ++ptrOut;// real
        //*ptrOut = ab_cd - ac - bd; ++ptrOut;// imag
    }

}

void ProgMLF2D::calculateFourierOffsets(const MultidimArray<double> &Mimg,
                                        const std::vector<double > &offsets,
                                        std::vector<double>  &out,
                                        MultidimArray<int> &Moffsets,
                                        MultidimArray<int> &Moffsets_mirror)
{

    int irefmir, ix, iy, iflip_start, iflip_stop, count;
    int nr_mir = (do_mirror) ? 2 : 1;

    double dxx, dyy;
    std::vector<double> Fimg_flip;
    MultidimArray<std::complex<double> > Fimg;
    MultidimArray<double> trans(2);
    MultidimArray<double> Maux2, Maux;

    Moffsets.resize(dim, dim);
    Moffsets.setXmippOrigin();
    Moffsets.initConstant(-1);
    Moffsets_mirror.resize(dim, dim);
    Moffsets_mirror.setXmippOrigin();
    Moffsets_mirror.initConstant(-1);
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Maux2.resize(dim, dim);
    Maux2.setXmippOrigin();

    // Flip images and store the precalculates Fourier Transforms in Fimg_flip
    out.clear();
    FOR_ALL_FLIPS()
    {
        Maux.setXmippOrigin();
        applyGeometry(LINEAR, Maux, Mimg, F[iflip], IS_INV, WRAP);

        FourierTransformHalf(Maux, Fimg);
        appendFTtoVector(Fimg,Fimg_flip);
    }

    // Now for all relevant offsets calculate the Fourier transforms
    // Matrices Moffsets & Moffsets_mirror contain pointers for the
    // out vector (access as count*4*dnr_points_2d + iflip*dnr_points_2d)
    count = 0;
    FOR_ALL_MODELS()
    {
        for (int imir = 0; imir < nr_mir; imir++)
        {
            irefmir = imir * model.n_ref + refno;
            iflip_start = imir * nr_nomirror_flips;
            iflip_stop = imir * nr_nomirror_flips + nr_nomirror_flips;
            FOR_ALL_LIMITED_TRANSLATIONS()
            {
                ix = ROUND(offsets[2*irefmir] + Vtrans[itrans](0));
                iy = ROUND(offsets[2*irefmir+1] + Vtrans[itrans](1));
                dxx = (double)intWRAP(ix, Moffsets.startingX(), Moffsets.finishingX());
                dyy = (double)intWRAP(iy, Moffsets.startingY(), Moffsets.finishingY());
                // For non-mirrors
                if (imir == 0 && A2D_ELEM(Moffsets, ROUND(dyy), ROUND(dxx)) < 0)
                {
                    for (int iflip = iflip_start; iflip < iflip_stop; iflip++)
                    {
                        Matrix2D<double> & refF_iflip = F[iflip];
                        dAi(trans, 0) = dxx * MAT_ELEM(refF_iflip, 0, 0) + dyy * MAT_ELEM(refF_iflip, 0, 1);
                        dAi(trans, 1) = dxx * MAT_ELEM(refF_iflip, 1, 0) + dyy * MAT_ELEM(refF_iflip, 1, 1);
                        fourierTranslate2D(Fimg_flip,trans,out,iflip*dnr_points_2d);
                    }
                    A2D_ELEM(Moffsets, ROUND(dyy), ROUND(dxx)) = count;
                    count++;
                }
                // For mirrors use a separate offset-matrix
                else if (imir == 1 && A2D_ELEM(Moffsets_mirror, ROUND(dyy), ROUND(dxx)) < 0)
                {
                    for (int iflip = iflip_start; iflip < iflip_stop; iflip++)
                    {
                        Matrix2D<double> & refF_iflip = F[iflip];
                        dAi(trans, 0) = dxx * MAT_ELEM(refF_iflip, 0, 0) + dyy * MAT_ELEM(refF_iflip, 0, 1);
                        dAi(trans, 1) = dxx * MAT_ELEM(refF_iflip, 1, 0) + dyy * MAT_ELEM(refF_iflip, 1, 1);
                        fourierTranslate2D(Fimg_flip,trans,out,iflip*dnr_points_2d);
                    }
                    A2D_ELEM(Moffsets_mirror, ROUND(dyy), ROUND(dxx)) = count;
                    count++;
                }
            }
        }
    }
}

double lcdf_tstudent_mlf1(double t)
{
    return cdf_tstudent(1,t);
}
double lcdf_tstudent_mlf3(double t)
{
    return cdf_tstudent(3,t);
}
double lcdf_tstudent_mlf6(double t)
{
    return cdf_tstudent(6,t);
}
double lcdf_tstudent_mlf9(double t)
{
    return cdf_tstudent(9,t);
}
double lcdf_tstudent_mlf30(double t)
{
    return cdf_tstudent(30,t);
}


// Exclude translations from the MLF_integration
// For significantly contributing refno+psi: re-calculate optimal shifts
void ProgMLF2D::processOneImage(const MultidimArray<double> &Mimg,
                                const size_t focus, bool apply_ctf,
                                double &fracweight, double &maxweight2,
                                double &sum_refw2, double &opt_scale,
                                size_t &opt_refno, double &opt_psi,
                                size_t &opt_ipsi, size_t &opt_iflip,
                                MultidimArray<double> &opt_offsets,
                                std::vector<double> &opt_offsets_ref,
                                std::vector<double > &pdf_directions,
                                bool do_kstest, bool write_histograms,
                                FileName fn_img, double &KSprob)
{
    std::vector<double> &Mwsum_sigma2_local = Mwsum_sigma2[focus];
    MultidimArray<double>                             Mweight;
    MultidimArray<int>                                Moffsets, Moffsets_mirror;
    std::vector<double>                               Fimg_trans;
    std::vector<MultidimArray<double> >               uniq_offsets;


    std::vector<double> refw(model.n_ref), refw2(model.n_ref), refw_mirror(model.n_ref), Pmax_refmir(2*model.n_ref);

    double aux, fracpdf, pdf, weight, weight2=0.;
    double tmpr, tmpi, sum_refw = 0.;
    double diff, maxweight = -99.e99, mindiff2 = 99.e99;
    double logsigma2, ldim, ref_scale = 1.;
    double scale_denom, scale_numer, wsum_sc = 0., wsum_sc2 = 0.;
    int    irot, irefmir, opt_irefmir=0, ix, iy;
    size_t point_trans=0;
    int    opt_itrans=0, iflip_start, iflip_stop, nr_mir;
    int    img_start, ref_start, wsum_start;

    if (!do_norm)
        opt_scale =1.;

    // Convert 1D Vsig vectors to the correct vector size for the 2D FTHalfs
    // Multiply by a factor of two because we consider both the real and the imaginary parts
    ldim = (double) 2 * nr_points_prob;
    // For t-student: df2 changes with effective resolution!
    if (do_student)
        df2 = - ( df +  ldim ) / 2. ;

    sum_refw2 = 0.;
    //change std::vectors by arrays
    double *sigma2=new double[nr_points_2d];
    double *inv_sigma2=new double [nr_points_2d];
    double *ctf=new double[nr_points_2d];
    double *decctf=new double[nr_points_2d];

    const int &ifocus = focus;
    FOR_ALL_POINTS()
    {
        int irr = DIRECT_MULTIDIM_ELEM(Mresol_int,pointer_2d[ipoint]);
        sigma2[ipoint] = 2. * VSIG_ITEM;
        inv_sigma2[ipoint] = 1./ sigma2[ipoint];
        decctf[ipoint] = VDEC_ITEM;
        ctf[ipoint] = VCTF_ITEM;
    }

    if (!apply_ctf)
    {
        FOR_ALL_POINTS()
        ctf[ipoint] = 1.;
    }

    // Precalculate normalization constant
    logsigma2 = 0.;
    double factor = (do_student) ? PI * df * 0.5 : PI;
    // Multiply by two because we treat real and imaginary parts!
    for (size_t ipoint = 0; ipoint < nr_points_prob; ipoint++)
        logsigma2 += 2 * log( sqrt(factor * sigma2[ipoint]));

    // Precalculate Fimg_trans, on pruned and expanded offset list
    calculateFourierOffsets(Mimg, opt_offsets_ref, Fimg_trans, Moffsets, Moffsets_mirror);

    Mweight.initZeros(nr_trans, model.n_ref, nr_flip*nr_psi);

    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            if (do_norm)
                ref_scale = opt_scale / refs_avgscale[refno];

            FOR_ALL_FLIPS()
            {
                irefmir = (int)floor(iflip / nr_nomirror_flips) * model.n_ref + refno;
                ix = (int)round(opt_offsets_ref[2*irefmir]);
                iy = (int)round(opt_offsets_ref[2*irefmir+1]);
                Pmax_refmir[irefmir] = 0.;
                point_trans = (iflip < nr_nomirror_flips) ? A2D_ELEM(Moffsets, iy, ix) : A2D_ELEM(Moffsets_mirror, iy, ix);
                if (point_trans < 0 || point_trans > dim2)
                {
                    std::cerr<<"point_trans = "<<point_trans<<" ix= "<<ix<<" iy= "<<iy<<std::endl;
                    REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,"mlf_align2d BUG: point_trans < 0");
                }

                fracpdf = (iflip < nr_nomirror_flips) ? 1. - mirror_fraction[refno] : mirror_fraction[refno];
                pdf =  alpha_k[refno] * fracpdf * A2D_ELEM(P_phi, iy, ix);
                // get the starting point in the Fimg_trans vector
                img_start = point_trans*4*dnr_points_2d + (iflip%nr_nomirror_flips)*dnr_points_2d;

                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    diff = 0.;
                    // get the starting point in the Fref vector
                    ref_start = refno * nr_psi * dnr_points_2d + ipsi * dnr_points_2d;
                    if (do_ctf_correction)
                    {
                        for (size_t ii = 0; ii < nr_points_prob; ii++)
                        {
                            tmpr = Fimg_trans[img_start + 2*ii] - ctf[ii] * ref_scale * Fref[ref_start + 2*ii];
                            tmpi = Fimg_trans[img_start + 2*ii+1] - ctf[ii] * ref_scale * Fref[ref_start + 2*ii+1];
                            tmpr = (tmpr * tmpr + tmpi * tmpi) * inv_sigma2[ii];
                            diff += tmpr;
                        }
                    }
                    else
                    {
                        for (size_t ii = 0; ii < nr_points_prob; ii++)
                        {
                            tmpr = Fimg_trans[img_start + 2*ii] - ref_scale * Fref[ref_start + 2*ii];
                            tmpi = Fimg_trans[img_start + 2*ii+1] - ref_scale * Fref[ref_start + 2*ii+1];
                            diff += (tmpr * tmpr + tmpi * tmpi) * inv_sigma2[ii];
                        }

                    }
                    // Multiply by two for t-student, because we divided by 2*sigma2 instead of sigma2!
                    if (do_student)
                        diff *= 2.;
                    dAkij(Mweight, zero_trans, refno, irot) = diff;
                    if (debug == 9)
                    {
                        std::cerr<<"refno= "<<refno<<" irot= "<<irot<<" diff= "<<diff<<" mindiff2= "<<mindiff2<<std::endl;
                    }
                    if (diff < mindiff2)
                    {
                        mindiff2 = diff;
                        opt_refno = refno;
                        opt_ipsi = ipsi;
                        opt_iflip = iflip;
                        opt_itrans = zero_trans;
                        opt_irefmir = irefmir;
                    }
                }
            }
        }
        if (debug == 9)
        {
            std::cerr<<">>>>> mindiff2= "<<mindiff2<<std::endl;
            exit(1);
        }

    }

    // Now that we have mindiff2 calculate all weights and maxweight
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            refw[refno] = 0.;
            refw_mirror[refno] = 0.;

            FOR_ALL_FLIPS()
            {
                irefmir = (int)floor(iflip / nr_nomirror_flips) * model.n_ref + refno;
                ix = (int)round(opt_offsets_ref[2*irefmir]);
                iy = (int)round(opt_offsets_ref[2*irefmir+1]);
                fracpdf = (iflip < nr_nomirror_flips) ? 1. - mirror_fraction[refno] : mirror_fraction[refno];
                pdf =  alpha_k[refno] * fracpdf * A2D_ELEM(P_phi, iy, ix);
                // get the starting point in the Fimg_trans vector
                img_start = point_trans*4*dnr_points_2d + (iflip%nr_nomirror_flips)*dnr_points_2d;

                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    diff = dAkij(Mweight, zero_trans, refno, irot);
                    // get the starting point in the Fref vector
                    ref_start = refno*nr_psi*dnr_points_2d + ipsi*dnr_points_2d;
                    if (!do_student)
                    {
                        // normal distribution
                        aux = diff - mindiff2;
                        // next line because of numerical precision of exp-function
                        if (aux > 100.)
                            weight = 0.;
                        else
                            weight = exp(-aux) * pdf;
                        // Store weight
                        dAkij(Mweight, zero_trans, refno, irot) = weight;
                    }
                    else
                    {
                        // t-student distribution
                        // pdf = (1 + diff2/df)^df2
                        // Correcting for mindiff2:
                        // pdfc = (1 + diff2/df)^df2 / (1 + mindiff2/df)^df2
                        //      = ( (1 + diff2/df)/(1 + mindiff2/df) )^df2
                        //      = ( (df + diff2) / (df + mindiff2) )^df2
                        // Extra factor two because we saved diff = | X - A |^2 / TWO * sigma
                        aux = (df + diff) / (df + mindiff2);
                        weight = pow(aux, df2) * pdf;
                        // Calculate extra weight acc. to Eq (10) Wang et al.
                        // Patt. Recognition Lett. 25, 701-710 (2004)
                        weight2 = ( df + ldim ) / ( df + diff  );
                        // Store probability weights
                        dAkij(Mweight, zero_trans, refno, irot) = weight * weight2;
                        refw2[refno] += weight * weight2;
                        sum_refw2 += weight * weight2;
                    }
                    // Accumulate sum weights
                    if (iflip < nr_nomirror_flips)
                        refw[refno] += weight;
                    else
                        refw_mirror[refno] += weight;
                    sum_refw += weight;
                    //std::cerr << "DEBUG_JM: refno << iflip << ipsi: " << refno << " " << iflip << " " << ipsi << std::endl;
                    //std::cerr << "DEBUG_JM:          weight: " << weight << " sum_refw: " << sum_refw << std::endl;

                    if (do_norm)
                    {
                        scale_numer = 0.;
                        scale_denom = 0.;
                        if (do_ctf_correction)
                        {
                            // NOTE: scale_denom could be precalculated in
                            // a much cheaper way!!!!
                            for (size_t ii = 0; ii < nr_points_prob; ii++)
                            {
                                scale_numer += Fimg_trans[img_start + 2*ii] * ctf[ii] * Fref[ref_start + 2*ii];
                                scale_numer += Fimg_trans[img_start + 2*ii+1] * ctf[ii] * Fref[ref_start + 2*ii+1];
                                scale_denom += ctf[ii] * Fref[ref_start + 2*ii] * ctf[ii] * Fref[ref_start + 2*ii];
                                scale_denom += ctf[ii] * Fref[ref_start + 2*ii+1] * ctf[ii] * Fref[ref_start + 2*ii+1];
                            }
                        }
                        else
                        {
                            for (size_t ii = 0; ii < nr_points_prob; ii++)
                            {
                                scale_numer += Fimg_trans[img_start + 2*ii] * Fref[ref_start + 2*ii];
                                scale_numer += Fimg_trans[img_start + 2*ii+1] * Fref[ref_start + 2*ii+1];
                                scale_denom += Fref[ref_start + 2*ii] * Fref[ref_start + 2*ii];
                                scale_denom += Fref[ref_start + 2*ii+1] * Fref[ref_start + 2*ii+1];
                            }
                        }
                        wsum_sc += dAkij(Mweight, zero_trans, refno, irot) * scale_numer;
                        wsum_sc2 += dAkij(Mweight, zero_trans, refno, irot) * scale_denom;
                    }
                    if (weight > Pmax_refmir[irefmir])
                        Pmax_refmir[irefmir] = weight;
                    if (weight > maxweight)
                    {
                        maxweight = weight;
                        if (do_student)
                            maxweight2 = weight2;
                        opt_refno = refno;
                        opt_ipsi = ipsi;
                        opt_iflip = iflip;
                        opt_itrans = zero_trans;
                        opt_irefmir = irefmir;
                    }
                }
            }
        }
    }

    //std::cerr << "DEBUG_JM: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<std::endl;
    // Now for all irefmir, check significant rotations...
    // and calculate their limited_translations probabilities
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            if (do_norm)
                ref_scale = opt_scale / refs_avgscale[refno];
            FOR_ALL_FLIPS()
            {
                irefmir = (int)floor(iflip / nr_nomirror_flips) * model.n_ref + refno;
                fracpdf = (iflip < nr_nomirror_flips) ? 1. - mirror_fraction[refno] : mirror_fraction[refno];

                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    // Note that for t-students distribution we now do something
                    // a little bit different compared to ML_integrate_complete!
                    // Instead of comparing "weight", we compare "weight*weight2"!!
                    if (dAkij(Mweight, zero_trans, refno, irot) > C_fast*Pmax_refmir[irefmir])
                    {
                        // get the starting point in the Fref vector
                        ref_start = refno*nr_psi*dnr_points_2d + ipsi*dnr_points_2d;
                        // expand for all limited translations
                        FOR_ALL_LIMITED_TRANSLATIONS()
                        {
                            if (itrans != zero_trans)
                            { // zero_trans has already been calculated!
                                ix = (int)round(opt_offsets_ref[2*irefmir] + Vtrans[itrans](0));
                                iy = (int)round(opt_offsets_ref[2*irefmir+1] + Vtrans[itrans](1));
                                ix = intWRAP(ix, Moffsets.startingX(), Moffsets.finishingX());
                                iy = intWRAP(iy, Moffsets.startingY(), Moffsets.finishingY());
                                if (iflip < nr_nomirror_flips)
                                    point_trans = A2D_ELEM(Moffsets, iy, ix);
                                else
                                    point_trans = A2D_ELEM(Moffsets_mirror, iy, ix);
                                if (point_trans < 0 || point_trans > dim2)
                                {
                                    std::cerr<<"point_trans = "<<point_trans<<" ix= "<<ix<<" iy= "<<iy<<std::endl;
                                    REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,"mlf_align2d BUG: point_trans < 0 or > dim2");
                                }
                                pdf =  alpha_k[refno] * fracpdf * A2D_ELEM(P_phi, iy, ix);
                                if (pdf > 0)
                                {
                                    // get the starting point in the Fimg_trans vector
                                    img_start = point_trans*4*dnr_points_2d + (iflip%nr_nomirror_flips)*dnr_points_2d;
                                    diff = 0.;
                                    if (do_ctf_correction)
                                    {
                                        for (size_t ii = 0; ii < nr_points_prob; ii++)
                                        {
                                            tmpr = Fimg_trans[img_start + 2*ii] - ctf[ii] * ref_scale * Fref[ref_start + 2*ii];
                                            tmpi = Fimg_trans[img_start + 2*ii+1] - ctf[ii] * ref_scale * Fref[ref_start + 2*ii+1];
                                            diff += (tmpr * tmpr + tmpi * tmpi) * inv_sigma2[ii];
                                        }
                                    }
                                    else
                                    {
                                        for (size_t ii = 0; ii < nr_points_prob; ii++)
                                        {
                                            tmpr = Fimg_trans[img_start + 2*ii] - ref_scale * Fref[ref_start + 2*ii];
                                            tmpi = Fimg_trans[img_start + 2*ii+1] - ref_scale * Fref[ref_start + 2*ii+1];
                                            diff += (tmpr * tmpr + tmpi * tmpi) * inv_sigma2[ii];
                                        }
                                    }
                                    if (!do_student)
                                    {
                                        // Normal distribution
                                        aux = diff - mindiff2;

                                        //FIXME: In this second loop the mindiff2 is not really the
                                        //minimum diff, that's why sometimes aux is less than zero
                                        //causing exponential to go inf, for now setting weight to zero

                                        // next line because of numerical precision of exp-function
                                        if (aux > 100. || aux < 0)
                                            weight = 0.;
                                        else
                                            weight = exp(-aux) * pdf;
                                        // Store weight
                                        dAkij(Mweight, itrans, refno, irot) = weight;
                                    }
                                    else
                                    {
                                        // t-student distribution
                                        // now diff2 and mindiff2 are already divided by sigma2!!
                                        // pdf = (1 + diff2/df)^df2
                                        // Correcting for mindiff2:
                                        // pdfc = (1 + diff2/df)^df2 / (1 + mindiff2/df)^df2
                                        //      = ( (1 + diff2/df)/(1 + mindiff2/df) )^df2
                                        //      = ( (df + diff2) / (df + mindiff2) )^df2
                                        // Multiply by two for t-student, because we divided by 2*sigma2
                                        diff *= 2.;
                                        // Extra factor two because we saved diff = | X - A |^2 / TWO * sigma
                                        aux = (df + diff) / (df + mindiff2);
                                        weight = pow(aux, df2) * pdf;
                                        // Calculate extra weight acc. to Eq (10) Wang et al.
                                        // Patt. Recognition Lett. 25, 701-710 (2004)
                                        weight2 = ( df + ldim ) / ( df + diff);
                                        // Store probability weights
                                        dAkij(Mweight, itrans, refno, irot) = weight * weight2;
                                        refw2[refno] += weight * weight2;
                                        sum_refw2 += weight * weight2;
                                    }
                                    // Accumulate sum weights
                                    if (iflip < nr_nomirror_flips)
                                        refw[refno] += weight;
                                    else
                                        refw_mirror[refno] += weight;
                                    sum_refw += weight;

                                    //std::cerr << "DEBUG_JM: LOOP2: refno << iflip << ipsi: " << refno << " " << iflip << " " << ipsi << std::endl;
                                    //std::cerr << "DEBUG_JM: LOOP2:         weight: " << weight << " sum_refw: " << sum_refw << std::endl;


                                    if (do_norm)
                                    {
                                        scale_numer = 0.;
                                        scale_denom = 0.;
                                        if (do_ctf_correction)
                                        {
                                            // NOTE: scale_denom could be precalculated in
                                            // a much cheaper way!!!!
                                            for (size_t ii = 0; ii < nr_points_prob; ii++)
                                            {
                                                scale_numer += Fimg_trans[img_start + 2*ii]
                                                               * ctf[ii] * Fref[ref_start + 2*ii];
                                                scale_numer += Fimg_trans[img_start + 2*ii+1]
                                                               * ctf[ii] * Fref[ref_start + 2*ii+1];
                                                scale_denom += ctf[ii] * Fref[ref_start + 2*ii]
                                                               * ctf[ii] * Fref[ref_start + 2*ii];
                                                scale_denom += ctf[ii] * Fref[ref_start + 2*ii+1]
                                                               * ctf[ii] * Fref[ref_start + 2*ii+1];
                                            }
                                        }
                                        else
                                        {
                                            for (size_t ii = 0; ii < nr_points_prob; ii++)
                                            {
                                                scale_numer += Fimg_trans[img_start + 2*ii]
                                                               * Fref[ref_start + 2*ii];
                                                scale_numer += Fimg_trans[img_start + 2*ii+1]
                                                               * Fref[ref_start + 2*ii+1];
                                                scale_denom += Fref[ref_start + 2*ii]
                                                               * Fref[ref_start + 2*ii];
                                                scale_denom += Fref[ref_start + 2*ii+1]
                                                               * Fref[ref_start + 2*ii+1];
                                            }
                                        }
                                        wsum_sc += dAkij(Mweight, itrans, refno, irot) * scale_numer;
                                        wsum_sc2 += dAkij(Mweight, itrans, refno, irot) * scale_denom;
                                    }
                                    if (weight > maxweight)
                                    {
                                        maxweight = weight;
                                        if (do_student)
                                            maxweight2 = weight2;
                                        opt_refno = refno;
                                        opt_ipsi = ipsi;
                                        opt_iflip = iflip;
                                        opt_itrans = itrans;
                                        opt_irefmir = irefmir;
                                    }
                                }
                                else
                                    dAkij(Mweight, itrans, refno, irot) = 0.;
                            }
                        }
                    }
                }
            }
        }
    }

    // Update optimal offsets
    opt_offsets(0) = opt_offsets_ref[2*opt_irefmir] + Vtrans[opt_itrans](0);
    opt_offsets(1) = opt_offsets_ref[2*opt_irefmir+1] + Vtrans[opt_itrans](1);
    opt_psi = -psi_step * (opt_iflip * nr_psi + opt_ipsi) - SMALLANGLE;

    // Perform KS-test on difference image of optimal hidden parameters
    // And calculate noise distribution histograms
    if (do_kstest)
    {
        MultidimArray<double> diff(2*nr_points_prob);
        std::vector<std::vector<double> > res_diff;
        res_diff.resize(hdim);
        if (opt_iflip < nr_nomirror_flips)
            point_trans = A2D_ELEM(Moffsets, (int)opt_offsets(1), (int)opt_offsets(0));
        else
            point_trans = A2D_ELEM(Moffsets_mirror, (int)opt_offsets(1), (int)opt_offsets(0));
        img_start = point_trans*4*dnr_points_2d + (opt_iflip%nr_nomirror_flips)*dnr_points_2d;
        ref_start = opt_refno*nr_psi*dnr_points_2d + opt_ipsi*dnr_points_2d;

        if (do_norm)
            ref_scale = opt_scale / refs_avgscale[opt_refno];
        for (size_t ii = 0; ii < nr_points_prob; ii++)
        {
            double inv_sqrt_sigma2 = 1. / sqrt(0.5*sigma2[ii]);
            if (do_ctf_correction)
            {
                dAi(diff,2*ii) = (Fimg_trans[img_start + 2*ii] -
                                  ctf[ii] * ref_scale * Fref[ref_start + 2*ii]) * inv_sqrt_sigma2;
                dAi(diff,2*ii+1) = (Fimg_trans[img_start + 2*ii+1] -
                                    ctf[ii] * ref_scale * Fref[ref_start + 2*ii+1]) * inv_sqrt_sigma2;
            }
            else
            {
                dAi(diff,2*ii) = (Fimg_trans[img_start + 2*ii] -
                                  ref_scale * Fref[ref_start + 2*ii])*  inv_sqrt_sigma2;
                dAi(diff,2*ii+1) = (Fimg_trans[img_start + 2*ii+1] -
                                    ref_scale * Fref[ref_start + 2*ii+1]) *  inv_sqrt_sigma2;
            }
            // Fill vectors for resolution-depedent histograms
            int ires = DIRECT_MULTIDIM_ELEM(Mresol_int, pointer_2d[ii]);
            res_diff[ires].push_back(dAi(diff,2*ii));
            res_diff[ires].push_back(dAi(diff,2*ii+1));
        }

        double * aux_array = MULTIDIM_ARRAY(diff) - 1;
        double KSD =0.;
        if (do_student)
        {
            if (df==1)
                ksone(aux_array, 2*nr_points_prob, &lcdf_tstudent_mlf1, &KSD, &KSprob);
            else if (df==3)
                ksone(aux_array, 2*nr_points_prob, &lcdf_tstudent_mlf3, &KSD, &KSprob);
            else if (df==6)
                ksone(aux_array, 2*nr_points_prob, &lcdf_tstudent_mlf6, &KSD, &KSprob);
            else if (df==9)
                ksone(aux_array, 2*nr_points_prob, &lcdf_tstudent_mlf9, &KSD, &KSprob);
            else if (df==30)
                ksone(aux_array, 2*nr_points_prob, &lcdf_tstudent_mlf30, &KSD, &KSprob);
            else
                REPORT_ERROR(ERR_VALUE_INCORRECT,"KS-test for t-distribution only implemented for df=1,3,6,9 or 30!");
        }
        else
        {
            ksone(aux_array, 2*nr_points_prob, &cdf_gauss, &KSD, &KSprob);
        }

        // Compute resolution-dependent histograms
        for (size_t ires = 0; ires < hdim; ires++)
        {
            for (size_t j=0; j<res_diff[ires].size(); j++)
                resolhist[ires].insert_value(res_diff[ires][j]);
        }

        // Overall average histogram
        if (sumhist.sampleNo()==0)
            sumhist.init(HISTMIN, HISTMAX, HISTSTEPS);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(diff)
        sumhist.insert_value(DIRECT_MULTIDIM_ELEM(diff, n));

        if (write_histograms)
        {
            Histogram1D     hist;
            double          val;
            hist.init(HISTMIN, HISTMAX, HISTSTEPS);
            compute_hist(diff,hist,HISTMIN, HISTMAX, HISTSTEPS);
            hist /= hist.sum()*hist.step_size;
            FileName fn_hist = fn_img + ".hist";
            std::ofstream fh_hist;
            fh_hist.open((fn_hist).c_str(), std::ios::out);
            if (!fh_hist)
                REPORT_ERROR(ERR_IO_NOTOPEN, (String)"Cannot write histogram file "+ fn_hist);
            FOR_ALL_ELEMENTS_IN_ARRAY1D(hist)
            {
                hist.index2val(i, val);
                val += 0.5*hist.step_size;
                fh_hist << val<<" "<<A1D_ELEM(hist, i)<<" ";
                if (do_student)
                    fh_hist << tstudent1D(val, df, 1., 0.)<<"\n";
                else
                    fh_hist << gaussian1D(val, 1., 0.)<<"\n";
            }
            fh_hist.close();
        }

    }

    // Update opt_scale
    if (do_norm)
    {
        opt_scale = wsum_sc / wsum_sc2;
    }

    // Acummulate all weighted sums
    // and normalize them by sum_refw, such that sum over all weights is one!
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            sumw[refno] += (refw[refno] + refw_mirror[refno]) / sum_refw;
            sumw2[refno] += refw2[refno] / sum_refw;
            if (do_student)
            {
                sumwsc[refno] += refw2[refno] * opt_scale / sum_refw;
                sumwsc2[refno] += refw2[refno] * opt_scale * opt_scale / sum_refw;
            }
            else
            {
                sumwsc[refno] += (refw[refno] + refw_mirror[refno]) * (opt_scale) / sum_refw;
                sumwsc2[refno] += (refw[refno] + refw_mirror[refno]) * (opt_scale * opt_scale) / sum_refw;
            }
            sumw_mirror[refno] += refw_mirror[refno] / sum_refw;
            FOR_ALL_FLIPS()
            {
                irefmir = (int)floor(iflip / nr_nomirror_flips) * model.n_ref + refno;
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    // get the starting point in the Fwsum_imgs vectors
                    wsum_start = refno*nr_psi*dnr_points_2d + ipsi*dnr_points_2d;
                    FOR_ALL_LIMITED_TRANSLATIONS()
                    {
                        ix = ROUND(opt_offsets_ref[2*irefmir] + Vtrans[itrans](0));
                        iy = ROUND(opt_offsets_ref[2*irefmir+1] + Vtrans[itrans](1));
                        ix = intWRAP(ix, Moffsets.startingX(), Moffsets.finishingX());
                        iy = intWRAP(iy, Moffsets.startingY(), Moffsets.finishingY());
                        if (iflip < nr_nomirror_flips)
                            point_trans = A2D_ELEM(Moffsets, iy, ix);
                        else
                            point_trans = A2D_ELEM(Moffsets_mirror, iy, ix);
                        weight = dAkij(Mweight, itrans, refno, irot);
                        if (weight > SIGNIFICANT_WEIGHT_LOW*maxweight)
                        {
                            weight /= sum_refw;
                            // get the starting point in the Fimg_trans vector
                            img_start = point_trans*4*dnr_points_2d + (iflip%nr_nomirror_flips)*dnr_points_2d;
                            wsum_sigma_offset += weight * (double)(ix * ix + iy * iy);
                            if (do_ctf_correction)
                            {
                                for (size_t ii = 0; ii < nr_points_2d; ii++)
                                {
                                    Fwsum_imgs[wsum_start + 2*ii] += weight
                                                                     * opt_scale * decctf[ii] * Fimg_trans[img_start + 2*ii];
                                    Fwsum_imgs[wsum_start + 2*ii+1] += weight
                                                                       * opt_scale * decctf[ii] * Fimg_trans[img_start + 2*ii+1];
                                    Fwsum_ctfimgs[wsum_start + 2*ii] += weight
                                                                        * opt_scale * Fimg_trans[img_start + 2*ii];
                                    Fwsum_ctfimgs[wsum_start + 2*ii+1] += weight
                                                                          * opt_scale * Fimg_trans[img_start + 2*ii+1];
                                    tmpr = Fimg_trans[img_start + 2*ii]
                                           - ctf[ii] * opt_scale * Fref[wsum_start + 2*ii];
                                    tmpi = Fimg_trans[img_start + 2*ii+1]
                                           - ctf[ii] * opt_scale * Fref[wsum_start + 2*ii+1];
                                    Mwsum_sigma2_local[ii] += weight * (tmpr * tmpr + tmpi * tmpi);
                                }
                            }
                            else
                            {
                                for (size_t ii = 0; ii < nr_points_2d; ii++)
                                {
                                    Fwsum_imgs[wsum_start + 2*ii] += weight
                                                                     * opt_scale * Fimg_trans[img_start + 2*ii];
                                    Fwsum_imgs[wsum_start + 2*ii+1] += weight
                                                                       * opt_scale * Fimg_trans[img_start + 2*ii+1];
                                    tmpr = Fimg_trans[img_start + 2*ii]
                                           - opt_scale * Fref[wsum_start + 2*ii];
                                    tmpi = Fimg_trans[img_start + 2*ii+1]
                                           - opt_scale * Fref[wsum_start + 2*ii+1];
                                    Mwsum_sigma2_local[ii] += weight * (tmpr * tmpr + tmpi * tmpi);
                                }
                            }
                        }// if weight > SIGNIFICANT_WEIGHT_LOW*maxweight
                    }//forall translation
                }//forall rotation
            }//forall flips
        }
    }//forall models

    // Update the optimal origin offsets
    nr_mir = (do_mirror) ? 2 : 1;

    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            for (int imir = 0; imir < nr_mir; imir++)
            {
                irefmir = imir * model.n_ref + refno;
                iflip_start = imir * nr_nomirror_flips;
                iflip_stop = imir * nr_nomirror_flips + nr_nomirror_flips;
                opt_itrans = zero_trans;
                for (int iflip = iflip_start; iflip < iflip_stop; iflip++)
                {
                    FOR_ALL_ROTATIONS()
                    {
                        irot = iflip * nr_psi + ipsi;
                        FOR_ALL_LIMITED_TRANSLATIONS()
                        {
                            weight = dAkij(Mweight, itrans, refno, irot);
                            if (weight > Pmax_refmir[irefmir])
                            {
                                Pmax_refmir[irefmir] = weight;
                                opt_itrans = itrans;
                            }
                        }
                    }
                }
                opt_offsets_ref[2*irefmir] += Vtrans[opt_itrans](0);
                opt_offsets_ref[2*irefmir+1] += Vtrans[opt_itrans](1);
            }
        }
    }

    // Distribution widths
    fracweight = maxweight / sum_refw;
    sum_refw2 /= sum_refw;


    // Compute Log Likelihood
    if (!do_student)
        // 1st term: log(refw_i)
        // 2nd term: for subtracting mindiff2
        // 3rd term: for missing normalization constant
        LL += log(sum_refw)
              - mindiff2
              - logsigma2;
    else
        // 1st term: log(refw_i)
        // 2nd term: for dividing by (1 + mindiff2/df)^df2
        // 3rd term: for sigma-dependent normalization term in t-student distribution
        // 4th&5th terms: gamma functions in t-distribution
        //LL += log(sum_refw) - log(1. + ( mindiff2 / df )) - log(sigma_noise * sigma_noise);
        LL += log(sum_refw)
              + df2 * log( 1. + ( mindiff2 / df ))
              - logsigma2
              + gammln(-df2) - gammln(df/2.);

    //    if (rank == 0)
    //    {
    //      std::cerr << "DEBUG_JM: current_image : " << current_image
    //          << " sum_refw: " << sum_refw
    //         << " LL: " << LL << std::endl;
    //    }
    delete []sigma2;
    delete []inv_sigma2;
    delete []ctf;
    delete []decctf;
}

void ProgMLF2D::iteration()
{
    // Integrate over all images
    expectation();
    // Update model with new estimates
    maximization();
}

void ProgMLF2D::expectation()
{

    Image<double> img;
    FileName fn_img, fn_trans;
    std::vector<double> allref_offsets, pdf_directions(model.n_ref);
    MultidimArray<double> trans(2);
    MultidimArray<double> opt_offsets(2);

    float old_phi = -999., old_theta = -999.;
    double opt_psi, opt_flip, opt_scale, maxcorr, maxweight2;
    double w2, KSprob = 0.;
    size_t opt_refno, opt_ipsi, opt_iflip;
    int focus = 0;
    bool apply_ctf;

    // Pre-calculate pdfs
    calculatePdfInplane();

    // Generate (FT of) each rotated version of all references
    rotateReference(Fref);

    // Initialize progress bar
    initProgress(nr_images_local);

    trans.initZeros();
    int n = model.n_ref;
    // Set all weighted sums to zero
    sumw.assign(n, 0.);
    sumw2.assign(n, 0.);
    sumwsc.assign(n, 0.);
    sumwsc2.assign(n, 0.);
    sumw_mirror.assign(n, 0.);

    sumw_defocus.clear();
    Mwsum_sigma2.clear();
    Fwsum_imgs.clear();
    Fwsum_ctfimgs.clear();
    LL = 0.;
    wsum_sigma_offset = 0.;
    sumcorr = 0.;
    std::vector<double> dum;
    dum.assign(nr_points_2d, 0.);
    Mwsum_sigma2.assign(nr_focus, dum);
    if (do_student && do_student_sigma_trick)
        sumw_defocus.assign(nr_focus, 0.);
    else
    {
        FOR_ALL_DEFOCUS_GROUPS()
        {
            sumw_defocus.push_back((double)count_defocus[ifocus]);
        }
    }

    if (dnr_points_2d > 0)
    {
        int nn = n * nr_psi * dnr_points_2d;
        Fwsum_imgs.assign(nn, 0.);
        if  (do_ctf_correction)
            Fwsum_ctfimgs.assign(nn, 0.);
    }

    std::stringstream ss;
    String s;
    // Loop over all images
    FOR_ALL_LOCAL_IMAGES()
    {
        size_t id = img_id[imgno];
        current_image = imgno;

        focus = 0;
        // Get defocus-group
        if (do_ctf_correction)
            MDimg.getValue(MDL_DEFGROUP, focus, id);

        //std::cerr << formatString("processing img: %lu with id: %lu\n", imgno, id);

        MDimg.getValue(MDL_IMAGE, fn_img, id);
        //std::cerr << formatString("   filename: %s\n", fn_img.c_str());
        //img.read(fn_img, false, false, false, false);

        img.read(fn_img);
        img().setXmippOrigin();

        // Get optimal offsets for all references
        allref_offsets = imgs_offsets[IMG_LOCAL_INDEX];

        // Read optimal orientations from memory
        if (limit_rot)
        {
            old_phi = imgs_oldphi[IMG_LOCAL_INDEX];
            old_theta = imgs_oldtheta[IMG_LOCAL_INDEX];
        }

        if (do_norm)
        {
            opt_scale = imgs_scale[IMG_LOCAL_INDEX];
        }

        // For limited orientational search: preselect relevant directions
        preselectDirections(old_phi, old_theta, pdf_directions);

        // Perform the actual expectation step for this image
        apply_ctf = !(iter == 1 && first_iter_noctf);

        processOneImage(img(), focus, apply_ctf, maxcorr, maxweight2, w2,
                        opt_scale, opt_refno, opt_psi, opt_ipsi, opt_iflip, opt_offsets,
                        allref_offsets, pdf_directions,
                        do_kstest, iter==iter_write_histograms, fn_img, KSprob);

        // for t-student, update sumw_defocus
        if (do_student && do_student_sigma_trick)
        {
            if (debug==8)
                std::cerr<<"sumw_defocus[focus]= "<<sumw_defocus[focus]<<" w2= "<<w2<<std::endl;
            sumw_defocus[focus] += w2;
        }

        // Store optimal scale in memory
        if (do_norm)
        {
            imgs_scale[IMG_LOCAL_INDEX] = opt_scale;
        }

        // Store optimal translations
        imgs_offsets[IMG_LOCAL_INDEX] = allref_offsets;

        // Store optimal phi and theta in memory
        if (limit_rot)
        {
            imgs_oldphi[IMG_LOCAL_INDEX] = model.Iref[opt_refno].rot();
            imgs_oldtheta[IMG_LOCAL_INDEX] = model.Iref[opt_refno].tilt();
        }

        // Output docfile
        sumcorr += maxcorr;
        opt_flip = 0.;
        if (-opt_psi > 360.)
        {
            opt_psi += 360.;
            opt_flip = 1.;
        }

        dAij(docfiledata,IMG_LOCAL_INDEX,0)
        = model.Iref[opt_refno].rot(); // rot
        dAij(docfiledata,IMG_LOCAL_INDEX,1)
        = model.Iref[opt_refno].tilt(); // tilt
        dAij(docfiledata,IMG_LOCAL_INDEX,2) = opt_psi + 360.; // psi
        dAij(docfiledata,IMG_LOCAL_INDEX,3) = trans(0) + opt_offsets(0); // Xoff
        dAij(docfiledata,IMG_LOCAL_INDEX,4) = trans(1) + opt_offsets(1); // Yoff
        dAij(docfiledata,IMG_LOCAL_INDEX,5) = (double) (opt_refno + 1); // Ref
        dAij(docfiledata,IMG_LOCAL_INDEX,6) = opt_flip; // Mirror
        dAij(docfiledata,IMG_LOCAL_INDEX,7) = maxcorr; // P_max/P_tot

        if (do_student)
        {
            dAij(docfiledata,IMG_LOCAL_INDEX,8) = maxweight2; // Robustness weight
        }
        if (do_norm)
        {
            dAij(docfiledata,IMG_LOCAL_INDEX,9) = opt_scale; // image scale
        }
        if (do_kstest)
        {
            dAij(docfiledata,IMG_LOCAL_INDEX,10) = KSprob;
        }

        // Report progress bar
        setProgress();
    }

    endProgress();

    if (do_ctf_correction)
    {
        reverseRotateReference(Fwsum_ctfimgs,wsum_ctfMref);
    }
    reverseRotateReference(Fwsum_imgs,wsum_Mref);

}

// Update all model parameters
void ProgMLF2D::maximization()
{

    MultidimArray<double> rmean_sigma2, rmean_signal2;
    MultidimArray<int>  radial_count;
    Matrix1D<int> center(2);
    MultidimArray<std::complex<double> > Faux, Faux2;
    MultidimArray<double> Maux;
    FileName fn_tmp;
    double aux;
    int c;

    // Pre-calculate sumw_allrefs & average Pmax/sumP or cross-correlation
    sumw_allrefs = 0.;

    // Update the reference images
    FOR_ALL_MODELS()
    {
        if (!do_student && sumw[refno] > 0.)
        {
            model.Iref[refno]() = wsum_Mref[refno];
            model.Iref[refno]() /= sumwsc2[refno];
            model.Iref[refno].setWeight(sumw[refno]);
            sumw_allrefs += sumw[refno];
            if (do_ctf_correction)
            {
                Ictf[refno]() = wsum_ctfMref[refno];
                Ictf[refno]() /= sumwsc2[refno];
                Ictf[refno].setWeight(sumw[refno]);
            }
            else
            {
                Ictf[refno]=model.Iref[refno];
            }
        }
        else if (do_student && sumw2[refno] > 0.)
        {
            model.Iref[refno]() = wsum_Mref[refno];
            model.Iref[refno]() /= sumwsc2[refno];
            model.Iref[refno].setWeight(sumw2[refno]);
            sumw_allrefs += sumw[refno];
            //sumw_allrefs2 += sumw2[refno];
            if (do_ctf_correction)
            {
                Ictf[refno]() = wsum_ctfMref[refno];
                Ictf[refno]() /= sumwsc2[refno];
                Ictf[refno].setWeight(sumw2[refno]);
            }
            else
            {
                Ictf[refno]=model.Iref[refno];
            }
        }
        else
        {
            model.Iref[refno].setWeight(0.);
            Ictf[refno].setWeight(0.);
            model.Iref[refno]().initZeros(dim, dim);
            Ictf[refno]().initZeros();
        }
    }

    // Adjust average scale (nr_classes will be smaller than n_ref for the 3D case!)
    if (do_norm)
    {
        int iclass, nr_classes = ROUND(model.n_ref / refs_per_class);
        std::vector<double> wsum_scale(nr_classes), sumw_scale(nr_classes);
        ldiv_t temp;
        average_scale = 0.;
        FOR_ALL_MODELS()
        {
            average_scale += sumwsc[refno];
            temp = ldiv( refno, refs_per_class );
            iclass = ROUND(temp.quot);
            wsum_scale[iclass] += sumwsc[refno];
            sumw_scale[iclass] += sumw[refno];
        }
        FOR_ALL_MODELS()
        {
            temp = ldiv( refno, refs_per_class );
            iclass = ROUND(temp.quot);
            if (sumw_scale[iclass]>0.)
            {
                refs_avgscale[refno] = wsum_scale[iclass]/sumw_scale[iclass];
                model.Iref[refno]() *= refs_avgscale[refno];
                Ictf[refno]() *= refs_avgscale[refno];
            }
            else
            {
                refs_avgscale[refno] = 1.;
            }
        }
        average_scale /= sumw_allrefs;
    }

    // Average corr
    sumcorr /= sumw_allrefs;

    // Update the model fractions
    if (!fix_fractions)
    {
        FOR_ALL_MODELS()
        {
            if (sumw[refno] > 0.)
            {
                alpha_k[refno] = sumw[refno] / sumw_allrefs;
                mirror_fraction[refno] = sumw_mirror[refno] / sumw[refno];
            }
            else
            {
                alpha_k[refno] = 0.;
                mirror_fraction[refno] = 0.;
            }
        }
    }

    // Update sigma of the origin offsets
    if (!fix_sigma_offset)
    {
        sigma_offset = sqrt(wsum_sigma_offset / (2. * sumw_allrefs));
    }

    // Update the noise parameters
    if (!fix_sigma_noise)
    {
        FOR_ALL_DEFOCUS_GROUPS()
        {
            getFTfromVector(Mwsum_sigma2[ifocus],0,Faux,true);
            Half2Whole(Faux, Faux2, dim);
            FFT_magnitude(Faux2, Maux);
            CenterFFT(Maux, true);
            Maux.setXmippOrigin();
            center.initZeros();
            rmean_sigma2.initZeros();
            radialAverage(Maux, center, rmean_sigma2, radial_count, true);
            // Factor 2 here, because the Gaussian distribution is 2D!
            for (size_t irr = 0; irr <= current_highres_limit; irr++)
            {
                aux = dAi(rmean_sigma2, irr) / (2. * sumw_defocus[ifocus]);
                if (aux > 0.)
                {
                    VSIG_ITEM = aux;
                }
            }
        }
    }

    // Calculate average spectral signal
    c = 0;
    FOR_ALL_MODELS()
    {
        if ( (!do_student && sumw[refno] > 0.) ||
             ( do_student && sumw2[refno] > 0.) )
        {
            FourierTransform(Ictf[refno](), Faux);
            FFT_magnitude(Faux, Maux);
            CenterFFT(Maux, true);
            Maux *= Maux;
            //if (!do_student)
            //Maux *= sumw[refno];
            //else
            //Maux *= sumw2[refno];
            Maux *= sumw[refno];
            center.initZeros();
            rmean_signal2.initZeros();
            Maux.setXmippOrigin();
            radialAverage(Maux, center, rmean_signal2, radial_count, true);
            if (c == 0)
                spectral_signal = rmean_signal2;
            else
                spectral_signal += rmean_signal2;
            c++;
            //            std::cerr << "DEBUG_JM: ====================" <<std::endl;
            //            std::cerr << "DEBUG_JM: refno: " << refno << std::endl;
            //            std::cerr << "DEBUG_JM: rmean_signal2: " << rmean_signal2 << std::endl;
            //            std::cerr << "DEBUG_JM: spectral_signal: " << spectral_signal << std::endl;
        }
    }
    double inv_nref = 1./(double)model.n_ref;
    spectral_signal *= inv_nref;
    //    std::cerr << "DEBUG_JM: inv_nref: " << inv_nref << std::endl;
    //    std::cerr << "DEBUG_JM, maximization: spectral_signal: " << spectral_signal << std::endl;
}

void ProgMLF2D::endIteration()
{
    ML2DBaseProgram::endIteration();
    //std::cerr << "DEBUG_JM, endIteration: spectral_signal: " << spectral_signal << std::endl;
    updateWienerFilters(spectral_signal, sumw_defocus, iter);
}

void writeImage(Image<double> &img, const FileName &fn, MDRow &row)
{
    img.write(fn);
    row.setValue(MDL_IMAGE, fn);
    row.setValue(MDL_ENABLED, 1);
    double weight = img.weight();
    double rot = 0., tilt = 0., psi = 0.;
    img.getEulerAngles(rot, tilt, psi);
    row.setValue(MDL_WEIGHT, weight);
    size_t count = ROUND(weight);
    row.setValue(MDL_CLASS_COUNT, count);
    //    row.setValue(MDL_ANGLE_ROT, rot);
    //    row.setValue(MDL_ANGLE_TILT, tilt);
    //    row.setValue(MDL_ANGLE_PSI, psi);
}

void ProgMLF2D::writeOutputFiles(const ModelML2D &model, OutputType outputType)
{

    FileName          fn_tmp, fn_base;
    Image<double>        Itmp;
    MetaData          MDo;
    String       comment;
    std::ofstream     fh;
    size_t id;
    bool write_conv = !do_ML3D;
    //Only override first time

    if (iter == 0)
        write_conv = false;

    fn_base = (outputType == OUT_ITER || outputType == OUT_REFS) ?
              getIterExtraPath(fn_root, iter) : fn_root;
    //    {
    //        fn_base = getIterExtraPath(fn_root, iter);
    //        fn_prefix = FN_ITER_PREFIX(iter);
    //
    //    }

    // Write out optimal image orientations
    fn_tmp = FN_IMAGES_MD(fn_base);
    MDimg.write(fn_tmp);

    // Write out current reference images and fill sel & log-file
    // First time for _ref, second time for _cref
    FileName fn_base_cref = FN_CREF_IMG;
    MDRow row;
    MDIterator mdIter(MDref);
    MDo = MDref;

    for (int refno = 0; refno < model.n_ref; ++refno, mdIter.moveNext())
    {
        //row.setValue(MDL_ITER, iter);
        row.setValue(MDL_REF, refno + 1);

        if (do_mirror)
            row.setValue(MDL_MIRRORFRAC, mirror_fraction[refno]);

        if (write_conv)
            row.setValue(MDL_SIGNALCHANGE, conv[refno]*1000);

        if (do_norm)
            row.setValue(MDL_INTSCALE, refs_avgscale[refno]);

        //write ctf
        fn_tmp.compose(refno + 1, fn_base_cref);
        writeImage(Ictf[refno], fn_tmp, row);
        MDo.setRow(row, mdIter.objId);

        //write image
        fn_tmp = FN_REF(fn_base, refno + 1);
        Itmp = model.Iref[refno];
        writeImage(Itmp, fn_tmp, row);
        MDref.setRow(row, mdIter.objId);
    }

    // Write out reference md file
    MDo.write(FN_CREF_IMG_MD);
    outRefsMd = FN_CLASSES_MD(fn_base);
    MDref.write(outRefsMd);

    // Write out log-file
    MetaData mdLog;
    //mdLog.setComment(cline);
    id = mdLog.addObject();
    mdLog.setValue(MDL_ITER, iter, id);
    mdLog.setValue(MDL_LL, LL, id);
    mdLog.setValue(MDL_PMAX, sumcorr, id);
    mdLog.setValue(MDL_SIGMAOFFSET, sigma_offset, id);
    mdLog.setValue(MDL_RANDOMSEED, seed, id);
    if (do_norm)
        mdLog.setValue(MDL_INTSCALE, average_scale, id);
    //    fn_tmp = FN_IMGMD(fn_prefix);
    //    mdLog.setValue(MDL_IMGMD, fn_tmp, id);
    //    mdLog.setValue(MDL_REFMD, outRefsMd, id);
    //    if (outputType == OUT_ITER)
    //      fn_prefix = fn_root + "_iter";

    mdLog.write(FN_LOGS_MD(fn_base), MD_APPEND);

    if (outputType != OUT_REFS)
    {
        MetaData mdImgs;
        int n = (int) MDref.size(); //Avoid int and size_t comparison warning
        for (int ref = 1; ref <= n; ++ref)
        {
            mdImgs.importObjects(MDimg, MDValueEQ(MDL_REF, ref));
            mdImgs.write(FN_CLASS_IMAGES_MD(fn_base, ref), MD_APPEND);
        }
    }
    //    fn_tmp = FN_LOGMD(fn_prefix);
    //    if (mode == MD_OVERWRITE)
    //      mdLog.write(fn_tmp);
    //    else
    //      mdLog.append(fn_tmp);
    //
    //    mode = MD_APPEND;


    // Write out average and resolution-dependent histograms
    if (do_kstest)
    {
        double          val;
        std::ofstream fh_hist;

        fn_tmp = fn_base + "avg.hist.log";
        fh_hist.open((fn_tmp).c_str(), std::ios::out);
        if (!fh_hist)
            REPORT_ERROR(ERR_IO_NOTOPEN, (String)"Cannot write histogram file "+ fn_tmp);
        sumhist /= (sumhist.sum()*sumhist.step_size);
        FOR_ALL_ELEMENTS_IN_ARRAY1D(sumhist)
        {
            sumhist.index2val(i, val);
            val += 0.5*sumhist.step_size;
            fh_hist << val<<" "<< A1D_ELEM(sumhist, i)<<" ";
            if (do_student)
                fh_hist << tstudent1D(val, df, 1., 0.)<<"\n";
            else
                fh_hist << gaussian1D(val, 1., 0.)<<"\n";
        }
        fh_hist.close();

        fn_tmp = fn_base + "_resol.hist";
        fh_hist.open((fn_tmp).c_str(), std::ios::out);
        if (!fh_hist)
            REPORT_ERROR(ERR_IO_NOTOPEN, (String)"Cannot write histogram file "+ fn_tmp);
        FOR_ALL_ELEMENTS_IN_ARRAY1D(sumhist)
        {
            sumhist.index2val(i, val);
            val += 0.5*sumhist.step_size;
            fh_hist << val<<" ";
            if (do_student)
                fh_hist << tstudent1D(val, df, 1., 0.)<<" ";
            else
                fh_hist << gaussian1D(val, 1., 0.)<<" ";
            for (size_t ires = 0; ires < hdim; ires++)
            {
                if (resolhist[ires].sampleNo() > 0)
                {
                    fh_hist <<A1D_ELEM(resolhist[ires], i)/(resolhist[ires].sum()*resolhist[ires].step_size)<<" ";
                }
            }
            fh_hist <<"\n";
        }
        fh_hist.close();

    }

    // Write out updated Vsig vectors
    if (!fix_sigma_noise)
    {
        FOR_ALL_DEFOCUS_GROUPS()
        {
            writeNoiseFile(fn_base, ifocus);
        }
    }
}//function writeOutputFiles

void ProgMLF2D::writeNoiseFile(const FileName &fn_base, int ifocus)
{
    LOG_LEVEL(writeNoiseFile);
    MetaData md;
    // Calculate resolution step
    double resol_step = 1. / (sampling * dim);
    md.addLabel(MDL_RESOLUTION_FREQ);

    FOR_ALL_DIGITAL_FREQS()
    md.setValue(MDL_MLF_NOISE, VSIG_ITEM, md.addObject());

    md.fillLinear(MDL_RESOLUTION_FREQ, 0., resol_step);
    // write Vsig vector to disc
    md.write(FN_VSIG(fn_base, ifocus, "noise.xmd"), (ifocus > 0 ? MD_APPEND : MD_OVERWRITE));

}//function writeNoiseFile

/// Add docfiledata to docfile
void ProgMLF2D::addPartialDocfileData(const MultidimArray<double> &data,
                                      size_t first, size_t last)
{
    LOG_LEVEL(addPartialDocfile);
    for (size_t imgno = first; imgno <= last; imgno++)
    {
        size_t index = imgno - first;
        size_t id = img_id[imgno];
        //FIXME now directly to MDimg
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
        //MDimg.setValue(MDL_ANGLE_PSI, dAij(data, index, 2), id);
        MDimg.setValue(MDL_SHIFT_X, dAij(data, index, 3), id);
        MDimg.setValue(MDL_SHIFT_Y, dAij(data, index, 4), id);
        MDimg.setValue(MDL_REF, ROUND(dAij(data, index, 5)), id);
        if (do_mirror)
        {
            MDimg.setValue(MDL_FLIP, dAij(data, index, 6) != 0., id);
        }
        MDimg.setValue(MDL_PMAX, dAij(data, index, 7), id);
        if (do_student)
        {
            MDimg.setValue(MDL_WROBUST, dAij(data, index, 8), id);
        }
        if (do_norm)
        {
            MDimg.setValue(MDL_INTSCALE, dAij(data, index, 9), id);
        }
        if (do_kstest)
        {
            MDimg.setValue(MDL_KSTEST, dAij(data, index, 10), id);
        }

    }
}//close function addDocfileData

