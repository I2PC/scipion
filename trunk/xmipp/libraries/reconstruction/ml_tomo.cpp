/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2007)
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
#include "ml_tomo.h"

//#define DEBUG
// For blocking of threads
pthread_mutex_t mltomo_weightedsum_update_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mltomo_selfile_access_mutex = PTHREAD_MUTEX_INITIALIZER;

// Read arguments ==========================================================
void Prog_ml_tomo_prm::read(int argc, char **argv, bool ML3D)
{
    // Generate new command line for restart procedure
    cline = "";
    int argc2 = 0;
    char ** argv2 = NULL;

    if (checkParameter(argc, argv, "-restart"))
    {
        std::string comment;
        FileName fn_sel;
        DocFile DFi;
        DFi.read(getParameter(argc, argv, "-restart"));
        DFi.go_beginning();
        comment = (DFi.get_current_line()).get_text();
        if (strstr(comment.c_str(), "ml_tomo-logfile") == NULL)
        {
            std::cerr << "Error!! Docfile is not of ml_tomo-logfile type. " << std::endl;
            exit(1);
        }
        else
        {
            char *copy;
            int n = 0;
            int nmax = DFi.dataLineNo();
            SFr.reserve(nmax);
            copy = NULL;
            DFi.next();
            comment = " -frac " + DFi.name();
            if (!ML3D)
            {
                fn_sel = DFi.name();
                fn_sel = fn_sel.without_extension() + "_restart.sel";
                comment += " -ref " + fn_sel;
            }
            comment += (DFi.get_current_line()).get_text();
            DFi.next();
            cline = (DFi.get_current_line()).get_text();
            comment = comment + cline;
            // regenerate command line
            generateCommandLine(comment, argc2, argv2, copy);
            if (!ML3D)
            {
                // Read images names from restart file
                DFi.next();
                while (n < nmax)
                {
                    n++;
                    DFi.next();
                    if (DFi.get_current_line().Is_comment()) fn_sel = ((DFi.get_current_line()).get_text()).erase(0, 3);
                    SFr.insert(fn_sel, SelLine::ACTIVE);
                    DFi.adjust_to_data_line();
                }
                fn_sel = DFi.name();
                fn_sel = fn_sel.without_extension() + "_restart.sel";
                SFr.write(fn_sel);
                SFr.clear();
            }
        }
    }
    else
    {
	// no restart, just copy argc to argc2 and argv to argv2
	argc2 = argc;
	argv2 = argv;	
        for (int i = 1; i < argc2; i++)
        {
            cline = cline + (std::string)argv2[i] + " ";
        }
    }

    // Read command line
    if (checkParameter(argc2, argv2, "-more_options"))
    {
	usage();
	extendedUsage();
    }
    nr_ref = textToInteger(getParameter(argc2, argv2, "-nref", "0"));
    fn_ref = getParameter(argc2, argv2, "-ref", "");
    fn_doc = getParameter(argc2, argv2, "-doc", "");
    fn_sel = getParameter(argc2, argv2, "-i");
    fn_root = getParameter(argc2, argv2, "-o", "mltomo");
    fn_sym = getParameter(argc2, argv2, "-sym", "c1");
    Niter = textToInteger(getParameter(argc2, argv2, "-iter", "100"));
    istart = textToInteger(getParameter(argc2, argv2, "-istart", "1"));
    sigma_noise = textToFloat(getParameter(argc2, argv2, "-noise", "1"));
    sigma_offset = textToFloat(getParameter(argc2, argv2, "-offset", "3"));
    fn_frac = getParameter(argc2, argv2, "-frac", "");
    fix_fractions = checkParameter(argc2, argv2, "-fix_fractions");
    fix_sigma_offset = checkParameter(argc2, argv2, "-fix_sigma_offset");
    fix_sigma_noise = checkParameter(argc2, argv2, "-fix_sigma_noise");
    eps = textToFloat(getParameter(argc2, argv2, "-eps", "5e-5"));
    verb = textToInteger(getParameter(argc2, argv2, "-verb", "1"));
    no_SMALLANGLE = checkParameter(argc2, argv2, "-no_SMALLANGLE");

    // Normalization 
    do_norm = checkParameter(argc2, argv2, "-norm");

    // regularization
    reg0=textToFloat(getParameter(argc2, argv2, "-reg0", "0."));
    regF=textToFloat(getParameter(argc2, argv2, "-regF", "0."));
    reg_step=textToFloat(getParameter(argc2, argv2, "-reg_step", "1."));

    // ML (with/without) imputation, or maxCC
    do_ml = !checkParameter(argc2, argv2, "-maxCC");
    do_impute = !checkParameter(argc2, argv2, "-dont_impute");
    noimp_threshold = textToFloat(getParameter(argc2, argv2, "-noimp_threshold", "1."));

    // Angular sampling
    angular_sampling = textToFloat(getParameter(argc2, argv2, "-ang", "10"));
    psi_sampling = textToFloat(getParameter(argc2, argv2, "-psi_sampling", "-1"));
    tilt_range0 = textToFloat(getParameter(argc2, argv2, "-tilt0", "-91."));
    tilt_rangeF = textToFloat(getParameter(argc2, argv2, "-tiltF", "91."));
    ang_search = textToFloat(getParameter(argc2, argv2, "-ang_search", "-1."));

    // Missing data structures 
    fn_missing = getParameter(argc2, argv2, "-missing","");
    max_resol = textToFloat(getParameter(argc2, argv2, "-max_resol", "0.5"));

    // Hidden arguments
    do_esthetics = checkParameter(argc2, argv2, "-esthetics");
    trymindiff_factor = textToFloat(getParameter(argc2, argv2, "-trymindiff_factor", "0.9"));

    // Number of threads
    threads = textToInteger(getParameter(argc2, argv2, "-thr","1"));

}

// Show ====================================================================
void Prog_ml_tomo_prm::show(bool ML3D)
{

    if (verb > 0)
    {
        // To screen
        std::cerr << " -----------------------------------------------------------------" << std::endl;
        std::cerr << " | Read more about this program in the following publication:    |" << std::endl;
        std::cerr << " |  Scheres ea. in preparation                                   |" << std::endl;
        std::cerr << " |                                                               |" << std::endl;
        std::cerr << " |   *** Please cite it if this program is of use to you! ***    |" << std::endl;
        std::cerr << " -----------------------------------------------------------------" << std::endl;
        std::cerr << "--> Maximum-likelihood multi-reference refinement " << std::endl;
        std::cerr << "  Input images            : " << fn_sel << " (" << nr_exp_images << ")" << std::endl;
        if (fn_ref != "")
            std::cerr << "  Reference image(s)      : " << fn_ref << std::endl;
        else
            std::cerr << "  Number of references:   : " << nr_ref << std::endl;
        std::cerr << "  Output rootname         : " << fn_root << std::endl;
        std::cerr << "  Angular sampling rate   : " << angular_sampling<< "degrees"<<std::endl;
        if (ang_search > 0.)
        {
            std::cerr << "  Local angular searches  : "<<ang_search<<" degrees"<<std::endl;
        }
        std::cerr << "  Symmetry group          : " << fn_sym<<std::endl;
        std::cerr << "  Stopping criterium      : " << eps << std::endl;
        std::cerr << "  initial sigma noise     : " << sigma_noise << std::endl;
        std::cerr << "  initial sigma offset    : " << sigma_offset << std::endl;
        if (reg0 > 0.)
        {
            std::cerr << "  Regularization from     : "<<reg0<<" to "<<regF<<std::endl;
            std::cerr << "  Regularization steps    : "<<reg_step<<std::endl;
        }
        if (fn_missing!="")
        {
            std::cerr << "  Missing data info       : "<<fn_missing <<std::endl;
            if (do_impute)
                std::cerr << "  Missing data treatment  : imputation "<<std::endl;
            else
                std::cerr << "  Missing data treatment  : conventional division "<<std::endl;
            std::cerr << "  Maximum resolution      : " << max_resol << std::endl;
        }
        if (fn_frac != "")
            std::cerr << "  Initial model fractions : " << fn_frac << std::endl;

        if (!do_ml) 
        {
            std::cerr << "  -> Use constrained correlation coefficient instead of ML-approach." << std::endl;
        }
        if (fix_fractions)
        {
            std::cerr << "  -> Do not update estimates of model fractions." << std::endl;
        }
        if (fix_sigma_offset)
        {
            std::cerr << "  -> Do not update sigma-estimate of origin offsets." << std::endl;
        }
        if (fix_sigma_noise)
        {
            std::cerr << "  -> Do not update sigma-estimate of noise." << std::endl;
        }
        if (do_norm)
        {
            std::cerr << "  -> Refine normalization for each experimental image"<<std::endl;
        }
        if (threads>1)
        {
            std::cerr << "  -> Using "<<threads<<" parallel threads"<<std::endl;
        }
        if (do_esthetics)
        {
            std::cerr << "  -> Perform esthetics on (0,0)-pixel artifacts" << std::endl;
        }
  

        std::cerr << " -----------------------------------------------------------------" << std::endl;

    }

}

// Usage ===================================================================
void Prog_ml_tomo_prm::usage()
{
    //TODO!!
    std::cerr << "Usage:  ml_tomo [options] " << std::endl;
    std::cerr << "   -i <selfile>                : Selfile with input images \n";
    std::cerr << "   -nref <int>                 : Number of references to generate automatically (recommended)\n";
    std::cerr << "   OR -ref <selfile/image>         OR selfile with initial references/single reference image \n";
    std::cerr << " [ -o <rootname> ]             : Output rootname (default = \"ml2d\")\n";
    std::cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_ml_tomo_prm::extendedUsage(bool ML3D)
{
    std::cerr << "Additional options: " << std::endl;
    std::cerr << " [ -eps <float=5e-5> ]         : Stopping criterium \n";
    std::cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
    std::cerr << " [ -noise <float=1> ]          : Expected standard deviation for pixel noise \n";
    std::cerr << " [ -offset <float=3> ]         : Expected standard deviation for origin offset [pix]\n";
    std::cerr << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n";
    if (!ML3D) std::cerr << " [ -restart <logfile> ]        : restart a run with all parameters as in the logfile \n";
    if (!ML3D) std::cerr << " [ -istart <int> ]             : number of initial iteration \n";
    std::cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise \n";
    std::cerr << " [ -fix_sigma_offset]          : Do not re-estimate the standard deviation in the origin offsets \n";
    std::cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
    std::cerr << " [ -doc <docfile=\"\"> ]         : Read initial angles and offsets from docfile \n";
    std::cerr << " [ -norm ]                     : Refined normalization parameters for each particle \n";
    std::cerr << std::endl;
    exit(1);
}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produce_Side_info2)
void Prog_ml_tomo_prm::produceSideInfo()
{

    FileName                    fn_img, fn_tmp, fn_base, fn_tmp2;
    VolumeXmipp                 vol;
    SelLine                     SL;
    SelFile                     SFtmp, SFpart;
    Matrix1D<double>            offsets(3), dum;
    Matrix3D<double>            Maux, Maux2;
    Matrix3D<std::complex<double> >  Faux;
    Matrix1D<int>               center(3), radial_count;
    double                      av, aux, Q0;
    int                         im, jm;

#ifdef  DEBUG
    std::cerr<<"Start produceSideInfo"<<std::endl;
#endif

    // Read selfile with experimental images
    SF.read(fn_sel);

    // Get image sizes and total number of images
    SF.go_beginning();
    vol.read(SF.NextImg());
    dim = XSIZE(vol());
    hdim = dim / 2;
    dim3 = dim * dim * dim;
    ddim3 = (double)dim3;
    nr_exp_images = SF.ImgNo();
    if (YSIZE(vol()) != dim || ZSIZE(vol()) != dim)
        REPORT_ERROR(1,"ml_tomo ERROR%: only cubic volumes are allowed");

    if (regF > reg0)
        REPORT_ERROR(1,"regF should be smaller than reg0!");
    reg_current = reg0;

    // Make fourier and real-space masks
    Matrix3D<int> int_mask(dim,dim,dim);
    int_mask.setXmippOrigin();
    real_mask.resize(dim, dim, dim);
    real_mask.setXmippOrigin();
    BinarySphericalMask(int_mask, hdim, INNER_MASK);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(int_mask)
    {
        DIRECT_MULTIDIM_ELEM(real_mask,n) = (double)DIRECT_MULTIDIM_ELEM(int_mask,n);
    }

    // Set-up fourier-space mask
    Matrix3D<double> cosine_mask(dim,dim,dim);
    cosine_mask.setXmippOrigin();
    fourier_mask.resize(dim,dim,hdim+1);
    int r_max_resol = XMIPP_MIN(hdim, FLOOR(max_resol * dim));
    RaisedCosineMask(cosine_mask, r_max_resol - 2, r_max_resol, INNER_MASK);
    // exclude origin voxel
    VOL_ELEM(cosine_mask,0,0,0) = 0.;
    int yy, zz;
    int dimb = (dim - 1)/2;
    for ( int z=0, ii=0; z<dim; z++ ) 
    {
        if ( z > dimb ) zz = z-dim;
        else zz = z;
        for ( int y=0; y<dim; y++ ) 
        {
            if ( y > dimb ) yy = y-dim;
            else yy = y;
            for ( int xx=0; xx<hdim + 1; xx++, ii++ ) 
            {
                if (xx <= FINISHINGX(cosine_mask))
                {
                    DIRECT_MULTIDIM_ELEM(fourier_mask,ii) = VOL_ELEM(cosine_mask,xx,yy,zz);
                }
            }
        }
    }

    // Set-up real-space mask
    real_mask.resize(dim,dim,dim);
    real_mask.setXmippOrigin();
    RaisedCosineMask(real_mask, hdim - 1, hdim, INNER_MASK);
    real_omask.resize(real_mask);
    real_omask = 1. - real_mask;

    // Get number of references
    if (fn_ref != "")
    {
        do_generate_refs = false;
        if (Is_VolumeXmipp(fn_ref)) nr_ref = 1;
        else
        {
            SFr.read(fn_ref);
            nr_ref = SFr.ImgNo();
        }
    }
    else
        do_generate_refs = true;

    // Precalculate sampling
    mysampling.SetSampling(angular_sampling);
    if (!mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
        REPORT_ERROR(3005, (std::string)"ml_refine3d::run Invalid symmetry" +  fn_sym);
    mysampling.SL.read_sym_file(fn_sym);
    mysampling.fill_L_R_repository();
    // by default max_tilt= +91., min_tilt= -91.
    mysampling.Compute_sampling_points(false, // half sphere?
                                       tilt_rangeF,
                                       tilt_range0);
    mysampling.remove_redundant_points_exhaustive(symmetry, 
                                                  sym_order, 
                                                  false, // half sphere?
                                                  0.75 * angular_sampling);
    if (psi_sampling < 0)
        psi_sampling = angular_sampling;

    // Read in docfile with information about the missing wedges
    nr_miss = 0;
    all_missing_info.clear();
    if (fn_missing=="")
    {
        do_missing = false;
    }
    else
    {
        do_missing = true;
        missing_info myinfo;
        DocFile DFm;
        FileName fn_tst;
        DFm.read(fn_missing);
        DFm.go_beginning();
        if (DFm.get_current_line().Is_comment()) 
        {
            fn_tst = (DFm.get_current_line()).get_text();
            if (strstr(fn_tst.c_str(), "Wedgeinfo") == NULL)
                REPORT_ERROR(1,"Missing wedge info file does not have a header starting with:\n ; Wedgeinfo");
            
            int n = 0;
            int nmax = DFm.dataLineNo();
            while (n < nmax)
            {
                n++;
                DFm.next();
                if (DFm.get_current_line().Is_comment()) 
                    fn_tst = (DFm.get_current_line()).get_text();
                DFm.adjust_to_data_line();

                if (strstr(fn_tst.c_str(), "wedge_y") != NULL)
                {
                    myinfo.type=MISSING_WEDGE_Y;
                    myinfo.thy0=DFm(0);
                    myinfo.thyF=DFm(1);
                    myinfo.thx0=-90.;
                    myinfo.thxF=90.;
                }
                else if (strstr(fn_tst.c_str(), "wedge_x") != NULL)
                {
                    myinfo.type=MISSING_WEDGE_X;
                    myinfo.thx0=DFm(0);
                    myinfo.thxF=DFm(1);
                    myinfo.thy0=-90.;
                    myinfo.thyF=90.;
                }
                else if (strstr(fn_tst.c_str(), "pyramid") != NULL)
                {
                    myinfo.type=MISSING_PYRAMID;
                    myinfo.thy0=DFm(0);
                    myinfo.thyF=DFm(1);
                    myinfo.thx0=DFm(2);
                    myinfo.thxF=DFm(3);
                }
                else if (strstr(fn_tst.c_str(), "cone") != NULL)
                {
                    myinfo.type=MISSING_CONE;
                    myinfo.thy0=DFm(0);
                    myinfo.thyF=0.;
                    myinfo.thx0=0.;
                    myinfo.thxF=0.;
                }
                else
                {
                    std::cerr<<" Missing type= "<<fn_tst<<std::endl;
                    REPORT_ERROR(1,"ERROR: unrecognized missing type!");
                }
                all_missing_info.push_back(myinfo);
                nr_miss++;
            }
        }
    }


    // Fill average scales 
    if (do_norm)
    {
        average_scale = 1.;
        refs_avgscale.clear();
        for (int refno=0;refno<nr_ref; refno++)
        {
            refs_avgscale.push_back(1.);
        }
    }

    // read in model fractions if given on command line
    if (fn_frac != "")
    {
        double sumfrac = 0.;
        DocFile DF;
        DocLine DL;
        DF.read(fn_frac);
        DF.go_first_data_line();
        for (int refno = 0; refno < nr_ref; refno++)
        {
            DL = DF.get_current_line();
            alpha_k(refno) = DL[0];
            if (do_norm)
            {
                refs_avgscale[refno] = DL[3];
            }
            sumfrac += alpha_k(refno);
            DF.next_data_line();
        }
        if (ABS(sumfrac - 1.) > 1e-3)
            if (verb > 0) std::cerr << " ->WARNING: Sum of all expected model fractions (" << sumfrac << ") is not one!" << std::endl;
        for (int refno = 0; refno < nr_ref; refno++)
        {
            alpha_k(refno) /= sumfrac;
        }
    }

}

// Generate initial references =============================================
void Prog_ml_tomo_prm::generateInitialReferences()
{

    SelFile SFtmp;
    VolumeXmipp Iave, Itmp;
    Matrix3D<double> Msumwedge, Mmissing;
    Matrix3D<std::complex<double> > Fave;
    double dummy;
    FileName fn_tmp;
    SelLine line;
    DocFile DF;

#ifdef  DEBUG
    std::cerr<<"Start generateInitialReferences"<<std::endl;
#endif

    if (verb > 0)
    {
        std::cerr << "  Generating initial references by averaging over random subsets" << std::endl;
        init_progress_bar(SF.ImgNo());
    }

    if (fn_doc!="") {
        DF.read(fn_doc);
    }
    // Make random subsets and calculate average images in random orientations
    // and correct using conventional division by the sum of the wedges
    randomize_random_generator();
    SFtmp = SF.randomize();
    int Nsub = ROUND((double)SFtmp.ImgNo() / nr_ref);
    Matrix2D<double> random_A;
    double random_rot, random_tilt, random_psi;
    int c = 0, cc = XMIPP_MAX(1, SFtmp.ImgNo() / 60);
    for (int refno = 0; refno < nr_ref; refno++)
    {
        SFtmp.go_beginning();
        SFtmp.jump_lines(Nsub*refno);
        if (refno == nr_ref - 1) Nsub = SFtmp.ImgNo() - refno * Nsub;
        for (int nn = 0; nn < Nsub; nn++)
        {
            fn_tmp=SFtmp.NextImg();
            Itmp.read(fn_tmp);
            Itmp().setXmippOrigin();
            random_rot = 360. * rnd_unif(0., 1.);
            random_tilt= ACOSD((2.*rnd_unif(0., 1.) - 1));
            random_psi = 360. * rnd_unif(0., 1.);
            random_A = Euler_rotation3DMatrix(random_rot, random_tilt, random_psi);
            Itmp().selfApplyGeometry( random_A, IS_NOT_INV, DONT_WRAP);

            // Going through the docfile again here is a bit dirty coding...
            // but let's first make it work...
            if (do_missing)
            {
                if (DF.search_comment(fn_tmp)) 
                {
                    int missno = (int)(DF(7)) - 1;
                    getMissingRegion(Mmissing, random_A, missno);
                }
                else 
                {
                    std::cerr << "ERROR% "<<fn_tmp
                              <<" not found in document file"
                              <<std::endl;
                    exit(0);
                }
            }
            if (nn == 0)
            {
                Iave() = Itmp();
                if (do_missing) 
                    Msumwedge = Mmissing;
            }
            else
            {
                Iave() += Itmp();
                if (do_missing) 
                    Msumwedge += Mmissing;
            }
            c++;
            if (verb > 0 && c % cc == 0) progress_bar(c);
        }

        if (do_missing)
        {
            // 1. Correct for missing wedge by division of sumwedge
            transformer.FourierTransform(Iave(),Fave,true);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fave)
            {
                if (DIRECT_MULTIDIM_ELEM(Msumwedge,n) > noimp_threshold)
                {
                    DIRECT_MULTIDIM_ELEM(Fave,n) /= 
                        DIRECT_MULTIDIM_ELEM(Msumwedge,n);
                }
                else
                {
                    DIRECT_MULTIDIM_ELEM(Fave,n) = 0.;
                }
            }
            transformer.inverseFourierTransform(Fave,Iave());
        }
        else
        {
            Iave() /= (double)Nsub;
        }

        // Enforce fourier_mask, symmetry and omask
        postProcessVolume(Iave);
        
        fn_tmp = fn_root + "_it";
        fn_tmp.compose(fn_tmp, 0, "");
        fn_tmp = fn_tmp + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "");
        fn_tmp = fn_tmp + ".vol";
        Iave.write(fn_tmp);
        SFr.insert(fn_tmp, SelLine::ACTIVE);
        
    }
    if (verb > 0) progress_bar(SF.ImgNo());
    fn_ref = fn_root + "_it";
    fn_ref.compose(fn_ref, 0, "sel");
    SFr.write(fn_ref);

#ifdef  DEBUG
    std::cerr<<"Finished generateInitialReferences"<<std::endl;
#endif
}

// Read reference images to memory and initialize offset vectors
// This side info is NOT general, i.e. in parallel mode it is NOT the
// same for all processors! (in contrast to produce_Side_info)
void Prog_ml_tomo_prm::produceSideInfo2(int nr_vols)
{

    int                       c, idum;
    DocFile                   DF, DFsub;
    DocLine                   DL;
    FileName                  fn_tmp;
    VolumeXmipp               img, Vaux;
    std::vector<Matrix1D<double> > Vdm;

#ifdef  DEBUG
    std::cerr<<"Start produceSideInfo2"<<std::endl;
    std::cerr<<"nr images= "<<SF.ImgNo()<<std::endl;
#endif

    // Read in all reference images in memory
    if (Is_VolumeXmipp(fn_ref))
    {
        SFr.reserve(1);
        SFr.insert(fn_ref);
    }
    else
    {
        SFr.read(fn_ref);
    }
    nr_ref = 0;
    SFr.go_beginning();
    while ((!SFr.eof()))
    {
        FileName fn_img=SFr.NextImg();
        img.read(fn_img);
        img().setXmippOrigin();

        // Enforce fourier_mask, symmetry and omask
        postProcessVolume(img);

        // Rotate some arbitrary (off-axis) angle and rotate back again to remove high frequencies
        // that will be affected by the interpolation due to rotation
        // This makes that the A2 values of the rotated references are much less sensitive to rotation
        img().selfApplyGeometry( Euler_rotation3DMatrix(32., 61., 53.), IS_NOT_INV, 
                                 DONT_WRAP, DIRECT_MULTIDIM_ELEM(img(),0) );
        img().selfApplyGeometry( Euler_rotation3DMatrix(32., 61., 53.), IS_INV, 
                                 DONT_WRAP, DIRECT_MULTIDIM_ELEM(img(),0) );

        Iref.push_back(img);
        Iold.push_back(img);
        nr_ref++;
    }

    // Store tomogram angles, offset vectors and missing wedge parameters
    imgs_missno.clear();
    imgs_optrefno.clear();
    imgs_optangno.clear();
    imgs_trymindiff.clear();
    imgs_scale.clear();
    imgs_bgmean.clear();

    for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
    {
        imgs_optrefno.push_back(0);
        imgs_optangno.push_back(0);
        imgs_trymindiff.push_back(-1.);
        if (do_missing)
            imgs_missno.push_back(-1);
        if (do_norm)
        {
            imgs_bgmean.push_back(0.);
            imgs_scale.push_back(1.);
        }
    }

    if (fn_doc!="") {
        DF.read(fn_doc);
        SF.go_beginning();
        int imgno = 0;
        DFsub.clear();
        while (!SF.eof()) 
        {
            fn_tmp=SF.NextImg();
            if (fn_tmp=="") break;
            DF.go_beginning();
            if (DF.search_comment(fn_tmp)) 
            {
                // Get missing wedge type
                if (do_missing)
                {
                    imgs_missno[imgno] = (int)(DF(7)) - 1;
                }
                // Get prior normalization parameters
                if (do_norm)
                {
                    imgs_bgmean[imgno] = DF(10);
                    imgs_scale[imgno] = DF(11);
                }
                // Generate a docfile from the (possible subset of) images in SF
                if (ang_search > 0.)
                {
                    DFsub.append_comment(fn_tmp);
                    DL=DF.get_current_line();
                    DFsub.append_line(DL);
                }
            } 
            else 
            {
                std::cerr << "ERROR% "<<fn_tmp
                          <<" not found in document file"
                          <<std::endl;
                exit(0);
            }
            imgno++;
        }

        // Set up local searches
        if (ang_search > 0.)
        {
            mysampling.SetNeighborhoodRadius(ang_search);
            mysampling.fill_exp_data_projection_direction_by_L_R(DFsub);
            mysampling.remove_points_far_away_from_experimental_data();
            mysampling.compute_neighbors();
        }

    } 

    // Now that we know where all experimental images are (for possible local searches)
    // we can finally calculate the sampling points for the references
    int nr_psi = CEIL(360. / psi_sampling);
    angle_info myinfo;
    all_angle_info.clear();
    nr_ang = 0;
    for (int i = 0; i < mysampling.no_redundant_sampling_points_angles.size(); i++)
    {
        double rot = XX(mysampling.no_redundant_sampling_points_angles[i]);
        double tilt = YY(mysampling.no_redundant_sampling_points_angles[i]);
        if (!no_SMALLANGLE)
        {
            rot  += SMALLANGLE;
            tilt += SMALLANGLE;
        }
        for (int ipsi = 0; ipsi < nr_psi; ipsi++)
        {
            double psi = (double)(ipsi * 360. / nr_psi);
            if (!no_SMALLANGLE)
            {
                psi  += SMALLANGLE;
            }
            myinfo.rot = rot;
            myinfo.tilt = tilt;
            myinfo.psi = psi;
            myinfo.direction = i;
            myinfo.A = Euler_rotation3DMatrix(rot, tilt, psi);
            all_angle_info.push_back(myinfo);
            nr_ang ++;
        }
    }


    // prepare priors
    alpha_k.resize(nr_ref);
    alpha_k_rot.resize(nr_ref*nr_ang);
    alpha_k.initConstant(1./(double)nr_ref);
    alpha_k_rot.initConstant(1./(double)(nr_ref*nr_ang));

    // Regularization
    regularize((double)SF.ImgNo());

//#define DEBUG_SAMPLING
#ifdef DEBUG_SAMPLING
    for (int angno = 0; angno < nr_ang; angno++)
    {
        double rot=all_angle_info[angno].rot;
        double tilt=all_angle_info[angno].tilt;
        double psi=all_angle_info[angno].psi;
        std::cerr<<"angno= "<<angno<<" rot= "<<rot<<" tilt= "<<tilt<<" psi= "<<psi<<std::endl;
    }
#endif

#define DEBUG_GENERAL
#ifdef DEBUG_GENERAL
    std::cerr<<"do_generate_refs ="<<do_generate_refs<<std::endl;
    std::cerr<<"nr_ref= "<<nr_ref<<std::endl; 
    std::cerr<<"nr_miss= "<<nr_miss<<std::endl;
    std::cerr<<"dim= "<<dim<<std::endl;
    std::cerr<<"Finished produceSideInfo"<<std::endl;
    std::cerr<<"nr_ang= "<<nr_ang<<std::endl; 
    std::cerr<<"nr_psi= "<<nr_psi<<std::endl; 
#endif

#ifdef DEBUG
    std::cerr<<"Finished produceSideInfo2"<<std::endl;
#endif

}

void Prog_ml_tomo_prm::getMissingRegion(Matrix3D<double> &Mmissing,
                                       Matrix2D<double> A,
                                       const int missno)
{

    if (missno < 0)
    {
        std::cerr<<" BUG: missno < 0"<<std::endl;
        exit(1);
    }

    Matrix2D<double> Ainv = A.inv();
    double xp, yp, zp;
    double tg0_y=0., tgF_y=0., tg0_x=0., tgF_x=0., tg=0., limx0=0., limxF=0., limy0=0., limyF=0., lim=0.;
    bool do_limit_x, do_limit_y, do_cone, is_observed;
    Mmissing.resize(dim,dim,hdim+1);

    if (all_missing_info[missno].type==MISSING_WEDGE_Y)
    {
        do_limit_x = true;
        do_limit_y = false;
        do_cone    = false;
         // TODO: Make this more logical!!
        tg0_y = -tan(PI * (-90. - all_missing_info[missno].thyF) / 180.);
        tgF_y = -tan(PI * (90. - all_missing_info[missno].thy0) / 180.);
    }
    else if (all_missing_info[missno].type==MISSING_WEDGE_X)
    {
        do_limit_x = false;
        do_limit_y = true;
        do_cone    = false;
        tg0_x = -tan(PI * (-90. - all_missing_info[missno].thxF) / 180.);
        tgF_x = -tan(PI * (90. - all_missing_info[missno].thx0) / 180.);
    }
    else if (all_missing_info[missno].type==MISSING_PYRAMID)
    {
        do_limit_x = true;
        do_limit_y = true;
        do_cone    = false;
        tg0_y = -tan(PI * (-90. - all_missing_info[missno].thyF) / 180.);
        tgF_y = -tan(PI * (90. - all_missing_info[missno].thy0) / 180.);
        tg0_x = -tan(PI * (-90. - all_missing_info[missno].thxF) / 180.);
        tgF_x = -tan(PI * (90. - all_missing_info[missno].thx0) / 180.);
    }
    else if (all_missing_info[missno].type==MISSING_CONE)
    {
        do_limit_x = false;
        do_limit_y = false;
        do_cone    = true;
        tg = -tan(PI * (-90. - all_missing_info[missno].thy0) / 180.);    }
    else
    {
        REPORT_ERROR(1,"bug: unrecognized type of missing region");
    }

//#define DEBUG_WEDGE
#ifdef DEBUG_WEDGE
    std::cerr<<"do_limit_x= "<<do_limit_x<<std::endl;
    std::cerr<<"do_limit_y= "<<do_limit_y<<std::endl;
    std::cerr<<"do_cone= "<<do_cone<<std::endl;
    std::cerr<<"tg0_y= "<<tg0_y<<std::endl;
    std::cerr<<"tgF_y= "<<tgF_y<<std::endl;
    std::cerr<<"tg0_x= "<<tg0_x<<std::endl;
    std::cerr<<"tgF_x= "<<tgF_x<<std::endl;
    std::cerr<<"XMIPP_EQUAL_ACCURACY= "<<XMIPP_EQUAL_ACCURACY<<std::endl;
    std::cerr<<"Ainv= "<<Ainv<<" A= "<<A<<std::endl;
#endif

    int zz, yy;
    int dimb = (dim - 1)/2;
    for ( int z=0, ii=0; z<dim; z++ ) {
        if ( z > dimb ) 
            zz = z-dim;
        else 
            zz = z;
        for ( int y=0; y<dim; y++ ) {
            if ( y > dimb ) yy = y-dim;
            else yy = y;
            for ( int xx=0; xx<hdim + 1; xx++, ii++ ) {

                double maskvalue= DIRECT_MULTIDIM_ELEM(fourier_mask,ii);
                if (maskvalue < XMIPP_EQUAL_ACCURACY)
                {
                    DIRECT_MULTIDIM_ELEM(Mmissing,ii) = 0.;
                }
                else
                {
                    // Rotate the wedge
                    xp = dMij(Ainv, 0, 0) * xx + dMij(Ainv, 0, 1) * yy + dMij(Ainv, 0, 2) * zz;
                    yp = dMij(Ainv, 1, 0) * xx + dMij(Ainv, 1, 1) * yy + dMij(Ainv, 1, 2) * zz;
                    zp = dMij(Ainv, 2, 0) * xx + dMij(Ainv, 2, 1) * yy + dMij(Ainv, 2, 2) * zz;

                    // Calculate the limits
                    if (do_cone)
                    {
                        lim = (tg * zp) * (tg * zp);
                        if (xp*xp + yp*yp >= lim)
                            is_observed = true;
                        else
                            is_observed = false;
                    }
                    else
                    {
                        is_observed = false; // for pyramid
                        if (do_limit_x)
                        {
                            limx0 = tg0_y * zp;
                            limxF = tgF_y * zp;
                            if (zp >= 0)
                            {
                                if (xp <= limx0 || xp >= limxF)
                                    is_observed = true;
                                else
                                    is_observed = false;
                            }
                            else
                            {
                                if (xp <= limxF || xp >= limx0)
                                    is_observed = true;
                                else
                                    is_observed = false;
                            }
                        }
                        if (do_limit_y && !is_observed)
                        {
                            limy0 = tg0_x * zp;
                            limyF = tgF_x * zp;
                            if (zp >= 0)
                            {
                                if (yp <= limy0 || yp >= limyF)
                                    is_observed = true;
                                else
                                    is_observed = false;
                            }
                            else
                            {
                                if (yp <= limyF || yp >= limy0)
                                    is_observed = true;
                                else
                                    is_observed = false;
                            }
                        }
                    }

                    if (is_observed)
                        DIRECT_MULTIDIM_ELEM(Mmissing,ii) = maskvalue;
                    else
                        DIRECT_MULTIDIM_ELEM(Mmissing,ii) = 0.;
                }
            }
        }
    }

#ifdef DEBUG_WEDGE
    VolumeXmipp test(dim,dim,dim), ori;
    ori()=Mmissing;
    ori.write("oriwedge.fft");
    test().initZeros();
    //test().setXmippOrigin();

    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(test())
        if (j<hdim+1)
            DIRECT_VOL_ELEM(test(),k,i,j)=
                DIRECT_VOL_ELEM(Mmissing,k,i,j);
        else
            DIRECT_VOL_ELEM(test(),k,i,j)=
                DIRECT_VOL_ELEM(Mmissing,
                                (dim-k)%dim,
                                (dim-i)%dim,
                                dim-j);

    test.write("wedge.ftt");
    CenterFFT(test(),true);
    test.write("Fwedge.vol");
    test().setXmippOrigin();
    Matrix3D<int> ress(dim,dim,dim);
    ress.setXmippOrigin();
    BinaryWedgeMask(test(),all_missing_info[missno].thy0, all_missing_info[missno].thyF, A);
    test.write("Mwedge.vol");
#endif

}

// Calculate probability density function of all in-plane transformations phi
void Prog_ml_tomo_prm::calculatePdfTranslations()
{

#ifdef  DEBUG
    std::cerr<<"start calculatePdfTranslations"<<std::endl;
#endif

    double r2, pdfpix;
    P_phi.resize(dim, dim, dim);
    P_phi.setXmippOrigin();
    Mr2.resize(dim, dim, dim);
    Mr2.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_MATRIX3D(P_phi)
    {
        r2 = (double)(j * j + i * i + k * k);
        if (sigma_offset > 0.)
        {
            pdfpix = exp(- r2 / (2. * sigma_offset * sigma_offset));
            pdfpix *= pow(2. * PI * sigma_offset * sigma_offset, -3./2.);
        }
        else
        {
            if (j == 0 && i == 0) pdfpix = 1.;
            else pdfpix = 0.;
        }
        VOL_ELEM(P_phi, k, i, j) = pdfpix;
        VOL_ELEM(Mr2, k, i, j) = (float)r2;
    }

//#define  DEBUG_PDF_SHIFT
#ifdef DEBUG_PDF_SHIFT
    std::cerr<<" Sum of translation pdfs (should be one) = "<<P_phi.sum()<<std::endl;
#endif 

#ifdef  DEBUG
    std::cerr<<"finished calculatePdfTranslations"<<std::endl;
#endif

}

void Prog_ml_tomo_prm::maskSphericalAverageOutside(Matrix3D<double> &Min)
{
    double outside_density = 0., sumdd = 0.;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(real_omask)
    {
        outside_density += DIRECT_MULTIDIM_ELEM(Min,n)*DIRECT_MULTIDIM_ELEM(real_omask,n);
        sumdd += DIRECT_MULTIDIM_ELEM(real_omask,n);
    }
    outside_density /= sumdd;

    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(real_omask)
    {
        DIRECT_MULTIDIM_ELEM(Min,n) *= DIRECT_MULTIDIM_ELEM(real_mask,n);
        DIRECT_MULTIDIM_ELEM(Min,n) += (outside_density)*DIRECT_MULTIDIM_ELEM(real_omask,n);
    }


}

void Prog_ml_tomo_prm::postProcessVolume(VolumeXmipp &Vin)
{

    Matrix3D<std::complex<double> > Faux;
#ifdef DEBUG
    std::cerr<<"start postProcessVolume"<<std::endl;
#endif

    // Fourier transform
    transformer.FourierTransform(Vin(),Faux,true);

    // Apply fourier mask
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
    {
        DIRECT_MULTIDIM_ELEM(Faux,n) *= DIRECT_MULTIDIM_ELEM(fourier_mask,n);
    }

    // Inverse Fourier Transform
    transformer.inverseFourierTransform(Faux,Vin());

    // Apply symmetry
    if (mysampling.SL.SymsNo() > 1)
    {
        // Note that for no-imputation this is not correct!
        // One would have to symmetrize the missing wedges and the sum of the images separately
        if (do_missing && !do_impute)
            REPORT_ERROR(1,"Symmetrization and dont_impute together are not implemented...");

        std::cerr<<"Symmetrizing... with SL.SymsNo()="<<mysampling.SL.SymsNo()<<std::endl;
        VolumeXmipp Vaux=Vin;
        Vaux().initZeros();
        symmetrize(mysampling.SL, Vin, Vaux);
        Vin = Vaux;
    }

    // Spherical mask with average outside density
    maskSphericalAverageOutside(Vin());

#ifdef DEBUG
    std::cerr<<"finished postProcessVolume"<<std::endl;
#endif

}


// Calculate FT of each reference and calculate A2 =============
void Prog_ml_tomo_prm::precalculateA2(std::vector< VolumeXmippT<double> > &Iref)
{

#ifdef DEBUG
    std::cerr<<"start precalculateA2"<<std::endl;
    TimeStamp t0; 
    time_config();
    annotate_time(&t0);
#endif

    double rot, tilt, psi, AA, stdAA, corr;
    Matrix2D<double>  A_rot_inv(4,4), I(4,4);
    Matrix3D<double> Maux(dim,dim,dim), Mmissing;
    Matrix3D<std::complex<double> > Faux, Faux2;

    A2.clear();
    corrA2.clear();
    I.initIdentity();
    Maux.setXmippOrigin();
    for (int refno = 0; refno < nr_ref; refno++)
    {
        // Calculate A2 for all different orientations
        for (int angno = 0; angno < nr_ang; angno++)
        {
            A_rot_inv = ((all_angle_info[angno]).A).inv();
            // use DONT_WRAP and put density of first element outside 
            // i.e. assume volume has been processed with omask
            applyGeometry(Maux, A_rot_inv, Iref[refno](), IS_NOT_INV, 
                          DONT_WRAP, DIRECT_MULTIDIM_ELEM(Iref[refno](),0));
//#define DEBUG_PRECALC_A2_ROTATE
#ifdef DEBUG_PRECALC_A2_ROTATE
            VolumeXmipp Vt;
            Vt()=Maux;
            Vt.write("rot.vol");
            Vt()=Iref[refno]();
            Vt.write("ref.vol");
            std::cerr<<" angles= "<<(all_angle_info[angno]).rot<<" "<<(all_angle_info[angno]).tilt<<" "<<(all_angle_info[angno]).psi<<std::endl;
            std::cerr<<" Written volume rot.vol and ref.vol, press any key to continue ..."<<std::endl; 
            char c;
            std::cin >> c;
#endif
            AA = Maux.sum2();
            if (angno==0) 
            {
                stdAA = AA;
            }
            if (AA > 0) corr = sqrt(stdAA / AA);
            corrA2.push_back(corr);
            Maux *= corr;
            if (do_missing)
            {
                transformer.FourierTransform(Maux,Faux,true);
                for (int missno = 0; missno < nr_miss; missno++)
                {
                    getMissingRegion(Mmissing,I,missno);
                    Faux2 = Faux;
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux2)
                    {
                        DIRECT_MULTIDIM_ELEM(Faux2,n) *= DIRECT_MULTIDIM_ELEM(Mmissing,n);
                    }
                    transformer.inverseFourierTransform(Faux2,Maux);
                    A2.push_back(Maux.sum2());
//#define DEBUG_PRECALC_A2
#ifdef DEBUG_PRECALC_A2
                    std::cerr<<"rot= "<<all_angle_info[angno].rot<<" tilt= "<<all_angle_info[angno].tilt<<" psi= "<<all_angle_info[angno].psi<<std::endl;
                    std::cerr<<"refno= "<<refno<<" angno= "<<angno<<" missno= "<<missno<<" A2= "<<Maux.sum2()<<" corrA2= "<<corr<<std::endl;
//#define DEBUG_PRECALC_A2_b
#ifdef DEBUG_PRECALC_A2_b
                    VolumeXmipp tt;
                    tt()=Maux;
                    tt.write("refrotwedge.vol");
                    std::cerr<<"press any key"<<std::endl;
                    char c;
                    std::cin >> c;
#endif
#endif
                }
            }
            else
            {
                A2.push_back(Maux.sum2());
#ifdef DEBUG_PRECALC_A2
                std::cerr<<"refno= "<<refno<<" angno= "<<angno<<" A2= "<<Maux.sum2()<<std::endl;
#endif
            }
        }
    }

#ifdef DEBUG
    std::cerr<<"finished precalculateA2"<<std::endl;
    print_elapsed_time(t0); 
#endif
}

// Maximum Likelihood calculation for one image ============================================
// Integration over all translation, given  model and in-plane rotation
void Prog_ml_tomo_prm::expectationSingleImage(
    Matrix3D<double> &Mimg, int imgno, int missno,
    std::vector< VolumeXmippT<double> > &Iref,
    std::vector<Matrix3D<double> > &wsumimgs,
    std::vector<Matrix3D<double> > &wsumweds,
    double &wsum_sigma_noise, double &wsum_sigma_offset,
    Matrix1D<double> &sumw, Matrix1D<double> &sumwsc, 
    Matrix1D<double> &sumwsc2, Matrix1D<double> &sumw_rot, 
    double &LL, double &dLL, double &fracweight, double &sumfracweight, 
    double &opt_scale, double &bgmean, double &trymindiff, 
    int &opt_refno, int &opt_angno, Matrix1D<double> &opt_offsets)
{

#ifdef DEBUG
    std::cerr<<"start expectationSingleImage"<<std::endl;
    TimeStamp t0, t1; 
    time_config();
    annotate_time(&t0);
    annotate_time(&t1);
#endif

    Matrix3D<double> Maux, Maux2, Mweight, Mmissing, Mzero(dim,dim,hdim+1), Mzero2(dim,dim,dim);
    Matrix3D<std::complex<double> > Faux, Faux2(dim,dim,hdim+1), Fimg, Fimg0, Fimg_rot;
    std::vector<Matrix3D<double> > mysumimgs;
    std::vector<Matrix3D<double> > mysumweds;
    Matrix1D<double> refw, refw_rot;
    double sigma_noise2, aux, pdf, fracpdf, myA2, mycorrAA, myXi2, A2_plus_Xi2;
    double mind, diff, mindiff, my_mindiff;
    double my_sumweight, weight, ref_scale = 1.;
    double wsum_sc, wsum_sc2, wsum_offset, old_bgmean;
    double wsum_corr, sum_refw, maxweight, my_maxweight;
    double rot, tilt, psi;
    int irot, irefmir, sigdim, xmax, ymax;
    int ioptpsi = 0, ioptlib = 0, ioptx = 0, iopty = 0, ioptz = 0, imax = 0;
    bool is_ok_trymindiff = false;
    int old_optrefno = opt_refno;
    int old_optangno = opt_angno;
    std::vector<double> all_Xi2;
    Matrix2D<double> A_rot(4,4), I(4,4), A_rot_inv(4,4);
    bool is_a_neighbor;
    XmippFftw local_transformer;

    // Only translations smaller than 6 sigma_offset are considered!
    // TODO: perhaps 3 sigma??
    I.initIdentity();
    sigdim = 2 * CEIL(sigma_offset * 6);
    sigdim++; // (to get uneven number)
    sigdim = XMIPP_MIN(dim, sigdim);
    // Setup matrices and constants
    Maux.resize(dim, dim, dim);
    Maux2.resize(dim, dim, dim);
    Maux.setXmippOrigin();
    Maux2.setXmippOrigin();
    Mweight.initZeros(sigdim, sigdim, sigdim);
    Mweight.setXmippOrigin();
    Mzero.initZeros();
    Mzero2.initZeros();
    Mzero2.setXmippOrigin();

    sigma_noise2 = sigma_noise * sigma_noise;
    if (!do_norm) 
        opt_scale =1.;

    // Calculate the unrotated Fourier transform with enforced wedge of Mimg (store in Fimg0)
    // also calculate myXi2;
    // Note that from here on local_transformer will act on Maux <-> Faux
    Maux=Mimg;
    local_transformer.FourierTransform(Maux, Faux, false);
    if (do_missing)
    {
        // Enforce missing wedge
        getMissingRegion(Mmissing,I,missno);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
        {
            DIRECT_MULTIDIM_ELEM(Faux,n) *= DIRECT_MULTIDIM_ELEM(Mmissing,n);
        }
    }
    Fimg0 = Faux;
    // BE CAREFUL: inverseFourierTransform messes up Faux!
    local_transformer.inverseFourierTransform();
    myXi2 = Maux.sum2();

    // To avoid numerical problems, subtract smallest difference from all differences.
    // That way: Pmax will be one and all other probabilities will be [0,1>
    // But to find mindiff I first have to loop over all hidden variables...
    // Fortunately, there is some flexibility as the reliable domain of the exp-functions goes from -700 to 700.
    // Therefore, make a first guess for mindiff (trymindiff) and check whether the difference 
    // with the real mindiff is not larger than 500 (to be on the save side).
    // If that is the case: OK; if not: do the entire loop again. 
    // The efficiency of this will depend on trymindiff. 
    // First cycle use trymindiff = trymindiff_factor * 0.5 * Xi2
    // From then on: use mindiff from the previous iteration
    if (trymindiff < 0.)
        // 90% of Xi2 may be a good idea (factor half because 0.5*diff is calculated)
        trymindiff = trymindiff_factor * 0.5 * myXi2;
    int redo_counter = 0;
    while (!is_ok_trymindiff)
    {
        // Initialize mindiff, weighted sums and maxweights
        mindiff = 99.e99;
        wsum_corr = wsum_offset = wsum_sc = wsum_sc2 = 0.;
        maxweight = sum_refw = 0.;
        mysumimgs.clear();
        mysumweds.clear();
        for (int refno = 0; refno < nr_ref; refno++)
        {
            mysumimgs.push_back(Mzero2); 
            if (do_missing)
                mysumweds.push_back(Mzero); 
        }
        refw.initZeros(nr_ref);
        refw_rot.initZeros(nr_ref*nr_ang);
        // The real stuff: now loop over all orientations, references and translations
        // Start the loop over all refno at old_optangno (=opt_angno from the previous iteration). 
        // This will speed-up things because we will find Pmax probably right away,
        // and this will make the if-statement that checks SIGNIFICANT_WEIGHT_LOW
        // effective right from the start
        for (int aa = old_optangno; aa < old_optangno+nr_ang; aa++)
        {
            int angno = aa;
            if (angno >= nr_ang) angno -= nr_ang;

            // See whether this image is in the neighborhoood for this imgno
            if (ang_search > 0.)
            {
                is_a_neighbor = false;
                for (int i = 0; i < mysampling.my_neighbors[imgno].size(); i++)
                {
                    if (mysampling.my_neighbors[imgno][i] == (all_angle_info[angno]).direction)
                    {
                        is_a_neighbor = true;
                        break;
                    }
                }
            }
            else
            {
                is_a_neighbor = true;
            }

            // If it is in the neighborhood: proceed
            if (is_a_neighbor)
            {

                A_rot = (all_angle_info[angno]).A;
                A_rot_inv = A_rot.inv();

                // Start the loop over all refno at old_optrefno (=opt_refno from the previous iteration). 
                // This will speed-up things because we will find Pmax probably right away,
                // and this will make the if-statement that checks SIGNIFICANT_WEIGHT_LOW
                // effective right from the start
                for (int rr = old_optrefno; rr < old_optrefno+nr_ref; rr++)
                {
                    int refno = rr;
                    if (refno >= nr_ref) refno-= nr_ref;
                    
                    fracpdf = alpha_k_rot(refno*nr_ang + angno);
                    // Now (inverse) rotate the reference and calculate its Fourier transform
                    // Use DONT_WRAP and assume map has been omasked
                    applyGeometry(Maux2, A_rot_inv, Iref[refno](), IS_NOT_INV, 
                                  WRAP, DIRECT_MULTIDIM_ELEM(Iref[refno](),0));
                    mycorrAA = corrA2[refno*nr_ang + angno];
                    Maux = Maux2 * mycorrAA;
                    local_transformer.FourierTransform();
                    if (do_norm) 
                        ref_scale = opt_scale / refs_avgscale[refno];
                    if (do_missing)
                        myA2 = A2[refno*nr_ang*nr_miss + angno*nr_miss + missno];
                    else
                        myA2 = A2[refno*nr_ang + angno];
                    A2_plus_Xi2 = 0.5 * ( ref_scale*ref_scale*myA2 + myXi2 );

//#define DEBUG_PRECALCULATE_A2
#ifdef DEBUG_PRECALCULATE_A2
                    std::cerr<<"A2_plus_Xi2 for refno= "<<refno<<" angno= "<<angno<<" = "<< A2_plus_Xi2<<" = 0.5 * ("<< ref_scale*ref_scale*myA2<<" + "<< myXi2<<")"<<std::endl;
#endif
                    
                    // A. Backward FFT to calculate weights in real-space
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
                    {
                        dVkij(Faux,k,i,j) = 
                            dVkij(Fimg0,k,i,j) * 
                            conj(dVkij(Faux,k,i,j));
                    }
                    local_transformer.inverseFourierTransform();
                    CenterFFT(Maux, true);
                    
//#define DEBUG_ALOT_SINGLEEXP
#ifdef DEBUG_ALOT_SINGLEEXP
                    std::cerr<<"rot= "<<all_angle_info[angno].rot<<" tilt= "<<all_angle_info[angno].tilt<<" psi= "<<all_angle_info[angno].psi<<" A2= "<<myA2<<" Xi2= "<<myXi2<<" corrA="<<mycorrAA <<" XA= "<<VOL_ELEM(Maux, 0,0,0) * ddim3<<" diff= "<< A2_plus_Xi2 - ref_scale * VOL_ELEM(Maux, 0, 0, 0) * ddim3<<" mindiff= "<<mindiff<<std::endl;
#endif

                    // B. Calculate weights for each pixel within sigdim (Mweight)
                    my_sumweight = my_maxweight = 0.;
                    FOR_ALL_ELEMENTS_IN_MATRIX3D(Mweight)
                    {
                        diff = A2_plus_Xi2 - ref_scale * VOL_ELEM(Maux, k, i, j) * ddim3;
                        mindiff = XMIPP_MIN(mindiff,diff);
#ifdef DEBUG_MINDIFF
                        if (mindiff < 0)
                        {
                            std::cerr<<"k= "<<k<<" j= "<<j<<" i= "<<i<<std::endl;
                            std::cerr<<"xaux="<<STARTINGX(Maux)<<" xw= "<<STARTINGX(Mweight)<<std::endl;
                            std::cerr<<"yaux="<<STARTINGY(Maux)<<" yw= "<<STARTINGY(Mweight)<<std::endl;
                            std::cerr<<"zaux="<<STARTINGZ(Maux)<<" zw= "<<STARTINGZ(Mweight)<<std::endl;
                            std::cerr<<"diff= "<<diff<<" A2_plus_Xi"<<std::endl;
                            std::cerr<<" mycorrAA= "<<mycorrAA<<" "<<std::endl;
                            std::cerr<<"debug mindiff= " <<mindiff<<" trymindiff= "<<trymindiff<< std::endl;
                            std::cerr<<"A2_plus_Xi2= "<<A2_plus_Xi2<<" myA2= "<<myA2<<" myXi2= "<<myXi2<<std::endl;
                            std::cerr<<"ref_scale= "<<ref_scale<<" volMaux="<<VOL_ELEM(Maux, k, i, j)<<std::endl;
                            std::cerr.flush();
                            VolumeXmipp tt;
                            tt()=Maux; tt.write("Maux.vol");
                            tt()=Mweight; tt.write("Mweight.vol");
                            std::cerr<<"Ainv= "<<A_rot_inv<<std::endl;
                            std::cerr<<"A= "<<A_rot_inv.inv()<<std::endl;                            
                            exit(1);
                        }
#endif
                        pdf = fracpdf * VOL_ELEM(P_phi, k, i, j);
                        // Normal distribution
                        aux = (diff - trymindiff) / sigma_noise2;
                        // next line because of numerical precision of exp-function
                        if (aux > 1000.) weight = 0.;
                        else weight = exp(-aux) * pdf;
                        VOL_ELEM(Mweight, k, i, j) = weight;
                        // Accumulate sum weights for this (my) matrix
                        my_sumweight += weight;
                        // calculate weighted sum of (X-A)^2 for sigma_noise update
                        wsum_corr += weight * diff;
                        // calculated weighted sum of offsets as well
                        wsum_offset += weight * VOL_ELEM(Mr2, k, i, j);
                        if (do_norm)
                        {
                            // weighted sum of Sum_j ( X_ij*A_kj )
                            wsum_sc += weight * (A2_plus_Xi2 - diff) / ref_scale;
                            // weighted sum of Sum_j ( A_kj*A_kj )
                            wsum_sc2 += weight * myA2;
                        }				
                        // keep track of optimal parameters
                        my_maxweight = XMIPP_MAX(my_maxweight, weight);
                        if (weight > maxweight)
                        {
                            maxweight = weight;
                            ioptz = k;
                            iopty = i;
                            ioptx = j;
                            opt_angno = angno;
                            opt_refno = refno;
                        }
                        
                    } // close for over all elements in Mweight

                    // C. only for significant settings, store weighted sums
                    if (my_maxweight > SIGNIFICANT_WEIGHT_LOW*maxweight )
                    {
                        sum_refw += my_sumweight;
                        refw(refno) += my_sumweight;
                        refw_rot(refno*nr_ang + angno) += my_sumweight;
                        // Back from smaller Mweight to original size of Maux
                        Maux.initZeros();
                        FOR_ALL_ELEMENTS_IN_MATRIX3D(Mweight)
                        {
                            VOL_ELEM(Maux, k, i, j) = VOL_ELEM(Mweight, k, i, j);
                        }
                        // Use forward FFT in convolution theorem again
                        CenterFFT(Maux, false);
                        local_transformer.FourierTransform();
                        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
                        {
                            dVkij(Faux, k, i, j) = conj(dVkij(Faux,k,i,j)) * dVkij(Fimg0,k,i,j);
                        }
                        local_transformer.inverseFourierTransform();
                        maskSphericalAverageOutside(Maux);
                        Maux.selfApplyGeometry(A_rot, IS_NOT_INV, 
                                               DONT_WRAP, DIRECT_MULTIDIM_ELEM(Maux,0));
                        mysumimgs[refno] += Maux * ddim3;
                        if (do_missing)
                        {
                            // Store sum of wedges!
                            getMissingRegion(Mmissing, A_rot, missno);
                            mysumweds[refno] += my_sumweight * Mmissing;
                        }
                    } // close if SIGNIFICANT_WEIGHT_LOW
                } // close for refno
            } // close if is_a_neighbor
        } // close for angno
            
        // Now check whether our trymindiff was OK.
        // The limit of the exp-function lies around 
        // exp(700)=1.01423e+304, exp(800)=inf; exp(-700) = 9.85968e-305; exp(-800) = 0
        // Use 500 to be on the save side?
        if (ABS((mindiff - trymindiff) / sigma_noise2) > 500.)
        {

//#define DEBUG_REDOCOUNTER
#ifdef DEBUG_REDOCOUNTER
            std::cerr<<"repeating mindiff "<<redo_counter<<"th time"<<std::endl;
            std::cerr<<"trymindiff= "<<trymindiff<<" mindiff= "<<mindiff<<std::endl;
            std::cerr<<"diff= "<<ABS((mindiff - trymindiff) / sigma_noise2)<<std::endl;
#endif
            // Re-do whole calculation now with the real mindiff
            trymindiff = mindiff;
            redo_counter++;
            // Never re-do more than once!
            if (redo_counter>1)
            {
                std::cerr<<"ml_tomo BUG% redo_counter > 1"<<std::endl;
                exit(1);
            }
        }
        else
        {
            is_ok_trymindiff = true;
            my_mindiff = trymindiff;
            trymindiff = mindiff;
        }
    }

    fracweight = maxweight / sum_refw;
    wsum_sc /= sum_refw;
    wsum_sc2 /= sum_refw;

    // Calculate remaining optimal parameters
    XX(opt_offsets) = -(double)ioptx;
    YY(opt_offsets) = -(double)iopty;
    ZZ(opt_offsets) = -(double)ioptz;
    if (do_norm)
    {
        // Note that bgmean will always be zero 
        // because we omitted the origin pixel in fourier_mas
        opt_scale = wsum_sc / wsum_sc2;
    }

    // From here on lock threads
    pthread_mutex_lock( &mltomo_weightedsum_update_mutex );

    // Update all global weighted sums after division by sum_refw
    wsum_sigma_noise += (2 * wsum_corr / sum_refw);
    wsum_sigma_offset += (wsum_offset / sum_refw);
    sumfracweight += fracweight;
    for (int refno = 0; refno < nr_ref; refno++)
    {
        sumw(refno) += refw(refno) / sum_refw;
        sumwsc(refno) += (refw(refno) * opt_scale) / sum_refw;
        sumwsc2(refno) += (refw(refno) * opt_scale * opt_scale) / sum_refw;
        // Sum mysumimgs to the global weighted sum
        wsumimgs[refno] += (opt_scale * mysumimgs[refno]) / sum_refw;
        if (do_missing)
        {
            wsumweds[refno] += mysumweds[refno] / sum_refw;
        }
        for (int angno = 0; angno < nr_ang; angno++)
            sumw_rot(refno*nr_ang + angno) += refw_rot(refno*nr_ang + angno) / sum_refw;
    }

    // 1st term: log(refw_i)
    // 2nd term: for subtracting mindiff
    // 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
    // TODO: check this!!
    dLL = log(sum_refw) 
        - my_mindiff / sigma_noise2 
        - ddim3 * log( sqrt(2. * PI * sigma_noise2));
    LL += dLL;

    pthread_mutex_unlock(  &mltomo_weightedsum_update_mutex );

#ifdef DEBUG
    std::cerr<<"finished expectationSingleImage"<<std::endl;
    print_elapsed_time(t0); 
#endif
}


void Prog_ml_tomo_prm::maxConstrainedCorrSingleImage(
    Matrix3D<double> &Mimg, int imgno, int missno,
    std::vector<VolumeXmippT<double> > &Iref,
    std::vector<Matrix3D<double> > &wsumimgs,
    std::vector<Matrix3D<double> > &wsumweds,
    Matrix1D<double> &sumw, double &maxCC, double &sumCC,
    int &opt_refno, int &opt_angno, Matrix1D<double> &opt_offsets)
{
#ifdef DEBUG
    std::cerr<<"start maxConstrainedCorrSingleImage"<<std::endl;
    TimeStamp t0, t1; 
    time_config();
    annotate_time(&t0);
    annotate_time(&t1);
#endif

    Matrix3D<double> Mimg0, Maux, Mref, Mmissing;
    Matrix3D<std::complex<double> > Faux, Fimg0, Fref;
    XmippFftw local_transformer;
    Matrix2D<double> A_rot(4,4), I(4,4), A_rot_inv(4,4);
    bool is_a_neighbor;
    double img_stddev, ref_stddev, corr, maxcorr=-9999.;
    int ioptx,iopty,ioptz;

    I.initIdentity();
    Maux.resize(dim, dim, dim);
    Maux.setXmippOrigin();

    // Calculate the unrotated Fourier transform with enforced wedge of Mimg (store in Fimg0)
    // also calculate myXi2;
    // Note that from here on local_transformer will act on Maux <-> Faux
    Maux=Mimg;
    local_transformer.FourierTransform(Maux, Faux, false);
    if (do_missing)
    {
        // Enforce missing wedge
        getMissingRegion(Mmissing,I,missno);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
        {
            DIRECT_MULTIDIM_ELEM(Faux,n) *= DIRECT_MULTIDIM_ELEM(Mmissing,n);
        }
        // BE CAREFUL: inverseFourierTransform messes up Faux!!
        Fimg0 = Faux;        
        local_transformer.inverseFourierTransform();
    }
    else
        Fimg0 = Faux;
    Mimg0 = Maux;
    
    // Calculate stddev of (wedge-inforced) image
    img_stddev = Mimg0.computeStddev();

    // Loop over all orientations
    for (int angno = 0; angno < nr_ang; angno++)
    {
        // See whether this image is in the neighborhoood for this imgno
        if (ang_search > 0.)
        {
            is_a_neighbor = false;
            for (int i = 0; i < mysampling.my_neighbors[imgno].size(); i++)
            {
                if (mysampling.my_neighbors[imgno][i] == (all_angle_info[angno]).direction)
                {
                    is_a_neighbor = true;
                    break;
                }
            }
        }
        else
        {
            is_a_neighbor = true;
        }
        
        // If it is in the neighborhoood: proceed
        if (is_a_neighbor)
        {
            A_rot = (all_angle_info[angno]).A;
            A_rot_inv = A_rot.inv();
            // Loop over all references
            for (int refno = 0; refno < nr_ref; refno++)
            {
                // Now (inverse) rotate the reference and calculate its Fourier transform
                // Use DONT_WRAP because the reference has been omasked
                applyGeometry(Maux, A_rot_inv, Iref[refno](), IS_NOT_INV, 
                              DONT_WRAP, DIRECT_MULTIDIM_ELEM(Iref[refno](),0));
                local_transformer.FourierTransform();
                if (do_missing)
                {
                    // Enforce wedge on the reference
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                    {
                        DIRECT_MULTIDIM_ELEM(Faux,n) *= DIRECT_MULTIDIM_ELEM(Mmissing,n);
                    }
                    // BE CAREFUL! inverseFourierTransform messes up Faux
                    Fref = Faux;
                    local_transformer.inverseFourierTransform();
                }
                else
                    Fref = Faux;
                // Calculate stddev of (wedge-inforced) reference
                Mref = Maux;
                ref_stddev = Mref.computeStddev();

                // Calculate correlation matrix via backward FFT
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
                {
                    dVkij(Faux,k,i,j) = 
                        dVkij(Fimg0,k,i,j) * 
                        conj(dVkij(Fref,k,i,j));
                }
                local_transformer.inverseFourierTransform();
                CenterFFT(Maux, true);

                FOR_ALL_ELEMENTS_IN_MATRIX3D(Maux)
                {
                    corr = VOL_ELEM(Maux,k,i,j) / (img_stddev * ref_stddev);
                    if (corr > maxcorr)
                    {
                        maxcorr = corr;
                        ioptz = k;
                        iopty = i;
                        ioptx = j;
                        opt_angno = angno;
                        opt_refno = refno;
                    }
                }
            }
        }
    }

    // Now that we have optimal hidden parameters, update output
    XX(opt_offsets) = -(double)ioptx;
    YY(opt_offsets) = -(double)iopty;
    ZZ(opt_offsets) = -(double)ioptz;
    Mimg0.selfTranslate(opt_offsets, DONT_WRAP);
    A_rot = (all_angle_info[opt_angno]).A;

    maskSphericalAverageOutside(Mimg0);
    Mimg0.selfApplyGeometry(A_rot, IS_NOT_INV, 
                            DONT_WRAP, DIRECT_MULTIDIM_ELEM(Mimg0,0));
    maxCC = maxcorr;

    // From here on lock threads
    pthread_mutex_lock( &mltomo_weightedsum_update_mutex );

    sumCC += maxCC;
    wsumimgs[opt_refno] += Mimg0;
    sumw(opt_refno)  += 1.;
    if (do_missing)
    {
        // Store sum of wedges
        getMissingRegion(Mmissing, A_rot, missno);
        wsumweds[opt_refno] += Mmissing;
    }

    pthread_mutex_unlock(  &mltomo_weightedsum_update_mutex );

#ifdef DEBUG
    std::cerr<<"finished expectationSingleImage"<<std::endl;
    print_elapsed_time(t0); 
#endif

}


void * threadMLTomoExpectationSingleImage( void * data )
{
    structThreadExpectationSingleImage * thread_data = (structThreadExpectationSingleImage *) data;

    // Variables from above
    int thread_id = thread_data->thread_id;
    int thread_num = thread_data->thread_num;
    Prog_ml_tomo_prm *prm = thread_data->prm;
    SelFile *SF = thread_data->SF;
    int *iter = thread_data->iter;
    double *wsum_sigma_noise = thread_data->wsum_sigma_noise;
    double *wsum_sigma_offset = thread_data->wsum_sigma_offset;
    double *sumfracweight = thread_data->sumfracweight;
    double *LL = thread_data->LL;
    std::vector<Matrix3D<double > > *wsumimgs = thread_data->wsumimgs;
    std::vector<Matrix3D<double > > *wsumweds = thread_data->wsumweds;
    std::vector< VolumeXmippT<double> > *Iref = thread_data->Iref;
    std::vector<Matrix1D<double > > *docfiledata = thread_data->docfiledata;
    Matrix1D<double> *sumw = thread_data->sumw;
    Matrix1D<double> *sumwsc = thread_data->sumwsc;
    Matrix1D<double> *sumwsc2 = thread_data->sumwsc2;
    Matrix1D<double> *sumw_rot = thread_data->sumw_rot;

//#define DEBUG_THREAD
#ifdef DEBUG_THREAD
    std::cerr<<"start threadMLTomoExpectationSingleImage"<<std::endl;
#endif

    // Local variables
    VolumeXmipp img;
    FileName fn_img, fn_trans;
    std::vector<Matrix1D<double> > allref_offsets;
    Matrix1D<double> opt_offsets(3);
    float old_phi = -999., old_theta = -999.;
    double fracweight, maxweight2, trymindiff, dLL;
    double opt_scale = 1., bgmean;
    int opt_refno, opt_angno, missno;

    // Calculate myFirst and myLast image for this thread
    int nn = (*SF).ImgNo();
    int remaining = nn % thread_num;
    int Npart = (int)(nn - remaining) / thread_num;
    int myFirst, myLast, myNum;
    if (thread_id < remaining)
    {
        myFirst = thread_id * (Npart + 1);
        myLast = myFirst + Npart;
    }
    else
    {
        myFirst = thread_id * Npart + remaining;
        myLast = myFirst + Npart - 1;
    }
    myNum = myLast - myFirst + 1;
            
    if (prm->verb > 0 && thread_id==0) init_progress_bar(myNum);

    // Loop over all images
    for (int imgno = myFirst; imgno <= myLast; imgno++)
    {
        pthread_mutex_lock(  &mltomo_selfile_access_mutex );
        (*SF).go_beginning();
        (*SF).jump(imgno, SelLine::ACTIVE);
        fn_img = (*SF).get_current_file();
        pthread_mutex_unlock(  &mltomo_selfile_access_mutex );

        img.read(fn_img);
        img().setXmippOrigin();
        
        if (prm->do_missing)
            missno = prm->imgs_missno[imgno];
        else
            missno = -1;

        if (prm->do_ml)
        {
            // A. Use maximum likelihood approach

            // These three parameters speed up expectationSingleImage
            trymindiff = prm->imgs_trymindiff[imgno];
            opt_refno = prm->imgs_optrefno[imgno];
            opt_angno = prm->imgs_optangno[imgno];
            if (prm->do_norm)
            {
                bgmean=prm->imgs_bgmean[imgno];
                opt_scale = prm->imgs_scale[imgno];
            }
            
            (*prm).expectationSingleImage(img(), imgno, missno, *Iref, *wsumimgs, *wsumweds, 
                                          *wsum_sigma_noise, *wsum_sigma_offset, 
                                          *sumw, *sumwsc, *sumwsc2, *sumw_rot, 
                                          *LL, dLL, fracweight, *sumfracweight, 
                                          opt_scale, bgmean, trymindiff, 
                                          opt_refno, opt_angno, opt_offsets);

            // Store for next iteration
            prm->imgs_trymindiff[imgno] = trymindiff;
            prm->imgs_optrefno[imgno] = opt_refno;
            prm->imgs_optangno[imgno] = opt_angno;
            if (prm->do_norm)
            {
                prm->imgs_scale[imgno] = opt_scale;
                prm->imgs_bgmean[imgno]  = bgmean;
            }
        }
        else
        {
            // B. Use constrained correlation coefficient approach
            (*prm).maxConstrainedCorrSingleImage(img(), imgno, missno, *Iref, 
                                                 *wsumimgs, *wsumweds, *sumw, fracweight, *sumfracweight, 
                                                 opt_refno, opt_angno, opt_offsets);
        }

        
        // Output docfile
        (*docfiledata)[imgno](0) = (prm->all_angle_info[opt_angno]).rot;// rot
        (*docfiledata)[imgno](1) = (prm->all_angle_info[opt_angno]).tilt;// tilt
        (*docfiledata)[imgno](2) = (prm->all_angle_info[opt_angno]).psi;// psi
        (*docfiledata)[imgno](3) = opt_offsets(0);            // Xoff
        (*docfiledata)[imgno](4) = opt_offsets(1);            // Yoff
        (*docfiledata)[imgno](5) = opt_offsets(2);            // Zoff
        (*docfiledata)[imgno](6) = (double)(opt_refno + 1);   // Ref
        (*docfiledata)[imgno](7) = (double)(missno + 1);      // Wedge number

        (*docfiledata)[imgno](8) = fracweight;                // P_max/P_tot
        (*docfiledata)[imgno](9) = dLL;                       // log-likelihood
        if (prm->do_norm)
        {
            (*docfiledata)[imgno](10)  = bgmean;               // background mean
            (*docfiledata)[imgno](11) = opt_scale;             // image scale 
        }
        else
        {
            (*docfiledata)[imgno](10)  = 0.;                   // background mean
            (*docfiledata)[imgno](11) = 1.;                    // image scale 

        }
        if (prm->verb > 0 && thread_id==0) progress_bar(imgno+1);

    }
    if (prm->verb > 0 && thread_id==0) progress_bar(myNum);

#ifdef DEBUG_THREAD
    std::cerr<<"finished threadMLTomoExpectationSingleImage"<<std::endl;
#endif

}

void Prog_ml_tomo_prm::expectation(
        SelFile &SF, std::vector< VolumeXmippT<double> > &Iref, int iter,
        double &LL, double &sumfracweight, DocFile &DFo,
        std::vector<Matrix3D<double> > &wsumimgs,
        std::vector<Matrix3D<double> > &wsumweds,
        double &wsum_sigma_noise, double &wsum_sigma_offset, 
	Matrix1D<double> &sumw, Matrix1D<double> &sumwsc, 
        Matrix1D<double> &sumwsc2, Matrix1D<double> &sumw_rot)
{

#ifdef DEBUG
    std::cerr<<"start expectation"<<std::endl;
#endif 

    Matrix1D<double> dataline(MLTOMODATALINELENGTH);
    Matrix3D<double> Mzero(dim,dim,hdim+1), Mzero2(dim,dim,dim);
    std::vector<Matrix1D<double> > docfiledata;
    bool fill_real_space;
    int num_img_tot;
    
    if (do_ml)
    {
        // Precalculate values outside real_mask and A2-values for all references
        precalculateA2(Iref);
        // Pre-calculate pdf of all in-plane transformations
        calculatePdfTranslations();
    }

    // Initialize weighted sums
    LL = 0.;
    wsum_sigma_noise = 0.;
    wsum_sigma_offset = 0.;
    sumfracweight = 0.;
    sumw.initZeros(nr_ref);
    sumwsc.initZeros(nr_ref);
    sumwsc2.initZeros(nr_ref);
    sumw_rot.initZeros(nr_ref*nr_ang);
    dataline.initZeros();
    for (int i = 0; i < SF.ImgNo(); i++)
        docfiledata.push_back(dataline);
    Mzero.initZeros();
    Mzero2.initZeros();
    Mzero2.setXmippOrigin();
    wsumimgs.clear();
    wsumweds.clear();
    for (int refno = 0; refno < nr_ref; refno++)
    {
        wsumimgs.push_back(Mzero2);
        if (do_missing)
            wsumweds.push_back(Mzero);
    }


    // Call threads to calculate the expectation of each image in the selfile
    pthread_t * th_ids = (pthread_t *)malloc( threads * sizeof( pthread_t));
    structThreadExpectationSingleImage * threads_d = (structThreadExpectationSingleImage *) malloc ( threads * sizeof( structThreadExpectationSingleImage ) );
    for( int c = 0 ; c < threads ; c++ )
    {
        threads_d[c].thread_id = c;
        threads_d[c].thread_num = threads;
        threads_d[c].prm = this;
        threads_d[c].SF=&SF;
        threads_d[c].iter=&iter;
        threads_d[c].wsum_sigma_noise=&wsum_sigma_noise;
        threads_d[c].wsum_sigma_offset=&wsum_sigma_offset;
        threads_d[c].sumfracweight=&sumfracweight;
        threads_d[c].LL=&LL;
        threads_d[c].wsumimgs=&wsumimgs;
        threads_d[c].wsumweds=&wsumweds;
        threads_d[c].Iref=&Iref;
        threads_d[c].docfiledata=&docfiledata;
        threads_d[c].sumw=&sumw;
        threads_d[c].sumwsc=&sumwsc;
        threads_d[c].sumwsc2=&sumwsc2;
        threads_d[c].sumw_rot=&sumw_rot;
        pthread_create( (th_ids+c), NULL, threadMLTomoExpectationSingleImage, (void *)(threads_d+c) );
    }

    // Wait for threads to finish and get joined DocFile
    for( int c = 0 ; c < threads ; c++ )
    {
        pthread_join(*(th_ids+c),NULL);
    }

    // Send back output in the form of a DocFile
    SF.go_beginning();
    for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
    {
        DFo.append_comment(SF.NextImg());
        DFo.append_data_line(docfiledata[imgno]);
    }

#ifdef DEBUG
    std::cerr<<"finished expectation"<<std::endl;
#endif 

}

// Update all model parameters
void Prog_ml_tomo_prm::maximization(std::vector<Matrix3D<double> > &wsumimgs,
                                    std::vector<Matrix3D<double> > &wsumweds,
                                    double &wsum_sigma_noise, double &wsum_sigma_offset,
                                    Matrix1D<double> &sumw, Matrix1D<double> &sumwsc, 
                                    Matrix1D<double> &sumwsc2, Matrix1D<double> &sumw_rot,
                                    double &sumfracweight, double &sumw_allrefs)
{
#ifdef DEBUG
    std::cerr<<"started maximization"<<std::endl;
#endif 

    Matrix1D<double> rmean_sigma2, rmean_signal2;
    Matrix1D<int> center(3), radial_count;
    Matrix3D<std::complex<double> > Faux(dim,dim,hdim+1), Fwsumimgs;
    Matrix3D<double> Maux, Msumallwedges(dim,dim,hdim+1);
    VolumeXmipp Vaux;
    FileName fn_tmp;
    double rr, thresh, aux, sumw_allrefs2 = 0., avg_origin=0.;
    int c;

    // Update the reference images
    sumw_allrefs = 0.;
    Msumallwedges.initZeros();
    for (int refno=0;refno<nr_ref; refno++)
    {
        if (!do_norm)
        {
            // This is needed for maxCC
            sumwsc2(refno) = sumw(refno);
        }
        sumw_allrefs += sumw(refno);
        transformer.FourierTransform(wsumimgs[refno],Fwsumimgs,true);
        if (do_missing)
        {
            Msumallwedges += wsumweds[refno];
            if (do_impute)
            {
                transformer.FourierTransform(Iref[refno](),Faux,true);
                if (sumw(refno) > 0.)
                {
//#define DEBUG_IMPUTE
#ifdef DEBUG_IMPUTE
                    VolumeXmipp tt;
                    FileName fnt;
                    FFT_magnitude(Faux,tt());
                    fnt.compose("Fref0",refno,"ampl");
                    tt.write(fnt);
                    tt()=wsumweds[refno];
                    fnt.compose("wsum0",refno,"ampl");
                    tt.write(fnt);
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                    {
                        // Impute old reference for missing pixels
                        DIRECT_MULTIDIM_ELEM(tt(),n) = 
                            (1. - DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) / sumw(refno));
                        DIRECT_MULTIDIM_ELEM(Faux,n) *= 
                            (1. - DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) / sumw(refno));
                    }
                    fnt.compose("inv_wsum",refno,"ampl");
                    tt.write(fnt);
                    FFT_magnitude(Faux,tt());
                    fnt.compose("Fref1",refno,"ampl");
                    tt.write(fnt);
                    FFT_magnitude(Fwsumimgs,tt());
                    fnt.compose("Fwsumimgs",refno,"ampl");
                    tt.write(fnt);
                    std::cerr<<" sumw(refno)= "<<sumw(refno)<<" sumwsc2(refno)= "<<sumwsc2(refno)<<std::endl;
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                    {
                        // And sum the weighted sum for observed pixels
                        DIRECT_MULTIDIM_ELEM(Faux,n) += 
                            (DIRECT_MULTIDIM_ELEM(Fwsumimgs,n) / sumwsc2(refno));
                    }
                    FFT_magnitude(Faux,tt());
                    fnt.compose("Fref2",refno,"ampl");
                    tt.write(fnt);
#else
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                    {
                        // Impute old reference for missing pixels
                        DIRECT_MULTIDIM_ELEM(Faux,n) *= 
                            (1. - DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) / sumw(refno));
                        // And sum the weighted sum for observed pixels
                        DIRECT_MULTIDIM_ELEM(Faux,n) += 
                            (DIRECT_MULTIDIM_ELEM(Fwsumimgs,n) / sumwsc2(refno));
                    }
#endif
                }
                // else do nothing (i.e. impute old reference completely
            }
            else 
            {
                // no imputation: divide by number of times a pixel has been observed 
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                {
                    if (DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) > noimp_threshold)
                    {
                        DIRECT_MULTIDIM_ELEM(Faux,n) = 
                            DIRECT_MULTIDIM_ELEM(Fwsumimgs,n) / 
                            DIRECT_MULTIDIM_ELEM(wsumweds[refno],n);
                        // TODO:  CHECK THE FOLLOWING LINE!!
                        DIRECT_MULTIDIM_ELEM(Faux,n) *= sumw(refno) / sumwsc2(refno);
                    }
                    else
                    {
                        DIRECT_MULTIDIM_ELEM(Faux,n) = 0.;
                    }
                }
            }
        }
        else
        {
            // No missing data
            if (sumw(refno) > 0.)
            {   
                // if no missing data: calculate straightforward average, 
                // i.e. just divide wsumimgs[refno] by sumwsc2(refno)
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                {
                    DIRECT_MULTIDIM_ELEM(Faux,n) = 
                        DIRECT_MULTIDIM_ELEM(Fwsumimgs,n) / sumwsc2(refno);
                }
            }
        }
        // Back to real space
        transformer.inverseFourierTransform(Faux,Iref[refno]());
    }  

    // Adjust average scale
    if (do_norm) {
        average_scale = 0.;
        for (int refno=0;refno<nr_ref; refno++)
        {
            average_scale += sumwsc(refno);
        }
        for (int refno=0;refno<nr_ref; refno++)
        {
            if (sumw(refno)>0.)
            {
                refs_avgscale[refno] =  sumwsc(refno)/sumw(refno);
                Iref[refno]() *= refs_avgscale[refno];
            }
            else
            {
                refs_avgscale[refno] = 1.;
            }
        }
        average_scale /= sumw_allrefs;
    }

    // Average fracweight
    sumfracweight /= sumw_allrefs;
    
    // Update the model fractions
    if (!fix_fractions)
    {
        for (int refno=0;refno<nr_ref; refno++)
        {
            if (sumw(refno) > 0.)
                alpha_k(refno) = sumw(refno) / sumw_allrefs;
            else
                alpha_k(refno) = 0.;
            for (int angno = 0; angno < nr_ang; angno++)
            {
                if (sumw_rot(refno*nr_ang + angno) > 0.)
                    alpha_k_rot(refno*nr_ang + angno) = sumw_rot(refno*nr_ang + angno) / sumw_allrefs;
                else
                    alpha_k_rot(refno*nr_ang + angno) = 0.;
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
        if (do_missing)
        {
            double sum_complete_wedge = 0.;
            double sum_complete_fourier = 0.;
            Matrix3D<double> Mcomplete(dim,dim,dim);
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Mcomplete)
            {
                if (j<XSIZE(Msumallwedges))
                {
                    sum_complete_wedge += DIRECT_VOL_ELEM(Msumallwedges,k,i,j);
                    sum_complete_fourier += DIRECT_VOL_ELEM(fourier_mask,k,i,j);
                }
                else
                {
                    sum_complete_wedge += DIRECT_VOL_ELEM(Msumallwedges,
                                                          (dim-k)%dim,(dim-i)%dim,dim-j);
                    sum_complete_fourier += DIRECT_VOL_ELEM(fourier_mask,
                                                            (dim-k)%dim,(dim-i)%dim,dim-j);
                }
            }
//#define DEBUG_UPDATE_SIGMA
#ifdef DEBUG_UPDATE_SIGMA
            std::cerr<<" sum_complete_wedge= "<<sum_complete_wedge<<" = "<<100*sum_complete_wedge/(sumw_allrefs*sum_complete_fourier)<<"%"<<std::endl;
            std::cerr<<" sum_complete_fourier= "<<sum_complete_fourier<<std::endl;
            std::cerr<<" sumw_allrefs= "<<sumw_allrefs<<std::endl;
            std::cerr<<" wsum_sigma_noise= "<<wsum_sigma_noise<<std::endl;
            std::cerr<<" sigma_noise_old= "<<sigma_noise<<std::endl;
            std::cerr<<" sigma_new_impute= "<<sqrt((sigma_noise*(sumw_allrefs*sum_complete_fourier -  sum_complete_wedge  ) + wsum_sigma_noise)/(sumw_allrefs * sum_complete_fourier))<<std::endl;
            std::cerr<<" sigma_new_noimpute= "<<sqrt(wsum_sigma_noise / (sum_complete_wedge))<<std::endl;
            std::cerr<<" sigma_new_nomissing= "<<sqrt(wsum_sigma_noise / (sumw_allrefs * sum_complete_fourier))<<std::endl;
#endif
            if (do_impute)
            {
                sigma_noise *= sigma_noise;
                sigma_noise *= sumw_allrefs*sum_complete_fourier - sum_complete_wedge ;
                sigma_noise += wsum_sigma_noise;
                sigma_noise =  sqrt(sigma_noise / (sumw_allrefs * sum_complete_fourier));
            }
            else
            {
                sigma_noise = sqrt(wsum_sigma_noise / (sum_complete_wedge));
            }
        }
        else
        {
            sigma_noise = sqrt(wsum_sigma_noise / (sumw_allrefs * ddim3));
        }
    }

    regularize(sumw_allrefs);
    // Update regularization constant
    reg_current -= reg_step;
    reg_current = XMIPP_MAX(reg_current, regF);

    // post-process reference volumes
    for (int refno=0; refno < nr_ref; refno++)
    {
        postProcessVolume(Iref[refno]);
        if (do_esthetics)
            avg_origin += VOL_ELEM(Iref[refno](), 0, 0, 0);
    }

    // Adjust origin pixel of all references
    if (do_esthetics)
        for (int refno=0;refno<nr_ref; refno++)
            VOL_ELEM(Iref[refno](), 0, 0, 0) = avg_origin/(double)nr_ref;

#ifdef DEBUG
    std::cerr<<"finished maximization"<<std::endl;
#endif 
}

// Apply regularization
bool Prog_ml_tomo_prm::regularize(double sumw_allrefs)
{
    if (reg_current > 0.)
    {
#ifdef DEBUG
        std::cerr<<"start regularize"<<std::endl;
#endif 

        // Normalized regularization (in N/K)
        double reg_norm = reg_current * sumw_allrefs / nr_ref;
        Matrix3D<std::complex<double> > Fref, Favg, Fzero(dim,dim,hdim+1);
        Matrix3D<double> Mavg, Mdiff;
        double sum_diff2=0.;

        // Calculate  FT of average of all references
        // and sum of squared differences between all references
        for (int refno = 0; refno < nr_ref; refno++)
        {
            if (refno==0) 
                Mavg = Iref[refno]();
            else
                Mavg += Iref[refno]();
            for (int refno2 = 0; refno2 < nr_ref; refno2++)
            {
                Mdiff = Iref[refno]() - Iref[refno2]();
                sum_diff2 += Mdiff.sum2();
            }
        }
        transformer.FourierTransform(Mavg,Favg,true);
        Favg *= reg_norm;

        // Update the regularized references
        for (int refno = 0; refno < nr_ref; refno++)
        {
            transformer.FourierTransform(Iref[refno](),Fref,true);
            double sumw = alpha_k(refno) * sumw_allrefs;
            // Fref = (sumw*Fref + reg_norm*Favg) /  (sumw + nr_ref * reg_norm)
#define DEBUG_REGULARISE
#ifdef DEBUG_REGULARISE
            std::cerr<<"refno= "<<refno<<" sumw = "<<sumw<<" reg_norm= "<<reg_norm<<" Fref1= "<<DIRECT_MULTIDIM_ELEM(Fref,1) <<" Favg1= "<<DIRECT_MULTIDIM_ELEM(Favg,1)<<" (sumw + nr_ref * reg_norm)= "<<(sumw + nr_ref * reg_norm)<<std::endl;
#endif
            Fref *= sumw;
            Fref += Favg;
            Fref /= (sumw + nr_ref * reg_norm);
            std::cerr<<" newFref1= "<<DIRECT_MULTIDIM_ELEM(Fref,1) <<std::endl;
        }

        // Update the regularized sigma_noise estimate
        if (!fix_sigma_noise)
        {
            double reg_sigma_noise2 = sigma_noise*sigma_noise*sumw_allrefs*ddim3;
#ifdef DEBUG_REGULARISE
            std::cerr<<"reg_sigma_noise2= "<<reg_sigma_noise2<<" sumw_allrefs="<<sumw_allrefs<<" ddim3= "<<ddim3<<" sigma_noise= "<<sigma_noise<<" sum_diff2= "<<sum_diff2<<" reg_norm= "<<reg_norm<<std::endl;
#endif
            reg_sigma_noise2 += reg_norm * sum_diff2;
            sigma_noise = sqrt(reg_sigma_noise2/(sumw_allrefs*ddim3));
            std::cerr<<"new sigma_noise= "<<sigma_noise<<std::endl;
        }

#ifdef DEBUG
        std::cerr<<"finished regularize"<<std::endl;
#endif 
    }
}
// Check convergence
bool Prog_ml_tomo_prm::checkConvergence(std::vector<double> &conv)
{

#ifdef DEBUG
    std::cerr<<"started checkConvergence"<<std::endl;
#endif 
    bool converged = true;
    double convv;
    Matrix3D<double> Maux;

    Maux.resize(dim, dim, dim);
    Maux.setXmippOrigin();

    conv.clear();
    for (int refno=0;refno<nr_ref; refno++)
    {
        if (alpha_k(refno) > 0.)
        {
            Maux = Iold[refno]() * Iold[refno]();
            convv = 1. / (Maux.computeAvg());
            Maux = Iold[refno]() - Iref[refno]();
            Maux = Maux * Maux;
            convv *= Maux.computeAvg();
            conv.push_back(convv);
            if (convv > eps) converged = false;
        }
        else
        {
            conv.push_back(-1.);
        }
    }

#ifdef DEBUG
    std::cerr<<"finished checkConvergence"<<std::endl;
#endif 
    return converged;
}

void Prog_ml_tomo_prm::writeOutputFiles(const int iter, DocFile &DFo,
                                        std::vector<Matrix3D<double> > &wsumweds,
                                        double &sumw_allrefs, double &LL, double &avefracweight, 
                                        std::vector<double> &conv)
{

    FileName          fn_tmp, fn_base, fn_tmp2;
    Matrix1D<double>  fracline(2);
    SelFile           SFo, SFc;
    DocFile           DFl;
    std::string       comment;
    std::ofstream     fh;
    VolumeXmipp       Vt;

    DFl.clear();
    SFo.clear();
    SFc.clear();

    fn_base = fn_root;
    if (iter >= 0)
    {
        fn_base += "_it";
        fn_base.compose(fn_base, iter, "");
    }

    if (do_norm) fracline.resize(4);
    // Write out current reference images and fill sel & log-file
    for (int refno=0;refno<nr_ref; refno++)
    {
        fn_tmp = fn_base + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "");
        fn_tmp = fn_tmp + ".vol";
        Iref[refno].write(fn_tmp);
        SFo.insert(fn_tmp, SelLine::ACTIVE);
        fracline(0) = alpha_k(refno);
        fracline(1) = 1000 * conv[refno]; // Output 1000x the change for precision
        DFl.insert_comment(fn_tmp);
        DFl.insert_data_line(fracline);

        if (iter >= 1 && do_missing)
        {

            double sumw = alpha_k(refno)*sumw_allrefs;
            Vt().resize(dim,dim,dim);
            Vt().setXmippOrigin();
            Vt().initZeros();

            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Vt())
                if (j<hdim+1)
                    DIRECT_VOL_ELEM(Vt(),k,i,j)=
                        DIRECT_VOL_ELEM(wsumweds[refno],k,i,j) / sumw_allrefs;
                else
                    DIRECT_VOL_ELEM(Vt(),k,i,j)=
                        DIRECT_VOL_ELEM(wsumweds[refno],
                                        (dim-k)%dim,
                                        (dim-i)%dim,
                                        dim-j) / sumw_allrefs;

            CenterFFT(Vt(),true);
            fn_tmp = fn_base + "_wedge";
            fn_tmp.compose(fn_tmp, refno + 1, "");
            fn_tmp = fn_tmp + ".vol";
            Vt.write(fn_tmp);
        }
    }


    // Write out sel & log-file
    fn_tmp = fn_base + ".sel";
    SFo.write(fn_tmp);

    DFl.go_beginning();
    comment = "ml_tomo-logfile: Number of images= " + floatToString(sumw_allrefs);
    comment += " LL= " + floatToString(LL, 15, 10) + " <Pmax/sumP>= " + floatToString(avefracweight, 10, 5);
    if (do_norm)
        comment+= " <scale>= " + floatToString(average_scale, 10, 5);
    DFl.insert_comment(comment);
    comment = "-noise " + floatToString(sigma_noise, 15, 12) + " -offset " + floatToString(sigma_offset, 15, 12) + " -istart " + integerToString(iter + 1) + " -doc " + fn_base + ".doc";
    DFl.insert_comment(comment);
    DFl.insert_comment(cline);
    DFl.insert_comment("columns: model fraction (1); 1000x signal change (2)");
    fn_tmp = fn_base + ".log";
    DFl.write(fn_tmp);

    // Write out docfile with optimal transformation & references
    fn_tmp = fn_base + ".doc";
    DFo.write(fn_tmp);

    // Also write out selfiles of all experimental images,
    // classified according to optimal reference image
    if (iter >= 1)
    {
        for (int refno = 0;refno < nr_ref; refno++)
        {
            DFo.go_beginning();
            SFo.clear();
            for (int n = 0; n < DFo.dataLineNo(); n++)
            {
                DFo.next();
                fn_tmp = ((DFo.get_current_line()).get_text()).erase(0, 3);
                DFo.adjust_to_data_line();
                if ((refno + 1) == (int)DFo(6)) SFo.insert(fn_tmp, SelLine::ACTIVE);
            }
            fn_tmp.compose(fn_base + "_ref", refno + 1, "sel");
            SFo.write(fn_tmp);
        }
    }

}

