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
//#define DEBUG

//Mutex for each thread update sums
pthread_mutex_t update_mutex = PTHREAD_MUTEX_INITIALIZER;
//Mutex for each thread get next work
pthread_mutex_t work_mutex = PTHREAD_MUTEX_INITIALIZER;

// Read arguments ==========================================================
void Prog_MLalign2D_prm::read(int argc, char **argv, bool ML3D)
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
        if (strstr(comment.c_str(), "MLalign2D-logfile") == NULL)
        {
            std::cerr << "Error!! Docfile is not of MLalign2D-logfile type. " << std::endl;
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
    n_ref = textToInteger(getParameter(argc2, argv2, "-nref", "0"));
    fn_ref = getParameter(argc2, argv2, "-ref", "");
    fn_sel = getParameter(argc2, argv2, "-i");
    fn_root = getParameter(argc2, argv2, "-o", "ml2d");
    psi_step = textToFloat(getParameter(argc2, argv2, "-psi_step", "5"));
    Niter = textToInteger(getParameter(argc2, argv2, "-iter", "100"));
    istart = textToInteger(getParameter(argc2, argv2, "-istart", "1"));
    sigma_noise = textToFloat(getParameter(argc2, argv2, "-noise", "1"));
    sigma_offset = textToFloat(getParameter(argc2, argv2, "-offset", "3"));
    do_mirror = checkParameter(argc2, argv2, "-mirror");
    eps = textToFloat(getParameter(argc2, argv2, "-eps", "5e-5"));
    fn_frac = getParameter(argc2, argv2, "-frac", "");
    fix_fractions = checkParameter(argc2, argv2, "-fix_fractions");
    fix_sigma_offset = checkParameter(argc2, argv2, "-fix_sigma_offset");
    fix_sigma_noise = checkParameter(argc2, argv2, "-fix_sigma_noise");
    verb = textToInteger(getParameter(argc2, argv2, "-verb", "1"));
    fast_mode = checkParameter(argc2, argv2, "-fast");
    C_fast = textToFloat(getParameter(argc2, argv2, "-C", "1e-12"));
    save_mem1 = checkParameter(argc2, argv2, "-save_memA");
    save_mem2 = checkParameter(argc2, argv2, "-save_memB");
    fn_doc = getParameter(argc2, argv2, "-doc", "");
    zero_offsets = checkParameter(argc2, argv2, "-zero_offsets");
    do_student = checkParameter(argc2, argv2, "-student");
    df = (double) textToInteger(getParameter(argc2, argv2, "-df", "6"));
    do_norm = checkParameter(argc2, argv2, "-norm");
    do_ML3D = ML3D;

    // Number of threads
    threads = textToInteger(getParameter(argc2, argv2, "-thr","1"));

    // Hidden arguments
    fn_scratch = getParameter(argc2, argv2, "-scratch", "");
    debug = textToInteger(getParameter(argc2, argv2, "-debug","0"));
    do_student_sigma_trick = !checkParameter(argc2, argv2, "-no_sigma_trick");
    trymindiff_factor = textToFloat(getParameter(argc2, argv2, "-trymindiff_factor", "0.9"));
    
    // Only for interaction with refine3d:
    search_rot = textToFloat(getParameter(argc2, argv2, "-search_rot", "999."));

}

// Show ====================================================================
void Prog_MLalign2D_prm::show(bool ML3D)
{

    if (verb > 0)
    {
        // To screen
        if (!ML3D)
        {
            std::cerr << " -----------------------------------------------------------------" << std::endl;
            std::cerr << " | Read more about this program in the following publications:   |" << std::endl;
	    std::cerr << " |  Scheres ea. (2005) J.Mol.Biol. 348(1), 139-49                |" << std::endl;
	    std::cerr << " |  Scheres ea. (2005) Bioinform. 21(suppl.2), ii243-4   (-fast) |" << std::endl;
            std::cerr << " |                                                               |" << std::endl;
            std::cerr << " |  *** Please cite them if this program is of use to you! ***   |" << std::endl;
            std::cerr << " -----------------------------------------------------------------" << std::endl;
        }
        std::cerr << "--> Maximum-likelihood multi-reference refinement " << std::endl;
        std::cerr << "  Input images            : " << fn_sel << " (" << nr_exp_images << ")" << std::endl;
        if (fn_ref != "")
            std::cerr << "  Reference image(s)      : " << fn_ref << std::endl;
        else
            std::cerr << "  Number of references:   : " << n_ref << std::endl;
        std::cerr << "  Output rootname         : " << fn_root << std::endl;
        std::cerr << "  Stopping criterium      : " << eps << std::endl;
        std::cerr << "  initial sigma noise     : " << sigma_noise << std::endl;
        std::cerr << "  initial sigma offset    : " << sigma_offset << std::endl;
        std::cerr << "  Psi sampling interval   : " << psi_step << std::endl;
        if (do_mirror)
            std::cerr << "  Check mirrors           : true" << std::endl;
        else
            std::cerr << "  Check mirrors           : false" << std::endl;
        if (fn_frac != "")
            std::cerr << "  Initial model fractions : " << fn_frac << std::endl;
        if (fast_mode)
        {
            std::cerr << "  -> Use fast, reduced search-space approach with C = " << C_fast << std::endl;
            if (zero_offsets)
                std::cerr << "    + Start from all-zero translations" << std::endl;
        }
        if (search_rot < 180.)
            std::cerr << "    + Limit orientational search to +/- " << search_rot << " degrees" << std::endl;
        if (save_mem1)
            std::cerr << "  -> Save_memory A: recalculate real-space rotations in -fast" << std::endl;
        if (save_mem2)
            std::cerr << "  -> Save_memory B: limit translations to 3 sigma_offset " << std::endl;
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
        if (do_student)
        {
            std::cerr << "  -> Use t-student distribution with df = " <<df<< std::endl;
            if (do_student_sigma_trick)
            {
                std::cerr << "  -> Use sigma-trick for t-student distributions" << std::endl;
            }
        }
        if (do_norm)
        {
            std::cerr << "  -> Refine normalization for each experimental image"<<std::endl;
        }
        if (threads>1)
        {
            std::cerr << "  -> Using "<<threads<<" parallel threads"<<std::endl;
        }
         

        std::cerr << " -----------------------------------------------------------------" << std::endl;

    }

}

// Usage ===================================================================
void Prog_MLalign2D_prm::usage()
{
    std::cerr << "Usage:  ml_align2d [options] " << std::endl;
    std::cerr << "   -i <selfile>                : Selfile with input images \n";
    std::cerr << "   -nref <int>                 : Number of references to generate automatically (recommended)\n";
    std::cerr << "   OR -ref <selfile/image>         OR selfile with initial references/single reference image \n";
    std::cerr << " [ -o <rootname> ]             : Output rootname (default = \"ml2d\")\n";
    std::cerr << " [ -mirror ]                   : Also check mirror image of each reference \n";
    std::cerr << " [ -fast ]                     : Use pre-centered images to pre-calculate significant orientations\n";
    std::cerr << " [ -thr <N=1> ]                : Use N parallel threads \n";
    std::cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_MLalign2D_prm::extendedUsage(bool ML3D)
{
    std::cerr << "Additional options: " << std::endl;
    std::cerr << " [ -eps <float=5e-5> ]         : Stopping criterium \n";
    std::cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
    std::cerr << " [ -psi_step <float=5> ]       : In-plane rotation sampling interval [deg]\n";
    std::cerr << " [ -noise <float=1> ]          : Expected standard deviation for pixel noise \n";
    std::cerr << " [ -offset <float=3> ]         : Expected standard deviation for origin offset [pix]\n";
    std::cerr << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n";
    std::cerr << " [ -C <double=1e-12> ]         : Significance criterion for fast approach \n";
    std::cerr << " [ -zero_offsets ]             : Kick-start the fast algorithm from all-zero offsets \n";
    if (!ML3D) std::cerr << " [ -restart <logfile> ]        : restart a run with all parameters as in the logfile \n";
    if (!ML3D) std::cerr << " [ -istart <int> ]             : number of initial iteration \n";
    std::cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise \n";
    std::cerr << " [ -fix_sigma_offset]          : Do not re-estimate the standard deviation in the origin offsets \n";
    std::cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
    std::cerr << " [ -doc <docfile=\"\"> ]         : Read initial angles and offsets from docfile \n";
    std::cerr << " [ -student ]                  : Use t-distributed instead of Gaussian model for the noise \n";
    std::cerr << " [ -df <int=6> ]               : Degrees of freedom for the t-distribution \n";
    std::cerr << " [ -norm ]                     : Refined normalization parameters for each particle \n";
    std::cerr << std::endl;
    exit(1);
}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produce_Side_info2)
void Prog_MLalign2D_prm::produceSideInfo()
{

    Matrix2D<double>            A(3, 3);

    // Read selfile with experimental images
    SF.read(fn_sel);
    nr_exp_images = SF.ImgNo();

    // Get original image size
    SF.ImgSize(dim, dim);
    hdim = dim / 2;
    dim2 = dim * dim;
    ddim2 = (double)dim2;
    sigma_noise2 = sigma_noise * sigma_noise;
    if (do_student) 
    {
        df2 = - ( df + ddim2 ) / 2. ;
        dfsigma2 = df * sigma_noise2;
    }


    // prepare masks for rotated references
    mask.resize(dim, dim);
    mask.setXmippOrigin();
    BinaryCircularMask(mask, hdim, INNER_MASK);
    omask.resize(dim, dim);
    omask.setXmippOrigin();
    BinaryCircularMask(omask, hdim, OUTSIDE_MASK);
    
    // Get number of references
    if (do_ML3D)
    {
        do_generate_refs = false;
    }
    else if (fn_ref != "")
    {
        do_generate_refs = false;
        if (Is_ImageXmipp(fn_ref)) n_ref = 1;
        else
        {
            SFr.read(fn_ref);
            n_ref = SFr.ImgNo();
        }
    }
    else
    {
        do_generate_refs = true;
    }

    // Construct matrices for 0, 90, 180 & 270 degree flipping and mirrors
    psi_max = 90.;
    nr_psi = CEIL(psi_max / psi_step);
    psi_step = psi_max / nr_psi;
    nr_flip = nr_nomirror_flips = 4;
    A.initIdentity();
    F.push_back(A);
    A(0, 0) = 0.;
    A(1, 1) = 0.;
    A(1, 0) = 1.;
    A(0, 1) = -1;
    F.push_back(A);
    A(0, 0) = -1.;
    A(1, 1) = -1.;
    A(1, 0) = 0.;
    A(0, 1) = 0;
    F.push_back(A);
    A(0, 0) = 0.;
    A(1, 1) = 0.;
    A(1, 0) = -1.;
    A(0, 1) = 1;
    F.push_back(A);
    if (do_mirror)
    {
        nr_flip = 8;
        A.initIdentity();
        A(0, 0) = -1;
        F.push_back(A);
        A(0, 0) = 0.;
        A(1, 1) = 0.;
        A(1, 0) = 1.;
        A(0, 1) = 1;
        F.push_back(A);
        A(0, 0) = 1.;
        A(1, 1) = -1.;
        A(1, 0) = 0.;
        A(0, 1) = 0;
        F.push_back(A);
        A(0, 0) = 0.;
        A(1, 1) = 0.;
        A(1, 0) = -1.;
        A(0, 1) = -1;
        F.push_back(A);
    }

    // Set limit_rot
    limit_rot = (search_rot < 180.);

    // Set sigdim, i.e. the number of pixels that will be considered in the translations
    if (save_mem2)
    	sigdim = 2 * CEIL(sigma_offset * 3);
    else
    	sigdim = 2 * CEIL(sigma_offset * 6);
    sigdim++; // (to get uneven number)
    sigdim = XMIPP_MIN(dim, sigdim);

    //Some vectors and matrixes initialization
    refw.resize(n_ref);
    refw2.resize(n_ref);
    refwsc2.resize(n_ref);
    refw_mirror.resize(n_ref);
    sumw_refpsi.resize(n_ref*nr_psi);
    A2.resize(n_ref);
    fref.resize(n_ref * nr_psi);
    mref.resize(n_ref * nr_psi);
    wsum_Mref.resize(n_ref);
    if (fast_mode)
    {
    	imax = n_ref * nr_flip / nr_nomirror_flips;
        ioptx_ref.resize(imax);
        iopty_ref.resize(imax);
        ioptflip_ref.resize(imax);
        maxw_ref.resize(imax);
    }
    else
    {
        imax = 0;
    }

}

// Generate initial references =============================================
void Prog_MLalign2D_prm::generateInitialReferences()
{

    SelFile SFtmp, SFout;
    ImageXmipp Iave, Itmp;
    double dummy;
    FileName fn_tmp;

    if (verb > 0)
    {
        std::cerr << "  Generating initial references by averaging over random subsets" << std::endl;
        init_progress_bar(n_ref);
    }

    // Make random subsets and calculate average images
    SFtmp = SF.randomize();
    int Nsub = ROUND((double)SFtmp.ImgNo() / n_ref);
    for (int refno = 0; refno < n_ref; refno++)
    {
        SFout.clear();
        SFout.reserve(Nsub);
        SFtmp.go_beginning();
        SFtmp.jump_lines(Nsub*refno);
        if (refno == n_ref - 1) Nsub = SFtmp.ImgNo() - refno * Nsub;
        for (int nn = 0; nn < Nsub; nn++)
        {
            SFout.insert(SFtmp.current());
            SFtmp.NextImg();
        }
        SFout.get_statistics(Iave, Itmp, dummy, dummy);
        fn_tmp = fn_root + "_it";
        fn_tmp.compose(fn_tmp, 0, "");
        fn_tmp = fn_tmp + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "");
        fn_tmp = fn_tmp + ".xmp";
        Iave.write(fn_tmp);
        SFr.insert(fn_tmp, SelLine::ACTIVE);
        if (verb > 0) progress_bar(refno);
    }
    if (verb > 0) progress_bar(n_ref);
    fn_ref = fn_root + "_it";
    fn_ref.compose(fn_ref, 0, "sel");
    SFr.write(fn_ref);

}

// Read reference images to memory and initialize offset vectors
// This side info is NOT general, i.e. in parallel mode it is NOT the
// same for all processors! (in contrast to produceSideInfo)
void Prog_MLalign2D_prm::produceSideInfo2(int nr_vols)
{

    int                       c, idum, refno = 0;
    DocFile                   DF;
    DocLine                   DL;
    double                    offx, offy, aux, sumfrac = 0.;
    FileName                  fn_tmp;
    ImageXmipp                img;
    std::vector<double>       Vdum;

    // Read in all reference images in memory
    if (Is_ImageXmipp(fn_ref))
    {
        SFr.reserve(1);
        SFr.insert(fn_ref);
    }
    else
    {
        SFr.read(fn_ref);
    }
    n_ref = 0;
    SFr.go_beginning();
    while ((!SFr.eof()))
    {
        fn_tmp=SFr.NextImg();
        img.read(fn_tmp, false, false, true, false);
        img().setXmippOrigin();
        Iref.push_back(img);
        Iold.push_back(img);
        // Default start is all equal model fractions
        alpha_k.push_back((double)1 / SFr.ImgNo());
        Iref[refno].set_weight(alpha_k[refno] * (double)nr_exp_images);
        // Default start is half-half mirrored images
        if (do_mirror) mirror_fraction.push_back(0.5);
        else mirror_fraction.push_back(0.);
        n_ref++;
        refno++;
    }

    // Initialize trymindiff for all images
    imgs_trymindiff.clear();
    for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
    {
        imgs_optrefno.push_back(0);
        imgs_trymindiff.push_back(-1.);
    }

    // Initialize scale and bgmean for all images 
    // (for now initialize to 1 and 0, below also include doc)
    if (do_norm)
    {
	imgs_scale.clear();
        imgs_bgmean.clear();
        for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
        {
	    imgs_bgmean.push_back(0.);
	    imgs_scale.push_back(1.);
	}
        average_scale = 1.;
        refs_avgscale.clear();
        for (int refno=0;refno<n_ref; refno++)
        {
            refs_avgscale.push_back(1.);
        }
    }

    // Initialize imgs_offsets vectors
    if (zero_offsets) offx = 0.;
    else offx = -999.;
    if (do_mirror) idum = 4 * n_ref;
    else idum = 2 * n_ref;
    for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
    {
        imgs_offsets.push_back(Vdum);
        for (int refno = 0; refno < idum; refno++)
        {
            imgs_offsets[imgno].push_back(offx);
        }
    }

    // For limited orientational search: initialize imgs_oldphi & imgs_oldtheta to -999.
    if (limit_rot)
    {
        imgs_oldphi.clear();
        imgs_oldtheta.clear();
        for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
        {
            imgs_oldphi.push_back(-999.);
            imgs_oldtheta.push_back(-999.);
        }
    }


    // Read optimal image-parameters from fn_doc
    if (fn_doc != "")
    {
        if (limit_rot || zero_offsets || do_norm)
        {
            DF.read(fn_doc);
            DF.go_beginning();
            SF.go_beginning();
            int imgno = 0;
            while (!SF.eof())
            {
                fn_tmp = SF.NextImg();
                if (DF.search_comment(fn_tmp))
                {
                    if (limit_rot)
                    {
                        imgs_oldphi[imgno] = DF(0);
                        imgs_oldtheta[imgno] = DF(1);
                    }
                    if (zero_offsets)
                    {
                        if (do_mirror) idum = 2 * n_ref;
                        else idum = n_ref;
                        for (int refno = 0; refno < idum; refno++)
                        {
                            imgs_offsets[imgno][2*refno] = DF(3);
                            imgs_offsets[imgno][2*refno+1] = DF(4);
                        }
                    }
                    if (do_norm)
                    {
                        imgs_bgmean[imgno] = DF(9);
                        imgs_scale[imgno] = DF(10);
                    }
                }
                else
                {
                    REPORT_ERROR(1, (std::string)"Prog_MLalign2D_prm: Cannot find " + fn_tmp + " in docfile " + fn_doc);
                }
                imgno++;
            }
            DF.clear();
        }
    }

    // read in model fractions if given on command line
    if (fn_frac != "")
    {
        DF.read(fn_frac);
        DF.go_first_data_line();
        for (refno = 0; refno < n_ref; refno++)
        {
            DL = DF.get_current_line();
            alpha_k[refno] = DL[0];
            if (do_mirror)
            {
                if (DL[1] > 1. || DL[1] < 0.)
                    REPORT_ERROR(1, "Prog_MLalign2D_prm: Mirror fraction (2nd column) should be [0,1]!");
                mirror_fraction[refno] = DL[1];
            }
            if (do_norm)
            {
                refs_avgscale[refno] = DL[3];
            }
            sumfrac += alpha_k[refno];
            DF.next_data_line();
        }
        if (ABS(sumfrac - 1.) > 1e-3)
            if (verb > 0) std::cerr << " ->WARNING: Sum of all expected model fractions (" << sumfrac << ") is not one!" << std::endl;
        for (refno = 0; refno < n_ref; refno++)
        {
            alpha_k[refno] /= sumfrac;
        }
    }

}

// Calculate probability density function of all in-plane transformations phi
void Prog_MLalign2D_prm::calculatePdfInplane()
{

    double x, y, r2, pdfpix, sum;
    P_phi.resize(dim, dim);
    P_phi.setXmippOrigin();
    Mr2.resize(dim, dim);
    Mr2.setXmippOrigin();

    sum=0.;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(P_phi)
    {
        x = (double)j;
        y = (double)i;
        r2 = x * x + y * y;
        if (sigma_offset > 0.)
        {
            pdfpix = exp(-r2 / (2 * sigma_offset * sigma_offset));
            pdfpix /= 2 * PI * sigma_offset * sigma_offset * nr_psi * nr_nomirror_flips;
        }
        else
        {
            if (j == 0 && i == 0) pdfpix = 1.;
            else pdfpix = 0.;
        }
        MAT_ELEM(P_phi, i, j) = pdfpix;
        MAT_ELEM(Mr2, i, j) = (float)r2;
	sum+=pdfpix;
    }
    // Normalization
    P_phi/=sum;

}

// Rotate reference for all models and rotations and fill Fref vectors =============
void Prog_MLalign2D_prm::rotateReference()
{

#ifdef DEBUG
    std::cerr<<"entering rotateReference"<<std::endl;
#endif

    awakeThreads(THREAD_ROTATE_REFERENCE_REFNO, 0);

#ifdef DEBUG
    std::cerr<<"leaving rotateReference"<<std::endl;
#endif
}

// Collect all rotations and sum to update Iref() for all models ==========
void Prog_MLalign2D_prm::reverseRotateReference()
{

#ifdef DEBUG
    std::cerr<<"entering reverseRotateReference"<<std::endl;
#endif

    awakeThreads(THREAD_REVERSE_ROTATE_REFERENCE_REFNO, 0);

#ifdef DEBUG
    std::cerr<<"leaving reverseRotateReference"<<std::endl;
#endif

}

void Prog_MLalign2D_prm::preselectLimitedDirections(float &phi, float &theta)
{

    float phi_ref, theta_ref, angle, angle2;
    Matrix1D<double> u, v;

    pdf_directions.clear();
    pdf_directions.resize(n_ref);

    for (int refno = 0; refno < n_ref; refno++)
    {
        if (!limit_rot || (phi == -999. && theta == -999.))
        	pdf_directions[refno] = 1.;
        else
        {
            phi_ref = Iref[refno].Phi();
            theta_ref = Iref[refno].Theta();
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
            if (fabs(angle) > search_rot) pdf_directions[refno] = 0.;
            else pdf_directions[refno] = 1.;
        }
    }

}

// Pre-selection of significant refno and ipsi, based on current optimal translation =======
void Prog_MLalign2D_prm::preselectFastSignificant()
{

#ifdef DEBUG
    std::cerr<<"entering preselectFastSignificant"<<std::endl;
#endif

    // Initialize Msignificant to all zeros
    // TODO: check whether this is strictly necessary? Probably not... 
    Msignificant.initZeros();
    awakeThreads(THREAD_PRESELECT_FAST_SIGNIFICANT_REFNO, 0);

#ifdef DEBUG
    std::cerr<<"leaving preselectFastSignificant"<<std::endl;
#endif


}

// Maximum Likelihood calculation for one image ============================================
// Integration over all translation, given  model and in-plane rotation
void Prog_MLalign2D_prm::expectationSingleImage(
    Matrix1D<double> &opt_offsets)
{

    Matrix2D<double> Maux, Mweight;
    Matrix2D<std::complex<double> > Faux, Fzero(dim, hdim + 1);
    double scale_dim2_sumw, my_mindiff;
    bool is_ok_trymindiff = false;
    XmippFftw local_transformer;
    ioptx = iopty = 0;

    // Setup matrices
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Mweight.initZeros(sigdim, sigdim);
    Mweight.setXmippOrigin();
    Fzero.initZeros();

    if (!do_norm) 
        opt_scale =1.;

    // precalculate all flipped versions of the image
    Fimg_flip.clear();
    for (int iflip = 0; iflip < nr_flip; iflip++)
    {
        Maux.setXmippOrigin();
        applyGeometry(Maux, F[iflip], Mimg, IS_INV, WRAP);
        local_transformer.FourierTransform(Maux,Faux,false);
        if (do_norm)
            dMij(Faux,0,0) -= bgmean;
         Fimg_flip.push_back(Faux);
    }

    // The real stuff: loop over all references, rotations and translations
    int redo_counter = 0;
    while (!is_ok_trymindiff)
    {
        // Initialize mindiff, weighted sums and maxweights
        mindiff = 99.e99;
        wsum_corr = wsum_offset = wsum_sc = wsum_sc2 = 0.;
        maxweight = maxweight2 = sum_refw = sum_refw2 = 0.;

        mysumimgs.clear();
        for (int i = 0; i < n_ref*nr_psi; i++)
        {
            mysumimgs.push_back(Fzero); 
            sumw_refpsi[i] = 0.;
        }

        for (int i = 0; i < n_ref; i++)
            refw[i] = refw2[i] = refw_mirror[i] = 0.;

        if (fast_mode)
        {
            for (int i = 0; i < imax; i++)
            {
                maxw_ref[i] = -99.e99;
                ioptx_ref[i] = 0;
                iopty_ref[i] = 0;
                ioptflip_ref[i] = 0;
            }
        }

        awakeThreads(THREAD_EXPECTATION_SINGLE_IMAGE_REFNO, opt_refno);

        // Now check whether our trymindiff was OK.
        // The limit of the exp-function lies around 
        // exp(700)=1.01423e+304, exp(800)=inf; exp(-700) = 9.85968e-305; exp(-88) = 0
        // Use 500 to be on the save side?
        if (ABS((mindiff - trymindiff) / sigma_noise2) > 500.)
        {
            // Re-do whole calculation now with the real mindiff
            trymindiff = mindiff;
            redo_counter++;
            // Never re-do more than once!
            if (redo_counter>1)
            {
                std::cerr<<"ml_align2d BUG% redo_counter > 1"<<std::endl;
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

    // Calculate optimal transformation parameters
    opt_psi = -psi_step * (iopt_flip * nr_psi + iopt_psi) - SMALLANGLE;
    if (fast_mode)
        for (int i = 0; i < imax;i++)
        {
            allref_offsets[2*i] = -(double)ioptx_ref[i] * DIRECT_MAT_ELEM(F[ioptflip_ref[i]], 0, 0) -
                                   (double)iopty_ref[i] * DIRECT_MAT_ELEM(F[ioptflip_ref[i]], 0, 1);
            allref_offsets[2*i+1] = -(double)ioptx_ref[i] * DIRECT_MAT_ELEM(F[ioptflip_ref[i]], 1, 0) -
                                     (double)iopty_ref[i] * DIRECT_MAT_ELEM(F[ioptflip_ref[i]], 1, 1);
        }
    opt_offsets(0) = -(double)ioptx * DIRECT_MAT_ELEM(F[iopt_flip], 0, 0) -
                      (double)iopty * DIRECT_MAT_ELEM(F[iopt_flip], 0, 1);
    opt_offsets(1) = -(double)ioptx * DIRECT_MAT_ELEM(F[iopt_flip], 1, 0) -
                      (double)iopty * DIRECT_MAT_ELEM(F[iopt_flip], 1, 1);

    // Update normalization parameters
    if (do_norm)
    {
        // 1. Calculate optimal setting of Mimg
        Matrix2D<double> Maux2 = Mimg;
        Maux2.selfTranslate(opt_offsets, true);
        Maux2.selfApplyGeometry(F[iopt_flip], IS_INV, WRAP);
        // 2. Calculate optimal setting of Mref
        int refnoipsi = opt_refno*nr_psi + iopt_psi;
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
        {
            dMij(Faux,i,j) = conj(dMij(fref[refnoipsi],i,j));
            dMij(Faux,i,j) *= opt_scale;
        }
        // Still take input from Faux and leave output in Maux
        local_transformer.inverseFourierTransform();
        Maux2 = Maux2 - Maux;
	if (debug==12) 
        {  
            std::cout<<std::endl;
            std::cout<<"scale= "<<opt_scale<<" changes to "<<wsum_sc/wsum_sc2<<std::endl;
            std::cout<<"bgmean= "<<bgmean<<" changes to "<<Maux2.computeAvg()<<std::endl;
        }
        // non-ML update of bgmean (this is much cheaper than true-ML update...)
        old_bgmean = bgmean;
        bgmean = Maux2.computeAvg();
        // ML-update of opt_scale
        opt_scale = wsum_sc / wsum_sc2;
    }

    // Update all global weighted sums after division by sum_refw
    wsum_sigma_noise += (2 * wsum_corr / sum_refw);
    wsum_sigma_offset += (wsum_offset / sum_refw);
    sumfracweight += fracweight;
    scale_dim2_sumw = (opt_scale * ddim2) / sum_refw;
    for (int refno = 0; refno < n_ref; refno++)
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            sumw[refno] += (refw[refno] + refw_mirror[refno]) / sum_refw;
	    sumw2[refno] += refw2[refno] / sum_refw;
            sumw_mirror[refno] += refw_mirror[refno] / sum_refw;
            if (do_student)
            {
                sumwsc[refno] += refw2[refno] * (opt_scale) / sum_refw;
                sumwsc2[refno] += refw2[refno] * (opt_scale * opt_scale) / sum_refw;
            }
            else
            {
                sumwsc[refno] += (refw[refno] + refw_mirror[refno]) * (opt_scale) / sum_refw;
                sumwsc2[refno] += (refw[refno] + refw_mirror[refno]) * (opt_scale * opt_scale) / sum_refw;
            }
            for (int ipsi = 0; ipsi < nr_psi; ipsi++)
            {
                int refnoipsi = refno*nr_psi + ipsi;
                // Correct weighted sum of images for new bgmean (only first element=origin in Fimg)
                if (do_norm)
                    dMij(mysumimgs[refnoipsi],0,0) -= sumw_refpsi[refnoipsi] * (bgmean - old_bgmean) / ddim2; 
                // Sum mysumimgs to the global weighted sum
                wsumimgs[refnoipsi] += (scale_dim2_sumw * mysumimgs[refnoipsi]);
            }
        }
    }

    if (!do_student)
	// 1st term: log(refw_i)
	// 2nd term: for subtracting mindiff
	// 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
        dLL = log(sum_refw) 
            - my_mindiff / sigma_noise2 
            - ddim2 * log( sqrt(2. * PI * sigma_noise2));
    else
	// 1st term: log(refw_i)
	// 2nd term: for dividing by (1 + 2. * mindiff/dfsigma2)^df2
	// 3rd term: for sigma-dependent normalization term in t-student distribution
	// 4th&5th terms: gamma functions in t-distribution
        dLL = log(sum_refw) 
            + df2 * log( 1. + (  2. * my_mindiff / dfsigma2 )) 
            - ddim2 * log( sqrt(PI * df * sigma_noise2)) 
            + gammln(-df2) - gammln(df/2.);
    LL += dLL;


}//close function expectationSingleImage

/** Function to create threads that will work later */
void Prog_MLalign2D_prm::createThreads()
{

	//Initialize some variables for using for threads

	th_ids = new pthread_t[threads];
	threads_d = new structThreadTasks[threads];

	barrier_init( &barrier, threads + 1);
	barrier_init( &barrier2, threads + 1);

	for (int i = 0; i < threads; i++)
	{
		threads_d[i].thread_id = i;
		threads_d[i].prm = this;

		int result = pthread_create( (th_ids + i), NULL, doThreadsTasks, (void *)(threads_d+i));

		if (result != 0)
		{
			REPORT_ERROR(result, "ERROR CREATING THREAD");
		}
	}

}//close function createThreads

/** Free threads memory and exit */
void Prog_MLalign2D_prm::destroyThreads()
{
	threadTask = THREAD_EXIT;
	barrier_wait(&barrier);
	delete [] th_ids;
	delete [] threads_d;
}

/// Function for threads do different tasks
void * doThreadsTasks(void * data)
{
    structThreadTasks * thread_data = (structThreadTasks *) data;
    
    int thread_id = thread_data->thread_id;
    Prog_MLalign2D_prm * prm = thread_data->prm;
    
    barrier_t & barrier = prm->barrier;
    barrier_t & barrier2 = prm->barrier2;
    
    //Loop until the threadTask become THREAD_EXIT
    do
    {
        //Wait until main threads order to start
        //debug_print("Waiting on barrier, th ", thread_id);
        barrier_wait(&barrier);
        
        //Check task to do
        switch (prm->threadTask)
        {
        case THREAD_PRESELECT_FAST_SIGNIFICANT_REFNO:
            prm->doThreadPreselectFastSignificantRefno();
            break;
 
        case THREAD_EXPECTATION_SINGLE_IMAGE_REFNO:
            prm->doThreadExpectationSingleImageRefno();
            break;
            
        case THREAD_ROTATE_REFERENCE_REFNO:
            prm->doThreadRotateReferenceRefno();
            break;
            
        case THREAD_REVERSE_ROTATE_REFERENCE_REFNO:
            prm->doThreadReverseRotateReferenceRefno();
            break;
            
        case THREAD_EXIT:
            pthread_exit(NULL);
            break;
            
        }
        
        barrier_wait(&barrier2);
        
    }while (1);
    
}//close function doThreadsTasks


/// Function to assign refno jobs to threads
/// return -1 when all jobs are done
int Prog_MLalign2D_prm::getThreadRefnoJob()
{
    int refno = -1;

    if (refno_count < n_ref)
    {
        refno = refno_index;
        refno_index = (refno_index + 1) % n_ref;
        refno_count++;
    }

    return refno;
}//close function getThreadRefnoJob

///Function for awake threads for different tasks
void Prog_MLalign2D_prm::awakeThreads(int task, int start_refno)
{
    threadTask = task;
    refno_index = start_refno;
    refno_count = 0;
    barrier_wait(&barrier);
    //Wait until done
    barrier_wait(&barrier2);
}//close function awakeThreads


void Prog_MLalign2D_prm::doThreadRotateReferenceRefno()
{
    double AA, stdAA, psi, dum, avg;
    Matrix2D<double> Maux(dim,dim);
    Matrix2D<std::complex<double> > Faux;
    XmippFftw local_transformer;
    int refnoipsi;
    
    Maux.setXmippOrigin();
    
    pthread_mutex_lock(&work_mutex);
    int refno = getThreadRefnoJob();
    pthread_mutex_unlock(&work_mutex);
    
    while (refno != -1)
    {
        computeStats_within_binary_mask(omask, Iref[refno](), dum, dum, avg, dum);
        for (int ipsi = 0; ipsi < nr_psi; ipsi++)
        {
            refnoipsi = refno*nr_psi + ipsi;
            // Add arbitrary number (small_angle) to avoid 0-degree rotation (lacking interpolation)
            psi = (double)(ipsi * psi_max / nr_psi) + SMALLANGLE;
            Iref[refno]().rotateBSpline(3, psi, Maux, WRAP);
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
            local_transformer.FourierTransform(Maux,Faux,false);
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
            {
                dMij(Faux, i, j) = conj(dMij(Faux, i, j));
            }
            
            fref[refnoipsi] = Faux;
        }
        // If we dont use save_mem1 Iref[refno] is useless from here on
        if (!save_mem1)
            Iref[refno]().resize(0, 0);
        
        pthread_mutex_lock(&work_mutex);
        refno = getThreadRefnoJob();
        pthread_mutex_unlock(&work_mutex);
    }//close while
    
}//close function doThreadRotateReferenceRefno

void Prog_MLalign2D_prm::doThreadReverseRotateReferenceRefno()
{
    double psi, dum, avg, ang;
    Matrix2D<double> Maux(dim, dim), Maux2(dim, dim), Maux3(dim, dim);
    Matrix2D<std::complex<double> > Faux;
    XmippFftw local_transformer;

    Maux.setXmippOrigin();
    Maux2.setXmippOrigin();
    Maux3.setXmippOrigin();

    pthread_mutex_lock(&work_mutex);
    int refno = getThreadRefnoJob();
    pthread_mutex_unlock(&work_mutex);
    
    while (refno != -1)
    {
        Maux.initZeros();
        wsum_Mref[refno] = Maux;
        
        for (int ipsi = 0; ipsi < nr_psi; ipsi++)
        {
            // Add arbitrary number to avoid 0-degree rotation without interpolation effects
            psi = (double)(ipsi * psi_max / nr_psi) + SMALLANGLE;
            int refnoipsi = refno * nr_psi + ipsi;
            // Do the backward FFT
            // The construction with Faux and Maux3 should perhaps not be necessary,
            // but I am having irreproducible segmentation faults for the iFFT
            Faux = wsumimgs[refnoipsi];
            local_transformer.inverseFourierTransform(Faux, Maux);
            Maux3 = Maux;
            CenterFFT(Maux3, true);
            computeStats_within_binary_mask(omask, Maux3, dum, dum, avg, dum);
            Maux3.rotateBSpline(3, -psi, Maux2, WRAP);
            apply_binary_mask(mask, Maux2, Maux2, avg);
            wsum_Mref[refno] += Maux2;
        }

        pthread_mutex_lock(&work_mutex);
        refno = getThreadRefnoJob();
        pthread_mutex_unlock(&work_mutex);
    }//close while refno
}//close function doThreadReverseRotateReference

void Prog_MLalign2D_prm::doThreadPreselectFastSignificantRefno()
{
    Matrix2D<double> Mtrans, Mflip;
    double ropt, aux, diff, pdf, fracpdf;
    double A2_plus_Xi2;
    int irot, irefmir, iiflip;
    Matrix1D<double> trans(2);
    Matrix1D<double> weight(nr_psi*nr_flip);
    double local_maxweight, local_mindiff;
    int nr_mirror = 1;

    if (do_mirror)
        nr_mirror++;

    Mtrans.resize(dim, dim);
    Mtrans.setXmippOrigin();
    Mflip.resize(dim, dim);
    Mflip.setXmippOrigin();

    pthread_mutex_lock(&work_mutex);
    int refno = getThreadRefnoJob();
    pthread_mutex_unlock(&work_mutex);
    
    while (refno != -1)
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            A2_plus_Xi2 = 0.5 * (A2[refno] + Xi2);
            for (int imirror = 0; imirror < nr_mirror; imirror++)
            {
                local_maxweight = -99.e99;
                local_mindiff = 99.e99;
                irefmir = imirror * n_ref + refno;
                // Get optimal offsets
                trans(0) = (double)allref_offsets[2*irefmir];
                trans(1) = (double)allref_offsets[2*irefmir+1];
                ropt = sqrt(trans(0) * trans(0) + trans(1) * trans(1));
                // Do not trust optimal offsets if they are larger than 3*sigma_offset:
                if (ropt > 3 * sigma_offset)
                {
                    for (int iflip = 0; iflip < nr_nomirror_flips; iflip++)
                    {
                        iiflip = imirror*nr_nomirror_flips + iflip;
                        for (int ipsi = 0; ipsi < nr_psi; ipsi++)
                        {
                            irot = iiflip * nr_psi + ipsi;
                            dMij(Msignificant, refno, irot) = 1;
                        }
                    }
                }
                else
                {
                    // Set priors
                    if (imirror == 0)
                        fracpdf = alpha_k[refno] * (1. - mirror_fraction[refno]);
                    else
                        fracpdf = alpha_k[refno] * mirror_fraction[refno];
                    pdf = fracpdf * MAT_ELEM(P_phi, 
                                             (int)trans(1),
                                             (int)trans(0));

                    // A. Translate image and calculate probabilities for every rotation
                    Mimg.translate(trans, Mtrans, true);
                    for (int iflip = 0; iflip < nr_nomirror_flips; iflip++)
                    {
                        iiflip = imirror*nr_nomirror_flips + iflip;
                        applyGeometry(Mflip, F[iiflip], Mtrans, IS_INV, WRAP);
                        for (int ipsi = 0; ipsi < nr_psi; ipsi++)
                        {
                            irot = iiflip * nr_psi + ipsi;
                            diff = A2_plus_Xi2;
                            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mflip)
                            {
                                diff -= dMij(Mflip, i, j) * dMij(mref[refno*nr_psi + ipsi], i, j);
                            }
                            weight(irot) = diff;
                            if (diff < local_mindiff)
                                local_mindiff = diff;
                        }
                    }

                    // B. Now that we have local_mindiff, calculate the weights
                    for (int iflip = 0; iflip < nr_nomirror_flips; iflip++)
                    {
                        iiflip = imirror*nr_nomirror_flips + iflip;
                        for (int ipsi = 0; ipsi < nr_psi; ipsi++)
                        {
                            irot = iiflip * nr_psi + ipsi;
                            diff = weight(irot);
                            if (!do_student)
                            {
                                // normal distribution
                                aux = (diff - local_mindiff) / sigma_noise2;
                                // next line because of numerical precision of exp-function
                                if (aux > 1000.) weight(irot) = 0.;
                                else weight(irot) = exp(-aux) * pdf;
                            }
                            else
                            {
                                // t-student distribution
                                aux = (dfsigma2 + 2. * diff) / (dfsigma2 + 2. * local_mindiff);
                                weight(irot) = pow(aux, df2) * pdf;
                            }
                            if (weight(irot) > local_maxweight)
                                local_maxweight = weight(irot);
                        } // close ipsi
                    } // close iflip

                    // C. Now that we know the weights for all rotations, set Msignificant
                    for (int iflip = 0; iflip < nr_nomirror_flips; iflip++)
                    {
                        iiflip = imirror*nr_nomirror_flips + iflip;
                        for (int ipsi = 0; ipsi < nr_psi; ipsi++)
                        {
                            irot = iiflip * nr_psi + ipsi;
                            if (weight(irot) >= C_fast * local_maxweight)
                                dMij(Msignificant, refno, irot) = 1; 
                            else
                                dMij(Msignificant, refno, irot) = 0; 
                        }
                    }
                } //endif ropt<3*sigma_offset
            } //end loop imirror
        } //endif limit_rot and pdf_directions
        pthread_mutex_lock(&work_mutex);
        refno = getThreadRefnoJob();
        pthread_mutex_unlock(&work_mutex);
    }//close while refno

}//close function doThreadPreselectFastSignificantRefno

void Prog_MLalign2D_prm::doThreadExpectationSingleImageRefno()
{
    double diff;
    int old_optrefno = opt_refno;
    double XiA, aux, pdf, fracpdf, A2_plus_Xi2;
    double weight, stored_weight, weight1, weight2, my_maxweight;
    double my_sumweight, my_sumstoredweight, ref_scale = 1.;
    int irot, irefmir, refnoipsi;
    //Some local variables to store partial sums of global sums variables
    double local_mindiff, local_wsum_corr, local_wsum_offset;
    double local_wsum_sc, local_wsum_sc2, local_maxweight;
    int local_iopty, local_ioptx, local_iopt_psi, local_iopt_flip;
    
    //FIXME: this will not be here
    Matrix2D<double> Maux, Mweight;
    Matrix2D<std::complex<double> > Faux, Fzero(dim, hdim + 1);
    XmippFftw local_transformer;

    // Setup matrices
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Mweight.initZeros(sigdim, sigdim);
    Mweight.setXmippOrigin();
    Fzero.initZeros();
    local_transformer.setReal(Maux);
    local_transformer.getFourierAlias(Faux);
    
    // Start the loop over all refno at old_optrefno (=opt_refno from the previous iteration).
    // This will speed-up things because we will find Pmax probably right away,
    // and this will make the if-statement that checks SIGNIFICANT_WEIGHT_LOW
    // effective right from the start
    pthread_mutex_lock(&work_mutex);
    int refno = getThreadRefnoJob();
    pthread_mutex_unlock(&work_mutex);
    
    while (refno != -1)
    {
        refw[refno] = refw_mirror[refno] = 0.;
        local_mindiff = 99.e99;
        local_wsum_sc = local_wsum_sc2 = local_wsum_corr = local_wsum_offset = 0;
        // This if is for limited rotation options
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            if (do_norm)
                ref_scale = opt_scale / refs_avgscale[refno];
            A2_plus_Xi2 = 0.5 * (ref_scale*ref_scale*A2[refno] + Xi2);

            for (int iflip = 0; iflip < nr_flip; iflip++)
            {
                for (int ipsi = 0; ipsi < nr_psi; ipsi++)
                {
                    refnoipsi = refno*nr_psi + ipsi;
                    irot = iflip * nr_psi + ipsi;
                    irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                    // This if is the speed-up caused by the -fast options
                    if (dMij(Msignificant, refno, irot))
                    {
                        if (iflip < nr_nomirror_flips)
                            fracpdf = alpha_k[refno] * (1. - mirror_fraction[refno]);
                        else
                            fracpdf = alpha_k[refno] * mirror_fraction[refno];

                        // A. Backward FFT to calculate weights in real-space
                        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
                        {
                            dMij(Faux,i,j) =
                                dMij(Fimg_flip[iflip],i,j) *
                                dMij(fref[refnoipsi],i,j);
                        }
                        // Takes the input from Faux, and leaves the output in Maux
                        local_transformer.inverseFourierTransform();
                        CenterFFT(Maux, true);

                        // B. Calculate weights for each pixel within sigdim (Mweight)
                        my_sumweight = my_sumstoredweight = my_maxweight = 0.;
                        FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight)
                        {
                            diff = A2_plus_Xi2 - ref_scale * MAT_ELEM(Maux, i, j) * ddim2;
                            local_mindiff = XMIPP_MIN(local_mindiff, diff);
                            pdf = fracpdf * MAT_ELEM(P_phi, i, j);

                            if (!do_student)
                            {
                                // Normal distribution
                                aux = (diff - trymindiff) / sigma_noise2;
                                // next line because of numerical precision of exp-function
                                if (aux > 1000.)
                                    weight = 0.;
                                else
                                    weight = exp(-aux) * pdf;
                                // store weight
                                stored_weight = weight;
                                MAT_ELEM(Mweight, i, j) = stored_weight;
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
                                aux = (dfsigma2 + 2. * diff) / (dfsigma2 + 2. * trymindiff);
                                weight = pow(aux, df2) * pdf;
                                // Calculate extra weight acc. to Eq (10) Wang et al.
                                // Patt. Recognition Lett. 25, 701-710 (2004)
                                weight2 = ( df + ddim2 ) / ( df + (2. * diff / sigma_noise2) );
                                // Store probability weights
                                stored_weight = weight * weight2;
                                MAT_ELEM(Mweight, i, j) = stored_weight;
                                // calculate weighted sum of (X-A)^2 for sigma_noise update
                                local_wsum_corr += stored_weight * diff;
                                refw2[refno] += stored_weight;
                            }
                            // Accumulate sum weights for this (my) matrix
                            my_sumweight += weight;
                            my_sumstoredweight += stored_weight;

                            // calculated weighted sum of offsets as well
                            local_wsum_offset += weight * MAT_ELEM(Mr2, i, j);

                            if (do_norm)
                            {
                                // weighted sum of Sum_j ( X_ij*A_kj )
                                local_wsum_sc += stored_weight * (A2_plus_Xi2 - diff) / ref_scale;
                                // weighted sum of Sum_j ( A_kj*A_kj )
                                local_wsum_sc2 += stored_weight * A2[refno];
                            }
                            // keep track of optimal parameters
                            my_maxweight = XMIPP_MAX(my_maxweight, weight);

                            //FIXME: This is only for test and obtaining same result
                            //		as sequential code for comparison
                            pthread_mutex_lock(&work_mutex);
                            if (weight > maxweight)
                            {
                                maxweight = weight;
                                if (do_student)
                                	maxweight2 = weight2;
                                iopty = i;
                                ioptx = j;
                                iopt_psi = ipsi;
                                iopt_flip = iflip;
                                opt_refno = refno;
                            }
                            pthread_mutex_unlock(&work_mutex);

                            if (fast_mode && weight > maxw_ref[irefmir])
                            {
                                maxw_ref[irefmir] = weight;
                                iopty_ref[irefmir] = i;
                                ioptx_ref[irefmir] = j;
                                ioptflip_ref[irefmir] = iflip;
                            }

                        } // close for over all elements in Mweight

                        // C. only for signifcant settings, store weighted sums
                        if (my_maxweight > SIGNIFICANT_WEIGHT_LOW * maxweight )
                        {
                            sumw_refpsi[refno*nr_psi + ipsi] += my_sumstoredweight;
                            //sum_refw += my_sumweight;
                            if (iflip < nr_nomirror_flips)
                                refw[refno] += my_sumweight;
                            else
                                refw_mirror[refno] += my_sumweight;

                            // Back from smaller Mweight to original size of Maux
                            Maux.initZeros();
                            FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight)
                            {
                                MAT_ELEM(Maux, i, j) = MAT_ELEM(Mweight, i, j);
                            }
                            // Use forward FFT in convolution theorem again
                            // Takes the input from Maux and leaves it in Faux
                            local_transformer.FourierTransform();
                            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
                            {
                                dMij(mysumimgs[refnoipsi],i,j) +=
                                    conj(dMij(Faux,i,j)) *
                                    dMij(Fimg_flip[iflip],i,j);
                            }
                        }
                    } // close if Msignificant
                } // close for ipsi
            } // close for iflip
        } // close if pdf_directions

        pthread_mutex_lock(&work_mutex);

        //Update sums
        sum_refw += refw[refno] + refw_mirror[refno];
        wsum_offset += local_wsum_offset;
        wsum_corr += local_wsum_corr;
        mindiff = XMIPP_MIN(mindiff, local_mindiff);

        if (do_norm)
        {
            wsum_sc += local_wsum_sc;
            wsum_sc2 += local_wsum_sc2;
        }

        //Ask for next job
        refno = getThreadRefnoJob();

        pthread_mutex_unlock(&work_mutex);
    } // close while refno

}//close function doThreadExpectationSingleImage


void Prog_MLalign2D_prm::expectation(int iter)
{

    Matrix1D<double> dataline(DATALINELENGTH);
    Matrix2D<std::complex<double> > Fdzero(dim,hdim+1); 
    std::vector<Matrix1D<double> > docfiledata;
    int num_img_tot;

#ifdef DEBUG
    std::cerr<<"entering expectation"<<std::endl;
#endif

    rotateReference();

    // Pre-calculate pdf of all in-plane transformations
    calculatePdfInplane();

    // Initialize weighted sums
    LL = 0.;
    wsumimgs.clear();
    sumw.clear();
    sumw2.clear();
    sumwsc.clear();
    sumwsc2.clear();
    sumw_mirror.clear();
    wsum_sigma_noise = 0.;
    wsum_sigma_offset = 0.;
    sumfracweight = 0.;
    dataline.initZeros();

    for (int i = 0; i < SF.ImgNo(); i++)
        docfiledata.push_back(dataline);

    Fdzero.initZeros();
    for (int i = 0; i < n_ref * nr_psi; i++)
    {
        wsumimgs.push_back(Fdzero);
    }

    for (int refno = 0; refno < n_ref; refno++)
    {
        sumw.push_back(0.);
        sumw2.push_back(0.);
        sumwsc.push_back(0.);
        sumw_mirror.push_back(0.);
        sumwsc2.push_back(0.);
    }

    // Local variables of old threadExpectationSingleImage
    ImageXmipp img;
    FileName fn_img, fn_trans;    
    Matrix1D<double> opt_offsets(2);
    float old_phi = -999., old_theta = -999.;
    double opt_flip;
    
    //Some initializations
    opt_scale = 1., bgmean = 0.;
    Msignificant.resize(n_ref, nr_psi * nr_flip);
    
    // Calculate myFirst and myLast image for this thread
    int nn = SF.ImgNo();
    int myNum = nn;

    if (verb > 0)
        init_progress_bar(myNum);

    // Loop over all images
    int c = XMIPP_MAX(1, myNum / 60);
    
    for (int imgno = 0; imgno < nn; imgno++)
    {
        SF.go_beginning();
        SF.jump(imgno, SelLine::ACTIVE);
        fn_img = SF.get_current_file();
        
        img.read(fn_img, false, false, false, false);
        img().setXmippOrigin();
        Xi2 = img().sum2();
        Mimg=img();
                    
        // These two parameters speed up expectationSingleImage
        opt_refno = imgs_optrefno[imgno];
        trymindiff = imgs_trymindiff[imgno];
        if (trymindiff < 0.)
            // 90% of Xi2 may be a good idea (factor half because 0.5*diff is calculated)
            trymindiff = trymindiff_factor * 0.5 * Xi2; 
        
        if (do_norm)
        {
            bgmean = imgs_bgmean[imgno];
            opt_scale = imgs_scale[imgno];
        }
        
        // Get optimal offsets for all references
        if (fast_mode)
        {
            allref_offsets = imgs_offsets[imgno];
        }
        
        // Read optimal orientations from memory
        if (limit_rot)
        {
            old_phi = imgs_oldphi[imgno];
            old_theta = imgs_oldtheta[imgno];
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
            imgs_offsets[imgno] = allref_offsets;
        }
        
        // Store mindiff for next iteration
        imgs_trymindiff[imgno] = trymindiff;
        // Store opt_refno for next iteration
        imgs_optrefno[imgno] = opt_refno;
        
        // Store optimal phi and theta in memory
        if (limit_rot)
        {
            imgs_oldphi[imgno] = Iref[opt_refno].Phi();
            imgs_oldtheta[imgno] = Iref[opt_refno].Theta();
        }
        
        // Store optimal normalization parameters in memory
        if (do_norm)
        {
            imgs_scale[imgno] = opt_scale;
            imgs_bgmean[imgno]  = bgmean;
        }
        
        // Output docfile
        if (-opt_psi > 360.)
        {
            opt_psi += 360.;
            opt_flip = 1.;
        }
        else
        {
            opt_flip = 0.;
        }

        docfiledata[imgno](0) = Iref[opt_refno].Phi();     // rot
        docfiledata[imgno](1) = Iref[opt_refno].Theta();   // tilt
        docfiledata[imgno](2) = opt_psi + 360.;            // psi
        docfiledata[imgno](3) = opt_offsets(0);            // Xoff
        docfiledata[imgno](4) = opt_offsets(1);            // Yoff
        docfiledata[imgno](5) = (double)(opt_refno + 1);   // Ref
        docfiledata[imgno](6) = opt_flip;                  // Mirror
        docfiledata[imgno](7) = fracweight;                // P_max/P_tot
        docfiledata[imgno](8) = dLL;                       // log-likelihood
        
        if (do_norm)
        {
            docfiledata[imgno](9)  = bgmean;               // background mean
            docfiledata[imgno](10) = opt_scale;            // image scale
        }
        
        if (do_student)
        {
            docfiledata[imgno](11) = maxweight2;           // Robustness weight
        }
        
        if (verb > 0 && imgno % c == 0)
            progress_bar(imgno);
    }
    
    if (verb > 0)
        progress_bar(myNum);
    
    // Rotate back and calculate weighted sums
    reverseRotateReference();

    // Send back output in the form of a DocFile
    SF.go_beginning();
    for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
    {
        DFo.append_comment(SF.NextImg());
        DFo.append_data_line(docfiledata[imgno]);
    }

#ifdef DEBUG
    std::cerr<<"leaving expectation"<<std::endl;
#endif

}

// Update all model parameters
void Prog_MLalign2D_prm::maximization(int refs_per_class)
{

#ifdef DEBUG
    std::cerr<<"entering maximization"<<std::endl;
#endif

    Matrix1D<double> rmean_sigma2, rmean_signal2;
    Matrix1D<int> center(2), radial_count;
    Matrix2D<std::complex<double> > Faux, Faux2;
    Matrix2D<double> Maux;
    FileName fn_tmp;
    double rr, thresh, aux, sumw_allrefs2 = 0.;
    int c;

    // Update the reference images
    sumw_allrefs = 0.;
    for (int refno=0; refno<n_ref; refno++)
    {
        if (!do_student && sumw[refno] > 0.)
        {
            Iref[refno]() = wsum_Mref[refno];
            Iref[refno]() /= sumwsc2[refno];
            Iref[refno].set_weight(sumw[refno]);
            sumw_allrefs += sumw[refno];
        }
 	else if (do_student && sumw2[refno] > 0.)
        {
	    Iref[refno]() = wsum_Mref[refno];
	    Iref[refno]() /= sumwsc2[refno];
	    Iref[refno].set_weight(sumw2[refno]);
	    sumw_allrefs += sumw[refno];
	    sumw_allrefs2 += sumw2[refno];
 	}
        else
        {
            Iref[refno].set_weight(0.);
            Iref[refno]().initZeros(dim, dim);
        }
    }

    // Adjust average scale (nr_classes will be smaller than n_ref for the 3D case!)
    if (do_norm) {
        int iclass, nr_classes = ROUND(n_ref / refs_per_class);
        std::vector<double> wsum_scale(nr_classes), sumw_scale(nr_classes);
        ldiv_t temp;
        average_scale = 0.;
        for (int refno=0;refno<n_ref; refno++)
        {
            average_scale += sumwsc[refno];
            temp = ldiv( refno, refs_per_class );
            iclass = ROUND(temp.quot);
            wsum_scale[iclass] += sumwsc[refno];
            sumw_scale[iclass] += sumw[refno];
        }
        for (int refno=0;refno<n_ref; refno++)
        {
            temp = ldiv( refno, refs_per_class );
            iclass = ROUND(temp.quot);
            if (sumw_scale[iclass]>0.)
            {
                refs_avgscale[refno] = wsum_scale[iclass]/sumw_scale[iclass];
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
        for (int refno=0;refno<n_ref; refno++)
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
        // Dont allow sigma_offset to drop below one pixel
        sigma_offset = XMIPP_MAX(1. , sigma_offset);

        // Update sigdim, i.e. the number of pixels that will be considered in the translations
        if (save_mem2)
            sigdim = 2 * CEIL(sigma_offset * 3);
        else
            sigdim = 2 * CEIL(sigma_offset * 6);
        sigdim++; // (to get uneven number)
        sigdim = XMIPP_MIN(dim, sigdim);
    }

    // Update the noise parameters
    if (!fix_sigma_noise)
    {
//      The following converges faster according to McLachlan&Peel (2000)
//	Finite Mixture Models, Wiley p. 228!
	if (do_student && do_student_sigma_trick)
	    sigma_noise = sqrt(wsum_sigma_noise / (sumw_allrefs2 * dim * dim));
	else
	    sigma_noise = sqrt(wsum_sigma_noise / (sumw_allrefs * dim * dim));

        sigma_noise2 = sigma_noise * sigma_noise;
        if (do_student)
            dfsigma2 = df * sigma_noise2;
    }

#ifdef DEBUG
    std::cerr<<"leaving maximization"<<std::endl;
#endif
}

// Check convergence
bool Prog_MLalign2D_prm::checkConvergence()
{

#ifdef DEBUG
    std::cerr<<"entering checkConvergence"<<std::endl;
#endif

    bool converged = true;
    double convv;
    Matrix2D<double> Maux;

    Maux.resize(dim, dim);
    Maux.setXmippOrigin();

    conv.clear();
    for (int refno=0;refno<n_ref; refno++)
    {
        if (Iref[refno].weight() > 0.)
        {
            multiplyElements(Iold[refno](), Iold[refno](), Maux);
            convv = 1. / (Maux.computeAvg());
            Maux = Iold[refno]() - Iref[refno]();
            multiplyElements(Maux, Maux, Maux);
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
    std::cerr<<"leaving checkConvergence"<<std::endl;
#endif
    return converged;
}

void Prog_MLalign2D_prm::writeOutputFiles(const int iter)
{

    FileName          fn_tmp, fn_base, fn_tmp2;
    ImageXmipp        Itmp;
    Matrix1D<double>  fracline(3);
    SelFile           SFo, SFc;
    DocFile           DFl;
    std::string       comment;
    std::ofstream     fh;

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
    for (int refno=0;refno<n_ref; refno++)
    {
        fn_tmp = fn_base + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "");
        fn_tmp = fn_tmp + ".xmp";
        Itmp=Iref[refno];
        Itmp.write(fn_tmp);
        SFo.insert(fn_tmp, SelLine::ACTIVE);
        fracline(0) = alpha_k[refno];
        fracline(1) = mirror_fraction[refno];
        fracline(2) = 1000 * conv[refno]; // Output 1000x the change for precision
        if (do_norm) fracline(3) = refs_avgscale[refno];
        DFl.insert_comment(fn_tmp);
        DFl.insert_data_line(fracline);
    }

    // Write out sel & log-file
    fn_tmp = fn_base + ".sel";
    SFo.write(fn_tmp);

    DFl.go_beginning();
    comment = "MLalign2D-logfile: Number of images= " + floatToString(sumw_allrefs);
    comment += " LL= " + floatToString(LL, 15, 10) + " <Pmax/sumP>= " + floatToString(sumfracweight, 10, 5);
    if (do_norm)
        comment+= " <scale>= " + floatToString(average_scale, 10, 5);
    DFl.insert_comment(comment);
    comment = "-noise " + floatToString(sigma_noise, 15, 12) + " -offset " + floatToString(sigma_offset, 15, 12) + " -istart " + integerToString(iter + 1) + " -doc " + fn_base + ".doc";
    DFl.insert_comment(comment);
    DFl.insert_comment(cline);
    if (do_norm) 
        DFl.insert_comment("columns: model fraction (1); mirror fraction (2); 1000x signal change (3); avg scale (4)");
    else
        DFl.insert_comment("columns: model fraction (1); mirror fraction (2); 1000x signal change (3)");
    fn_tmp = fn_base + ".log";
    DFl.write(fn_tmp);

    // Write out docfile with optimal transformation & references
    fn_tmp = fn_base + ".doc";
    DFo.write(fn_tmp);

    // Also write out selfiles of all experimental images,
    // classified according to optimal reference image
    for (int refno = 0;refno < n_ref; refno++)
    {
        DFo.go_beginning();
        SFo.clear();
        for (int n = 0; n < DFo.dataLineNo(); n++)
        {
            DFo.next();
            fn_tmp = ((DFo.get_current_line()).get_text()).erase(0, 3);
            DFo.adjust_to_data_line();
            if ((refno + 1) == (int)DFo(5)) SFo.insert(fn_tmp, SelLine::ACTIVE);
        }
        fn_tmp = fn_root + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "sel");
        SFo.write(fn_tmp);
    }

}

