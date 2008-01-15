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
#include "ml_align2d.h"

// Read arguments ==========================================================
void Prog_MLalign2D_prm::read(int argc, char **argv, bool ML3D)
{

    // Generate new command line for restart procedure
    cline = "";
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
            generateCommandLine(comment, argc, argv, copy);
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
        for (int i = 1; i < argc; i++)
        {
            cline = cline + (std::string)argv[i] + " ";
        }
    }

    // Read command line
    if (checkParameter(argc, argv, "-more_options"))
    {
	usage();
	extended_usage();
    }
    n_ref = textToInteger(getParameter(argc, argv, "-nref", "0"));
    fn_ref = getParameter(argc, argv, "-ref", "");
    fn_sel = getParameter(argc, argv, "-i");
    fn_root = getParameter(argc, argv, "-o", "ml2d");
    psi_step = textToFloat(getParameter(argc, argv, "-psi_step", "5"));
    Niter = textToInteger(getParameter(argc, argv, "-iter", "100"));
    istart = textToInteger(getParameter(argc, argv, "-istart", "1"));
    sigma_noise = textToFloat(getParameter(argc, argv, "-noise", "1"));
    sigma_offset = textToFloat(getParameter(argc, argv, "-offset", "3"));
    do_mirror = checkParameter(argc, argv, "-mirror");
    eps = textToFloat(getParameter(argc, argv, "-eps", "5e-5"));
    fn_frac = getParameter(argc, argv, "-frac", "");
    write_docfile = !checkParameter(argc, argv, "-dont_output_docfile");
    write_selfiles = !checkParameter(argc, argv, "-dont_output_selfiles");
    write_intermediate = !checkParameter(argc, argv, "-dont_output_intermediate");
    fix_fractions = checkParameter(argc, argv, "-fix_fractions");
    fix_sigma_offset = checkParameter(argc, argv, "-fix_sigma_offset");
    fix_sigma_noise = checkParameter(argc, argv, "-fix_sigma_noise");
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    maxCC_rather_than_ML = checkParameter(argc, argv, "-maxCC");
    fast_mode = checkParameter(argc, argv, "-fast");
    C_fast = textToFloat(getParameter(argc, argv, "-C", "1e-12"));
    max_shift = textToFloat(getParameter(argc, argv, "-max_shift", "-1"));
    save_mem1 = checkParameter(argc, argv, "-save_memA");
    save_mem2 = checkParameter(argc, argv, "-save_memB");
    save_mem3 = checkParameter(argc, argv, "-save_memC");
    search_shift = textToFloat(getParameter(argc, argv, "-search_shift", "999."));
    fn_doc = getParameter(argc, argv, "-doc", "");
    do_ML3D = ML3D;

    // Hidden arguments
    do_esthetics = checkParameter(argc, argv, "-esthetics");
    anneal = textToFloat(getParameter(argc, argv, "-anneal", "1"));
    anneal_step = textToFloat(getParameter(argc, argv, "-anneal_step", "1"));
    do_write_offsets = checkParameter(argc, argv, "-write_offsets");
    fn_scratch = getParameter(argc, argv, "-scratch", "");
    zero_offsets = checkParameter(argc, argv, "-zero_offsets");
    debug = textToInteger(getParameter(argc, argv, "-debug","0"));
    do_student = checkParameter(argc, argv, "-robust");
    df = (double) textToInteger(getParameter(argc, argv, "-df", "3"));

    //only for interaction with Refine3D:
    search_rot = textToFloat(getParameter(argc, argv, "-search_rot", "999."));

    // For improved control of MPI jobs
    fn_control = getParameter(argc, argv, "-control", "");
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
        if (maxCC_rather_than_ML)
        {
            std::cerr << "  -> Use a maxCC instead of a maximum likelihood target." << std::endl;
        }
        if (fast_mode)
        {
            std::cerr << "  -> Use fast, reduced search-space approach with C = " << C_fast << std::endl;
            if (zero_offsets)
                std::cerr << "    + Start from all-zero translations" << std::endl;
        }
        if (search_shift < 999.)
            std::cerr << "    + Limit translational search to +/- " << search_shift << " pixels" << std::endl;
        if (search_rot < 180.)
            std::cerr << "    + Limit orientational search to +/- " << search_rot << " degrees" << std::endl;
        if (save_mem1)
            std::cerr << "  -> Save_memory A: recalculate real-space rotations in -fast" << std::endl;
        if (save_mem2)
            std::cerr << "  -> Save_memory B: limit translations to 3 sigma_offset " << std::endl;
        if (save_mem3)
            std::cerr << "  -> Save_memory C: do not store rotated references; rotate experimental image instead " << std::endl;


        // Hidden stuff
        if (!write_intermediate)
        {
            std::cerr << "  -> Do not write out images after each iteration." << std::endl;
        }
        if (fix_fractions && !maxCC_rather_than_ML)
        {
            std::cerr << "  -> Do not update estimates of model fractions." << std::endl;
        }
        if (fix_sigma_offset && !maxCC_rather_than_ML)
        {
            std::cerr << "  -> Do not update sigma-estimate of origin offsets." << std::endl;
        }
        if (fix_sigma_noise && !maxCC_rather_than_ML)
        {
            std::cerr << "  -> Do not update sigma-estimate of noise." << std::endl;
        }
        if (do_esthetics)
        {
            std::cerr << "  -> Perform esthetics on (0,0)-pixel artifacts" << std::endl;
        }
        if (do_student)
        {
            std::cerr << "  -> Use t-student distribution for improved robustness " << endl;
        }
        std::cerr << " -----------------------------------------------------------------" << endl;

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
    std::cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_MLalign2D_prm::extended_usage(bool ML3D)
{
    std::cerr << "Additional options: " << std::endl;
    std::cerr << " [ -eps <float=5e-5> ]         : Stopping criterium \n";
    std::cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
    std::cerr << " [ -psi_step <float=5> ]       : In-plane rotation sampling interval [deg]\n";
    std::cerr << " [ -noise <float=1> ]          : Expected standard deviation for pixel noise \n";
    std::cerr << " [ -offset <float=3> ]         : Expected standard deviation for origin offset [pix]\n";
    std::cerr << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n";
    std::cerr << " [ -C <double=1e-12> ]         : Significance criterion for fast approach \n";
    if (!ML3D) std::cerr << " [ -restart <logfile> ]        : restart a run with all parameters as in the logfile \n";
    if (!ML3D) std::cerr << " [ -istart <int> ]             : number of initial iteration \n";
    std::cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise \n";
    std::cerr << " [ -fix_sigma_offset]          : Do not re-estimate the standard deviation in the origin offsets \n";
    std::cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
    std::cerr << " [ -dont_output_docfile ]      : Do not write out docfile with most likely angles & translations \n";
    std::cerr << " [ -dont_output_selfiles ]     : Do not write out selfiles with most likely class assignments \n";
    std::cerr << " [ -maxCC ]                    : Use maximum cross-correlation instead of maximum likelihood target \n";
    std::cerr << " [ -search_shift <float=999>]  : Limit of translational search [pix] (does NOT use FFT) \n";
    std::cerr << " [ -max_shift <float=dim/4>]   : Dont trust shifts larger than max_shift \n";
    std::cerr << " [ -doc <docfile=\"\"> ]         : Read initial angles and offsets from docfile \n";
    std::cerr << " [ -write_offsets ]            : Save memory by writing optimal offsets to disc (disc-access intensive) \n";
    std::cerr << std::endl;
    exit(1);
}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produce_Side_info2)
void Prog_MLalign2D_prm::produce_Side_info()
{

    FileName                    fn_img, fn_tmp, fn_base, fn_tmp2;
    ImageXmipp                  img;
    FourierImageXmipp           fourimg;
    SelLine                     SL;
    SelFile                     SFtmp, SFpart;
    Matrix1D<double>            offsets(2), dum;
    Matrix2D<double>            A(3, 3), Maux, Maux2;
    Matrix2D<std::complex<double> >  Faux;
    Matrix1D<int>               center(2), radial_count;
    std::vector<int>                 tmppointp, tmppointp_nolow, tmppointi, tmppointj;
    bool                        uniqname, nomoredirs;
    float                       xx, yy;
    double                      av, psi, aux, Q0;
    int                         im, jm;

    // Read selfile with experimental images
    SF.read(fn_sel);

    // Get image sizes and total number of images
    SF.ImgSize(dim, dim);
    hdim = dim / 2;
    dim2 = dim * dim;
    nr_exp_images = SF.ImgNo();
    if (do_student)
    {
	df2 = - ( df + dim2 ) / 2. ;
    }

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

    // Check deterministic annealing parameters
    if (anneal < 1.) REPORT_ERROR(1, "Initial annealing parameter should be larger than 1 (recommended 2-5)");
    if (anneal_step < 0.) REPORT_ERROR(1, "Anneal_step should be positive (recommended 0.5-1)");

    // Check the uniqueness of all filenames (for names of temporary offsets files)
    uniqname = false;
    nomoredirs = false;
    offsets_keepdir = 0;
    while ((uniqname == false) && (nomoredirs == false))
    {
        SF.go_beginning();
        SFtmp.clear();
        while (!SF.eof())
        {
            fn_tmp = SF.NextImg();
            fn_tmp2 = fn_tmp.remove_directories(offsets_keepdir);
            if (fn_tmp == fn_tmp2) nomoredirs = true;
            SFtmp.insert(fn_tmp2);
        }
        SFtmp.sort_by_filenames();
        SFtmp.go_beginning();
        uniqname = true;
        while (!SFtmp.eof())
        {
            fn_tmp = SFtmp.NextImg();
            fn_tmp2 = SFtmp.NextImg();
            if (fn_tmp == fn_tmp2)
            {
                uniqname = false;
                offsets_keepdir++;
                break;
            }
        }
    }
    SFtmp.clear();
    if (!uniqname)
        REPORT_ERROR(1, "Prog_MLalign2D_prm: Provide a selfile with unique image names (preferably all in one directory)");

    // Set nr_psi & nr_flip and construct flipping matrices
    if (save_mem3)
    {
        nr_psi = 1;
        nr_flip = nr_nomirror_flips = CEIL(360. / psi_step);
        psi_step = 360. / nr_psi;
        // store all rotation (and mirror) matrices
        FOR_ALL_FLIPS()
        {
            double ang = (double)(iflip * 360. / nr_flip) + SMALLANGLE;
            A = rotation2DMatrix(ang);
            F.push_back(A);
        }
        if (do_mirror)
        {
            FOR_ALL_FLIPS()
            {
                double ang = (double)(iflip * 360. / nr_flip);
                A = rotation2DMatrix(ang);
                A(0, 0) *= -1.;
                A(0, 1) *= -1.;
                F.push_back(A);
            }
            nr_flip *= 2;
        }
    }
    else
    {
        psi_max = 90.;
        nr_psi = CEIL(psi_max / psi_step);
        psi_step = psi_max / nr_psi;
        nr_flip = nr_nomirror_flips = 4;
        // 0, 90, 180 & 270 degree flipping, as well as mirror
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
    }

    // Set some stuff for maxCC-mode
    if (maxCC_rather_than_ML)
    {
        fix_sigma_noise = true;
        fix_sigma_offset = true;
        fix_fractions = true;
    }

    // Set & check max_shift
    if (max_shift < 0) max_shift = dim / 4.;

    // Set limit_rot & limit_trans
    if (search_rot < 180.) limit_rot = true;
    else limit_rot = false;
    if (search_shift < 999) limit_trans = true;
    else limit_trans = false;

    // Fill limited translation shift-vectors
    if (limit_trans)
    {
        nr_trans = 0;
        Maux.resize(dim, dim);
        Maux.setXmippOrigin();
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Maux)
        {
            double rr = (double)i * i + (double)j * j;
            if (rr <= search_shift*search_shift)
            {
                offsets(0) = (double)j;
                offsets(1) = (double)i;
                Vtrans.push_back(offsets);
                if (i == 0 && j == 0) zero_trans = nr_trans;
                nr_trans++;
            }
        }
    }


}

// Generate initial references =============================================
void Prog_MLalign2D_prm::generate_initial_references()
{

    SelFile SFtmp, SFout;
    ImageXmipp Iave, Itmp;
    double dummy;
    FileName fn_tmp;
    SelLine line;

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
// same for all processors! (in contrast to produce_Side_info)
void Prog_MLalign2D_prm::produce_Side_info2(int nr_vols)
{

    int                       c, idum, refno = 0;
    DocFile                   DF, DF2;
    DocLine                   DL;
    double                    offx, offy, aux, sumfrac = 0.;
    FileName                  fn_tmp;
    ImageXmipp                img;
    std::vector<double>            Vdum;

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
        img.read(SFr.NextImg(), false, false, true, false);
        img().setXmippOrigin();
        Iref.push_back(img);
        if (!save_mem3) Iold.push_back(img);
        // Default start is all equal model fractions
        alpha_k.push_back((double)1 / SFr.ImgNo());
        Iref[refno].weight() = alpha_k[refno] * (double)nr_exp_images;
        // Default start is half-half mirrored images
        if (do_mirror) mirror_fraction.push_back(0.5);
        else mirror_fraction.push_back(0.);
        n_ref++;
        refno++;
    }

    // Make scratch directory for the temporary origin offsets
    if (fn_scratch!="")
    {
        fn_scratch += "/ml_align2d_offsets";
	// Clean it if already existing
    	system(((std::string)"rm -rf "+fn_scratch).c_str());
 	// Generate new one
    	system(((std::string)"mkdir -p " + fn_scratch).c_str());
    }

    // Read optimal origin offsets from fn_doc
    if (fn_doc != "")
    {
        DF.read(fn_doc);
        SF.go_beginning();
        while (!SF.eof())
        {
            fn_tmp = SF.NextImg();
            if (DF.search_comment(fn_tmp))
            {
                imgs_oldxoff.push_back(DF(3));
                imgs_oldyoff.push_back(DF(4));
            }
            else
            {
                REPORT_ERROR(1, (std::string)"Prog_MLalign2D_prm: Cannot find " + fn_tmp + " in docfile " + fn_doc);
            }
        }
        DF.go_beginning();
    }

    // If we don't write offsets to disc, fill imgs_offsets vectors
    if (!do_write_offsets)
    {
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
    }

    // For limited orientational search: fill imgs_oldphi & imgs_oldtheta
    // (either read from fn_doc or initialize to -999.)
    if (limit_rot)
    {
        imgs_oldphi.clear();
        imgs_oldtheta.clear();
        SF.go_beginning();
        while (!SF.eof())
        {
            fn_tmp = SF.NextImg();
            if (fn_doc != "")
            {
                if (DF.search_comment(fn_tmp))
                {
                    imgs_oldphi.push_back(DF(0));
                    imgs_oldtheta.push_back(DF(1));
                }
                else
                {
                    REPORT_ERROR(1, (std::string)"Prog_MLalign2D_prm: Cannot find " + fn_tmp + " in docfile " + fn_doc);
                }
            }
            else
            {
                imgs_oldphi.push_back(-999.);
                imgs_oldtheta.push_back(-999.);
            }
        }
    }
    DF.clear();

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

void Prog_MLalign2D_prm::write_offsets(FileName fn, std::vector<double> &data)
{

    std::ofstream fh;
    int itot;

    if (fn_scratch!="")
    {
    	// Write to scratch disc
    	fn = fn_scratch + "/" + fn;
    }

    fh.open((fn).c_str(), std::ios::out);
    if (!fh)
    {
        fh.clear();
        // Create the directory if it does not exist yet, and try again
        std::string dirname;
        int last_slash = ((std::string)fn).rfind("/");
        dirname = ((std::string)fn).erase(last_slash);
        if (!exists(dirname)) system(((std::string)"mkdir -p " + dirname).c_str());
        fh.open((fn).c_str(), std::ios::out);
        if (!fh)
            REPORT_ERROR(3008, (std::string)"Prog_MLalign2D_prm: Cannot write file: " + fn);
    }

    itot = data.size();
    fh << itot << "\n";
    for (int i = 0; i < itot; i += 2)
    {
        fh << data[i] << " " << data[i+1] << "\n";
    }
    fh.close();
    data.clear();

}

bool Prog_MLalign2D_prm::read_offsets(FileName fn, std::vector<double> &data)
{

    std::ifstream fh;
    int ii, itot, nr_off, itoth, nr_offh;
    double remain;
    std::vector<double> data1;

    if (fn_scratch!="")
    {
    	// Read from scratch disc
    	fn = fn_scratch + "/" + fn;
    }

    if (!exists(fn)) return false;
    else
    {
        fh.open((fn).c_str(), std::ios::in);
        if (!fh) return false;
        else
        {
            fh >> itot;
            if (do_mirror) nr_off = n_ref * 4;
            else nr_off = n_ref * 2;
            if (itot != nr_off)
            {
                fh.close();
                return false;
            }
            else
            {
                data.clear();
                data.resize(itot);
                for (int i = 0; i < itot; i += 2)
                {
                    fh >> data[i];
                    fh >> data[i+1];
                }
                fh.close();
                return true;
            }
        }
    }

}


// Calculate probability density function of all in-plane transformations phi
void Prog_MLalign2D_prm::calculate_pdf_phi()
{

    double r2, pdfpix, sum;
    P_phi.resize(dim, dim);
    P_phi.setXmippOrigin();
    Mr2.resize(dim, dim);
    Mr2.setXmippOrigin();

    sum=0.;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(P_phi)
    {
        r2 = (double)(j * j + i * i);
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
void Prog_MLalign2D_prm::rotate_reference(std::vector< ImageXmippT<double> > &Iref,
        bool fill_real_space,
        bool fill_fourier_space,
        std::vector <std::vector< Matrix2D<double> > > &Mref,
        std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref)
{

    double AA, stdAA, psi, dum, avg, mean_ref, stddev_ref, dummy;
    Matrix2D<double> Maux;
    Matrix2D<std::complex<double> > Faux;
    std::vector<Matrix2D<std::complex <double> > > dumF;
    std::vector<Matrix2D<double> > dumM;
    Matrix2D<int> mask, omask;
    Matrix2D<double> cmask;

    Maux.initZeros(dim, dim);
    Maux.setXmippOrigin();
    Fref.clear();
    Mref.clear();
    A2.clear();

    // prepare masks
    mask.resize(dim, dim);
    mask.setXmippOrigin();
    BinaryCircularMask(mask, hdim, INNER_MASK);
    omask.resize(dim, dim);
    omask.setXmippOrigin();
    BinaryCircularMask(omask, hdim, OUTSIDE_MASK);

    FOR_ALL_MODELS()
    {
        Mref.push_back(dumM);
        Fref.push_back(dumF);
        computeStats_within_binary_mask(omask, Iref[refno](), dum, dum, avg, dum);
        FOR_ALL_ROTATIONS()
        {
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
                A2.push_back(AA);
            }
            // Subtract mean_ref from image prior to FFT for maxCC
            if (maxCC_rather_than_ML)
            {
                Maux.computeStats(mean_ref, stddev_ref, dummy, dummy);
                Maux -= mean_ref;
                if (ipsi == 0) A2[refno] = stddev_ref;
            }
            if (fill_real_space)
            {
                Mref[refno].push_back(Maux);
                if (AA > 0) Mref[refno][ipsi] *= sqrt(stdAA / AA);
            }
            if (fill_fourier_space)
            {
                FourierTransformHalf(Maux, Faux);
		Faux *= dim * dim;
		FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
                {
		    dMij(Faux, i, j) = conj(dMij(Faux, i, j));
		}
                Fref[refno].push_back(Faux);
                if (AA > 0) Fref[refno][ipsi] *= sqrt(stdAA / AA);
            }
        }
        // If we dont use save_mem1 Iref[refno] is useless from here on
        if (!save_mem1) Iref[refno]().resize(0, 0);
    }

}

// Collect all rotations and sum to update Iref() for all models ==========
void Prog_MLalign2D_prm::reverse_rotate_reference(
    std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref,
    std::vector <std::vector< Matrix2D<double> > > &Mref, bool real_space,
    std::vector<Matrix2D<double > > &Mnew)
{

    double psi, dum, avg, ang;
    Matrix2D<double> Maux, Maux2;
    Matrix2D<std::complex<double> > Faux;
    Matrix2D<int> mask, omask;
    Maux.resize(dim, dim);
    Maux2.resize(dim, dim);
    Maux.setXmippOrigin();
    Maux2.setXmippOrigin();

    Mnew.clear();
    mask.resize(dim, dim);
    mask.setXmippOrigin();
    BinaryCircularMask(mask, hdim, INNER_MASK);
    omask.resize(dim, dim);
    omask.setXmippOrigin();
    BinaryCircularMask(omask, hdim, OUTSIDE_MASK);

    FOR_ALL_MODELS()
    {
        Maux.initZeros();
        Mnew.push_back(Maux);
        FOR_ALL_ROTATIONS()
        {
            // Add arbitrary number to avoid 0-degree rotation without interpolation effects
            psi = (double)(ipsi * psi_max / nr_psi) + SMALLANGLE;
            if (real_space)
            {
                Maux = Mref[refno][ipsi];
            }
            else
            {
                InverseFourierTransformHalf(Fref[refno][ipsi], Maux, dim);
		Maux /= dim * dim;
		CenterFFT(Maux, true);
            }
            computeStats_within_binary_mask(omask, Maux, dum, dum, avg, dum);
            Maux.rotateBSpline(3, -psi, Maux2, WRAP);
            apply_binary_mask(mask, Maux2, Maux2, avg);
            Mnew[refno] += Maux2;
        }

        // perform correction of origin-pixel artifacts
        if (do_esthetics)
            MAT_ELEM(Mnew[refno], 0, 0) = (MAT_ELEM(Mnew[refno], 1, 0) + MAT_ELEM(Mnew[refno], 0, 1) +
                                           MAT_ELEM(Mnew[refno], -1, 0) + MAT_ELEM(Mnew[refno], 0, -1)) / 4;

    }

}

void Prog_MLalign2D_prm::preselect_directions(float &phi, float &theta,
        std::vector<double> &pdf_directions)
{

    float phi_ref, theta_ref, angle, angle2;
    Matrix1D<double> u, v;

    pdf_directions.clear();
    pdf_directions.resize(n_ref);
    FOR_ALL_MODELS()
    {
        if (!limit_rot || (phi == -999. && theta == -999.)) pdf_directions[refno] = 1.;
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
void Prog_MLalign2D_prm::preselect_significant_model_phi(
    Matrix2D<double> &Mimg, std::vector<double > &offsets,
    std::vector <std::vector< Matrix2D<double > > > &Mref,
    Matrix2D<int> &Msignificant,
    std::vector<double> &pdf_directions)
{


    Matrix2D<double> Maux, Maux2, Mrefl, Mdsig(n_ref, nr_psi*nr_flip);
    double ropt, sigma_noise2, aux, weight, diff, pdf, fracpdf;
    double Xi2, A2_plus_Xi2, dfsigma2, CC;
    double mindiff = 99.e99;
    double maxweight = -99.e99;
    int irot, irefmir;
    std::vector<double> maxw_ref(2*n_ref);
    Matrix1D<double> trans(2);
    std::vector<Matrix2D<double> > Mrot;

    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Maux2.resize(dim, dim);
    Maux2.setXmippOrigin();
    sigma_noise2 = sigma_noise * sigma_noise;
    Xi2 = Mimg.sum2();
    Msignificant.initZeros();
    dfsigma2 = df * sigma_noise2;

    // Flip images and calculate correlations and maximum correlation
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            A2_plus_Xi2 = 0.5 * (A2[refno] + Xi2);
            if (save_mem1)
            {
                Mrot.clear();
                FOR_ALL_ROTATIONS()
                {
                    double psi = (double)(ipsi * psi_max / nr_psi) + SMALLANGLE;
                    Iref[refno]().rotateBSpline(3, psi, Maux, WRAP);
                    Mrot.push_back(Maux);
                }
            }
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                // Do not trust optimal offsets if they are larger than 3*sigma_offset:
                ropt = sqrt(offsets[2*irefmir] * offsets[2*irefmir] + offsets[2*irefmir+1] * offsets[2*irefmir+1]);
                if (ropt > 3*sigma_offset)
                {
                    FOR_ALL_ROTATIONS()
                    {
                        irot = iflip * nr_psi + ipsi;
                        dMij(Msignificant, refno, irot) = 1;
                    }
                }
                else
                {
                    trans(0) = offsets[2*irefmir];
                    trans(1) = offsets[2*irefmir+1];
                    Mimg.translate(trans, Maux, true);
                    applyGeometry(Maux2, F[iflip], Maux, IS_INV, WRAP);
                    FOR_ALL_ROTATIONS()
                    {
                        irot = iflip * nr_psi + ipsi;
                        dMij(Msignificant, refno, irot) = 0;
                        CC = A2_plus_Xi2;
                        if (save_mem1) Mrefl = Mrot[ipsi];
                        else Mrefl = Mref[refno][ipsi];
                        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Maux2)
                        {
                            CC -= dMij(Maux2, i, j) * dMij(Mrefl, i, j);
                        }
                        dMij(Mdsig, refno, irot) = CC;
                        if (CC < mindiff) mindiff = CC;
                    }
                }
            }
        }
    }

    // Now that we have mindiff calculate the weighting matrices and maxweight
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    if (!dMij(Msignificant, refno, irot))
                    {
                        if (iflip < nr_nomirror_flips) fracpdf = alpha_k[refno] * (1. - mirror_fraction[refno]);
                        else fracpdf = alpha_k[refno] * mirror_fraction[refno];
			diff = dMij(Mdsig, refno, irot);
			pdf = fracpdf * MAT_ELEM(P_phi, (int)offsets[2*irefmir+1], (int)offsets[2*irefmir]);
			if (!do_student)
			{
			    // normal distribution
			    aux = (diff - mindiff) / sigma_noise2;
			    // next line because of numerical precision of exp-function
			    if (aux > 1000.) weight = 0.;
			    else weight = exp(-aux) * pdf;
			}
			else
			{
			    // t-student distribution
			    aux = (dfsigma2 + 2. * diff) / (dfsigma2 + 2. * mindiff);
			    weight = pow(aux, df2) * pdf;
			}
                        dMij(Mdsig, refno, irot) = weight;
                        if (weight > maxw_ref[irefmir]) maxw_ref[irefmir] = weight;
                    }
                }
            }
        }
    }

    // Now that we have maxweight calculate which weighting matrices are significant
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    if (!dMij(Msignificant, refno, irot))
                    {
                        if (dMij(Mdsig, refno, irot) >= C_fast*maxw_ref[irefmir]) dMij(Msignificant, refno, irot) = 1;
                        else dMij(Msignificant, refno, irot) = 0;
                    }
                }
            }
        }
    }

}

// Calculate translated matrices for all limited translations
// for each of the flipped variants
void Prog_MLalign2D_prm::calculate_realspace_offsets(
    Matrix2D<double> &Mimg, std::vector<double > &offsets,
    std::vector<double > &pdf_directions,
    std::vector<std::vector<Matrix2D<double> > > &Mimg_trans,
    Matrix2D<int> &Moffsets, Matrix2D<int> &Moffsets_mirror)
{

    int irefmir, ix, iy, opt_ix, opt_iy, iflip, irot, opt_iflip, nr_mir, iflip_start, iflip_stop, count;
    double ropt2, maxCC, CC, dxx, dyy;
    std::vector<Matrix2D<double> > Mflip, dum;
    Matrix1D<double> trans(2);
    Matrix2D<double> Maux2, Maux;
    std::vector<Matrix2D<double> > Finv;

    if (do_mirror) nr_mir = 2;
    else nr_mir = 1;

    Moffsets.resize(dim, dim);
    Moffsets.setXmippOrigin();
    Moffsets.init_constant(-1);
    Moffsets_mirror.resize(dim, dim);
    Moffsets_mirror.setXmippOrigin();
    Moffsets_mirror.init_constant(-1);
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Mimg_trans.clear();

    // Pre-calculate flipped image matrices
    Mflip.clear();
    FOR_ALL_FLIPS()
    {
        Maux.setXmippOrigin();
        applyGeometry(Maux, F[iflip], Mimg, IS_INV, WRAP);
        Mflip.push_back(Maux);
    }

    // Calculate inverted flipping matrices
    // This is to have the offsets consistent with those in fast_mode
    FOR_ALL_FLIPS()
    {
        Maux = F[iflip].inv();
        Finv.push_back(Maux);
    }

    // If offsets > max_shift: reset to zero...
    count = 0;
    FOR_ALL_MODELS()
    {
        for (int imir = 0; imir < nr_mir; imir++)
        {
            irefmir = imir * n_ref + refno;
            iflip_start = imir * nr_nomirror_flips;
            iflip_stop = imir * nr_nomirror_flips + nr_nomirror_flips;
            ropt2 = offsets[2*irefmir] * offsets[2*irefmir] + offsets[2*irefmir+1] * offsets[2*irefmir+1];
            if (ropt2 > max_shift*max_shift)
            {
                offsets[2*irefmir] = 0.;
                offsets[2*irefmir+1] = 0.;
            }
            if (!limit_rot || pdf_directions[refno] > 0.)
            {
                FOR_ALL_LIMITED_TRANSLATIONS()
                {
                    ix = ROUND(offsets[2*irefmir] + Vtrans[itrans](0));
                    iy = ROUND(offsets[2*irefmir+1] + Vtrans[itrans](1));
                    dxx = (double)intWRAP(ix, Moffsets.startingX(), Moffsets.finishingX());
                    dyy = (double)intWRAP(iy, Moffsets.startingY(), Moffsets.finishingY());
                    // For non-mirrors
                    if (imir == 0 && MAT_ELEM(Moffsets, ROUND(dyy), ROUND(dxx)) < 0)
                    {
                        Mimg_trans.push_back(dum);
                        for (int iflip = iflip_start; iflip < iflip_stop; iflip++)
                        {
                            trans(0) = dxx * DIRECT_MAT_ELEM(Finv[iflip], 0, 0) + dyy * DIRECT_MAT_ELEM(Finv[iflip], 0, 1);
                            trans(1) = dxx * DIRECT_MAT_ELEM(Finv[iflip], 1, 0) + dyy * DIRECT_MAT_ELEM(Finv[iflip], 1, 1);
                            Mflip[iflip].translate(trans, Maux, WRAP);
                            Mimg_trans[count].push_back(Maux);
                        }
                        MAT_ELEM(Moffsets, ROUND(dyy), ROUND(dxx)) = count;
                        count++;
                    }
                    // For mirrors use a separate offset-matrix
                    else if (imir == 1 && MAT_ELEM(Moffsets_mirror, ROUND(dyy), ROUND(dxx)) < 0)
                    {
                        Mimg_trans.push_back(dum);
                        for (int iflip = iflip_start; iflip < iflip_stop; iflip++)
                        {
                            trans(0) = dxx * DIRECT_MAT_ELEM(Finv[iflip], 0, 0) + dyy * DIRECT_MAT_ELEM(Finv[iflip], 0, 1);
                            trans(1) = dxx * DIRECT_MAT_ELEM(Finv[iflip], 1, 0) + dyy * DIRECT_MAT_ELEM(Finv[iflip], 1, 1);
                            Mflip[iflip].translate(trans, Maux, WRAP);
                            Mimg_trans[count].push_back(Maux);
                        }
                        MAT_ELEM(Moffsets_mirror, ROUND(dyy), ROUND(dxx)) = count;
                        count++;
                    }
                }
            }
        }
    }

}

void Prog_MLalign2D_prm::ML_integrate_locally(
    Matrix2D<double> &Mimg, std::vector <std::vector< Matrix2D<double> > > &Mref,
    std::vector <std::vector< Matrix2D<double> > > &Mwsum_imgs,
    double &wsum_sigma_noise, double &wsum_sigma_offset,
    std::vector<double> &sumw, std::vector<double> &sumw2, std::vector<double> &sumw_mirror,
    double &LL, double &fracweight,
    int &opt_refno, double &opt_psi,
    Matrix1D<double> &opt_offsets, std::vector<double> &opt_offsets_ref,
    std::vector<double> &pdf_directions)
{

    Matrix3D<double> Mweight;
    Matrix2D<double> Maux, Mdzero, Mtrans;
    std::vector<Matrix2D<double> > Mflip;
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Mtrans.resize(Maux);
    Matrix2D<int> Moffsets, Moffsets_mirror;
    std::vector<std::vector<Matrix2D<double> > > Mimg_trans;
    std::vector<double> refw(n_ref), refw2(n_ref), refw_mirror(n_ref);
    double sigma_noise2, Xi2, aux, fracpdf, A2_plus_Xi2, weight, weight1, weight2;
    double maxw, pdf, diff, dfsigma2, mindiff2 = 99.e99;
    double sum, wsum_corr = 0., sum_refw = 0., sum_refw2 = 0., maxweight = -99.e99;
    int point_trans, nr_mir, irefmir, iflip_start, iflip_stop, ixx, iyy, irot;
    int opt_ipsi = 0, opt_iflip = 0, opt_irefmir = 0, opt_itrans = 0;
    int ii, imax = n_ref * nr_flip / nr_nomirror_flips;
    std::vector<double> Pmax_refmir(imax);
    for (int i = 0; i < imax; i++) Pmax_refmir[i] = -99.e99;
    sigma_noise2 = sigma_noise * sigma_noise;
    dfsigma2 = df * sigma_noise2;

    // Calculate all flipped and translated versions of Mimg
    calculate_realspace_offsets(Mimg, opt_offsets_ref, pdf_directions,
                                Mimg_trans, Moffsets, Moffsets_mirror);

    // Calculate all squared differences & mindiff2 (for optimal offsets only)
    Mweight.initZeros(nr_trans, n_ref, nr_flip*nr_psi);
    Xi2 = Mimg.sum2();
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            A2_plus_Xi2 = 0.5 * (A2[refno] + Xi2);
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                ixx = ROUND(opt_offsets_ref[2*irefmir]);
                iyy = ROUND(opt_offsets_ref[2*irefmir+1]);
                if (iflip < nr_nomirror_flips) point_trans = MAT_ELEM(Moffsets, iyy, ixx);
                else point_trans = MAT_ELEM(Moffsets_mirror, iyy, ixx);
                Mtrans = Mimg_trans[point_trans][iflip%nr_nomirror_flips];
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    Maux = Mref[refno][ipsi];
                    diff = A2_plus_Xi2;
                    for (int ii = 0; ii < dim2; ii++)
                    {
                        diff -= (Maux).data[ii] * (Mtrans).data[ii];
                    }
                    dVkij(Mweight, zero_trans, refno, irot) = diff;
		    if (debug) {
			std::cout << 360. + (-psi_step * (iflip * nr_psi + ipsi) - SMALLANGLE)<<" "<<diff<<" "<<A2_plus_Xi2<<" "<<diff-A2_plus_Xi2<<" "<<A2[refno]<<" "<<Xi2<<" "<<-ixx<<" "<< -iyy<<std::endl;
		    }
                    if (diff < mindiff2) mindiff2 = diff;
                }
            }
        }
    }

    // Now that we have mindiff2 calculate the weighting matrices and maxweight
    FOR_ALL_MODELS()
    {
        refw[refno] = 0.;
        refw_mirror[refno] = 0.;
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                ixx = ROUND(opt_offsets_ref[2*irefmir]);
                iyy = ROUND(opt_offsets_ref[2*irefmir+1]);
                if (iflip < nr_nomirror_flips) fracpdf = (1. - mirror_fraction[refno]);
                else fracpdf = mirror_fraction[refno];
                pdf = fracpdf * alpha_k[refno] * pdf_directions[refno] * MAT_ELEM(P_phi, iyy, ixx);
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    diff = dVkij(Mweight, zero_trans, refno, irot);
		    if (!do_student)
		    {
			// normal distribution
			aux = (diff - mindiff2) / (anneal * sigma_noise2);
			// next line because of numerical precision of exp-function
			if (aux > 1000.) weight = 0.;
			else weight = exp(-aux) * pdf;
			// Store running sum for new sigma2
			wsum_corr += weight * diff;
			// Store weight
			dVkij(Mweight, zero_trans, refno, irot) = weight;
			// Accumulate sum weights
			if (iflip < nr_nomirror_flips) refw[refno] += weight;
			else refw_mirror[refno] += weight;
			sum_refw += weight;
		    }
		    else
		    {
			// t-student distribution
			// pdf = (1 + diff2/sigma2*df)^df2
			// Correcting for mindiff2:
			// pdfc = (1 + diff2/sigma2*df)^df2 / (1 + mindiff2/sigma2*df)^df2
			//      = ( (1 + diff2/sigma2*df)/(1 + mindiff2/sigma2*df) )^df2
			//      = ( (sigma2*df + diff2) / (sigma2*df + mindiff2) )^df2
			// Extra factor two because we had saved 0.5*diff2!!
			aux = (dfsigma2 + 2. * diff) / (dfsigma2 + 2. * mindiff2);
			weight1 = pow(aux, df2) * pdf;
			// Calculate extra weight acc. to Eq (10) Wang et al.
			// Patt. Recognition Lett. 25, 701-710 (2004)
			weight2 = ( df + dim2 ) / ( df + (2. * diff / sigma_noise2) );
			// total weight
			weight = weight1 * weight2;
			// Store running sum for new sigma2
			wsum_corr += weight * diff;
			// Store weight
			dVkij(Mweight, zero_trans, refno, irot) = weight;
			// Accumulate sum weights
			if (iflip < nr_nomirror_flips) refw[refno] += weight1;
			else refw_mirror[refno] += weight1;
			sum_refw += weight1;
			refw2[refno] += weight;
			sum_refw2 += weight;
		    }
                    if (weight > Pmax_refmir[irefmir]) Pmax_refmir[irefmir] = weight;
                    if (weight > maxweight)
                    {
                        maxweight = weight;
			if (do_student) maxweight2 = weight2;
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

    // Now for all irefmir, check significant rotations
    // and calculate their limited_translations probabilities
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            A2_plus_Xi2 = 0.5 * (A2[refno] + Xi2);
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                if (iflip < nr_nomirror_flips) fracpdf = (1. - mirror_fraction[refno]);
                else fracpdf = mirror_fraction[refno];
                fracpdf *= alpha_k[refno] * pdf_directions[refno];
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
		    // Note that for t-students distribution we now do something 
		    // a little bit different compared to ML_integrate_complete! 
		    // Instead of comparing "weight1", we compare "weight1*weight2"!!
                    if (dVkij(Mweight, zero_trans, refno, irot) > C_fast*Pmax_refmir[irefmir])
                    {
                        Maux = Mref[refno][ipsi];
                        // expand for all limited translations
                        FOR_ALL_LIMITED_TRANSLATIONS()
                        {
                            if (itrans != zero_trans)
                            {
                                ixx = ROUND(opt_offsets_ref[2*irefmir] + Vtrans[itrans](0));
                                iyy = ROUND(opt_offsets_ref[2*irefmir+1] + Vtrans[itrans](1));
                                if (iflip < nr_nomirror_flips) point_trans = MAT_ELEM(Moffsets, iyy, ixx);
                                else point_trans = MAT_ELEM(Moffsets_mirror, iyy, ixx);
                                Mtrans = Mimg_trans[point_trans][iflip%nr_nomirror_flips];
                                diff = A2_plus_Xi2;
				pdf = fracpdf * MAT_ELEM(P_phi, iyy, ixx);
                                for (int ii = 0; ii < dim2; ii++)
                                {
                                    diff -= (Maux).data[ii] * (Mtrans).data[ii];
                                }
				if (!do_student)
				{
				    // Normal distribution
				    aux = (diff - mindiff2) / (anneal * sigma_noise2);
				    // next line because of numerical precision of exp-function
				    if (aux > 1000.) weight = 0.;
				    else weight = exp(-aux) * pdf;
				    wsum_corr += weight * diff;
				    dVkij(Mweight, itrans, refno, irot) = weight;
				    // Accumulate sum weights
				    if (iflip < nr_nomirror_flips) refw[refno] += weight;
				    else refw_mirror[refno] += weight;
				    sum_refw += weight;
				}
				else
				{
				    // t-student distribution
				    // pdf = (1 + diff2/sigma2*df)^df2
				    // Correcting for mindiff2:
				    // pdfc = (1 + diff2/sigma2*df)^df2 / (1 + mindiff2/sigma2*df)^df2
				    //      = ( (1 + diff2/sigma2*df)/(1 + mindiff2/sigma2*df) )^df2
				    //      = ( (sigma2*df + diff2) / (sigma2*df + mindiff2) )^df2
				    // Extra factor two because we saved 0.5*diff2!!
				    aux = (dfsigma2 + 2. * diff) / (dfsigma2 + 2. * mindiff2);
				    weight1 = pow(aux, df2) * pdf;
				    // Calculate extra weight acc. to Eq (10) Wang et al.
				    // Patt. Recognition Lett. 25, 701-710 (2004)
				    weight2 = ( df + dim2 ) / ( df + (2. * diff / sigma_noise2) );
				    // Total weight
				    weight = weight1 * weight2;
				    // Store running sum for new sigma2
				    wsum_corr += weight * diff;
				    dVkij(Mweight, itrans, refno, irot) = weight;
				    // Accumulate sum weights
				    if (iflip < nr_nomirror_flips) refw[refno] += weight1;
				    else refw_mirror[refno] += weight1;
				    sum_refw += weight1;
				    refw2[refno] += weight;
				    sum_refw2 += weight;
				}
                                if (weight > maxweight)
                                {
                                    maxweight = weight;
				    if (do_student) maxweight2 = weight2;
                                    opt_refno = refno;
                                    opt_ipsi = ipsi;
                                    opt_iflip = iflip;
                                    opt_itrans = itrans;
                                    opt_irefmir = irefmir;
                                }
				if (debug) {
				    std::cout << 360. + (-psi_step * (iflip * nr_psi + ipsi) - SMALLANGLE)<<" "<<diff<<" "<<A2_plus_Xi2<<" "<<diff-A2_plus_Xi2<<" "<<A2[refno]<<" "<<Xi2<<" "<<-ixx<<" "<< -iyy<<std::endl;
				}
                            }
                        }
                    }
                }
            }
        }
    }

    // Normalize all weighted sums by sum_refw such that sum over all weights is one!
    // And accumulate weighted sums
    wsum_sigma_noise += (2 * wsum_corr / sum_refw);
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            sumw[refno] += (refw[refno] + refw_mirror[refno]) / sum_refw;
	    sumw2[refno] += refw2[refno] / sum_refw;
            sumw_mirror[refno] += refw_mirror[refno] / sum_refw;
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    FOR_ALL_LIMITED_TRANSLATIONS()
                    {
                        weight = dVkij(Mweight, itrans, refno, irot);
                        if (weight > SIGNIFICANT_WEIGHT_LOW*maxweight)
                        {
                            weight /= sum_refw;
                            ixx = ROUND(opt_offsets_ref[2*irefmir] + Vtrans[itrans](0));
                            iyy = ROUND(opt_offsets_ref[2*irefmir+1] + Vtrans[itrans](1));
                            if (iflip < nr_nomirror_flips) point_trans = MAT_ELEM(Moffsets, iyy, ixx);
                            else point_trans = MAT_ELEM(Moffsets_mirror, iyy, ixx);
			    wsum_sigma_offset += weight * (double)(ixx * ixx + iyy * iyy);
			    Mwsum_imgs[refno][ipsi] += weight * Mimg_trans[point_trans][iflip%nr_nomirror_flips];
                        }
                    }
                }
            }
        }
    }

    // Update the optimal origin in-plane transformations
    opt_offsets(0) = opt_offsets_ref[2*opt_irefmir] + Vtrans[opt_itrans](0);
    opt_offsets(1) = opt_offsets_ref[2*opt_irefmir+1] + Vtrans[opt_itrans](1);
    // Sjors 16jul07: repaired bug for -save_memC?
    if (save_mem3)
        opt_psi = -opt_iflip * 360. / nr_flip - SMALLANGLE;
    else
	opt_psi = -psi_step * (opt_iflip * nr_psi + opt_ipsi) - SMALLANGLE;
    if (!do_student)
	fracweight = maxweight / sum_refw;
    else
	fracweight = maxweight / sum_refw2;

    if (do_mirror) nr_mir = 2;
    else nr_mir = 1;
    for (int i = 0; i < imax; i++) Pmax_refmir[i] = -99.e99;
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            for (int imir = 0; imir < nr_mir; imir++)
            {
                irefmir = imir * n_ref + refno;
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
                            weight = dVkij(Mweight, itrans, refno, irot);
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

    // Compute Log Likelihood
    if (!do_student)
	// 1st term: log(refw_i)
	// 2nd term: for subtracting mindiff2
	// 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
	LL += log(sum_refw) - mindiff2 / sigma_noise2 - dim * dim * log(2.50663 * sigma_noise);
    else
	// 1st term: log(refw_i)
	// 2nd term: for dividing by mindiff2
	// 3rd term: for constant term in t-student distribution
	// (ignoring constants due to gamma functions!)
	LL += log(sum_refw) - log(1. + ( mindiff2 / df )) - log(sigma_noise * sigma_noise);

}


// Maximum Likelihood calculation for one image ============================================
// Integration over all translation, given  model and in-plane rotation
void Prog_MLalign2D_prm::ML_integrate_complete(
    Matrix2D<double> &Mimg, std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref,
    Matrix2D<int> &Msignificant,
    std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fwsum_imgs,
    double &wsum_sigma_noise, double &wsum_sigma_offset,
    std::vector<double> &sumw, std::vector<double> &sumw2, std::vector<double> &sumw_mirror,
    double &LL, double &fracweight, double &maxweight2,
    int &opt_refno, double &opt_psi, Matrix1D<double> &opt_offsets,
    std::vector<double> &opt_offsets_ref, std::vector<double> &pdf_directions)
{

    Matrix2D<double> Maux, Mdzero;
    Matrix2D<std::complex<double> > Fimg, Faux, Faux2;
    std::vector<Matrix2D<std::complex<double> > > Fimg_flip;
    std::vector<Matrix2D<double> > dumM;
    std::vector<bool> dumb;
    std::vector <std::vector< Matrix2D<double> > > Mweight;
    std::vector<double> refw(n_ref), refw2(n_ref), refw_mirror(n_ref);
    double sigma_noise2, XiA, Xi2, aux, pdf, fracpdf, A2_plus_Xi2, maxw, mind, dfsigma2, mindiff2 = 99.e99;
    double weight, weight1, weight2, diff, sum, sum2, wsum_corr = 0., sum_refw = 0., sum_refw2 = 0., maxweight = -99.e99;
    int irot, irefmir, sigdim, xmax, ymax;
    int ioptx = 0, iopty = 0, ioptpsi = 0, ioptflip = 0, imax = 0;
    if (fast_mode) imax = n_ref * nr_flip / nr_nomirror_flips;
    std::vector<int> ioptx_ref(imax), iopty_ref(imax), ioptflip_ref(imax);
    std::vector<double> maxw_ref(imax);


    /* Not to store all 360-degrees rotations of the references (and pdf, Fwsum_imgs etc.) in memory,
       the experimental image is rotated over 0, 90, 180 & 270 degrees (called FLIPS), and only
       rotations over 90 degrees of the references are stored.
       This save a factor of 4 in memory requirements, and these FLIPS do not require interpolation
       of the experiemental image and thereby deterioration of the process.
       If do_mirror there are 8 flips, now also including all mirrored versions.
       Rotation is done in real space, and then the Fourier Transform is calculated four/eight times.
       In principle, rotation could be done in Fourier space, but then there is a problem with
       even-sized images, where the origin is not exactly in the center, and 1 pixel wrapping is required.
       Anyway, the total number of (I)FFT's is determined in much greater extent by n_ref and n_rot!
    */


    // Only translations smaller than save_mem2 (default=6) sigma_offset are considered!
    // This saves a lot of memory and CPU! (typically a factor 2, depending on sigma_offset vs. dim)
    if (save_mem2) sigdim = 2 * CEIL(sigma_offset * 3);
    else sigdim = 2 * CEIL(sigma_offset * 6);
    sigdim++; // (to get uneven number)
    sigdim = XMIPP_MIN(dim, sigdim);
    maxweight2 = 0.;

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
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Fimg_flip.clear();
    Mweight.clear();
    sigma_noise2 = sigma_noise * sigma_noise;
    dfsigma2 = df * sigma_noise2;
    Xi2 = Mimg.sum2();

    // Flip images and calculate correlation matrices and maximum correlation
    FOR_ALL_FLIPS()
    {
        Maux.setXmippOrigin();
        applyGeometry(Maux, F[iflip], Mimg, IS_INV, WRAP);
        FourierTransformHalf(Maux, Fimg);
        Fimg *= dim * dim;
        Fimg_flip.push_back(Fimg);
        FOR_ALL_MODELS()
        {
            Mweight.push_back(dumM);
            if (!limit_rot || pdf_directions[refno] > 0.)
            {
                A2_plus_Xi2 = 0.5 * (A2[refno] + Xi2);
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
		    Mweight[refno].push_back(Mdzero);
                    if (dMij(Msignificant, refno, irot))
                    {
                        multiplyElements(Fimg, Fref[refno][ipsi], Faux);
                        InverseFourierTransformHalf(Faux, Maux, dim);
                        Maux /= dim * dim;
                        CenterFFT(Maux, true);
                        Mweight[refno][irot].resize(sigdim, sigdim);
                        Mweight[refno][irot].setXmippOrigin();
                        FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot])
                        {
                            MAT_ELEM(Mweight[refno][irot], i, j) = A2_plus_Xi2 - MAT_ELEM(Maux, i, j);
                        }
                        mind = Mweight[refno][irot].compute_min();
                        if (mind < mindiff2) mindiff2 = mind;
                    }
                    else
                    {
                        Mweight[refno].push_back(Mdzero);
                    }
                }
            }
        }
    }

    // Now that we have mindiff2 calculate the weighting matrices and maxweight
    FOR_ALL_MODELS()
    {
        refw[refno] = 0.;
        refw_mirror[refno] = 0.;
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            FOR_ALL_ROTATIONS()
            {
                FOR_ALL_FLIPS()
                {
                    irot = iflip * nr_psi + ipsi;
                    irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                    if (dMij(Msignificant, refno, irot))
                    {
                        if (iflip < nr_nomirror_flips) fracpdf = alpha_k[refno] * (1. - mirror_fraction[refno]);
                        else fracpdf = alpha_k[refno] * mirror_fraction[refno];
                        sum = sum2 = 0.;
			FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot])
			{
			    diff = MAT_ELEM(Mweight[refno][irot], i, j);
			    pdf = fracpdf * MAT_ELEM(P_phi, i, j);
			    if (!do_student)
			    {
				// Normal distribution
				aux = (diff - mindiff2) / (anneal * sigma_noise2);
				// next line because of numerical precision of exp-function
				if (aux > 1000.) weight = 0.;
				else weight = exp(-aux) * pdf;
				// calculate weighted sum of (X-A)^2 for sigma_noise update
				wsum_corr += weight * diff;
				// store weight
				MAT_ELEM(Mweight[refno][irot], i, j) = weight;
				// Accumulate sum weights
				if (iflip < nr_nomirror_flips) refw[refno] += weight;
				else refw_mirror[refno] += weight;
				sum_refw += weight;
			    }
			    else
			    {
				// t-student distribution
				aux = (dfsigma2 + 2. * diff) / (dfsigma2 + 2. * mindiff2);
				weight1 = pow(aux, df2) * pdf;
				// Calculate extra weight acc. to Eq (10) Wang et al.
				// Patt. Recognition Lett. 25, 701-710 (2004)
				weight2 = ( df + dim2 ) / ( df + (2. * diff / sigma_noise2) );
				// Total weight
				weight = weight1 * weight2;
				// Store running sum for new sigma_noise
				wsum_corr += weight * diff;
				// Store probability weights
				MAT_ELEM(Mweight[refno][irot], i, j) = weight;
				// Accumulate sum weights
				if (iflip < nr_nomirror_flips) refw[refno] += weight1;
				else refw_mirror[refno] += weight1;
				sum_refw += weight1;
				refw2[refno] += weight;
				sum_refw2+= weight;
			    }
			    if (weight > maxweight)
			    {
				maxweight = weight;
				if (do_student) maxweight2 = weight2;
				iopty = i;
				ioptx = j;
				ioptpsi = ipsi;
				ioptflip = iflip;
				opt_refno = refno;
			    }
			    if (fast_mode && weight > maxw_ref[irefmir])
			    {
				maxw_ref[irefmir] = weight;
				iopty_ref[irefmir] = i;
				ioptx_ref[irefmir] = j;
				ioptflip_ref[irefmir] = iflip;
			    }
			}
		    }
                }
            }
        }
    }

    // Normalize all weighted sums by sum_refw such that sum over all weights is one!
    // And accumulate the FT of the weighted, shifted images.
    wsum_sigma_noise += (2 * wsum_corr / sum_refw);
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            sumw[refno] += (refw[refno] + refw_mirror[refno]) / sum_refw;
	    sumw2[refno] += refw2[refno] / sum_refw;
            sumw_mirror[refno] += refw_mirror[refno] / sum_refw;
            FOR_ALL_ROTATIONS()
            {
                FOR_ALL_FLIPS()
                {
                    irot = iflip * nr_psi + ipsi;
                    if (dMij(Msignificant, refno, irot))
                    {
                        if (Mweight[refno][irot].compute_max() > SIGNIFICANT_WEIGHT_LOW*maxweight)
                        {
                            Mweight[refno][irot] /= sum_refw;
                            // Use Maux, because Mweight is smaller than dim x dim!
                            Maux.initZeros();
			    FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot])
			    {
				MAT_ELEM(Maux, i, j) = MAT_ELEM(Mweight[refno][irot], i, j);
				wsum_sigma_offset += MAT_ELEM(Maux, i, j) * MAT_ELEM(Mr2, i, j);
			    }
			    FourierTransformHalf(Maux, Faux);
                            Faux *= dim * dim;
                            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
                            {
                                dMij(Fwsum_imgs[refno][ipsi], i, j) += conj(dMij(Faux, i, j)) * dMij(Fimg_flip[iflip], i, j);
                            }
                        }
                    }
                }
            }
        }
    }

    // Calculate optimal transformation parameters
    if (fast_mode)
    {
        for (int i = 0; i < imax;i++)
        {
            opt_offsets_ref[2*i] = -(double)ioptx_ref[i] * DIRECT_MAT_ELEM(F[ioptflip_ref[i]], 0, 0) -
                                   (double)iopty_ref[i] * DIRECT_MAT_ELEM(F[ioptflip_ref[i]], 0, 1);
            opt_offsets_ref[2*i+1] = -(double)ioptx_ref[i] * DIRECT_MAT_ELEM(F[ioptflip_ref[i]], 1, 0) -
                                     (double)iopty_ref[i] * DIRECT_MAT_ELEM(F[ioptflip_ref[i]], 1, 1);
        }
    }
    opt_offsets(0) = -(double)ioptx * DIRECT_MAT_ELEM(F[ioptflip], 0, 0) -
                     (double)iopty * DIRECT_MAT_ELEM(F[ioptflip], 0, 1);
    opt_offsets(1) = -(double)ioptx * DIRECT_MAT_ELEM(F[ioptflip], 1, 0) -
                     (double)iopty * DIRECT_MAT_ELEM(F[ioptflip], 1, 1);

    // Sjors 16jul07: repaired bug for -save_memC?
    if (save_mem3)
        opt_psi = -ioptflip * 360. / nr_flip - SMALLANGLE;
    else
	opt_psi = -psi_step * (ioptflip * nr_psi + ioptpsi) - SMALLANGLE;
    if (!do_student)
	fracweight = maxweight / sum_refw;
    else
	fracweight = maxweight / sum_refw2;

    if (!do_student)
	// 1st term: log(refw_i)
	// 2nd term: for subtracting mindiff2
	// 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
	LL += log(sum_refw) - mindiff2 / sigma_noise2 - dim * dim * log(2.50663 * sigma_noise);
    else
	// 1st term: log(refw_i)
	// 2nd term: for dividing by mindiff2
	// 3rd term: for constant term in t-student distribution
	// (ignoring constants due to gamma functions!)
	LL += log(sum_refw) - log(1. + ( mindiff2 / df )) - log(sigma_noise * sigma_noise);

}

void Prog_MLalign2D_prm::maxCC_search_complete(Matrix2D<double> &Mimg,
        std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref,
        std::vector <std::vector< Matrix2D<double> > > &Mref,
        double &search_shift, std::vector <std::vector< Matrix2D<double> > > &Msum_imgs,
        std::vector<double> &sumw, std::vector<double> &sumw_mirror,
        double &maxCC, int &opt_refno, double &opt_psi, Matrix1D<double> &opt_offsets,
        std::vector<double> &pdf_directions)
{

    Matrix2D<double> Maux, Maux2;
    Matrix2D<std::complex<double> > Fimg, Faux;
    double sigma_noise2, aux, avg, std, CC;
    int irot, sigdim, xmax = 0, ymax = 0;
    int ioptx = 0, iopty = 0, ioptpsi = 0, ioptflip = 0, imax = 0;
    double stddev_img, mean_img, dummy;

    maxCC = -99.e99;
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Maux2.resize(dim, dim);
    Maux2.setXmippOrigin();
    Faux.resize(hdim + 1, dim);
    sigma_noise2 = sigma_noise * sigma_noise;
    Matrix2D<int> shiftmask;

    if (search_shift > 0.)
    {
        shiftmask.resize(dim, dim);
        shiftmask.setXmippOrigin();
        BinaryCircularMask(shiftmask, search_shift, INNER_MASK);
    }

    Mimg.computeStats(mean_img, stddev_img, dummy, dummy);
    Maux2 = Mimg;
    Maux2 -= mean_img;

    // Flip images and calculate correlation matrices and maximum correlation
    FOR_ALL_FLIPS()
    {
        applyGeometry(Maux, F[iflip], Maux2, IS_INV, WRAP);
        if (search_shift > 0.)
        {
            FourierTransformHalf(Maux, Fimg);
            Fimg *= dim * dim;
        }
        FOR_ALL_MODELS()
        {
            if (!limit_rot || pdf_directions[refno] > 0.)
            {
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    if (search_shift > 0.)
                    {
                        multiplyElements(Fimg, Fref[refno][ipsi], Faux);
                        InverseFourierTransformHalf(Faux, Maux, dim);
                        Maux /= dim * dim;
                        CenterFFT(Maux, true);
                        apply_binary_mask(shiftmask, Maux, Maux, 0.);
                        Maux.maxIndex(ymax, xmax);
                        CC = MAT_ELEM(Maux, ymax, xmax);
                    }
                    else
                    {
                        CC = 0.;
                        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Maux)
                        {
                            CC += dMij(Maux, i, j) * dMij(Mref[refno][ipsi], i, j);
                        }
                    }
                    CC /= A2[refno] * stddev_img; // For maxCC-mode, A2[refno] holds stddev_ref!
                    if (CC > maxCC)
                    {
                        maxCC = CC;
                        iopty = ymax;
                        ioptx = xmax;
                        ioptpsi = ipsi;
                        ioptflip = iflip;
                        opt_refno = refno;
                    }
                }
            }
        }
    }
    maxCC /= dim * dim;

    // Calculate optimal transformation parameters
    opt_offsets(0) = -(double)ioptx * DIRECT_MAT_ELEM(F[ioptflip], 0, 0)
                     - (double)iopty * DIRECT_MAT_ELEM(F[ioptflip], 0, 1);
    opt_offsets(1) = -(double)ioptx * DIRECT_MAT_ELEM(F[ioptflip], 1, 0)
                     - (double)iopty * DIRECT_MAT_ELEM(F[ioptflip], 1, 1);
    if (save_mem3)
        opt_psi = -ioptflip * 360. / nr_flip - SMALLANGLE;
    else
        opt_psi = -psi_step * (ioptflip * nr_psi + ioptpsi) - SMALLANGLE;

    // Store sums of the aligned images
    Mimg.translate(opt_offsets, Maux, true);
    applyGeometry(Maux2, F[ioptflip], Maux, IS_INV, WRAP);
    Msum_imgs[opt_refno][ioptpsi] += Maux2;
    sumw[opt_refno] += 1.;
    if (ioptflip > 3) sumw_mirror[opt_refno] += 1.;

}


void Prog_MLalign2D_prm::ML_sum_over_all_images(SelFile &SF, std::vector< ImageXmippT<double> > &Iref, int iter,
        double &LL, double &sumcorr, DocFile &DFo,
        std::vector<Matrix2D<double> > &wsum_Mref,
        double &wsum_sigma_noise, double &wsum_sigma_offset, 
	std::vector<double> &sumw, std::vector<double> &sumw2, std::vector<double> &sumw_mirror)
{

    ImageXmipp img;
    SelLine line;
    FileName fn_img, fn_trans;
    std::vector <std::vector< Matrix2D<double> > > Mref, Msum_imgs;
    std::vector <std::vector< Matrix2D<std::complex<double> > > > Fref, Fwsum_imgs;
    std::vector<Matrix2D<std::complex <double> > > dum;
    std::vector<Matrix2D<double> > dum2;
    std::vector<double> allref_offsets, pdf_directions(n_ref);
    Matrix2D<double>  Mdzero;
    Matrix2D<std::complex<double> >  Fdzero;
    Matrix2D<int> Msignificant;
    Msignificant.resize(n_ref, nr_psi*nr_flip);
    Matrix1D<double> dataline(9), opt_offsets(2), trans(2);

    float old_phi = -999., old_theta = -999.;
    double opt_psi, opt_flip, maxcorr, maxweight2;
    double opt_xoff, opt_yoff;
    int c, nn, imgno, opt_refno;
    bool fill_real_space, fill_fourier_space;


    // Generate (FT of) each rotated version of all references
    if (limit_trans || (maxCC_rather_than_ML && !(search_shift > 0.)))
    {
        fill_real_space = true;
        fill_fourier_space = false;
    }
    else if (fast_mode)
    {
        fill_fourier_space = true;
        if (save_mem1) fill_real_space = false;
        else fill_real_space = true;
    }
    else
    {
        fill_real_space = false;
        fill_fourier_space = true;
    }
    rotate_reference(Iref, fill_real_space, fill_fourier_space, Mref, Fref);

    // Initialize
    LL = 0.;
    nn = SF.ImgNo();
    if (verb > 0) init_progress_bar(nn);
    c = XMIPP_MAX(1, nn / 60);
    Fwsum_imgs.clear();
    Msum_imgs.clear();
    sumw.clear();
    sumw2.clear();
    sumw_mirror.clear();
    Fdzero.resize(hdim + 1, dim);
    Mdzero.resize(dim, dim);
    Mdzero.setXmippOrigin();
    LL = 0.;
    wsum_sigma_noise = 0.;
    wsum_sigma_offset = 0.;
    sumcorr = 0.;
    trans.initZeros();
    FOR_ALL_MODELS()
    {
        sumw.push_back(0.);
        sumw2.push_back(0.);
        sumw_mirror.push_back(0.);
        if ((maxCC_rather_than_ML || limit_trans))
        {
            Msum_imgs.push_back(dum2);
            FOR_ALL_ROTATIONS()
            {
                Msum_imgs[refno].push_back(Mdzero);
            }
        }
	else
	{
	    Fwsum_imgs.push_back(dum);
	    FOR_ALL_ROTATIONS()
	    {
		Fwsum_imgs[refno].push_back(Fdzero);
	    }
	}
    }   

    // Loop over all images
    imgno = 0;
    SF.go_beginning();
    while ((!SF.eof()))
    {

	// Check whether to kill job
	exit_if_not_exists(fn_control);

        fn_img = SF.NextImg();
        fn_trans = fn_img.remove_directories(offsets_keepdir);
        fn_trans = fn_root + "_offsets/" + fn_trans + ".off";

        img.read(fn_img, false, false, false, false);
        img().setXmippOrigin();
        if (fn_doc != "")
        {
            trans(0) = (double)ROUND(imgs_oldxoff[imgno]);
            trans(1) = (double)ROUND(imgs_oldyoff[imgno]);
            img().selfTranslate(trans, true);
        }

        // Get optimal offsets for all references
        if (fast_mode || limit_trans)
        {
            if (do_write_offsets)
            {
                // read from disc
                if (!read_offsets(fn_trans, allref_offsets))
                {
                    int itot = n_ref * 2;
                    if (do_mirror) itot *= 2;
                    allref_offsets.clear();
                    allref_offsets.resize(itot);
                    if (zero_offsets) for (int i = 0; i < itot; i++) allref_offsets[i] = 0.;
                    else for (int i = 0; i < itot; i++) allref_offsets[i] = -999.;
                }
            }
            else
            {
                // get from memory
                allref_offsets = imgs_offsets[imgno];
            }
        }

        // Read optimal orientations from memory
        if (limit_rot)
        {
            old_phi = imgs_oldphi[imgno];
            old_theta = imgs_oldtheta[imgno];
        }

        // For limited orientational search: preselect relevant directions
        preselect_directions(old_phi, old_theta, pdf_directions);

        if (maxCC_rather_than_ML)
        {
            // A. Use a maximum cross-correlation target function

            maxCC_search_complete(img(), Fref, Mref, search_shift, Msum_imgs, sumw, sumw_mirror,
                                  maxcorr, opt_refno, opt_psi, opt_offsets, pdf_directions);

        }
	else if (limit_trans)
        {
            // B. Use a maximum-likelihood target function in real space
            //    with limited translational searches

            ML_integrate_locally(img(), Mref, Msum_imgs, wsum_sigma_noise, wsum_sigma_offset,
                                 sumw, sumw2, sumw_mirror, LL, maxcorr, maxweight2,
                                 opt_refno, opt_psi, opt_offsets,
                                 allref_offsets, pdf_directions);

        }
        else
        {
            // C. Use a maximum-likelihood target function in real space
            //    with complete or reduced-space translational searches (-fast)

            if (fast_mode) preselect_significant_model_phi(img(), allref_offsets, Mref,
                        Msignificant, pdf_directions);
            else Msignificant.init_constant(1);
            ML_integrate_complete(img(), Fref, Msignificant,
                                  Fwsum_imgs, wsum_sigma_noise, wsum_sigma_offset, 
                                  sumw, sumw2, sumw_mirror, LL, maxcorr, maxweight2,
                                  opt_refno, opt_psi, opt_offsets,
                                  allref_offsets, pdf_directions);

        }

        // Write optimal offsets for all references to disc
        if (fast_mode || limit_trans)
        {
            if (do_write_offsets) write_offsets(fn_trans, allref_offsets);
            else imgs_offsets[imgno] = allref_offsets;
        }

        // Store optimal phi and theta in memory
        if (limit_rot)
        {
            imgs_oldphi[imgno] = Iref[opt_refno].Phi();
            imgs_oldtheta[imgno] = Iref[opt_refno].Theta();
        }

        // Output odcfile
        sumcorr += maxcorr;
        if (write_docfile)
        {
            opt_flip = 0.;
            if (-opt_psi > 360.)
            {
                opt_psi += 360.;
                opt_flip = 1.;
            }
            dataline(0) = Iref[opt_refno].Phi(); // rot
            dataline(1) = Iref[opt_refno].Theta(); // tilt
            dataline(2) = opt_psi + 360.;        // psi
            dataline(3) = trans(0) + opt_offsets(0); // Xoff
            dataline(4) = trans(1) + opt_offsets(1); // Yoff
            dataline(5) = (double)(opt_refno + 1);   // Ref
            dataline(6) = opt_flip;              // Mirror
            dataline(7) = maxcorr;               // P_max/P_tot or Corr
	    if (do_student)
	    {
		dataline(8) = maxweight2; 
	    }
            DFo.append_comment(img.name());
            DFo.append_data_line(dataline);
        }

        // Output docfile
        if (verb > 0) if (imgno % c == 0) progress_bar(imgno);
        imgno++;
    }
    if (verb > 0) progress_bar(nn)
;
    // reverse rotation of the weighted sums
    if (maxCC_rather_than_ML || limit_trans)
        reverse_rotate_reference(Fwsum_imgs, Msum_imgs, true, wsum_Mref);
    else
        reverse_rotate_reference(Fwsum_imgs, Msum_imgs, false, wsum_Mref);

}

// Update all model parameters
void Prog_MLalign2D_prm::update_parameters(std::vector<Matrix2D<double> > &wsum_Mref,
        double &wsum_sigma_noise, double &wsum_sigma_offset,
        std::vector<double> &sumw, std::vector<double> &sumw2, std::vector<double> &sumw_mirror,
        double &sumcorr, double &sumw_allrefs)
{

    Matrix1D<double> rmean_sigma2, rmean_signal2;
    Matrix1D<int> center(2), radial_count;
    Matrix2D<std::complex<double> > Faux, Faux2;
    Matrix2D<double> Maux;
    FileName fn_tmp;
    double rr, thresh, aux, sumw_allrefs2 = 0.;
    int c;

    // Update the reference images
    sumw_allrefs = 0.;
    FOR_ALL_MODELS()
    {
        if (!do_student && sumw[refno] > 0.)
        {
            Iref[refno]() = wsum_Mref[refno];
	    Iref[refno]() /= sumw[refno];
	    Iref[refno].weight() = sumw[refno];
	    sumw_allrefs += sumw[refno];
        }
	else if (do_student && sumw2[refno] > 0.)
	{
//	    std::cerr<<" refno= "<<refno<<" sumw= "<<sumw[refno]<<" sumw2= "<<sumw2[refno]<<std::endl;
	    Iref[refno]() = wsum_Mref[refno];
	    Iref[refno]() /= sumw2[refno];
	    Iref[refno].weight() = sumw2[refno];
	    sumw_allrefs += sumw[refno];
	    sumw_allrefs2 += sumw2[refno];
 	}
        else
        {
            Iref[refno].weight() = 0.;
            Iref[refno]().initZeros(dim, dim);
        }
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

    std::cerr<<" wsum_allrefs2= "<<sumw_allrefs2<<" sumw_allrefs= "<<sumw_allrefs<<std::endl;
    // Update sigma of the origin offsets
    if (!fix_sigma_offset) 
    {
	if (do_student)
	    sigma_offset = sqrt(wsum_sigma_offset / (2. * sumw_allrefs2));
	else
	    sigma_offset = sqrt(wsum_sigma_offset / (2. * sumw_allrefs));

    }
    // Update the noise parameters
    if (!fix_sigma_noise)
    {
//      The following converges faster according to McLachlan&Peel (2000)
//	Finite Mixture MOdels, Wiley pp. 228!
	if (do_student)
	    sigma_noise = sqrt(wsum_sigma_noise / (sumw_allrefs2 * dim * dim));
	else
	    sigma_noise = sqrt(wsum_sigma_noise / (sumw_allrefs * dim * dim));
    }

    // Update annealing parameter
    //anneal*=anneal_step;
    //anneal=MAX(1.,anneal);

}

// Check convergence
bool Prog_MLalign2D_prm::check_convergence(std::vector<double> &conv)
{

    bool converged = true;
    double convv;
    Matrix2D<double> Maux;

    Maux.resize(dim, dim);
    Maux.setXmippOrigin();

    conv.clear();
    FOR_ALL_MODELS()
    {
        if (Iref[refno].weight() > 0.)
        {
            Maux = multiplyElements(Iold[refno](), Iold[refno]());
            convv = 1. / (Maux.compute_avg());
            Maux = Iold[refno]() - Iref[refno]();
            Maux = multiplyElements(Maux, Maux);
            convv *= Maux.compute_avg();
            conv.push_back(convv);
            if (convv > eps) converged = false;
        }
        else
        {
            conv.push_back(-1.);
        }
    }

    return converged;
}

// Output to screen
void Prog_MLalign2D_prm::output_to_screen(int &iter, double &sumcorr, double &LL)
{
    if (verb > 0)
    {
        if (maxCC_rather_than_ML) std::cout << "  iter " << iter << " <CC>= " + floatToString(sumcorr, 10, 5);
        else
        {
            std::cout << "  iter " << iter << " noise= " << floatToString(sigma_noise, 10, 7) << " offset= " << floatToString(sigma_offset, 10, 7);
            std::cout << "  LL= " << LL << " <Pmax/sumP>= " << sumcorr << std::endl;
            std::cout << "  Model  fraction  mirror-fraction " << std::endl;
            FOR_ALL_MODELS()
            {
                std::cout << "  " << integerToString(refno + 1, 5) << " " << floatToString(alpha_k[refno], 10, 7) << " " << floatToString(mirror_fraction[refno], 10, 7) << std::endl;
            }
        }
    }

}

void Prog_MLalign2D_prm::write_output_files(const int iter, DocFile &DFo,
        double &sumw_allrefs, double &LL, double &avecorr,
        std::vector<double> &conv)
{

    FileName          fn_tmp, fn_base, fn_tmp2;
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
    else
    {
    	if (fn_scratch!="") 
	{
	    // Clean scrath disc
	    system(((std::string)"rm -rf "+fn_scratch).c_str());
	}
    
    }

    // Write out current reference images and fill sel & log-file
    FOR_ALL_MODELS()
    {
        fn_tmp = fn_base + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "");
        fn_tmp = fn_tmp + ".xmp";
        Iref[refno].write(fn_tmp);
        SFo.insert(fn_tmp, SelLine::ACTIVE);
        fracline(0) = alpha_k[refno];
        fracline(1) = mirror_fraction[refno];
        fracline(2) = 1000 * conv[refno]; // Output 1000x the change for precision
        DFl.insert_comment(fn_tmp);
        DFl.insert_data_line(fracline);
    }

    // Write out sel & log-file
    fn_tmp = fn_base + ".sel";
    SFo.write(fn_tmp);

    DFl.go_beginning();
    comment = "MLalign2D-logfile: Number of images= " + floatToString(sumw_allrefs);
    if (maxCC_rather_than_ML) comment += " <CC>= " + floatToString(avecorr, 10, 5);
    else
    {
        comment += " LL= " + floatToString(LL, 15, 10) + " <Pmax/sumP>= " + floatToString(avecorr, 10, 5);
        DFl.insert_comment(comment);
        comment = "-noise " + floatToString(sigma_noise, 15, 12) + " -offset " + floatToString(sigma_offset, 15, 12) + " -istart " + integerToString(iter + 1);
        if (anneal > 1.) comment += " -anneal " + floatToString(anneal, 10, 7);
    }
    DFl.insert_comment(comment);
    DFl.insert_comment(cline);
    DFl.insert_comment("columns: model fraction (1); mirror fraction (2); 1000x signal change (3)");
    fn_tmp = fn_base + ".log";
    DFl.write(fn_tmp);

    if (write_docfile)
    {
        // Write out docfile with optimal transformation & references
        fn_tmp = fn_base + ".doc";
        DFo.write(fn_tmp);
    }

    if (write_selfiles)
    {
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

}

