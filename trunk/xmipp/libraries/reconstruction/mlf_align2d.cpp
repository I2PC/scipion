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
#include "mlf_align2d.h"

// Read arguments ==========================================================
void Prog_MLFalign2D_prm::read(int argc, char **argv, bool ML3D)
{

    // Generate new command line for restart procedure
    cline = "";
    if (checkParameter(argc, argv, "-restart"))
    {
        string comment;
        FileName fn_sel;
        DocFile DFi;
        DFi.read(getParameter(argc, argv, "-restart"));
        DFi.go_beginning();
        comment = (DFi.get_current_line()).get_text();
        if (strstr(comment.c_str(), "MLFalign2D-logfile") == NULL)
        {
            cerr << "Error!! Docfile is not of MLFalign2D-logfile type. " << endl;
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
            cline = cline + (string)argv[i] + " ";
        }
    }

    // Read command line
    if (checkParameter(argc, argv, "-more_options"))
    {
	usage();
	extendedUsage();
    }

    // Main parameters
    n_ref = textToInteger(getParameter(argc, argv, "-nref", "0"));
    fn_ref = getParameter(argc, argv, "-ref", "");
    fn_sel = getParameter(argc, argv, "-i");
    do_ctf_correction = !checkParameter(argc, argv, "-no_ctf");
    if (do_ctf_correction)
    {
	fn_ctfdat = getParameter(argc, argv, "-ctfdat","");
    }
    else
    {
	sampling = textToFloat(getParameter(argc, argv, "-pixel_size","1"));
    }
    fn_root = getParameter(argc, argv, "-o", "mlf2d");
    search_shift = textToInteger(getParameter(argc, argv, "-search_shift", "3"));
    psi_step = textToFloat(getParameter(argc, argv, "-psi_step", "5"));
    do_mirror = checkParameter(argc, argv, "-mirror");
    lowres_limit = textToInteger(getParameter(argc, argv, "-low", "999"));
    highres_limit = textToInteger(getParameter(argc, argv, "-high", "0"));
    ini_highres_limit = textToInteger(getParameter(argc, argv, "-ini_high", "0"));
    phase_flipped = !checkParameter(argc, argv, "-not_phase_flipped");
    reduce_snr = textToFloat(getParameter(argc, argv, "-reduce_snr", "1"));
    first_iter_noctf = checkParameter(argc, argv, "-ctf_affected_refs");

    // Less common stuff
    Niter = textToInteger(getParameter(argc, argv, "-iter", "100"));
    istart = textToInteger(getParameter(argc, argv, "-istart", "1"));
    sigma_offset = textToFloat(getParameter(argc, argv, "-offset", "3"));
    eps = textToFloat(getParameter(argc, argv, "-eps", "5e-5"));
    fn_frac = getParameter(argc, argv, "-frac", "");
    write_docfile = !checkParameter(argc, argv, "-dont_output_docfile");
    write_selfiles = !checkParameter(argc, argv, "-dont_output_selfiles");
    fix_fractions = checkParameter(argc, argv, "-fix_fractions");
    fix_sigma_offset = checkParameter(argc, argv, "-fix_sigma_offset");
    fix_sigma_noise = checkParameter(argc, argv, "-fix_sigma_noise");
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    C_fast = textToFloat(getParameter(argc, argv, "-C", "1e-12"));
    fn_doc = getParameter(argc, argv, "-doc", "");
    do_include_allfreqs=checkParameter(argc,argv,"-include_allfreqs");

    // Only for interaction with Refine3D:
    do_ML3D = ML3D;
    search_rot = textToFloat(getParameter(argc, argv, "-search_rot", "999."));

    // Hidden arguments
    do_write_offsets = checkParameter(argc, argv, "write_offsets");
    fn_scratch = getParameter(argc, argv, "-scratch", "");
    debug = textToInteger(getParameter(argc, argv, "-debug","0"));
    do_variable_psi = checkParameter(argc, argv, "-var_psi");
    do_variable_trans = checkParameter(argc, argv, "-var_trans");
}

// Show ====================================================================
void Prog_MLFalign2D_prm::show(bool ML3D)
{

    if (verb > 0)
    {
        // To screen
        if (!ML3D)
        {
            cerr << " -----------------------------------------------------------------" << endl;
            cerr << " | Read more about this program in the following publication:    |" << endl;
	    cerr << " |  Scheres ea. (2007) Structure, 15, 1167-1177                  |" << endl;
            cerr << " |                                                               |" << endl;
            cerr << " |   *** Please cite it if this program is of use to you! ***    |" << endl;
            cerr << " -----------------------------------------------------------------" << endl;
        }
        cerr << "--> Multi-reference refinement " << endl;
        cerr << "--> using a maximum-likelihood in Fourier-space (MLF) target " <<endl;
	if (do_ctf_correction)
	{
	    cerr << "--> with CTF correction "<<endl;
	}
	else
	{
	    cerr << "--> ignoring CTF effects "<<endl;
	}
        cerr << "  Input images            : " << fn_sel << " (" << nr_exp_images << ")" << endl;
        if (fn_ref != "")
	{
            cerr << "  Reference image(s)      : " << fn_ref << endl;
	}
        else
	{
            cerr << "  Number of references:   : " << n_ref << endl;
	}
	if (do_ctf_correction)
	{
	    cerr << "  CTF-parameters file     : " << fn_ctfdat << endl;
	}
        cerr << "  Output rootname         : " << fn_root << endl;
        cerr << "  Stopping criterium      : " << eps << endl;
        cerr << "  initial sigma offset    : " << sigma_offset << endl;
        cerr << "  Psi sampling interval   : " << psi_step << " degrees" << endl;
	cerr << "  Translational searches  : " << search_shift << " pixels" << endl;
	cerr << "  Low resolution limit    : " << lowres_limit << " Ang" << endl;
	cerr << "  High resolution limit   : " << highres_limit << " Ang" << endl;
	if (reduce_snr != 1.)
	    cerr << "  Multiply estimated SNR  : " << reduce_snr << endl;
	if (reduce_snr > 1.)
	{
	    cerr << "  --> WARNING!! With reduce_snr>1 you may likely overfit the noise!" << endl;
	}
        if (do_mirror)
            cerr << "  Check mirrors           : true" << endl;
        else
            cerr << "  Check mirrors           : false" << endl;
        if (fn_frac != "")
            cerr << "  Initial model fractions : " << fn_frac << endl;
	if (do_ctf_correction)
	{
	    FOR_ALL_DEFOCUS_GROUPS()
	    {
		cerr << "    + defocus group "<<ifocus+1<<" contains "<<count_defocus[ifocus]<<" images"<<endl;
	    }
	}
	if (ini_highres_limit > 0.)
	    cerr << "    + High resolution limit for 1st iteration set to " << highres_limit << "Ang"<<endl;
        if (search_rot < 180.)
            cerr << "    + Limit orientational search to +/- " << search_rot << " degrees" << endl;
	if (do_variable_psi)
            cerr << "    + Vary in-plane rotation angle sampling with resolution " << endl;
	if (do_variable_trans)
            cerr << "    + Vary in-plane translation sampling with resolution " << endl;

        // Hidden stuff
        if (fix_fractions)
        {
            cerr << "    + Do not update estimates of model fractions." << endl;
        }
        if (fix_sigma_offset)
        {
            cerr << "    + Do not update sigma-estimate of origin offsets." << endl;
        }
        if (fix_sigma_noise)
        {
            cerr << "    + Do not update estimated noise spectra." << endl;
        }
	cerr << " -----------------------------------------------------------------" << endl;

    }

}

// Fourier mode usage ==============================================================
void Prog_MLFalign2D_prm::usage()
{
    cerr << "Usage:  mlf_align2d [options] " << endl;
    cerr << "   -i <selfile>                : Selfile with all input images \n";
    cerr << "   -ctfdat <ctfdatfile>        : Two-column ASCII file with filenames and CTF parameter files of all images \n";
    cerr << "      OR -no_ctf                   OR do not use any CTF correction \n";
    cerr << "   -ref <selfile/image>        : Selfile with initial references/single reference image \n";
    cerr << "      OR -nref <int>               OR number of references to generate automatically (bias-free)\n";
    cerr << " [ -o <rootname> ]             : Output rootname (default = \"mlf2d\")\n";
    cerr << " [ -mirror ]                   : Also check mirror image of each reference \n";
    cerr << " [ -search_shift <int=3>]      : Limited translational searches (in pixels) \n";
    cerr << " [ -reduce_noise <factor=1> ]  : Use a value smaller than one to decrease the estimated SSNRs \n";
    cerr << " [ -not_phase_flipped ]        : Use this if the experimental images have not been phase flipped \n";
    cerr << " [ -ctf_affected_refs ]        : Use this if the references (-ref) are not CTF-deconvoluted \n";
    cerr << " [ -low <Ang=999> ]            : Exclude lowest frequencies from P-calculations (in Ang) \n";
    cerr << " [ -high <Ang=0> ]             : Exclude highest frequencies from P-calculations (in Ang) \n";
    cerr << " [ -ini_high <Ang=0> ]         : Exclude highest frequencies during first iteration (in Ang) \n";
    cerr << " [ -pixel_size <Ang=1> ]       : Pixel size in Angstrom (only necessary for -no_ctf mode) \n";
    cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_MLFalign2D_prm::extendedUsage(bool ML3D)
{
    cerr << "Additional options: " << endl;
    cerr << " [ -eps <float=5e-5> ]         : Stopping criterium \n";
    cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
    cerr << " [ -psi_step <float=5> ]       : In-plane rotation sampling interval [deg]\n";
    cerr << " [ -offset <float=3> ]         : Expected standard deviation for origin offset [pix]\n";
    cerr << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n";
    cerr << " [ -C <double=1e-12> ]         : Significance criterion for fast approach \n";
    if (!ML3D) cerr << " [ -restart <logfile> ]        : restart a run with all parameters as in the logfile \n";
    if (!ML3D) cerr << " [ -istart <int> ]             : number of initial iteration \n";
    cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the noise spectra \n";
    cerr << " [ -fix_sigma_offset]          : Do not re-estimate the standard deviation in the origin offsets \n";
    cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
    cerr << " [ -dont_output_docfile ]      : Do not write out docfile with most likely angles & translations \n";
    cerr << " [ -dont_output_selfiles ]     : Do not write out selfiles with most likely class assignments \n";
    cerr << " [ -doc <docfile=\"\"> ]         : Read initial angles and offsets from docfile \n";
    cerr << endl;
    exit(1);
}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produceSideInfo2)
void Prog_MLFalign2D_prm::produceSideInfo()
{

    FileName                    fn_img, fn_tmp, fn_base, fn_tmp2;
    ImageXmipp                  img;
    FourierImageXmipp           fourimg;
    XmippCTF                    ctf;
    SelLine                     SL;
    SelFile                     SFtmp, SFpart;
    Matrix1D<double>            offsets(2), dum, rmean_ctf;
    Matrix2D<double>            A(3, 3), Maux, Maux2;
    Matrix2D<complex<double> >  Faux, ctfmask;
    Matrix1D<int>               center(2), radial_count;
    vector<int>                 tmppointp, tmppointp_nolow, tmppointi, tmppointj;
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
        first_iter_noctf = true;
    }

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
        REPORT_ERROR(1, "Prog_MLFalign2D_prm: Provide a selfile with unique image names (preferably all in one directory)");

    psi_max = 90.;
    nr_psi = CEIL(psi_max / psi_step);
    max_nr_psi = nr_psi;
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

    // Set limit_rot & limit_trans
    if (search_rot < 180.) limit_rot = true;
    else limit_rot = false;
    

    // Fill limited translation shift-vectors
    nr_trans = 0;
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Maux)
    {
	double rr = (double)i * i + (double)j * j;
	if (rr <= (double)search_shift*search_shift)
	{
	    offsets(0) = (double)j;
	    offsets(1) = (double)i;
	    Vtrans.push_back(offsets);
	    if (i == 0 && j == 0) zero_trans = nr_trans;
	    nr_trans++;
	}
    }

    FileName fnt_img, fnt_ctf, fnt;
    vector<FileName> all_fn_ctfs;
    bool found, is_unique;
    int iifocus;

    count_defocus.clear();
    Vctf.clear();
    Vdec.clear();
    Vsig.clear();
    if (!do_ctf_correction)
    {
	nr_focus = 1;
	count_defocus.push_back(nr_exp_images);
	Vctf.push_back(dum);
	Vdec.push_back(dum);
	Vsig.push_back(dum);
	Vctf[0].resize(hdim);
	Vsig[0].resize(hdim);
	Vdec[0].resize(hdim);
	for (int irr = 0; irr < hdim; irr++)
	{
	    dVi(Vctf[0], irr) = 1.;
	    dVi(Vdec[0], irr) = 1.;
	}
    }
    else
    {
	// Read ctfdat and determine the number of CTF groups
	nr_focus = 0;
	all_fn_ctfs.clear();
	ctfdat.read(fn_ctfdat);
	ctfdat.goFirstLine();
	while (!ctfdat.eof())
	{
	    ctfdat.getCurrentLine(fnt_img,fnt_ctf);
	    is_unique=true;
	    for (int i = 0; i< all_fn_ctfs.size(); i++)
	    {
		if (fnt_ctf == all_fn_ctfs[i])
		{
		    is_unique=false;
		    break;
		}
	    }
	    if (is_unique)
	    {
		all_fn_ctfs.push_back(fnt_ctf);
		count_defocus.push_back(0);
		nr_focus++;
	    }
	    ctfdat.nextLine();
	}
	
	// Fill count_defocus vectors and add CTF-group labels to selfile
	SF.go_beginning();
	SFtmp.clear();
	while (!SF.eof())
	{
	    SL=SF.current();
	    fnt=SF.NextImg();
	    // Find which CTF group it belongs to
	    found = false;
	    ctfdat.goFirstLine();
	    while (!ctfdat.eof())
	    {
		ctfdat.getCurrentLine(fnt_img,fnt_ctf);
		if (fnt_img == fnt)
		{
		    found = true;
		    for (iifocus = 0; iifocus< all_fn_ctfs.size(); iifocus++)
			if (fnt_ctf == all_fn_ctfs[iifocus])
			    break;
		    break;
		}
		ctfdat.nextLine();
	    }
	    if (!found)
		REPORT_ERROR(1, "Prog_MLFalign2D_prm: Did not find image "+fnt+" in the CTFdat file");
	    count_defocus[iifocus]++;
	    SL.set_number(iifocus + 1);
	    SFtmp.insert(SL);
	}
	SF = SFtmp;
	
	// Check the number of images in each defocus group
	// and read CTF-parameters from disc
	FOR_ALL_DEFOCUS_GROUPS()
	{
	    if (count_defocus[ifocus] < 50 && verb > 0)
		cerr << "WARNING%% CTF group " << ifocus + 1 << " contains less than 50 images!" << endl;
	    
	    Vctf.push_back(dum);
	    Vdec.push_back(dum);
	    Vsig.push_back(dum);
	    ctf.read(all_fn_ctfs[ifocus]);
	    if (ABS(ctf.DeltafV - ctf.DeltafU) >1.) 
	    {
		REPORT_ERROR(1, "Prog_MLFalign2D-ERROR%% Only non-astigmatic CTFs are allowed!");
	    }
	    ctf.K = 1.;
	    ctf.enable_CTF = true;
	    ctf.Produce_Side_Info();
	    ctf.Generate_CTF(dim, dim, ctfmask);
	    if (ifocus == 0)
	    {
		sampling = ctf.Tm;
		Q0 = ctf.Q0;
	    }
	    else
	    {
		if (sampling != ctf.Tm)
		    REPORT_ERROR(1, "Prog_MLFalign2D-ERROR%% Different sampling rates in CTF parameter files!");
		if (Q0 != ctf.Q0 )
		    REPORT_ERROR(1, "Prog_MLFalign2D-ERROR%% Avoid different Q0 values in the CTF parameter files!");
	    }
	    Maux.resize(dim, dim);
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Maux)
	    {
		if (phase_flipped) dMij(Maux, i, j) = fabs(dMij(ctfmask, i, j).real());
		else dMij(Maux, i, j) = dMij(ctfmask, i, j).real();
	    }
	    CenterFFT(Maux, true);
	    center.initZeros();
	    rmean_ctf.initZeros();
	    Maux.setXmippOrigin();
	    radialAverage(Maux, center, rmean_ctf, radial_count, true);
	    Vctf[ifocus].resize(hdim);
	    for (int irr = 0; irr < hdim; irr++)
	    {
		dVi(Vctf[ifocus], irr) = dVi(rmean_ctf, irr);
	    }
	    Vdec[ifocus].resize(Vctf[ifocus]);
	    Vsig[ifocus].resize(Vctf[ifocus]);
	}
    }

    // Get a resolution pointer in Fourier-space
    Mresol.resize(dim, dim);
    Mresol.setXmippOrigin();
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Mresol)
    {
	MAT_ELEM(Maux, i, j) = XMIPP_MIN((double)hdim, sqrt((double)(i * i + j * j)));
	if (i==0 && j==0) 
	{
	    MAT_ELEM(Mresol, i, j) = 99999.;
	}
	else
	{
	    MAT_ELEM(Mresol, i, j) = (sampling * dim) / sqrt((double)(i * i + j * j));
	}
    }
    CenterFFT(Mresol, false);
    CenterFFT(Maux, false);
    Mresol_int.resize(dim, dim);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mresol_int)
    {
	dMij(Mresol_int, i, j) = ROUND(dMij(Maux, i, j));
    }
    
}

// For initial noise variances
void Prog_MLFalign2D_prm::estimateInitialNoiseSpectra()
{

    Matrix2D<double>            Maux, Mallave;
    Matrix2D<complex<double> >  Fimg, Faux, Fave;
    Matrix1D<double>            rmean_noise, rmean_signal, rmean_avesignal;
    vector<Matrix2D<double> >   Msigma2, Mave;
    Matrix1D<int>               center(2), radial_count;
    double                      rr;
    ImageXmipp                  img;
    FileName                    fn_tmp;
    SelLine                     SL;
    int                         focus, c, nn;
    ofstream                    fh;

    center.initZeros();

    // For first iteration only: calculate sigma2 (& spectral noise) from power
    // spectra of all images and subtract power spectrum of average image
    // (to take away low-res frequencies where the signal dominates!)
    if (istart == 1)
    {

        if (verb > 0) cerr << "--> Estimating initial noise models from average power spectra ..." << endl;
        if (verb > 0)
        {
            nn = SF.ImgNo();
            init_progress_bar(nn);
            c = XMIPP_MAX(1, nn / 60);
        }
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
        SF.go_beginning();
        int imgno = 0;
        while (!SF.eof())
        {
            SL = SF.current();
            focus = SL.get_number() - 1;
            img.read(SF.NextImg(), false, false, false, false);
            img().setXmippOrigin();
            FourierTransform(img(), Fimg);
            FFT_magnitude(Fimg, Maux);
            Maux.setXmippOrigin();
            Maux *= Maux;
            Msigma2[focus] += Maux;
            Mave[focus] += img();
            imgno++;
            if (verb > 0) if (imgno % c == 0) progress_bar(imgno);
        }
        if (verb > 0) progress_bar(nn);

        // Calculate Vsig vectors and write them to disc
        FOR_ALL_DEFOCUS_GROUPS()
        {
            Mave[ifocus] /= count_defocus[ifocus];
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
            rmean_noise /= count_defocus[ifocus];
            // Subtract signal terms
            // Divide by factor 2 because of the 2D-Gaussian distribution!
            for (int irr = 0; irr < hdim; irr++)
            {
                dVi(Vsig[ifocus], irr) = (dVi(rmean_noise, irr) - dVi(rmean_signal, irr)) / 2.;
            }

            // write Vsig vector to disc
            fn_tmp = fn_root + "_it";
            fn_tmp.compose(fn_tmp, istart - 1, "");
            fn_tmp += "_sig";
            if (nr_focus > 1) fn_tmp.compose(fn_tmp, ifocus + 1, "");
            fn_tmp += ".dat";
            fh.open((fn_tmp).c_str(), ios::out);
            if (!fh) REPORT_ERROR(1, (string)"Prog_MLFalign2D_prm: Cannot write file: " + fn_tmp);
	    for (int irr = 0; irr < hdim; irr++)
            {
                fh << (double)irr/(sampling*dim) << " " << dVi(Vsig[ifocus], irr) << "\n";
            }
            fh.close();
        }
        Msigma2.clear();
        Mave.clear();
    }



}

void Prog_MLFalign2D_prm::updateWienerFilters(Matrix1D<double> &spectral_signal,
					      int iter)
{

    // Use formula 2.32b on p60 from Frank's book 2nd ed.,
    // Assume that Vctf, Vdec and Vsig exist (with their right sizes)
    // and that Vctf, Vsig are already filled with the right values
    vector<Matrix1D<double> >  Vsnr;
    Matrix1D<double>           Vzero, Vdenom, Vavgctf2;
    Matrix2D<double>           Maux;
    Matrix2D<complex<double> > Faux;
    ofstream                   fh;
    int                        maxres = 0;
    double                     noise, ssnr, current_probres_limit, current_highres_limit;
    FileName                   fn_base, fn_tmp;

    // Pre-calculate average CTF^2 and initialize Vsnr
    Vavgctf2.initZeros(hdim);
    Vzero.initZeros(hdim);
    FOR_ALL_DEFOCUS_GROUPS()
    {
        Vsnr.push_back(Vzero);
        for (int irr = 0; irr < hdim; irr++)
        {
            dVi(Vavgctf2, irr) += dVi(Vctf[ifocus], irr) * dVi(Vctf[ifocus], irr) * count_defocus[ifocus];
        }
    }
    Vavgctf2 /= nr_exp_images;

    // Calculate SSNR for all defocus groups
    // For each group the spectral noise is estimated via (2*Vsig)/(count_defocus-1)
    // The spectral signal is not split into defocus groups
    // Therefore, affect the average spectral_signal for each defocu
    // group with its CTF^2 and divide by the average CTF^2
    FOR_ALL_DEFOCUS_GROUPS()
    {
	for (int irr = 0; irr < hdim; irr++)
	{
	    noise = 2. * dVi(Vsig[ifocus], irr) * (double)count_defocus[ifocus];
	    noise /= (double)(count_defocus[ifocus] - 1);
	    if (noise < 1e-20) dVi(Vsnr[ifocus], irr) = 0.;
	    else
	    {
		if (dVi(Vavgctf2, irr) > 0.)
		{
		    dVi(Vsnr[ifocus], irr) = reduce_snr * dVi(Vctf[ifocus], irr) * dVi(Vctf[ifocus], irr) *
			dVi(spectral_signal, irr) / (dVi(Vavgctf2, irr) * noise);
		}
		else dVi(Vsnr[ifocus], irr) = 0.;
	    }
	    // For start from already CTF-deconvoluted references:
	    if ((iter == istart - 1) && !first_iter_noctf)
	    {
		dVi(Vsnr[ifocus], irr) *= dVi(Vavgctf2, irr);
	    }
	    // Take ini_highres_limit into account
	    if ( iter == istart - 1 && ini_highres_limit > 0. 
		&& sampling*dim/(double)irr < ini_highres_limit )
	    {
		dVi(Vsnr[ifocus], irr) = 0.;
	    }
	    // Subtract 1 according Unser et al.
	    dVi(Vsnr[ifocus], irr) = XMIPP_MAX(0., dVi(Vsnr[ifocus], irr) - 1.);
	    // Prevent spurious high-frequency significant SNRs from random averages
	    if (iter == 0 && do_generate_refs)
	    {
		dVi(Vsnr[ifocus], irr) = XMIPP_MAX(0., dVi(Vsnr[ifocus], irr) - 2.);
	    }
	    if (dVi(Vsnr[ifocus], irr) > 0. && irr > maxres) 
	    {
		maxres = irr;
	    }
	}
    }
    
    // Check that at least some frequencies have non-zero SSNR...
    if (maxres == 0)
	REPORT_ERROR(1, "Prog_MLFalign2D_prm: All frequencies have zero spectral SNRs... (increase -reduce_snr) ");

    if (do_ctf_correction)
    {
	// Pre-calculate denominator of eq 2.32b of Frank's book (2nd ed.)
	Vdenom.initZeros(hdim);
	for (int irr = 0; irr < hdim; irr++)
	{
	    FOR_ALL_DEFOCUS_GROUPS()
	    {
		dVi(Vdenom, irr) += count_defocus[ifocus] * dVi(Vsnr[ifocus], irr) *
		    dVi(Vctf[ifocus], irr) * dVi(Vctf[ifocus], irr);
	    }
	    dVi(Vdenom, irr) += 1.;
	    dVi(Vdenom, irr) /= nr_exp_images;
	}

	// Calculate Wiener filters
	FOR_ALL_DEFOCUS_GROUPS()
	{
	    for (int irr = 0; irr < hdim; irr++)
	    {
		if (dVi(Vsnr[ifocus], irr) > 0.)
		{
		    dVi(Vdec[ifocus], irr) = dVi(Vsnr[ifocus], irr) * dVi(Vctf[ifocus], irr) / dVi(Vdenom, irr);
		    // Prevent too strong Wiener filter artefacts
		    dVi(Vdec[ifocus], irr) = XMIPP_MIN(10.,dVi(Vdec[ifocus], irr) );
		}
		else
		{
		    dVi(Vdec[ifocus], irr) = 0.;
		}
	    }
	}
    }

    // Write Wiener filters and spectral SNR to text files
    if (verb > 0)
    {
        fn_base = fn_root;
        fn_base += "_it";
        fn_base.compose(fn_base, iter, "");
        // Defocus group-specific Wiener filter files
        FOR_ALL_DEFOCUS_GROUPS()
        {
            fn_tmp = fn_base + "_wien";
            if (nr_focus > 1) fn_tmp.compose(fn_tmp, ifocus + 1, "");
            fn_tmp += ".txt";
            fh.open((fn_tmp).c_str(), ios::out);
            if (!fh)
                REPORT_ERROR(3008, (string)"Prog_MLFalign2D_prm: Cannot write file: " + fn_tmp);
            fh  << "#  Resol      SSNR       CTF    Wiener    signal     noise       Ang" << endl;
            for (int irr = 0; irr < hdim; irr++)
            {
                fh << floatToString((double)irr / (sampling*dim));
                fh.width(10);
                fh << floatToString(XMIPP_MIN(25., dVi(Vsnr[ifocus], irr)));
                fh.width(10);
                fh << floatToString(dVi(Vctf[ifocus], irr));
                fh.width(10);
                fh << floatToString(dVi(Vdec[ifocus], irr));
                fh.width(10);
                fh << floatToString(dVi(spectral_signal, irr));
                fh.width(10);
                noise = 2. * dVi(Vsig[ifocus], irr) * (double)count_defocus[ifocus];
                noise /= (double)(count_defocus[ifocus] - 1);
                fh << floatToString(noise);
		fh.width(10);
                if (irr>0) 
		{
		    fh << floatToString((sampling*dim)/(double)irr);
		}
		else
		{
		    fh << floatToString(999.);
		}
		fh << endl;
            }
            fh.close();
        }
    }

    // Set the current resolution limits
    if (do_include_allfreqs) 
    {
	maxres = hdim;
	current_highres_shell = maxres;
	current_probres_limit = sampling * dim / maxres;
	current_highres_limit = sampling * dim / maxres;
    }
    else
    {
	maxres = XMIPP_MIN(maxres, hdim);
	current_probres_limit = sampling * dim / maxres;
	maxres += 5; // hard-code increase_highres_limit to 5
	maxres = XMIPP_MIN(maxres, hdim);
	current_highres_shell = maxres;
 	current_highres_limit = sampling * dim / maxres;
    }
    current_probres_limit = XMIPP_MAX(current_probres_limit, highres_limit);
    current_highres_limit = XMIPP_MAX(current_highres_limit, highres_limit);

    if (debug) cerr<<"resolution limits: low= "<<lowres_limit<<" prob= "<<current_probres_limit
		   <<" high= "<<current_highres_limit<<" overall high= "<<highres_limit<<" Ang"<<endl;

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
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
    {
        double resol = dMij(Mresol, i, j);
	if (resol <= lowres_limit &&
	    resol >= current_probres_limit &&
	    !(i == 0 && j > hdim) ) // exclude first half row in FourierTransformHalf
	{
            pointer_2d.push_back(i*XSIZE(Maux) + j);
            pointer_i.push_back(i);
            pointer_j.push_back(j);
            nr_points_prob++;
	}
    }
    // Second, get the rest of the currently relevant pixels:
    // These are within <current_probres_limit,current_highres_limit]
    nr_points_2d = nr_points_prob;
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
    {
        double resol = dMij(Mresol, i, j);
	if (resol < current_probres_limit &&
	    resol >= current_highres_limit &&
	    !(i == 0 && j > hdim) ) // exclude first half row in FourierTransformHalf
	{
            pointer_2d.push_back(i*XSIZE(Maux) + j);
            pointer_i.push_back(i);
            pointer_j.push_back(j);
            nr_points_2d++;
	}
    }
    dnr_points_2d = 2 * nr_points_2d;
    if (debug) cerr<<"nr_points_2d= "<<nr_points_2d<<" nr_points_prob= "<<nr_points_prob<<endl;

    setCurrentSamplingRates(current_probres_limit);
}


// Vary psi_step and trans_step with resolution =============================================
void Prog_MLFalign2D_prm::setCurrentSamplingRates(double current_probres_limit)
{

    if (do_variable_psi)
    {
	// Let's try a resolution-dependent in-plane rotational sampling
	// 2x oversample for current_probres_limit
	// Use user-provided psi_step as minimum sampling
	// SIND(0.5*psi_step) = D / dim;
	psi_step = ASIND(current_probres_limit/(dim*sampling));
	nr_psi = CEIL(psi_max / psi_step);
	if (verb>0)
	{
	    cerr<<" Current resolution= "<<current_probres_limit<<" Ang; current psi_step = "<<psi_step<<endl;
	}    
	nr_psi = XMIPP_MIN(nr_psi, max_nr_psi);
	psi_step = psi_max / nr_psi;
    }

    if (do_variable_trans)
    {
	Matrix1D<double>  offsets(2);   
	nr_trans = 0;
	Vtrans.clear();
	// Lets sample the in-plane translations 4x the current resolution
	int trans_step = ROUND(current_probres_limit/4.*sampling);
	if (verb>0)
	{
	    cerr<<" Current resolution= "<<current_probres_limit<<" Ang; current trans_step = "<<trans_step<<endl;
	}    
	trans_step = XMIPP_MAX(1,trans_step);
	for (int ix = -search_shift*trans_step; ix <= search_shift*trans_step; ix+=trans_step)
	{
	    for (int iy = -search_shift*trans_step; iy <= search_shift*trans_step; iy+=trans_step)
	    {
		double rr = (double)ix * ix + (double)iy * iy;
		if (rr <= (double)trans_step*trans_step*search_shift*search_shift)
		{
		    offsets(0) = ix;
		    offsets(1) = iy;
		    cerr<<"trans: "<<ix<<" "<<iy<<endl;
		    Vtrans.push_back(offsets);
		    if (ix == 0 && iy == 0) zero_trans = nr_trans;
		    nr_trans++;
		}
	    }
	}
    }

}

// Generate initial references =============================================
void Prog_MLFalign2D_prm::generateInitialReferences()
{

    SelFile SFtmp, SFout;
    ImageXmipp Iave, Itmp;
    double dummy;
    FileName fn_tmp;
    SelLine line;

    if (verb > 0)
    {
        cerr << "  Generating initial references by averaging over random subsets" << endl;
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
void Prog_MLFalign2D_prm::produceSideInfo2(int nr_vols)
{

    int                       c, idum, refno = 0;
    DocFile                   DF, DF2;
    DocLine                   DL;
    double                    offx, offy, aux, sumfrac = 0.;
    FileName                  fn_tmp;
    ImageXmipp                img;
    vector<double>            Vdum;

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
        Iold.push_back(img);
        Ictf.push_back(img);
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
    	system(((string)"rm -rf "+fn_scratch).c_str());
 	// Generate new one
    	system(((string)"mkdir -p " + fn_scratch).c_str());
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
                REPORT_ERROR(1, (string)"Prog_MLFalign2D_prm: Cannot find " + fn_tmp + " in docfile " + fn_doc);
            }
        }
        DF.go_beginning();
    }

    // If we don't write offsets to disc, fill imgs_offsets vectors
    // with zeros
    if (!do_write_offsets)
    {
        if (do_mirror) idum = 4 * n_ref;
        else idum = 2 * n_ref;
        for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
        {
            imgs_offsets.push_back(Vdum);
            for (int refno = 0; refno < idum; refno++)
            {
                imgs_offsets[imgno].push_back(0.);
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
                    REPORT_ERROR(1, (string)"Prog_MLFalign2D_prm: Cannot find " + fn_tmp + " in docfile " + fn_doc);
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
                    REPORT_ERROR(1, "Prog_MLFalign2D_prm: Mirror fraction (2nd column) should be [0,1]!");
                mirror_fraction[refno] = DL[1];
            }
            sumfrac += alpha_k[refno];
            DF.next_data_line();
        }
        if (ABS(sumfrac - 1.) > 1e-3)
            if (verb > 0) cerr << " ->WARNING: Sum of all expected model fractions (" << sumfrac << ") is not one!" << endl;
        for (refno = 0; refno < n_ref; refno++)
        {
            alpha_k[refno] /= sumfrac;
        }
    }


    Matrix2D<complex<double> >   Faux;
    Matrix2D<double>             Maux;
    Matrix1D<double>             dum, rmean_signal2, spectral_signal;
    Matrix1D<int>                center(2), radial_count;
    ifstream                     fh;

    center.initZeros();
    // Calculate average spectral signal
    c = 0;
    FOR_ALL_MODELS()
    {
	if (alpha_k[refno] > 0.)
	{
	    FourierTransform(Iref[refno](), Faux);
	    FFT_magnitude(Faux, Maux);
	    CenterFFT(Maux, true);
	    Maux *= Maux;
	    Maux *= alpha_k[refno] * (double)nr_exp_images;
	    Maux.setXmippOrigin();
	    rmean_signal2.initZeros();
	    radialAverage(Maux, center, rmean_signal2, radial_count, true);
	    if (c == 0) spectral_signal = rmean_signal2;
	    else spectral_signal += rmean_signal2;
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
	spectral_signal /= (double)nr_vols * 2.;
    }
    else
    {
	// Divide by the number of reference images
	spectral_signal /= (double)n_ref;
    }

    // Read in Vsig-vectors with fixed file names
    FOR_ALL_DEFOCUS_GROUPS()
    {
	fn_tmp = fn_root + "_it";
	fn_tmp.compose(fn_tmp, istart - 1, "");
	fn_tmp += "_sig";
	if (nr_focus > 1) fn_tmp.compose(fn_tmp, ifocus + 1, "");
	fn_tmp += ".dat";
	fh.open((fn_tmp).c_str(), ios::in);
	if (!fh) REPORT_ERROR(1, (string)"Prog_MLFalign2D_prm: Cannot read file: " + fn_tmp);
	else
	{
	    for (int irr = 0; irr < hdim; irr++)
	    {
		fh >> aux;
		if (ABS(aux - ((double)irr/(sampling*dim)) ) > 0.01 ) 
		{
		    cerr<<"aux= "<<aux<<" resol= "<<(double)irr/(sampling*dim)<<endl;
		    REPORT_ERROR(1, (string)"Prog_MLFalign2D_prm: Wrong format: " + fn_tmp);
		}
		fh >> aux;
		dVi(Vsig[ifocus], irr) = aux;
	    }
	}
	fh.close();
    }
    
    // Calculate all Wiener filters
    updateWienerFilters(spectral_signal, istart - 1);

}

void Prog_MLFalign2D_prm::writeOffsets(FileName fn, vector<double> &data)
{

    ofstream fh;
    int itot;

    if (fn_scratch!="")
    {
    	// Write to scratch disc
    	fn = fn_scratch + "/" + fn;
    }

    fh.open((fn).c_str(), ios::out);
    if (!fh)
    {
        fh.clear();
        // Create the directory if it does not exist yet, and try again
        string dirname;
        int last_slash = ((string)fn).rfind("/");
        dirname = ((string)fn).erase(last_slash);
        if (!exists(dirname)) system(((string)"mkdir -p " + dirname).c_str());
        fh.open((fn).c_str(), ios::out);
        if (!fh)
            REPORT_ERROR(3008, (string)"Prog_MLFalign2D_prm: Cannot write file: " + fn);
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

bool Prog_MLFalign2D_prm::readOffsets(FileName fn, vector<double> &data)
{

    ifstream fh;
    int ii, itot, nr_off, itoth, nr_offh;
    double remain;
    vector<double> data1;

    if (fn_scratch!="")
    {
    	// Read from scratch disc
    	fn = fn_scratch + "/" + fn;
    }

    if (!exists(fn)) return false;
    else
    {
        fh.open((fn).c_str(), ios::in);
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
void Prog_MLFalign2D_prm::calculateInPlanePDF()
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

void Prog_MLFalign2D_prm::appendFTtoVector(const Matrix2D<complex<double> > &in,
					   vector<double> &out)
{

    // First, store the points used in the probability calculations
    for (int ipoint = 0; ipoint < nr_points_2d; ipoint++)
    {
	int ii = pointer_2d[ipoint];
	out.push_back((in.data[ii]).real());
	out.push_back((in.data[ii]).imag());
    }

}

void Prog_MLFalign2D_prm::getFTfromVector(const vector<double> &in, 
					  const int start_point,
					  Matrix2D<complex<double> > &out, 
					  bool only_real)
{

    out.resize(hdim + 1, dim);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(out)
    {
	dMij(out,i,j) = 0.;
    }
    if (only_real)
    {
	for (int ipoint = 0; ipoint < nr_points_2d; ipoint++)
	{
	    int ii = pointer_2d[ipoint];
	    complex<double> aux(in[start_point+ipoint],0.);
	    out.data[ii] = aux;
	}
    }
    else
    {
	for (int ipoint = 0; ipoint < nr_points_2d; ipoint++)
	{
	    int ii = pointer_2d[ipoint];
	    complex<double> aux(in[start_point+2*ipoint],in[start_point+2*ipoint+1]);
	    out.data[ii] = aux;
	}
    }
}

// Rotate reference for all models and rotations and fill Fref vectors =============
void Prog_MLFalign2D_prm::rotateReference(vector<double> &out)
{
    out.clear();

    double AA, stdAA, psi;
    Matrix2D<double> Maux;
    Matrix2D<complex<double> > Faux;

    Maux.initZeros(dim, dim);
    Maux.setXmippOrigin();

    FOR_ALL_MODELS()
    {
        FOR_ALL_ROTATIONS()
	{
	    // Add arbitrary number (small_angle) to avoid 0-degree rotation (lacking interpolation)
	    psi = (double)(ipsi * psi_max / nr_psi) + SMALLANGLE;
	    Maux = Iref[refno]().rotateBSpline(3, psi, WRAP);
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
		Faux *= sqrt(stdAA / AA);
	    }
	    // Add all points as doubles to the vector
	    appendFTtoVector(Faux,out);
	}

	// Free memory
	Iref[refno]().resize(0,0);
    }
}


// Collect all rotations and sum to update Iref() for all models ==========
void Prog_MLFalign2D_prm::reverseRotateReference(const vector<double> &in,
						 vector<Matrix2D<double > > &out)
{

    double psi, dum, avg, ang;
    Matrix2D<double> Maux, Maux2;
    Matrix2D<complex<double> > Faux;
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
	    Maux2 = Maux.rotateBSpline(3, -psi, WRAP);
            out[refno] += Maux2;
        }
    }

}

void Prog_MLFalign2D_prm::preselectDirections(float &phi, float &theta,
        vector<double> &pdf_directions)
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
            u.normalize();
            v.normalize();
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

void Prog_MLFalign2D_prm::fourierTranslate2D(const vector<double> &in,
					     Matrix1D<double> &trans,
					     vector<double> &out,
					     int point_start)
{
    double xx, yy, xxshift, yyshift, dotp;
    double a, b, c, d, ac, bd, ab_cd;
    xxshift = -trans(0) / (double)dim;
    yyshift = -trans(1) / (double)dim;

    for (int i = 0; i < nr_points_2d; i++)
    {
        xx = (double)pointer_j[i];
        yy = (double)pointer_i[i];
	dotp = 2 * PI * (xx * xxshift + yyshift * yy);
        a = cos(dotp);
        b = sin(dotp);
	c = in[point_start + 2*i];
        d = in[point_start + 2*i+1];
        ac = a * c;
        bd = b * d;
        ab_cd = (a + b) * (c + d);
	out.push_back(ac - bd); // real
	out.push_back(ab_cd - ac - bd); // imag
    }

}

void Prog_MLFalign2D_prm::calculateFourierOffsets(const Matrix2D<double> &Mimg,
						  const vector<double > &offsets,
						  vector<double>  &out,
						  Matrix2D<int> &Moffsets, 
						  Matrix2D<int> &Moffsets_mirror)
{

    int irefmir, ix, iy, opt_ix, opt_iy, iflip, nr_mir, iflip_start, iflip_stop, count;
    double dxx, dyy;
    vector<double> Fimg_flip;
    Matrix2D<complex<double> > Fimg;
    Matrix1D<double> trans(2);
    Matrix2D<double> Maux2, Maux;

    Moffsets.resize(dim, dim);
    Moffsets.setXmippOrigin();
    Moffsets.init_constant(-1);
    Moffsets_mirror.resize(dim, dim);
    Moffsets_mirror.setXmippOrigin();
    Moffsets_mirror.init_constant(-1);
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();
    Maux2.resize(dim, dim);
    Maux2.setXmippOrigin();

    if (do_mirror) nr_mir = 2;
    else nr_mir = 1;

    // Flip images and store the precalculates Fourier Transforms in Fimg_flip
    out.clear();
    FOR_ALL_FLIPS()
    {
        Maux.setXmippOrigin();
        applyGeometry(Maux, F[iflip], Mimg, IS_INV, WRAP);
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
            irefmir = imir * n_ref + refno;
            iflip_start = imir * nr_nomirror_flips;
            iflip_stop = imir * nr_nomirror_flips + nr_nomirror_flips;
            FOR_ALL_LIMITED_TRANSLATIONS()
            {
                ix = ROUND(offsets[2*irefmir] + Vtrans[itrans](0));
                iy = ROUND(offsets[2*irefmir+1] + Vtrans[itrans](1));
                dxx = (double)intWRAP(ix, Moffsets.startingX(), Moffsets.finishingX());
                dyy = (double)intWRAP(iy, Moffsets.startingY(), Moffsets.finishingY());
                // For non-mirrors
                if (imir == 0 && MAT_ELEM(Moffsets, ROUND(dyy), ROUND(dxx)) < 0)
                {
                    for (int iflip = iflip_start; iflip < iflip_stop; iflip++)
                    {
                        trans(0) = dxx * DIRECT_MAT_ELEM(F[iflip], 0, 0) + dyy * DIRECT_MAT_ELEM(F[iflip], 0, 1);
                        trans(1) = dxx * DIRECT_MAT_ELEM(F[iflip], 1, 0) + dyy * DIRECT_MAT_ELEM(F[iflip], 1, 1);
                        fourierTranslate2D(Fimg_flip,trans,out,iflip*dnr_points_2d);
                    }
                    MAT_ELEM(Moffsets, ROUND(dyy), ROUND(dxx)) = count;
                    count++;
                }
                // For mirrors use a separate offset-matrix
                else if (imir == 1 && MAT_ELEM(Moffsets_mirror, ROUND(dyy), ROUND(dxx)) < 0)
                {
                    for (int iflip = iflip_start; iflip < iflip_stop; iflip++)
                    {
                        trans(0) = dxx * DIRECT_MAT_ELEM(F[iflip], 0, 0) + dyy * DIRECT_MAT_ELEM(F[iflip], 0, 1);
                        trans(1) = dxx * DIRECT_MAT_ELEM(F[iflip], 1, 0) + dyy * DIRECT_MAT_ELEM(F[iflip], 1, 1);
                        fourierTranslate2D(Fimg_flip,trans,out,iflip*dnr_points_2d);
                    }
                    MAT_ELEM(Moffsets_mirror, ROUND(dyy), ROUND(dxx)) = count;
                    count++;
                }
            }
        }
    }
}


// Exclude translations from the MLF_integration
// For significantly contributing refno+psi: re-calculate optimal shifts
void Prog_MLFalign2D_prm::processOneImage(const Matrix2D<double> &Mimg, 
					  const int focus, bool apply_ctf,
					  const vector<double> &Fref,
					  vector<double> &Fwsum_imgs,
					  vector<double> &Fwsum_ctfimgs,
					  vector<double> &Mwsum_sigma2,
					  double &wsum_sigma_offset, 
					  vector<double> &sumw,
					  vector<double> &sumw_mirror,
					  double &LL, double &fracweight, 
					  int &opt_refno, double &opt_psi,
					  Matrix1D<double> &opt_offsets, 
					  vector<double> &opt_offsets_ref,
					  vector<double > &pdf_directions)
{

    Matrix3D<double>                             Mweight;
    Matrix2D<int>                                Moffsets, Moffsets_mirror;
    vector<double>                               Fimg_trans;
    vector<Matrix1D<double> >                    uniq_offsets;

    vector<double> refw(n_ref), refw_mirror(n_ref), Pmax_refmir(2*n_ref);
    vector<double> sigma2, ctf, decctf;
    double sigma_noise2, XiA, Xi2, aux, fracpdf, pdf, weight;
    double tmpr, tmpi, sum_refw = 0.;
    double diff, maxweight = -99.e99, mindiff2 = 99.e99;
    double logsigma2;
    int    irot, opt_ipsi, opt_iflip, opt_irot, irefmir, opt_irefmir, ix, iy, ii;
    int    nr_uniq, point_trans, zscore;
    int    opt_itrans, iflip_start, iflip_stop, nr_mir;

    TimeStamp t0; 
    time_config();

    annotate_time(&t0);

    // Convert 1D Vsig vectors to the correct vector size for the 2D FTHalfs
    sigma2.clear();
    ctf.clear();
    decctf.clear();
    for (int ipoint = 0; ipoint < nr_points_2d; ipoint++)
    {
	int ii = pointer_2d[ipoint];
	int ires = Mresol_int.data[ii];
	sigma2.push_back(2.* Vsig[focus].data[ires]);
	decctf.push_back(Vdec[focus].data[ires]);
	ctf.push_back(Vctf[focus].data[ires]);
    }
    if (!apply_ctf)
    {
	for (int ipoint = 0; ipoint < nr_points_2d; ipoint++)
	{
	    ctf[ipoint] = 1.;
	}
    }

    // Precalculate normalization constant
    logsigma2 = 0.;
    for (int ipoint = 0; ipoint < nr_points_prob; ipoint++)
    {
	logsigma2 += log(PI*sigma2[ipoint]);
    }

    // Precalculate Fimg_trans, on pruned and expanded offset list
    calculateFourierOffsets(Mimg, opt_offsets_ref, Fimg_trans, Moffsets, Moffsets_mirror);

    if (debug==1) {cout<<"processOneImage 1 "; print_elapsed_time(t0); annotate_time(&t0);}

    Mweight.initZeros(nr_trans, n_ref, nr_flip*nr_psi);
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                ix = ROUND(opt_offsets_ref[2*irefmir]);
                iy = ROUND(opt_offsets_ref[2*irefmir+1]);
                Pmax_refmir[irefmir] = 0.;
                if (iflip < nr_nomirror_flips) 
		{ 
		    point_trans = MAT_ELEM(Moffsets, iy, ix);
                }
		else 
		{
		    point_trans = MAT_ELEM(Moffsets_mirror, iy, ix);
		}
		if (point_trans < 0 || point_trans > dim2)
		{
		    cerr<<"point_trans = "<<point_trans<<" ix= "<<ix<<" iy= "<<iy<<endl;
		    REPORT_ERROR(1,"mlf_align2d BUG: point_trans < 0");
		}
                if (iflip < nr_nomirror_flips) 
		{
		    pdf = alpha_k[refno] * (1. - mirror_fraction[refno]) * MAT_ELEM(P_phi, iy, ix);
		}
                else 
		{
		    pdf = alpha_k[refno] * mirror_fraction[refno] * MAT_ELEM(P_phi, iy, ix);
		}
		// get the starting point in the Fimg_trans vector
		int img_start = point_trans*4*dnr_points_2d + (iflip%nr_nomirror_flips)*dnr_points_2d;
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    diff = 0.;
		    // get the starting point in the Fref vector
		    int ref_start = refno*nr_psi*dnr_points_2d + ipsi*dnr_points_2d;
		    if (do_ctf_correction)
		    {
			for (int ii = 0; ii < nr_points_prob; ii++)
			{
			    tmpr = Fimg_trans[img_start + 2*ii] - ctf[ii] * Fref[ref_start + 2*ii];	
			    tmpi = Fimg_trans[img_start + 2*ii+1] - ctf[ii] * Fref[ref_start + 2*ii+1];	
			    tmpr = (tmpr * tmpr + tmpi * tmpi) / sigma2[ii];
			    diff += tmpr;
			}
		    }
		    else
		    {
			for (int ii = 0; ii < nr_points_prob; ii++)
			{
			    tmpr = Fimg_trans[img_start + 2*ii] - Fref[ref_start + 2*ii];	
			    tmpi = Fimg_trans[img_start + 2*ii+1] - Fref[ref_start + 2*ii+1];	
			    tmpr = (tmpr * tmpr + tmpi * tmpi) / sigma2[ii];
			    diff += tmpr;
			}
		    }
                    dVkij(Mweight, zero_trans, refno, irot) = diff;
		    if (debug == 9)
		    {
			cerr<<"refno= "<<refno<<" irot= "<<irot<<" diff= "<<diff<<" mindiff2= "<<mindiff2<<endl;
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
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                ix = ROUND(opt_offsets_ref[2*irefmir]);
                iy = ROUND(opt_offsets_ref[2*irefmir+1]);
                if (iflip < nr_nomirror_flips) 
		{
		    pdf = alpha_k[refno] * (1. - mirror_fraction[refno]) * MAT_ELEM(P_phi, iy, ix);
		}
                else 
		{
		    pdf = alpha_k[refno] * mirror_fraction[refno] * MAT_ELEM(P_phi, iy, ix);
		}
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    aux = dVkij(Mweight, zero_trans, refno, irot) - mindiff2;
                    // next line because of numerical precision of exp-function
                    if (aux > 1000.) weight = 0.;
                    else weight = exp(-aux) * pdf;
                    dVkij(Mweight, zero_trans, refno, irot) = weight;
                    if (weight > Pmax_refmir[irefmir]) Pmax_refmir[irefmir] = weight;
                    if (weight > maxweight)
                    {
                        maxweight = weight;
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

    if (debug==1) {cout<<"processOneImage 2 "; print_elapsed_time(t0); annotate_time(&t0);}

    // Now for all irefmir, check significant rotations...
    // and calculate their limited_translations probabilities
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            FOR_ALL_FLIPS()
            {
                irefmir = FLOOR(iflip / nr_nomirror_flips) * n_ref + refno;
                if (iflip < nr_nomirror_flips) fracpdf = alpha_k[refno] * (1. - mirror_fraction[refno]);
                else fracpdf = alpha_k[refno] * mirror_fraction[refno];
                FOR_ALL_ROTATIONS()
                {
                    irot = iflip * nr_psi + ipsi;
                    if (dVkij(Mweight, zero_trans, refno, irot) > C_fast*Pmax_refmir[irefmir])
                    {
			// get the starting point in the Fref vector
			int ref_start = refno*nr_psi*dnr_points_2d + ipsi*dnr_points_2d;
                        // expand for all limited translations
                        FOR_ALL_LIMITED_TRANSLATIONS()
                        {
                            if (itrans != zero_trans)
                            { // zero_trans has already been calculated!
                                ix = ROUND(opt_offsets_ref[2*irefmir] + Vtrans[itrans](0));
                                iy = ROUND(opt_offsets_ref[2*irefmir+1] + Vtrans[itrans](1));
                                ix = intWRAP(ix, Moffsets.startingX(), Moffsets.finishingX());
                                iy = intWRAP(iy, Moffsets.startingY(), Moffsets.finishingY());
                                if (iflip < nr_nomirror_flips) point_trans = MAT_ELEM(Moffsets, iy, ix);
                                else point_trans = MAT_ELEM(Moffsets_mirror, iy, ix);
				if (point_trans < 0 || point_trans > dim2)
				{
				    cerr<<"point_trans = "<<point_trans<<" ix= "<<ix<<" iy= "<<iy<<endl;
				    REPORT_ERROR(1,"mlf_align2d BUG: point_trans < 0");
				}
				pdf = fracpdf * MAT_ELEM(P_phi, iy, ix);
                                if (pdf > 0)
                                {
				    // get the starting point in the Fimg_trans vector
				    int img_start = point_trans*4*dnr_points_2d + (iflip%nr_nomirror_flips)*dnr_points_2d;
                                    diff = 0.;
				    if (do_ctf_correction)
				    {
					for (int ii = 0; ii < nr_points_prob; ii++)
					{
					    tmpr = Fimg_trans[img_start + 2*ii] - ctf[ii] * Fref[ref_start + 2*ii];	
					    tmpi = Fimg_trans[img_start + 2*ii+1] - ctf[ii] * Fref[ref_start + 2*ii+1];	
					    tmpr = (tmpr * tmpr + tmpi * tmpi) / sigma2[ii];
					    diff += tmpr;
					}
				    }
				    else
				    {
					for (int ii = 0; ii < nr_points_prob; ii++)
					{
					    tmpr = Fimg_trans[img_start + 2*ii] - Fref[ref_start + 2*ii];	
					    tmpi = Fimg_trans[img_start + 2*ii+1] - Fref[ref_start + 2*ii+1];	
					    tmpr = (tmpr * tmpr + tmpi * tmpi) / sigma2[ii];
					    diff += tmpr;
					}
				    }
                                    aux = diff - mindiff2;
                                    // next line because of numerical precision of exp-function
                                    if (aux > 1000.) weight = 0.;
                                    else weight = exp(-aux) * fracpdf * MAT_ELEM(P_phi, iy, ix);
                                    dVkij(Mweight, itrans, refno, irot) = weight;
                                    if (weight > maxweight)
                                    {
                                        maxweight = weight;
                                        opt_refno = refno;
                                        opt_ipsi = ipsi;
                                        opt_iflip = iflip;
                                        opt_itrans = itrans;
                                        opt_irefmir = irefmir;
                                    }
                                }
                                else dVkij(Mweight, itrans, refno, irot) = 0.;
				if (debug == 9)
				{
				    cerr<<"refno= "<<refno<<" irot= "<<irot<<" itrans= "<<itrans<<" weight= "<<dVkij(Mweight, itrans, refno, irot)<<" maxweight= "<<maxweight<<endl;
				}
                            }
                        }
                    }
                }
            }
        }
    }

    opt_offsets(0) = opt_offsets_ref[2*opt_irefmir] + Vtrans[opt_itrans](0);
    opt_offsets(1) = opt_offsets_ref[2*opt_irefmir+1] + Vtrans[opt_itrans](1);
    opt_psi = -psi_step * (opt_iflip * nr_psi + opt_ipsi) - SMALLANGLE;

    // Now that we have maxweight, determine the sums over all weights
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
                    FOR_ALL_LIMITED_TRANSLATIONS()
                    {
                        weight = dVkij(Mweight, itrans, refno, irot);
                        if (weight > SIGNIFICANT_WEIGHT_LOW*maxweight)
                        {
                            if (iflip < nr_nomirror_flips) refw[refno] += weight;
                            else refw_mirror[refno] += weight;
                            sum_refw += weight;
                        }
                        else dVkij(Mweight, itrans, refno, irot) = 0.;
                    }
                }
            }
        }
    }
    fracweight = maxweight / sum_refw;


    // Accumulate weighted sums of fraction-parameters
    FOR_ALL_MODELS()
    {
        if (!limit_rot || pdf_directions[refno] > 0.)
        {
            sumw[refno] += (refw[refno] + refw_mirror[refno]) / sum_refw;
            sumw_mirror[refno] += refw_mirror[refno] / sum_refw;
        }
    }
    
    if (debug==1) {cout<<"processOneImage 3 "; print_elapsed_time(t0); annotate_time(&t0);}

    // Acummulate all weighted sums
    // and normalize them by sum_refw, such that sum over all weights is one!
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
		    // get the starting point in the Fwsum_imgs vectors
		    int wsum_start = refno*nr_psi*dnr_points_2d + ipsi*dnr_points_2d;
                    FOR_ALL_LIMITED_TRANSLATIONS()
                    {
                        ix = ROUND(opt_offsets_ref[2*irefmir] + Vtrans[itrans](0));
                        iy = ROUND(opt_offsets_ref[2*irefmir+1] + Vtrans[itrans](1));
                        ix = intWRAP(ix, Moffsets.startingX(), Moffsets.finishingX());
                        iy = intWRAP(iy, Moffsets.startingY(), Moffsets.finishingY());
                        if (iflip < nr_nomirror_flips) point_trans = MAT_ELEM(Moffsets, iy, ix);
                        else point_trans = MAT_ELEM(Moffsets_mirror, iy, ix);
                        weight = dVkij(Mweight, itrans, refno, irot) / sum_refw;
                        if (weight > 0.)
                        {
			    // get the starting point in the Fimg_trans vector
			    int img_start = point_trans*4*dnr_points_2d + (iflip%nr_nomirror_flips)*dnr_points_2d;
                            wsum_sigma_offset += weight * (double)(ix * ix + iy * iy);
			    if (do_ctf_correction)
			    {
				for (int ii = 0; ii < nr_points_2d; ii++)
				{
				    Fwsum_imgs[wsum_start + 2*ii] += weight * decctf[ii] * Fimg_trans[img_start + 2*ii];
				    Fwsum_imgs[wsum_start + 2*ii+1] += weight * decctf[ii] * Fimg_trans[img_start + 2*ii+1];
				    Fwsum_ctfimgs[wsum_start + 2*ii] += weight * Fimg_trans[img_start + 2*ii];
				    Fwsum_ctfimgs[wsum_start + 2*ii+1] += weight * Fimg_trans[img_start + 2*ii+1];
				    tmpr = Fimg_trans[img_start + 2*ii] - ctf[ii] * Fref[wsum_start + 2*ii];	
				    tmpi = Fimg_trans[img_start + 2*ii+1] - ctf[ii] * Fref[wsum_start + 2*ii+1];
				    Mwsum_sigma2[ii] += weight * (tmpr * tmpr + tmpi * tmpi);
				}
			    }
			    else
			    {
				for (int ii = 0; ii < nr_points_2d; ii++)
				{
				    Fwsum_imgs[wsum_start + 2*ii] += weight * Fimg_trans[img_start + 2*ii];
				    Fwsum_imgs[wsum_start + 2*ii+1] += weight * Fimg_trans[img_start + 2*ii+1];
				    tmpr = Fimg_trans[img_start + 2*ii] -  Fref[wsum_start + 2*ii];	
				    tmpi = Fimg_trans[img_start + 2*ii+1] - Fref[wsum_start + 2*ii+1];
				    Mwsum_sigma2[ii] += weight * (tmpr * tmpr + tmpi * tmpi);
				}
			    }
                        }
                    }
                }
            }
        }
    }

    // Update the optimal origin offsets
    if (do_mirror) nr_mir = 2;
    else nr_mir = 1;
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

    if (debug==1) {cout<<"processOneImage 4 "; print_elapsed_time(t0); annotate_time(&t0);}

    // Compute Log Likelihood
    // 1st term: log(refw_i)
    // 2nd term: for subtracting mindiff2
    // 3rd term: for missing normalization constant
    LL += log(sum_refw) - mindiff2 - logsigma2;

}

void Prog_MLFalign2D_prm::sumOverAllImages(SelFile &SF, vector<ImageXmipp> &Iref, int iter,
					   double &LL, double &sumcorr, DocFile &DFo,
					   vector<Matrix2D<double> > &wsum_Mref,
					   vector<Matrix2D<double> > &wsum_ctfMref,
					   vector<vector<double> > &Mwsum_sigma2,
					   double &wsum_sigma_offset, 
					   vector<double> &sumw, 
					   vector<double> &sumw_mirror)
{

    ImageXmipp img;
    SelLine line;
    FileName fn_img, fn_trans;
    vector<double> Fref, Fwsum_imgs, Fwsum_ctfimgs, dum;
    vector<double> allref_offsets, pdf_directions(n_ref);
    Matrix1D<double> dataline(9), opt_offsets(2), trans(2);

    float old_phi = -999., old_theta = -999.;
    double opt_psi, opt_flip, maxcorr;
    double opt_xoff, opt_yoff;
    int c, nn, imgno, opt_refno, focus;
    bool apply_ctf;

    // Generate (FT of) each rotated version of all references
    rotateReference(Fref);

    // Initialize
    nn = SF.ImgNo();
    if (verb > 0) init_progress_bar(nn);
    c = XMIPP_MAX(1, nn / 60);
    trans.initZeros();

    // Set all weighted sums to zero
    sumw.clear();
    sumw_mirror.clear();
    Mwsum_sigma2.clear();
    Fwsum_imgs.clear();
    Fwsum_ctfimgs.clear();
    LL = 0.;
    wsum_sigma_offset = 0.;
    sumcorr = 0.;
    FOR_ALL_DEFOCUS_GROUPS()
    {
	Mwsum_sigma2.push_back(dum);
	for (int ii = 0; ii < nr_points_2d; ii++)
	{
	    Mwsum_sigma2[ifocus].push_back(0.);
	}
    }
    FOR_ALL_MODELS()
    {
        sumw.push_back(0.);
        sumw_mirror.push_back(0.);
        FOR_ALL_ROTATIONS()
	{
	    for (int ii = 0; ii < dnr_points_2d; ii++)
	    {
		Fwsum_imgs.push_back(0.);
		if (do_ctf_correction) Fwsum_ctfimgs.push_back(0.);
	    }
	}
    }

    // Loop over all images
    imgno = 0;
    SF.go_beginning();
    while ((!SF.eof()))
    {
        // For defocus-groups in Fourier-mode
        line = SF.current();
        focus = line.get_number() - 1;

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
	if (do_write_offsets)
	{
	    // A. read them from disc
	    if (!readOffsets(fn_trans, allref_offsets))
	    {
		int itot = n_ref * 2;
		if (do_mirror) itot *= 2;
		allref_offsets.clear();
		allref_offsets.resize(itot);
		for (int i = 0; i < itot; i++) allref_offsets[i] = 0.;
	    }
	}
	else
	{
	    // B. read them from memory
	    allref_offsets = imgs_offsets[imgno];
	}

        // Read optimal orientations from memory
        if (limit_rot)
        {
            old_phi = imgs_oldphi[imgno];
            old_theta = imgs_oldtheta[imgno];
        }

        // For limited orientational search: preselect relevant directions
        preselectDirections(old_phi, old_theta, pdf_directions);

	// Perform the actual expectation step for this image
	if (iter == 1 && first_iter_noctf) apply_ctf = false;
	else apply_ctf = true;
	processOneImage(img(), focus, apply_ctf, Fref,
			Fwsum_imgs, Fwsum_ctfimgs, Mwsum_sigma2[focus], wsum_sigma_offset,
			sumw, sumw_mirror, LL, maxcorr,
			opt_refno, opt_psi, opt_offsets,
			allref_offsets, pdf_directions);

	// Store optimal translations
	if (do_write_offsets) writeOffsets(fn_trans, allref_offsets);
	else imgs_offsets[imgno] = allref_offsets;

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
            dataline(0) = Iref[opt_refno].Phi();     // rot
            dataline(1) = Iref[opt_refno].Theta();   // tilt
            dataline(2) = opt_psi + 360.;            // psi
            dataline(3) = trans(0) + opt_offsets(0); // Xoff
            dataline(4) = trans(1) + opt_offsets(1); // Yoff
            dataline(5) = (double)(opt_refno + 1);   // Ref
            dataline(6) = opt_flip;                  // Mirror
            dataline(7) = maxcorr;                   // P_max/P_tot or Corr
            DFo.append_comment(img.name());
            DFo.append_data_line(dataline);
        }

        // Output docfile
        if (verb > 0) if (imgno % c == 0) progress_bar(imgno);
        imgno++;
    }
    if (verb > 0) progress_bar(nn);

    if (do_ctf_correction) reverseRotateReference(Fwsum_ctfimgs,wsum_ctfMref);
    reverseRotateReference(Fwsum_imgs,wsum_Mref);

}

// Update all model parameters
void Prog_MLFalign2D_prm::updateParameters(vector<Matrix2D<double> > &wsum_Mref,
					   vector<Matrix2D<double> > &wsum_ctfMref,
					   vector<vector<double> > &Mwsum_sigma2,
					   double &wsum_sigma_offset,
					   vector<double> &sumw, vector<double> &sumw_mirror,
					   double &sumcorr, double &sumw_allrefs,
					   Matrix1D<double> &spectral_signal)
{

    Matrix1D<double> rmean_sigma2, rmean_signal2;
    Matrix1D<int> center(2), radial_count;
    Matrix2D<complex<double> > Faux, Faux2;
    Matrix2D<double> Maux;
    FileName fn_tmp;
    double rr, thresh, aux;
    int c;

    // Pre-calculate sumw_allrefs & average Pmax/sumP or cross-correlation
    sumw_allrefs = 0.;
    FOR_ALL_MODELS()
    {
        sumw_allrefs += sumw[refno];
    }
    sumcorr /= sumw_allrefs;

    // Update the reference images
    FOR_ALL_MODELS()
    {
        if (sumw[refno] > 0.)
        {
            Iref[refno]() = wsum_Mref[refno];
            Iref[refno]() /= sumw[refno];
            Iref[refno].weight() = sumw[refno];
	    if (do_ctf_correction)
	    {
		Ictf[refno]() = wsum_ctfMref[refno];
		Ictf[refno]() /= sumw[refno];
		Ictf[refno].weight() = sumw[refno];
	    }
	    else
	    {
		Ictf[refno]=Iref[refno];
	    }
        }
        else
        {
            Iref[refno].weight() = 0.;
            Iref[refno]().initZeros(dim, dim);
	    Ictf[refno]().initZeros();
	    Ictf[refno].weight() = 0.;
        }
    }

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
    if (!fix_sigma_offset) sigma_offset = sqrt(wsum_sigma_offset / (2 * sumw_allrefs));

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
	    // Use current_highres_shell-1 to prevent artefacts from
	    for (int irr = 0; irr < current_highres_shell; irr++)
	    {
		dVi(Vsig[ifocus], irr) = dVi(rmean_sigma2, irr) / (double)(2 * count_defocus[ifocus]);
	    }
	}
    }

    // Calculate average spectral signal
    c = 0;
    FOR_ALL_MODELS()
    {
	if (sumw[refno] > 0.)
	{
	    FourierTransform(Ictf[refno](), Faux);
	    FFT_magnitude(Faux, Maux);
	    CenterFFT(Maux, true);
	    Maux *= Maux;
	    Maux *= sumw[refno];
	    center.initZeros();
	    rmean_signal2.initZeros();
	    Maux.setXmippOrigin();
	    radialAverage(Maux, center, rmean_signal2, radial_count, true);
	    if (c == 0) spectral_signal = rmean_signal2;
	    else spectral_signal += rmean_signal2;
	    c++;
	}
    }
    spectral_signal /= (double)n_ref;

}

// Check convergence
bool Prog_MLFalign2D_prm::checkConvergence(vector<double> &conv)
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

void Prog_MLFalign2D_prm::writeOutputFiles(const int iter, DocFile &DFo,
					   double &sumw_allrefs, double &LL, 
					   double &avecorr, vector<double> &conv)
{

    FileName          fn_tmp, fn_base, fn_tmp2;
    Matrix1D<double>  fracline(3);
    SelFile           SFo, SFc;
    DocFile           DFl;
    string            comment;
    ofstream          fh;

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
	    system(((string)"rm -rf "+fn_scratch).c_str());
	}
    
    }

    // Write out current reference images and fill sel & log-file
    FOR_ALL_MODELS()
    {
	fn_tmp = fn_root + "_cref";
	fn_tmp.compose(fn_tmp, refno + 1, "");
	fn_tmp = fn_tmp + ".xmp";
	Ictf[refno].write(fn_tmp);
	SFc.insert(fn_tmp, SelLine::ACTIVE);
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
    fn_tmp = fn_root + "_cref.sel";
    SFc.write(fn_tmp);

    DFl.go_beginning();
    comment = "MLFalign2D-logfile: Number of images= " + floatToString(sumw_allrefs);
    comment += " LL= " + floatToString(LL, 15, 10) + " <Pmax/sumP>= " + floatToString(avecorr, 10, 5);
    DFl.insert_comment(comment);
    comment = " -offset " + floatToString(sigma_offset, 15, 12) + " -istart " + integerToString(iter + 1);
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

    // Write out updated Vsig vectors
    if (!fix_sigma_noise)
    {
	FOR_ALL_DEFOCUS_GROUPS()
        {
	    fn_tmp = fn_base + "_sig";
	    if (nr_focus > 1) fn_tmp.compose(fn_tmp, ifocus + 1, "");
	    fn_tmp += ".dat";
	    fh.open((fn_tmp).c_str(), ios::out);
	    if (!fh) REPORT_ERROR(1, (string)"Prog_MLFalign2D_prm: Cannot write file: " + fn_tmp);
            for (int irr = 0; irr < hdim; irr++)
            {
                fh << irr/(sampling*dim) << " " << dVi(Vsig[ifocus], irr) << "\n";
            }
	    fh.close();
	}
    }

}

