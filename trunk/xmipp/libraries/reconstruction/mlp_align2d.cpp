/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
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
#include "mlp_align2d.h"

// Read arguments ==========================================================
void Prog_MLPalign2D_prm::read(int argc, char **argv, bool ML3D)
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
	extendedUsage();
    }
    nr_ref = textToInteger(getParameter(argc, argv, "-nref", "0"));
    fn_ref = getParameter(argc, argv, "-ref", "");
    fn_sel = getParameter(argc, argv, "-i");
    fn_root = getParameter(argc, argv, "-o", "mlp2d");
    Niter = textToInteger(getParameter(argc, argv, "-iter", "100"));
    istart = textToInteger(getParameter(argc, argv, "-istart", "1"));
    sigma_noise = textToFloat(getParameter(argc, argv, "-noise", "1"));
    sigma_offset = textToFloat(getParameter(argc, argv, "-offset", "3"));
    do_mirror = checkParameter(argc, argv, "-mirror");
    eps = textToFloat(getParameter(argc, argv, "-eps", "5e-5"));
    fn_frac = getParameter(argc, argv, "-frac", "");
    fix_fractions = checkParameter(argc, argv, "-fix_fractions");
    fix_sigma_noise = checkParameter(argc, argv, "-fix_sigma_noise");
    fix_sigma_offset = checkParameter(argc, argv, "-fix_sigma_offset");
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    search_shift = textToFloat(getParameter(argc, argv, "-search_shift", "3."));
    fn_doc = getParameter(argc, argv, "-doc", "");
    Ri = textToInteger(getParameter(argc,argv,"-Ri","-1"));
    Ro = textToInteger(getParameter(argc,argv,"-Ro","-1"));
    do_ML3D = ML3D;

    // Only for interaction with refine3D:
    search_rot = textToFloat(getParameter(argc, argv, "-search_rot", "999."));

}

// Show ====================================================================
void Prog_MLPalign2D_prm::show(bool ML3D)
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
        std::cerr << "--> Using polar coordinates and limited translational searches " << std::endl;
        std::cerr << "  Input images            : " << fn_sel << " (" << nr_exp_images << ")" << std::endl;
        if (fn_ref != "")
            std::cerr << "  Reference image(s)      : " << fn_ref << std::endl;
        else
            std::cerr << "  Number of references:   : " << nr_ref << std::endl;
        std::cerr << "  Output rootname         : " << fn_root << std::endl;
        std::cerr << "  Stopping criterium      : " << eps << std::endl;
        std::cerr << "  initial sigma noise     : " << sigma_noise << std::endl;
        std::cerr << "  initial sigma offset    : " << sigma_offset << std::endl;
        if (do_mirror)
            std::cerr << "  Check mirrors           : true" << std::endl;
        else
            std::cerr << "  Check mirrors           : false" << std::endl;
        if (fn_frac != "")
            std::cerr << "  Initial model fractions : " << fn_frac << std::endl;

	std::cerr << "    + Limit polar coordinates to radii between "<<Ri<<" and "<<Ro<<" pixels"<<std::endl; 
        if (search_shift < 999.)
            std::cerr << "    + Limit translational search to +/- " << search_shift << " pixels" << std::endl;
        if (search_rot < 180.)
            std::cerr << "    + Limit orientational search to +/- " << search_rot << " degrees" << std::endl;

        // Hidden stuff
        if (fix_fractions)
        {
            std::cerr << "  -> Do not update estimates of model fractions." << std::endl;
        }
        if (fix_sigma_noise)
        {
            std::cerr << "  -> Do not update sigma-estimate of noise." << std::endl;
        }
        std::cerr << " -----------------------------------------------------------------" << std::endl;

    }

}

// Usage ===================================================================
void Prog_MLPalign2D_prm::usage()
{
    std::cerr << "Usage:  mlp_align2d [options] " << std::endl;
    std::cerr << "   -i <selfile>                : Selfile with input images \n";
    std::cerr << "   -ref <selfile/image>        : Selfile with initial references/single reference image \n";
    std::cerr << "      OR -nref <int>               OR number of references to generate automatically (bias-free)\n";
    std::cerr << " [ -o <rootname> ]             : Output rootname (default = \"ml2d\")\n";
    std::cerr << " [ -mirror ]                   : Also check mirror image of each reference \n";
    std::cerr << " [ -Ri <float=1> ]             : Inner radius to limit rotational search \n";
    std::cerr << " [ -Ro <float=dim/2 - 1> ]     : Outer radius to limit rotational search \n";
    std::cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_MLPalign2D_prm::extendedUsage(bool ML3D)
{
    std::cerr << "Additional options: " << std::endl;
    std::cerr << " [ -eps <float=5e-5> ]         : Stopping criterium \n";
    std::cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
    std::cerr << " [ -noise <float=1> ]          : Expected standard deviation for pixel noise \n";
    std::cerr << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n";
    if (!ML3D) std::cerr << " [ -restart <logfile> ]        : restart a run with all parameters as in the logfile \n";
    if (!ML3D) std::cerr << " [ -istart <int> ]             : number of initial iteration \n";
    std::cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise \n";
    std::cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
    std::cerr << " [ -search_shift <float=999>]  : Limit of translational search [pix] (does NOT use FFT) \n";
    std::cerr << " [ -doc <docfile=\"\"> ]         : Read initial angles and offsets from docfile \n";
    std::cerr << std::endl;
    exit(1);
}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produce_Side_info2)
void Prog_MLPalign2D_prm::produceSideInfo()
{

    // Read selfile with experimental images
    SF.read(fn_sel);
    // Get image sizes and total number of images
    SF.ImgSize(dim, dim);
    nr_exp_images = SF.ImgNo();

    // Automatically set maximum ring radii
    if (Ri<0) Ri=1;
    if (Ro<0) Ro=(dim/2)-1;

    // Radii limits: between 1 and half dimension
    Ri=MAX(1.,Ri);
    Ro=MIN((dim/2)-1, Ro);

    // Get number of references
    if (do_ML3D)
    {
        do_generate_refs = false;
    }
    else if (fn_ref != "")
    {
        do_generate_refs = false;
        if (Is_ImageXmipp(fn_ref)) nr_ref = 1;
        else
        {
            SFr.read(fn_ref);
            nr_ref = SFr.ImgNo();
        }
    }
    else
    {
        do_generate_refs = true;
    }

    // Set-up limit_rot & limit_trans
    if (search_rot < 180.) limit_rot = true;
    else limit_rot = false;
    nr_trans = 0;
    Vxtrans.clear();
    Vytrans.clear();
    for (int ix = -search_shift; ix <= search_shift; ix++)
    {
	for (int iy = -search_shift; iy <= search_shift; iy++)
	{
	    double x = (double)ix;
	    double y = (double)iy;
	    if (x*x + y*y <= search_shift * search_shift)
	    {
		Vxtrans.push_back(x);
		Vytrans.push_back(y);
		nr_trans++;
	    }
	}
    }

}

// Generate initial references =============================================
void Prog_MLPalign2D_prm::generateInitialReferences()
{

    SelFile SFtmp, SFout;
    ImageXmipp Iave, Itmp;
    double dummy;
    FileName fn_tmp;
    SelLine line;

    if (verb > 0)
    {
        std::cerr << "  Generating initial references by averaging over random subsets" << std::endl;
        init_progress_bar(nr_ref);
    }

    // Make random subsets and calculate average images
    SFtmp = SF.randomize();
    int Nsub = ROUND((double)SFtmp.ImgNo() / nr_ref);
    for (int refno = 0; refno < nr_ref; refno++)
    {
        SFout.clear();
        SFout.reserve(Nsub);
        SFtmp.go_beginning();
        SFtmp.jump_lines(Nsub*refno);
        if (refno == nr_ref - 1) Nsub = SFtmp.ImgNo() - refno * Nsub;
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
    if (verb > 0) progress_bar(nr_ref);
    fn_ref = fn_root + "_it";
    fn_ref.compose(fn_ref, 0, "sel");
    SFr.write(fn_ref);

}

// Read reference images to memory
// This side info is NOT general, i.e. in parallel mode it is NOT the
// same for all processors! (in contrast to produce_Side_info)
void Prog_MLPalign2D_prm::produceSideInfo2(int nr_vols)
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
    nr_ref = 0;
    SFr.go_beginning();
    while ((!SFr.eof()))
    {
        img.read(SFr.NextImg(), false, false, true, false);
        img().setXmippOrigin();
        Iref.push_back(img);
	Iold.push_back(img);
        // Default start is all equal model fractions
        alpha_k.push_back((double)1 / SFr.ImgNo());
        Iref[refno].weight() = alpha_k[refno] * (double)nr_exp_images;
        // Default start is half-half mirrored images
        if (do_mirror) mirror_fraction.push_back(0.5);
        else mirror_fraction.push_back(0.);
        nr_ref++;
        refno++;
    }

    // Calculate FT of rings of polars for all references
    calculateFtRingsAllRefs(Iref,fP_refs,fP_zero,sum2_refs,Ri,Ro);

    // Read optimal origin offsets from fn_doc
    if (fn_doc != "")
        DF.read(fn_doc);

    // Fill vectors with optimal positions of all images
    SF.go_beginning();
    imgs_oldphi.clear();
    imgs_oldtheta.clear();
    imgs_oldxoff.clear();
    imgs_oldyoff.clear();
    while (!SF.eof())
    {
	fn_tmp = SF.NextImg();
	if (fn_doc != "")
	{
	    if (DF.search_comment(fn_tmp))
            {
		imgs_oldphi.push_back(DF(0));
		imgs_oldtheta.push_back(DF(1));
                imgs_oldxoff.push_back(DF(3));
                imgs_oldyoff.push_back(DF(4));
            }
            else
            {
                REPORT_ERROR(1, (std::string)"Prog_MLPalign2D_prm: Cannot find " + fn_tmp + " in docfile " + fn_doc);
            }
	    DF.go_beginning();
        }
	else
	{
	    imgs_oldphi.push_back(-999.);
	    imgs_oldtheta.push_back(-999.);
	    imgs_oldxoff.push_back(0.);
	    imgs_oldyoff.push_back(0.);
	}
    }
    DF.clear();

    // read in model fractions if given on command line
    if (fn_frac != "")
    {
        DF.read(fn_frac);
        DF.go_first_data_line();
        for (refno = 0; refno < nr_ref; refno++)
        {
            DL = DF.get_current_line();
            alpha_k[refno] = DL[0];
            if (do_mirror)
            {
                if (DL[1] > 1. || DL[1] < 0.)
                    REPORT_ERROR(1, "Prog_MLPalign2D_prm: Mirror fraction (2nd column) should be [0,1]!");
                mirror_fraction[refno] = DL[1];
            }
            sumfrac += alpha_k[refno];
            DF.next_data_line();
        }
        if (ABS(sumfrac - 1.) > 1e-3)
            if (verb > 0) std::cerr << " ->WARNING: Sum of all expected model fractions (" << sumfrac << ") is not one!" << std::endl;
        for (refno = 0; refno < nr_ref; refno++)
        {
            alpha_k[refno] /= sumfrac;
        }
    }

}

void Prog_MLPalign2D_prm::preselectDirections(float &phi, float &theta,
					      std::vector<double> &pdf_directions)
{

    float phi_ref, theta_ref, angle, angle2;
    Matrix1D<double> u, v;

    pdf_directions.clear();
    pdf_directions.resize(nr_ref);
    for (int iref = 0; iref < nr_ref; iref++)
    {
        if (!limit_rot || (phi == -999. && theta == -999.)) pdf_directions[iref] = 1.;
        else
        {
            phi_ref = Iref[iref].Phi();
            theta_ref = Iref[iref].Theta();
            Euler_direction(phi, theta, 0., u);
            Euler_direction(phi_ref, theta_ref, 0., v);
            u.normalize();
            v.normalize();
            angle = RAD2DEG(acos(dotProduct(u, v)));
            angle = fabs(realWRAP(angle, -180, 180));
            // also check mirror
            angle2 = 180. + angle;
            angle2 = fabs(realWRAP(angle2, -180, 180));
            angle = MIN(angle, angle2);
            if (fabs(angle) > search_rot) pdf_directions[iref] = 0.;
            else pdf_directions[iref] = 1.;
        }
    }

}

void Prog_MLPalign2D_prm::updatePdfTranslations()
{
    double r2, pdfpix, sum;
    Mpdf_trans.resize(dim, dim);
    Mpdf_trans.setXmippOrigin();
    Mr2.resize(dim, dim);
    Mr2.setXmippOrigin();

    sum=0.;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Mpdf_trans)
    {
        r2 = (double)(j * j + i * i);
        if (sigma_offset > 0.)
        {
            pdfpix = exp(-r2 / (2 * sigma_offset * sigma_offset));
        }
        else
        {
            if (j == 0 && i == 0) pdfpix = 1.;
            else pdfpix = 0.;
        }
        MAT_ELEM(Mpdf_trans, i, j) = pdfpix;
        MAT_ELEM(Mr2, i, j) = (double)r2;
	sum += pdfpix;
    }

    // Normalization (for limited searches this is not completely correct!)
    Mpdf_trans /= sum;

}

void Prog_MLPalign2D_prm::calculateFtRingsAllRefs(const std::vector< ImageXmippT<double> > &Iref,
						  std::vector< Polar< std::complex <double> > > &fP_refs,
						  Polar< std::complex <double> > &fP_zero,
						  std::vector< double > &sum2_refs,
						  const int &first, const int &last)
{
	Matrix2D<double>         Maux;
	Polar<double>            P;
	Polar<std::complex <double> > fP;
	std::vector<double>           x,y,data;

	fP_refs.clear();
	sum2_refs.clear();
	for (int iref = 0; iref < nr_ref; iref++)
	{
	    // Calculate polars using gridding interpolation
	    produceReverseGriddingMatrix2D(Iref[iref](),Maux,kb);
	    P.getPolarFromCartesian(Maux,kb,first,last);
	    if (iref == 0)
	    {
		// Pre-calculate voronoi weights for the polar structure
		P.getCartesianCoordinates(x,y,data);
		approximateVoronoiArea(voronoi_area,x,y);
	    }
	    // Store complex conjugated
	    fP = P.fourierTransformRings(true); 
	    fP_refs.push_back(fP);
	    // Also store A2 (power of the reference)
	    sum2_refs.push_back(P.computeSum2());

	    // nr_psi is number of samples in last ring
	    // dnr_psi is 2 x nr_psi for mirror searches
	    nr_psi = P.getSampleNo(P.getRingNo()-1);
	    if (do_mirror) 
		dnr_psi = 2 * nr_psi;
	    else
		dnr_psi = nr_psi;

	    // Calculate an all-zero Polar
	    if (iref == 0)
	    {
		fP_zero = fP;
		fP_zero *= 0.;
	    }
	}

}

void Prog_MLPalign2D_prm::calculateFtRingsAllTransImg(const ImageXmippT<double> &img,
						      std::vector< Polar< std::complex <double> > > &fP_trans,
						      std::vector< Polar< std::complex <double> > > &fPm_trans,
						      double &Xi2, const int &first, const int &last)
{
    Matrix2D<double>         Maux;
    Polar<double>            P;
    Polar<std::complex <double> > fP;

    produceReverseGriddingMatrix2D(img(),Maux,kb);
    fP_trans.clear();
    fPm_trans.clear();
    for (int itrans = 0; itrans < nr_trans; itrans++)
    {
	P.getPolarFromCartesian(Maux,kb,first,last,Vxtrans[itrans],Vytrans[itrans]);
	fP = P.fourierTransformRings(false);
	fP_trans.push_back(fP);
	if (ABS(Vxtrans[itrans]) < XMIPP_EQUAL_ACCURACY &&
	    ABS(Vytrans[itrans]) < XMIPP_EQUAL_ACCURACY)
	    Xi2 = P.computeSum2();
	
	if (do_mirror)
	{
	    fP = P.fourierTransformRings(true);
	    fPm_trans.push_back(fP);
	}
    }

}

void Prog_MLPalign2D_prm::processOneImage(const ImageXmippT<double> &img,
					  const std::vector < Polar <std::complex <double> >  > &fP_ref,
					  const std::vector < double > &sum2_refs,
					  const std::vector < double > &pdf_directions,
					  std::vector < Polar <std::complex <double> > > &fP_wsum_imgs,
					  double &wsum_sigma_noise, double &wsum_sigma_offset, 
					  std::vector < double > &sumw, std::vector < double > &sumw_mirror,
					  double &LL, double &fracweight,
					  int &opt_iref, double &opt_psi, double &opt_flip, 
					  double &opt_xoff, double &opt_yoff)
{

    std::vector<double>                      refw, refw_mirror;
    Matrix1D<double>                    ang,corr;
    std::vector< Polar< std::complex <double> > > fP_trans, fPm_trans;
    double                              A2, Xi2, A2_plus_Xi2;
    double                              aux, weight, maxweight, diff2, mindiff2;
    double                              wsum_corr, pdf, pdfm, pdft, sumw_allrefs;
    int                                 ipos;
    double                              sigma_noise2 = sigma_noise * sigma_noise;
    int                                 opt_ipsi, opt_ipos, opt_itrans, nr_pos;
    Matrix1D<double>                    Mweight(nr_psi);
    Matrix1D<std::complex <double> >         Fweight(nr_psi);
    std::vector<int>                         img_xoffs, img_yoffs;

    //TimeStamp t0,t1; 
    //time_config();
    //annotate_time(&t0);

    // reserve memory for the weights vector
    nr_pos = nr_trans * nr_ref * dnr_psi;
    std::vector<double> weights;

    // Store actual offsets for pdf calculations
    img_xoffs.clear();
    img_yoffs.clear();
    for (int itrans = 0; itrans < nr_trans; itrans++)
    {
	img_xoffs.push_back(img.Xoff() - Vxtrans[itrans]);
	img_yoffs.push_back(img.Yoff() - Vytrans[itrans]);
    }

    // Calculate all FT-rings for all translated images
    calculateFtRingsAllTransImg(img,fP_trans,fPm_trans,Xi2,Ri,Ro);
    //std::cerr<<"generate polars: "; print_elapsed_time(t0);
    //annotate_time(&t1);

   // Search over all references
    mindiff2 = 99.e99;
    ipos = 0;
    for (int iref = 0; iref < nr_ref; iref++)
    {
	A2 = sum2_refs[iref];
	// If within limited rot-search criteria
        if (!limit_rot || pdf_directions[iref] > 0.)
	{
	    // Search over all translations
	    for (int itrans = 0; itrans < nr_trans; itrans++)
	    {
		A2_plus_Xi2 = 0.5 * (A2 + Xi2);
		// A. Check straight image
		// Search all in-plane rotations via convolution theorem
		rotationalCorrelation(fP_trans[itrans],fP_ref[iref],ang,corr);
		for (int ipsi = 0; ipsi < nr_psi; ipsi++)
		{
		    ipos++;
		    diff2 = A2_plus_Xi2 - corr(ipsi);
		    weights.push_back(diff2);
                    if (diff2 < mindiff2) mindiff2 = diff2;
		}
		// B. Check mirror image
		if (do_mirror) 
		{
		    rotationalCorrelation(fPm_trans[itrans],fP_ref[iref],ang,corr);
		    for (int ipsi = 0; ipsi < nr_psi; ipsi++)
		    {
			ipos++;
			diff2 = A2_plus_Xi2 - corr(ipsi);
			weights.push_back(diff2);
			if (diff2 < mindiff2) mindiff2 = diff2;
		    }
		}
	    }
	}
    }

    // Now that we have mindiff2 calculate all probabilities and maxweight
    ipos = 0;
    wsum_corr = 0.;
    sumw_allrefs = 0.;
    maxweight = 0.;
    for (int iref = 0; iref < nr_ref; iref++)
    {
	refw.push_back(0.);
        refw_mirror.push_back(0.);
	pdf = alpha_k[iref] * pdf_directions[iref] * (1. - mirror_fraction[iref]);
	pdfm = alpha_k[iref] * pdf_directions[iref] * mirror_fraction[iref];
        if (!limit_rot || pdf_directions[iref] > 0.)
	{
	    for (int itrans = 0; itrans < nr_trans; itrans++)
	    {
		pdft = MAT_ELEM(Mpdf_trans, img_yoffs[itrans], img_xoffs[itrans]);
		for (int ipsi = 0; ipsi < dnr_psi; ipsi++)
		{
		    ipos++;
		    diff2 = weights[ipos - 1];
                    aux = (diff2 - mindiff2) / (sigma_noise2);
                    // next line because of numerical precision of exp-function
                    if (aux > 1000.) weight = 0.;
		    else 
		    {
			if (ipsi<nr_psi) 
			{   // straight image
			    weight = exp(-aux) * pdf * pdft;
			    refw[iref] += weight;
			}
			else
			{   // mirror image
			    weight = exp(-aux) * pdfm * pdft;
			    refw_mirror[iref] += weight;
			}
		    }
		    weights[ipos - 1] = weight;
		    sumw_allrefs += weight;
		    wsum_corr += weight * diff2;
		    if (weight > maxweight)
		    {
			maxweight = weight;
                        opt_iref = iref;
                        opt_itrans = itrans;
                        opt_ipsi = ipsi;
			opt_ipos = ipos;
		    }
		}
	    }
	}
    }
    
    // Normalize all weighted sums by sumw such that sum over all weights is one!
    // And accumulate all weighted sums
    wsum_sigma_noise += (2. * wsum_corr / sumw_allrefs);
    ipos = 0;
    for (int iref = 0; iref < nr_ref; iref++)
    {
	sumw[iref] += (refw[iref] + refw_mirror[iref]) / sumw_allrefs;
	sumw_mirror[iref] += refw_mirror[iref] / sumw_allrefs;
        if (!limit_rot || pdf_directions[iref] > 0.)
	{
	    for (int itrans = 0; itrans < nr_trans; itrans++)
	    {
		// A. Store sum of straight images (via convolution theorem)
		aux = 0.;
		for (int ipsi = 0; ipsi < nr_psi; ipsi++)
		{
		    ipos++;
		    weight = weights[ipos - 1] / sumw_allrefs;
		    Mweight(ipsi) = weight;
		    if (weight > aux) aux = weight; 
		    wsum_sigma_offset += weight * MAT_ELEM(Mr2, img_yoffs[itrans], img_xoffs[itrans]); 
		}
		// Only store sum for non-zero weights
		if (aux > SIGNIFICANT_WEIGHT_LOW*maxweight)
		{
		    FourierTransformHalf(Mweight,Fweight);
		    for (int iring = 0; iring < fP_wsum_imgs[iref].getRingNo(); iring++)
			for (int iphi = 0; iphi < fP_wsum_imgs[iref].getSampleNo(iring); iphi++)
			    fP_wsum_imgs[iref](iring,iphi) += (double)(nr_psi) * conj(Fweight(iphi)) *
				fP_trans[itrans](iring,iphi);
		}
		if (do_mirror)
		{
		    // B. Store sum of mirror images 
		    aux = 0.;
		    for (int ipsi = 0; ipsi < nr_psi; ipsi++)
		    {
			ipos++;
			weight = weights[ipos - 1] / sumw_allrefs;
			Mweight(ipsi) = weight;
			if (weight > aux) aux = weight; 
			wsum_sigma_offset += weight * MAT_ELEM(Mr2, img_yoffs[itrans], img_xoffs[itrans]); 
		    }
		    if (aux > SIGNIFICANT_WEIGHT_LOW*maxweight)
		    {
			FourierTransformHalf(Mweight,Fweight);
			for (int iring = 0; iring < fP_wsum_imgs[iref].getRingNo(); iring++)
			    for (int iphi = 0; iphi < fP_wsum_imgs[iref].getSampleNo(iring); iphi++)
				fP_wsum_imgs[iref](iring,iphi) += (double)(nr_psi) * conj(Fweight(iphi)) *
				    fPm_trans[itrans](iring,iphi);
		    }
		}
	    }
	}
    }

    // Compute Log Likelihood
    // 1st term: log(refw_i)
    // 2nd term: for subtracting mindiff2
    // 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
    LL += log(sumw_allrefs) - mindiff2 / sigma_noise2 - dim * dim * log(2.50663 * sigma_noise);

    // Return values
    fracweight = maxweight / sumw_allrefs;
    opt_xoff = - Vxtrans[opt_itrans];
    opt_yoff = - Vytrans[opt_itrans];
    if (opt_ipsi>nr_psi-1)
    {
	opt_flip = 1.;
	opt_psi = ang(opt_ipsi-nr_psi);
    }
    else
    {
	opt_flip = 0.;
	opt_psi = ang(opt_ipsi);
    }
    //std::cerr<<"all the rest: "; print_elapsed_time(t1);

}


void Prog_MLPalign2D_prm::sumOverAllImages(SelFile &SF, const std::vector< ImageXmippT<double> > &Iref,
					   double &LL, double &sumcorr, DocFile &DFo,
					   std::vector < Polar <std::complex <double> > > &fP_wsum_imgs,
					   double &wsum_sigma_noise, double &wsum_sigma_offset,  
					   std::vector<double> &sumw, std::vector<double> &sumw_mirror)
{

    ImageXmipp        img;
    FileName          fn_img;
    std::vector<double>    pdf_directions(nr_ref);
    Matrix1D<double>  dataline(8), opt_offsets(2), trans(2);
    float             old_phi = -999., old_theta = -999.;
    double            opt_psi, opt_flip, opt_xoff, opt_yoff, maxcorr;
    int               c, nn, imgno, opt_iref;
    //std::vector<Polar<std::complex<double> > > fP_refs;
    //Polar<std::complex<double> >          fP_zero;

    // Calculate FTs of polar rings of all references
    //calculateFtRingsAllRefs(Iref,fP_refs,fP_zero,sum2_refs,Ri,Ro);

    // Initialize
    LL = 0.;
    nn = SF.ImgNo();
    if (verb > 0) init_progress_bar(nn);
    c = MAX(1, nn / 60);

    sumw.clear();
    sumw_mirror.clear();
    LL = 0.;
    wsum_sigma_noise = 0.;
    wsum_sigma_offset = 0.;
    sumcorr = 0.;
    trans.initZeros();
    fP_wsum_imgs.clear();
    for (int iref = 0; iref < nr_ref; iref++)
    {
        sumw.push_back(0.);
        sumw_mirror.push_back(0.);
	fP_wsum_imgs.push_back(fP_zero);
    }
    
    // Loop over all images
    imgno = 0;
    SF.go_beginning();
    while ((!SF.eof()))
    {
        fn_img = SF.NextImg();
        img.read(fn_img, false, false, false, false);
        img().setXmippOrigin();

	// Set Xoff and Yoff from docfile or previous iteration
	img.Xoff() = trans(0) = (double)ROUND(imgs_oldxoff[imgno]);
	img.Yoff() = trans(1) = (double)ROUND(imgs_oldyoff[imgno]);
	img().selfTranslate(trans, true);
	
       // Read optimal orientations from memory
        if (limit_rot)
        {
            old_phi = imgs_oldphi[imgno];
            old_theta = imgs_oldtheta[imgno];
        }
	
        // For limited orientational search: preselect relevant directions
        preselectDirections(old_phi, old_theta, pdf_directions);
	
	// Perform the actual expectation step for the current image
	processOneImage(img, fP_refs, sum2_refs, pdf_directions, fP_wsum_imgs, 
			wsum_sigma_noise, wsum_sigma_offset, sumw, sumw_mirror, LL, maxcorr,
			opt_iref, opt_psi, opt_flip, opt_xoff, opt_yoff);
	
        // Store optimal position in memory
	imgs_oldxoff[imgno] += opt_xoff;
	imgs_oldyoff[imgno] += opt_yoff;
        if (limit_rot)
        {
            imgs_oldphi[imgno] = Iref[opt_iref].Phi();
            imgs_oldtheta[imgno] = Iref[opt_iref].Theta();
        }
	
        // Output docfile
        sumcorr += maxcorr;
	dataline(0) = Iref[opt_iref].Phi();      // rot
	dataline(1) = Iref[opt_iref].Theta();    // tilt
	dataline(2) = opt_psi;                   // psi
	dataline(3) = trans(0) + opt_xoff;       // Xoff
	dataline(4) = trans(1) + opt_yoff;       // Yoff
	dataline(5) = (double)(opt_iref + 1);    // Ref
	dataline(6) = opt_flip;                  // Mirror
	dataline(7) = maxcorr;                   // Pmax/sumP
	DFo.append_comment(img.name());
	DFo.append_data_line(dataline);

	// Output docfile
	if (verb > 0) if (imgno % c == 0) progress_bar(imgno);
	imgno++;
    }
    if (verb > 0) progress_bar(nn);
    
}

// Update all model parameters
void Prog_MLPalign2D_prm::updateParameters(std::vector < Polar <std::complex <double> > > &fP_wsum_imgs,
					   double &wsum_sigma_noise, double &wsum_sigma_offset,
					   std::vector<double> &sumw, std::vector<double> &sumw_mirror,
					   double &sumcorr, double &sumw_allrefs)
{


    // Pre-calculate sumw_allrefs & average Pmax/sumP or cross-correlation
    sumw_allrefs = 0.;
    for (int iref = 0; iref < nr_ref; iref++)
	sumw_allrefs += sumw[iref];

    // Average corr/fracweight
    sumcorr /= sumw_allrefs;

    // Update the reference images
    Polar<double> P;
    std::vector<double> x, y, data;
    Matrix2D<double> Maux;
    for (int iref = 0; iref < nr_ref; iref++)
    {
        if (sumw[iref] > 0.)
        {
	    // For now, keep everything in polar coordinates
	    // For references: store the complex conjugated FTs!
	    for (int iring = 0; iring < fP_refs[iref].getRingNo(); iring++)
	    {
		for (int iphi = 0; iphi < fP_refs[iref].getSampleNo(iring); iphi++)
		{
		    fP_refs[iref](iring,iphi) = conj(fP_wsum_imgs[iref](iring,iphi));
		}
	    }
	    fP_refs[iref] /= sumw[iref];

	    inverseFourierTransformRings(fP_refs[iref],P,true);
	    // Store updated A2 (power of the reference)
	    sum2_refs[iref] = P.computeSum2();
	    // Gridding to get back to cartesian coordinates
	    P.getCartesianCoordinates(x,y,data);
	    Maux = interpolateCartesianFromArbitrarySampling(dim,dim,x,y,data,voronoi_area,kb);
	    produceForwardGriddingMatrix2D(Maux,Iref[iref](),kb);
            Iref[iref].weight() = sumw[iref];
	}
        else
        {
            fP_refs[iref] *= 0.;
            Iref[iref]().initZeros(dim, dim);	    
	    Iref[iref].weight() = 0.;
        }
    }

    // Update sigma of the origin offsets
    if (!fix_sigma_offset) sigma_offset = sqrt(wsum_sigma_offset / (2. * sumw_allrefs));

    // Update the model fractions
    if (!fix_fractions)
    {
	for (int iref = 0; iref < nr_ref; iref++)
        {
            if (sumw[iref] > 0.)
            {
                alpha_k[iref] = sumw[iref] / sumw_allrefs;
                mirror_fraction[iref] = sumw_mirror[iref] / sumw[iref];
            }
            else
            {
                alpha_k[iref] = 0.;
                mirror_fraction[iref] = 0.;
            }
        }
    }

    // Update the noise parameters
    if (!fix_sigma_noise)
    {
	double Rop = MIN(dim/2., Ro + 0.5);
	double Rip = MAX(0.,     Ri - 0.5);
	double nr_pixels = PI * (Rop * Rop - Rip * Rip);
	sigma_noise = sqrt(wsum_sigma_noise / (sumw_allrefs * nr_pixels));
	// prevent division by zero
	sigma_noise = MAX(sigma_noise, SIGNIFICANT_WEIGHT_LOW);
    }

}

// Check convergence
bool Prog_MLPalign2D_prm::checkConvergence(std::vector<double> &conv)
{

    bool converged = true;
    double convv;
    Matrix2D<double> Maux;

    Maux.resize(dim, dim);
    Maux.setXmippOrigin();

    conv.clear();
    for (int iref = 0; iref < nr_ref; iref++)
    {
        if (Iref[iref].weight() > 0.)
        {
            Maux = multiplyElements(Iold[iref](), Iold[iref]());
            convv = 1. / (Maux.compute_avg());
            Maux = Iold[iref]() - Iref[iref]();
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

void Prog_MLPalign2D_prm::writeOutputFiles(const int iter, DocFile &DFo,
        double &sumw_allrefs, double &LL, double &avecorr,
        std::vector<double> &conv)
{

    FileName          fn_tmp, fn_base, fn_tmp2;
    Matrix1D<double>  fracline(3);
    SelFile           SFo, SFc;
    DocFile           DFl;
    string            comment;
    std::ofstream          fh;

    DFl.clear();
    SFo.clear();
    SFc.clear();

    fn_base = fn_root;
    if (iter >= 0)
    {
        fn_base += "_it";
        fn_base.compose(fn_base, iter, "");
    }

    // Write out current reference images and fill sel & log-file
    for (int iref = 0; iref < nr_ref; iref++)
    {
        fn_tmp = fn_base + "_ref";
        fn_tmp.compose(fn_tmp, iref + 1, "");
        fn_tmp = fn_tmp + ".xmp";
        Iref[iref].write(fn_tmp);
        SFo.insert(fn_tmp, SelLine::ACTIVE);
        fracline(0) = alpha_k[iref];
        fracline(1) = mirror_fraction[iref];
        fracline(2) = 1000 * conv[iref]; // Output 1000x the change for precision
        DFl.insert_comment(fn_tmp);
        DFl.insert_data_line(fracline);
    }

    // Write out sel & log-file
    fn_tmp = fn_base + ".sel";
    SFo.write(fn_tmp);

    DFl.go_beginning();
    comment = "MLalign2D-logfile: Number of images= " + floatToString(sumw_allrefs);
    comment += " LL= " + floatToString(LL, 15, 10) + " <Pmax/sumP>= " + floatToString(avecorr, 10, 5);
    DFl.insert_comment(comment);
    comment = "-noise " + floatToString(sigma_noise, 15, 12) + " -offset " + floatToString(sigma_offset, 15, 12) + " -istart " + integerToString(iter + 1);

    DFl.insert_comment(comment);
    DFl.insert_comment(cline);
    DFl.insert_comment("columns: model fraction (1); mirror fraction (2); 1000x signal change (3)");
    fn_tmp = fn_base + ".log";
    DFl.write(fn_tmp);

    // Write out docfile with optimal transformation & references
    fn_tmp = fn_base + ".doc";
    DFo.write(fn_tmp);

    // Also write out selfiles of all experimental images,
    // classified according to optimal reference image
    for (int iref = 0; iref < nr_ref; iref++)
    {
	DFo.go_beginning();
	SFo.clear();
	for (int n = 0; n < DFo.dataLineNo(); n++)
	{
	    DFo.next();
	    fn_tmp = ((DFo.get_current_line()).get_text()).erase(0, 3);
	    DFo.adjust_to_data_line();
	    if ((iref + 1) == (int)DFo(5)) SFo.insert(fn_tmp, SelLine::ACTIVE);
	}
	fn_tmp = fn_root + "_ref";
	fn_tmp.compose(fn_tmp, iref + 1, "sel");
	SFo.write(fn_tmp);
    }

}

