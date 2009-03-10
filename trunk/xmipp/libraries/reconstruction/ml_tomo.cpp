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

    // Imputation
    do_impute = !checkParameter(argc2, argv2, "-dont_impute");
    noimp_threshold = textToFloat(getParameter(argc2, argv2, "-noimp_threshold", "1."));

    // Angular sampling
    angular_sampling = textToFloat(getParameter(argc2, argv2, "-ang", "10"));
    psi_sampling = textToFloat(getParameter(argc2, argv2, "-psi_sampling", "-1"));
    tilt_range0 = textToFloat(getParameter(argc2, argv2, "-tilt0", "-91."));
    tilt_rangeF = textToFloat(getParameter(argc2, argv2, "-tiltF", "91."));
    //fn_sym = getParameter(argc2, argv2, "-sym", "c1");
    fn_sym="c1";

    // Missing data structures 
    do_wedge = checkParameter(argc2, argv2, "-wedges");
    do_pyramid = checkParameter(argc2, argv2, "-pyramids");
    do_cone = checkParameter(argc2, argv2, "-cones");
    if (do_wedge) 
        fn_missing = getParameter(argc2, argv2, "-wedges");
    else if (do_pyramid) 
        fn_missing = getParameter(argc2, argv2, "-pyramids");
    else if (do_cone) 
        fn_missing = getParameter(argc2, argv2, "-cones");
    max_resol = textToFloat(getParameter(argc2, argv2, "-max_resol", "0.5"));

    // Hidden arguments
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
        std::cerr << " |  *** Please cite them if this program is of use to you! ***   |" << std::endl;
        std::cerr << " -----------------------------------------------------------------" << std::endl;
        std::cerr << "--> Maximum-likelihood multi-reference refinement " << std::endl;
        std::cerr << "  Input images            : " << fn_sel << " (" << nr_exp_images << ")" << std::endl;
        if (fn_ref != "")
            std::cerr << "  Reference image(s)      : " << fn_ref << std::endl;
        else
            std::cerr << "  Number of references:   : " << nr_ref << std::endl;
        std::cerr << "  Output rootname         : " << fn_root << std::endl;
        std::cerr << "  Angular sampling rate   : " << angular_sampling<<std::endl;
        std::cerr << "  Stopping criterium      : " << eps << std::endl;
        std::cerr << "  initial sigma noise     : " << sigma_noise << std::endl;
        std::cerr << "  initial sigma offset    : " << sigma_offset << std::endl;
        if (do_missing)
        {
            std::cerr << "  Maximum resolution      : " << max_resol << std::endl;
            if (do_impute)
                std::cerr << "  Use imputation for data with missing regions "<<std::endl;
            else
                std::cerr << "  Divide by number of observations for data with missing regions "<<std::endl;
        }
        if (fn_frac != "")
            std::cerr << "  Initial model fractions : " << fn_frac << std::endl;
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
    XmippSampling mysampling;
    mysampling.SetSampling(angular_sampling);
    if (!mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
        REPORT_ERROR(3005, (std::string)"ml_refine3d::run Invalid symmetry" +  fn_sym);
    // by default max_tilt= +91., min_tilt= -91.
    mysampling.Compute_sampling_points(false,tilt_rangeF,tilt_range0);
    mysampling.remove_redundant_points(symmetry, sym_order);
    if (psi_sampling < 0)
        psi_sampling = angular_sampling;

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
            myinfo.A = Euler_rotation3DMatrix(rot, tilt, psi);
            all_angle_info.push_back(myinfo);
            nr_ang ++;
        }
    }

//#define DEBUG_SAMPLING
#ifdef DEBUG_SAMPLING
    for (int angno = 0; angno < nr_ang; angno++)
    {
        double rot=all_angle_info[angno].rot;
        double tilt=all_angle_info[angno].tilt;
        double psi=all_angle_info[angno].psi;
        std::cerr<<" rot= "<<rot<<" tilt= "<<tilt<<" psi= "<<psi<<std::endl;
    }
#endif


    // Read in docfile with information about the missing wedges
    nr_miss = 0;
    miss_thx0.clear();
    miss_thxF.clear();
    miss_thy0.clear();
    miss_thyF.clear();
    do_missing = (do_wedge || do_pyramid || do_cone);
    if (do_missing)
    {
        DocFile DFm;
        DFm.read(fn_missing);
        DFm.go_first_data_line();
        while (!DFm.eof())
        {
            if (do_wedge)
            {
                miss_thy0.push_back(DFm(0));
                miss_thyF.push_back(DFm(1));
            }
            else if (do_pyramid)
            {
                miss_thy0.push_back(DFm(0));
                miss_thyF.push_back(DFm(1));
                miss_thx0.push_back(DFm(2));
                miss_thxF.push_back(DFm(3));
            }
            else if (do_cone)
            {
                miss_thy0.push_back(DFm(0));
            }
            else
                REPORT_ERROR(1,"BUG: resolved type of missing data structure!");
            DFm.next_data_line();
            nr_miss++;
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
            alpha_k[refno] = DL[0];
            if (do_norm)
            {
                refs_avgscale[refno] = DL[3];
            }
            sumfrac += alpha_k[refno];
            DF.next_data_line();
        }
        if (ABS(sumfrac - 1.) > 1e-3)
            if (verb > 0) std::cerr << " ->WARNING: Sum of all expected model fractions (" << sumfrac << ") is not one!" << std::endl;
        for (int refno = 0; refno < nr_ref; refno++)
        {
            alpha_k[refno] /= sumfrac;
        }
    }

#define DEBUG_GENERAL
#ifdef DEBUG_GENERAL
    std::cerr<<"do_generate_refs ="<<do_generate_refs<<std::endl;
    std::cerr<<"nr_ref= "<<nr_ref<<std::endl; 
    std::cerr<<"nr_miss= "<<nr_miss<<std::endl;
    std::cerr<<"dim= "<<dim<<std::endl;
    std::cerr<<"nr_ang= "<<nr_ang<<std::endl;
    std::cerr<<"nr_psi= "<<nr_psi<<std::endl;
    std::cerr<<"Finished produceSideInfo"<<std::endl;
#endif
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
        init_progress_bar(nr_ref);
    }

    if (fn_doc!="") {
        DF.read(fn_doc);
    }
    // Make random subsets and calculate average images in random orientations
    // and correct using conventional division by the sum of the wedges
    randomize_random_generator();
    SFtmp = SF.randomize();
    int Nsub = ROUND((double)SFtmp.ImgNo() / nr_ref);
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
            int iran= ROUND(rnd_unif(0., (float)(nr_ang-1)));
            Itmp().selfApplyGeometry( all_angle_info[iran].A, IS_NOT_INV, DONT_WRAP);

            // This is dirty coding...
            // but let's first make it work...
            if (do_wedge || do_pyramid || do_cone)
            {
                if (DF.search_comment(fn_tmp)) 
                {
                    int missno = (int)(DF(0)) - 1;
                    getMissingWedge(Mmissing, all_angle_info[iran].A, missno);
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
                if (do_wedge || do_pyramid || do_cone) Msumwedge = Mmissing;
            }
            else
            {
                Iave() += Itmp();
                if (do_wedge || do_pyramid || do_cone) Msumwedge += Mmissing;
            }
        }
        // Correct missing wedge by division by sumwedge in Fourier space
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
        fn_tmp = fn_root + "_it";
        fn_tmp.compose(fn_tmp, 0, "");
        fn_tmp = fn_tmp + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "");
        fn_tmp = fn_tmp + ".vol";
        Iave.write(fn_tmp);
        SFr.insert(fn_tmp, SelLine::ACTIVE);
        if (verb > 0) progress_bar(refno);
    }
    if (verb > 0) progress_bar(nr_ref);
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
    DocFile                   DF;
    FileName                  fn_tmp;
    VolumeXmipp               img;
    std::vector<Matrix1D<double> > Vdum;

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

        // Rotate some arbitrary (off-axis) angle and rotate back again to remove high frequencies
        // that will be affected by the interpolation due to rotation
        // This makes that the A2 values of the rotated references are much less sensitive to rotation
        img().selfApplyGeometry( Euler_rotation3DMatrix(32., 61., 53.), IS_NOT_INV, DONT_WRAP);
        img().selfApplyGeometry( Euler_rotation3DMatrix(32., 61., 53.), IS_INV, DONT_WRAP);

        Iref.push_back(img);
        Iold.push_back(img);
        // Default start is all equal model fractions
        alpha_k.push_back((double)1 / SFr.ImgNo());
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
        while (!SF.eof()) 
        {
            fn_tmp=SF.NextImg();
            if (fn_tmp=="") break;
            DF.go_beginning();
            if (DF.search_comment(fn_tmp)) 
            {
                // TODO: If the only thing I want to get from the docfile is the missno, then
                // at a later stage move this if upwards to save time!
                if (do_wedge || do_pyramid || do_cone)
                {
                    imgs_missno[imgno] = (int)(DF(0)) - 1;
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
    } 


#ifdef  DEBUG
    std::cerr<<"Finished produceSideInfo2"<<std::endl;
#endif

}

void Prog_ml_tomo_prm::getMissingWedge(Matrix3D<double> &Mmissing,
                                       Matrix2D<double> A,
                                       const int missno)
{

/*
#define NEVER
#ifdef NEVER

    std::cerr<<"in"<<std::endl;
    Matrix3D<double> mask(dim,dim,dim);
    mask.setXmippOrigin();
    double xp, yp, zp;
    double tg0, tgF, limx0, limxF;

    tg0 = -tan(PI * (-90. - miss_thyF[missno]) / 180.);
    tgF = -tan(PI * (90. - miss_thy0[missno]) / 180.);
    std::cerr<<"old tg0= "<<tg0<<" tgF= "<<tgF<<std::endl;

    A = A.inv();
    FOR_ALL_ELEMENTS_IN_MATRIX3D(mask)
    {
        xp = dMij(A, 0, 0) * (double)j + dMij(A, 0, 1) * (double)i + dMij(A, 0, 2) * (double)k;
        zp = dMij(A, 2, 0) * (double)j + dMij(A, 2, 1) * (double)i + dMij(A, 2, 2) * (double)k;
        limx0 = tg0 * zp;
        limxF = tgF * zp;
        if (zp >= 0)
        {
            if (xp <= limx0 || xp >= limxF)
                VOL_ELEM(mask, k, i, j) = 1.;
            else
                VOL_ELEM(mask, k, i, j) = 0.;
        }
        else
        {
            if (xp <= limxF || xp >= limx0)
                VOL_ELEM(mask, k, i, j) = 1.;
            else
                VOL_ELEM(mask, k, i, j) = 0.;
        }
    }
    VolumeXmipp Vt;
    Vt()=mask;
    Vt.write("mask.vol");


    tg0 = -tan(PI * miss_thy0[missno] / 180.);
    tgF = -tan(PI * miss_thyF[missno] / 180.);
    std::cerr<<"tg0= "<<tg0<<" tgF= "<<tgF<<std::endl;

    A = A.inv();
    FOR_ALL_ELEMENTS_IN_MATRIX3D(mask)
    {
        xp = dMij(A, 0, 0) * (double)j + dMij(A, 0, 1) * (double)i + dMij(A, 0, 2) * (double)k;
        zp = dMij(A, 2, 0) * (double)j + dMij(A, 2, 1) * (double)i + dMij(A, 2, 2) * (double)k;
        if (xp >= 0)
        {
            if (zp >= tgF * xp && zp <= tg0 * xp)
               VOL_ELEM(mask, k, i, j) = 1.;
            else
               VOL_ELEM(mask, k, i, j) = 0.;
        }
        else
        {
            if (zp >= tg0 * xp && zp <= tgF * xp)
                VOL_ELEM(mask, k, i, j) = 1.;
            else
                VOL_ELEM(mask, k, i, j) = 0.;
        }

    }
    Vt()=mask;
    Vt.write("masknew.vol");


    exit(1);
#else
*/

    double theta0_alongy, thetaF_alongy, theta0_alongx, thetaF_alongx;
    Mmissing.resize(dim,dim,hdim+1);

    if (do_wedge) 
    {
        theta0_alongy=miss_thy0[missno];
        thetaF_alongy=miss_thyF[missno];
    }
    else if (do_pyramid)
    {
        REPORT_ERROR(1,"Missing pyramides not implemented yet...");
        theta0_alongy=miss_thy0[missno];
        thetaF_alongy=miss_thyF[missno];
        theta0_alongx=miss_thx0[missno];
        thetaF_alongx=miss_thxF[missno];
    }
    else if (do_cone)
    {
        REPORT_ERROR(1,"Missing cones not implemented yet...");
    }
    else
    {
        REPORT_ERROR(1,"bug: not do_wedge, nor do_pyramid, nor do_cone; but inside getMissingWedge");
    }

    Matrix2D<double> Ainv = A.inv();
    double xp, yp, zp;
    double tg0_y, tgF_y, tg0_x, tgF_x, limx0, limxF, limy0, limyF;

    tg0_y = -tan(PI * (-90. - thetaF_alongy) / 180.);
    tgF_y = -tan(PI * (90. - theta0_alongy) / 180.);

//#define DEBUG_WEDGE
#ifdef DEBUG_WEDGE
    std::cerr<<"tg0_y= "<<tg0_y<<std::endl;
    std::cerr<<"tgF_y= "<<tgF_y<<std::endl;
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
                    zp = dMij(Ainv, 2, 0) * xx + dMij(Ainv, 2, 1) * yy + dMij(Ainv, 2, 2) * zz;
                    // Calculate the limits
                    limx0 = tg0_y * zp;
                    limxF = tgF_y * zp;
                    if (zp >= 0)
                    {
                        if (xp <= limx0 || xp >= limxF)
                            DIRECT_MULTIDIM_ELEM(Mmissing,ii) = maskvalue;
                        else
                            DIRECT_MULTIDIM_ELEM(Mmissing,ii) = 0.;
                    }
                    else
                    {
                        if (xp <= limxF || xp >= limx0)
                            DIRECT_MULTIDIM_ELEM(Mmissing,ii) = maskvalue;
                        else
                            DIRECT_MULTIDIM_ELEM(Mmissing,ii) = 0.;
                    }
                }
            }
        }
    }
//#define DEBUG_WEDGE
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
    Matrix1D<double> off(3);
    off.initConstant(hdim);
    test().selfTranslate(off);
    test.write("Fwedge.vol");
    test().setXmippOrigin();
    Matrix3D<int> ress(dim,dim,dim);
    ress.setXmippOrigin();
    // TODO: CHECK A or Ainv ...??!!
    BinaryWedgeMask(test(),theta0_alongy,thetaF_alongy, A);
    test.write("Mwedge.vol");
#endif

}

// Calculate probability density function of all in-plane transformations phi
void Prog_ml_tomo_prm::calculatePdfInplane()
{

#ifdef  DEBUG
    std::cerr<<"start calculatePdfInplane"<<std::endl;
#endif

    double r2, pdfpix, sum;
    P_phi.resize(dim, dim, dim);
    P_phi.setXmippOrigin();
    Mr2.resize(dim, dim, dim);
    Mr2.setXmippOrigin();

    sum=0.;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(P_phi)
    {
        r2 = (double)(j * j + i * i + k * k);
        if (sigma_offset > 0.)
        {
            pdfpix = exp(-r2 / (2 * sigma_offset * sigma_offset));
            pdfpix /= 2 * PI * sigma_offset * sigma_offset * nr_ang;
        }
        else
        {
            if (j == 0 && i == 0) pdfpix = 1.;
            else pdfpix = 0.;
        }
        VOL_ELEM(P_phi, k, i, j) = pdfpix;
        VOL_ELEM(Mr2, k, i, j) = (float)r2;
	sum+=pdfpix;
    }
    // Normalization
    P_phi/=sum;

#ifdef  DEBUG
    std::cerr<<"finished calculatePdfInplane"<<std::endl;
#endif

}

/// Mask references 
void Prog_ml_tomo_prm::maskReferences(std::vector< VolumeXmippT<double> > &Iref)
{
    // This is now obsolete because afterwards real_mask is applied to the references
    /* 
    Matrix3D<int> mask, omask;
    double dum, avg;

    // prepare masks
    mask.resize(dim, dim, dim);
    mask.setXmippOrigin();
    BinarySphericalMask(mask, hdim, INNER_MASK);
    omask.resize(dim, dim, dim);
    omask.setXmippOrigin();
    BinarySphericalMask(omask, hdim, OUTSIDE_MASK);
    for (int refno = 0; refno < nr_ref; refno++)
    {
        computeStats_within_binary_mask(omask, Iref[refno](), dum, dum, avg, dum);
        apply_binary_mask(mask, Iref[refno](), Iref[refno](), avg);
    }
    */

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
            Maux.initZeros(); // This somehow seems necessary before applyGeometry!!
            applyGeometry(Maux, A_rot_inv, Iref[refno](), IS_NOT_INV, DONT_WRAP);
            Maux *= real_mask;
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
                    getMissingWedge(Mmissing,I,missno);
                    Faux2 = Faux;
                    double sumw =0.;
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux2)
                    {
                        DIRECT_MULTIDIM_ELEM(Faux2,n) *= DIRECT_MULTIDIM_ELEM(Mmissing,n);
                        sumw += DIRECT_MULTIDIM_ELEM(Mmissing,n);
                    }
                    transformer.inverseFourierTransform(Faux2,Maux);
                    A2.push_back(Maux.sum2());
//#define DEBUG_PRECALC_A2
#ifdef DEBUG_PRECALC_A2
                    std::cerr<<"rot= "<<all_angle_info[angno].rot<<" tilt= "<<all_angle_info[angno].tilt<<" psi= "<<all_angle_info[angno].psi<<std::endl;
                    std::cerr<<"refno= "<<refno<<" angno= "<<angno<<" missno= "<<missno<<" A2= "<<Maux.sum2()<<" corrA2= "<<corr<<" sumw= "<<sumw<<std::endl;
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
    Matrix3D<double> &Mimg, int &missno,
    std::vector< VolumeXmippT<double> > &Iref,
    std::vector<Matrix3D<double> > &wsumimgs,
    std::vector<Matrix3D<double> > &wsumweds,
    double &wsum_sigma_noise, double &wsum_sigma_offset,
    std::vector<double> &sumw, std::vector<double> &sumwsc, std::vector<double> &sumwsc2, 
    double &LL, double &dLL, double &fracweight, double &sumfracweight, 
    double &opt_scale, double &bgmean, double &trymindiff, 
    int &opt_refno, int &opt_angno, Matrix1D<double> &opt_offsets,
    std::vector<double> &pdf_directions)
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
    std::vector<double> refw(nr_ref), refwsc2(nr_ref);
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

    // Calculate Fimg0: the unrotated Fourier transform (with enforced wedge) of Mimg
    // and calculate myXi2;
    local_transformer.FourierTransform(Mimg,Fimg0,true);
    if (do_missing)
    {
        // Enforce missing wedge
        getMissingWedge(Mmissing,I,missno);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg0)
        {
            DIRECT_MULTIDIM_ELEM(Fimg0,n) *= DIRECT_MULTIDIM_ELEM(Mmissing,n);
        }
    }
    if (do_norm) dVkij(Fimg0,0,0,0) -= bgmean;
    local_transformer.inverseFourierTransform(Fimg0,Maux);
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
            refw[refno] = 0.;
        }
 
        // The real stuff: now loop over all orientations, references and translations
        // Start the loop over all refno at old_optangno (=opt_angno from the previous iteration). 
        // This will speed-up things because we will find Pmax probably right away,
        // and this will make the if-statement that checks SIGNIFICANT_WEIGHT_LOW
        // effective right from the start
        for (int aa = old_optangno; aa < old_optangno+nr_ang; aa++)
        {
            int angno = aa;
            if (angno >= nr_ang) angno -= nr_ang;

            // Precalculate rotated version of Mimg for its storage in the weighted sums
            // At this point, no longer worry about the interpolation effects on the missing wedges
            A_rot = (all_angle_info[angno]).A;
            A_rot_inv = A_rot.inv();
            Maux.setXmippOrigin();
            Maux.initZeros();
            applyGeometry(Maux, A_rot, Mimg, IS_NOT_INV, WRAP);
            // From here on local_transformer will act on Maux and Faux
            local_transformer.FourierTransform(Maux,Faux,false);
            Fimg_rot=Faux;

            // Start the loop over all refno at old_optrefno (=opt_refno from the previous iteration). 
            // This will speed-up things because we will find Pmax probably right away,
            // and this will make the if-statement that checks SIGNIFICANT_WEIGHT_LOW
            // effective right from the start
            for (int rr = old_optrefno; rr < old_optrefno+nr_ref; rr++)
            {
                int refno = rr;
                if (refno >= nr_ref) refno-= nr_ref;
                
                // Now (inverse) rotate the reference and calculate its Fourier transform
                // TODO: check the inverse rotation and its A2!!
                // Do that by reenforcing the wedge, and recalculating A2
                
                Maux2.initZeros();
                applyGeometry(Maux2, A_rot_inv, Iref[refno](), IS_NOT_INV, DONT_WRAP);
                Maux2 *= real_mask;
                mycorrAA = corrA2[refno*nr_ang + angno];
                Maux = Maux2 * mycorrAA;

//#define DEBUG_MINDIFF
#ifdef DEBUG_MINDIFF
                Matrix3D<double> Maux_ori2=Maux2;
                Matrix3D<double> Maux_ori=Iref[refno]();
                Matrix3D<double> Maux_ori3=Maux;
#endif
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
                Matrix3D<std::complex<double> > Faux_ori=Faux;
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
                {
                    dVkij(Faux,k,i,j) = 
                        dVkij(Fimg0,k,i,j) * 
                        conj(dVkij(Faux,k,i,j));
                }
                // Takes the input from Faux, and leaves the output in Maux
                local_transformer.inverseFourierTransform();
                CenterFFT(Maux, true);

                // B. Calculate weights for each pixel within sigdim (Mweight)
                my_sumweight = my_maxweight = 0.;
//#define DEBUG_ALOT_SINGLEEXP
#ifdef DEBUG_ALOT_SINGLEEXP
                std::cerr<<"rot= "<<all_angle_info[angno].rot<<" tilt= "<<all_angle_info[angno].tilt<<" psi= "<<all_angle_info[angno].psi<<" A2= "<<myA2<<" Xi2= "<<myXi2<<" corrA="<<mycorrAA <<" XA= "<<VOL_ELEM(Maux, 0,0,0) * ddim3<<" diff= "<< A2_plus_Xi2 - ref_scale * VOL_ELEM(Maux, 0, 0, 0) * ddim3<<" mindiff= "<<mindiff<<std::endl;
#endif
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
                        
                        local_transformer.inverseFourierTransform();
                        tt()=Maux; tt.write("Faux.vol");
                        Faux=Faux_ori;
                        local_transformer.inverseFourierTransform();
                        tt()=Maux; tt.write("Faux_ori.vol");
                        tt()=Maux_ori; tt.write("Maux_ori.vol");
                        tt()=Maux_ori2; tt.write("Maux_ori2.vol");
                        tt()=Maux_ori3; tt.write("Maux_ori3.vol");
                        std::cerr<<"Ainv= "<<A_rot_inv<<std::endl;
                        std::cerr<<"A= "<<A_rot_inv.inv()<<std::endl;

                        exit(1);
                    }
#endif
                    pdf = alpha_k[refno] * VOL_ELEM(P_phi, k, i, j);
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
                    refw[refno] += my_sumweight;
                        
                    // Back from smaller Mweight to original size of Maux
                    Maux.initZeros();
                    FOR_ALL_ELEMENTS_IN_MATRIX3D(Mweight)
                    {
                        VOL_ELEM(Maux, k, i, j) = VOL_ELEM(Mweight, k, i, j);
                    }
                    // Use forward FFT in convolution theorem again
                    // Takes the input from Maux and leaves it in Faux
                    CenterFFT(Maux, false);
                    local_transformer.FourierTransform();
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
                    {
                        dVkij(Faux, k, i, j) = conj(dVkij(Faux,k,i,j)) * dVkij(Fimg0,k,i,j);
                    }
                    local_transformer.inverseFourierTransform();

                    Maux.selfApplyGeometry(A_rot, IS_NOT_INV, DONT_WRAP);
                    mysumimgs[refno] += Maux * ddim3;
                    if (do_missing)
                    {
                        // Store sum of wedges!
                        // TODO: check A_rot or A_rot_inv...
                        getMissingWedge(Mmissing, A_rot, missno);
                        mysumweds[refno] += my_sumweight * Mmissing;
                    }
                } // close if SIGNIFICANT_WEIGHT_LOW
            } // close for refno
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

    // Calculate optimal transformation parameters
    A_rot = all_angle_info[opt_angno].A;
    A_rot = A_rot.inv();
    /*
    XX(opt_offsets) = -(double)ioptx * dMij(A_rot, 0, 0) 
                      -(double)iopty * dMij(A_rot, 0, 1)
                      -(double)ioptz * dMij(A_rot, 0, 2);
    YY(opt_offsets) = -(double)ioptx * dMij(A_rot, 1, 0) 
                      -(double)iopty * dMij(A_rot, 1, 1)
                      -(double)ioptz * dMij(A_rot, 1, 2);
    ZZ(opt_offsets) = -(double)ioptx * dMij(A_rot, 2, 0) 
                      -(double)iopty * dMij(A_rot, 2, 1)
                      -(double)ioptz * dMij(A_rot, 2, 2);
    */
    XX(opt_offsets) = -(double)ioptx;
    YY(opt_offsets) = -(double)iopty;
    ZZ(opt_offsets) = -(double)ioptz;

//#define DEBUG_TRANSFORMATIONS
#ifdef DEBUG_TRANSFORMATIONS
    std::cerr<<"exp mindiff= "<<mindiff<<" opt angno= "<<opt_angno<<" "<< all_angle_info[opt_angno].rot<<" "<< all_angle_info[opt_angno].tilt<<" "<< all_angle_info[opt_angno].psi<<std::endl;
    std::cerr<<"opt_ioffs= "<<ioptx<<" "<<iopty<<" "<<ioptz<<std::endl;
    std::cerr<<"opt_offsets= "<<opt_offsets<<std::endl;
#endif

    // Update normalization parameters
    if (do_norm)
    {
        // 1. Calculate optimal setting of Mimg
        Maux2 = Mimg;
        Maux2.selfTranslate(opt_offsets, true);

        // 2. Calculate optimal setting of Iref[opt_refno]
        A_rot_inv = (all_angle_info[opt_angno].A).inv();
        Maux.initZeros();
        applyGeometry(Maux, A_rot_inv, Iref[opt_refno](), IS_NOT_INV, WRAP);
        Maux *= opt_scale;

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

    // From here on lock threads
    pthread_mutex_lock( &mltomo_weightedsum_update_mutex );

    // Update all global weighted sums after division by sum_refw
    wsum_sigma_noise += (2 * wsum_corr / sum_refw);
    wsum_sigma_offset += (wsum_offset / sum_refw);
    sumfracweight += fracweight;
    for (int refno = 0; refno < nr_ref; refno++)
    {
        //if (!limit_rot || pdf_directions[refno] > 0.)
        sumw[refno] += refw[refno] / sum_refw;
        sumwsc[refno] += (refw[refno] * opt_scale) / sum_refw;
        sumwsc2[refno] += (refw[refno] * opt_scale * opt_scale) / sum_refw;
        // Correct weighted sum of images for new bgmean (only first element=origin in Fimg)
        if (do_norm)
            // TODO: check factor ddim3!
            //mysumimgs[refno] -= refw[refno] * (bgmean - old_bgmean) / ddim3; 
            mysumimgs[refno] -= refw[refno] * (bgmean - old_bgmean); 
        // Sum mysumimgs to the global weighted sum
        wsumimgs[refno] += (opt_scale * mysumimgs[refno]) / sum_refw;
        if (do_missing)
        {
            wsumweds[refno] += mysumweds[refno] / sum_refw;
        }
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
    std::vector<double> *sumw = thread_data->sumw;
    std::vector<double> *sumwsc = thread_data->sumwsc;
    std::vector<double> *sumwsc2 = thread_data->sumwsc2;

//#define DEBUG_THREAD
#ifdef DEBUG_THREAD
    std::cerr<<"start threadMLTomoExpectationSingleImage"<<std::endl;
#endif

    // Local variables
    VolumeXmipp img;
    FileName fn_img, fn_trans;
    std::vector<double> pdf_directions(prm->nr_ref);
    std::vector<Matrix1D<double> > allref_offsets;
    Matrix1D<double> opt_offsets(3);
    float old_phi = -999., old_theta = -999.;
    double fracweight, maxweight2, trymindiff, dLL;
    double opt_scale = 1., bgmean = 0.;
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
        
        // These three parameters speed up expectationSingleImage
        trymindiff = prm->imgs_trymindiff[imgno];
        opt_refno = prm->imgs_optrefno[imgno];
        opt_angno = prm->imgs_optangno[imgno];
        if (prm->do_wedge || prm->do_pyramid || prm->do_cone)
            missno = prm->imgs_missno[imgno];
        else
            missno = -1;

        if (prm->do_norm)
        {
            bgmean=prm->imgs_bgmean[imgno];
            opt_scale = prm->imgs_scale[imgno];
        }
            
        (*prm).expectationSingleImage(img(), missno, *Iref, *wsumimgs, *wsumweds, 
                                      *wsum_sigma_noise, *wsum_sigma_offset, 
                                      *sumw, *sumwsc, *sumwsc2, 
                                      *LL, dLL, fracweight, *sumfracweight, 
                                      opt_scale, bgmean, trymindiff, opt_refno, opt_angno, 
                                      opt_offsets,pdf_directions);

        // Store mindiff for next iteration
        prm->imgs_trymindiff[imgno] = trymindiff;
        // Store opt_refno for next iteration
        prm->imgs_optrefno[imgno] = opt_refno;
        // Store opt_angno for next iteration
        prm->imgs_optangno[imgno] = opt_angno;

        // Store optimal normalization parameters in memory
        if (prm->do_norm)
        {
            prm->imgs_scale[imgno] = opt_scale;
            prm->imgs_bgmean[imgno]  = bgmean;
        }
        
        // Output docfile
        (*docfiledata)[imgno](0) = (prm->all_angle_info[opt_angno]).rot;// rot
        (*docfiledata)[imgno](1) = (prm->all_angle_info[opt_angno]).tilt;// tilt
        (*docfiledata)[imgno](2) = (prm->all_angle_info[opt_angno]).psi;// psi
        (*docfiledata)[imgno](3) = opt_offsets(0);            // Xoff
        (*docfiledata)[imgno](4) = opt_offsets(1);            // Yoff
        (*docfiledata)[imgno](5) = opt_offsets(2);            // Zoff
        (*docfiledata)[imgno](6) = (double)(opt_refno + 1);   // Ref
        (*docfiledata)[imgno](7) = fracweight;                // P_max/P_tot
        (*docfiledata)[imgno](8) = dLL;                       // log-likelihood
        if (prm->do_norm)
        {
            (*docfiledata)[imgno](9)  = bgmean;               // background mean
            (*docfiledata)[imgno](10) = opt_scale;            // image scale 
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
	std::vector<double> &sumw, std::vector<double> &sumwsc, std::vector<double> &sumwsc2)
{

#ifdef DEBUG
    std::cerr<<"start expectation"<<std::endl;
#endif 

    Matrix1D<double> dataline(MLTOMODATALINELENGTH);
    Matrix3D<double> Mzero(dim,dim,hdim+1), Mzero2(dim,dim,dim);
    std::vector<Matrix1D<double> > docfiledata;
    bool fill_real_space;
    int num_img_tot;
    
    // Precalculate A2-values for all references
    precalculateA2(Iref);

    // Pre-calculate pdf of all in-plane transformations
    calculatePdfInplane();

    // Initialize weighted sums
    LL = 0.;
    sumw.clear();
    sumwsc.clear();
    sumwsc2.clear();
    wsum_sigma_noise = 0.;
    wsum_sigma_offset = 0.;
    sumfracweight = 0.;
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
        sumw.push_back(0.);
        sumwsc.push_back(0.);
        sumwsc2.push_back(0.);
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
                                    std::vector<double> &sumw, std::vector<double> &sumwsc, 
                                    std::vector<double> &sumwsc2,
                                    double &sumfracweight, double &sumw_allrefs)
{
#ifdef DEBUG
    std::cerr<<"started maximization"<<std::endl;
#endif 

    Matrix1D<double> rmean_sigma2, rmean_signal2;
    Matrix1D<int> center(3), radial_count;
    Matrix3D<std::complex<double> > Faux(dim,dim,hdim+1), Fwsumimgs;
    Matrix3D<double> Maux, Msumallwedges(dim,dim,hdim+1);
    FileName fn_tmp;
    double rr, thresh, aux, sumw_allrefs2 = 0.;
    int c;

    // Update the reference images
    sumw_allrefs = 0.;
    Msumallwedges.initZeros();
    for (int refno=0;refno<nr_ref; refno++)
    {
        sumw_allrefs += sumw[refno];
        transformer.FourierTransform(wsumimgs[refno],Fwsumimgs,true);
        if (do_missing)
        {
            Msumallwedges += wsumweds[refno];
            if (do_impute)
            {
                transformer.FourierTransform(Iref[refno](),Faux,true);
                if (sumw[refno] > 0.)
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
                            (1. - DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) / sumw[refno]);
                        DIRECT_MULTIDIM_ELEM(Faux,n) *= 
                            (1. - DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) / sumw[refno]);
                    }
                    fnt.compose("inv_wsum",refno,"ampl");
                    tt.write(fnt);
                    FFT_magnitude(Faux,tt());
                    fnt.compose("Fref1",refno,"ampl");
                    tt.write(fnt);
                    FFT_magnitude(Fwsumimgs,tt());
                    fnt.compose("Fwsumimgs",refno,"ampl");
                    tt.write(fnt);
                    std::cerr<<" sumw[refno]= "<<sumw[refno]<<" sumwsc2[refno]= "<<sumwsc2[refno]<<std::endl;
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                    {
                        // And sum the weighted sum for observed pixels
                        DIRECT_MULTIDIM_ELEM(Faux,n) += 
                            (DIRECT_MULTIDIM_ELEM(Fwsumimgs,n) / sumwsc2[refno]);
                    }
                    FFT_magnitude(Faux,tt());
                    fnt.compose("Fref2",refno,"ampl");
                    tt.write(fnt);
#else
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                    {
                        // Impute old reference for missing pixels
                        DIRECT_MULTIDIM_ELEM(Faux,n) *= 
                            (1. - DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) / sumw[refno]);
                        // And sum the weighted sum for observed pixels
                        DIRECT_MULTIDIM_ELEM(Faux,n) += 
                            (DIRECT_MULTIDIM_ELEM(Fwsumimgs,n) / sumwsc2[refno]);
                    }
#endif
                }
                // else do nothing (i.e. impute old reference completely
            }
            else // no imputation: divide by number of times a pixel has been observed 
            {
                std::cerr<< " "<<std::endl;
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                {
                    if (DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) > noimp_threshold)
                    {
                        DIRECT_MULTIDIM_ELEM(Faux,n) = 
                            DIRECT_MULTIDIM_ELEM(Fwsumimgs,n) / 
                            DIRECT_MULTIDIM_ELEM(wsumweds[refno],n);
                        // TODO:  CHECK THE FOLLOWING LINE!!
                        DIRECT_MULTIDIM_ELEM(Faux,n) *= sumw[refno] / sumwsc2[refno];
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
            // if no missing data: calculate straightforward average, 
            // i.e. just divide wsumimgs[refno] by sumwsc2[refno]
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
            {
                DIRECT_MULTIDIM_ELEM(Faux,n) = 
                    DIRECT_MULTIDIM_ELEM(Fwsumimgs,n) / sumwsc2[refno];
            }
        }
        transformer.inverseFourierTransform(Faux,Iref[refno]());

    }

    // Adjust average scale (nr_classes will be smaller than nr_ref for the 3D case!)
    if (do_norm) {
        std::vector<double> wsum_scale(nr_ref), sumw_scale(nr_ref);
        average_scale = 0.;
        for (int refno=0;refno<nr_ref; refno++)
        {
            average_scale += sumwsc[refno];
            wsum_scale[refno] += sumwsc[refno];
            sumw_scale[refno] += sumw[refno];
        }
        for (int refno=0;refno<nr_ref; refno++)
        {
            if (sumw_scale[refno]>0.)
            {
                refs_avgscale[refno] = wsum_scale[refno]/sumw_scale[refno];
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
            if (sumw[refno] > 0.)
                alpha_k[refno] = sumw[refno] / sumw_allrefs;
            else
                alpha_k[refno] = 0.;
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

#ifdef DEBUG
    std::cerr<<"finished maximization"<<std::endl;
#endif 
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
        if (alpha_k[refno] > 0.)
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
        fracline(0) = alpha_k[refno];
        fracline(1) = 1000 * conv[refno]; // Output 1000x the change for precision
        DFl.insert_comment(fn_tmp);
        DFl.insert_data_line(fracline);

        if (iter >= 1 && do_missing)
        {

            double sumw = alpha_k[refno]*sumw_allrefs;
            Vt().resize(dim,dim,dim);
            Vt().setXmippOrigin();
            Vt().initZeros();
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
                        if (xx <= FINISHINGX(Vt()))
                        {
                            VOL_ELEM(Vt(),zz,yy,xx)=
                                DIRECT_MULTIDIM_ELEM(wsumweds[refno],ii) / sumw;
                            VOL_ELEM(Vt(),zz,yy,-xx)=
                                DIRECT_MULTIDIM_ELEM(wsumweds[refno],ii) / sumw;
                        }
                    } 
                }
            }
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
    for (int refno = 0;refno < nr_ref; refno++)
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

