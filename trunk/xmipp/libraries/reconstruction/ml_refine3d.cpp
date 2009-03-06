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

#include "ml_refine3d.h"


// Read ===================================================================
void Prog_Refine3d_prm::read(int argc, char ** argv, int &argc2, char ** &argv2)
{

    bool do_restart = false;

    if (checkParameter(argc, argv, "-more_options"))
    {
        if (fourier_mode)
        {
            MLF_usage();
            extended_usage();
        }
        else
        {
            usage();
            extended_usage();
        }
    }
    if (checkParameter(argc, argv, "-show_all_ML_options"))
    {
	if (fourier_mode)

	{
	    Prog_MLFalign2D_prm MLF_prm;
	    MLF_prm.extendedUsage(true);
	}
	else
	{
	    Prog_MLalign2D_prm ML_prm;
	    ML_prm.extendedUsage(true);
	}
    }
    if (checkParameter(argc, argv, "-show_all_ART_options"))
    {
        Basic_ART_Parameters   art_prm;
        art_prm.usage_more();
        exit(0);
    }

    // Generate new command line for restart procedure
    if (checkParameter(argc, argv, "-restart"))
    {
        std::string   comment, cline = "";
        DocFile  DFi;
        FileName fn_tmp;

        do_restart = true;
        DFi.read(getParameter(argc, argv, "-restart"));
        DFi.go_beginning();
        comment = (DFi.get_current_line()).get_text();
        if (fourier_mode && strstr(comment.c_str(), "MLFalign2D-logfile") == NULL)
        {
            std::cerr << "Error!! Docfile is not of MLFalign2D-logfile type. " << std::endl;
            exit(1);
        }
        else if (!fourier_mode && strstr(comment.c_str(), "MLalign2D-logfile") == NULL)
        {
            std::cerr << "Error!! Docfile is not of MLalign2D-logfile type. " << std::endl;
            exit(1);
        }
        else
        {
            char *copy;
            copy = NULL;
            DFi.next();
            // get new part of command line (with -istart)
            comment = (DFi.get_current_line()).get_text();
            DFi.next();
            // get original command line
            cline = (DFi.get_current_line()).get_text();
            comment = comment + cline;
            // regenerate command line
	    argv2 = NULL;
	    argc2 = 0;
            generateCommandLine(comment, argc2, argv2, copy);
            // Get number of volumes and names to generate SFvol
            if (fourier_mode) fn_root = getParameter(argc2, argv2, "-o", "mlf3d");
            else fn_root = getParameter(argc2, argv2, "-o", "ml3d");
            fn_vol = getParameter(argc2, argv2, "-vol");
            istart = textToInteger(getParameter(argc2, argv2, "-istart"));
            if (Is_VolumeXmipp(fn_vol))
            {
                SFvol.reserve(1);
                SFvol.insert(fn_vol);
            }
            else
            {
                SFvol.read(fn_vol);
            }
            Nvols = SFvol.ImgNo();
            SFvol.clear();
            SFvol.go_beginning();
            for (int ivol = 0; ivol < Nvols; ivol++)
            {
                fn_tmp = fn_root + "_it";
                fn_tmp.compose(fn_tmp, istart - 1, "");
                if (Nvols > 1)
                {
                    fn_tmp += "_vol";
                    fn_tmp.compose(fn_tmp, ivol + 1, "");
                }
                fn_tmp += ".vol";
                SFvol.insert(fn_tmp);
            }
            fn_vol = fn_root + "_it";
            fn_vol.compose(fn_vol, istart - 1, "");
            fn_vol += "_restart.sel";
            SFvol.write(fn_vol);
        }
    }
    else
    {
	// no restart, just copy argc to argc2 and argv to argv2
	argc2 = argc;
	argv2 = argv;
    }


    //Read Refine3d parameters
    fn_sel = getParameter(argc2, argv2, "-i");
    if (fourier_mode) fn_root = getParameter(argc2, argv2, "-o", "mlf3d");
    else fn_root = getParameter(argc2, argv2, "-o", "ml3d");
    if (!do_restart)
    {
        // Fill volume selfile
        fn_vol = getParameter(argc2, argv2, "-vol");
        if (Is_VolumeXmipp(fn_vol))
        {
            SFvol.reserve(1);
            SFvol.insert(fn_vol);
        }
        else
        {
            SFvol.read(fn_vol);
        }
        Nvols = SFvol.ImgNo();
    }

    angular = textToFloat(getParameter(argc2, argv2, "-ang", "10"));
    fn_sym = getParameter(argc2, argv2, "-sym", "c1");
    eps = textToFloat(getParameter(argc2, argv2, "-eps", "5e-5"));
    verb = textToInteger(getParameter(argc2, argv2, "-verb", "1"));
    Niter = textToInteger(getParameter(argc2, argv2, "-iter", "25"));
    istart = textToInteger(getParameter(argc2, argv2, "-istart", "1"));
    tilt_range0 = textToFloat(getParameter(argc2, argv2, "-tilt0", "-91."));
    tilt_rangeF = textToFloat(getParameter(argc2, argv2, "-tiltF", "91."));
    fn_symmask = getParameter(argc2, argv2, "-sym_mask", "");
    lowpass = textToFloat(getParameter(argc2, argv2, "-filter", "-1"));
    wlsart_no_start = checkParameter(argc2, argv2, "-nostart");
    do_perturb = checkParameter(argc2, argv2, "-perturb");

    // Hidden for now
    fn_solv = getParameter(argc2, argv2, "-solvent", "");
    reconstruct_wbp = checkParameter(argc2, argv2, "-WBP");
    reconstruct_fourier = checkParameter(argc2, argv2, "-fourier");
    do_prob_solvent = checkParameter(argc2, argv2, "-prob_solvent");
    threshold_solvent = textToFloat(getParameter(argc2, argv2, "-threshold_solvent", "999"));
    do_deblob_solvent = checkParameter(argc2, argv2, "-deblob_solvent");
    dilate_solvent = textToInteger(getParameter(argc2, argv2, "-dilate_solvent", "0"));
    skip_reconstruction = checkParameter(argc2, argv2, "-skip_reconstruction");

    // Checks
    if (lowpass > 0.5) REPORT_ERROR(1, "Digital frequency for low-pass filter should be smaller than 0.5");

}

// Usage ===================================================================
void Prog_Refine3d_prm::usage()
{
    std::cerr << "Usage:  ml_refine3d [options] " << std::endl;
    std::cerr << "   -i <selfile>                : Selfile with input images \n"
    << "   -vol <volume/selfile>       : Initial reference volume \n"
    << "                               :  OR selfile with multiple reference volumes\n"
    << " [ -o <root=\"ml3d\"> ]          : Output rootname \n"
    << " [ -ang <float=10> ]           : Angular sampling (degrees) \n"
    << " [ -iter <int=25> ]            : Maximum number of iterations \n"
    << " [ -more_options ]             : Show additional parameters for 3D-refinement\n";

}

// MLF Usage =================================================================
void Prog_Refine3d_prm::MLF_usage()
{
    std::cerr << "Usage:  mlf_refine3d [options] " << std::endl;
    std::cerr << "   -i <selfile>                : Selfile with all input images \n";
    std::cerr << "   -ctfdat <ctfdatfile>        : Two-column ASCII file with filenames and CTF parameter files of all images \n";
    std::cerr << "      OR -no_ctf                   OR do not use any CTF correction \n";
    std::cerr << "   -vol <volume/selfile>        : Initial reference volume \n";
    std::cerr << "                               :  OR selfile with multiple reference volumes\n";
    std::cerr << " [ -o <rootname> ]             : Output rootname (default = \"mlf2d\")\n";
    std::cerr << " [ -ang <float=10> ]           : Angular sampling (degrees) \n";
    std::cerr << " [ -iter <int=25>  ]           : Maximum number of iterations \n";
    std::cerr << " [ -search_shift <float=3>]    : Limited translational searches (in pixels) \n";
    std::cerr << " [ -reduce_noise <factor=1> ]  : Use a value smaller than one to decrease the estimated SSNRs \n";
    std::cerr << " [ -not_phase_flipped ]        : Use this if the experimental images have not been phase flipped \n";
    std::cerr << " [ -ctf_affected_refs ]        : Use this if the references (-ref) are not CTF-deconvoluted \n";
    std::cerr << " [ -low <Ang=999> ]            : Exclude lowest frequencies from P-calculations (in Ang) \n";
    std::cerr << " [ -high <Ang=0> ]             : Exclude highest frequencies from P-calculations (in Ang) \n";
    std::cerr << " [ -ini_high <Ang=0> ]         : Exclude highest frequencies during first iteration (in Ang) \n";
    std::cerr << " [ -pixel_size <Ang=1> ]       : Pixel size in Angstrom (only necessary for -no_ctf mode) \n";
    std::cerr << " [ -more_options ]             : Show additional parameters for 3D-refinement \n";

}

// Extended usage =============================================================
void Prog_Refine3d_prm::extended_usage()
{
    std::cerr << "Additional options: " << std::endl;
    std::cerr << " [ -l <float=0.2> ]            : wlsART-relaxation parameter (lambda)  \n"
    << " [ -k <float=0.5> ]            : wlsART-relaxation parameter for residual (kappa)\n"
    << " [ -n <int=10> ]               : Number of wlsART-iterations \n"
    << " [ -nostart ]                  : Start wlsART reconstructions from all-zero volumes \n"
    << " [ -sym <symfile> ]            : Enforce symmetry \n"
    << " [ -filter <dig.freq.=-1> ]    : Low-pass filter volume every iteration \n"
    << " [ -sym_mask <maskfile> ]      : Local symmetry (only inside mask) \n"
    << " [ -tilt0 <float=-91.> ]       : Lower-value for restricted tilt angle search \n"
    << " [ -tiltF <float=91.> ]        : Higher-value for restricted tilt angle search \n"
    << " [ -perturb ]                  : Randomly perturb reference projection directions \n"
    << " [ -show_all_ML_options ]      : Show all parameters for the ML-refinement\n"
    << " [ -show_all_ART_options ]     : Show all parameters for the wlsART reconstruction \n";
    std::cerr << std::endl;
    exit(1);
}

// Show ======================================================================
void Prog_Refine3d_prm::show()
{

    if (verb > 0)
    {
        // To screen
        std::cerr << " -----------------------------------------------------------------" << std::endl;
        std::cerr << " | Read more about this program in the following publication:    |" << std::endl;
        if (fourier_mode)
            std::cerr << " |  Scheres ea. (2007)  Structure, 15, 1167-1177                 |" << std::endl;
        else
            std::cerr << " |  Scheres ea. (2007)  Nature Methods, 4, 27-29                 |" << std::endl;
        std::cerr << " |                                                               |" << std::endl;
        std::cerr << " |    *** Please cite it if this program is of use to you! ***   |" << std::endl;
        std::cerr << " -----------------------------------------------------------------" << std::endl;
        std::cerr << "--> Maximum-likelihood multi-reference 3D-refinement" << std::endl;
        if (Nvols == 1)
            std::cerr << "  Initial reference volume : " << fn_vol << std::endl;
        else
        {
            std::cerr << "  Selfile with references  : " << fn_vol << std::endl;
            std::cerr << "    with # of volumes      : " << Nvols << std::endl;
        }
        std::cerr << "  Experimental images:     : " << fn_sel << std::endl;
        std::cerr << "  Angular sampling rate    : " << angular << std::endl;
        if (fn_sym != "")
            std::cerr << "  Symmetry file:           : " << fn_sym << std::endl;
        if (fn_symmask != "")
            std::cerr << "  Local symmetry mask      : " << fn_symmask << std::endl;
        std::cerr << "  Output rootname          : " << fn_root << std::endl;
        std::cerr << "  Convergence criterion    : " << eps << std::endl;
        if (lowpass > 0)
            std::cerr << "  Low-pass filter          : " << lowpass << std::endl;
        if (tilt_range0 > -91. || tilt_rangeF < 91.)
            std::cerr << "  Limited tilt range       : " << tilt_range0 << "  " << tilt_rangeF << std::endl;
        if (wlsart_no_start)
            std::cerr << "  -> Start wlsART reconstructions from all-zero volumes " << std::endl;
        if (reconstruct_wbp)
            std::cerr << "  -> Use weighted back-projection instead of wlsART for reconstruction" << std::endl;
        else if (reconstruct_fourier)
            std::cerr << "  -> Use fourier-interpolation instead of wlsART for reconstruction" << std::endl;
        if (do_prob_solvent)
            std::cerr << "  -> Perform probabilistic solvent flattening" << std::endl;
        std::cerr << " -----------------------------------------------------------------" << std::endl;

        // Also open and fill history file
        fh_hist.open((fn_root + ".hist").c_str(), std::ios::app);
        if (!fh_hist)
            REPORT_ERROR(3008, (std::string)"Prog_Refine3d: Cannot open file " + fn_root + ".hist");

        fh_hist << " -----------------------------------------------------------------" << std::endl;
        fh_hist << " | Read more about this program in the following publication:    |" << std::endl;
        if (fourier_mode)
            fh_hist << " |  Scheres ea. (2007)  Structure, 15, 1167-1177                 |" << std::endl;
        else
            fh_hist << " |  Scheres ea. (2007)  Nature Methods, 4, 27-29                 |" << std::endl;
        fh_hist << " |                                                               |" << std::endl;
        fh_hist << " |    *** Please cite it if this program is of use to you! ***   |" << std::endl;
        fh_hist << " -----------------------------------------------------------------" << std::endl;
        fh_hist << "--> Maximum-likelihood multi-reference 3D-refinement" << std::endl;
        fh_hist << "  Initial reference volume : " << fn_vol << std::endl;
        fh_hist << "  Experimental images:     : " << fn_sel << std::endl;
        fh_hist << "  Angular sampling rate    : " << angular << std::endl;
        if (fn_sym != "")
            fh_hist << "  Symmetry file:           : " << fn_sym << std::endl;
        fh_hist << "  Output rootname          : " << fn_root << std::endl;
        fh_hist << "  Convergence criterion    : " << eps << std::endl;
        if (lowpass > 0)
            fh_hist << "  Low-pass filter          : " << lowpass << std::endl;
        if (tilt_range0 > -91. || tilt_rangeF < 91.)
            fh_hist << "  Limited tilt range       : " << tilt_range0 << "  " << tilt_rangeF << std::endl;
        if (wlsart_no_start)
            fh_hist << "  -> Start wlsART reconstructions from all-zero volumes " << std::endl;
        if (reconstruct_wbp)
            fh_hist << "  -> Use weighted back-projection instead of wlsART for reconstruction" << std::endl;
        else if (reconstruct_fourier)
            fh_hist << "  -> Use fourier-interpolation instead of wlsART for reconstruction" << std::endl;
        if (do_prob_solvent)
            fh_hist << "  -> Perform probabilistic solvent flattening" << std::endl;
        fh_hist << " -----------------------------------------------------------------" << std::endl;

    }

}

// Fill sampling and create DFlib
void Prog_Refine3d_prm::produceSideInfo(int rank)
{
    // Precalculate sampling
    mysampling.SetSampling(angular);
    if (!mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
        REPORT_ERROR(3005, (std::string)"ml_refine3d::run Invalid symmetry" +  fn_sym);
    // by default max_tilt= +91., min_tilt= -91.
    mysampling.Compute_sampling_points(true,tilt_rangeF,tilt_range0);
    mysampling.remove_redundant_points(symmetry, sym_order);

    //Only the master creates a docfile with the library angles
    if (rank == 0)
    {
        DocFile  DFlib;
        for (int ilib = 0; ilib < mysampling.no_redundant_sampling_points_angles.size(); ilib++)
        {
            double rot=XX(mysampling.no_redundant_sampling_points_angles[ilib]);
            double tilt=YY(mysampling.no_redundant_sampling_points_angles[ilib]);
            double psi = 0.;
            DFlib.append_angles(rot,tilt,psi,"rot","tilt","psi");
        }
        FileName fn_tmp = fn_root + "_lib.doc";
        DFlib.write(fn_tmp);
    }

}

// Projection of the reference (blob) volume =================================
void Prog_Refine3d_prm::project_reference_volume(SelFile &SFlib, int rank, int size)
{

    VolumeXmipp                   vol;
    FileName                      fn_proj, fn_tmp;
    Projection                    proj;
    double                        rot, tilt, psi = 0.;
    int                           nvol, nl, nr_dir, my_rank;


    // Here all nodes fill SFlib and DFlib, but each node actually projects 
    // only a part of the projections. In this way parallellization is obtained

    // Total number of projections
    nl = Nvols * mysampling.no_redundant_sampling_points_angles.size();

    // Initialize
    SFlib.clear();
    eachvol_start.clear();
    eachvol_end.clear();

    if (verb > 0)
    {
        std::cerr << "--> projecting reference library ..." << std::endl;
        init_progress_bar(nl);
    }

    // Loop over all reference volumes
    nvol = 0;
    nr_dir = 0;
    fn_tmp = fn_root + "_lib";
    SFvol.go_beginning();
    while (!SFvol.eof())
    {
        FileName fn_img=SFvol.NextImg();
        if (fn_img=="") break;
        eachvol_start.push_back(nr_dir);
        vol.read(fn_img);
        vol().setXmippOrigin();
        
        for (int ilib = 0; ilib < mysampling.no_redundant_sampling_points_angles.size(); ilib++)
        {
            fn_proj.compose(fn_tmp, nr_dir + 1, "proj");
            rot=XX(mysampling.no_redundant_sampling_points_angles[ilib]);
            tilt=YY(mysampling.no_redundant_sampling_points_angles[ilib]);

            // Parallellization: each rank projects and writes a different direction
            my_rank = nr_dir % size;
            if (rank == my_rank)
            {
                project_Volume(vol(), proj, vol().rowNumber(), vol().colNumber(), rot, tilt, psi);
                proj.set_eulerAngles(rot, tilt, psi);
                proj.write(fn_proj);
            }

            // But all ranks gather the information in SFlib (and in eachvol_end and eachvol_start)
            SFlib.insert(fn_proj, SelLine::ACTIVE);
            nr_dir++;
            if (verb > 0 && (nr_dir % XMIPP_MAX(1, nl / 60) == 0)) progress_bar(nr_dir);
        }
        eachvol_end.push_back(nr_dir - 1);
        nvol++;
    }
    if (verb > 0)
    {
        progress_bar(nl);
        std::cerr << " -----------------------------------------------------------------" << std::endl;
    }

    // Only the master write the complete SFlib
    if (rank == 0)
    {
        fn_tmp = fn_root + "_lib.sel";
        SFlib.write(fn_tmp);
    }

}

// Make noise images for 3D SSNR calculation ===================================
void Prog_Refine3d_prm::make_noise_images(std::vector<ImageXmipp> &Iref)
{

    ImageXmipp img;
    FileName   fn_img;
    SelFile    SFt;

    SFt.clear();
    for (int i = 0; i < Iref.size(); i++)
    {
        img = Iref[i];
        img().initZeros();
        img().addNoise(0, 1, "gaussian");
        if (Iref[i].weight() > 1.) img() /= sqrt(Iref[i].weight());
        fn_img = fn_root + "_noise";
        fn_img.compose(fn_img, i, "xmp");
        img.write(fn_img);
        SFt.insert(fn_img);
    }
    fn_img = fn_root + "_noise.sel";
    SFt.write(fn_img);

}

// Reconstruction using the ML-weights ==========================================
void Prog_Refine3d_prm::reconstruction(int argc, char **argv,
                                       int iter, int volno, int noise)
{

    VolumeXmipp            new_vol;
    FileName               fn_tmp, fn_insel, fn_blob;
    SelFile                SFall, SFone;

    if (noise == 1) fn_tmp = fn_root + "_noise";
    else if (noise == 2) fn_tmp = fn_root + "_cref";
    else
    {
        fn_tmp = fn_root + "_it";
        fn_tmp.compose(fn_tmp, iter, "");
        if (iter > 1)
        {
            fn_blob = fn_root + "_it";
            fn_blob.compose(fn_blob, iter - 1, "");
        }
        else fn_blob = "";
    }

    // Setup selfile for reconstruction
    fn_insel = fn_tmp + ".sel";
    if (Nvols > 1)
    {
        fn_tmp += "_vol";
        fn_tmp.compose(fn_tmp, volno + 1, "");
        if (fn_blob != "")
        {
            fn_blob += "_vol";
            fn_blob.compose(fn_blob, volno + 1, "basis");
        }
        // Select only relevant projections to reconstruct
        SFall.read(fn_insel);
        SFall.go_beginning();
        for (int nr = eachvol_start[volno]; nr <= eachvol_end[volno]; nr++)
        {
            SFall.go_beginning();
            SFall.jump_lines(nr);
            SFone.insert(SFall.current());
        }
        fn_insel = fn_tmp + ".sel";
        SFone.write(fn_insel);
    }
    else
    {
        if (fn_blob != "") fn_blob += ".basis";
    }

    if (reconstruct_wbp)
    {
        Prog_WBP_prm           wbp_prm;
        if (verb > 0) std::cerr << "--> WBP reconstruction " << std::endl;

        // read command line (fn_sym, angular etc.)
        wbp_prm.read(argc, argv);
        wbp_prm.fn_sel = fn_insel;
        wbp_prm.do_weights = true;
        wbp_prm.do_all_matrices = true;
        wbp_prm.show();
        wbp_prm.verb = verb;
        if (volno > 0) wbp_prm.verb = 0;
        wbp_prm.fn_out = fn_tmp + ".vol";
        wbp_prm.produce_Side_info();
        wbp_prm.apply_2Dfilter_arbitrary_geometry(wbp_prm.SF, new_vol);
        new_vol.write(wbp_prm.fn_out);
    }
    else if (reconstruct_fourier)
    {
        // read command line (fn_sym, angular etc.)
        Prog_RecFourier_prm   fourier_prm;
        if (verb > 0) std::cerr << "--> Fourier-interpolation reconstruction " << std::endl;
        fourier_prm.read(argc, argv);
        fourier_prm.fn_sel = fn_insel;
        fourier_prm.fn_doc="";
        fourier_prm.do_weights = true;
        fourier_prm.fn_out = fn_tmp + ".vol";
        fourier_prm.verb = verb;
        if (volno > 0) fourier_prm.verb = 0;
        fourier_prm.show();
        fourier_prm.produce_Side_info();
        fourier_prm.run();
        new_vol=fourier_prm.Vout;
    }
    else // use wlsART
    {
        Basic_ART_Parameters   art_prm;
        Plain_ART_Parameters   dummy;
        GridVolume             new_blobs;
        GridVolume             start_blobs;
        if (verb > 0) std::cerr << "--> weighted least-squares ART reconstruction " << std::endl;

        // Read ART parameters from command line & I/O with outer loop of Refine3d
        art_prm.read(argc, argv);
        art_prm.WLS = true;
        if (fn_symmask != "") art_prm.fn_sym = "";
        if (!checkParameter(argc, argv, "-n"))
            art_prm.no_it = 10;
        if (!checkParameter(argc, argv, "-l"))
        {
            art_prm.lambda_list.resize(1);
            art_prm.lambda_list.initConstant(0.2);
        }
        if (!checkParameter(argc, argv, "-k"))
        {
            art_prm.kappa_list.resize(1);
            art_prm.kappa_list.initConstant(0.5);
        }
        art_prm.fn_sel = fn_insel;
        art_prm.fn_root = fn_tmp;
        if (noise == 1 || noise == 2)
        {
            art_prm.fn_start = "";
            art_prm.tell = false;
        }
        else if (!wlsart_no_start)
        {
            art_prm.tell = TELL_SAVE_BASIS;
            art_prm.fn_start = fn_blob;
        }
        // Reconstruct using weighted least-squares ART
        Basic_ROUT_Art(art_prm, dummy, new_vol, new_blobs);

    }

    if (verb > 0) std::cerr << " -----------------------------------------------------------------" << std::endl;

}

void Prog_Refine3d_prm::calculate_3DSSNR(Matrix1D<double> &spectral_signal, int iter)
{

    SelFile                     SFnoise;
    Matrix2D<std::complex<double> >  Faux;
    headerXmipp                 head;
    VolumeXmipp                 vol, nvol;
    FileName                    fn_tmp, fn_tmp2;
    Matrix1D<double>            alpha_signal, alpha_noise, input_signal, avg_alphaS, avg_alphaN;
    Matrix2D<double>            alpha_T, alpha_N, Msignal, Maux, Mone, mask;
    Projection                  proj;
    int                         c, dim;
    double                      ssnr, issnr, alpha, resol, volweight, sum;
    Matrix1D<int>               center(2), radial_count;

    // Read in noise reconstruction and calculate alpha's
    SFnoise.read(fn_root + "_noise.sel");
    SFnoise.ImgSize(dim, dim);

    center.initZeros();
    proj().resize(dim, dim);
    proj().setXmippOrigin();
    mask.resize(dim, dim);
    mask.setXmippOrigin();
    RaisedCosineMask(mask, dim / 2 - 2, dim / 2);

    if (verb > 0)
    {
        std::cerr << "--> calculating 3D-SSNR ..." << std::endl;
        init_progress_bar(SFnoise.ImgNo());
    }

    for (int volno = 0; volno < Nvols; volno++)
    {

        fn_tmp = fn_root + "_noise";
        fn_tmp2 = fn_root + "_cref";
        if (Nvols > 1)
        {
            fn_tmp += "_vol";
            fn_tmp.compose(fn_tmp, volno + 1, "");
            fn_tmp2 += "_vol";
            fn_tmp2.compose(fn_tmp2, volno + 1, "");
        }
        fn_tmp += ".vol";
        fn_tmp2 += ".vol";
        nvol.read(fn_tmp);
        vol.read(fn_tmp2);
        nvol().setXmippOrigin();
        vol().setXmippOrigin();
        Mone.resize(dim, dim);
        Mone.initConstant(1. / (double)(dim*dim));
        Mone.setXmippOrigin();
        SFnoise.go_beginning();

        c = 0;
        volweight = 0.;
        for (int nr = eachvol_start[volno]; nr <= eachvol_end[volno]; nr++)
        {
            SFnoise.go_beginning();
            SFnoise.jump_lines(nr);
            head.read(SFnoise.get_current_file());
            // alpha denominator
            if (c == 0) alpha_N = Mone * head.Weight();
            else alpha_N += Mone * head.Weight();
            // alpha nominator
            project_Volume(nvol(), proj, dim, dim, head.Phi(), head.Theta(), head.Psi());
            apply_cont_mask(mask, proj(), proj());
            FourierTransform(proj(), Faux);
            FFT_magnitude(Faux, Maux);
            CenterFFT(Maux, true);
            Maux.setXmippOrigin();
            Maux *= Maux;
            if (c == 0) alpha_T = Maux * head.Weight();
            else alpha_T += Maux * head.Weight();
            // input signal
            project_Volume(vol(), proj, dim, dim, head.Phi(), head.Theta(), head.Psi());
            apply_cont_mask(mask, proj(), proj());
            FourierTransform(proj(), Faux);
            FFT_magnitude(Faux, Maux);
            CenterFFT(Maux, true);
            Maux.setXmippOrigin();
            Maux *= Maux;
            if (c == 0) Msignal = Maux * head.Weight();
            else Msignal += Maux * head.Weight();
            volweight += head.Weight();
            c++;
            if (c % XMIPP_MAX(1, SFnoise.ImgNo() / 60) == 0 && verb > 0) progress_bar(c);
        }

        alpha_T.setXmippOrigin();
        alpha_N.setXmippOrigin();
        Msignal.setXmippOrigin();
        alpha_signal.initZeros();
        alpha_noise.initZeros();
        input_signal.initZeros();
        radialAverage(alpha_T, center, alpha_signal, radial_count, true);
        radialAverage(alpha_N, center, alpha_noise, radial_count, true);
        radialAverage(Msignal, center, input_signal, radial_count, true);
        input_signal /= volweight;

        // Calculate spectral_signal =input_signal/alpha!!
        // Also store averages of alphaN and alphaS for output
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(input_signal)
        {
            dVi(input_signal, i) = dVi(input_signal, i) * dVi(alpha_noise, i) / dVi(alpha_signal, i);
        }
        if (volno == 0)
        {
            spectral_signal = input_signal;
            avg_alphaN = alpha_noise;
            avg_alphaS = alpha_signal;
        }
        else
        {
            spectral_signal += input_signal;
            avg_alphaN += alpha_noise;
            avg_alphaS += alpha_signal;
        }
    }
    if (verb > 0) progress_bar(SFnoise.ImgNo());
    spectral_signal /= (double)Nvols;
    avg_alphaN /= (double)Nvols;
    avg_alphaS /= (double)Nvols;

    if (verb > 0)
    {
        fn_tmp = fn_root + "_it";
        fn_tmp.compose(fn_tmp, iter, "3dssnr");
        std::ofstream out(fn_tmp.c_str(), std::ios::out);
        out  << "#        signal    1/alpha    alpha-S    alpha-N" << std::endl;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(spectral_signal)
        {
            if (i > 0 && i < dim / 2)
            {
                out.width(5);
                out  << integerToString(i);
                out.width(10);
                out <<  floatToString(VEC_ELEM(spectral_signal, i));
                out.width(10);
                out <<  floatToString(VEC_ELEM(avg_alphaN, i) / VEC_ELEM(avg_alphaS, i));
                out.width(10);
                out <<  floatToString(VEC_ELEM(avg_alphaS, i));
                out.width(10);
                out <<  floatToString(VEC_ELEM(avg_alphaN, i));
                out << std::endl;
            }
        }
        out.close();
    }

}

void Prog_Refine3d_prm::remake_SFvol(int iter, bool rewrite, bool include_noise)
{

    FileName               fn_tmp, fn_tmp2;
    int                    volno = 0;
    VolumeXmipp            ref_vol;

    fn_tmp = fn_root + "_it";
    fn_tmp.compose(fn_tmp, iter, "");

    // Initial iteration: copy volumes to correct name for iteration
    // loop, and rewrite with this name to disc
    if (rewrite)
    {
        SFvol.go_beginning();
        while (!SFvol.eof())
        {
            FileName fn_img=SFvol.NextImg();
            if (fn_img=="") break;
            ref_vol.read(fn_img);
            ref_vol().setXmippOrigin();
            if (Nvols > 1)
            {
                fn_tmp2 = fn_tmp + "_vol";
                fn_tmp2.compose(fn_tmp2, volno + 1, "vol");
            }
            else fn_tmp2 = fn_tmp + ".vol";
            ref_vol.write(fn_tmp2);
            volno++;
        }
    }

    // Update selection file for reference volumes
    SFvol.clear();
    if (Nvols > 1)
    {
        fn_tmp += "_vol";
        volno = 0;
        while (volno < Nvols)
        {
            fn_tmp2.compose(fn_tmp, volno + 1, "vol");
            SFvol.insert(fn_tmp2);
            volno++;
        }
    }
    else
    {
        SFvol.insert(fn_tmp + ".vol");
    }
    if (include_noise)
    {
        fn_tmp = fn_root + "_noise";
        if (Nvols > 1)
        {
            fn_tmp += "_vol";
            volno = 0;
            while (volno < Nvols)
            {
                fn_tmp2.compose(fn_tmp, volno + 1, "vol");
                SFvol.insert(fn_tmp2);
                volno++;
            }
        }
        else
        {
            SFvol.insert(fn_tmp + ".vol");
        }
        // Besides noise volumes, also include cref volumes
        fn_tmp = fn_root + "_cref";
        if (Nvols > 1)
        {
            fn_tmp += "_vol";
            volno = 0;
            while (volno < Nvols)
            {
                fn_tmp2.compose(fn_tmp, volno + 1, "vol");
                SFvol.insert(fn_tmp2);
                volno++;
            }
        }
        else
        {
            SFvol.insert(fn_tmp + ".vol");
        }
    }

}

// Concatenate MLalign2D selfiles ==============================================
void Prog_Refine3d_prm::concatenate_selfiles(int iter)
{

    FileName fn_tmp, fn_class;

    // Concatenate all hard-classification selfiles
    // Only after an iteration has been performed, and thus rewrite==false...
    for (int volno = 0; volno < Nvols; volno++)
    {
        fn_class = fn_root + "_it";
        fn_class.compose(fn_class, iter, "");
        fn_class += "_class_vol";
        fn_class.compose(fn_class, volno + 1, "sel");
        system(((std::string)"rm -f " + fn_class).c_str());
        for (int nr = eachvol_start[volno]; nr <= eachvol_end[volno]; nr++)
        {
            fn_tmp = fn_root + "_ref";
            fn_tmp.compose(fn_tmp, nr + 1, "sel");
            system(((std::string)"cat " + fn_tmp + " >> " + fn_class).c_str());
            system(((std::string)"rm -f " + fn_tmp).c_str());
        }
    }

}

// Modify reference volume ======================================================
void Prog_Refine3d_prm::post_process_volumes(int argc, char **argv)
{

    Prog_segment_prm       segm_prm;
    FileName               fn_vol, fn_tmp;
    VolumeXmipp            vol, Vaux, Vsymmask, Vsolv;
    SymList                SL;
    Matrix3D<int>          mask3D;
    double                 avg, dummy, in, out;
    int                    dim;

    if ((fn_sym != "") || (lowpass > 0) ||
        (fn_solv != "") || (do_prob_solvent) || (threshold_solvent != 999))
    {

        SFvol.go_beginning();
        while (!SFvol.eof())
        {
            // Read corresponding volume from disc
            fn_vol = SFvol.NextImg();
            if (fn_vol=="") break;
            vol.read(fn_vol);
            vol().setXmippOrigin();
            dim = vol().rowNumber();
            // Store the original volume on disc
            fn_tmp = fn_vol + ".original";
            vol.write(fn_tmp);

            // Symmetrize if requested
            if (fn_sym != "")
            {
                Vaux().resize(vol());
                SL.read_sym_file(fn_sym);
                symmetrize(SL, vol, Vaux);
                // Read local symmetry mask if requested
                if (fn_symmask != "")
                {
                    Vsymmask.read(fn_symmask);
                    Vsymmask().setXmippOrigin();
                    if (Vsymmask().computeMax() > 1. || Vsymmask().computeMin() < 0.)
                        REPORT_ERROR(1, "ERROR: sym_mask should have values between 0 and 1!");
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Vsymmask())
                    {
                        in = dVkij(Vsymmask(), k, i, j);
                        out = 1. - in;
                        dVkij(vol(), k, i, j) = out * dVkij(vol(), k, i, j) + in * dVkij(Vaux(), k, i, j);
                    }
                    Vsymmask.clear();
                }
                else
                {
                    vol = Vaux;
                }
                Vaux.clear();
            }

            // Filtering the volume
            if (lowpass > 0)
            {
                FourierMask fmask;
                fmask.raised_w = 0.02;
                fmask.FilterShape = RAISED_COSINE;
                fmask.FilterBand = LOWPASS;
                fmask.w1 = lowpass;
                fmask.apply_mask_Space(vol());
            }

            // Different types of solvent flattening
            if (do_prob_solvent || (fn_solv != "") || (threshold_solvent != 999))
            {
                if (do_prob_solvent)
                {
                    // A. Probabilistic solvent flattening
                    // Write already processed volume to disc (for segment program)
                    vol.write(fn_vol);
                    segm_prm.read(argc, argv);
                    segm_prm.fn_vol = fn_vol;
                    segm_prm.fn_mask = fn_vol + ".solv";
                    segm_prm.do_prob = true;
                    std::cerr << segm_prm;
                    fh_hist << segm_prm;
                    segm_prm.produce_side_info();
                    segm_prm.segment(Vsolv);
                }
                else if (threshold_solvent != 999)
                {
                    // B. Perform flooding and separate_objects-like solvent mask
                    Vsolv = vol;
                    Vsolv().threshold("below", threshold_solvent, 0.);
                    // The following is because binarize() seems buggy
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Vsolv())
                    {
                        if (dVkij(Vsolv(), k, i, j) != 0.) dVkij(Vsolv(), k, i, j) = 1.;
                    }
                }
                else if (fn_solv != "")
                {
                    // C. Read user-provided solvent mask from disc
                    Vsolv.read(fn_solv);
                    if (Vsolv().computeMax() > 1. || Vsolv().computeMin() < 0.)
                        REPORT_ERROR(1, "ERROR: solvent mask should have values between 0 and 1!");
                }
                // Binarize Vsolv, avoiding buggy Vsolv().binarize()
                if (do_deblob_solvent || dilate_solvent > 0)
                {
                    Vsolv().threshold("below", 0.5, 0.);
                    // The following is because binarize() seems buggy
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Vsolv())
                    {
                        if (dVkij(Vsolv(), k, i, j) != 0.) dVkij(Vsolv(), k, i, j) = 1.;
                    }
                }
                if (do_deblob_solvent)
                {
                    int object_no, maxo;
                    double nr_vox, max_vox = 0.;
                    VolumeXmipp label;
                    object_no = label_volume(Vsolv(), label());
                    max_vox = 0;
                    for (int o = 0; o <= object_no; o++)
                    {
                        Vaux() = label();
                        FOR_ALL_ELEMENTS_IN_MATRIX3D(Vaux())
                        {
                            Vaux(k, i, j) = Vaux(k, i, j) == o;
                        }
                        nr_vox = Vaux().sum();
                        if (o != 0 && (nr_vox > max_vox))
                        {
                            max_vox = nr_vox;
                            maxo = o;
                            Vsolv() = Vaux();
                        }
                    }
                    label.clear();
                }
                // Dilate solvent mask (only for binary masks)
                // Dilate several times, result is summed iteratively
                if (dilate_solvent > 0)
                {
                    VolumeXmipp Vsum;
                    Vsum() = Vsolv();
                    for (int i = 0; i < dilate_solvent; i++)
                    {
                        dilate3D(Vsolv(), Vaux(), 18, 0, 1);
                        Vsum() = Vsum() + Vaux();
                        Vsolv() = Vaux();
                    }
                    Vsum() /= (double)(dilate_solvent + 1);
                    Vsolv() = Vsum();
                    Vsum.clear();
                }
                // Apply solvent mask
                Vsolv() = 1. - Vsolv();
                double solvavg = 0., sumsolv = 0.;
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(vol())
                {
                    solvavg += dVkij(Vsolv(), k, i, j) * dVkij(vol(), k, i, j);
                    sumsolv += dVkij(Vsolv(), k, i, j);
                }
                solvavg /= sumsolv;
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Vsolv())
                {
                    dVkij(vol(), k, i, j) -= dVkij(Vsolv(), k, i, j) * (dVkij(vol(), k, i, j) - solvavg);
                }
            }

            // (Re-) write post-processed volume to disc
            vol.write(fn_vol);

        }
        if (verb > 0) std::cerr << " -----------------------------------------------------------------" << std::endl;
    }

}

// Convergence check ===============================================================
bool Prog_Refine3d_prm::check_convergence(int iter)
{

    VolumeXmipp            vol, old_vol, diff_vol;
    FileName               fn_tmp;
    Mask_Params            mask_prm;
    Matrix3D<int>          mask3D;
    double                 signal, change;
    int                    dim;
    bool                   converged = true;

    if (verb > 0) std::cerr << "--> checking convergence " << std::endl;

    for (int volno = 0; volno < Nvols; volno++)
    {
        // Read corresponding volume from disc
        fn_vol = fn_root + "_it";
        fn_vol.compose(fn_vol, iter, "");
        if (Nvols > 1)
        {
            fn_vol += "_vol";
            fn_vol.compose(fn_vol, volno + 1, "vol");
        }
        else fn_vol += ".vol";
        vol.read(fn_vol);
        vol().setXmippOrigin();
        dim = vol().rowNumber();
        old_vol().initZeros(vol());
        diff_vol().initZeros(vol());

        // Only consider voxels within the spherical mask
        mask_prm.R1 = dim / 2;
        mask_prm.type = BINARY_CIRCULAR_MASK;
        mask_prm.mode = INNER_MASK;
        mask_prm.generate_3Dmask(vol());
        fn_tmp = fn_root + "_it";
        fn_tmp.compose(fn_tmp, iter - 1, "");
        if (Nvols > 1)
        {
            fn_tmp += "_vol";
            fn_tmp.compose(fn_tmp, volno + 1, "vol");
        }
        else fn_tmp += ".vol";
        old_vol.read(fn_tmp);
        diff_vol() = vol() - old_vol();
        mask_prm.apply_mask(old_vol(), old_vol());
        mask_prm.apply_mask(diff_vol(), diff_vol());
        change = diff_vol().sum2();
        signal = old_vol().sum2();
        if (change / signal > eps) converged = false;
        if (verb > 0)
        {
            if (Nvols > 1)
            {
                std::cerr << "Relative signal change volume " << volno + 1 << " = " << change / signal << std::endl;
                fh_hist << "Relative signal change volume " << volno + 1 << " = " << change / signal << std::endl;
            }
            else
            {
                std::cerr << "Relative signal change volume = " << change / signal << std::endl;
                fh_hist << "Relative signal change volume = " << change / signal << std::endl;
            }
        }
    }

    return converged;

}
