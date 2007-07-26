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
void Prog_Refine3d_prm::read(int &argc, char ** &argv)
{

    // This flag is set with scripts, so that for the user the
    // mlf_align2d and the ml_align2d are distinct programs
    fourier_mode = checkParameter(argc, argv, "-MLF");

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
        Prog_MLalign2D_prm ML_prm;
        ML_prm.extended_usage(true);
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
        string   comment, cline = "";
        DocFile  DFi;
        FileName fn_tmp;

        do_restart = true;
        DFi.read(getParameter(argc, argv, "-restart"));
        DFi.go_beginning();
        comment = (DFi.get_current_line()).get_text();
        if (strstr(comment.c_str(), "MLalign2D-logfile") == NULL)
        {
            cerr << "Error!! Docfile is not of MLalign2D-logfile type. " << endl;
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
            generateCommandLine(comment, argc, argv, copy);
            // Get number of volumes and names to generate SFvol
            if (fourier_mode) fn_root = getParameter(argc, argv, "-o", "mlf3d");
            else fn_root = getParameter(argc, argv, "-o", "ml3d");
            fn_vol = getParameter(argc, argv, "-vol");
            istart = textToInteger(getParameter(argc, argv, "-istart"));
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

    //Read Refine3d parameters
    fn_sel = getParameter(argc, argv, "-i");
    fn_root = getParameter(argc, argv, "-o", "MLrefine3D");
    if (!do_restart)
    {
        // Fill volume selfile
        fn_vol = getParameter(argc, argv, "-vol");
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

    angular = textToFloat(getParameter(argc, argv, "-ang", "10"));
    fn_sym = getParameter(argc, argv, "-sym", "");
    eps = textToFloat(getParameter(argc, argv, "-eps", "5e-5"));
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    Niter = textToInteger(getParameter(argc, argv, "-iter", "100"));
    istart = textToInteger(getParameter(argc, argv, "-istart", "1"));
    tilt_range0 = textToFloat(getParameter(argc, argv, "-tilt0", "0."));
    tilt_rangeF = textToFloat(getParameter(argc, argv, "-tiltF", "90."));
    fn_symmask = getParameter(argc, argv, "-sym_mask", "");
    lowpass = textToFloat(getParameter(argc, argv, "-filter", "-1"));
    wlsart_no_start = checkParameter(argc, argv, "-nostart");

    // Hidden for now
    fn_solv = getParameter(argc, argv, "-solvent", "");
    do_wbp = checkParameter(argc, argv, "-WBP");
    do_prob_solvent = checkParameter(argc, argv, "-prob_solvent");
    threshold_solvent = textToFloat(getParameter(argc, argv, "-threshold_solvent", "999"));
    do_deblob_solvent = checkParameter(argc, argv, "-deblob_solvent");
    dilate_solvent = textToInteger(getParameter(argc, argv, "-dilate_solvent", "0"));
    skip_reconstruction = checkParameter(argc, argv, "-skip_reconstruction");

    // Checks
    if (lowpass > 0.5) REPORT_ERROR(1, "Digital frequency for low-pass filter should be smaller than 0.5");

}

// Usage ===================================================================
void Prog_Refine3d_prm::usage()
{
    cerr << "Usage:  ml_refine3d [options] " << endl;
    cerr << "   -i <selfile>                : Selfile with input images \n"
    << "   -vol <volume/selfile>       : Initial reference volume \n"
    << "                               :  OR selfile with multiple reference volumes\n"
    << " [ -o <root=\"ml3d\"> ]          : Output rootname \n"
    << " [ -ang <float=10> ]           : Angular sampling (degrees) \n"
    << " [ -iter <int=100> ]           : Maximum number of iterations \n"
    << " [ -more_options ]             : Show additional parameters for 3D-refinement\n";

}

// MLF Usage =================================================================
void Prog_Refine3d_prm::MLF_usage()
{
    cerr << "Usage:  mlf_refine3d [options] " << endl;
    cerr << "   -i <selfile>                : Selfile of selfiles with input images per defocus group \n"
    << "   -ctfs <selfile>             : Selfile of CTF parameters files for each defocus group \n"
    << "   -vol <volume/selfile>       : Initial reference volume \n"
    << "                               :  OR selfile with multiple reference volumes\n"
    << " [ -o <root=\"mlf\"> ]           : Output rootname \n"
    << " [ -ang <float=10> ]           : Angular sampling (degrees) \n"
    << " [ -iter <int=100> ]           : Maximum number of iterations \n"
    << " [ -search_shift <float=3>]    : Limited translational searches (in pixels) \n"
    << " [ -not_phase_flipped ]        : Use this if the experimental images have not been phase flipped \n"
    << " [ -ctf_affected_refs ]        : Use this if the references are not CTF-deconvoluted \n"
    << " [ -low <pix=0> ]              : Exclude lowest freq. Fourier pixels from P-calculations (in pixels) \n"
    << " [ -more_options ]             : Show additional parameters for 3D-refinement\n";

}

// Extended usage =============================================================
void Prog_Refine3d_prm::extended_usage()
{
    cerr << "Additional options: " << endl;
    cerr << " [ -l <float=0.2> ]            : wlsART-relaxation parameter (lambda)  \n"
    << " [ -k <float=0.5> ]            : wlsART-relaxation parameter for residual (kappa)\n"
    << " [ -n <int=10> ]               : Number of wlsART-iterations \n"
    << " [ -nostart ]                  : Start wlsART reconstructions from all-zero volumes \n"
    << " [ -sym <symfile> ]            : Enforce symmetry \n"
    << " [ -filter <dig.freq.=-1> ]    : Low-pass filter volume every iteration \n"
    << " [ -sym_mask <maskfile> ]      : Local symmetry (only inside mask) \n"
    << " [ -tilt0 <float=0.> ]         : Lower-value for restricted tilt angle search \n"
    << " [ -tiltF <float=90.> ]        : Higher-value for restricted tilt angle search \n"
    << " [ -show_all_ML_options ]      : Show all parameters for the ML-refinement\n"
    << " [ -show_all_ART_options ]     : Show all parameters for the wlsART reconstruction \n";
    cerr << endl;
    exit(1);
}

// Show ======================================================================
void Prog_Refine3d_prm::show()
{

    if (verb > 0)
    {
        // To screen
        cerr << " -----------------------------------------------------------------" << endl;
        cerr << " | Read more about this program in the following publication:    |" << endl;
        if (fourier_mode)
            cerr << " |  Scheres ea. (2007)  in preparation                           |" << endl;
        else
            cerr << " |  Scheres ea. (2007)  Nature Methods, 4, 27-29                 |" << endl;
        cerr << " |                                                               |" << endl;
        cerr << " |    *** Please cite it if this program is of use to you! ***   |" << endl;
        cerr << " -----------------------------------------------------------------" << endl;
        cerr << "--> Maximum-likelihood multi-reference 3D-refinement" << endl;
        if (Nvols == 1)
            cerr << "  Initial reference volume : " << fn_vol << endl;
        else
        {
            cerr << "  Selfile with references  : " << fn_vol << endl;
            cerr << "    with # of volumes      : " << Nvols << endl;
        }
        cerr << "  Experimental images:     : " << fn_sel << endl;
        cerr << "  Angular sampling rate    : " << angular << endl;
        if (fn_sym != "")
            cerr << "  Symmetry file:           : " << fn_sym << endl;
        if (fn_symmask != "")
            cerr << "  Local symmetry mask      : " << fn_symmask << endl;
        cerr << "  Output rootname          : " << fn_root << endl;
        cerr << "  Convergence criterion    : " << eps << endl;
        if (lowpass > 0)
            cerr << "  Low-pass filter          : " << lowpass << endl;
        if (tilt_range0 > 0. || tilt_rangeF < 90.)
            cerr << "  Limited tilt range       : " << tilt_range0 << "  " << tilt_rangeF << endl;
        if (wlsart_no_start)
            cerr << "  -> Start wlsART reconstructions from all-zero volumes " << endl;
        if (do_wbp)
            cerr << "  -> Use weighted back-projection instead of ART for reconstruction" << endl;
        if (do_prob_solvent)
            cerr << "  -> Perform probabilistic solvent flattening" << endl;
        cerr << " -----------------------------------------------------------------" << endl;

        // Also open and fill history file
        fh_hist.open((fn_root + ".hist").c_str(), ios::app);
        if (!fh_hist)
            REPORT_ERROR(3008, (string)"Prog_Refine3d: Cannot open file " + fn_root + ".hist");

        fh_hist << " -----------------------------------------------------------------" << endl;
        fh_hist << " | Read more about this program in the following publication:    |" << endl;
        if (fourier_mode)
            fh_hist << " |  Scheres ea. (2007)  in preparation                           |" << endl;
        else
            fh_hist << " |  Scheres ea. (2007)  Nature Methods, 4, 27-29                 |" << endl;
        fh_hist << " |                                                               |" << endl;
        fh_hist << " |    *** Please cite it if this program is of use to you! ***   |" << endl;
        fh_hist << " -----------------------------------------------------------------" << endl;
        fh_hist << "--> Maximum-likelihood multi-reference 3D-refinement" << endl;
        fh_hist << "  Initial reference volume : " << fn_vol << endl;
        fh_hist << "  Experimental images:     : " << fn_sel << endl;
        fh_hist << "  Angular sampling rate    : " << angular << endl;
        if (fn_sym != "")
            fh_hist << "  Symmetry file:           : " << fn_sym << endl;
        fh_hist << "  Output rootname          : " << fn_root << endl;
        fh_hist << "  Convergence criterion    : " << eps << endl;
        if (lowpass > 0)
            fh_hist << "  Low-pass filter          : " << lowpass << endl;
        if (tilt_range0 > 0. || tilt_rangeF < 90.)
            fh_hist << "  Limited tilt range       : " << tilt_range0 << "  " << tilt_rangeF << endl;
        if (wlsart_no_start)
            fh_hist << "  -> Start wlsART reconstructions from all-zero volumes " << endl;
        if (do_wbp)
            fh_hist << "  -> Use weighted back-projection instead of wlsART for reconstruction" << endl;
        if (do_prob_solvent)
            fh_hist << "  -> Perform probabilistic solvent flattening" << endl;
        fh_hist << " -----------------------------------------------------------------" << endl;

    }

}

// Projection of the reference (blob) volume =================================
void Prog_Refine3d_prm::project_reference_volume(SelFile &SFlib, int rank)
{

    VolumeXmipp                   vol;
    SymList                       SL;
    DocFile                       DFlib;
    FileName                      fn_proj, fn_tmp;
    Projection                    proj;
    double                        rot, tilt, psi = 0.;
    int                           nvol, nl, nr_dir;

    if (fn_sym != "" && fn_symmask == "") SL.read_sym_file(fn_sym);
    make_even_distribution(DFlib, angular, SL, false);
    if (tilt_range0 > 0. || tilt_rangeF < 90.)
        limit_tilt_range(DFlib, tilt_range0, tilt_rangeF);
    // Select use-provided tilt range
    if (tilt_range0 > 0. || tilt_rangeF < 90.)
    {
        DocLine DL;
        DocFile Dt;
        DFlib.go_first_data_line();
        while (!DFlib.eof())
        {
            DL = DFlib.get_current_line();
            tilt = DFlib(1);
            if (tilt >= tilt_range0 && tilt <= tilt_rangeF) Dt.append_line(DL);
            DFlib.next_data_line();
        }
        DFlib = Dt;
        Dt.clear();
    }
    nl = Nvols * DFlib.dataLineNo();

    SFlib.clear();
    eachvol_start.clear();
    eachvol_end.clear();

    if (verb > 0 && rank == 0)
    {
        cerr << "--> projecting reference library ..." << endl;
        init_progress_bar(nl);
    }

    // Loop over all reference volumes
    nvol = 0;
    nr_dir = 0;
    fn_tmp = fn_root + "_lib";
    SFvol.go_beginning();
    while (!SFvol.eof())
    {
        eachvol_start.push_back(nr_dir);
        vol.read(SFvol.NextImg());
        vol().setXmippOrigin();
        DFlib.go_beginning();
        DFlib.adjust_to_data_line();
        while (!DFlib.eof())
        {
            fn_proj.compose(fn_tmp, nr_dir + 1, "proj");
            rot  = DFlib(0);
            tilt = DFlib(1);
            // In parallel case: only master actually projects and writes to disc
            if (rank == 0)
            {
                project_Volume(vol(), proj, vol().rowNumber(), vol().colNumber(), rot, tilt, psi);
                proj.set_eulerAngles(rot, tilt, psi);
                proj.write(fn_proj);
            }
            SFlib.insert(fn_proj, SelLine::ACTIVE);
            DFlib.next_data_line();
            nr_dir++;
            if (verb > 0  && rank == 0 && (nr_dir % MAX(1, nl / 60) == 0)) progress_bar(nr_dir);
        }
        eachvol_end.push_back(nr_dir - 1);
        nvol++;
    }
    if (verb > 0 && rank == 0)
    {
        progress_bar(nl);
        cerr << " -----------------------------------------------------------------" << endl;
    }

    if (rank == 0)
    {
        fn_tmp = fn_root + "_lib.sel";
        SFlib.write(fn_tmp);
        fn_tmp = fn_root + "_lib.doc";
        DFlib.write(fn_tmp);
    }

    // Free memory
    vol.clear();

}

// Make noise images for 3D SSNR calculation ===================================
void Prog_Refine3d_prm::make_noise_images(vector<ImageXmipp> &Iref)
{

    ImageXmipp img;
    FileName   fn_img;
    SelFile    SFt;

    SFt.clear();
    for (int i = 0; i < Iref.size(); i++)
    {
        img = Iref[i];
        img().initZeros();
        img().add_noise(0, 1, "gaussian");
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

    if (!do_wbp)
    {
        Basic_ART_Parameters   art_prm;
        Plain_ART_Parameters   dummy;
        GridVolume             new_blobs;
        GridVolume             start_blobs;
        if (verb > 0) cerr << "--> weighted least-squares ART reconstruction " << endl;

        // Read ART parameters from command line & I/O with outer loop of Refine3d
        art_prm.read(argc, argv);
        art_prm.WLS = true;
        if (fn_symmask != "") art_prm.fn_sym = "";
        if (!checkParameter(argc, argv, "-n"))
            art_prm.no_it = 10;
        if (!checkParameter(argc, argv, "-l"))
        {
            art_prm.lambda_list.resize(1);
            art_prm.lambda_list.init_constant(0.2);
        }
        if (!checkParameter(argc, argv, "-k"))
        {
            art_prm.kappa_list.resize(1);
            art_prm.kappa_list.init_constant(0.5);
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
    else
    {

        Prog_WBP_prm           wbp_prm;
        if (verb > 0) cerr << "--> WBP reconstruction " << endl;

        // read command line (fn_sym, angular etc.)
        wbp_prm.read(argc, argv);
        wbp_prm.fn_sel = fn_insel;
        wbp_prm.do_weights = true;
        wbp_prm.do_all_matrices = true;
        wbp_prm.show();
        wbp_prm.verb = verb;
        wbp_prm.fn_out = fn_tmp + ".vol";
        wbp_prm.produce_Side_info();
        wbp_prm.apply_2Dfilter_arbitrary_geometry(wbp_prm.SF, new_vol);
        new_vol.write(wbp_prm.fn_out);
    }

    if (verb > 0) cerr << " -----------------------------------------------------------------" << endl;

}

void Prog_Refine3d_prm::calculate_3DSSNR(Matrix1D<double> &spectral_signal, int iter)
{

    SelFile                     SFnoise;
    Matrix2D<complex<double> >  Faux;
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
        cerr << "--> calculating 3D-SSNR ..." << endl;
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
        Mone.init_constant(1. / (double)(dim*dim));
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
            if (c % MAX(1, SFnoise.ImgNo() / 60) == 0 && verb > 0) progress_bar(c);
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
        ofstream out(fn_tmp.c_str(), ios::out);
        out  << "#        signal    1/alpha    alpha-S    alpha-N" << endl;
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
                out << endl;
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
            ref_vol.read(SFvol.NextImg());
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
        system(((string)"rm -f " + fn_class).c_str());
        for (int nr = eachvol_start[volno]; nr <= eachvol_end[volno]; nr++)
        {
            fn_tmp = fn_root + "_ref";
            fn_tmp.compose(fn_tmp, nr + 1, "sel");
            system(((string)"cat " + fn_tmp + " >> " + fn_class).c_str());
            system(((string)"rm -f " + fn_tmp).c_str());
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
                    if (Vsymmask().compute_max() > 1. || Vsymmask().compute_min() < 0.)
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
                    cerr << segm_prm;
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
                    if (Vsolv().compute_max() > 1. || Vsolv().compute_min() < 0.)
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
        if (verb > 0) cerr << " -----------------------------------------------------------------" << endl;
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

    if (verb > 0) cerr << "--> checking convergence " << endl;

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
                cerr << "Relative signal change volume " << volno + 1 << " = " << change / signal << endl;
                fh_hist << "Relative signal change volume " << volno + 1 << " = " << change / signal << endl;
            }
            else
            {
                cerr << "Relative signal change volume = " << change / signal << endl;
                fh_hist << "Relative signal change volume = " << change / signal << endl;
            }
        }
    }

    return converged;

}
