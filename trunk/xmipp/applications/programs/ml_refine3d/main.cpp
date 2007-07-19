/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

#include <reconstruction/ml_refine3d.h>

int main(int argc, char **argv)
{

    int                         c, iter, volno, converged = 0;
    ;
    double                      LL, sumw_allrefs, convv, sumcorr, wsum_sigma_noise, wsum_sigma_offset;
    vector<double>              conv;
    vector<Matrix2D<double> >   wsum_Mref, wsum_ctfMref, Mwsum_sigma2;
    vector<double>              sumw, sumw_cv, sumw_mirror;
    Matrix1D<double>            spectral_signal;
    DocFile                     DFo;

    Prog_Refine3d_prm           prm;
    Prog_MLalign2D_prm          ML2D_prm;

    // Get input parameters
    try
    {

        // Read command line
        prm.read(argc, argv);
        prm.show();
        // Write starting volume(s) to disc with correct name for iteration loop
        prm.remake_SFvol(prm.istart - 1, true);

        // Read MLalign2D-stuff
        ML2D_prm.read(argc, argv, true);
        if (!checkParameter(argc, argv, "-psi_step")) ML2D_prm.psi_step = prm.angular;
        ML2D_prm.fn_root = prm.fn_root;
        ML2D_prm.fast_mode = true;
        ML2D_prm.do_mirror = true;
        ML2D_prm.save_mem2 = true;
        ML2D_prm.write_docfile = true;
        ML2D_prm.write_selfiles = true;
        ML2D_prm.write_intermediate = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.sel";
        // Project volume and read lots of stuff into memory
        prm.project_reference_volume(ML2D_prm.SFr);
        ML2D_prm.produce_Side_info();
        if (ML2D_prm.fourier_mode) ML2D_prm.estimate_initial_sigma2();
        ML2D_prm.produce_Side_info2(prm.Nvols);
        ML2D_prm.show(true);

        // Initialize some stuff
        for (int refno = 0; refno < ML2D_prm.n_ref; refno++) conv.push_back(-1.);
        ML2D_prm.Iold.clear(); // To save memory

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        if (prm.fourier_mode) prm.MLF_usage();
        else prm.usage();
        exit(0);
    }

    try
    {

        // Loop over all iterations
        iter = prm.istart;
        while (!converged && iter <= prm.Niter)
        {

            if (prm.verb > 0)
            {
                cerr        << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << endl;
                prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << endl;
            }

            DFo.clear();
            if (ML2D_prm.maxCC_rather_than_ML)
                DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Corr (8)");
            else
                DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");

            // Pre-calculate pdfs
            if (!ML2D_prm.maxCC_rather_than_ML) ML2D_prm.calculate_pdf_phi();

            // Integrate over all images
            ML2D_prm.ML_sum_over_all_images(ML2D_prm.SF, ML2D_prm.Iref, iter,
                                            LL, sumcorr, DFo, wsum_Mref, wsum_ctfMref,
                                            wsum_sigma_noise, Mwsum_sigma2,
                                            wsum_sigma_offset, sumw, sumw_mirror);

            // Update model parameters
            ML2D_prm.update_parameters(wsum_Mref, wsum_ctfMref,
                                       wsum_sigma_noise, Mwsum_sigma2,
                                       wsum_sigma_offset, sumw,
                                       sumw_mirror, sumcorr, sumw_allrefs,
                                       spectral_signal);

            // Write intermediate output files
            ML2D_prm.write_output_files(iter, DFo, sumw_allrefs, LL, sumcorr, conv);
            prm.concatenate_selfiles(iter);

            // Write noise images to disc
            if (ML2D_prm.fourier_mode) prm.make_noise_images(ML2D_prm.Iref);

            // Reconstruct new volumes from the reference images
            for (volno = 0; volno < prm.Nvols; volno++)
                prm.reconstruction(argc, argv, iter, volno, 0);
            if (ML2D_prm.fourier_mode)
            {
                prm.make_noise_images(ML2D_prm.Iref);
                for (volno = 0; volno < prm.Nvols; volno++)
                    prm.reconstruction(argc, argv, iter, volno, 1);
            }

            // Update the reference volume selection file
            // and post-process the volumes (for -FS also the noise volumes!)
            prm.remake_SFvol(iter, false, ML2D_prm.fourier_mode);
            prm.post_process_volumes(argc, argv);
            prm.remake_SFvol(iter, false, false);

            // Calculate 3D-SSNR and new Wiener filters
            if (ML2D_prm.fourier_mode)
            {
                prm.calculate_3DSSNR(spectral_signal, iter);
                if (!ML2D_prm.do_divide_ctf)
		    ML2D_prm.calculate_wiener_defocus_series(spectral_signal, iter);
            }

            // Check convergence
            if (prm.check_convergence(iter))
            {
                converged = 1;
                if (prm.verb > 0) cerr << "--> Optimization converged!" << endl;
            }

            // Re-project volumes
            if (!converged && iter + 1 <= prm.Niter)
            {
                prm.project_reference_volume(ML2D_prm.SFr);
                // Read new references from disc (I could just as well keep them in memory, maybe...)
                ML2D_prm.SFr.go_beginning();
                c = 0;
                while (!ML2D_prm.SFr.eof())
                {
                    ML2D_prm.Iref[c].read(ML2D_prm.SFr.NextImg(), false, false, false, false);
                    ML2D_prm.Iref[c]().setXmippOrigin();
                    c++;
                }
            }

            iter++;
        } // end loop iterations

	// Write out converged doc and logfiles
	ML2D_prm.write_output_files(-1, DFo, sumw_allrefs, LL, sumcorr, conv);

        if (!converged && prm.verb > 0) cerr << "--> Optimization was stopped before convergence was reached!" << endl;
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        if (prm.fourier_mode) prm.MLF_usage();
        else prm.usage();
        exit(0);
    }

}




