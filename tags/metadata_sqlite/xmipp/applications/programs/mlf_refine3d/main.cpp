/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <reconstruction/ml_refine3d.h>

int main(int argc, char **argv)
{

    int                         c, iter, volno, converged = 0, argc2 = 0;
    char                        **argv2=NULL;
    double                      LL, sumw_allrefs, convv, sumcorr, wsum_sigma_noise, wsum_sigma_offset;
    std::vector<double>              conv;
    std::vector<Matrix2D<double> >   wsum_Mref, wsum_ctfMref;
    std::vector<std::vector<double> >     Mwsum_sigma2;
    std::vector<double>              sumw, sumw2, sumwsc, sumwsc2, sumw_cv, sumw_mirror, sumw_defocus;
    Matrix1D<double>            spectral_signal;
    DocFile                     DFo;

    Prog_Refine3d_prm           prm;
    Prog_MLFalign2D_prm         ML2D_prm;

    // Set to true for MLF!
    prm.fourier_mode = true;

    // Get input parameters
    try
    {

        // Read command line
        prm.read(argc, argv, argc2, argv2);
        prm.produceSideInfo();

        prm.show();
        // Write starting volume(s) to disc with correct name for iteration loop
        prm.remake_SFvol(prm.istart - 1, true);

        // Read MLalign2D-stuff
        ML2D_prm.read(argc2, argv2, true);
        if (!checkParameter(argc2, argv2, "-psi_step")) ML2D_prm.psi_step = prm.angular;
        ML2D_prm.fn_root = prm.fn_root;
        ML2D_prm.do_mirror = true;
        ML2D_prm.write_docfile = true;
        ML2D_prm.write_selfiles = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.sel";
        // Project volume and read lots of stuff into memory
        prm.project_reference_volume(ML2D_prm.SFr);
        ML2D_prm.produceSideInfo();
        ML2D_prm.estimateInitialNoiseSpectra();
        ML2D_prm.produceSideInfo2(prm.Nvols);
        ML2D_prm.show(true);

        // Initialize some stuff
        for (int refno = 0; refno < ML2D_prm.n_ref; refno++) conv.push_back(-1.);
        ML2D_prm.Iold.clear(); // To save memory

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.MLF_usage();
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
                std::cerr        << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << std::endl;
                prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << std::endl;
            }

            DFo.clear();
            DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8), w_robust (9), scale (10), KSprob (11)");

            // Pre-calculate pdfs
            ML2D_prm.calculateInPlanePDF();

            // Integrate over all images
            ML2D_prm.sumOverAllImages(ML2D_prm.SF, ML2D_prm.Iref, iter,
				      LL, sumcorr, DFo, 
				      wsum_Mref, wsum_ctfMref,
				      Mwsum_sigma2, wsum_sigma_offset, 
				      sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus);

            // Update model parameters
            ML2D_prm.updateParameters(wsum_Mref, wsum_ctfMref,
				      Mwsum_sigma2, wsum_sigma_offset, 
				      sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus, 
				      sumcorr, sumw_allrefs,
				      spectral_signal, prm.eachvol_end[0]+1);

            // Write intermediate output files
            ML2D_prm.writeOutputFiles(iter, DFo, sumw_allrefs, LL, sumcorr, conv);
            prm.concatenate_selfiles(iter);

	    // Jump out before 3D reconstruction 
	    // (Useful for some parallelization protocols)
	    if (prm.skip_reconstruction)
		exit(1);

            // Write noise images to disc
            prm.make_noise_images(ML2D_prm.Iref);

            // Reconstruct new volumes from the reference images
            for (volno = 0; volno < prm.Nvols; volno++)
	    {
                prm.reconstruction(argc2, argv2, iter, volno, 0);
	        prm.reconstruction(argc2, argv2, iter, volno, 1);
	        prm.reconstruction(argc2, argv2, iter, volno, 2);
            }

            // Update the reference volume selection file
            // and post-process the volumes (for -FS also the noise volumes!)
            prm.remake_SFvol(iter, false, true);
            prm.post_process_volumes(argc2, argv2);
            prm.remake_SFvol(iter, false, false);

            // Calculate 3D-SSNR and new Wiener filters
	    prm.calculate_3DSSNR(spectral_signal, iter);
	    ML2D_prm.updateWienerFilters(spectral_signal, sumw_defocus, iter);

            // Check convergence
            if (prm.check_convergence(iter))
            {
                converged = 1;
                if (prm.verb > 0) std::cerr << "--> Optimization converged!" << std::endl;
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
                    FileName fn_img=ML2D_prm.SFr.NextImg();
                    if (fn_img=="") break;
                    ML2D_prm.Iref[c].read(fn_img, false, false, false, false);
                    ML2D_prm.Iref[c]().setXmippOrigin();
                    c++;
                }
            }

            iter++;
        } // end loop iterations

	// Write out converged doc and logfiles
	ML2D_prm.writeOutputFiles(-1, DFo, sumw_allrefs, LL, sumcorr, conv);

        if (!converged && prm.verb > 0) std::cerr << "--> Optimization was stopped before convergence was reached!" << std::endl;
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.MLF_usage();
        exit(0);
    }

}




