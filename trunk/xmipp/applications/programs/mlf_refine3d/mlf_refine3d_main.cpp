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
    std::vector<MultidimArray<double> >   wsum_Mref, wsum_ctfMref;
    std::vector<std::vector<double> >     Mwsum_sigma2;
    std::vector<double>              sumw, sumw2, sumwsc, sumwsc2, sumw_cv, sumw_mirror, sumw_defocus;
    MultidimArray<double>            spectral_signal;
    MetaData                     DFo;

    ProgRefine3D           prm;
    ProgMLF2D         ML2D_prm(true); //set 3d flag to true

    // Set to true for MLF!
    prm.fourier_mode = true;

    // Get input parameters
    try
    {

        // Read command line
        prm.read(argc, argv);
        prm.produceSideInfo();
        prm.show();
        // Write starting volume(s) to disc with correct name for iteration loop
        prm.remakeSFvol(prm.istart - 1, true);

        // Read MLalign2D-stuff
        ML2D_prm.read(argc, argv);
        if (!checkParameter(argc2, argv2, "-psi_step"))
            ML2D_prm.psi_step = prm.angular;
        ML2D_prm.fn_root = prm.fn_root;
        ML2D_prm.do_mirror = true;
        //ML2D_prm.write_docfile = true;
        //ML2D_prm.write_selfiles = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.xmd";
        // Project volume and read lots of stuff into memory
        prm.projectReferenceVolume(ML2D_prm.MDref);
        //ML2D_prm.MDref.write(std::cerr);
        ML2D_prm.produceSideInfo();
        ML2D_prm.estimateInitialNoiseSpectra();
        ML2D_prm.produceSideInfo2(prm.Nvols);
        ML2D_prm.show(true);

        // Initialize some stuff
        for (int refno = 0; refno < ML2D_prm.n_ref; refno++)
            conv.push_back(-1.);
        ML2D_prm.Iold.clear(); // To save memory

    }
    catch (XmippError XE)
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

            if (prm.verbose > 0)
            {
                std::cerr        << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << std::endl;
                prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << std::endl;
            }

            DFo.clear();
            //DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8), w_robust (9), scale (10), KSprob (11)");

            // Pre-calculate pdfs
            ML2D_prm.calculateInPlanePDF();

            // Integrate over all images
            ML2D_prm.expectation(ML2D_prm.Iref, iter,
                                 LL, sumcorr, wsum_Mref, wsum_ctfMref,
                                 Mwsum_sigma2, wsum_sigma_offset,
                                 sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus);

            // Update model parameters
            ML2D_prm.maximization(wsum_Mref, wsum_ctfMref,
                                  Mwsum_sigma2,wsum_sigma_offset,
                                  sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus,
                                  sumcorr, sumw_allrefs,
                                  spectral_signal);
            // Write intermediate output files
            ML2D_prm.addPartialDocfileData(ML2D_prm.docfiledata, ML2D_prm.myFirstImg, ML2D_prm.myLastImg);
            ML2D_prm.writeOutputFiles(iter, sumw_allrefs, LL, sumcorr, conv);
            prm.concatenateSelfiles(iter);

            // Jump out before 3D reconstruction
            // (Useful for some parallelization protocols)
            if (prm.skip_reconstruction)
                exit(1);

            // Write noise images to disc
            prm.makeNoiseImages(ML2D_prm.Iref);

            // Reconstruct new volumes from the reference images
            for (volno = 0; volno < prm.Nvols; volno++)
            {
                prm.reconstruction(argc2, argv2, iter, volno, 0);
                prm.reconstruction(argc2, argv2, iter, volno, 1);
                prm.reconstruction(argc2, argv2, iter, volno, 2);
            }

            // Update the reference volume selection file
            // and post-process the volumes (for -FS also the noise volumes!)
            prm.remakeSFvol(iter, false, true);
            prm.postProcessVolumes(argc2, argv2);
            prm.remakeSFvol(iter, false, false);

            // Calculate 3D-SSNR and new Wiener filters
            prm.calculate3DSSNR(spectral_signal, iter);
            ML2D_prm.updateWienerFilters(spectral_signal, sumw_defocus, iter);

            // Check convergence
            if (prm.checkConvergence(iter))
            {
                converged = 1;
                if (prm.verbose > 0)
                    std::cerr << "--> Optimization converged!" << std::endl;
            }

            // Re-project volumes
            if (!converged && iter + 1 <= prm.Niter)
            {
                prm.projectReferenceVolume(ML2D_prm.MDref);
                // Read new references from disc (I could just as well keep them in memory, maybe...)
                FileName fn_img;
                int c = 0;
                FOR_ALL_OBJECTS_IN_METADATA(ML2D_prm.MDref)
                {
                    ML2D_prm.MDref.getValue(MDL_IMAGE, fn_img,__iter.objId);
                    ML2D_prm.Iref[c].read(fn_img);
                    ML2D_prm.Iref[c]().setXmippOrigin();
                    c++;
                }
            }

            iter++;
        } // end loop iterations

        // Write out converged doc and logfiles
        ML2D_prm.writeOutputFiles(-1, sumw_allrefs, LL, sumcorr, conv);

        if (!converged && prm.verbose > 0)
            std::cerr << "--> Optimization was stopped before convergence was reached!" << std::endl;
    }
    catch (XmippError XE)
    {
        std::cout << XE;
        prm.MLF_usage();
        exit(0);
    }

}




