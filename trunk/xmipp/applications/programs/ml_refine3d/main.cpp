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
    int                         c, volno, converged = 0, argc2 = 0;
    char                        **argv2=NULL;
    Prog_Refine3d_prm           prm;
    Prog_MLalign2D_prm          ML2D_prm;

    // Set to false for ML3D
    prm.fourier_mode = false;

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
        ML2D_prm.fast_mode = true;
        ML2D_prm.do_mirror = true;
        ML2D_prm.save_mem2 = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.sel";
        // Project volume and read lots of stuff into memory
        prm.project_reference_volume(ML2D_prm.SFr);

        //FIXME: TAKE CARE ABOUT SHOW
        ML2D_prm.newProduceSideInfo();
        //ML2D_prm.produceSideInfo();
        //ML2D_prm.produceSideInfo2(prm.Nvols);
        //ML2D_prm.show(true);

        // Initialize some stuff
        ML2D_prm.Iold.clear(); // To save memory
        ML2D_prm.createThreads();

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

    try
    {

        // Loop over all iterations
        ML2D_prm.iter = prm.istart;
        while (!converged && ML2D_prm.iter <= prm.Niter)
        {

            if (prm.verb > 0)
            {
                std::cerr        << "--> 3D-EM volume refinement:  iteration " << ML2D_prm.iter << " of " << prm.Niter << std::endl;
                prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << ML2D_prm.iter << " of " << prm.Niter << std::endl;
            }

            //ML2D_prm.DFo.clear();
            //ML2D_prm.DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8), LL (9), bgmean (10), scale (11), w_robust (12)");

            // Integrate over all images
            ML2D_prm.expectation();

            // Update model parameters
            ML2D_prm.maximization(ML2D_prm.model, prm.eachvol_end[0]+1);

            // Write intermediate output files
            ML2D_prm.writeOutputFiles(ML2D_prm.model, OUT_ITER);
            prm.concatenate_selfiles(ML2D_prm.iter);

            // Jump out before 3D reconstruction
            // (Useful for some parallelization protocols)
            if (prm.skip_reconstruction)
                exit(1);

            // Reconstruct new volumes from the reference images
            for (volno = 0; volno < prm.Nvols; volno++)
                prm.reconstruction(argc2, argv2, ML2D_prm.iter, volno, 0);

            // Update the reference volume selection file
            // and post-process the volumes (for -FS also the noise volumes!)
            prm.remake_SFvol(ML2D_prm.iter, false, false);
            prm.post_process_volumes(argc2, argv2);
            prm.remake_SFvol(ML2D_prm.iter, false, false);

            // Check convergence
            if (prm.check_convergence(ML2D_prm.iter))
            {
                converged = 1;
                if (prm.verb > 0) std::cerr << "--> Optimization converged!" << std::endl;
            }

            // Re-project volumes
            if (!converged && ML2D_prm.iter + 1 <= prm.Niter)
            {
                prm.project_reference_volume(ML2D_prm.SFr);
                // Read new references from disc (I could just as well keep them in memory, maybe...)
                ML2D_prm.SFr.go_beginning();
                c = 0;
                while (!ML2D_prm.SFr.eof())
                {
                    ML2D_prm.model.Iref[c].read(ML2D_prm.SFr.NextImg(), false, false, false, false);
                    ML2D_prm.model.Iref[c]().setXmippOrigin();
                    c++;
                }
            }

            ML2D_prm.iter++;
        } // end loop iterations

        // Write out converged doc and logfiles
        ML2D_prm.writeOutputFiles(ML2D_prm.model);
        ML2D_prm.destroyThreads();

        if (!converged && prm.verb > 0) std::cerr << "--> Optimization was stopped before convergence was reached!" << std::endl;
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}
