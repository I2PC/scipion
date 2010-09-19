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
    ProgML2D          ML2D_prm;
    FileName                    fnt;

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
        if (!checkParameter(argc2, argv2, "-psi_step"))
            ML2D_prm.psi_step = prm.angular;
        ML2D_prm.fn_root = prm.fn_root;
        ML2D_prm.fast_mode = true;
        ML2D_prm.do_mirror = true;
        ML2D_prm.save_mem2 = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.xmd";
        // Project volume and read lots of stuff into memory
        prm.project_reference_volume(ML2D_prm.MDref);
        ML2D_prm.produceSideInfo();
        ML2D_prm.produceSideInfo2();
        prm.Nvols *= ML2D_prm.factor_nref;
        ML2D_prm.Iold.clear(); // To save memory
        ML2D_prm.createThreads();

    }
    catch (XmippError XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

    try
    {
        ModelML2D block_model(ML2D_prm.model.n_ref);

        // Loop over all iterations
        for (ML2D_prm.iter = ML2D_prm.istart; !converged && ML2D_prm.iter <= ML2D_prm.Niter; ML2D_prm.iter++)
        {
            if (prm.verb > 0)
            {
                std::cerr << "--> 3D-EM volume refinement:  iteration " << ML2D_prm.iter << " of " << prm.Niter << std::endl;
                prm.fh_hist  << "--> 3D-EM volume refinement:  iteration " << ML2D_prm.iter << " of " << prm.Niter << std::endl;
            }

            for (ML2D_prm.current_block = 0; ML2D_prm.current_block < ML2D_prm.blocks; ML2D_prm.current_block++)
            {
                // Project volumes
                if (ML2D_prm.iter > ML2D_prm.istart || ML2D_prm.current_block > 0)
                {
                    prm.project_reference_volume(ML2D_prm.MDref);
                    c = 0;
                    // Read new references from disc (I could just as well keep them in memory, maybe...)
                    FOR_ALL_OBJECTS_IN_METADATA(ML2D_prm.MDref)
                    {
                        ML2D_prm.MDref.getValue(MDL_IMAGE, fnt);
                        ML2D_prm.model.Iref[c].read(fnt);
                        ML2D_prm.model.Iref[c]().setXmippOrigin();
                        c++;
                    }
                }

                // Integrate over all images
                ML2D_prm.expectation();

                ML2D_prm.maximizationBlocks(prm.nr_projections);

                // Write out 2D reference images (to be used in reconstruction)
                ML2D_prm.writeOutputFiles(ML2D_prm.model, OUT_REFS);

                // Jump out before 3D reconstruction
                // (Useful for some parallelization protocols)
                if (prm.skip_reconstruction)
                    exit(1);

                // Reconstruct new volumes from the reference images
                for (volno = 0; volno < prm.Nvols; volno++)
                    prm.reconstruction(argc2, argv2, ML2D_prm.iter, volno, 0);

                // Update the reference volume selection file
                // and post-process the volumes
                prm.remake_SFvol(ML2D_prm.iter, false, false);
                prm.post_process_volumes(argc2, argv2);

            } // end loop blocks

            // Check convergence
            converged = (prm.check_convergence(ML2D_prm.iter)) ? 1 :0;

            // Write output ML2D files
            ML2D_prm.addPartialDocfileData(ML2D_prm.docfiledata, ML2D_prm.myFirstImg, ML2D_prm.myLastImg);
            ML2D_prm.writeOutputFiles(ML2D_prm.model, OUT_IMGS);
            prm.concatenate_selfiles(ML2D_prm.iter);

        } // end loop iterations

        if (prm.verb > 0)
        {
        	if (converged)
            	std::cerr << "--> Optimization converged!" << std::endl;
        	else
                std::cerr << "--> Optimization was stopped before convergence was reached!" << std::endl;
        }

        // Write converged output ML2D files
        ML2D_prm.writeOutputFiles(ML2D_prm.model);
        ML2D_prm.destroyThreads();

    }
    catch (XmippError XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}
