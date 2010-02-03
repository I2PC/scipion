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

#include <reconstruction/ml_align2d.h>

int main(int argc, char **argv)
{
    int c, nn, imgno, opt_refno;
    bool converged;
    double aux;
    Matrix2D<double> Maux;
    FileName fn_img, fn_tmp;

    Prog_MLalign2D_prm prm;

    // Get input parameters
    try
    {
        prm.read(argc, argv);
        prm.produceSideInfo();
        prm.show();

        if (prm.fn_ref == "")
        {
            if (prm.model.n_ref != 0)
            {
                prm.generateInitialReferences();
            }
            else
            {
                REPORT_ERROR(1, "Please provide -ref or -nref");
            }
        }

        prm.produceSideInfo2();
        prm.createThreads();


    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

    try
    {
        Maux.resize(prm.dim, prm.dim);
        Maux.setXmippOrigin();

        // Loop over all iterations
        for (prm.iter = prm.istart; prm.iter <= prm.Niter; prm.iter++)
        {
#ifdef TIMING
            prm.timer.tic(ITER);
#endif
            if (prm.verb > 0)
                std::cerr << "  Multi-reference refinement:  iteration " << prm.iter << " of " << prm.Niter << std::endl;

            for (int refno = 0;refno < prm.model.n_ref; refno++)
                prm.Iold[refno]() = prm.model.Iref[refno]();

            Model_MLalign2D tmp_model(prm.model.n_ref);

            for (prm.current_block = 0; prm.current_block < prm.blocks; prm.current_block++)
            {
#ifdef TIMING
            prm.timer.tic(ITER_E);
#endif
            std::cout << "------- ITER: " << prm.iter << " ------ BLOCK: " << prm.current_block << "-------------" << std::endl;
            // Integrate over all images
            prm.expectation();
#ifdef TIMING
            prm.timer.toc(ITER_E);
            prm.timer.tic(ITER_M);
#endif
            if (prm.blocks == 1) //ie not IEM
            {
                prm.maximization(prm.model);
            }
            else //do IEM
            {
                Model_MLalign2D block_model(prm.model.n_ref);
                prm.maximization(block_model);
                if (prm.current_block == 0)
                    tmp_model = block_model;
                else
                    tmp_model.addModel(block_model);
            }

#ifdef TIMING
            prm.timer.toc(ITER_M);
#endif

            }//close for blocks

            if (prm.blocks > 1)
                prm.model = tmp_model;

            if (prm.do_norm)
                prm.correctScaleAverage();

            std::cout << "------- AFTER ITER: " << prm.iter << " ------" << std::endl;
            std::cerr << "sumw_allrefs: " << prm.model.sumw_allrefs << std::endl;
               std::cerr << "wsum_sigma_offset: " << prm.model.get_wsum_sigma_offset() << std::endl;
               std::cerr << "wsum_sigma_noise: " << prm.model.get_wsum_sigma_noise() << std::endl;
               std::cerr << "sigma_offset: " << prm.model.sigma_offset << std::endl;
               std::cerr << "sigma_noise: " << prm.model.sigma_noise << std::endl;

               for (int refno = 0; refno < prm.model.n_ref; refno++)
               {
                   std::cerr << "refno:       " << refno << std::endl;
                   std::cerr << "sumw:        " << prm.model.get_sumw(refno) << std::endl;
                   std::cerr << "sumw_mirror: " << prm.model.get_sumw_mirror(refno) << std::endl;
                   std::cerr << "alpha_k:        " << prm.model.alpha_k[refno] << std::endl;
                  std::cerr << "mirror_fraction: " << prm.model.mirror_fraction[refno] << std::endl;

               }

            // Check convergence
            converged = prm.checkConvergence();

            // Write output files
            FileName fn_base = prm.getBaseName("_it", prm.iter);
            prm.writeDocfile(fn_base);
            prm.writeModel(prm.model, fn_base);

            if (converged)
            {
                if (prm.verb > 0)
                    std::cerr << " Optimization converged!" << std::endl;

                break;
            }

#ifdef TIMING
            std::cout << "-------------------- ITER: " << prm.iter << " ----------------------" << std::endl;
            prm.timer.toc(ITER);
            prm.timer.printTimes(true);

#endif

        } // end loop iterations

        //Write final output files
        prm.writeDocfile(prm.fn_root);
        prm.writeModel(prm.model, prm.fn_root);
        prm.destroyThreads();

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}
