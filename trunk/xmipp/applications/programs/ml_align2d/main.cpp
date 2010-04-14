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

        Model_MLalign2D block_model(prm.model.n_ref);

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

            for (prm.current_block = 0; prm.current_block < prm.blocks; prm.current_block++)
            {
#ifdef TIMING
                prm.timer.tic(ITER_E);
#endif
                std::cout << "------- ITER: " << prm.iter << " ------ BLOCK: " << prm.current_block + 1 << "-------------" << std::endl;
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
                    if (prm.iter > 1)
                    {
                        prm.readModel(block_model, prm.getBaseName("_block", prm.current_block + 1));
                        //block_model.print();
                        prm.model.substractModel(block_model);
                    }
                    prm.maximization(block_model);
                    prm.writeOutputFiles(block_model, OUT_BLOCK);
                    if (prm.iter > 1)
                    {
                        prm.model.addModel(block_model);
                    }
                }

#ifdef TIMING
            prm.timer.toc(ITER_M);
#endif

            }//close for blocks

            if (prm.blocks > 1 && prm.iter == 1)
            {
                for (prm.current_block = 0; prm.current_block < prm.blocks; prm.current_block++)
                {
                    prm.readModel(block_model, prm.getBaseName("_block", prm.current_block + 1));

                    if (prm.current_block == 0)
                        prm.model = block_model;
                    else
                        prm.model.addModel(block_model);
                }
            }

            if (prm.do_norm)
                prm.correctScaleAverage();
            // Check convergence
            converged = prm.checkConvergence();
            // Write output files
            prm.addDocfileHeaderComment();
            prm.addDocfileData(prm.docfiledata, prm.myFirstImg, prm.myLastImg);
            prm.writeOutputFiles(prm.model, OUT_ITER);

            if (converged)
            {
                if (prm.verb > 0) std::cerr << " Optimization converged!"
                        << std::endl;

                break;
            }
            else
            {
                // reset DFo
                prm.DFo.clear();
            }

#ifdef TIMING
            std::cout << "-------------------- ITER: " << prm.iter << " ----------------------" << std::endl;
            prm.timer.toc(ITER);
            prm.timer.printTimes(true);

#endif

        } // end loop iterations

        prm.writeOutputFiles(prm.model);
        prm.destroyThreads();

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}
