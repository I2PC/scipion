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
    bool converged = false;
    double aux;
    MultidimArray<double> Maux;
    FileName fn_img, fn_tmp;


    Prog_MLalign2D_prm prm;

    // Get input parameters
    try
    {
        //Read arguments
        prm.read(argc, argv);
        //Generate initial references if not provided
        prm.produceSideInfo();
        //Do some initialization work
        prm.produceSideInfo2();
        //Create threads to be ready for work
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
        for (prm.iter = prm.istart; !converged && prm.iter <= prm.Niter; prm.iter++)
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
                // Integrate over all images
                prm.expectation();
#ifdef TIMING

                prm.timer.toc(ITER_E);
                prm.timer.tic(ITER_M);
#endif

                prm.maximizationBlocks();
#ifdef TIMING
                prm.timer.toc(ITER_M);
#endif

            }//close for blocks
            //std::cerr << "======End of blocks, MODEL: =========" <<std::endl;
            //prm.model.print();

            // Check convergence
            converged = prm.checkConvergence();

            // Write output files
            prm.addPartialDocfileData(prm.docfiledata, prm.myFirstImg, prm.myLastImg);
            prm.writeOutputFiles(prm.model, OUT_ITER);


#ifdef TIMING
            std::cout << "-------------------- ITER: " << prm.iter << " ----------------------" << std::endl;
            prm.timer.toc(ITER);
            prm.timer.printTimes(true);

#endif

        } // end loop iterations

        if (converged && prm.verb > 0)
            std::cerr << " Optimization converged!" << std::endl;

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
