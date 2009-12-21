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
            if (prm.n_ref != 0)
            {
                prm.generateInitialReferences();
            }
            else
            {
                REPORT_ERROR(1, "Please provide -ref or -nref");
            }
        }
        prm.produceSideInfo2();

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
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {

            if (prm.verb > 0) std::cerr << "  Multi-reference refinement:  iteration " << iter << " of " << prm.Niter << std::endl;

            for (int refno = 0;refno < prm.n_ref; refno++) prm.Iold[refno]() = prm.Iref[refno]();

            prm.DFo.clear();
            prm.DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8), LL (9), bgmean (10), scale (11), w_robust (12)");

            // Integrate over all images
            prm.expectation(iter);

            // Update model parameters
            prm.maximization();

            // Check convergence
            converged = prm.checkConvergence();

            // Write output files
            prm.writeOutputFiles(iter);

            if (converged)
            {
                if (prm.verb > 0) std::cerr << " Optimization converged!" << std::endl;
                break;
            }

        } // end loop iterations
        prm.writeOutputFiles(-1);

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}
