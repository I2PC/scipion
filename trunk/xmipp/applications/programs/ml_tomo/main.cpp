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

#include <reconstruction/ml_tomo.h>

int main(int argc, char **argv)
{
    int c, nn, imgno, opt_refno;
    double LL, sumw_allrefs, convv, sumcorr;
    bool converged;
    std::vector<double> conv;
    double aux, wsum_sigma_noise, wsum_sigma_offset;
    std::vector<Matrix3D<double > > wsumimgs;
    std::vector<Matrix3D<double > > wsumweds;
    std::vector<Matrix1D<double > > fsc;
    Matrix1D<double> sumw;
    Matrix3D<double> P_phi, Mr2, Maux;
    FileName fn_img, fn_tmp;
    Matrix1D<double> oneline(0);
    DocFile DFo;

    Prog_ml_tomo_prm prm;

    // Get input parameters
    try
    {
        prm.read(argc, argv);
        prm.produceSideInfo();
        prm.show();
        if (prm.fn_ref == "")
        {
            if (prm.nr_ref != 0)
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
        Maux.resize(prm.dim, prm.dim, prm.dim);
        Maux.setXmippOrigin();

        // Loop over all iterations
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {

            if (prm.verb > 0) std::cerr << "  Multi-reference refinement:  iteration " << iter << " of " << prm.Niter << std::endl;

            // Save old reference images
            for (int refno = 0;refno < prm.nr_ref; refno++) 
                prm.Iold[refno]() = prm.Iref[refno]();

            DFo.clear();
            DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Zoff (6), Ref (7), Wedge (8), Pmax/sumP (9), LL (10), bgmean (11), scale (12)");

            // Integrate over all images
            prm.expectation(prm.SF, prm.Iref, iter,
                            LL, sumcorr, DFo, wsumimgs, wsumweds,
                            wsum_sigma_noise, wsum_sigma_offset, sumw);

            // Update model parameters
            prm.maximization(wsumimgs, wsumweds,
                             wsum_sigma_noise, wsum_sigma_offset, 
                             sumw, sumcorr, sumw_allrefs, fsc, iter);

            // Check convergence
            converged = prm.checkConvergence(conv);

            prm.writeOutputFiles(iter, DFo, wsumweds, sumw_allrefs, LL, sumcorr, conv, fsc);

            if (converged)
            {
                if (prm.verb > 0) std::cerr << " Optimization converged!" << std::endl;
                break;
            }

        } // end loop iterations
        prm.writeOutputFiles(-1, DFo, wsumweds, sumw_allrefs, LL, sumcorr, conv, fsc);

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}
