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

#include <reconstruction/mlf_align2d.h>

int main(int argc, char **argv)
{

    int c, nn, imgno, opt_refno;
    double LL, sumw_allrefs, convv, sumcorr;
    bool converged;
    std::vector<double> conv;
    double aux, wsum_sigma_noise, wsum_sigma_offset;
    std::vector<Matrix2D<double > > wsum_Mref, wsum_ctfMref;
    std::vector<double> sumw, sumw_mirror;
    Matrix2D<double> P_phi, Mr2, Maux;
    std::vector<std::vector<double> > Mwsum_sigma2;
    FileName fn_img, fn_tmp;
    Matrix1D<double> oneline(0), spectral_signal;
    DocFile DFo;

    Prog_MLFalign2D_prm prm;

    // Get input parameters
    try
    {
        prm.read(argc, argv);
        prm.produceSideInfo();
        prm.show();
        prm.estimateInitialNoiseSpectra();
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

            DFo.clear();
	    DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");

            // Pre-calculate pdfs
            prm.calculateInPlanePDF();

            // Integrate over all images
            prm.sumOverAllImages(prm.SF, prm.Iref, iter,
				 LL, sumcorr, DFo, wsum_Mref, wsum_ctfMref,
				 Mwsum_sigma2, wsum_sigma_offset, 
				 sumw, sumw_mirror);

            // Update model parameters
            prm.updateParameters(wsum_Mref, wsum_ctfMref,
                                  Mwsum_sigma2,wsum_sigma_offset, 
                                  sumw, sumw_mirror, 
                                  sumcorr, sumw_allrefs,
                                  spectral_signal);

            // Check convergence
            converged = prm.checkConvergence(conv);

            prm.writeOutputFiles(iter, DFo, sumw_allrefs, LL, sumcorr, conv);

            // Calculate new wiener filters
	    prm.updateWienerFilters(spectral_signal, iter);

            if (converged)
            {
		if (prm.verb > 0) std::cerr << " Optimization converged!" << std::endl;
		break;
	    }

        } // end loop iterations
        prm.writeOutputFiles(-1, DFo, sumw_allrefs, LL, sumcorr, conv);

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

}




