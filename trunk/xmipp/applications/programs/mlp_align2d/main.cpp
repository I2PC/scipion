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

#include <reconstruction/mlp_align2d.h>

int main(int argc, char **argv)
{

    int c, nn, imgno, opt_refno;
    bool converged;
    std::vector<double> conv;
    double aux, LL, sumw_allrefs, sumcorr, wsum_sigma_noise, wsum_sigma_offset;
    std::vector< Polar <complex<double> > > fP_wsum_imgs;
    std::vector<double> sumw, sumw_mirror;
    FileName fn_img;
    DocFile DFo;

    Prog_MLPalign2D_prm prm;

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
        // Loop over all iterations
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {

            if (prm.verb > 0) std::cerr << "  Multi-reference refinement:  iteration " << iter << " of " << prm.Niter << std::endl;

            for (int iref = 0;iref < prm.nr_ref; iref++) prm.Iold[iref]() = prm.Iref[iref]();

            DFo.clear();
	    DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");
	    
	    // Update pdf of the translations
	    prm.updatePdfTranslations();

	    // Integrate over all images
            prm.sumOverAllImages(prm.SF, prm.Iref,
				 LL, sumcorr, DFo, fP_wsum_imgs,
				 wsum_sigma_noise, wsum_sigma_offset, sumw, sumw_mirror);

            // Update model parameters
            prm.updateParameters(fP_wsum_imgs, wsum_sigma_noise, wsum_sigma_offset,
                                  sumw, sumw_mirror, sumcorr, sumw_allrefs);

            // Check convergence
            converged = prm.checkConvergence(conv);

	    // Output files
	    prm.writeOutputFiles(iter, DFo, sumw_allrefs, LL, sumcorr, conv);

            if (converged)
            {
		if (prm.verb > 0) std::cerr << " Optimization converged!" << std::endl;
		break;
	    }

        } // end loop iterations

	// Write out final files
        prm.writeOutputFiles(-1, DFo, sumw_allrefs, LL, sumcorr, conv);

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

}




