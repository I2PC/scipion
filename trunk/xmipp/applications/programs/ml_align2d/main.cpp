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

#include <reconstruction/ml_align2d.h>

int main(int argc, char **argv)
{

    int c, nn, imgno, opt_refno;
    double LL, sumw_allrefs, convv, sumcorr;
    bool converged;
    std::vector<double> conv;
    double aux, wsum_sigma_noise, wsum_sigma_offset;
    std::vector<Matrix2D<double > > wsum_Mref;
    std::vector<double> sumw, sumw2, sumwsc, sumwsc2, sumw_mirror;
    Matrix2D<double> P_phi, Mr2, Maux;
    FileName fn_img, fn_tmp;
    Matrix1D<double> oneline(0);
    DocFile DFo;

    Prog_MLalign2D_prm prm;

    // Get input parameters
    try
    {
        prm.read(argc, argv);
        prm.produce_Side_info();
        prm.show();

        if (prm.fn_ref == "")
        {
            if (prm.n_ref != 0)
            {
                prm.generate_initial_references();
            }
            else
            {
                REPORT_ERROR(1, "Please provide -ref or -nref");
            }
        }

        prm.produce_Side_info2();

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
            if (prm.maxCC_rather_than_ML)
                DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Corr (8)");
            else
		DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8), w_robust (9), bgmean (10), scale (11), sigma (12), KSprob (13)");

            // Pre-calculate pdfs
            if (!prm.maxCC_rather_than_ML) prm.calculate_pdf_phi();

            // Integrate over all images
            prm.ML_sum_over_all_images(prm.SF, prm.Iref, iter,
                                       LL, sumcorr, DFo, wsum_Mref,
                                       wsum_sigma_noise, wsum_sigma_offset, 
				       sumw, sumw2, sumwsc, sumwsc2, sumw_mirror);

            // Update model parameters
            prm.update_parameters(wsum_Mref, 
                                  wsum_sigma_noise, wsum_sigma_offset, 
				  sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, 
				  sumcorr, sumw_allrefs);

            // Check convergence
            converged = prm.check_convergence(conv);

            if (prm.write_intermediate)
                prm.write_output_files(iter, DFo, sumw_allrefs, LL, sumcorr, conv);
            else prm.output_to_screen(iter, sumcorr, LL);

            if (converged)
            {
                if (prm.anneal - 1. < 1e-5)
                {
                    if (prm.verb > 0) std::cerr << " Optimization converged!" << std::endl;
                    break;
                }
                else
                {
                    prm.anneal -= prm.anneal_step;
                    if (prm.verb > 0) std::cerr << " Lowering annealing parameter to " << prm.anneal << std::endl;
                }

            }

        } // end loop iterations
        prm.write_output_files(-1, DFo, sumw_allrefs, LL, sumcorr, conv);

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

}




