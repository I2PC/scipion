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

/* INCLUDES ---------------------------------------------------------------- */
#include <reconstruction/mlf_tomo.h>

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv)
{

    int c, nn, imgno, opt_refno;
    double LL, sumw_allrefs, sumcorr;
    double aux, wsum_sigma_offset, wsum_sigma_noise2;
    vector<matrix3D<double > > wsum_Mref, Mref;
    vector<matrix3D<double > > wsum_Mwedge;
    vector<double> sumw, wsum_sigma2, sum_nonzero_pixels;
    matrix3D<double> Maux;
    vector<int> count_defocus;
    FileName fn_img, fn_tmp;
    Matrix1D<double> oneline(0);
    DocFile DFo, DFf;
    SelFile SFa;

    Prog_mlf_tomo_prm prm;

    // Get input parameters
    try
    {
        prm.read(argc, argv);
        prm.produce_Side_info();
        prm.show();
        if (prm.fn_sig == "") prm.estimate_initial_sigma2();
        prm.produce_Side_info2();

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        prm.usage();
        exit(0);
    }

    try
    {
        Maux.resize(prm.dim, prm.dim, prm.dim);
        Maux.set_Xmipp_origin();
        DFo.reserve(2*prm.SF.ImgNo() + 1);
        DFf.reserve(2*prm.SFr.ImgNo() + 4);
        SFa.reserve(prm.Niter*prm.nr_ref);
        SFa.clear();

        // Loop over all iterations
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {

            if (prm.verb > 0) cerr << "  Sub-tomogram refinement:  iteration " << iter << " of " << prm.Niter << endl;

            DFo.clear();
            DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Zoff (6), WedNo (7) Ref (8), Pmax/sumP (9)");

            // Integrate over all images
            prm.sum_over_all_images(prm.SF, wsum_Mref, wsum_Mwedge, sum_nonzero_pixels, wsum_sigma2,
                                    wsum_sigma_offset, sumw, LL, sumcorr, DFo);

            // Update model parameters
            prm.update_parameters(wsum_Mref, wsum_Mwedge, sum_nonzero_pixels, wsum_sigma2,
                                  wsum_sigma_offset, sumw, sumcorr, sumw_allrefs, iter);

            // Do some post-processing and calculate real-space references
            prm.post_process_references(Mref);

            prm.write_output_files(iter, SFa, DFf, Mref, sumw_allrefs, sumw, LL, sumcorr);

            // Write out docfile with optimal transformation & references
            fn_tmp = prm.fn_root + "_it";
            fn_tmp.compose(fn_tmp, iter, "doc");
            DFo.write(fn_tmp);

        } // end loop iterations

        // Write out converged structures
        prm.write_output_files(-1, SFa, DFf, Mref, sumw_allrefs, sumw, LL, sumcorr);

        // Write out docfile with optimal transformation & references
        fn_img = prm.fn_root + ".doc";
        DFo.write(fn_img);

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        prm.usage();
        exit(0);
    }
}



