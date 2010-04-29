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

/* INCLUDES ---------------------------------------------------------------- */
#include <mpi.h>
#include <reconstruction/mlf_tomo.h>

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv)
{

    Prog_mlf_tomo_prm prm;

    int c, iaux, nn, imgno, opt_refno;
    double LL, sumw_allrefs, sumcorr;
    double aux, wsum_sigma_offset, wsum_sigma_noise2;
    std::vector<Matrix3D<double > > wsum_Mref, Mref;
    std::vector<Matrix3D<double > > wsum_Mwedge;
    std::vector<double> sumw, wsum_sigma2, sum_nonzero_pixels;
    Matrix3D<double> Maux;
    std::vector<int> count_defocus;
    FileName fn_img, fn_tmp;
    Matrix1D<double> oneline(0);
    DocFile DFo, DFf;
    SelFile SFa;
    // For parallelization
    int rank, size, num_img_tot;

    // Init Parallel interface
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get input parameters
    try
    {
        prm.read(argc, argv);

        // Silence all slaves
        if (rank != 0) prm.verb = 0;

        prm.produce_Side_info();

        // Some output to screen
        if (rank == 0) prm.show();

        if (prm.fn_sig == "") prm.estimate_initial_sigma2();
        prm.produce_Side_info2();

    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    try
    {
        Maux.resize(prm.dim, prm.dim, prm.dim);
        Maux.setXmippOrigin();
        DFo.reserve(2*prm.SF.ImgNo() + 1);
        DFf.reserve(2*prm.SFr.ImgNo() + 4);
        SFa.reserve(prm.Niter*prm.nr_ref);
        SFa.clear();

        // Loop over all iterations
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {

            if (prm.verb > 0) std::cerr << "  Sub-tomogram refinement:  iteration " << iter << " of " << prm.Niter << std::endl;

            DFo.clear();
            if (rank == 0)
            {
                DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Zoff (6), WedNo (7) Ref (8), Pmax/sumP (9)");
            }

            // Integrate over all images
            prm.sum_over_all_images(prm.SF, wsum_Mref, wsum_Mwedge, sum_nonzero_pixels, wsum_sigma2,
                                    wsum_sigma_offset, sumw, LL, sumcorr, DFo);

            // Here MPI_allreduce of all wsums,LL and sumcorr !!!
            MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            LL = aux;
            MPI_Allreduce(&sumcorr, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sumcorr = aux;
            MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_offset = aux;
            for (int ires = 0; ires < prm.resol_max; ires++)
            {
                MPI_Allreduce(&wsum_sigma2[ires], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                wsum_sigma2[ires] = aux;
                MPI_Allreduce(&sum_nonzero_pixels[ires], &iaux, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                sum_nonzero_pixel[ires] = iaux;
            }
            for (int refno = 0; refno < prm.n_ref; refno++)
            {
                MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[2*refno]), MULTIDIM_ARRAY(Maux),
                              MULTIDIM_SIZE(wsum_Mref[2*refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                wsum_Mref[2*refno] = Maux;
                MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[2*refno+1]), MULTIDIM_ARRAY(Maux),
                              MULTIDIM_SIZE(wsum_Mref[2*refno+1]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                wsum_Mref[2*refno+1] = Maux;
                MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mwedge[refno]), MULTIDIM_ARRAY(Maux),
                              MULTIDIM_SIZE(wsum_Mwedge[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                wsum_Mwedge[refno] = Maux;
                MPI_Allreduce(&sumw[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw[refno] = aux;
            }

            // Update model parameters
            prm.update_parameters(wsum_Mref, wsum_Mwedge, sum_nonzero_pixels, wsum_sigma2,
                                  wsum_sigma_offset, sumw, sumcorr, sumw_allrefs, iter);


            // Do some post-processing and calculate real-space references
            prm.post_process_references(Mref);

            // All nodes write out temporary DFo
            fn_img.compose(prm.fn_root, rank, "tmpdoc");
            DFo.write(fn_img);
            MPI_Barrier(MPI_COMM_WORLD);

            if (rank == 0)
            {
                // Write out docfile with optimal transformation & references
                DFo.clear();
                for (int rank2 = 0; rank2 < size; rank2++)
                {
                    fn_img.compose(prm.fn_root, rank2, "tmpdoc");
                    int ln = DFo.LineNo();
                    DFo.append(fn_img);
                    DFo.locate(DFo.get_last_key());
                    DFo.next();
                    DFo.remove_current();
                    system(((std::string)"rm -f " + fn_img).c_str());
                }
                fn_tmp = prm.fn_root + "_it";
                fn_tmp.compose(fn_tmp, iter, "doc");
                DFo.write(fn_tmp);

                // Write rest of output files
                prm.write_output_files(iter, SFa, DFf, Mref, sumw_allrefs, sumw, LL, sumcorr);
            }


        } // end loop iterations

        if (rank == 0)
        {
            // Write out converged structures
            prm.write_output_files(-1, SFa, DFf, Mref, sumw_allrefs, sumw, LL, sumcorr);

            // Write out docfile with optimal transformation & references
            fn_img = prm.fn_root + ".doc";
            DFo.write(fn_img);
        }

    }

    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(0);
    }

    MPI_Finalize();
    return 0;

}



