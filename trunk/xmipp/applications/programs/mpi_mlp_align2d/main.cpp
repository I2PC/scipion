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

#include <mpi.h>

#include <reconstruction/mlp_align2d.h>

int main(int argc, char **argv)
{

    int c, nn, imgno, opt_refno;
    bool converged;
    std::vector<double> conv;
    double aux, LL, sumw_allrefs, sumcorr, wsum_sigma_noise, wsum_sigma_offset;
    std::vector< Polar <complex<double> > > fP_wsum_imgs;
    std::vector<double> sumw, sumw_mirror;
    Matrix1D<double> Vaux, Vsum;
    FileName fn_img;
    DocFile DFo;
    // For parallelization
    int rank, size, num_img_tot;

    // Init Parallel interface
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Prog_MLPalign2D_prm prm;

    // Get input parameters
    try
    {
        // Read subsequently to avoid problems in restart procedure
        for (int proc = 0; proc < size; proc++)
        {
            if (proc == rank) prm.read(argc, argv);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (rank != 0) prm.verb = 0;

        // All nodes produce general side-info
        prm.produceSideInfo();

        // Some output to screen
        if (rank == 0) prm.show();

        // Create references from random subset averages, or read them from selfile
        if (prm.fn_ref == "")
        {
            if (prm.nr_ref != 0)
            {
                if (rank == 0)
                {
                    prm.generateInitialReferences();
                }
                else
                {
                    prm.fn_ref = prm.fn_root + "_it";
                    prm.fn_ref.compose(prm.fn_ref, 0, "sel");
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            else
            {
                REPORT_ERROR(1, "Please provide -ref or -nref");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Select only relevant part of selfile for this rank
        prm.SF.mpi_select_part(rank, size, num_img_tot);

        // And produce selfile-specific side-info
        prm.produceSideInfo2();
        MPI_Barrier(MPI_COMM_WORLD);

    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.usage();
            MPI_Finalize();
            exit(1);
        }
    }

    try
    {
        // Loop over all iterations
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {

            if (prm.verb > 0) std::cerr << "  Multi-reference refinement:  iteration " << iter << " of " << prm.Niter << std::endl;

            // Save old reference images
            //for (int iref = 0;iref < prm.nr_ref; iref++) prm.Iold[iref]() = prm.Iref[iref]();

            DFo.clear();
            if (rank == 0)
            {
                DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");
            }

            // Update pdf of the translations
            prm.updatePdfTranslations();

            // Integrate over all images
            prm.sumOverAllImages(prm.SF, prm.Iref,
                                 LL, sumcorr, DFo, fP_wsum_imgs,
                                 wsum_sigma_noise, wsum_sigma_offset, sumw, sumw_mirror);

            // Here MPI_allreduce of all wsums,LL and sumcorr !!!
            MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            LL = aux;

            MPI_Allreduce(&sumcorr, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sumcorr = aux;

            MPI_Allreduce(&wsum_sigma_noise, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_noise = aux;

            MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_offset = aux;

            for (int refno = 0; refno < prm.nr_ref; refno++)
            {
                convertPolarToSingleArray(fP_wsum_imgs[refno], Vaux);
                Vsum.initZeros(Vaux);
                MPI_Allreduce(MULTIDIM_ARRAY(Vaux), MULTIDIM_ARRAY(Vsum),
                              MULTIDIM_SIZE(Vaux), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                convertSingleArrayToPolar(Vsum,fP_wsum_imgs[refno]);

                MPI_Allreduce(&sumw[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw[refno] = aux;

                MPI_Allreduce(&sumw_mirror[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw_mirror[refno] = aux;
            }

            // Update model parameters
            prm.updateParameters(fP_wsum_imgs, wsum_sigma_noise, wsum_sigma_offset,
                                 sumw, sumw_mirror, sumcorr, sumw_allrefs);

            // Check convergence
            converged = prm.checkConvergence(conv);

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
                prm.writeOutputFiles(iter, DFo, sumw_allrefs, LL, sumcorr, conv);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            if (converged)
            {
                if (prm.verb > 0) std::cerr << " Optimization converged!" << std::endl;
                break;
            }
            MPI_Barrier(MPI_COMM_WORLD);

        } // end loop iterations

        // Write out final files
        if (rank == 0)
            prm.writeOutputFiles(-1, DFo, sumw_allrefs, LL, sumcorr, conv);

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

    MPI_Finalize();
    return 0;

}




