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

#include <reconstruction/mlf_align2d.h>
#define TAG_DOCFILE 12
#define TAG_DOCFILESIZE 13

int main(int argc, char **argv)
{

    int c, nn, imgno, opt_refno;
    double LL, sumw_allrefs, convv, sumcorr;
    std::vector<double> conv;
    double aux, wsum_sigma_noise, wsum_sigma_offset;
    std::vector<Matrix2D<double > > wsum_Mref, wsum_ctfMref;
    std::vector<double> sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus;
    Matrix2D<double> P_phi, Mr2, Maux;
    std::vector<std::vector<double> > Mwsum_sigma2;
    FileName fn_img;
    Matrix1D<double> oneline(0), spectral_signal, Vaux;
    // For parallelization
    int rank, size, num_img_tot;
    bool converged;

    // Init Parallel interface
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    Prog_MLFalign2D_prm prm;

    // Get input parameters
    try
    {

        // Read subsequently to avoid problems in restart procedure
        for (int proc = 0; proc < size; proc++)
        {
            if (proc == rank)
                prm.read(argc, argv);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        //Send "master" seed to slaves for same randomization
        if (rank==0)
        {
            for (int slave = 1; slave < size; slave++)
                MPI_Send(&prm.seed, 1, MPI_INT, slave, TAG_DOCFILE,
                         MPI_COMM_WORLD);
        }
        else
        {
            prm.verb = 0;
            MPI_Recv(&prm.seed, 1, MPI_INT, 0, TAG_DOCFILE, MPI_COMM_WORLD,
                     &status);
        }

        // All nodes produce general side-info
        prm.produceSideInfo(rank);
        MPI_Barrier(MPI_COMM_WORLD);

        // Some output to screen
        if (rank == 0)
        {
            prm.show();
            prm.estimateInitialNoiseSpectra();
        }
        // And produce selfile-specific side-info
        prm.produceSideInfo2(1, size, rank);
		MPI_Barrier(MPI_COMM_WORLD);

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
        Maux.resize(prm.dim, prm.dim);
        Maux.setXmippOrigin();

        // Loop over all iterations
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {
            if (prm.verb > 0)
                std::cerr << "  Multi-reference refinement:  iteration " << iter << " of " << prm.Niter << std::endl;

            // Save old reference images
            for (int refno = 0;refno < prm.n_ref; refno++)
                prm.Iold[refno]() = prm.Iref[refno]();

            // Pre-calculate pdfs
            prm.calculateInPlanePDF();

            // Integrate over all images
            prm.expectation(prm.Iref, iter,
                                 LL, sumcorr, wsum_Mref, wsum_ctfMref,
                                 Mwsum_sigma2, wsum_sigma_offset,
                                 sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus);

            // Here MPI_allreduce of all wsums,LL and sumcorr !!!
            MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            LL = aux;
            MPI_Allreduce(&sumcorr, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sumcorr = aux;
            MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_offset = aux;
            if (prm.do_kstest)
            {
                Vaux.resize(prm.sumhist);
                MPI_Allreduce(MULTIDIM_ARRAY(prm.sumhist), MULTIDIM_ARRAY(Vaux),
                              MULTIDIM_SIZE(prm.sumhist), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vaux)
                DIRECT_MULTIDIM_ELEM(prm.sumhist, n) = DIRECT_MULTIDIM_ELEM(Vaux, n);
                for (int ires = 0; ires < prm.hdim; ires++)
                {
                    Vaux.resize(prm.sumhist);
                    MPI_Allreduce(MULTIDIM_ARRAY(prm.resolhist[ires]), MULTIDIM_ARRAY(Vaux),
                                  MULTIDIM_SIZE(prm.resolhist[ires]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vaux)
                    DIRECT_MULTIDIM_ELEM(prm.resolhist[ires], n) = DIRECT_MULTIDIM_ELEM(Vaux, n);
                }
            }
            for (int refno = 0;refno < prm.n_ref; refno++)
            {
                MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[refno]), MULTIDIM_ARRAY(Maux),
                              MULTIDIM_SIZE(wsum_Mref[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                wsum_Mref[refno] = Maux;
                if (prm.do_ctf_correction)
                {
                    MPI_Allreduce(MULTIDIM_ARRAY(wsum_ctfMref[refno]), MULTIDIM_ARRAY(Maux),
                                  MULTIDIM_SIZE(wsum_ctfMref[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    wsum_ctfMref[refno] = Maux;
                }
                MPI_Allreduce(&sumw[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw[refno] = aux;
                MPI_Allreduce(&sumw2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw2[refno] = aux;
                MPI_Allreduce(&sumwsc2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumwsc2[refno] = aux;
                MPI_Allreduce(&sumwsc[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumwsc[refno] = aux;
                MPI_Allreduce(&sumw_mirror[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw_mirror[refno] = aux;
            }
            for (int ifocus = 0;ifocus < prm.nr_focus;ifocus++)
            {
                for (int ii = 0; ii <  Mwsum_sigma2[ifocus].size(); ii++)
                {
                    MPI_Allreduce(&Mwsum_sigma2[ifocus][ii], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    Mwsum_sigma2[ifocus][ii] = aux;
                }
                if (prm.do_student)
                {
                    MPI_Allreduce(&sumw_defocus[ifocus], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    sumw_defocus[ifocus] = aux;
                }
            }

            // Update model parameters
            prm.maximization(wsum_Mref, wsum_ctfMref,
                                 Mwsum_sigma2, wsum_sigma_offset,
                                 sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus,
                                 sumcorr, sumw_allrefs,
                                 spectral_signal);

            // Check convergence
            converged = prm.checkConvergence(conv);

            // Write intermediate files
            // Write intermediate files
            if (rank != 0)
            {
                // All slaves send docfile data to the master
                int s_size = MULTIDIM_SIZE(prm.docfiledata);
                MPI_Send(&s_size, 1, MPI_INT, 0, TAG_DOCFILESIZE,
                         MPI_COMM_WORLD);
                MPI_Send(MULTIDIM_ARRAY(prm.docfiledata), s_size, MPI_DOUBLE,
                         0, TAG_DOCFILE, MPI_COMM_WORLD);
            }
            else
            {
                // Master fills docfile
                // Master's own contribution
                prm.addPartialDocfileData(prm.docfiledata, prm.myFirstImg,
                                          prm.myLastImg);
                int s_size, first_img, last_img;
                int docCounter = 1;

                while (docCounter < size)
                {
                    // receive in order
                    MPI_Recv(&s_size, 1, MPI_INT, docCounter, TAG_DOCFILESIZE,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(MULTIDIM_ARRAY(prm.docfiledata), s_size,
                             MPI_DOUBLE, docCounter, TAG_DOCFILE,
                             MPI_COMM_WORLD, &status);
                    divide_equally(prm.nr_images_global, size, docCounter,
                                   first_img, last_img);
                    prm.addPartialDocfileData(prm.docfiledata, first_img,
                                              last_img);
                    docCounter++;

                }

                // Output all intermediate files
                prm.writeOutputFiles(iter, sumw_allrefs, LL, sumcorr, conv);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // Calculate new wiener filters
            prm.updateWienerFilters(spectral_signal, sumw_defocus, iter);

            if (converged)
            {
                if (prm.verb > 0)
                    std::cerr << " Optimization converged!" << std::endl;
                break;
            }
            MPI_Barrier(MPI_COMM_WORLD);

        } // end loop iterations

        if (rank == 0)
            prm.writeOutputFiles(-1, sumw_allrefs, LL, sumcorr, conv);

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




