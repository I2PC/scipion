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

#include <reconstruction/ml_refine3d.h>
#define TAG_DOCFILE 12
#define TAG_DOCFILESIZE 13

int main(int argc, char **argv)
{

    int                         c, iter, volno, converged = 0, argc2 = 0;
    char                        **argv2=NULL;
    double                      LL, sumw_allrefs, convv, sumcorr, wsum_sigma_offset;
    std::vector<double>              conv;
    std::vector<MultidimArray<double> >   wsum_Mref, wsum_ctfMref;
    std::vector<std::vector<double> >     Mwsum_sigma2;
    std::vector<double>              sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus;
    MultidimArray<double>            spectral_signal, Vaux;
    MetaData                     DFo;

    // For parallelization
    int rank, size, num_img_tot;
    double                      aux;
    MultidimArray<double>       Maux;
    FileName                    fn_tmp;
    MetaData                     SFo;

    Prog_Refine3d_prm           prm;
    Prog_MLFalign2D_prm         ML2D_prm;

    // Set to true for MLF!
    prm.fourier_mode = true;

    // Init Parallel interface
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    // Get input parameters
    try
    {

        // Read command line
        prm.read(argc, argv, argc2, argv2);
        prm.produceSideInfo(rank);

        // Write starting volumes to disc with correct name for iteration loop
        if (rank == 0)
        {
            prm.show();
            prm.remake_SFvol(prm.istart - 1, true, false);
        }
        else
            prm.remake_SFvol(prm.istart - 1, false, false);
        MPI_Barrier(MPI_COMM_WORLD);

        // Read and set general MLalign2D-stuff
        ML2D_prm.read(argc2, argv2, true);
        if (rank != 0)
            ML2D_prm.verb = prm.verb = 0;
        if (!checkParameter(argc2, argv2, "-psi_step"))
            ML2D_prm.psi_step = prm.angular;
        ML2D_prm.fn_root = prm.fn_root;
        ML2D_prm.fast_mode = true;
        ML2D_prm.do_mirror = true;
        //ML2D_prm.write_docfile = true;
        //ML2D_prm.write_selfiles = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.xmd";

        // Project the reference volume
        prm.project_reference_volume(ML2D_prm.MDref, rank, size);
        MPI_Barrier(MPI_COMM_WORLD);

        // All nodes produce general side-info
        ML2D_prm.produceSideInfo();
        if (rank == 0)
            ML2D_prm.estimateInitialNoiseSpectra();
        MPI_Barrier(MPI_COMM_WORLD);

        // Select only relevant part of selfile for this rank
        mpiSelectPart(ML2D_prm.MDref, rank, size, num_img_tot);

        // All nodes read node-specific side-info into memory
        ML2D_prm.produceSideInfo2(prm.Nvols);
        ML2D_prm.Iold.clear(); // To save memory

        // Some output to screen
        if (rank == 0)
            ML2D_prm.show(true);

    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.MLF_usage();
        }
        MPI_Finalize();
        exit(1);
    }

    try
    {
        // Initialize some additional stuff
        Maux.resize(ML2D_prm.dim, ML2D_prm.dim);
        Maux.setXmippOrigin();
        for (int refno = 0; refno < ML2D_prm.n_ref; refno++)
            conv.push_back(-1.);

        // Loop over all iterations
        iter = prm.istart;
        while (!converged && iter <= prm.Niter)
        {

            if (prm.verb > 0)
            {
                std::cerr        << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << std::endl;
                prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << std::endl;
            }

            // initialize
            DFo.clear();
            conv.clear();

            // Pre-calculate pdfs
            ML2D_prm.calculateInPlanePDF();

            // Integrate over all images
            ML2D_prm.expectation(ML2D_prm.Iref, iter,
                                 LL, sumcorr,
                                 wsum_Mref, wsum_ctfMref,
                                 Mwsum_sigma2, wsum_sigma_offset,
                                 sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus);

            // Here MPI_allreduce of all weighted sums, LL, etc.
            // All nodes need the answer to calculate internally
            // sigma_noise etc. for the next iteration!
            MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            LL = aux;
            MPI_Allreduce(&sumcorr, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sumcorr = aux;
            MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_offset = aux;
            if (ML2D_prm.do_kstest)
            {
                Vaux.resize(ML2D_prm.sumhist);
                MPI_Allreduce(MULTIDIM_ARRAY(ML2D_prm.sumhist), MULTIDIM_ARRAY(Vaux),
                              MULTIDIM_SIZE(ML2D_prm.sumhist), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vaux)
                DIRECT_MULTIDIM_ELEM(ML2D_prm.sumhist, n) = DIRECT_MULTIDIM_ELEM(Vaux, n);
                for (int ires = 0; ires < ML2D_prm.hdim; ires++)
                {
                    Vaux.resize(ML2D_prm.sumhist);
                    MPI_Allreduce(MULTIDIM_ARRAY(ML2D_prm.resolhist[ires]), MULTIDIM_ARRAY(Vaux),
                                  MULTIDIM_SIZE(ML2D_prm.resolhist[ires]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vaux)
                    DIRECT_MULTIDIM_ELEM(ML2D_prm.resolhist[ires], n) = DIRECT_MULTIDIM_ELEM(Vaux, n);
                }
            }
            for (int ifocus = 0; ifocus < ML2D_prm.nr_focus; ifocus++)
            {
                for (int ii = 0; ii <  Mwsum_sigma2[ifocus].size(); ii++)
                {
                    MPI_Allreduce(&Mwsum_sigma2[ifocus][ii], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    Mwsum_sigma2[ifocus][ii] = aux;
                }
                if (ML2D_prm.do_student)
                {
                    MPI_Allreduce(&sumw_defocus[ifocus], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    sumw_defocus[ifocus] = aux;
                }
            }
            for (int refno = 0;refno < ML2D_prm.n_ref; refno++)
            {
                if (ML2D_prm.do_ctf_correction)
                {
                    MPI_Allreduce(MULTIDIM_ARRAY(wsum_ctfMref[refno]), MULTIDIM_ARRAY(Maux),
                                  MULTIDIM_SIZE(wsum_ctfMref[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    wsum_ctfMref[refno] = Maux;
                }
                MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[refno]), MULTIDIM_ARRAY(Maux),
                              MULTIDIM_SIZE(wsum_Mref[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                wsum_Mref[refno] = Maux;
                MPI_Allreduce(&sumw[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw[refno] = aux;
                MPI_Allreduce(&sumw2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw2[refno] = aux;
                MPI_Allreduce(&sumwsc[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumwsc[refno] = aux;
                MPI_Allreduce(&sumwsc2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumwsc2[refno] = aux;
                MPI_Allreduce(&sumw_mirror[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw_mirror[refno] = aux;
            }

            // Update model parameters
            ML2D_prm.maximization(wsum_Mref, wsum_ctfMref,
                                  Mwsum_sigma2, wsum_sigma_offset,
                                  sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus,
                                  sumcorr, sumw_allrefs,
                                  spectral_signal);

            // Write intermediate files
            if (rank != 0)
            {
                // All slaves send docfile to the master
                int s_size = MULTIDIM_SIZE(ML2D_prm.docfiledata);
                MPI_Send(&s_size, 1, MPI_INT, 0, TAG_DOCFILESIZE,
                         MPI_COMM_WORLD);
                MPI_Send(MULTIDIM_ARRAY(ML2D_prm.docfiledata), s_size, MPI_DOUBLE,
                         0, TAG_DOCFILE, MPI_COMM_WORLD);
            }
            else
            {
                // Master fills docfile
                // Master's own contribution
                ML2D_prm.addPartialDocfileData(ML2D_prm.docfiledata, ML2D_prm.myFirstImg, ML2D_prm.myLastImg);

                int s_size, first_img, last_img;
                int docCounter=1;

                while (docCounter < size)
                {
                    // receive in order
                    MPI_Recv(&s_size, 1, MPI_INT, docCounter, TAG_DOCFILESIZE,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(MULTIDIM_ARRAY(ML2D_prm.docfiledata), s_size,
                             MPI_DOUBLE, docCounter, TAG_DOCFILE,
                             MPI_COMM_WORLD, &status);
                    divide_equally(ML2D_prm.nr_images_global, size, docCounter,
                                   first_img, last_img);
                    ML2D_prm.addPartialDocfileData(ML2D_prm.docfiledata, first_img,
                                              last_img);
                    docCounter++;
                }

                // Output all intermediate files
                ML2D_prm.writeOutputFiles(iter, sumw_allrefs, LL, sumcorr, conv);
                prm.concatenate_selfiles(iter);
                prm.make_noise_images(ML2D_prm.Iref);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // Reconstruct the new reference volumes also in parallel
            // Assume that the number of processors is larger than the
            // number of volumes to reconstruct ...
            for (volno = 0; volno < prm.Nvols; volno++)
                for (int nn=0; nn < 3; nn++)
                {
                    // nn=0: new reference reconstruction
                    // nn=1: noise reconstruction
                    // nn=2: ctf-corrupted reconstruction
                    if (rank == (volno * 3 + nn) % size)
                        prm.reconstruction(argc2, argv2, iter, volno, nn);
                }
            MPI_Barrier(MPI_COMM_WORLD);

            // Only the master does post-processing & convergence check (i.e. sequentially)
            if (rank == 0)
            {

                // Solvent flattening and/or symmetrizing (if requested)
                prm.remake_SFvol(iter, false, true);
                prm.post_process_volumes(argc2, argv2);

                // Calculate 3D-SSNR
                prm.calculate_3DSSNR(spectral_signal, iter);

                // Check convergence
                if (prm.check_convergence(iter))
                {
                    converged = 1;
                    if (prm.verb > 0)
                        std::cerr << "--> Optimization converged!" << std::endl;
                }

            }

            // Broadcast new spectral_signal and converged to all nodes
            MPI_Bcast(&converged, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(MULTIDIM_ARRAY(spectral_signal), MULTIDIM_SIZE(spectral_signal),
                      MPI_DOUBLE, 0, MPI_COMM_WORLD);

            // Update filenames in SFvol (now without noise volumes!)
            prm.remake_SFvol(iter, false, false);
            ML2D_prm.updateWienerFilters(spectral_signal, sumw_defocus, iter);

            if (!converged && iter + 1 <= prm.Niter)
            {
                // All nodes again: project and read new references from disc
                prm.project_reference_volume(ML2D_prm.MDref, rank, size);
                MPI_Barrier(MPI_COMM_WORLD);
                FileName fn_img;
                int c = 0;
                FOR_ALL_OBJECTS_IN_METADATA(ML2D_prm.MDref)
                {
                    ML2D_prm.MDref.getValue(MDL_IMAGE, fn_img);
                    ML2D_prm.Iref[c].read(fn_img);
                    ML2D_prm.Iref[c]().setXmippOrigin();
                    c++;
                }
            }

            iter++;
        } // end loop iterations

        if (rank == 0)
            ML2D_prm.writeOutputFiles(-1, sumw_allrefs, LL, sumcorr, conv);

        if (!converged && prm.verb > 0)
            std::cerr << "--> Optimization was stopped before convergence was reached!" << std::endl;

    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.MLF_usage();
        }
        MPI_Finalize();
        exit(1);
    }

    MPI_Finalize();
    return 0;

}




