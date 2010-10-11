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
    int                         c, volno, converged = 0, argc2 = 0;
    char                        **argv2=NULL;
    double                      convv, aux;
    FileName                     fnt;

    // For parallelization
    int rank, size, num_img_tot;
    MultidimArray<double>            Maux;
    MultidimArray<double>            Vaux;

    Prog_Refine3d_prm            prm;
    ProgML2D           ML2D_prm(true);

    // Set to false for ML3D
    prm.fourier_mode = false;

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
        if (IS_MASTER)
        {
            prm.show();
            prm.remake_SFvol(prm.istart - 1, true, false);
        }
        else
            prm.remake_SFvol(prm.istart - 1, false, false);
        MPI_Barrier(MPI_COMM_WORLD);

        // Read and set general MLalign2D-stuff
        ML2D_prm.read(argc2, argv2);
        if (!IS_MASTER)
            ML2D_prm.verbose = prm.verb = 0;
        if (!checkParameter(argc2, argv2, "-psi_step"))
            ML2D_prm.psi_step = prm.angular;
        ML2D_prm.fn_root = prm.fn_root;
        ML2D_prm.fast_mode = true;
        ML2D_prm.do_mirror = true;
        ML2D_prm.save_mem2 = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.xmd";
        // Project the reference volume
        prm.project_reference_volume(ML2D_prm.MDref, rank, size);
        MPI_Barrier(MPI_COMM_WORLD);

        // All nodes produce general side-info
        ML2D_prm.produceSideInfo();
        //Syncronize all before read references (slaves will wait for master to terminate generateInitialReferences)
        MPI_Barrier(MPI_COMM_WORLD);
        ML2D_prm.produceSideInfo2();
        prm.Nvols *= ML2D_prm.factor_nref;
        ML2D_prm.Iold.clear(); // To save memory
        ML2D_prm.createThreads();
        MPI_Barrier(MPI_COMM_WORLD);

    }
    catch (XmippError XE)
    {
        if (IS_MASTER)
        {
            std::cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    try
    {
        // Initialize some additional stuff
        Maux.resize(ML2D_prm.dim, ML2D_prm.dim);
        Maux.setXmippOrigin();
        ModelML2D block_model(ML2D_prm.model.n_ref);

        // Loop over all iterations
        for (ML2D_prm.iter = ML2D_prm.istart; !converged && ML2D_prm.iter <= ML2D_prm.Niter; ML2D_prm.iter++)
        {

            if (prm.verb > 0)
            {
                std::cerr << "--> 3D-EM volume refinement:  iteration " << ML2D_prm.iter << " of " << prm.Niter << std::endl;
                prm.fh_hist  << "--> 3D-EM volume refinement:  iteration " << ML2D_prm.iter << " of " << prm.Niter << std::endl;
            }

            for (ML2D_prm.current_block = 0; ML2D_prm.current_block < ML2D_prm.blocks; ML2D_prm.current_block++)
            {
                // Project volumes
                if (ML2D_prm.iter > ML2D_prm.istart || ML2D_prm.current_block > 0)
                {
                    prm.project_reference_volume(ML2D_prm.MDref, rank, size);
                    c = 0;
                    // Read new references from disc (I could just as well keep them in memory, maybe...)
                    FOR_ALL_OBJECTS_IN_METADATA(ML2D_prm.MDref)
                    {
                        ML2D_prm.MDref.getValue(MDL_IMAGE, fnt);
                        ML2D_prm.model.Iref[c].read(fnt);
                        ML2D_prm.model.Iref[c]().setXmippOrigin();
                        c++;
                    }
                }

                // Integrate over all images
                ML2D_prm.expectation();

                // Here MPI_allreduce of all weighted sums, LL, etc.
                // All nodes need the answer to calculate internally
                // sigma_noise etc. for the next iteration!
                MPI_Allreduce(&ML2D_prm.LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                ML2D_prm.LL = aux;
                MPI_Allreduce(&ML2D_prm.sumfracweight, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                ML2D_prm.sumfracweight = aux;
                MPI_Allreduce(&ML2D_prm.wsum_sigma_noise, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                ML2D_prm.wsum_sigma_noise = aux;
                MPI_Allreduce(&ML2D_prm.wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                ML2D_prm.wsum_sigma_offset = aux;
                for (int refno = 0;refno < ML2D_prm.model.n_ref; refno++)
                {
                    MPI_Allreduce(MULTIDIM_ARRAY(ML2D_prm.wsum_Mref[refno]), MULTIDIM_ARRAY(Maux),
                                  MULTIDIM_SIZE(ML2D_prm.wsum_Mref[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    ML2D_prm.wsum_Mref[refno] = Maux;
                    MPI_Allreduce(&ML2D_prm.sumw[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    ML2D_prm.sumw[refno] = aux;
                    MPI_Allreduce(&ML2D_prm.sumw_mirror[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    ML2D_prm.sumw_mirror[refno] = aux;
                    MPI_Allreduce(&ML2D_prm.sumw2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    ML2D_prm.sumw2[refno] = aux;
                    MPI_Allreduce(&ML2D_prm.sumwsc[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    ML2D_prm.sumwsc[refno] = aux;
                    MPI_Allreduce(&ML2D_prm.sumwsc2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    ML2D_prm.sumwsc2[refno] = aux;
                }

                ML2D_prm.maximizationBlocks(prm.nr_projections);

                // Write out 2D reference images (to be used in reconstruction)
                if (IS_MASTER)
                    ML2D_prm.writeOutputFiles(ML2D_prm.model, OUT_REFS);
                MPI_Barrier(MPI_COMM_WORLD);

                // Jump out before 3D reconstruction
                // (Useful for some parallelization protocols)
                if (prm.skip_reconstruction)
                    exit(1);

                // Reconstruct the new reference volumes also in parallel
                // Assume that the number of processors is larger than the
                // number of volumes to reconstruct ...
                // Reconstruct new volumes from the reference images
                for (volno = 0; volno < prm.Nvols; volno++)
                    if (rank == volno % size)
                        prm.reconstruction(argc2, argv2, ML2D_prm.iter, volno, 0);

                if (rank < prm.Nvols)
                    // new reference reconstruction
                    prm.reconstruction(argc2, argv2, ML2D_prm.iter, rank, 0);
                MPI_Barrier(MPI_COMM_WORLD);

                // Update filenames in SFvol (now without noise volumes!)
                prm.remake_SFvol(ML2D_prm.iter, false, false);

                // Only the master does post-processing & convergence check (i.e. sequentially)
                if (IS_MASTER)
                {
                    // Solvent flattening and/or symmetrize (if requested)
                    prm.post_process_volumes(argc2, argv2);
                }


            } //end loop blocks

            // Check convergence
            converged = (prm.check_convergence(ML2D_prm.iter)) ? 1 :0;

            // Write intermediate files
            if (!IS_MASTER)
            {
                // All slaves send docfile data to the master
                int s_size = MULTIDIM_SIZE(ML2D_prm.docfiledata);
                MPI_Send(&s_size, 1, MPI_INT, 0, TAG_DOCFILESIZE,
                         MPI_COMM_WORLD);
                MPI_Send(MULTIDIM_ARRAY(ML2D_prm.docfiledata), s_size, MPI_DOUBLE,
                         0, TAG_DOCFILE, MPI_COMM_WORLD);
            }
            else
            {
                // Master fills docfile
                //ML2D_prm.addDocfileHeaderComment();
                // Master's own contribution
                ML2D_prm.addPartialDocfileData(ML2D_prm.docfiledata, ML2D_prm.myFirstImg,
                                               ML2D_prm.myLastImg);
                int s_size, first_img, last_img;
                int docCounter = 1;

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
                ML2D_prm.writeOutputFiles(ML2D_prm.model, OUT_IMGS);
                prm.concatenate_selfiles(ML2D_prm.iter);
            }
            MPI_Barrier(MPI_COMM_WORLD);

        } // end loop iterations

        if (IS_MASTER)
        {
            if (converged && prm.verb > 0)
                std::cerr << "--> Optimization converged!" << std::endl;
            else if (prm.verb > 0)
                std::cerr << "--> Optimization was stopped before convergence was reached!" << std::endl;

            ML2D_prm.writeOutputFiles(ML2D_prm.model);
        }

        ML2D_prm.destroyThreads();

    }
    catch (XmippError XE)
    {
        if (IS_MASTER)
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
