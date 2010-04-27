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

#include <reconstruction/ml_align2d.h>
#define TAG_DOCFILE 12
#define TAG_DOCFILESIZE 13

int main(int argc, char **argv)
{
    Prog_MLalign2D_prm prm;

    int c, nn, imgno, opt_refno, iaux;
    double aux, convv;
    Matrix2D<double> Maux;
    Matrix1D<double> Vaux;
    FileName fn_img;
    // For parallelization
    int rank, size, num_img_tot;
    bool converged;

    // Init Parallel interface
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    // Get input parameters
    try
    {
        // Read subsequently to avoid problems in restart procedure
        for (int proc = 0; proc < size; proc++)
        {
            if (proc == rank) prm.read(argc, argv);
            MPI_Barrier( MPI_COMM_WORLD);
        }

        //Send "master" seed to slaves for same randomization
        if (IS_MASTER)
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

        // Create references from random subset averages
        prm.produceSideInfo(rank);
        //Syncronize all before read references (slaves will wait for master to terminate generateInitialReferences)
        MPI_Barrier(MPI_COMM_WORLD);
        // All nodes produce general side-info
        prm.produceSideInfo2(size, rank);
        prm.createThreads();
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
        Model_MLalign2D block_model(prm.model.n_ref);

        // Loop over all iterations
        for (prm.iter = prm.istart; prm.iter <= prm.Niter; prm.iter++)
        {
            if (prm.verb > 0) std::cerr
                    << "  Multi-reference refinement:  iteration " << prm.iter
                    << " of " << prm.Niter << std::endl;

            // Save old reference images
            for (int refno = 0; refno < prm.model.n_ref; refno++)
                prm.Iold[refno]() = prm.model.Iref[refno]();

            for (prm.current_block = 0; prm.current_block < prm.blocks; prm.current_block++)
            {
                // Integrate over all images
                prm.expectation();

                // Here MPI_allreduce of all wsums,LL and sumfracweight !!!
                MPI_Allreduce(&prm.LL, &aux, 1, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD);
                prm.LL = aux;
                MPI_Allreduce(&prm.sumfracweight, &aux, 1, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD);
                prm.sumfracweight = aux;
                MPI_Allreduce(&prm.wsum_sigma_noise, &aux, 1, MPI_DOUBLE,
                        MPI_SUM, MPI_COMM_WORLD);
                prm.wsum_sigma_noise = aux;
                MPI_Allreduce(&prm.wsum_sigma_offset, &aux, 1, MPI_DOUBLE,
                        MPI_SUM, MPI_COMM_WORLD);
                prm.wsum_sigma_offset = aux;
                for (int refno = 0; refno < prm.model.n_ref; refno++)
                {
                    MPI_Allreduce(MULTIDIM_ARRAY(prm.wsum_Mref[refno]),
                            MULTIDIM_ARRAY(Maux),
                            MULTIDIM_SIZE(prm.wsum_Mref[refno]), MPI_DOUBLE,
                            MPI_SUM, MPI_COMM_WORLD);
                    prm.wsum_Mref[refno] = Maux;
                    MPI_Allreduce(&prm.sumw[refno], &aux, 1, MPI_DOUBLE,
                            MPI_SUM, MPI_COMM_WORLD);
                    prm.sumw[refno] = aux;
                    MPI_Allreduce(&prm.sumwsc2[refno], &aux, 1, MPI_DOUBLE,
                            MPI_SUM, MPI_COMM_WORLD);
                    prm.sumwsc2[refno] = aux;
                    MPI_Allreduce(&prm.sumw_mirror[refno], &aux, 1, MPI_DOUBLE,
                            MPI_SUM, MPI_COMM_WORLD);
                    prm.sumw_mirror[refno] = aux;
                    MPI_Allreduce(&prm.sumw2[refno], &aux, 1, MPI_DOUBLE,
                            MPI_SUM, MPI_COMM_WORLD);
                    prm.sumw2[refno] = aux;
                    MPI_Allreduce(&prm.sumwsc[refno], &aux, 1, MPI_DOUBLE,
                            MPI_SUM, MPI_COMM_WORLD);
                    prm.sumwsc[refno] = aux;
                }

                if (prm.blocks == 1) //ie not IEM
                {
                    prm.maximization(prm.model);
                }
                else //do IEM
                {
                    if (prm.iter > 1)
                    {
                        prm.readModel(block_model, prm.getBaseName("_block",
                                prm.current_block + 1));
                        //block_model.print();
                        prm.model.substractModel(block_model);
                    }
                    prm.maximization(block_model);

                    if (rank == 0) prm.writeOutputFiles(block_model, OUT_BLOCK);

                    if (prm.iter > 1)
                    {
                        prm.model.addModel(block_model);
                    }
                }
            }//close for blocks

            if (prm.blocks > 1 && prm.iter == 1)
            {
                for (prm.current_block = 0; prm.current_block < prm.blocks; prm.current_block++)
                {
                    prm.readModel(block_model, prm.getBaseName("_block",
                            prm.current_block + 1));

                    if (prm.current_block == 0)
                        prm.model = block_model;
                    else
                        prm.model.addModel(block_model);
                }
            }

            if (prm.do_norm) prm.correctScaleAverage();
            // Check convergence
            converged = prm.checkConvergence();

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
                //prm.addDocfileHeaderComment();
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

                //Renumber docfile
                prm.DFo.renum();
                // Output all intermediate files
                prm.writeOutputFiles(prm.model, OUT_ITER);
            }
            MPI_Barrier( MPI_COMM_WORLD);

            if (converged)
            {
                if (prm.verb > 0) std::cerr << " Optimization converged!"
                        << std::endl;
                break;
            }
            MPI_Barrier(MPI_COMM_WORLD);

        } // end loop iterations

        if (rank == 0)
        {
            //Write final output files
            prm.writeOutputFiles(prm.model);
        }
        prm.destroyThreads();

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
