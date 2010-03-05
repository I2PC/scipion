/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao Castro (jrbcast@cnb.csic.es)
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

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <mpi.h>

#include <reconstruction/reconstruct_art.h>
#include <reconstruction/art_crystal.h>

/* ------------------------------------------------------------------------- */
/* Time managing stuff                                                       */
/* ------------------------------------------------------------------------- */

typedef struct
{
    double user;  /* User time. */
    double sys;   /* System time. */
    double cpu;   /* CPU time = User + System. */
    double wall;  /* Wall time. */
}
USWtime_t;

// Gets User and System times for use with MPI
void uswtime(USWtime_t *tm)
{
    struct rusage buffer;

    tm->wall = MPI_Wtime();
    getrusage(RUSAGE_SELF, &buffer);
    tm->user = (double) buffer.ru_utime.tv_sec + 1.0e-6 * buffer.ru_utime.tv_usec;
    tm->sys  = (double) buffer.ru_stime.tv_sec + 1.0e-6 * buffer.ru_stime.tv_usec;
}

/* ------------------------------------------------------------------------- */
/* Main code                                                                 */
/* ------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
    Basic_ART_Parameters   art_prm;
    Plain_ART_Parameters   eprm;
    Crystal_ART_Parameters crystal_art_prm;
    VolumeXmipp            vol_voxels, vol_voxels_aux; // Volume to reconstruct
    GridVolume             vol_basis;
    int                    crystal_mode;
    MPI_Status    status;             // Stores MPI directives status
    int      num_img_tot;  // The total amount of images
    int     num_img_node;  // Stores the total amount of images each node deals with
    int      remaining;          // Number of projections still to compute
    int      Npart;              // Number of projection to process
    int      myFirst, myLast; // Limits projections to process
    double    comms_t, aux_comm_t; // Comunications time
    double    it_t;   // iteration time
    double     cav_t;  // time for CAV weights calculation
    double               cavk_it_t;          // BiCAV weights calculation time (1 iter.)
    double     cavk_total_t;       // Sum( cavk_it_t )
    USWtime_t    recons_t;  // Reconstruction time
    double    total_t;  // Program execution time
    double     comms_t_it, aux_t; // communicattions time in one iteration
    GridVolumeT<int> GVNeq_aux;    // This is a buffer for the communication
    Matrix1D<int> Ordered_aux;

    int      rank, size;       // MPI number of proccess and number of proccesses


    // Init Parallel interface
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) total_t = MPI_Wtime();   // to perform overall execution time

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Read Art Parameters
    try
    {
        art_prm.read(argc, argv);
        // Crystal
        crystal_mode = checkParameter(argc, argv, "-crystal");
        if (crystal_mode) crystal_art_prm.read(argc, argv, art_prm);
    }
    catch (Xmipp_error &XE)
    {
        std::cout << XE;
        art_prm.usage();
        crystal_art_prm.usage_more();
        exit(1);
    }

    /*************************** PARAMETERS INITIALIZATION ***************************/
    int num_iter = art_prm.no_it;
    cavk_total_t = 0.0;
    if (rank == 0)  // Master code
    {
        comms_t = 0.0; // Initializes time

        art_prm.produce_Side_Info(vol_basis, FULL, rank);
        eprm.produce_Side_Info(art_prm, vol_basis);

        Basic_ART_Init_history(art_prm, eprm, vol_basis);

        // ordered list must be the same in all nodes
        aux_comm_t = MPI_Wtime();
        MPI_Bcast(MULTIDIM_ARRAY(art_prm.ordered_list), MULTIDIM_SIZE(art_prm.ordered_list), MPI_INT, 0, MPI_COMM_WORLD);
        comms_t += MPI_Wtime() - aux_comm_t;

        std::cout << "Filename root: " << art_prm.fn_root << std::endl;
    }
    else
    {
        // each proccess can handle its own history file
        // so we add the id number to the root filename
        FileName aux = art_prm.fn_root;

        art_prm.fn_root = art_prm.fn_root + integerToString(rank);
        art_prm.produce_Side_Info(vol_basis, FULL, rank);

        // Restore original filename.
        art_prm.fn_root = aux;
        eprm.produce_Side_Info(art_prm, vol_basis);

        // ordered list must be the same in all nodes
        MPI_Bcast(MULTIDIM_ARRAY(art_prm.ordered_list), MULTIDIM_SIZE(art_prm.ordered_list), MPI_INT, 0, MPI_COMM_WORLD);
    }

    // How many images processes each node. It is done in such a way that every node receives the
    // same amount of them

    num_img_tot = art_prm.numIMG;

    Npart = (int)((float)num_img_tot / (float)size);

    num_img_node = Npart;

    remaining = num_img_tot % size;

    // each process will only see the images it is iterested in.
    if (rank < remaining)
    {
        num_img_node++;
        myFirst = rank * Npart + rank;
    }
    else
    {
        myFirst = Npart * rank + remaining;
    }

    myLast = myFirst + num_img_node - 1;

    STARTINGX(art_prm.ordered_list) = -myFirst;

    art_prm.no_it = 1;

    GridVolume   vol_basis_aux = vol_basis;

    GridVolume   vol_aux2 = vol_basis;

    // Print some data
    if (rank == 0)
    {
        if (art_prm.parallel_mode == Basic_ART_Parameters::pSART)
        {
            if (art_prm.block_size < size) art_prm.block_size = size;  // Each processor will have at least one projection
            if (art_prm.block_size > num_img_tot) art_prm.block_size = num_img_tot;  // block_size is as much equal to num_img_tot
            std::cout << "pSART " << "TB = " << art_prm.block_size << std::endl;
        }
        else if (art_prm.parallel_mode == Basic_ART_Parameters::pCAV)
            std::cout << "pCAV" << std::endl;
        else if (art_prm.parallel_mode == Basic_ART_Parameters::pBiCAV)
        {
            if (art_prm.block_size < size) art_prm.block_size = size;  // Each processor will have at least one projection
            if (art_prm.block_size > num_img_tot) art_prm.block_size = num_img_tot;  // block_size is as much equal to num_img_tot
            std::cout << "pBiCAV " << "TB = " << art_prm.block_size << std::endl;
        }
        else if (art_prm.parallel_mode == Basic_ART_Parameters::pSIRT)
            std::cout << "pSIRT" << std::endl;
        else if (art_prm.parallel_mode == Basic_ART_Parameters::pfSIRT)
            std::cout << "pfSIRT" << std::endl;
        else std::cout << "pAVSP" << std::endl;

        std::cout << "Number of processors: " << size << std::endl;
        std::cout << "Lambda: " << art_prm.lambda_list(0) << std::endl;
        std::cout << "Iterations: " << num_iter << std::endl;
        std::cout << "Projections: " << num_img_tot << std::endl;
        std::cout << "________________________________________________________________________________\n\n" << std::endl;
    }

    art_prm.block_size /= size; // Now this variable stores how many projs. from each block belong to each node.

    /*************************** CAV weights precalculation *************************/
    if (art_prm.parallel_mode == Basic_ART_Parameters::pCAV)
    {

        // Creates and initializes special variables needed to CAV weights computation.

        /*
        EXTRA CALCULATIONS FOR CAV WEIGHTS: Each node computes its own part related to its images
        and after that they send each other their results and sum up them. This part of code has been taken and
        modified from Basic_art.cc produce_side_info().
        */

        cav_t = MPI_Wtime();
        art_prm.compute_CAV_weights(vol_basis, num_img_node);
        GVNeq_aux = *(art_prm.GVNeq);
        for (int n = 0 ; n < (art_prm.GVNeq)->VolumesNo(); n++)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            aux_comm_t = MPI_Wtime();
            MPI_Allreduce(MULTIDIM_ARRAY((*(art_prm.GVNeq))(n)()), MULTIDIM_ARRAY(GVNeq_aux(n)()), MULTIDIM_SIZE((*(art_prm.GVNeq))(n)()), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            comms_t += MPI_Wtime() - aux_comm_t;
            (*(art_prm.GVNeq))(n)() = GVNeq_aux(n)();
        }
        if (rank == 0) std::cout << "Elapsed time for pCAV weights computation: " << MPI_Wtime() - cav_t << std::endl;
    }

    /*************************** PARALLEL ART ALGORITHM ***************************/

    for (int i = 0; i < num_iter; i++)
    {
        comms_t_it = 0.0;
        cavk_it_t  = 0.0;

        it_t = MPI_Wtime();

        // A bit tricky. Allows us to use the sequential Basic_ART_iterations as
        // in parallel it runs internally only one iteration.

        art_prm.lambda_list(0) = art_prm.lambda(i);

        if (art_prm.parallel_mode == Basic_ART_Parameters::pSART)
        {

            int numsteps = Npart / art_prm.block_size;

            // could be necessary another step for remaining projections
            if ((Npart % art_prm.block_size) != 0)
                numsteps++;

            int processed = 0;

            art_prm.numIMG = art_prm.block_size;

            STARTINGX(art_prm.ordered_list) = -myFirst;

            for (int ns = 0; ns < numsteps ; ns++)
            {

                if (ns == numsteps - 1)
                {
                    art_prm.numIMG = num_img_node - processed;
                }

                for (int j = 0 ; j < vol_basis.VolumesNo() ; j++)
                    vol_aux2(j)() = vol_basis(j)();

                Basic_ART_iterations(art_prm, eprm, vol_basis, rank);

                STARTINGX(art_prm.ordered_list) -= art_prm.numIMG;

                processed += art_prm.numIMG;

                int blocksize = art_prm.numIMG * size;

                // All processors send their result and get the other's so all of them
                // have the same volume for the next step.
                for (int j = 0 ; j < vol_basis.VolumesNo() ; j++)
                {
                    vol_basis(j)() = vol_basis(j)() - vol_aux2(j)(); // Adapt result to parallel ennvironment from sequential routine
                    MPI_Barrier(MPI_COMM_WORLD);
                    aux_comm_t = MPI_Wtime();
                    MPI_Allreduce(MULTIDIM_ARRAY(vol_basis(j)()), MULTIDIM_ARRAY(vol_basis_aux(j)()), MULTIDIM_SIZE(vol_basis(j)()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    aux_t = MPI_Wtime() - aux_comm_t;
                    comms_t += aux_t;
                    comms_t_it += aux_t;
                    vol_basis(j)() = vol_aux2(j)() + (vol_basis_aux(j)() / (double) blocksize);
                }
            }
        }
        else if (art_prm.parallel_mode == Basic_ART_Parameters::pBiCAV)
        {
            // Creates and initializes special variables needed to BICAV weights computation.
            GridVolumeT<int> GVNeq_aux; // This is a buffer for communications

            int numsteps = num_img_node / art_prm.block_size;

            int processed = 0;

            if ((num_img_node % art_prm.block_size) != 0)
                numsteps++;

            art_prm.numIMG = art_prm.block_size;

            art_prm.parallel_mode = Basic_ART_Parameters::pCAV; // Another trick

            STARTINGX(art_prm.ordered_list) = -myFirst;

            for (int ns = 0; ns < numsteps ; ns++)
            {
                if (ns == numsteps - 1)
                {
                    art_prm.numIMG = num_img_node - processed;
                }

                /*
                EXTRA CALCULATIONS FOR BICAV WEIGHTS: Each node computes its own part related to its images
                and after that send each others their results and sum up them. This part of code has been taken and
                modified from Basic_art.cc produce_side_info().
                */

                cav_t = MPI_Wtime();

                art_prm.compute_CAV_weights(vol_basis, art_prm.numIMG);
                GVNeq_aux = *(art_prm.GVNeq);

                // All processors send their result and get the other's so all of them
                // have the weights.

                for (int n = 0 ; n < (art_prm.GVNeq)->VolumesNo(); n++)
                {
                    MPI_Barrier(MPI_COMM_WORLD);
                    aux_comm_t = MPI_Wtime();
                    MPI_Allreduce(MULTIDIM_ARRAY((*(art_prm.GVNeq))(n)()), MULTIDIM_ARRAY(GVNeq_aux(n)()), MULTIDIM_SIZE((*(art_prm.GVNeq))(n)()), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                    aux_t = MPI_Wtime() - aux_comm_t;
                    comms_t += aux_t;
                    comms_t_it += aux_t;
                    (*(art_prm.GVNeq))(n)() = GVNeq_aux(n)();
                }
                cavk_it_t += MPI_Wtime() - cav_t;

                for (int jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
                    vol_aux2(jj)() = vol_basis(jj)();

                Basic_ART_iterations(art_prm, eprm, vol_basis, rank);

                if (ns < (numsteps - 1))
                {
                    STARTINGX(art_prm.ordered_list) -= art_prm.numIMG;
                    processed += art_prm.numIMG;
                }

                // All processors send their result and get the other's so all of them
                // have the same volume for the next step.

                for (int jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
                {
                    vol_basis(jj)() = vol_basis(jj)() - vol_aux2(jj)(); // Adapt result to parallel ennvironment from sequential routine
                    MPI_Barrier(MPI_COMM_WORLD);
                    aux_comm_t = MPI_Wtime();
                    MPI_Allreduce(MULTIDIM_ARRAY(vol_basis(jj)()),
                        MULTIDIM_ARRAY(vol_basis_aux(jj)()),
                        MULTIDIM_SIZE(vol_basis(jj)()),
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    aux_t = MPI_Wtime() - aux_comm_t;
                    comms_t += aux_t;
                    comms_t_it += aux_t;

                    FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_basis(jj)())
                        if (VOL_ELEM(GVNeq_aux(jj)(), k, i, j) == 0)
                        {
                            if (VOL_ELEM(vol_basis_aux(jj)(), k, i, j) == 0)
                            {
                                VOL_ELEM(vol_basis(jj)(), k, i, j) = 0;
                            }
                            else
                            {
                                // This case should not happen as this element
                                // is not affected by actual projections

                                std::cerr << "Error with weights, contact developers!" << std::endl;
                                exit(0);
                            }
                        }
                        else
                        {
                            VOL_ELEM(GVNeq_aux(jj)(), k, i, j) =
                                VOL_ELEM(vol_aux2(jj)(), k, i, j) +
                                VOL_ELEM(vol_basis_aux(jj)(), k, i, j) /
                                VOL_ELEM(GVNeq_aux(jj)(), k, i, j);
                        }
                }
            }
            art_prm.parallel_mode = Basic_ART_Parameters::pBiCAV; // trick undone
        }
        else if (art_prm.parallel_mode == Basic_ART_Parameters::pCAV)
        {

            // CAV weights calculations have been done before iterations begin in order to avoid recalculate them
            for (int j = 0 ; j < vol_basis.VolumesNo() ; j++)
                vol_aux2(j)() = vol_basis(j)();

            art_prm.numIMG = num_img_node;

            STARTINGX(art_prm.ordered_list) = -myFirst;

            Basic_ART_iterations(art_prm, eprm, vol_basis, rank);

            // All processors send their result and get the other's so all of them
            // have the same volume for the next step.
            for (int jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
            {
                vol_basis(jj)() = vol_basis(jj)() - vol_aux2(jj)(); // Adapt result to parallel ennvironment from sequential routine
                MPI_Barrier(MPI_COMM_WORLD);
                aux_comm_t = MPI_Wtime();
                MPI_Allreduce(MULTIDIM_ARRAY(vol_basis(jj)()),
                    MULTIDIM_ARRAY(vol_basis_aux(jj)()),
                    MULTIDIM_SIZE(vol_basis(jj)()),
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                aux_t = MPI_Wtime() - aux_comm_t;
                comms_t += aux_t;
                comms_t_it += aux_t;

                FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_basis(jj)())
                    if (VOL_ELEM(GVNeq_aux(jj)(), k, i, j) == 0)  // Division by 0
                    {
                        if (VOL_ELEM(vol_basis_aux(jj)(), k, i, j) == 0)
                        {
                            VOL_ELEM(vol_basis(jj)(), k, i, j) = 0;
                        }
                        else
                        {
                            // This case should not happen as this element
                            // is not affected by actual projections

                            std::cerr << "Error with weights, contact developers!" << std::endl;
                            exit(0);
                        }
                    }
                    else
                    {
                        VOL_ELEM(vol_basis(jj)(), k, i, j) =
                            VOL_ELEM(vol_aux2(jj)(), k, i, j) +
                            VOL_ELEM(vol_basis_aux(jj)(), k, i, j) /
                            VOL_ELEM(GVNeq_aux(jj)(), k, i, j);
                    }
            }
        }
        else
        {   // SIRT AND ASS

            art_prm.numIMG = num_img_node;

            STARTINGX(art_prm.ordered_list) = -myFirst;

            if (art_prm.parallel_mode == Basic_ART_Parameters::pSIRT ||
                art_prm.parallel_mode == Basic_ART_Parameters::pfSIRT)
            {
                for (int jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
                    vol_aux2(jj)() = vol_basis(jj)();
            }

            Basic_ART_iterations(art_prm, eprm, vol_basis, rank);

            // All processors send their result and get the other's so all of them
            // have the same volume for the next step.
            for (int jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
            {
                // SIRT Alg. needs to store previous results but AVSP doesn't
                if (art_prm.parallel_mode == Basic_ART_Parameters::pSIRT ||
                    art_prm.parallel_mode == Basic_ART_Parameters::pfSIRT)
                {
                    vol_basis(jj)() = vol_basis(jj)() - vol_aux2(jj)(); // Adapt result to parallel ennvironment from sequential routine
                }

                MPI_Barrier(MPI_COMM_WORLD);
                aux_comm_t = MPI_Wtime();

                MPI_Allreduce(MULTIDIM_ARRAY(vol_basis(jj)()),
                    MULTIDIM_ARRAY(vol_basis_aux(jj)()),
                    MULTIDIM_SIZE(vol_basis(jj)()),
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                aux_t = MPI_Wtime() - aux_comm_t;
                comms_t += aux_t;
                comms_t_it += aux_t;

                if (art_prm.parallel_mode == Basic_ART_Parameters::pfSIRT)
                {
                    double norm_value = (double) num_img_tot;
                    vol_basis(jj)() = vol_aux2(jj)() + (vol_basis_aux(jj)() / norm_value);
                }
                else if (art_prm.parallel_mode == Basic_ART_Parameters::pSIRT)
                {
                    double norm_value = (double) num_img_tot * (double)(art_prm.ProjXdim() * art_prm.ProjYdim());
                    vol_basis(jj)() = vol_aux2(jj)() + (vol_basis_aux(jj)() / norm_value);
                }
                else // ASS
                {
                    vol_basis(jj)() = vol_basis_aux(jj)();
                    vol_basis(jj)() /= size; // Non-SIRT Normalization
                }
            }
        }
        cavk_total_t += cavk_it_t;

        // Only the proccess with rank=0 (Master) will output the results
        if (rank == 0)
        {
            std::cout << "\nIteration " << i << std::endl;
            std::cout << "Time: " << MPI_Wtime() - it_t << " secs." << std::endl;
            std::cout << "Comms. time: " << comms_t_it << " secs." << std::endl;
            if (art_prm.parallel_mode == Basic_ART_Parameters::pBiCAV)
                std::cout << "BiCAV weighting time: " << cavk_it_t << std::endl;

            if (i < num_iter - 1)
            {
                if ((art_prm.tell&TELL_SAVE_BASIS) && (i < num_iter - 1))
                    vol_basis.write(art_prm.fn_root + "it" + integerToString(i + 1) + ".basis");
            }
        }
    }

    /*************************** FINISHING AND STORING VALUES ***************************/

    if (rank > 0)
    {
        art_prm.fh_hist->close();
        MPI_Finalize();   // Must exist for each proccess on MPI evironment
        return 0;
    }

    int Xoutput_volume_size=(art_prm.Xoutput_volume_size==0) ?
       art_prm.projXdim:art_prm.Xoutput_volume_size;
    int Youtput_volume_size=(art_prm.Youtput_volume_size==0) ?
       art_prm.projYdim:art_prm.Youtput_volume_size;
    int Zoutput_volume_size=(art_prm.Zoutput_volume_size==0) ?
       art_prm.projXdim:art_prm.Zoutput_volume_size;
    art_prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
        Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
    vol_voxels.write(art_prm.fn_root+".vol");

    if (art_prm.tell&TELL_SAVE_BASIS)
        vol_basis.write(art_prm.fn_root + ".vol");
    art_prm.fh_hist->close();

    uswtime(&recons_t);
    std::cout << "\n\n------ FINAL STATISTICS ------" << std::endl;
    std::cout << "\nTOTAL EXECUTION TIME: " << recons_t.wall - total_t << std::endl;
    std::cout << "Communications time: " << comms_t << " secs." << std::endl;
    std::cout << "CPU time: " << recons_t.user + recons_t.sys << " secs." << std::endl;
    std::cout << "USER: " << recons_t.user << " SYSTEM: " << recons_t.sys << "\n\n" << std::endl;
    if (art_prm.parallel_mode == Basic_ART_Parameters::pBiCAV)
        std::cout << "total pBiCAV Weighting time: " << cavk_total_t << std::endl;
    MPI_Finalize();
    return 0 ;
}






