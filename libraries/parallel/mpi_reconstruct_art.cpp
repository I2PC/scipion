/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include "mpi_reconstruct_art.h"
#include <sys/resource.h>

/** Empty constructor */
ProgMPIReconsArt::ProgMPIReconsArt()
{
    isMpi = true;
}

/*  constructor ------------------------------------------------------- */
ProgMPIReconsArt::ProgMPIReconsArt(int argc, char *argv[])
{
    this->read(argc, argv);
}

/* constructor providing an MpiNode
 * this is useful for using this program from others
 */
ProgMPIReconsArt::ProgMPIReconsArt(MpiNode * node)
{
    this->setNode(node);
}


/* Run --------------------------------------------------------------------- */
void ProgMPIReconsArt::run()
{
    BasicARTParameters &artPrm = artRecons->artPrm;

    Image<double> vol_voxels, vol_voxels_aux; // Volume to reconstruct
    GridVolume    vol_basis;

    int     num_img_tot;  // The total amount of images
    int     num_img_node;  // Stores the total amount of images each node deals with
    size_t  remaining;          // Number of projections still to compute
    int      Npart;              // Number of projection to process
    int      myFirst;     // Limits projections to process
    double    comms_t, aux_comm_t; // Communications time
    double    it_t;   // iteration time
    double     cav_t;  // time for CAV weights calculation
    double     cavk_it_t;          // BiCAV weights calculation time (1 iter.)
    double     cavk_total_t;       // Sum( cavk_it_t )
    USWtime_t  recons_t;  // Reconstruction time
    double    total_t;  // Program execution time
    double     comms_t_it, aux_t; // Communications time in one iteration
    GridVolumeT<int> GVNeq_aux;    // This is a buffer for the communication
    Matrix1D<int> Ordered_aux;


    /*************************** PARAMETERS INITIALIZATION ***************************/
    cavk_total_t = 0.0;

    if (node->isMaster())  // Master code
    {
        total_t = MPI_Wtime();   // to perform overall execution time
        comms_t = 0.0; // Initializes time

        show();
        artRecons->preProcess(vol_basis, FULL, node->rank);
        if (verbose > 0)
        {
            std::cout << " ---------------------------------------------------------------------" << std::endl;
            std::cout << " Projections                  : "  << artRecons->artPrm.numIMG << std::endl;
        }
        //Initialize history
        artRecons->initHistory(vol_basis);

        // ordered list must be the same in all nodes
        aux_comm_t = MPI_Wtime();
        MPI_Bcast(MULTIDIM_ARRAY(artPrm.ordered_list), MULTIDIM_SIZE(artPrm.ordered_list), MPI_INT, 0, MPI_COMM_WORLD);
        comms_t += MPI_Wtime() - aux_comm_t;
    }
    else
    {
        // each process can handle its own history file
        // so we add the id number to the root filename
        FileName aux = artPrm.fn_root;

        artPrm.fn_root = artPrm.fn_root + integerToString(node->rank);
        artRecons->preProcess(vol_basis, FULL, node->rank);

        // Restore original filename.
        artPrm.fn_root = aux;

        // ordered list must be the same in all nodes
        MPI_Bcast(MULTIDIM_ARRAY(artPrm.ordered_list),
                  MULTIDIM_SIZE(artPrm.ordered_list), MPI_INT, 0, MPI_COMM_WORLD);
    }

    // How many images processes each node. It is done in such a way that every node receives the
    // same amount of them

    num_img_tot = artPrm.numIMG;

    Npart = (int)((float)num_img_tot / (float)nProcs);

    num_img_node = Npart;

    remaining = num_img_tot % nProcs;

    // each process will only see the images it is interested in.
    if (node->rank < remaining)
    {
        num_img_node++;
        myFirst = node->rank * Npart + node->rank;
    }
    else
        myFirst = Npart * node->rank + remaining;

    // Shift the starting position for each node
//    myLast = myFirst + num_img_node - 1; // Actually it isn't used

    STARTINGX(artPrm.ordered_list) = -myFirst;

    GridVolume   vol_basis_aux = vol_basis;
    GridVolume   vol_aux2 = vol_basis;

    if (artPrm.parallel_mode == BasicARTParameters::pSART ||
        artPrm.parallel_mode == BasicARTParameters::pBiCAV)
    {
        if (artPrm.block_size < 1)
            artPrm.block_size = 1;  // Each processor will have at least one projection
        else if (artPrm.block_size > num_img_node)
            artPrm.block_size = num_img_node;  // block_size is as much equal to num_img_node
    }

    // Print some data
    if (verbose > 0)
    {
        std::cout << " Parallel mode                : "  ;
        if (artPrm.parallel_mode == BasicARTParameters::pSART)
            std::cout << "pSART " << "BlockSize = " << artPrm.block_size << std::endl;
        else if (artPrm.parallel_mode == BasicARTParameters::pSIRT)
            std::cout << "pSIRT" << std::endl;
        else if (artPrm.parallel_mode == BasicARTParameters::pfSIRT)
            std::cout << "pfSIRT" << std::endl;
        else if (artPrm.parallel_mode == BasicARTParameters::pBiCAV)
            std::cout << "pBiCAV " << "BlockSize = " << artPrm.block_size << std::endl;
        else if (artPrm.parallel_mode == BasicARTParameters::pAVSP)
            std::cout << "pAVSP" << std::endl;
        else if (artPrm.parallel_mode == BasicARTParameters::pCAV)
            std::cout << "pCAV" << std::endl;
        else if (artPrm.parallel_mode == BasicARTParameters::ART)
            std::cout << "ART" << std::endl;
        else if (artPrm.parallel_mode == BasicARTParameters::SIRT)
            std::cout << "SIRT" << std::endl;
        std::cout << " Number of processors         : " << nProcs << std::endl;
        std::cout << " ---------------------------------------------------------------------" << std::endl;
    }


    /*************************** CAV weights precalculation *************************/
    if (artPrm.parallel_mode == BasicARTParameters::pCAV)
    {
        // Creates and initializes special variables needed to CAV weights computation.

        /* EXTRA CALCULATIONS FOR CAV WEIGHTS: Each node computes its own part related to its images
        and after that they send each other their results and sum up them. This part of code has been taken and
        modified from Basic_art.cc produceSideInfo().*/

        cav_t = MPI_Wtime();
        artPrm.computeCAVWeights(vol_basis, num_img_node, verbose-1);
        GVNeq_aux = *(artPrm.GVNeq);

        MPI_Barrier(MPI_COMM_WORLD); // Actually, this isn't necessary.
        for (size_t n = 0 ; n < (artPrm.GVNeq)->VolumesNo(); n++)
        {
            aux_comm_t = MPI_Wtime();
            MPI_Allreduce(MULTIDIM_ARRAY((*(artPrm.GVNeq))(n)()),
                          MULTIDIM_ARRAY(GVNeq_aux(n)()),
                          MULTIDIM_SIZE((*(artPrm.GVNeq))(n)()), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            comms_t += MPI_Wtime() - aux_comm_t;
            (*(artPrm.GVNeq))(n)() = GVNeq_aux(n)();
        }
        if (verbose > 0)
            std::cout << "Elapsed time for pCAV weights computation : " << MPI_Wtime() - cav_t << std::endl;
    }

    /*************************** PARALLEL ART ALGORITHM ***************************/
    // A bit tricky. Allows us to use the sequential Basic_ART_iterations as
    // in parallel it runs internally only one iteration.
    int num_iter = artPrm.no_it;
    artPrm.no_it = 1;

    for (int i = 0; i < num_iter; i++)
    {
        comms_t_it = 0.0;
        cavk_it_t  = 0.0;

        it_t = MPI_Wtime();

        artPrm.lambda_list(0) = artPrm.lambda(i);

        if (artPrm.parallel_mode == BasicARTParameters::pSART)
        {
            int numsteps = num_img_node / artPrm.block_size;

            // could be necessary another step for remaining projections
            int stepR;
            bool oneMoreStep = (stepR = num_img_node % artPrm.block_size) != 0;
            if (oneMoreStep)
                numsteps++;

            artPrm.numIMG = artPrm.block_size;
            int blockSizeTot = artPrm.block_size * nProcs;

            //            STARTINGX(artPrm.ordered_list) = -myFirst; //Already assigned
            for (int ns = 0; ns < numsteps ; ns++)
            {
                // Calculus of blockSizeTot to normalize
                if (ns == numsteps - 2 && oneMoreStep)
                    blockSizeTot = (artPrm.block_size - 1) * nProcs + stepR;
                else if (ns == numsteps - 1)
                {
                    artPrm.numIMG = num_img_node - artPrm.block_size*ns;
                    if (oneMoreStep)
                        blockSizeTot = artPrm.numIMG * stepR;
                    else
                        blockSizeTot = artPrm.block_size*stepR + artPrm.numIMG*(nProcs - stepR);
                }

                for (size_t j = 0 ; j < vol_basis.VolumesNo() ; j++)
                    vol_aux2(j)() = vol_basis(j)();

                artRecons->iterations(vol_basis, node->rank);

                STARTINGX(artPrm.ordered_list) -= artPrm.numIMG;

//                node->updateComm();//Update communicator to avoid already finished nodes

                // All processors send their result and get the other's so all of them
                // have the same volume for the next step.
                for (size_t j = 0 ; j < vol_basis.VolumesNo() ; j++)
                {
                    vol_basis(j)() = vol_basis(j)() - vol_aux2(j)(); // Adapt result to parallel environment from sequential routine
                    aux_comm_t = MPI_Wtime();
                    MPI_Allreduce(MULTIDIM_ARRAY(vol_basis(j)()),
                                  MULTIDIM_ARRAY(vol_basis_aux(j)()),
                                  MULTIDIM_SIZE(vol_basis(j)()),
                                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    aux_t = MPI_Wtime() - aux_comm_t;
                    comms_t += aux_t;
                    comms_t_it += aux_t;
                    vol_basis(j)() = vol_aux2(j)() + (vol_basis_aux(j)() / (double) blockSizeTot);
                }
            }
        }
        else if (artPrm.parallel_mode == BasicARTParameters::pBiCAV)
        {
            // Creates and initializes special variables needed to BICAV weights computation.
            GridVolumeT<int> GVNeq_aux; // This is a buffer for communications

            int numsteps = num_img_node / artPrm.block_size;

            int processed = 0;

            if ((num_img_node % artPrm.block_size) != 0)
                numsteps++;

            artPrm.numIMG = artPrm.block_size;

            artPrm.parallel_mode = BasicARTParameters::pCAV; // Another trick

            STARTINGX(artPrm.ordered_list) = -myFirst;

            for (int ns = 0; ns < numsteps ; ns++)
            {
                if (ns == numsteps - 1)
                    artPrm.numIMG = num_img_node - processed;

                /*
                EXTRA CALCULATIONS FOR BICAV WEIGHTS: Each node computes its own part related to its images
                and after that send each others their results and sum up them. This part of code has been taken and
                modified from Basic_art.cc produceSideInfo().
                */

                cav_t = MPI_Wtime();

                artPrm.computeCAVWeights(vol_basis, artPrm.numIMG, verbose-1);
                GVNeq_aux = *(artPrm.GVNeq);

//                node->updateComm();
                // All processors send their result and get the other's so all of them
                // have the weights.

                for (size_t n = 0 ; n < (artPrm.GVNeq)->VolumesNo(); n++)
                {
                    aux_comm_t = MPI_Wtime();
                    MPI_Allreduce(MULTIDIM_ARRAY((*(artPrm.GVNeq))(n)()),
                                  MULTIDIM_ARRAY(GVNeq_aux(n)()),
                                  MULTIDIM_SIZE((*(artPrm.GVNeq))(n)()), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                    aux_t = MPI_Wtime() - aux_comm_t;
                    comms_t += aux_t;
                    comms_t_it += aux_t;
                    (*(artPrm.GVNeq))(n)() = GVNeq_aux(n)();
                }
                cavk_it_t += MPI_Wtime() - cav_t;

                for (size_t jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
                    vol_aux2(jj)() = vol_basis(jj)();

                artRecons->iterations(vol_basis, node->rank);

                if (ns < (numsteps - 1))
                {
                    STARTINGX(artPrm.ordered_list) -= artPrm.numIMG;
                    processed += artPrm.numIMG;
                }

                // All processors send their result and get the other's so all of them
                // have the same volume for the next step.

                for (size_t jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
                {
                    vol_basis(jj)() = vol_basis(jj)() - vol_aux2(jj)(); // Adapt result to parallel environment from sequential routine
                    aux_comm_t = MPI_Wtime();
                    MPI_Allreduce(MULTIDIM_ARRAY(vol_basis(jj)()),
                                  MULTIDIM_ARRAY(vol_basis_aux(jj)()),
                                  MULTIDIM_SIZE(vol_basis(jj)()),
                                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    aux_t = MPI_Wtime() - aux_comm_t;
                    comms_t += aux_t;
                    comms_t_it += aux_t;

                    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_basis(jj)())
                    if (A3D_ELEM(GVNeq_aux(jj)(), k, i, j) == 0)
                    {
                        if (A3D_ELEM(vol_basis_aux(jj)(), k, i, j) == 0)
                            A3D_ELEM(vol_basis(jj)(), k, i, j) = 0;
                        else // This case should not happen as this element is not affected by actual projections
                            REPORT_ERROR(ERR_VALUE_INCORRECT,"Error with weights, contact developers!");
                    }
                    else
                    {
                        A3D_ELEM(vol_basis(jj)(), k, i, j) =
                            A3D_ELEM(vol_aux2(jj)(), k, i, j) +
                            A3D_ELEM(vol_basis_aux(jj)(), k, i, j) /
                            A3D_ELEM(GVNeq_aux(jj)(), k, i, j);
                    }
                }
            }
            artPrm.parallel_mode = BasicARTParameters::pBiCAV; // trick undone
        }
        else if (artPrm.parallel_mode == BasicARTParameters::pCAV)
        {
            // CAV weights calculations have been done before iterations begin in order to avoid recalculate them
            for (size_t j = 0 ; j < vol_basis.VolumesNo() ; j++)
                vol_aux2(j)() = vol_basis(j)();

            artPrm.numIMG = num_img_node;

            STARTINGX(artPrm.ordered_list) = -myFirst;

            artRecons->iterations(vol_basis, node->rank);

            // All processors send their result and get the other's so all of them
            // have the same volume for the next step.
            MPI_Barrier(MPI_COMM_WORLD);
            for (size_t jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
            {
                vol_basis(jj)() = vol_basis(jj)() - vol_aux2(jj)(); // Adapt result to parallel ennvironment from sequential routine
                aux_comm_t = MPI_Wtime();
                MPI_Allreduce(MULTIDIM_ARRAY(vol_basis(jj)()),
                              MULTIDIM_ARRAY(vol_basis_aux(jj)()),
                              MULTIDIM_SIZE(vol_basis(jj)()),
                              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                aux_t = MPI_Wtime() - aux_comm_t;
                comms_t += aux_t;
                comms_t_it += aux_t;

                FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_basis(jj)())
                if (A3D_ELEM(GVNeq_aux(jj)(), k, i, j) == 0)  // Division by 0
                {
                    if (A3D_ELEM(vol_basis_aux(jj)(), k, i, j) == 0)
                        A3D_ELEM(vol_basis(jj)(), k, i, j) = 0;
                    else // This case should not happen as this element is not affected by actual projections
                        REPORT_ERROR(ERR_VALUE_INCORRECT,"Error with weights, contact developers!");
                }
                else
                {
                    A3D_ELEM(vol_basis(jj)(), k, i, j) =
                        A3D_ELEM(vol_aux2(jj)(), k, i, j) +
                        A3D_ELEM(vol_basis_aux(jj)(), k, i, j) /
                        A3D_ELEM(GVNeq_aux(jj)(), k, i, j);
                }
            }
        }
        else
        {   // SIRT AND ASS

            artPrm.numIMG = num_img_node;

            STARTINGX(artPrm.ordered_list) = -myFirst;

            if (artPrm.parallel_mode == BasicARTParameters::pSIRT ||
                artPrm.parallel_mode == BasicARTParameters::pfSIRT)
            {
                for (size_t jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
                    vol_aux2(jj)() = vol_basis(jj)();
            }

            artRecons->iterations(vol_basis, node->rank);

            // All processors send their result and get the other's so all of them
            // have the same volume for the next step.
            MPI_Barrier(MPI_COMM_WORLD);
            for (size_t jj = 0 ; jj < vol_basis.VolumesNo() ; jj++)
            {
                // SIRT Alg. needs to store previous results but AVSP doesn't
                if (artPrm.parallel_mode == BasicARTParameters::pSIRT ||
                    artPrm.parallel_mode == BasicARTParameters::pfSIRT)
                    vol_basis(jj)() = vol_basis(jj)() - vol_aux2(jj)(); // Adapt result to parallel ennvironment from sequential routine

                aux_comm_t = MPI_Wtime();

                MPI_Allreduce(MULTIDIM_ARRAY(vol_basis(jj)()),
                              MULTIDIM_ARRAY(vol_basis_aux(jj)()),
                              MULTIDIM_SIZE(vol_basis(jj)()),
                              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                aux_t = MPI_Wtime() - aux_comm_t;
                comms_t += aux_t;
                comms_t_it += aux_t;

                if (artPrm.parallel_mode == BasicARTParameters::pfSIRT)
                {
                    double norm_value = (double) num_img_tot;
                    vol_basis(jj)() = vol_aux2(jj)() + (vol_basis_aux(jj)() / norm_value);
                }
                else if (artPrm.parallel_mode == BasicARTParameters::pSIRT)
                {
                    double norm_value = (double) num_img_tot * (double)(artPrm.ProjXdim() * artPrm.ProjYdim());
                    vol_basis(jj)() = vol_aux2(jj)() + (vol_basis_aux(jj)() / norm_value);
                }
                else // ASS
                {
                    vol_basis(jj)() = vol_basis_aux(jj)();
                    vol_basis(jj)() /= nProcs; // Non-SIRT Normalization
                }
            }
        }
        cavk_total_t += cavk_it_t;

        // Only the Master process will output the results (it's the only with verbose!=0)
        if (verbose)
        {
            std::cout << "\nIteration " << i << std::endl;
            std::cout << "Time: " << MPI_Wtime() - it_t << " secs." << std::endl;
            std::cout << "Comms. time: " << comms_t_it << " secs." << std::endl;
            if (artPrm.parallel_mode == BasicARTParameters::pBiCAV)
                std::cout << "BiCAV weighting time: " << cavk_it_t << std::endl;
        }
        if (node->isMaster() && (artPrm.tell&TELL_SAVE_BASIS) && (i < num_iter - 1))
            vol_basis.write(artPrm.fn_root + "_it" + integerToString(i + 1) + ".basis");
    }

    /*************************** FINISHING AND STORING VALUES ***************************/

    if (node->isMaster())
    {
        int Xoutput_volume_size=(artPrm.Xoutput_volume_size==0) ?
                                artPrm.projXdim:artPrm.Xoutput_volume_size;
        int Youtput_volume_size=(artPrm.Youtput_volume_size==0) ?
                                artPrm.projYdim:artPrm.Youtput_volume_size;
        int Zoutput_volume_size=(artPrm.Zoutput_volume_size==0) ?
                                artPrm.projXdim:artPrm.Zoutput_volume_size;
        artPrm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                    Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
        vol_voxels.write(artPrm.fn_root+".vol");

        if (artPrm.tell&TELL_SAVE_BASIS)
            vol_basis.write(artPrm.fn_root + ".vol");

        uswtime(&recons_t);
        if (verbose)
        {
            std::cout << "\n\n------ FINAL STATISTICS ------" << std::endl;
            std::cout << "\nTOTAL EXECUTION TIME: " << recons_t.wall - total_t << std::endl;
            std::cout << "Communications time: " << comms_t << " secs." << std::endl;
            std::cout << "CPU time: " << recons_t.user + recons_t.sys << " secs." << std::endl;
            std::cout << "USER: " << recons_t.user << " SYSTEM: " << recons_t.sys << std::endl;
            if (artPrm.parallel_mode == BasicARTParameters::pBiCAV)
                std::cout << "\nTotal pBiCAV Weighting time: " << cavk_total_t << std::endl;
        }
    }
    artPrm.fh_hist->close();
}


/* ------------------------------------------------------------------------- */
/* Time managing stuff                                                       */
/* ------------------------------------------------------------------------- */
// Gets User and System times for use with MPI
void uswtime(USWtime_t *tm)
{
    struct rusage buffer;

    tm->wall = MPI_Wtime();
    getrusage(RUSAGE_SELF, &buffer);
    tm->user = (double) buffer.ru_utime.tv_sec + 1.0e-6 * buffer.ru_utime.tv_usec;
    tm->sys  = (double) buffer.ru_stime.tv_sec + 1.0e-6 * buffer.ru_stime.tv_usec;
}


