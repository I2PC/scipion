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

#include <reconstruction/ml_refine3d.h>

#include <mpi.h>
#define TAG_DOCFILE 12
#define TAG_DOCFILESIZE 13
 

int main(int argc, char **argv)
{

    int                         c, iter, volno, converged = 0, argc2 = 0;
    char                        **argv2=NULL;
    double                      LL, sumw_allrefs, convv, sumcorr, wsum_sigma_noise, wsum_sigma_offset;
    std::vector<double>              conv;
    std::vector<Matrix2D<double> >   wsum_Mref;
    std::vector<double>              sumw, sumw2, sumwsc, sumwsc2, sumw_cv, sumw_mirror;
    DocFile                     DFo;

    // For parallelization
    int rank, size, num_img_tot;
    double                      aux;
    Matrix2D<double>            Maux;
    Matrix1D<double>            Vaux;
    SelFile                     SFo;

    Prog_Refine3d_prm           prm;
    Prog_MLalign2D_prm          ML2D_prm;

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
        if (rank == 0)
        {
            prm.show();
            prm.remake_SFvol(prm.istart - 1, true, false);
        }
        else prm.remake_SFvol(prm.istart - 1, false, false);
        MPI_Barrier(MPI_COMM_WORLD);

        // Read and set general MLalign2D-stuff
        ML2D_prm.read(argc2, argv2, true);
        if (rank != 0) ML2D_prm.verb = prm.verb = 0;
        if (!checkParameter(argc2, argv2, "-psi_step")) ML2D_prm.psi_step = prm.angular;
        ML2D_prm.fn_root = prm.fn_root;
        ML2D_prm.fast_mode = true;
        ML2D_prm.do_mirror = true;
        ML2D_prm.save_mem2 = true;
        ML2D_prm.write_docfile = true;
        ML2D_prm.write_selfiles = true;
        ML2D_prm.write_intermediate = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.sel";

        // Check that there are enough computing nodes
        if (prm.Nvols > size)
            REPORT_ERROR(1, "mpi_MLrefine3D requires that you use more CPUs than reference volumes");

        // Project the reference volume
        prm.project_reference_volume(ML2D_prm.SFr, rank, size);
        MPI_Barrier(MPI_COMM_WORLD);

        // All nodes produce general side-info
        ML2D_prm.produce_Side_info();
        MPI_Barrier(MPI_COMM_WORLD);

        // Select only relevant part of selfile for this rank
        ML2D_prm.SF.mpi_select_part(rank, size, num_img_tot);

        // All nodes read node-specific side-info into memory
        ML2D_prm.produce_Side_info2(prm.Nvols);
        ML2D_prm.Iold.clear(); // To save memory

        // Some output to screen
        if (rank == 0) ML2D_prm.show(true);

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
        // Initialize some additional stuff
        Maux.resize(ML2D_prm.dim, ML2D_prm.dim);
        Maux.setXmippOrigin();
        for (int refno = 0; refno < ML2D_prm.n_ref; refno++) conv.push_back(-1.);

        // Loop over all iterations
        iter = prm.istart;
        while (!converged && iter <= prm.Niter)
        {

            if (prm.verb > 0)
            {
                std::cerr        << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << std::endl;
                prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << std::endl;
            }

            // Initialize
            DFo.clear();
            conv.clear();

            if (!ML2D_prm.maxCC_rather_than_ML) ML2D_prm.calculate_pdf_phi();

            // Integrate over all images
            ML2D_prm.ML_sum_over_all_images(ML2D_prm.SF, ML2D_prm.Iref, iter,
                                            LL, sumcorr, DFo, wsum_Mref,
                                            wsum_sigma_noise, wsum_sigma_offset, 
					    sumw, sumw2, sumwsc, sumwsc2, sumw_mirror);

            // Here MPI_allreduce of all weighted sums, LL, etc.
            // All nodes need the answer to calculate internally
            // sigma_noise etc. for the next iteration!
            MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            LL = aux;
            MPI_Allreduce(&sumcorr, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sumcorr = aux;
            MPI_Allreduce(&wsum_sigma_noise, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_noise = aux;
            MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_offset = aux;
            if (ML2D_prm.do_kstest)
            {
                Vaux.resize(ML2D_prm.sumhist);
                MPI_Allreduce(MULTIDIM_ARRAY(ML2D_prm.sumhist), MULTIDIM_ARRAY(Vaux),
                              MULTIDIM_SIZE(ML2D_prm.sumhist), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                ML2D_prm.sumhist = Vaux;
            }
            for (int refno = 0;refno < ML2D_prm.n_ref; refno++)
            {
                MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[refno]), MULTIDIM_ARRAY(Maux),
                              MULTIDIM_SIZE(wsum_Mref[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                wsum_Mref[refno] = Maux;
                MPI_Allreduce(&sumw[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw[refno] = aux;
                MPI_Allreduce(&sumw_mirror[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw_mirror[refno] = aux;
                MPI_Allreduce(&sumw2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumw2[refno] = aux;
                MPI_Allreduce(&sumwsc[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumwsc[refno] = aux;
                MPI_Allreduce(&sumwsc2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                sumwsc2[refno] = aux;
            }

            // Update model parameters
            ML2D_prm.update_parameters(wsum_Mref, 
                                       wsum_sigma_noise, wsum_sigma_offset, 
				       sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, 
				       sumcorr, sumw_allrefs, prm.eachvol_end[0]+1);


            // Write intermediate files 
            if (rank != 0)
            {
                // All slaves send docfile to the master
                std::ostringstream doc;
                doc << DFo;
                int s_size=  doc.str().size();
                char results[s_size];
                strncpy(results,doc.str().c_str(),s_size);
                results[s_size]='\0';
                MPI_Send(&s_size, 1, MPI_INT, 0, TAG_DOCFILESIZE, MPI_COMM_WORLD);
                MPI_Send(results, s_size, MPI_CHAR, 0, TAG_DOCFILE, MPI_COMM_WORLD);
            }
            else
            {
                // Master fills docfile 
                std::ofstream myDocFile;
                FileName fn_tmp;
                fn_tmp.compose(prm.fn_root + "_it",iter,"doc");
                myDocFile.open (fn_tmp.c_str());
                if (ML2D_prm.maxCC_rather_than_ML)
                    myDocFile << " ; Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Corr (8)\n";
                else
                    myDocFile << " ; Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8), w_robust (9), bgmean (10), scale (11), sigma (12), KSprob (13)\n";

                // Master's own contribution
                myDocFile << DFo;
                int docCounter=1;
                while (docCounter < size)
                {
                    // receive in order
                    int iNumber, s_size;
                    MPI_Recv(&s_size, 1, MPI_INT, docCounter, TAG_DOCFILESIZE, MPI_COMM_WORLD, &status);
                    char results[s_size];
                    MPI_Recv(results, s_size, MPI_CHAR, docCounter, TAG_DOCFILE, MPI_COMM_WORLD, &status);
                    results[s_size]='\0';
                    myDocFile<<results ;
                    docCounter++;
                }

                //save doc_file and renumber it
                myDocFile.close();
                DFo.clear();
                DFo.read(fn_tmp);
                DFo.renum();

                // Output all intermediate files
                ML2D_prm.write_output_files(iter, DFo, sumw_allrefs, LL, sumcorr, conv);
                prm.concatenate_selfiles(iter);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // Reconstruct the new reference volumes also in parallel
            // Assume that the number of processors is larger than the
            // number of volumes to reconstruct ...
            if (rank < prm.Nvols)
                // new reference reconstruction
                prm.reconstruction(argc2, argv2, iter, rank, 0);
            MPI_Barrier(MPI_COMM_WORLD);

            // Only the master does post-processing & convergence check (i.e. sequentially)
            if (rank == 0)
            {

                // Solvent flattening and/or symmetrization (if requested)
                prm.remake_SFvol(iter, false, false);
                prm.post_process_volumes(argc2, argv2);

                // Check convergence
                if (prm.check_convergence(iter))
                {
                    converged = 1;
                    if (prm.verb > 0) std::cerr << "--> Optimization converged!" << std::endl;
                }

            }

            // Broadcast new spectral_signal and converged to all nodes
            MPI_Bcast(&converged, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // Update filenames in SFvol (now without noise volumes!)
            prm.remake_SFvol(iter, false, false);

            if (!converged && iter + 1 <= prm.Niter)
            {
                // All nodes again: project and read new references from disc
                prm.project_reference_volume(ML2D_prm.SFr, rank, size);
                MPI_Barrier(MPI_COMM_WORLD);
                ML2D_prm.SFr.go_beginning();
                c = 0;
                while (!ML2D_prm.SFr.eof())
                {
                    FileName fn_img=ML2D_prm.SFr.NextImg();
                    if (ML2D_prm.SFr.eof()) break;
                    ML2D_prm.Iref[c].read(fn_img, false, false, false, false);
                    ML2D_prm.Iref[c]().setXmippOrigin();
                    c++;
                }
            }

            iter++;
        } // end loop iterations

	if (rank == 0)
	    ML2D_prm.write_output_files(-1, DFo, sumw_allrefs, LL, sumcorr, conv);

        if (!converged && prm.verb > 0)
            std::cerr << "--> Optimization was stopped before convergence was reached!" << std::endl;

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




