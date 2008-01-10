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
#include <reconstruction/mlp_align2d.h>

#include <mpi.h>

int main(int argc, char **argv)
{

    int                         c, iter, volno, converged = 0;
    double                      LL, sumw_allrefs, convv, sumcorr, wsum_sigma_noise, wsum_sigma_offset;
    vector<double>              conv;
    vector< Polar <complex<double> > > fP_wsum_imgs;
    vector<double>              sumw, sumw_cv, sumw_mirror;
    Matrix1D<double>            Vaux, Vsum;
    DocFile                     DFo;

    // For parallelization
    int rank, size, num_img_tot;
    double                      aux;
    Matrix2D<double>            Maux;
    FileName                    fn_tmp;
    SelFile                     SFo;

    Prog_Refine3d_prm           prm;
    Prog_MLPalign2D_prm         ML2D_prm;

    // Init Parallel interface
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    // Get input parameters
    try
    {

        // Read command line
        prm.read(argc, argv);

        // Write starting volumes to disc with correct name for iteration loop
        if (rank == 0)
        {
            prm.show();
            prm.remake_SFvol(prm.istart - 1, true, false);
        }
        else prm.remake_SFvol(prm.istart - 1, false, false);
        MPI_Barrier(MPI_COMM_WORLD);

        // Read and set general MLalign2D-stuff
        ML2D_prm.read(argc, argv, true);
        if (rank != 0) ML2D_prm.verb = prm.verb = 0;
        ML2D_prm.fn_root = prm.fn_root;
        ML2D_prm.do_mirror = true;
        ML2D_prm.fn_ref = prm.fn_root + "_lib.sel";

        // Check that there are enough computing nodes
        if (prm.Nvols > size)
            REPORT_ERROR(1, "mpi_MLrefine3D requires that you use more CPUs than reference volumes");

        // Project the reference volume
        prm.project_reference_volume(ML2D_prm.SFr, rank);
        MPI_Barrier(MPI_COMM_WORLD);

        // All nodes produce general side-info
        ML2D_prm.produceSideInfo();

        // Select only relevant part of selfile for this rank
        ML2D_prm.SF.mpi_select_part(rank, size, num_img_tot);

        // All nodes read node-specific side-info into memory
        ML2D_prm.produceSideInfo2(prm.Nvols);
        ML2D_prm.Iold.clear(); // To save memory

        // Some output to screen
        if (rank == 0) ML2D_prm.show(true);

    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            cout << XE;
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
        for (int refno = 0; refno < ML2D_prm.nr_ref; refno++) conv.push_back(-1.);

        // Loop over all iterations
        iter = prm.istart;
        while (!converged && iter <= prm.Niter)
        {

            if (prm.verb > 0)
            {
                cerr        << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << endl;
                prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << iter << " of " << prm.Niter << endl;
            }

            // Prepare DFo header
            DFo.clear();
            conv.clear();
            if (rank == 0)
            {
		DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");
            }

            // Pre-calculate pdfs
            ML2D_prm.updatePdfTranslations();

            // Integrate over all images
            ML2D_prm.sumOverAllImages(ML2D_prm.SF, ML2D_prm.Iref,
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

            for (int refno = 0;refno < ML2D_prm.nr_ref; refno++)
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
            ML2D_prm.updateParameters(fP_wsum_imgs, wsum_sigma_noise, wsum_sigma_offset,
				      sumw, sumw_mirror, sumcorr, sumw_allrefs);

            // All nodes write out temporary DFo
            fn_tmp.compose(prm.fn_root, rank, "tmpdoc");
            DFo.write(fn_tmp);
            MPI_Barrier(MPI_COMM_WORLD);

            // Only the master outputs intermediate files
            if (rank == 0)
            {
                DFo.clear();
                for (int rank2 = 0; rank2 < size; rank2++)
                {
                    fn_tmp.compose(prm.fn_root, rank2, "tmpdoc");
                    DFo.append(fn_tmp);
                    DFo.locate(DFo.get_last_key());
                    DFo.next();
                    DFo.remove_current();
                    system(((string)"rm -f " + fn_tmp).c_str());
                }
		ML2D_prm.writeOutputFiles(iter, DFo, sumw_allrefs, LL, sumcorr, conv);
                prm.concatenate_selfiles(iter);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // Reconstruct the new reference volumes also in parallel
            // Assume that the number of processors is larger than the
            // number of volumes to reconstruct ...
            if (rank < prm.Nvols)
                // new reference reconstruction
                prm.reconstruction(argc, argv, iter, rank, 0);
            MPI_Barrier(MPI_COMM_WORLD);

            // Only the master does post-processing & convergence check (i.e. sequentially)
            if (rank == 0)
            {

                // Solvent flattening and/or symmetrization (if requested)
                prm.remake_SFvol(iter, false, false);
                prm.post_process_volumes(argc, argv);

                // Check convergence
                if (prm.check_convergence(iter))
                {
                    converged = 1;
                    if (prm.verb > 0) cerr << "--> Optimization converged!" << endl;
                }

            }

            // Broadcast new spectral_signal and converged to all nodes
            MPI_Bcast(&converged, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // Update filenames in all SFvol
            prm.remake_SFvol(iter, false, false);

            if (!converged && iter + 1 <= prm.Niter)
            {
                // All nodes again: project and read new references from disc
                prm.project_reference_volume(ML2D_prm.SFr, rank);
                MPI_Barrier(MPI_COMM_WORLD);
                ML2D_prm.SFr.go_beginning();
                c = 0;
                while (!ML2D_prm.SFr.eof())
                {
                    ML2D_prm.Iref[c].read(ML2D_prm.SFr.NextImg(), false, false, false, false);
                    ML2D_prm.Iref[c]().setXmippOrigin();
                    c++;
                }
            }

            iter++;
        } // end loop iterations

	if (rank == 0)
	    ML2D_prm.writeOutputFiles(-1, DFo, sumw_allrefs, LL, sumcorr, conv);

        if (!converged && prm.verb > 0)
            cerr << "--> Optimization was stopped before convergence was reached!" << endl;

    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    MPI_Finalize();
    return 0;

}




