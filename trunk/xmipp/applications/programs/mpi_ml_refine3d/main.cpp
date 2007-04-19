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

int main(int argc, char **argv) {

  int                         c,iter,volno,converged=0;;
  double                      LL,sumw_allrefs,convv,sumcorr,wsum_sigma_noise, wsum_sigma_offset;
  vector<double>              conv;
  vector<matrix2D<double> >   wsum_Mref,wsum_ctfMref,Mwsum_sigma2;
  vector<double>              sumw,sumw_cv,sumw_mirror;
  matrix1D<double>            spectral_signal;
  DocFile                     DFo;

  // For parallelization
  int rank, size, num_img_tot;
  double                      aux;
  matrix2D<double>            Maux;
  FileName                    fn_tmp;
  SelFile                     SFo;

  Prog_Refine3d_prm           prm;
  Prog_MLalign2D_prm          ML2D_prm;

  // Init Parallel interface		
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;

  // Get input parameters
  try {

    // Read command line
    prm.read(argc,argv);

    // Write starting volumes to disc with correct name for iteration loop
    if (rank==0) {
      prm.show();
      prm.remake_SFvol(prm.istart-1,true,false);
    } else prm.remake_SFvol(prm.istart-1,false,false);
    MPI_Barrier(MPI_COMM_WORLD);

    // Read and set general MLalign2D-stuff
    ML2D_prm.read(argc,argv,true);
    if (rank!=0) ML2D_prm.verb=prm.verb=0;
    if (!check_param(argc,argv,"-psi_step")) ML2D_prm.psi_step=prm.angular;
    ML2D_prm.fn_root=prm.fn_root;
    ML2D_prm.fast_mode=true;
    ML2D_prm.do_mirror=true;
    ML2D_prm.save_mem2=true;
    ML2D_prm.write_docfile=true;
    ML2D_prm.write_selfiles=true;
    ML2D_prm.write_intermediate=true;
    ML2D_prm.fn_ref=prm.fn_root+"_lib.sel";

    // Check that there are enough computing nodes
    if (prm.Nvols>size)
      REPORT_ERROR(1,"mpi_MLrefine3D requires that you use more CPUs than reference volumes");
    if (ML2D_prm.fourier_mode && 3*prm.Nvols>size)
      REPORT_ERROR(1,"mpi_mlf_refine3d requires that you use three times more CPUs than reference volumes");

    // Project the reference volume
    prm.project_reference_volume(ML2D_prm.SFr,rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // All nodes produce general side-info
    ML2D_prm.produce_Side_info();
    if (ML2D_prm.fourier_mode && rank==0) ML2D_prm.estimate_initial_sigma2();
    MPI_Barrier(MPI_COMM_WORLD);

    // Select only relevant part of selfile for this rank
    ML2D_prm.SF.mpi_select_part(rank,size,num_img_tot);

    // All nodes read node-specific side-info into memory
    ML2D_prm.produce_Side_info2(prm.Nvols);
    ML2D_prm.Iold.clear(); // To save memory

    // Some output to screen
    if (rank==0) ML2D_prm.show(true);

  } catch (Xmipp_error XE) {
    if (rank==0) {
      cout << XE;
      if (prm.fourier_mode) prm.MLF_usage();
      else prm.usage();
    }
    MPI_Finalize(); exit(1);
  }

  try {
    // Initialize some additional stuff
    Maux.resize(ML2D_prm.dim,ML2D_prm.dim);
    Maux.set_Xmipp_origin();
    for (int refno=0; refno<ML2D_prm.n_ref; refno++) conv.push_back(-1.);

    // Loop over all iterations
    iter=prm.istart;
    while (!converged && iter<=prm.Niter ) {

      if (prm.verb>0) {
	cerr        << "--> 3D-EM volume refinement:  iteration " << iter <<" of "<< prm.Niter<<endl;
	prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << iter <<" of "<< prm.Niter<<endl;
      }

      // Prepare DFo header
      DFo.clear();
      conv.clear();
      if (rank==0) {
	if (ML2D_prm.maxCC_rather_than_ML)
	  DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Corr (8)");
	else
	  DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");
      }

      // Pre-calculate pdfs
      if (!ML2D_prm.maxCC_rather_than_ML) ML2D_prm.calculate_pdf_phi();

      // Integrate over all images
      ML2D_prm.ML_sum_over_all_images(ML2D_prm.SF,ML2D_prm.Iref,iter,
				      LL,sumcorr,DFo,wsum_Mref,wsum_ctfMref,
				      wsum_sigma_noise,Mwsum_sigma2,
				      wsum_sigma_offset,sumw,sumw_mirror);

      // Here MPI_allreduce of all weighted sums, LL, etc.
      // All nodes need the answer to calculate internally
      // sigma_noise etc. for the next iteration!
      MPI_Allreduce(&LL,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      LL=aux;
      MPI_Allreduce(&sumcorr,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      sumcorr=aux;
      MPI_Allreduce(&wsum_sigma_noise,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      wsum_sigma_noise=aux;
      MPI_Allreduce(&wsum_sigma_offset,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      wsum_sigma_offset=aux;
      if (ML2D_prm.fourier_mode) {
        for (int ifocus=0;ifocus<ML2D_prm.nr_focus;ifocus++) {
          MPI_Allreduce(MULTIDIM_ARRAY(Mwsum_sigma2[ifocus]),MULTIDIM_ARRAY(Maux),
                        MULTIDIM_SIZE(Mwsum_sigma2[ifocus]),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
          Mwsum_sigma2[ifocus]=Maux;
        }
	for (int refno=0;refno<ML2D_prm.n_ref; refno++) {
	  MPI_Allreduce(MULTIDIM_ARRAY(wsum_ctfMref[refno]),MULTIDIM_ARRAY(Maux),
			MULTIDIM_SIZE(wsum_ctfMref[refno]),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  wsum_ctfMref[refno]=Maux;
	}
      }
      for (int refno=0;refno<ML2D_prm.n_ref; refno++) {
	MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[refno]),MULTIDIM_ARRAY(Maux),
		      MULTIDIM_SIZE(wsum_Mref[refno]),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	wsum_Mref[refno]=Maux;
	MPI_Allreduce(&sumw[refno],&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	sumw[refno]=aux;
	MPI_Allreduce(&sumw_mirror[refno],&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	sumw_mirror[refno]=aux;
      }

      // Update model parameters
      ML2D_prm.update_parameters(wsum_Mref,wsum_ctfMref,
				 wsum_sigma_noise,Mwsum_sigma2,
				 wsum_sigma_offset,sumw,
				 sumw_mirror,sumcorr,sumw_allrefs,
				 spectral_signal);

      // All nodes write out temporary DFo
      fn_tmp.compose(prm.fn_root,rank,"tmpdoc");
      DFo.write(fn_tmp);
      MPI_Barrier(MPI_COMM_WORLD);

      // Only the master outputs intermediate files
      if (rank==0) {
	DFo.clear();
	for (int rank2=0; rank2<size; rank2++) {
	  fn_tmp.compose(prm.fn_root,rank2,"tmpdoc");
	  DFo.append(fn_tmp);
	  DFo.locate(DFo.get_last_key());
	  DFo.next();
	  DFo.remove_current();
	  system(((string)"rm -f "+fn_tmp).c_str());
	}
	ML2D_prm.write_output_files(iter,DFo,sumw_allrefs,LL,sumcorr,conv);
	prm.concatenate_selfiles(iter);

	// Write noise images to disc
	if (ML2D_prm.fourier_mode) prm.make_noise_images(ML2D_prm.Iref);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      // Reconstruct the new reference volumes also in parallel
      // Assume that the number of processors is larger than the
      // number of volumes to reconstruct ...
      if (rank<prm.Nvols)
        // new reference reconstruction
        prm.reconstruction(argc, argv, iter, rank, 0);
      else if (ML2D_prm.fourier_mode && rank>=prm.Nvols && rank<2*prm.Nvols)
	// noise reconstruction
	prm.reconstruction(argc, argv, iter, rank%prm.Nvols, 1);
      else if (ML2D_prm.fourier_mode && rank>=2*prm.Nvols && rank<3*prm.Nvols)
	// ctf-corrupted reconstruction
	prm.reconstruction(argc, argv, iter, rank%prm.Nvols, 2);
      MPI_Barrier(MPI_COMM_WORLD);

      // Only the master does post-processing & convergence check (i.e. sequentially)
      if (rank==0) {

	// Solvent flattening and/or symmetrization (if requested)
	prm.remake_SFvol(iter,false,ML2D_prm.fourier_mode);
	prm.post_process_volumes(argc, argv);

	// Calculate 3D-SSNR
	if (ML2D_prm.fourier_mode) prm.calculate_3DSSNR(spectral_signal,iter);

	// Check convergence
	if (prm.check_convergence(iter)) {
	  converged=1;
	  if (prm.verb>0) cerr <<"--> Optimization converged!"<<endl;
	}

      }

      // Broadcast new spectral_signal and converged to all nodes
      MPI_Bcast(&converged, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(MULTIDIM_ARRAY(spectral_signal),MULTIDIM_SIZE(spectral_signal),
                MPI_DOUBLE,0,MPI_COMM_WORLD);

      // Update filenames in SFvol (now without noise volumes!)
      prm.remake_SFvol(iter,false,false);
      if (ML2D_prm.fourier_mode)
	ML2D_prm.calculate_wiener_defocus_series(spectral_signal,iter);

      if (!converged && iter+1<=prm.Niter) {
	// All nodes again: project and read new references from disc
	prm.project_reference_volume(ML2D_prm.SFr,rank);
	MPI_Barrier(MPI_COMM_WORLD);
	ML2D_prm.SFr.go_beginning();
	c=0;
	while(!ML2D_prm.SFr.eof()) {
	  ML2D_prm.Iref[c].read(ML2D_prm.SFr.NextImg(),false,false,false,false);
	  ML2D_prm.Iref[c]().set_Xmipp_origin();
	  c++;
	}
      }

      iter++;
    } // end loop iterations

    if (!converged && prm.verb>0)
      cerr <<"--> Optimization was stopped before convergence was reached!"<<endl;

  } catch (Xmipp_error XE) {
    if (rank==0) {
      cout << XE;
      if (prm.fourier_mode) prm.MLF_usage();
      else prm.usage();
    }
    MPI_Finalize(); exit(1);
  }

  MPI_Finalize();	
  return 0;

}




