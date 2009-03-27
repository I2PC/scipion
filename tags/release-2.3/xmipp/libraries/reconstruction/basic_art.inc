/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include "recons_misc.h"
#include "symmetrize.h"

#include <data/docfile.h>
#include <data/mask.h>
#include <data/wavelet.h>
#include <data/filters.h>
#include <sys/time.h>
#include <algorithm>

class POCSClass;
class VariabilityClass;
void update_residual_vector(Basic_ART_Parameters &prm, GridVolume &vol_basis,
   double &kappa, double &pow_residual_vol, double &pow_residual_imgs);

/* ------------------------------------------------------------------------- */
/* Init history                                                              */
/* ------------------------------------------------------------------------- */
template <class Extra_ART_Parameters>
void Basic_ART_Init_history(Basic_ART_Parameters &prm,
   const Extra_ART_Parameters &eprm, const GridVolume &vol_basis0) {

// Show general information ................................................
   *prm.fh_hist << " RECONSTRUCTION HISTORY: \n\n";
   *prm.fh_hist << " Parameters -------------------------------------------------\n";
   *prm.fh_hist << " Projections file: " << prm.fn_sel << std::endl;
   *prm.fh_hist << " CTF file:         " << prm.fn_ctf << std::endl;
   *prm.fh_hist << " Unmatched projectors:" << prm.unmatched << std::endl;
   *prm.fh_hist << " Denoise: " << prm.denoise << std::endl;
   *prm.fh_hist << " Sampling= " << prm.sampling << std::endl;
   *prm.fh_hist << prm.basis << std::endl;
   switch (prm.grid_type) {
      case FCC: *prm.fh_hist << " Grid Type: FCC "; break;
      case BCC: *prm.fh_hist << " Grid Type: BCC "; break;
      case CC:  *prm.fh_hist << " Grid Type: CC ";  break;
   }
   *prm.fh_hist << " Ray length: " << prm.ray_length << std::endl;
   *prm.fh_hist << "\n Radius of the interest sphere= " << prm.R
               << " pixels" << std::endl;
   *prm.fh_hist << " Grid unit=" << prm.grid_relative_size
               << " pixels" << std::endl;
   if (prm.proj_ext==0) *prm.fh_hist << " No Projection extension\n";
   else *prm.fh_hist << " Projection extension frame: " << prm.proj_ext
      << " pixels\n";
   if (prm.Zoutput_volume_size==0) *prm.fh_hist << " Output volume size as input images\n";
   else *prm.fh_hist << " Output volume size (ZxYxX): "
      << prm.Zoutput_volume_size << "x" << prm.Youtput_volume_size
      << "x" << prm.Xoutput_volume_size << std::endl;
   *prm.fh_hist << " Iterations:    lambda=" << prm.lambda_list.transpose() << std::endl
               << "                No. Iter=" << prm.no_it << std::endl;
   if (prm.WLS) {
     *prm.fh_hist << " Perform weighted least-squares ART "<< std::endl;
     *prm.fh_hist << " Iterations:    kappa=" << prm.kappa_list.transpose() << std::endl
                 << "                No. Iter=" << prm.no_it << std::endl;
   }
   *prm.fh_hist << " Parallel mode: ";
   switch (prm.parallel_mode) {
      case Basic_ART_Parameters::ART:   *prm.fh_hist << "ART\n"; break;
      case Basic_ART_Parameters::SIRT:  *prm.fh_hist << "SIRT\n"; break;
      case Basic_ART_Parameters::pSIRT:  *prm.fh_hist << "Parallel SIRT\n"; break;
      case Basic_ART_Parameters::pfSIRT:  *prm.fh_hist << "Parallel False SIRT\n"; break;
      case Basic_ART_Parameters::pSART:
          *prm.fh_hist << "SART, block size=" << prm.block_size << std::endl;
          break;
      case Basic_ART_Parameters::pAVSP:   *prm.fh_hist << "AVSP\n"; break;
      case Basic_ART_Parameters::pBiCAV:
          *prm.fh_hist << "BiCAV, block size=" << prm.block_size << std::endl;
          break;
      case Basic_ART_Parameters::pCAV:  *prm.fh_hist << "CAV (global algorithm)\n";break;
   }
   *prm.fh_hist << " Equation mode: ";
   if (prm.eq_mode==CAV) *prm.fh_hist << "CAV (Projection update)\n";
   else if (prm.eq_mode==CAVK) *prm.fh_hist << "CAVK\n";
   else if (prm.eq_mode==CAVARTK) *prm.fh_hist << "CAVARTK\n";
   else *prm.fh_hist << "ARTK\n";
   *prm.fh_hist << " Projections: Ydim x Xdim = " << prm.projYdim
      << " x " << prm.projXdim << std::endl;
   *prm.fh_hist << " Surface mask: " << prm.fn_surface_mask << std::endl;
   *prm.fh_hist << " POCS freq: " << prm.POCS_freq << std::endl;
   *prm.fh_hist << " Known volume: " << prm.known_volume << std::endl;
   *prm.fh_hist << " Positivity: " <<prm.positivity << std::endl;
   *prm.fh_hist << " Apply shifts in images headers: " << prm.apply_shifts << std::endl;
   *prm.fh_hist << " Symmetry file: " << prm.fn_sym << std::endl;
   *prm.fh_hist << " Force Symmetry each: " << prm.sym_each << " projections"
                <<std::endl;
   *prm.fh_hist << " Maximun absolute tilt angle: " << prm.max_tilt << " degrees"
                <<std::endl;
   *prm.fh_hist << " Generating symmetry group: " << !prm.do_not_generate_subgroup << std::endl;
   *prm.fh_hist << " Do not use symmetrized projections: "
       << !prm.do_not_use_symproj << std::endl;
   *prm.fh_hist << " Number of total projections (including symmetrized): "
      << prm.numIMG << std::endl;
   *prm.fh_hist << " Number of different projections: " << prm.trueIMG << std::endl;
   *prm.fh_hist << " Forstd::cing symmetry: " << prm.force_sym << std::endl;
   *prm.fh_hist << " Stop at: " << prm.stop_at << std::endl;
   if (prm.random_sort)
       *prm.fh_hist << "Random sort" << std::endl;
   else
       *prm.fh_hist << "Sort with last " << prm.sort_last_N << " images\n";
   *prm.fh_hist << " Variability analysis: " << prm.variability_analysis << std::endl;
   *prm.fh_hist << " Sparsity epsilon: " << prm.sparseEps << std::endl;
   *prm.fh_hist << " Diffusion weight: " << prm.diffusionWeight << std::endl;
   if (prm.SL.SymsNo()!=0) {
      Matrix2D<double> L(4,4),R(4,4); // A matrix from the list
      *prm.fh_hist << "Symmetry matrices -------\n";
      for (int j=0; j<prm.SL.SymsNo(); j++) {
         prm.SL.get_matrices(j,L,R);
         *prm.fh_hist << "Left  Symmetry matrix " << j << std::endl << L;
         *prm.fh_hist << "Right Symmetry matrix " << j << std::endl << R << std::endl;
      }
   }
   *prm.fh_hist << " Saving intermidiate at every "
      << prm.save_intermidiate_every << " projections\n";
   if(prm.ref_trans_step > 0.)
     {
     *prm.fh_hist << " Refine translational alignement after "
                 << prm.ref_trans_after << " projection presentations\n"
		 << " Maximun allowed shift is: " << prm.ref_trans_step << "\n\n";
     }

// Show angles .............................................................
   // Prepare info structure for showing
   DocFile DF;
   Matrix1D<double> v(4);
   for (int i=0; i<prm.numIMG; i++) {
      v(0)=prm.IMG_Inf[i].rot;
      v(1)=prm.IMG_Inf[i].tilt;
      v(2)=prm.IMG_Inf[i].psi;
      v(3)=prm.IMG_Inf[i].sym;
      DF.append_data_line(v);
   }

   // Now show
   *prm.fh_hist << " Projection angles -----------------------------------------\n";
   *prm.fh_hist << DF;
   *prm.fh_hist << " -----------------------------------------------------------\n";

// Show Initial volume and volume structure ................................
   if (prm.fn_start != "")
      *prm.fh_hist << "Starting from file: " << prm.fn_start << std::endl;
   else
      *prm.fh_hist << "Starting from a zero volume\n";

   *prm.fh_hist << "Grid structure ------\n" << vol_basis0.grid();

// Show extra information ..................................................
   *prm.fh_hist << "Extra information ------------------------------------------\n";
   *prm.fh_hist << eprm;
}

/* ------------------------------------------------------------------------- */
/* ART Several Iterations                                                    */
/* ------------------------------------------------------------------------- */
//#define DEBUG_POCS
template <class Extra_ART_Parameters>
void Basic_ART_iterations(Basic_ART_Parameters &prm,
   Extra_ART_Parameters &eprm, GridVolume &vol_basis, int rank) {
   // Some variables .......................................................
   int        ART_numIMG;              // Normalizing factor
   GridVolume *ptr_vol_out;            // Pointer to output volume
   GridVolume vol_basis_out;           // Output volume (only useful in SIRT)
   VolumeXmipp vol_voxels;             // This one is useful only in the
                                       // case of saving intermidiate volumes
   int Xoutput_volume_size, Youtput_volume_size, Zoutput_volume_size;
   double aux_tilt;

   // Projection related ...................................................
   Projection      read_proj;          // Projection read from file
   Projection      theo_proj;          // Projection from the
                                       // reconstruction
   Projection      alig_proj;          // Projection with the correlation
                                       // maps needed for translation aligment				
   Projection      corr_proj;          // Image with the correction
                                       // factors for unitary basis
   Projection      diff_proj;          // Difference between the
                                       // theoretical and real image
									   
   // As this is a threaded implementation, create structures for threads, and
   // create threads
   pthread_t * th_ids;
   
   if( prm.threads > 1 )
   {
		th_ids = (pthread_t *)malloc( prm.threads * sizeof( pthread_t));
	
		// Initialize the structures which will contain the parameters passed to different
		// threads
		project_threads= (project_thread_params *) malloc ( prm.threads * sizeof( project_thread_params ) );
	
		// Initialize barrier to wait for working threads and the master thread.
		barrier_init( &project_barrier, (prm.threads+1) );

		// Threas are created in a waiting state. They can only run when master thread unlocks them by calling
		// barrier_wait() 
		
		for( int c = 0 ; c < prm.threads ; c++ )
		{
			project_threads[c].thread_id = c;
			project_threads[c].threads_count = prm.threads;
			project_threads[c].destroy = false;
		
			pthread_create( (th_ids+c), NULL, project_SimpleGridThread<double>, (void *)(project_threads+c) );
		}
	}
   
   // Reconstruction results ...............................................
   double          mean_error,mean_error_1stblock;
   double          global_mean_error,global_mean_error_1stblock;

   // Initialize residual image vector for wlsART ..........................
   if (prm.WLS) {
     prm.residual_imgs.clear();
     for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++) {
       read_proj.read(prm.IMG_Inf[iact_proj].fn_proj);
       read_proj().setXmippOrigin();
       read_proj().initZeros();
       prm.residual_imgs.push_back(read_proj);
     }
   }

   // Some initialisation ..................................................
   Xoutput_volume_size=(prm.Xoutput_volume_size==0) ?
      prm.projXdim:prm.Xoutput_volume_size;
   Youtput_volume_size=(prm.Youtput_volume_size==0) ?
      prm.projYdim:prm.Youtput_volume_size;
   Zoutput_volume_size=(prm.Zoutput_volume_size==0) ?
      prm.projXdim:prm.Zoutput_volume_size;

   // POCS constraints .....................................................
   POCSClass POCS(&prm, Zoutput_volume_size, Youtput_volume_size,
      Xoutput_volume_size);

   // Variability analysis .................................................
   VariabilityClass VC(&prm, Zoutput_volume_size, Youtput_volume_size,
      Xoutput_volume_size);

   // Noisy reconstruction .................................................
   GridVolume vol_basis_noisy; // Output volume (only for ART)
   Projection noisy_projection;
   SelFile SF_noise, SF_signal;
   if (prm.noisy_reconstruction) {
     vol_basis_noisy=vol_basis;
     vol_basis_noisy.initZeros();
   }

   // If SIRT set normalizing factor and create output volume ..............
   if (prm.parallel_mode==Basic_ART_Parameters::SIRT ) {
  	ART_numIMG=prm.numIMG;
      	vol_basis_out=vol_basis;         // Copy the structure of vol_basis
	ptr_vol_out=&vol_basis_out;      // Pointer to output volume
	if (prm.variability_analysis) ptr_vol_out->initZeros();
   } else if ( prm.parallel_mode==Basic_ART_Parameters::pfSIRT ||
       prm.parallel_mode==Basic_ART_Parameters::pSIRT ||
       prm.parallel_mode==Basic_ART_Parameters::pSART ||
       prm.parallel_mode==Basic_ART_Parameters::pCAV )
   {
   // In those cases, normalization must be done at the top level program once we
   // have the reconstruction values. Have a look ar /Applications/Src/MPIArt/mpi_art.cc for
   // an example
      ART_numIMG= 1; 		
      ptr_vol_out=&vol_basis_out;      // Pointer to output volume
      vol_basis_out=vol_basis;         // Copy the structure of vol_basis
                                       // and possible initial values
   } else if (prm.eq_mode==CAV) {
      ART_numIMG=1;                    // Normalizing factor = Total no. images
      ptr_vol_out=&vol_basis_out;      // Pointer to output volume
      vol_basis_out=vol_basis;         // Copy the structure of vol_basis
                                       // and possible initial values
   } else {
      ART_numIMG=1;                    // No normalizing factor
      ptr_vol_out=&vol_basis;          // Output volume is the same as
                                       // input one
   }
   // Now iterate ..........................................................
   TimeStamp time0;                    // For measuring the elapsed time
   annotate_time(&time0);
   int images=0;
   double mean_error_2ndblock,pow_residual_imgs;
   bool iv_launched=false;
   for (int it = 0; it < prm.no_it; it++) {
      // Initialisation of some variables
      global_mean_error = 0;
      global_mean_error_1stblock = 0;
      POCS.newIteration();
      VC.newIteration();
      if (rank==-1) {
	 std::cerr << "Running iteration " << it << " with lambda= " << prm.lambda(it)<< "\n" << std::endl;
	 if (!(prm.tell&TELL_SHOW_ERROR)) init_progress_bar(prm.numIMG);
      }

      // For each projection -----------------------------------------------
      for (int act_proj = 0; act_proj < prm.numIMG ; act_proj++) {
          POCS.newProjection();

          // Select next projection ........................................
          int iact_proj; // Index inside the sorting information for act_proj
          if (prm.tell&TELL_MANUAL_ORDER) {
             int proj_number;
             std::cerr << "Introduce next projection to study: ";
             std::cin  >> proj_number;
             int sym_number=-1;
             if (prm.SL.SymsNo()!=0) {
                std::cerr << "Introduce symmetry to study: ";
                std::cin  >> sym_number;
             }
             iact_proj=0;
             while (iact_proj<prm.numIMG) {
                if (prm.IMG_Inf[iact_proj].fn_proj.get_number()==proj_number &&
                    prm.IMG_Inf[iact_proj].sym==sym_number) break;
                iact_proj++;
             }
          } else {
             iact_proj = prm.ordered_list(act_proj);
          }
	  read_proj.read(prm.IMG_Inf[iact_proj].fn_proj,prm.apply_shifts);
          read_proj.moveOriginTo_center();
          read_proj.set_rot (prm.IMG_Inf[iact_proj].rot);
          read_proj.set_tilt(prm.IMG_Inf[iact_proj].tilt);
          read_proj.set_psi (prm.IMG_Inf[iact_proj].psi);

      	  // If noisy reconstruction
	  if (prm.noisy_reconstruction) {
	     init_random_generator(prm.IMG_Inf[iact_proj].seed);
	     noisy_projection().resize(read_proj());
	     noisy_projection().initRandom(0,1,"gaussian");
             noisy_projection.moveOriginTo_center();
             noisy_projection.set_rot (prm.IMG_Inf[iact_proj].rot);
             noisy_projection.set_tilt(prm.IMG_Inf[iact_proj].tilt);
             noisy_projection.set_psi (prm.IMG_Inf[iact_proj].psi);
	     if (prm.IMG_Inf[iact_proj].sym==-1 && it==0) {
	        FileName fn_noise=prm.fn_root+"_noise_proj"+
		   integerToString(read_proj.name().get_number(),5)+".xmp";
	        noisy_projection.write(fn_noise);
	        SF_noise.insert(fn_noise);
		SF_signal.insert(read_proj.name());
	     }
	  }

          //skipping if  tilt greater than max_tilt
          //tilt is inbetween 0 and 360
          aux_tilt=read_proj.tilt();
          if((aux_tilt > prm.max_tilt && aux_tilt < 180.-prm.max_tilt) ||
             (aux_tilt > prm.max_tilt + 180 && aux_tilt < 360.-prm.max_tilt)) {
              std::cerr << "Skipping Proj no: " << iact_proj
		             << "tilt=" << read_proj.tilt()  << std::endl;
		        continue;
          }

          // Projection extension? .........................................
          if (prm.proj_ext!=0) {
             read_proj().window(
                STARTINGY (read_proj())-prm.proj_ext,
                STARTINGX (read_proj())-prm.proj_ext,
                FINISHINGY(read_proj())+prm.proj_ext,
                FINISHINGX(read_proj())+prm.proj_ext);
	     noisy_projection().resize(read_proj());
	  }

         //Skip if desired	
          if (prm.tell&TELL_ONLY_SYM)
	     if( prm.IMG_Inf[iact_proj].sym != -1) {
		std::cerr << "Skipping Proj no: " << iact_proj
		     << " with symmetry no: " << prm.IMG_Inf[iact_proj].sym
		     <<  std::endl;
		continue;
	     } else
		std::cerr << "NO Skipping Proj no: " << iact_proj
		     << " with symmetry no: " << prm.IMG_Inf[iact_proj].sym
		     <<  std::endl;
		
          // For wlsART: use alig_proj for residual image!!
	  if (prm.WLS) alig_proj=prm.residual_imgs[iact_proj];

          // Apply the reconstruction algorithm ............................
      	  // Notice that the following function is extern
	  ART_single_step(vol_basis, ptr_vol_out,
      	      prm, eprm,
              theo_proj, read_proj, prm.IMG_Inf[iact_proj].sym , diff_proj,
	      corr_proj, alig_proj,
              mean_error, ART_numIMG, prm.lambda(it),
	      images, prm.IMG_Inf[iact_proj].fn_ctf);

          if (prm.WLS) {
            prm.residual_imgs[iact_proj]=alig_proj;
	    global_mean_error_1stblock += diff_proj().sum2()/(XSIZE(diff_proj())*YSIZE(diff_proj()));
          }

          global_mean_error += mean_error;
	  if (prm.noisy_reconstruction) {
	     double noise_mean_error;
	     ART_single_step(vol_basis_noisy, &vol_basis_noisy,
      		 prm, eprm,
        	 theo_proj, noisy_projection, prm.IMG_Inf[iact_proj].sym,
		 diff_proj,  corr_proj, alig_proj,
        	 noise_mean_error, ART_numIMG, prm.lambda(it),
		 images, prm.IMG_Inf[iact_proj].fn_ctf);
	  }

          // Force symmetry in the volume ..................................
          // so far only crystallographic
          if (prm.sym_each &&
             act_proj%prm.sym_each==0 &&
             (act_proj!=0 || prm.sym_each==1)) {
             apply_symmetry(vol_basis,ptr_vol_out,eprm, prm.grid_type);
	     if (prm.noisy_reconstruction)
		apply_symmetry(vol_basis_noisy,&vol_basis_noisy,eprm, prm.grid_type);
	  }

          // Apply POCS ....................................................
	  POCS.apply(vol_basis,it,images);
	  if (prm.noisy_reconstruction)
	     POCS.apply(vol_basis_noisy,it,images);

          // Variability analysis ..........................................
	  if (prm.variability_analysis)
	     VC.newUpdateVolume(ptr_vol_out,read_proj);

          // Apply anisotropic diffusion ...................................
          if (it>=1 && prm.diffusionWeight>-1)
          {
             prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
             Matrix1D<double> alpha;
             alpha=vectorR3(1.0,1.0,prm.diffusionWeight);
             double regError=tomographicDiffusion(vol_voxels(),alpha,
                prm.lambda(it)/100);
             if (prm.tell&TELL_SHOW_ERROR)
                std::cout << "Regularization error = " << regError << std::endl;
             *prm.fh_hist << "Regularization error = " << regError << std::endl;
             prm.basis.changeFromVoxels(vol_voxels(),vol_basis,prm.grid_type,
                prm.grid_relative_size, NULL, NULL, prm.R, prm.threads);
          }

          // Apply sparsity constraint .....................................
          if (it>=1 && prm.sparseEps>0)
          {
             prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
             VolumeXmipp vol_wavelets, vol_wavelets_abs;
             set_DWT_type(DAUB12);
             DWT(vol_voxels(),vol_wavelets());
             vol_wavelets_abs()=vol_wavelets();
             vol_wavelets_abs().selfABSnD();
             double *begin=MULTIDIM_ARRAY(vol_wavelets_abs());
             double *end=MULTIDIM_ARRAY(vol_wavelets_abs())+
                 MULTIDIM_SIZE(vol_wavelets_abs());
             std::sort(begin,end);
             double threshold1=DIRECT_MULTIDIM_ELEM(vol_wavelets_abs(),
                (long int)((1-prm.sparseEps)*MULTIDIM_SIZE(vol_wavelets_abs())));
             std::cout << "Threshold=" << threshold1 << std::endl;
             vol_wavelets().threshold("abs_below", threshold1, 0.0);
             IDWT(vol_wavelets(),vol_voxels());
             prm.basis.changeFromVoxels(vol_voxels(),vol_basis,prm.grid_type,
                prm.grid_relative_size, NULL, NULL, prm.R, prm.threads);
          }

          // Show results ..................................................
          *prm.fh_hist << prm.IMG_Inf[iact_proj].fn_proj << ", sym="
             << prm.IMG_Inf[iact_proj].sym << "\t\t" << mean_error;
          if (POCS.apply_POCS)
             *prm.fh_hist << "\tPOCS:" << POCS.POCS_mean_error;
          *prm.fh_hist << std::endl;
          if (prm.tell&TELL_SHOW_ERROR) {
             std::cout << prm.IMG_Inf[iact_proj].fn_proj << ", sym="
                << prm.IMG_Inf[iact_proj].sym << "\t\t" << mean_error;
             if (POCS.apply_POCS)
                std::cout        << "\tPOCS:" << POCS.POCS_mean_error;
             std::cout << std::endl;
          } else if (act_proj%XMIPP_MAX(1,prm.numIMG/60)==0) progress_bar(act_proj);

          if (prm.tell&TELL_STATS) {
             std::cout << "   read      ";   read_proj().printStats();
             std::cout << "\n   theo      "; theo_proj().printStats();
             std::cout << "\n   corr      "; corr_proj().printStats();
             std::cout << "\n   alig      "; alig_proj().printStats();
             std::cout << "\n   diff      "; diff_proj().printStats();
             std::cout << "\n   subvol(0) "; (*ptr_vol_out)(0)().printStats();
             std::cout << "\n";
          }

          if (prm.tell&TELL_SAVE_AT_EACH_STEP) {
      	     std::cout << "Stats PPPdiff.xmp: "; diff_proj().printStats(); std::cout << std::endl;
      	     std::cout << "Stats PPPtheo.xmp: "; theo_proj().printStats(); std::cout << std::endl;
      	     std::cout << "Stats PPPread.xmp: "; read_proj().printStats(); std::cout << std::endl;
      	     std::cout << "Stats PPPcorr.xmp: "; corr_proj().printStats(); std::cout << std::endl;
             diff_proj.write("PPPdiff.xmp");
             theo_proj.write("PPPtheo.xmp");
             read_proj.write("PPPread.xmp");
             corr_proj.write("PPPcorr.xmp");
	     if(act_proj!=0)
                alig_proj.write("PPPalign.xmp");
             vol_basis.write("PPPbasis.basis");
             prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
      	     std::cout << "Stats PPPvol.vol : "; vol_voxels().printStats(); std::cout << std::endl;
             vol_voxels.write("PPPvol.vol");

	     if (!iv_launched) {
	        system("xmipp_show -img PPPdiff.xmp PPPtheo.xmp PPPread.xmp PPPcorr.xmp -dont_apply_geo -poll &");
		system("xmipp_show -vol PPPvol.vol -poll &");
		iv_launched=true;
	     }
             std::cout << "\nHit any key and enter\n"; char c; std::cin >> c;
         }

      	 // Save intermidiate
	 if ((prm.tell&TELL_SAVE_INTERMIDIATE | prm.tell&TELL_IV) &&
             prm.save_intermidiate_every!=0 &&
	     act_proj%prm.save_intermidiate_every==0) {
	    if (prm.tell&TELL_SAVE_INTERMIDIATE)
               std::cerr << "\nSaving intermidiate ...\n"
                    << "Converting basis volume to voxels ...\n";
            // Save reconstructed volume
            prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
               Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
      	    if (prm.tell&TELL_SAVE_INTERMIDIATE)
               vol_voxels.write(prm.fn_root+"it"+integerToString(it)+"proj"+
		  integerToString(act_proj,5)+".vol");
	    else
               vol_voxels.write("PPPvol.vol");

            // Launch viewer
	    if (!iv_launched) {
	       system("xmipp_show -vol PPPvol.vol -poll &");
	       iv_launched=true;
	    }
	 }

         // Check if algorithm must stop via stop_at
         if (++images==prm.stop_at) break;
      }

      if (!(prm.tell&TELL_SHOW_ERROR)) progress_bar(prm.numIMG);

      // Update residual images for WLS
      if (prm.WLS) {
        double kappa=prm.kappa(it);
        update_residual_vector( prm, vol_basis, kappa,
           mean_error_2ndblock, pow_residual_imgs);
      }

      // Prepare for next global iteration ---------------------------------
      // Calculate norm and print
      if (rank==-1) {
	 *prm.fh_hist << "Finished - Iteration " << it << std::endl;
         if (prm.WLS)
	    std::cerr        << "   Weighted error: " << global_mean_error
	                << " 1st block: "        << global_mean_error_1stblock
                        << " 2nd block: "        << mean_error_2ndblock
	                << " residual: "         << pow_residual_imgs << std::endl;
         else
            std::cerr << "   Global mean squared error: "
                 <<  global_mean_error/prm.numIMG << std::endl;

         if (prm.WLS)
	    *prm.fh_hist << "   Weighted error: " << global_mean_error
	                << " 1st block: "        << global_mean_error_1stblock
                        << " 2nd block: "        << mean_error_2ndblock
	                << " residual: "         << pow_residual_imgs << std::endl;
         else
	    *prm.fh_hist << "   Global mean squared error: "
                        << global_mean_error/prm.numIMG << std::endl;

	 if (POCS.apply_POCS) {
            std::cerr        << "   POCS Global mean squared error: "
                	<< POCS.POCS_global_mean_error/POCS.POCS_N << std::endl;
            *prm.fh_hist << "   POCS Global mean squared error: "
                	<< POCS.POCS_global_mean_error/POCS.POCS_N << std::endl;
	 }
      }

      // Convert volume and write if not last iteration
      // If in SIRT mode move the volume to the reference one
      if ((prm.parallel_mode==Basic_ART_Parameters::SIRT
             && !prm.variability_analysis)||
          prm.parallel_mode==Basic_ART_Parameters::pSIRT ||
          prm.parallel_mode==Basic_ART_Parameters::pfSIRT ||
          prm.parallel_mode==Basic_ART_Parameters::pSART ||
	  prm.parallel_mode==Basic_ART_Parameters::pCAV ||
	  prm.eq_mode==CAV )
          vol_basis=vol_basis_out;

      if (prm.tell&TELL_SAVE_INTERMIDIATE && it!=prm.no_it-1) {
         if (rank==-1) std::cerr << "Converting basis volume to voxels ...\n";
         prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
            Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
         vol_voxels.write(prm.fn_root+"it"+integerToString(it)+".vol");

         if (prm.tell&TELL_SAVE_BASIS)
            vol_basis.write(prm.fn_root+"it"+integerToString(it)+".basis");
      }
      
      // Check if algorithm must stop via stop_at
      if (images==prm.stop_at) break;
   }

   // Times on screen
   if (rank==-1) {
      std::cout << "\nTime of " << prm.no_it << " iterations: \n";
      print_elapsed_time(time0);
   }

   // Finish variability analysis
   VC.finishAnalysis();

   // Save the noisy reconstruction
   if (prm.noisy_reconstruction) {
      VolumeXmipp vol_voxels_noisy;
      prm.basis.changeToVoxels(vol_basis_noisy, &(vol_voxels_noisy()),
         Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
      vol_voxels_noisy.write(prm.fn_root+"_noise.vol");
      SF_noise.write(prm.fn_root+"_noise_proj.sel");
      SF_signal.write(prm.fn_root+"_signal_proj.sel");
   }
   
   // Destroy created threads. This is done in a tricky way. At this point, threads
   // are "slept" waiting for a barrier to be reached by the master thread to continue
   // projecting/backprojecting a new projection. Here we set the flag destroy=true so
   // thre threads won't process a projections but will return.
   if( prm.threads > 1 )
   {
		for( int c = 0 ; c < prm.threads ; c++ )
		{	
			project_threads[c].destroy = true;
		}	

		// Trigger threads "self-destruction"
		barrier_wait( &project_barrier );

		// Wait for effective threads death
		for( int c = 0 ; c < prm.threads ; c++ )
		{
			pthread_join(*(th_ids+c),NULL);
		}

		// Destroy barrier and mutex, as they are no longer needed.
                // Sjors, 27march2009
                // NO, dont do this, because running a second threaded-ART
                // in a single program will not have mutexes anymore...
		//pthread_mutex_destroy( &project_mutex );
		//barrier_destroy( &project_barrier );
	}
}

/* ------------------------------------------------------------------------- */
/* ROUTINE ART                                                               */
/* ------------------------------------------------------------------------- */
template <class Extra_ART_Parameters>
void Basic_ROUT_Art(Basic_ART_Parameters &prm,
   Extra_ART_Parameters &eprm, VolumeXmipp &vol_voxels,
   GridVolume &vol_basis) {
// Configure time clock
   time_config();

struct timeval start_time, end_time;
long int init_usecs, process_usecs, finish_usecs;

gettimeofday(&start_time, NULL);

// Produce side information and initial volume
   prm.produce_Side_Info(vol_basis);
//calculate symmetry in produce_Side_Info and exit if wrong sampling
//printting the right sampling
   eprm.produce_Side_Info(prm,vol_basis);

// Show parameters and initiate history
   Basic_ART_Init_history(prm, eprm, vol_basis);

gettimeofday(&end_time, NULL);

init_usecs = (end_time.tv_sec-start_time.tv_sec)*1000000+(end_time.tv_usec-start_time.tv_usec);
 
gettimeofday(&start_time,NULL);

// Iterations
   Basic_ART_iterations(prm, eprm, vol_basis);

gettimeofday(&end_time,NULL);
   
process_usecs = (end_time.tv_sec-start_time.tv_sec)*1000000+(end_time.tv_usec-start_time.tv_usec);

gettimeofday(&start_time,NULL);

// Finish iterations
   finish_ART_iterations(prm, eprm, vol_basis);

// Write final volume
   int Xoutput_volume_size=(prm.Xoutput_volume_size==0) ?
      prm.projXdim:prm.Xoutput_volume_size;
   int Youtput_volume_size=(prm.Youtput_volume_size==0) ?
      prm.projYdim:prm.Youtput_volume_size;
   int Zoutput_volume_size=(prm.Zoutput_volume_size==0) ?
      prm.projXdim:prm.Zoutput_volume_size;

//   int min_distance = ceil( ( 2 * prm.grid_relative_size ) / prm.basis.blob.radius) + 1;
   
   prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
      Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size,prm.threads);

   vol_voxels.write(prm.fn_root+".vol");
   if (prm.tell&TELL_SAVE_BASIS) vol_basis.write(prm.fn_root+".basis");
   prm.fh_hist->close();

gettimeofday(&end_time,NULL);

finish_usecs = (end_time.tv_sec-start_time.tv_sec)*1000000+(end_time.tv_usec-start_time.tv_usec);
      
std::cout << "INIT_TIME: " << (double)init_usecs/(double)1000000 << std::endl;
std::cout << "PROCESS_TIME: " << (double)process_usecs/(double)1000000 << std::endl;
std::cout << "FINISH_TIME: " << (double)finish_usecs/(double)1000000 << std::endl;

}
