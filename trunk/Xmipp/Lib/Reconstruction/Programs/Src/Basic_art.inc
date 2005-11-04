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

#include "../recons_misc.hh"
#include "../Prog_symmetrize.hh"
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippMasks.hh>

class POCSClass;
class VariabilityClass;

/* ------------------------------------------------------------------------- */
/* Init history                                                              */
/* ------------------------------------------------------------------------- */
template <class Extra_ART_Parameters>
void Basic_ART_Init_history(Basic_ART_Parameters &prm,
   const Extra_ART_Parameters &eprm, const GridVolume &vol_basis0) {
   
// Show general information ................................................
   prm.fh_hist << " RECONSTRUCTION HISTORY: \n\n";
   prm.fh_hist << " Parameters -------------------------------------------------\n";
   prm.fh_hist << " Projections file: " << prm.fn_sel << endl;
   prm.fh_hist << " CTF file:         " << prm.fn_ctf << endl;
   prm.fh_hist << " Unmatched projectors:" << prm.unmatched << endl;
   prm.fh_hist << " Sampling= " << prm.sampling << endl;
   prm.fh_hist << prm.basis << endl;
   switch (prm.grid_type) {
      case FCC: prm.fh_hist << " Grid Type: FCC "; break;
      case BCC: prm.fh_hist << " Grid Type: BCC "; break;
      case CC:  prm.fh_hist << " Grid Type: CC ";  break;
   }
   prm.fh_hist << " Ray length: " << prm.ray_length << endl;
   prm.fh_hist << "\n Radius of the interest sphere= " << prm.R 
               << " pixels" << endl;
   prm.fh_hist << " Grid unit=" << prm.grid_relative_size 
               << " pixels" << endl;
   if (prm.proj_ext==0) prm.fh_hist << " No Projection extension\n";
   else prm.fh_hist << " Projection extension frame: " << prm.proj_ext
      << " pixels\n";
   if (prm.Zoutput_volume_size==0) prm.fh_hist << " Output volume size as input images\n";
   else prm.fh_hist << " Output volume size (ZxYxX): "
      << prm.Zoutput_volume_size << "x" << prm.Youtput_volume_size
      << "x" << prm.Xoutput_volume_size << endl;
   prm.fh_hist << " Iterations:    lambda=" << prm.lambda_list.transpose() << endl
               << "                No. Iter=" << prm.no_it << endl;
   if (prm.WLS) {
     prm.fh_hist << " Perform weighted least-squares ART "<< endl;
     prm.fh_hist << " Iterations:    kappa=" << prm.kappa_list.transpose() << endl
                 << "                No. Iter=" << prm.no_it << endl;
   }
   prm.fh_hist << " Parallel mode: ";
   switch (prm.parallel_mode) {
      case Basic_ART_Parameters::ART:   prm.fh_hist << "ART\n"; break;
      case Basic_ART_Parameters::SIRT:  prm.fh_hist << "SIRT\n"; break;
      case Basic_ART_Parameters::pSIRT:  prm.fh_hist << "Parallel SIRT\n"; break;
      case Basic_ART_Parameters::pfSIRT:  prm.fh_hist << "Parallel False SIRT\n"; break;
      case Basic_ART_Parameters::pSART:
          prm.fh_hist << "SART, block size=" << prm.block_size << endl;
          break;
      case Basic_ART_Parameters::pAVSP:   prm.fh_hist << "AVSP\n"; break;
      case Basic_ART_Parameters::pBiCAV:
          prm.fh_hist << "BiCAV, block size=" << prm.block_size << endl;
          break;
      case Basic_ART_Parameters::pCAV:  prm.fh_hist << "CAV (global algorithm)\n";break;
   }
   prm.fh_hist << " Equation mode: ";
   if (prm.eq_mode==CAV) prm.fh_hist << "CAV (Projection update)\n";
   else if (prm.eq_mode==CAVK) prm.fh_hist << "CAVK\n";
   else if (prm.eq_mode==CAVARTK) prm.fh_hist << "CAVARTK\n";
   else prm.fh_hist << "ARTK\n";
   prm.fh_hist << " Projections: Ydim x Xdim = " << prm.projYdim
      << " x " << prm.projXdim << endl;
   prm.fh_hist << " Surface mask: " << prm.fn_surface_mask << endl;
   prm.fh_hist << " POCS freq: " << prm.POCS_freq << endl;
   prm.fh_hist << " Known volume: " << prm.known_volume << endl;
   prm.fh_hist << " Positivity: "; print(prm.fh_hist,prm.positivity);
   prm.fh_hist << " Apply shifts in images headers: "; print(prm.fh_hist,prm.apply_shifts);
      prm.fh_hist << endl;
   prm.fh_hist << " Symmetry file: " << prm.fn_sym << endl;
   prm.fh_hist << " Force Symmetry each: " << prm.sym_each << " projections" 
                                                           <<endl;
   prm.fh_hist << " Maximun absolute tilt angle: " << prm.max_tilt << " degrees" 
                                                           <<endl;
   prm.fh_hist << " Generating symmetry group: ";
       print(prm.fh_hist,!prm.do_not_generate_subgroup); prm.fh_hist << endl;
   prm.fh_hist << " Do not use symmetrized projections: ";
       print(prm.fh_hist,!prm.do_not_use_symproj); prm.fh_hist << endl;
   prm.fh_hist << " Number of total projections (including symmetrized): "
      << prm.numIMG << endl;
   prm.fh_hist << " Number of different projections: " << prm.trueIMG << endl;
   prm.fh_hist << " Forcing symmetry: " << prm.force_sym << endl;
   prm.fh_hist << " Stop at: " << prm.stop_at << endl;
   if (prm.random_sort)
       prm.fh_hist << "Random sort" << endl;
   else
       prm.fh_hist << "Sort with last " << prm.sort_last_N << " images\n";
   prm.fh_hist << " Variability analysis: " << prm.variability_analysis << endl;
   if (prm.SL.SymsNo()!=0) {
      matrix2D<double> L(4,4),R(4,4); // A matrix from the list
      prm.fh_hist << "Symmetry matrices -------\n";
      for (int j=0; j<prm.SL.SymsNo(); j++) {
         prm.SL.get_matrices(j,L,R);
         prm.fh_hist << "Left  Symmetry matrix " << j << endl << L;
         prm.fh_hist << "Right Symmetry matrix " << j << endl << R << endl;
      }
   }
   prm.fh_hist << " Saving intermidiate at every "
      << prm.save_intermidiate_every << " projections\n";
   if(prm.ref_trans_step > 0.)
     {
     prm.fh_hist << " Refine translational alignement after "
                 << prm.ref_trans_after << " projection presentations\n"
		 << " Maximun allowed shift is: " << prm.ref_trans_step << "\n\n";
     }

// Show angles .............................................................
   // Prepare info structure for showing
   DocFile DF;
   matrix1D<double> v(4);
   for (int i=0; i<prm.numIMG; i++) {
      v(0)=prm.IMG_Inf[i].rot;
      v(1)=prm.IMG_Inf[i].tilt;
      v(2)=prm.IMG_Inf[i].psi;
      v(3)=prm.IMG_Inf[i].sym;
      DF.append_data_line(v);
   }

   // Now show
   prm.fh_hist << " Projection angles -----------------------------------------\n";
   prm.fh_hist << DF;
   prm.fh_hist << " -----------------------------------------------------------\n";

// Show Initial volume and volume structure ................................
   if (prm.fn_start != "")
      prm.fh_hist << "Starting from file: " << prm.fn_start << endl;
   else
      prm.fh_hist << "Starting from a zero volume\n";

   prm.fh_hist << "Grid structure ------\n" << vol_basis0.grid();

// Show extra information ..................................................
   prm.fh_hist << "Extra information ------------------------------------------\n";
   prm.fh_hist << eprm;
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

   // Reconstruction results ...............................................
   double          mean_error,mean_error_1stblock;
   double          global_mean_error,global_mean_error_1stblock;

   // Initialize residual image vector for wlsART ..........................
   if (prm.WLS) {
     prm.residual_imgs.clear();
     for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++) {
       read_proj.read(prm.IMG_Inf[iact_proj].fn_proj);
       read_proj().set_Xmipp_origin();
       read_proj().init_zeros();
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
     vol_basis_noisy.init_zeros();
   }

   // If SIRT set normalizing factor and create output volume ..............
   if (prm.parallel_mode==Basic_ART_Parameters::SIRT ) {
  	ART_numIMG=prm.numIMG;
      	vol_basis_out=vol_basis;         // Copy the structure of vol_basis 
	ptr_vol_out=&vol_basis_out;      // Pointer to output volume
	if (prm.variability_analysis) ptr_vol_out->init_zeros();
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
	 cerr << "Running iteration " << it << " with lambda= " << prm.lambda(it)<< "\n" << endl;
	 if (!(prm.tell&TELL_SHOW_ERROR)) init_progress_bar(prm.numIMG);
      }
      
      // For each projection -----------------------------------------------
      for (int act_proj = 0; act_proj < prm.numIMG ; act_proj++) {
          POCS.newProjection();
  
          // Select next projection ........................................
          int iact_proj; // Index inside the sorting information for act_proj
          if (prm.tell&TELL_MANUAL_ORDER) {
             int proj_number;
             cerr << "Introduce next projection to study: ";
             cin  >> proj_number;
             int sym_number=-1;
             if (prm.SL.SymsNo()!=0) {
                cerr << "Introduce symmetry to study: ";
                cin  >> sym_number;
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
          read_proj.move_origin_to_center();
          read_proj.rot ()=prm.IMG_Inf[iact_proj].rot;
          read_proj.tilt()=prm.IMG_Inf[iact_proj].tilt;
          read_proj.psi ()=prm.IMG_Inf[iact_proj].psi;
          
      	  // If noisy reconstruction
	  if (prm.noisy_reconstruction) {
	     init_random_generator(prm.IMG_Inf[iact_proj].seed);
	     noisy_projection().resize(read_proj());
	     noisy_projection().init_random(0,1,"gaussian");
             noisy_projection.move_origin_to_center();
             noisy_projection.rot ()=prm.IMG_Inf[iact_proj].rot;
             noisy_projection.tilt()=prm.IMG_Inf[iact_proj].tilt;
             noisy_projection.psi ()=prm.IMG_Inf[iact_proj].psi;
	     if (prm.IMG_Inf[iact_proj].sym==-1 && it==0) {
	        FileName fn_noise=prm.fn_root+"_noise_proj"+
		   ItoA(read_proj.name().get_number(),5)+".xmp";
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
              cerr << "Skipping Proj no: " << iact_proj 
		             << "tilt=" << read_proj.tilt()  << endl;
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
		cerr << "Skipping Proj no: " << iact_proj 
		     << " with symmetry no: " << prm.IMG_Inf[iact_proj].sym 
		     <<  endl;
		continue;
	     } else
		cerr << "NO Skipping Proj no: " << iact_proj 
		     << " with symmetry no: " << prm.IMG_Inf[iact_proj].sym 
		     <<  endl;
		
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

          // Show results ..................................................
          prm.fh_hist << prm.IMG_Inf[iact_proj].fn_proj << ", sym="
             << prm.IMG_Inf[iact_proj].sym << "\t\t" << mean_error;
          if (POCS.apply_POCS)
             prm.fh_hist << "\tPOCS:" << POCS.POCS_mean_error;
          prm.fh_hist << endl;
          if (prm.tell&TELL_SHOW_ERROR) {
             cout << prm.IMG_Inf[iact_proj].fn_proj << ", sym="
                << prm.IMG_Inf[iact_proj].sym << "\t\t" << mean_error;
             if (POCS.apply_POCS)
                cout        << "\tPOCS:" << POCS.POCS_mean_error;
             cout << endl;
          } else if (act_proj%MAX(1,prm.numIMG/60)==0) progress_bar(act_proj);

          if (prm.tell&TELL_STATS) {
             cout << "   read      ";   read_proj().print_stats();
             cout << "\n   theo      "; theo_proj().print_stats();
             cout << "\n   corr      "; corr_proj().print_stats();
             cout << "\n   alig      "; alig_proj().print_stats();
             cout << "\n   diff      "; diff_proj().print_stats();
             cout << "\n   subvol(0) "; (*ptr_vol_out)(0)().print_stats();
             cout << "\n";
          }

          if (prm.tell&TELL_SAVE_AT_EACH_STEP) {
      	     cout << "Stats PPPdiff.xmp: "; diff_proj().print_stats(); cout << endl;
      	     cout << "Stats PPPtheo.xmp: "; theo_proj().print_stats(); cout << endl;
      	     cout << "Stats PPPread.xmp: "; read_proj().print_stats(); cout << endl;
      	     cout << "Stats PPPcorr.xmp: "; corr_proj().print_stats(); cout << endl;
             diff_proj.write("PPPdiff.xmp"); 
             theo_proj.write("PPPtheo.xmp");
             read_proj.write("PPPread.xmp");
             corr_proj.write("PPPcorr.xmp");
	     if(act_proj!=0)
                alig_proj.write("PPPalign.xmp");
             vol_basis.write("PPPbasis.basis");
             prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
      	     cout << "Stats PPPvol.vol : "; vol_voxels().print_stats(); cout << endl;
             vol_voxels.write("PPPvol.vol");

	     if (!iv_launched) {
	        system("xmipp_show -img PPPdiff.xmp PPPtheo.xmp PPPread.xmp PPPcorr.xmp -dont_apply_geo -poll &");
		system("xmipp_show -vol PPPvol.vol -poll &");
		iv_launched=true;
	     }
             cout << "\nHit any key and enter\n"; char c; cin >> c;
         }

      	 // Save intermidiate
	 if ((prm.tell&TELL_SAVE_INTERMIDIATE | prm.tell&TELL_IV) &&
             prm.save_intermidiate_every!=0 &&
	     act_proj%prm.save_intermidiate_every==0) {
	    if (prm.tell&TELL_SAVE_INTERMIDIATE) 
               cerr << "\nSaving intermidiate ...\n"
                    << "Converting basis volume to voxels ...\n";   
            // Save reconstructed volume
            prm.basis.changeToVoxels(vol_basis, &(vol_voxels()), 
               Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
      	    if (prm.tell&TELL_SAVE_INTERMIDIATE)
               vol_voxels.write(prm.fn_root+"it"+ItoA(it)+"proj"+
		  ItoA(act_proj,5)+".vol");
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
	 prm.fh_hist << "Finished - Iteration " << it << endl;
         if (prm.WLS)
	    cerr        << "   Weighted error: " << global_mean_error
	                << " 1st block: "        << global_mean_error_1stblock
                        << " 2nd block: "        << mean_error_2ndblock
	                << " residual: "         << pow_residual_imgs << endl;
         else
            cerr << "   Global mean squared error: "
                 <<  global_mean_error/prm.numIMG << endl;

         if (prm.WLS)
	    prm.fh_hist << "   Weighted error: " << global_mean_error
	                << " 1st block: "        << global_mean_error_1stblock
                        << " 2nd block: "        << mean_error_2ndblock
	                << " residual: "         << pow_residual_imgs << endl;
         else
	    prm.fh_hist << "   Global mean squared error: "
                        << global_mean_error/prm.numIMG << endl;

	 if (POCS.apply_POCS) {
            cerr        << "   POCS Global mean squared error: "
                	<< POCS.POCS_global_mean_error/POCS.POCS_N << endl;
            prm.fh_hist << "   POCS Global mean squared error: "
                	<< POCS.POCS_global_mean_error/POCS.POCS_N << endl;
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
         if (rank==-1) cerr << "Converting basis volume to voxels ...\n";
         prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
            Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
         vol_voxels.write(prm.fn_root+"it"+ItoA(it)+".vol");

         if (prm.tell&TELL_SAVE_BASIS)
            vol_basis.write(prm.fn_root+"it"+ItoA(it)+".basis");
      }
      
      // Check if algorithm must stop via stop_at
      if (images==prm.stop_at) break;
   }

   // Times on screen
   if (rank==-1) {
      cout << "\nTime of " << prm.no_it << " iterations: \n";
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

// Produce side information and initial volume
   prm.produce_Side_Info(vol_basis);
//calculate symmetry in produce_Side_Info and exit if wrong sampling
//printting the right sampling
   eprm.produce_Side_Info(prm,vol_basis);

// Show parameters and initiate history
   Basic_ART_Init_history(prm, eprm, vol_basis);
   
// Iterations
   Basic_ART_iterations(prm, eprm, vol_basis);
   
// Finish iterations
   finish_ART_iterations(prm, eprm, vol_basis);

// Write final volume
   int Xoutput_volume_size=(prm.Xoutput_volume_size==0) ?
      prm.projXdim:prm.Xoutput_volume_size;
   int Youtput_volume_size=(prm.Youtput_volume_size==0) ?
      prm.projYdim:prm.Youtput_volume_size;
   int Zoutput_volume_size=(prm.Zoutput_volume_size==0) ?
      prm.projXdim:prm.Zoutput_volume_size;
   prm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
      Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
   vol_voxels.write(prm.fn_root+".vol");
   if (prm.tell&TELL_SAVE_BASIS) vol_basis.write(prm.fn_root+".basis");
   prm.fh_hist.close();

}
