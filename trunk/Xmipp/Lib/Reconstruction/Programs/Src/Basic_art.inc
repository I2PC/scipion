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

#include "../Basic_art.hh"
#include "../Prog_symmetrize.hh"
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippHistograms.hh>
#include <XmippData/xmippMasks.hh>

/* ------------------------------------------------------------------------- */
/* Init history                                                              */
/* ------------------------------------------------------------------------- */
template <class Extra_ART_Parameters>
void Basic_ART_Init_history(Basic_ART_Parameters &prm,
   const Extra_ART_Parameters &eprm, const GridVolume &vol_blobs0) {
   
// Show general information ................................................
   prm.fh_hist << " RECONSTRUCTION HISTORY: \n\n";
   prm.fh_hist << " Parameters -------------------------------------------------\n";
   prm.fh_hist << " Projections file: " << prm.fn_sel << endl;
   prm.fh_hist << " CTF file:         " << prm.fn_ctf << endl;
   prm.fh_hist << " Sampling= " << prm.sampling << endl;   
   prm.fh_hist << " Blobs:         radius=" << prm.blob.radius 
      << " pixels" << " alpha="
      << prm.blob.alpha << " order=" << prm.blob.order << endl;
   switch (prm.grid_type) {
      case FCC: prm.fh_hist << " Grid Type: FCC "; break;
      case BCC: prm.fh_hist << " Grid Type: BCC "; break;
      case CC:  prm.fh_hist << " Grid Type: CC ";  break;
   }
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
   prm.fh_hist << " Parallel mode: ";
   switch (prm.parallel_mode) {
      case Basic_ART_Parameters::ART:   prm.fh_hist << "ART\n"; break;
      case Basic_ART_Parameters::SIRT:  prm.fh_hist << "SIRT\n"; break;
      case Basic_ART_Parameters::SART:
          prm.fh_hist << "SART, block size=" << prm.block_size << endl;
          break;
      case Basic_ART_Parameters::AVSP:  prm.fh_hist << "AVSP\n"; break;
      case Basic_ART_Parameters::BiCAV:
          prm.fh_hist << "BiCAV, block size=" << prm.block_size << endl;
          break;
   }
   prm.fh_hist << " Equation mode: ";
   if (prm.eq_mode==CAV) prm.fh_hist << "CAV\n";
   else if (prm.eq_mode==CAVK) prm.fh_hist << "CAVK\n";
   else prm.fh_hist << "ARTK\n";
   prm.fh_hist << " Projections: Ydim x Xdim = " << prm.projYdim
      << " x " << prm.projXdim << endl;
   prm.fh_hist << " Surface mask: " << prm.fn_surface_mask << endl;
   prm.fh_hist << " POCS freq: " << prm.POCS_freq << endl;
   prm.fh_hist << " Known volume: " << prm.known_volume << endl;
   prm.fh_hist << " Positivity: "; print(prm.fh_hist,prm.positivity);
      prm.fh_hist << endl;
   prm.fh_hist << " Symmetry file: " << prm.fn_sym << endl;
   prm.fh_hist << " Force Symmetry each: " << prm.sym_each << " projections" 
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

   prm.fh_hist << "Grid structure ------\n" << vol_blobs0.grid();

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
   const Extra_ART_Parameters &eprm, GridVolume &vol_blobs,
   int rank) {
   // Some variables .......................................................
   int        ART_numIMG;              // Normalizing factor
   GridVolume *ptr_vol_out;            // Pointer to output volume
   GridVolume vol_blobs_out;           // Output volume (only useful in SIRT)
   VolumeXmipp vol_voxels;             // This one is useful only in the
                                       // case of saving intermidiate volumes
   VolumeXmipp vol_POCS;               // Used for applying POCS
   bool        apply_POCS=FALSE;
   int Xoutput_volume_size, Youtput_volume_size, Zoutput_volume_size;
   
   // Projection related ...................................................
   Projection      read_proj;          // Projection read from file
   Projection      theo_proj;          // Projection from the
                                       // reconstruction
   Projection      alig_proj;          // Projection with the correlation
                                       // maps needed for translation aligment				       
   Projection      corr_proj;          // Image with the correction
                                       // factors for unitary blobs
   Projection      diff_proj;          // Difference between the
                                       // theoretical and real image
   vtkImageData    *FFTtheo=NULL;      // Fourier transform of a projection

   // Reconstruction results ...............................................
   double          mean_error;
   double          global_mean_error;

   // POCS constraints .....................................................
   VolumeXmipp     theo_POCS_vol, corr_POCS_vol;
   double          POCS_mean_error, POCS_max_error, POCS_global_mean_error;
   #define POCS_measuring  1
   #define POCS_use        2
   #define POCS_lowering   3
   #define POCS_N_measure  8
   #define POCS_N_use      2
   int             POCS_state=POCS_measuring;
   int             POCS_freq=prm.POCS_freq, POCS_i=0;
   int             POCS_vec_i=0, POCS_used=0, POCS_N=0;
   matrix1D<double> POCS_errors(POCS_N_measure);
   double          POCS_avg, POCS_stddev, POCS_min, POCS_max;

   // Some initialisation ..................................................
   Xoutput_volume_size=(prm.Xoutput_volume_size==0) ?
      prm.projXdim:prm.Xoutput_volume_size;
   Youtput_volume_size=(prm.Youtput_volume_size==0) ?
      prm.projYdim:prm.Youtput_volume_size;
   Zoutput_volume_size=(prm.Zoutput_volume_size==0) ?
      prm.projXdim:prm.Zoutput_volume_size;
   POCS_errors.init_zeros();

   // If SIRT set normalizing factor and create output volume ..............
   if (prm.parallel_mode==Basic_ART_Parameters::SIRT ||
       prm.parallel_mode==Basic_ART_Parameters::SART) {
      ART_numIMG=prm.numIMG;           // Normalizing factor = Total no. images
      ptr_vol_out=&vol_blobs_out;      // Pointer to output volume
      vol_blobs_out=vol_blobs;         // Copy the structure of vol_blobs
                                       // and possible initial values
   } else if (prm.eq_mode==CAV) {
      ART_numIMG=1;                    // Normalizing factor = Total no. images
      ptr_vol_out=&vol_blobs_out;      // Pointer to output volume
      vol_blobs_out=vol_blobs;         // Copy the structure of vol_blobs
                                       // and possible initial values
   } else {
      ART_numIMG=1;                    // No normalizing factor
      ptr_vol_out=&vol_blobs;          // Output volume is the same as
                                       // input one
   }
   // Now iterate ..........................................................
   TimeStamp time0;                    // For measuring the elapsed time
   annotate_time(&time0);
   int images=0;
    for (int it = 0; it < prm.no_it; it++) {
      // Initialisation of some variables
      global_mean_error = 0;
      POCS_global_mean_error = 0;

      // For each projection -----------------------------------------------
      if (rank==-1) {
	 cerr << "Running iteration " << it << " with lambda= " << prm.lambda(it)<< "\n" << endl;
	 if (!(prm.tell&TELL_SHOW_ERROR)) init_progress_bar(prm.numIMG);
      }
      
      for (int act_proj = 0; act_proj < prm.numIMG ; act_proj++) {
          POCS_N=0;
  
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
	  read_proj.read(prm.IMG_Inf[iact_proj].fn_proj);
          read_proj.move_origin_to_center();
          read_proj.rot ()=prm.IMG_Inf[iact_proj].rot;
          read_proj.tilt()=prm.IMG_Inf[iact_proj].tilt;
          read_proj.psi ()=prm.IMG_Inf[iact_proj].psi;

          // Projection extension? .........................................
          if (prm.proj_ext!=0)
             read_proj().window(
                STARTINGY (read_proj())-prm.proj_ext,
                STARTINGX (read_proj())-prm.proj_ext,
                FINISHINGY(read_proj())+prm.proj_ext,
                FINISHINGX(read_proj())+prm.proj_ext);

         //Skip if desired	  
          if (prm.tell&TELL_ONLY_SYM) 
	     {
	     if( prm.IMG_Inf[iact_proj].sym != -1)
	        {
		cerr << "Skipping Proj no: " << iact_proj 
		     << " with symmetry no: " << prm.IMG_Inf[iact_proj].sym 
		     <<  endl;
		continue;
		}
	     else
	        {
		cerr << "NO Skipping Proj no: " << iact_proj 
		     << " with symmetry no: " << prm.IMG_Inf[iact_proj].sym 
		     <<  endl;
		}
	     	
	      }
		
          // Apply the reconstruction algorithm ............................
      	  // Notice that the following function is extern
	  ART_single_step(vol_blobs, ptr_vol_out,
      	      prm, eprm,
              theo_proj, read_proj, prm.IMG_Inf[iact_proj].sym , diff_proj, 
	      corr_proj, alig_proj,
              mean_error, ART_numIMG, prm.lambda(it),
	      images, prm.IMG_Inf[iact_proj].fn_ctf, prm.print_system_matrix);

          global_mean_error += mean_error;
          
              // Force symmetry, so far only crystallographic

              //Force symmetry in the volume each
              if (prm.sym_each &&
	          act_proj%prm.sym_each==0 &&
		  (act_proj!=0 || prm.sym_each==1) ) {
                 apply_symmetry(vol_blobs,ptr_vol_out,eprm, prm.grid_type);
              }
          // Apply POCS ....................................................
          if ((prm.surface_mask!=NULL || 
	      prm.positivity || prm.force_sym!=0 ||
	      prm.known_volume!=-1) && POCS_i%POCS_freq==0) {
	      apply_POCS=TRUE;
              VolumeXmipp vol_aux;
              VolumeXmipp *desired_volume=NULL;

	      // Compute the corresponding voxel volume
              blobs2voxels(vol_blobs, prm.blob, &vol_voxels, prm.D, 
                 Zoutput_volume_size, Youtput_volume_size,
		 Xoutput_volume_size);
              if (prm.tell&TELL_SAVE_AT_EACH_STEP) {
		 vol_voxels.write("PPPvolPOCS0.vol");
		 cout << "Stats PPPvolPOCS0.vol: "; vol_voxels().print_stats(); cout << endl;
	      }
              // Apply surface restriction
	      if (prm.surface_mask!=NULL) {
      	         vol_POCS()=(*(prm.surface_mask))();
	      } else {
		 vol_POCS().resize(vol_voxels());
		 vol_POCS().init_zeros();
              }
              if (prm.tell&TELL_SAVE_AT_EACH_STEP) {
		 vol_POCS.write("PPPvolPOCS1.vol");
		 cout << "Stats PPPvolPOCS1.vol: "; vol_POCS().print_stats(); cout << endl;
      	      }
              // Force symmetry
              if (prm.force_sym!=0) {
                 symmetrize(prm.SL,vol_voxels,vol_aux,FALSE);
                 desired_volume=&vol_aux;
                 if (prm.tell&TELL_SAVE_AT_EACH_STEP) {
		    vol_aux.write("PPPvolPOCS2.vol");
                    cout << "Stats PPPvolPOCS2.vol: "; vol_aux().print_stats(); cout << endl;
	         }
              }
	      // Apply volume constraint
	      if (prm.known_volume!=-1) {
	         histogram1D hist;
                 matrix3D<int> aux_mask;
                 aux_mask.resize(vol_POCS());
                 FOR_ALL_ELEMENTS_IN_MATRIX3D(aux_mask)
                    aux_mask(k,i,j)=1-(int)vol_POCS(k,i,j);
                 long mask_voxels=vol_POCS().count_threshold("below",0.5,0);
	         compute_hist_within_binary_mask(
                    aux_mask,vol_voxels(),hist,300);
	         double known_percentage;
	         known_percentage=MIN(100,100*prm.known_volume/mask_voxels);
	         double threshold;
	         threshold=hist.percentil(100-known_percentage);
	         FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_voxels())
	            if (vol_voxels(k,i,j)<threshold) vol_POCS(k,i,j)=1;
	      }
              if (prm.tell&TELL_SAVE_AT_EACH_STEP) {
		 vol_POCS.write("PPPvolPOCS3.vol");
		 cout << "Stats PPPvolPOCS3.vol: "; vol_POCS().print_stats(); cout << endl;
	      }

              // Do not allow positivity outside interest region
              // and do not allow negativity inside the interest region
              // if positivity restrictions are to be applied
              int bg=(int) vol_POCS().sum();
              int fg=MULTIDIM_SIZE(vol_POCS())-bg;
              int relax=0, posi=0;
	      FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_voxels())
	         if (vol_POCS(k,i,j)==1 && vol_voxels(k,i,j)<0)
                    {vol_POCS(k,i,j)=0; relax++;}
                 else if (vol_POCS(k,i,j)==0 && vol_voxels(k,i,j)<0 &&
                    prm.positivity)
                    {vol_POCS(k,i,j)=1; posi++;}
               // Debugging messages
               //cerr << "Relaxation/Positivity " << (double)relax/(double)bg << " "
               //     << (double)posi/(double)fg << " " << endl;

              // Solve volumetric equations
              if (desired_volume==NULL)
                 ART_voxels2blobs_single_step(vol_blobs, &vol_blobs,
                    prm.blob, prm.D, prm.lambda(it),
	            (Volume *) &theo_POCS_vol, NULL,
	            (Volume *) &corr_POCS_vol,
                    (Volume *) &vol_POCS,
                    POCS_mean_error, POCS_max_error, VARTK);
              else {
                 FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_POCS())
                    if (vol_POCS(k,i,j)==1) (*desired_volume)(k,i,j)=0;
                 for (int i=0; i<prm.force_sym; i++) {
                    ART_voxels2blobs_single_step(vol_blobs, &vol_blobs,
                       prm.blob, prm.D, prm.lambda(it),
	               (Volume *) &theo_POCS_vol, desired_volume,
	               (Volume *) &corr_POCS_vol,
                       NULL,
                       POCS_mean_error, POCS_max_error, VARTK);
                    if (prm.tell&TELL_SAVE_AT_EACH_STEP)
		       cout << "    POCS Iteration " << i
		            << " POCS Error=" <<  POCS_mean_error << endl;
                 }
              }
              POCS_i=1;
              POCS_global_mean_error+=POCS_mean_error;
              POCS_N++;

              // Now some control logic
              if (prm.numIMG-images<100 || images%100==0 ||
	          desired_volume!=NULL) {
                 POCS_freq=1;
                 POCS_state=POCS_measuring;
                 POCS_vec_i=0;
              } else {
                 double dummy;
                 switch (POCS_state) {
                    case POCS_measuring:
                       #ifdef DEBUG_POCS
                          cout << "M:" << POCS_vec_i << " " << POCS_mean_error << endl;
                       #endif
                       POCS_errors(POCS_vec_i++)=POCS_mean_error;
                       if (POCS_vec_i==POCS_N_measure) {
                          POCS_vec_i=0;
                          // Change to use state
                          POCS_used=0;
                          POCS_freq++;
                          POCS_state=POCS_use;
                          #ifdef DEBUG_POCS
                             cerr << "1: Changing to " << POCS_freq << endl; 
                          #endif
                       }
                       break;
                    case POCS_use:
                       POCS_used++;
                       POCS_errors.compute_stats(POCS_avg,
                          POCS_stddev, dummy, POCS_min);
                       #ifdef DEBUG_POCS
                          cout << "Reference errors: " << POCS_errors.transpose() << endl;
                          cout << "Checking " << ABS(POCS_mean_error-POCS_avg) << " " << 1.2*1.96*POCS_stddev << endl;
                       #endif
                       if (ABS(POCS_mean_error-POCS_avg)<1.2*1.96*POCS_stddev) {
                          if (POCS_mean_error<POCS_avg) {
                             double max_error=POCS_errors(0);
                             POCS_vec_i=0;
                             for (int i=1; i<POCS_N_measure; i++)
                                 if (POCS_errors(i)>max_error) {
                                    max_error=POCS_errors(i);
                                    POCS_vec_i=i;
                                 }
                             POCS_errors(POCS_vec_i)=POCS_mean_error;
                          }
                          if (POCS_used<POCS_N_use) { // While not enough uses
                          } else if (POCS_freq<3) { // increase frequency
                             POCS_freq++;
                             #ifdef DEBUG_POCS
                                cerr << "2: Changing to " << POCS_freq << endl;
                             #endif
                             POCS_used=0;
                          }
                       } else {
                          // It is behaving worse
                          if (POCS_freq>prm.POCS_freq+1) {
                             POCS_freq=prm.POCS_freq+1;
                             POCS_used=0;
                             #ifdef DEBUG_POCS
                                cerr << "3: Changing to " << POCS_freq << endl;
                             #endif
                          } else if (POCS_used>2) {
                             POCS_freq=prm.POCS_freq;
                             // Change status
                             POCS_used=0;
                             POCS_state=POCS_lowering;
                             #ifdef DEBUG_POCS
                                cerr << "Lowering\n";
                             #endif
                          }
                       }
                       break;
                    case POCS_lowering:
                       // Lower the POCS error before measuring again
                       POCS_errors.compute_stats(POCS_avg,
                          POCS_stddev, POCS_max, dummy);
                       POCS_used++;
                       if (POCS_mean_error<POCS_max || POCS_used>2*POCS_N_measure) {
                             // Change status
                             POCS_vec_i=0;
                             POCS_state=POCS_measuring;
                       }
                       break;
                 }
               }
          } else {POCS_i++; POCS_mean_error=-1;}

          // Show results ..................................................
          prm.fh_hist << prm.IMG_Inf[iact_proj].fn_proj << ", sym="
             << prm.IMG_Inf[iact_proj].sym << "\t\t" << mean_error;
          if (apply_POCS)
             prm.fh_hist << "\tPOCS:" << POCS_mean_error;
          prm.fh_hist << endl;
          if (prm.tell&TELL_SHOW_ERROR) {
             cout << prm.IMG_Inf[iact_proj].fn_proj << ", sym="
                << prm.IMG_Inf[iact_proj].sym << "\t\t" << mean_error;
             if (apply_POCS)
                cout        << "\tPOCS:" << POCS_mean_error;
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
             vol_blobs.write("PPPblobs.blob");
             blobs2voxels(vol_blobs, prm.blob, &vol_voxels, prm.D,
                Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
      	     cout << "Stats PPPvol.vol : "; vol_voxels().print_stats(); cout << endl;
             vol_voxels.write("PPPvol.vol");
	     //VolumeXmipp Vcoefs;
	     //blobs2space_coefficients(vol_blobs, prm.blob, &Vcoefs);
	     //Vcoefs.write((string)"PPPblob_coefs.vol");
             cout << "\nHit any key and enter\n"; char c; cin >> c;
         }

	 if (prm.tell&TELL_SAVE_INTERMIDIATE &&
	     act_proj%prm.save_intermidiate_every==0) {
            cerr << "\nSaving intermidiate ...\n"
                 << "Converting blob volume to voxels ...\n";   
            blobs2voxels(vol_blobs, prm.blob, &vol_voxels, prm.D, 
               Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
            vol_voxels.write(prm.fn_root+"it"+ItoA(it)+"proj"+
	       ItoA(act_proj,5)+".vol");	       
	 }

         // Check if algorithm must stop via stop_at
         if (++images==prm.stop_at) break;
      }

      if (!(prm.tell&TELL_SHOW_ERROR)) progress_bar(prm.numIMG);

      // Prepare for next global iteration ---------------------------------
      // Calculate norm and print
      if (rank==-1) {
	 prm.fh_hist << "Finished - Iteration " << it << endl;
	 cerr        <<"   Global mean squared error: "
                     << global_mean_error/prm.numIMG << endl;
	 prm.fh_hist << "   Global mean squared error: "
                     << global_mean_error/prm.numIMG << endl;
	 if (apply_POCS) {
            cerr        << "   POCS Global mean squared error: "
                	<< POCS_global_mean_error/POCS_N << endl;
            prm.fh_hist << "   POCS Global mean squared error: "
                	<< POCS_global_mean_error/POCS_N << endl;
	 }
      }

      // Convert volume and write if not last iteration
      // If in SIRT mode move the volume to the reference one
      if (prm.parallel_mode==Basic_ART_Parameters::SIRT || prm.eq_mode==CAV ||
          prm.parallel_mode==Basic_ART_Parameters::SART)
	  vol_blobs=vol_blobs_out;

      if (prm.tell&TELL_SAVE_INTERMIDIATE && it!=prm.no_it-1) {
         if (rank==-1) cerr << "Converting blob volume to voxels ...\n";   
         blobs2voxels(vol_blobs, prm.blob, &vol_voxels, prm.D, 
            Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
         vol_voxels.write(prm.fn_root+"it"+ItoA(it)+".vol");

         if (prm.tell&TELL_SAVE_BLOBS)
            vol_blobs.write(prm.fn_root+"it"+ItoA(it)+".blob");
      }

      // Check if algorithm must stop via stop_at
      if (images==prm.stop_at) break;
   }
   // Times on screen
   if (rank==-1) {
      cout << "\nTime of " << prm.no_it << " iterations: \n";
      print_elapsed_time(time0);
   }
}

/* ------------------------------------------------------------------------- */
/* ROUTINE ART                                                               */
/* ------------------------------------------------------------------------- */
template <class Extra_ART_Parameters>
void Basic_ROUT_Art(Basic_ART_Parameters &prm,
   Extra_ART_Parameters &eprm, VolumeXmipp &vol_voxels,
   GridVolume &vol_blobs) {
// Configure time clock
   time_config();

// Produce side information and initial volume
   prm.produce_Side_Info(vol_blobs);
//calculate symmetry in produce_Side_Info and exit if wrong sampling
//printting the right sampling
   eprm.produce_Side_Info(prm,vol_blobs);

// Show parameters and initiate history
   Basic_ART_Init_history(prm, eprm, vol_blobs);
   
// Iterations
   Basic_ART_iterations(prm, eprm, vol_blobs);
   
// Finish iterations
   finish_ART_iterations(prm, eprm, vol_blobs);

// Write final volume
   int Xoutput_volume_size=(prm.Xoutput_volume_size==0) ?
      prm.projXdim:prm.Xoutput_volume_size;
   int Youtput_volume_size=(prm.Youtput_volume_size==0) ?
      prm.projYdim:prm.Youtput_volume_size;
   int Zoutput_volume_size=(prm.Zoutput_volume_size==0) ?
      prm.projXdim:prm.Zoutput_volume_size;
   blobs2voxels(vol_blobs, prm.blob, &vol_voxels, prm.D,
      Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
   vol_voxels.write(prm.fn_root+".vol");
   if (prm.tell&TELL_SAVE_BLOBS) vol_blobs.write(prm.fn_root+".blob");
      
   prm.fh_hist.close();

}
