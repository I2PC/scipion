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

#include "../Prog_art.hh"
#include "Basic_art.inc"
#include "../Prog_FourierFilter.hh"
#include <XmippData/xmippWavelets.hh>
#include <XmippData/Programs/Prog_denoising.hh>

/* ------------------------------------------------------------------------- */
/* Plain ART Parameters                                                      */
/* ------------------------------------------------------------------------- */
/* Produce Side information ------------------------------------------------ */
void Plain_ART_Parameters::produce_Side_Info(const Basic_ART_Parameters &prm,
   GridVolume &vol_blobs0) {}

/* Cout -------------------------------------------------------------------- */
ostream & operator << (ostream &o, const Plain_ART_Parameters &eprm) {
   return o;
}

/* ------------------------------------------------------------------------- */
/* Process correction                                                        */
/* ------------------------------------------------------------------------- */
void process_correction(Projection &corr_proj) {
   // Mask correction
   int Rmin=CEIL(MIN(XSIZE(corr_proj()),YSIZE(corr_proj()))/2);
   int R2=Rmin*Rmin;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(corr_proj()) 
      if (i*i+j*j>R2) corr_proj(i,j)=0;

   // Denoise the corrections
   Denoising_parameters denoiser;
   denoiser.denoising_type=Denoising_parameters::ADAPTIVE_SOFT;
   denoiser.scale=3;
   denoiser.output_scale=0;
   denoiser.produce_side_info();
   matrix2D<double> denoised=corr_proj();
   denoiser.denoise(denoised);
   ImageXmipp save; save()=corr_proj(); save.write("PPPbefore.xmp");
   if (denoised(0,0)==denoised(0,0))
      // This condition is not true if there are NaNs
      corr_proj()=denoised;
   save()=corr_proj(); save.write("PPPafter.xmp");
   cout << "Press\n";
   char c; cin >> c;
}

/* ------------------------------------------------------------------------- */
/* Update residual vector for WLS                                            */
/* ------------------------------------------------------------------------- */
void update_residual_vector(Basic_ART_Parameters &prm, GridVolume &vol_blobs, 
			    double &kappa, double &pow_residual_vol, double &pow_residual_imgs) {
  GridVolume       residual_vol;
  Projection       read_proj,dummy_proj,new_proj;
  FileName         fn_resi,fn_tmp;
  double           sqrtweight,dim2,norma,normb,apply_kappa;
  ImageOver        *footprint=(ImageOver *)&prm.blobprint;
  ImageOver        *footprint2=(ImageOver *)&prm.blobprint2;
  matrix2D<double> *A=NULL;
  vector<matrix2D<double> > newres_imgs;
  matrix2D<int>    mask;

  residual_vol.resize(vol_blobs);
  residual_vol.init_zeros();

  // Calculate volume from all backprojected residual images
  cerr << "Backprojection of residual images " <<endl;
  if (!(prm.tell&TELL_SHOW_ERROR)) init_progress_bar(prm.numIMG);

  for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++) {

    // backprojection of the weighted residual image
    sqrtweight=sqrt(prm.residual_imgs[iact_proj].weight()/prm.sum_weight);

    read_proj=prm.residual_imgs[iact_proj];
    read_proj()*=sqrtweight;
    dummy_proj().resize(read_proj());
    dummy_proj.set_angles(prm.IMG_Inf[iact_proj].rot,prm.IMG_Inf[iact_proj].tilt,prm.IMG_Inf[iact_proj].psi);

    project_Volume(residual_vol,prm.blob,*footprint,*footprint2,dummy_proj,
		   read_proj,YSIZE(read_proj()),XSIZE(read_proj()),
		   prm.IMG_Inf[iact_proj].rot,prm.IMG_Inf[iact_proj].tilt,prm.IMG_Inf[iact_proj].psi,BACKWARD,prm.eq_mode,
		   prm.GVNeq,NULL,prm.ray_length);

   if (!(prm.tell&TELL_SHOW_ERROR)) if (iact_proj%MAX(1,prm.numIMG/60)==0) progress_bar(iact_proj);
  }
  if (!(prm.tell&TELL_SHOW_ERROR)) progress_bar(prm.numIMG);

  // Convert to voxels: solely for output of power of residual volume
  VolumeXmipp      residual_vox;
  int Xoutput_volume_size=(prm.Xoutput_volume_size==0) ?
    prm.projXdim:prm.Xoutput_volume_size;
  int Youtput_volume_size=(prm.Youtput_volume_size==0) ?
    prm.projYdim:prm.Youtput_volume_size;
  int Zoutput_volume_size=(prm.Zoutput_volume_size==0) ?
    prm.projXdim:prm.Zoutput_volume_size;
  blobs2voxels(residual_vol, prm.blob, &residual_vox, prm.D,
               Zoutput_volume_size, Youtput_volume_size,
               Xoutput_volume_size);
  pow_residual_vol=residual_vox().sum2()/(Xoutput_volume_size*Youtput_volume_size*Zoutput_volume_size);
  residual_vox.clear();

  cerr << "Projection of residual volume; kappa = " << kappa<< endl;
  if (!(prm.tell&TELL_SHOW_ERROR)) init_progress_bar(prm.numIMG);

  // Now that we have the residual volume: project in all directions
  pow_residual_imgs=0.;
  new_proj().resize(read_proj()); 
  mask.resize(read_proj());
  BinaryCircularMask(mask,YSIZE(read_proj())/2,INNER_MASK);

  dim2=(double)YSIZE(read_proj())*XSIZE(read_proj());
  for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++) {

    project_Volume(residual_vol,prm.blob,*footprint,*footprint2,new_proj,
		   dummy_proj,YSIZE(read_proj()),XSIZE(read_proj()),
		   prm.IMG_Inf[iact_proj].rot,prm.IMG_Inf[iact_proj].tilt,prm.IMG_Inf[iact_proj].psi,FORWARD,prm.eq_mode,
		   prm.GVNeq,A,prm.ray_length);

    sqrtweight=sqrt(prm.residual_imgs[iact_proj].weight()/prm.sum_weight);

    // Next lines like normalization in [EHL] (2.18)?
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(new_proj()) {
      dMij(dummy_proj(),i,j)=MAX(1.,dMij(dummy_proj(),i,j)); // to avoid division by zero
      dMij(new_proj(),i,j)/=dMij(dummy_proj(),i,j);
    }
    new_proj()*=sqrtweight*kappa;

    /*
    fn_tmp="residual_"+ItoA(iact_proj);
    dummy_proj()=1000*prm.residual_imgs[iact_proj]();
    dummy_proj.write(fn_tmp+".old");
    */

    prm.residual_imgs[iact_proj]()-=new_proj();
    pow_residual_imgs+=prm.residual_imgs[iact_proj]().sum2();

    // Mask out edges of the images
    apply_binary_mask(mask,prm.residual_imgs[iact_proj](),prm.residual_imgs[iact_proj](),0.);

    /*
    dummy_proj()=1000*new_proj();
    dummy_proj.write(fn_tmp+".change");
    dummy_proj()=1000*prm.residual_imgs[iact_proj]();
    dummy_proj.write(fn_tmp+".new");
    */

    if (!(prm.tell&TELL_SHOW_ERROR)) if (iact_proj%MAX(1,prm.numIMG/60)==0) progress_bar(iact_proj);
  }

  pow_residual_imgs/=dim2;
  newres_imgs.clear();

  if (!(prm.tell&TELL_SHOW_ERROR)) progress_bar(prm.numIMG);

}

/* ------------------------------------------------------------------------- */
/* ART Single step                                                           */
/* ------------------------------------------------------------------------- */
#define blob       prm.blob
void ART_single_step(
   GridVolume        	   &vol_in,          // Input Reconstructed volume
   GridVolume        	   *vol_out,         // Output Reconstructed volume
   const Basic_ART_Parameters &prm,          // blob, lambda
   const Plain_ART_Parameters &eprm,         // In this case, nothing
   Projection        	   &theo_proj,       // Projection of the reconstruction
                                             // It is outside to make it visible
                                             // just if it needed for any
                                             // other purpose
   /*const*/ Projection  	   &read_proj,       // Real projection
   int sym_no,                               // Symmetry matrix index
   Projection        	   &diff_proj,       // Difference between read and
                                             // theoretical projection
   Projection        	   &corr_proj,       // Correcting projection
   Projection              &alig_proj,       // Translation alignement aux proj
   double             	   &mean_error,      // Mean error over the pixels
   int               	   numIMG,           // number of images in the set
                                             // in SIRT the correction must
                                             // be divided by this number
   double                  lambda,           // Lambda to be used
   int                     act_proj,         // Projection number
   const FileName          &fn_ctf,          // CTF to apply
   bool                    unmatched,        // Apply unmatched projectors
   double                  ray_length,       // Ray length for the projection
   bool                    print_system_matrix) // Print matrix (A in Ax=b) of
                                             // the equation system, as well as
					     // the independent vector (b)
{
// Prepare to work with CTF ................................................
   FourierMask ctf;
   ImageOver *footprint=(ImageOver *)&prm.blobprint;
   ImageOver *footprint2=(ImageOver *)&prm.blobprint2;
   bool remove_footprints=false;
   double weight,sqrtweight;

   if (fn_ctf!="" && !unmatched)
      if (Is_FourierImageXmipp(fn_ctf)) ctf.read_mask(fn_ctf);
      else {
	 // It is a description of the CTF
	 ctf.FilterShape=ctf.FilterBand=CTF;
	 ctf.ctf.read(fn_ctf);
	 ctf.ctf.Tm/=BLOB_SUBSAMPLING;
	 ctf.ctf.Produce_Side_Info();

	 // Create new footprints
	 footprint=new ImageOver;
	 footprint2=new ImageOver;
	 remove_footprints=true;

	 // Enlarge footprint, bigger than necessary to avoid
	 // aliasing
	 *footprint=prm.blobprint;
	 (*footprint)().set_Xmipp_origin();
	 int finalsize=2*CEIL(30+blob.radius)+1;
	 footprint->window(
            FIRST_XMIPP_INDEX(finalsize),FIRST_XMIPP_INDEX(finalsize),
            LAST_XMIPP_INDEX(finalsize), LAST_XMIPP_INDEX(finalsize));

	 // Apply CTF
	 ctf.apply_mask_Space((*footprint)());

	 // Remove unnecessary regions
	 finalsize=2*CEIL(15+blob.radius)+1;
	 footprint->window(
            FIRST_XMIPP_INDEX(finalsize),FIRST_XMIPP_INDEX(finalsize),
            LAST_XMIPP_INDEX(finalsize), LAST_XMIPP_INDEX(finalsize));
	 #ifdef DEBUG
            ImageXmipp save; save()=(*footprint)();
            save.write("PPPfootprint.xmp");
	 #endif

	 // Create footprint2
	 *footprint2=*footprint;
	 (*footprint2)()*=(*footprint2)();
      }

// Project structure .......................................................
   // The correction image is reused in this call to store the normalising
   // projection, ie, the projection of an all-1 volume
   matrix2D<double> *A=NULL;
   if (print_system_matrix) A=new matrix2D<double>;
   corr_proj().init_zeros();
   project_Volume(vol_in,blob,*footprint,*footprint2,theo_proj,
      corr_proj,YSIZE(read_proj()),XSIZE(read_proj()),
      read_proj.rot(),read_proj.tilt(),read_proj.psi(),FORWARD,prm.eq_mode,
      prm.GVNeq,A,prm.ray_length);

   if (fn_ctf!="" && unmatched)
      if (Is_FourierImageXmipp(fn_ctf)) {
         // If it is a fft file
         ctf.read_mask(fn_ctf);
         ctf.apply_mask_Space(theo_proj());
      }

   // Print system matrix
   if (print_system_matrix) {
      cout << "Equation system (Ax=b) ----------------------\n";
      cout << "Size: "; A->print_shape(); cout << endl;
      for (int i=0; i<YSIZE(*A); i++) {
         bool null_row=true;
         for (int j=0; j<YSIZE(*A); j++)
	    if (DIRECT_MAT_ELEM(*A,i,j)!=0) {null_row=false; break;}
	 if (!null_row) {
	    cout << "pixel=" << ItoA(i,3) << " --> "
	         << MULTIDIM_ELEM(read_proj(),i) << " = ";
            for (int j=0; j<XSIZE(*A); j++)
	       cout << DIRECT_MAT_ELEM(*A,i,j) << " ";
	    cout << endl;
	 }
      }
      cout << "---------------------------------------------\n";
      delete A;
   }

   // Now compute differences .................................................
   double applied_lambda=lambda/numIMG; // In ART mode, numIMG=1 

   mean_error=0;
   diff_proj().resize(read_proj());
   
   // Weighted least-squares ART for Maximum-Likelihood refinement
   if (prm.WLS) { 
     weight=read_proj.weight()/prm.sum_weight;
     sqrtweight=sqrt(weight);

     FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(read_proj)) {
       // Compute difference image and error
       IMGPIXEL(diff_proj,i,j)=IMGPIXEL(read_proj,i,j)-IMGPIXEL(theo_proj,i,j);
       mean_error += IMGPIXEL(diff_proj,i,j) * IMGPIXEL(diff_proj,i,j);
       
       // Subtract the residual image (stored in alig_proj!)
       IMGPIXEL(diff_proj,i,j)=sqrtweight*IMGPIXEL(diff_proj,i,j)-IMGPIXEL(alig_proj,i,j);
       
       // Calculate the correction and the updated residual images
       IMGPIXEL(corr_proj,i,j)=
         applied_lambda*IMGPIXEL(diff_proj,i,j)/(weight*IMGPIXEL(corr_proj,i,j) + 1.);
       IMGPIXEL(alig_proj,i,j)+=IMGPIXEL(corr_proj,i,j);
       IMGPIXEL(corr_proj,i,j)*=sqrtweight;

     }
     mean_error /= XSIZE(diff_proj())*YSIZE(diff_proj());
     mean_error*=weight;

   } else {

     FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(read_proj)) {
       // Compute difference image and error
       IMGPIXEL(diff_proj,i,j)=IMGPIXEL(read_proj,i,j)-IMGPIXEL(theo_proj,i,j);
       mean_error += IMGPIXEL(diff_proj,i,j) * IMGPIXEL(diff_proj,i,j);

       // Compute the correction image
       IMGPIXEL(corr_proj,i,j)=MAX(IMGPIXEL(corr_proj,i,j),1);
       IMGPIXEL(corr_proj,i,j)=
         applied_lambda*IMGPIXEL(diff_proj,i,j)/IMGPIXEL(corr_proj,i,j);
     }
     mean_error /= XSIZE(diff_proj())*YSIZE(diff_proj());
   }
   
   // Denoising of the correction image
   //process_correction(corr_proj);

   // Backprojection of correction plane ......................................
   project_Volume(*vol_out,blob,*footprint,*footprint2,theo_proj,
      corr_proj,YSIZE(read_proj()),XSIZE(read_proj()),
      read_proj.rot(),read_proj.tilt(),read_proj.psi(),BACKWARD,prm.eq_mode,
      prm.GVNeq,NULL,prm.ray_length);

   // Remove footprints if necessary
   if (remove_footprints) {
      delete footprint;
      delete footprint2;
   }
}
#undef blob

/* Instantiation of the ART process ---------------------------------------- */
void instantiate_Plain_ART() {
   Basic_ART_Parameters prm;
   Plain_ART_Parameters eprm;
   VolumeXmipp vol_voxels;
   GridVolume  vol_blobs;
   Basic_ROUT_Art(prm,eprm,vol_voxels, vol_blobs);
}

/* Finish iterations ------------------------------------------------------- */
void finish_ART_iterations(const Basic_ART_Parameters &prm,
   const Plain_ART_Parameters &eprm, GridVolume &vol_blobs) { }

/* Apply_symmetry ------------------------------------------------------- */
void apply_symmetry(GridVolume &vol_in, GridVolume *vol_out,
                    const Plain_ART_Parameters &eprm, int grid_type)
    {cout <<"\nERROR: Function not implemented for single particles"<<endl;
     }
