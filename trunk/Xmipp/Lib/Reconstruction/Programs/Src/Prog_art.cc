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
void process_correction(Projection &corr_proj, 
   Projection &theo_proj, GridVolume *vol_var) {
   // Mask correction
   double R2=vol_var->grid(0).R2;
   if (R2!=-1)
      FOR_ALL_ELEMENTS_IN_MATRIX2D(corr_proj())
         if (i*i+j*j>R2) corr_proj(i,j)=0;

   // Denoise correction
   set_DWT_type(DAUB12);
   DWT(corr_proj(),corr_proj());
   bayesian_wiener_filtering(corr_proj(),1);
   clean_quadrant(corr_proj(),0,"01");
   clean_quadrant(corr_proj(),0,"10");
   clean_quadrant(corr_proj(),0,"11");
   IDWT(corr_proj(),corr_proj());

   // Robust estimation of the variance
   double min, max, avg, stddev, stddev_ant;
   corr_proj().compute_stats(avg, stddev, min, max);
   int N=1;
   do {
      stddev_ant=stddev;
      max=3*stddev;
      min=-max;

      double sum=0, sum2=0;
      int    N_accounted=0;

      FOR_ALL_ELEMENTS_IN_MATRIX2D(corr_proj()) 
         if (corr_proj(i,j)>=min && corr_proj(i,j)<=max) {
            N_accounted++;
	    sum+=corr_proj(i,j);
	    sum2+=corr_proj(i,j)*corr_proj(i,j);
         }

      if (N_accounted!=0) {
         sum2/=N_accounted;
         sum /=N_accounted;
         stddev=sqrt(sum2-sum*sum);
      } else stddev=0;

      N++;
   } while (ABS(stddev-stddev_ant)/stddev>0.01 && N<10);

   // Remove outliers
   min=avg-2*stddev;
   max=avg+2*stddev;
   double avg_theo=theo_proj().compute_avg();
   double min_theo=theo_proj().compute_avg();
   double th_theo=0.5*(avg_theo+min_theo);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(corr_proj())
      if (corr_proj(i,j)<min || corr_proj(i,j)>max)
          corr_proj(i,j)=0;
      //else if (theo_proj(i,j)<th_theo) corr_proj(i,j)=0;

   // Low pass filter the result
   FourierMask Filter;
   Filter.FilterShape=RAISED_COSINE;
   Filter.FilterBand=LOWPASS;
   Filter.w1=0.25;
   Filter.raised_w=0.02;
   Filter.apply_mask_Space(corr_proj());
   
   // High pass filter the variance tracking
   FourierMask Filter2;
   Filter2.FilterShape=RAISED_COSINE;
   Filter2.FilterBand=HIGHPASS;
   Filter2.w1=0.04;
   Filter2.raised_w=0.02;
   Filter2.apply_mask_Space((*vol_var)(0)());
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
   GridVolume              *vol_var,         // Keep track of the variance
   bool                    print_system_matrix) // Print matrix (A in Ax=b) of
                                             // the equation system, as well as
					     // the independent vector (b)
{
// Prepare to work with CTF ................................................
   FourierMask ctf;
   ImageOver *footprint=(ImageOver *)&prm.blobprint;
   ImageOver *footprint2=(ImageOver *)&prm.blobprint2;
   bool remove_footprints=false;

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
   project_Volume(vol_in,*footprint,*footprint2,theo_proj,
      corr_proj,YSIZE(read_proj()),XSIZE(read_proj()),
      read_proj.rot(),read_proj.tilt(),read_proj.psi(),FORWARD,prm.eq_mode,
      prm.GVNeq,A,vol_var);

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
   
   FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(read_proj)) {
   // Compute difference image and error
      IMGPIXEL(diff_proj,i,j)=IMGPIXEL(read_proj,i,j)-IMGPIXEL(theo_proj,i,j);
      mean_error += IMGPIXEL(diff_proj,i,j) * IMGPIXEL(diff_proj,i,j);

   // Compute the correction image
      if (ABS(IMGPIXEL(corr_proj,i,j))<1)
         IMGPIXEL(corr_proj,i,j)=SGN(IMGPIXEL(corr_proj,i,j));
      IMGPIXEL(corr_proj,i,j)=
         applied_lambda*IMGPIXEL(diff_proj,i,j)/IMGPIXEL(corr_proj,i,j);
   }
   mean_error /= XSIZE(diff_proj())*YSIZE(diff_proj());
   
   // Denoising of the correction image
   if (vol_var!=NULL) process_correction(corr_proj,theo_proj,vol_var);

   // Backprojection of correction plane ......................................
   project_Volume(*vol_out,*footprint,*footprint2,theo_proj,
      corr_proj,YSIZE(read_proj()),XSIZE(read_proj()),
      read_proj.rot(),read_proj.tilt(),read_proj.psi(),BACKWARD,prm.eq_mode,
      prm.GVNeq,NULL,vol_var);

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
   GridVolume  vol_blobs, *vol_blobs_var=NULL;
   Basic_ROUT_Art(prm,eprm,vol_voxels, vol_blobs, vol_blobs_var);
}

/* Finish iterations ------------------------------------------------------- */
void finish_ART_iterations(const Basic_ART_Parameters &prm,
   const Plain_ART_Parameters &eprm, GridVolume &vol_blobs) {}

/* Apply_symmetry ------------------------------------------------------- */
void apply_symmetry(GridVolume &vol_in, GridVolume *vol_out,
                    const Plain_ART_Parameters &eprm, int grid_type)
    {cout <<"\nERROR: Function not implemented for single particles"<<endl;
     }
