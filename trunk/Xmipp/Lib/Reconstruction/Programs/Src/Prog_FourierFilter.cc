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

#ifdef _HAVE_VTK
   #include "../Prog_FourierFilter.hh"
   #include <XmippData/xmippArgs.hh>
   #include <XmippData/xmippMasks.hh>

/* Clear ------------------------------------------------------------------- */
void FourierMask::clear() {
   FilterShape=RAISED_COSINE;
   FilterBand=LOWPASS;
   w2=w1=0;
   raised_w=0;
   ctf.clear();
   ctf.enable_CTFnoise=TRUE;
   if (mask!=NULL) {mask->Delete(); mask=NULL;}
}

/* Assignment -------------------------------------------------------------- */
FourierMask & FourierMask::operator = (const FourierMask &F) {
   if (this!=&F) {
      clear();

      FilterShape=F.FilterShape;
      FilterBand=F.FilterBand;
      w2=F.w2;
      w1=F.w1;
      raised_w=F.raised_w;
      ctf=F.ctf;
      if (F.mask!=NULL) VTK2VTK(F.mask,mask);
   }
   return *this;
}

/* Read parameters from command line. -------------------------------------- */
void FourierMask::read(int argc, char **argv) _THROW {
   clear();

   // Filter shape .........................................................
   int i=position_param(argc,argv,"-fourier_mask");
   if (i==-1) return;
   if (i+1>=argc) REPORT_ERROR(3000,"FourierMask: -mask with no mask_type");
   if (strcmp(argv[i+1],"raised_cosine")==0) {
      if (i+2>=argc)
         REPORT_ERROR(3000,"FourierMask: Raised cosine needs a number of pixels");
      raised_w=AtoF(argv[i+2]);
      FilterShape=RAISED_COSINE;
   } else if (strcmp(argv[i+1],"ctf")==0) {
      if (i+2>=argc)
         REPORT_ERROR(3000,"FourierMask: CTF needs a CTF file");
      ctf.read(argv[i+2]);
      FilterShape=FilterBand=CTF;
      ctf.Produce_Side_Info();
   } else {
      if (i+1>=argc)
         REPORT_ERROR(3000,"FourierMask: you haven't supplied a file");
      fn_mask=argv[i+1];
      FilterShape=FilterBand=FROM_FILE;
   }

   // Filter band ..........................................................
   if (check_param(argc,argv,"-low_pass")) {
      w1=AtoF(get_param(argc,argv,"-low_pass")); FilterBand=LOWPASS;
   } else if (check_param(argc,argv,"-high_pass")) {
      w1=AtoF(get_param(argc,argv,"-high_pass")); FilterBand=HIGHPASS;
   } else if (check_param(argc,argv,"-band_pass")) {
      if (!get_2_double_params(argc,argv,"-band_pass",w1,w2,0,0))
         REPORT_ERROR(3000,"FourierMask: Not enough parameters after -band_pass");
      FilterBand=BANDPASS;
   } else if (check_param(argc,argv,"-stop_band")) {
      if (!get_2_double_params(argc,argv,"-stop_band",w1,w2,0,0))
         REPORT_ERROR(3000,"FourierMask: Not enough parameters after -stop_band");
      FilterBand=STOPBAND;
   }
   if (check_param(argc,argv,"-sampling")) {
      double sampling_rate=AtoF(get_param(argc,argv,"-sampling"));
      if (w1!=0)       w1=sampling_rate/w1;
      if (w2!=0)       w2=sampling_rate/w2;
      /*CO: I think it is more confusing */
      //if (raised_w!=0) raised_w=sampling_rate/raised_w;
   }
}

/* Show -------------------------------------------------------------------- */
void FourierMask::show() {
   cout << "Filter Band: ";
   switch (FilterBand) {
      case LOWPASS:  cout << "Lowpass before " << w1 << endl;  break;
      case HIGHPASS: cout << "Highpass after " << w1 << endl; break;
      case BANDPASS: cout << "Bandpass between " << w1 << " and " << w2 << endl;
         break;
      case STOPBAND: cout << "Stopband between " << w1 << " and " << w2 << endl;
         break;
      case CTF:      cout << "CTF\n";      break;
      case FROM_FILE: cout <<"From file " << fn_mask << endl; break;
   }
   cout << "Filter Shape: ";
   switch (FilterShape) {
      case RAISED_COSINE:
         cout << "Raised cosine with " << raised_w
              << " raised frequencies\n";
         break;
      case CTF:
         cout << "CTF\n" << ctf;
         break;
	  case FROM_FILE:
	     cout << "From file " << fn_mask << endl;
		 break;
   }
}

/* Usage ------------------------------------------------------------------- */
void FourierMask::usage() {
   cerr << "   -low_pass  <w1>                   : Cutoff freq (<1/2 or A)\n"
        << "   -high_pass <w1>                   : Cutoff freq (<1/2 or A)\n"
        << "   -band_pass <w1> <w2>              : Cutoff freq (<1/2 or A)\n"
        << "   -stop_band <w1> <w2>              : Cutoff freq (<1/2 or A)\n"
        << "   -fourier_mask <file>              : Provide a Fourier file\n"
        << "   -fourier_mask raised_cosine <raisedw>: Use raised cosine edges (in dig.freq.)\n"
        << "   -fourier_mask ctf                 : In that case the following\n"
	<< "                                       parameters apply\n"
        << "  [-sampling <sampling_rate>]        : If provided pass frequencies\n"
        << "                                       are taken in Angstroms\n"
   ;
   cerr << "CTF parameters -----------------------------\n";
   ctf.Usage();
}

/* Generate mask for a resized image --------------------------------------- */
void FourierMask::generate_mask(vtkImageData *v) _THROW {
   if (FilterBand!=FROM_FILE) {
      // Resize Xmipp real mask
      bool copy_from_Xmipp_real_mask=TRUE;
      Mask_Params real_mask;
      int dim[3]; v->GetDimensions(dim);
      double N1=w1*dim[0];
      double N2=w2*dim[0];
      double raised_pixels=raised_w*dim[0];

      // Generate mask
      switch (FilterBand) {
	 case LOWPASS:
            switch (FilterShape) {
	       case RAISED_COSINE:
		  real_mask.type=RAISED_COSINE_MASK;
		  real_mask.mode=INNER_MASK;
		  real_mask.R1=N1;
		  real_mask.R2=N1+raised_pixels;
		  real_mask.x0=real_mask.y0=real_mask.z0=0;
		  break;
	    }
            break;
	 case HIGHPASS:
            switch (FilterShape) {
	       case RAISED_COSINE:
		  real_mask.type=RAISED_COSINE_MASK;
		  real_mask.mode=OUTSIDE_MASK;
		  real_mask.R1=N1-raised_pixels;
		  real_mask.R2=N1;
		  real_mask.x0=real_mask.y0=real_mask.z0=0;
		  break;
	    }
            break;
	 case BANDPASS:
            switch (FilterShape) {
	       case RAISED_COSINE:
		  real_mask.type=RAISED_CROWN_MASK;
		  real_mask.mode=INNER_MASK;
		  real_mask.R1=N1;
		  real_mask.R2=N2;
		  real_mask.x0=real_mask.y0=real_mask.z0=0;
		  real_mask.pix_width=raised_pixels;
		  break;
	    }
            break;
	 case STOPBAND:
            switch (FilterShape) {
	       case RAISED_COSINE:
		  real_mask.type=RAISED_CROWN_MASK;
		  real_mask.mode=OUTSIDE_MASK;
		  real_mask.R1=N1;
		  real_mask.R2=N2;
		  real_mask.x0=real_mask.y0=real_mask.z0=0;
		  real_mask.pix_width=raised_pixels;
		  break;
	    }
            break;
	 case CTF:
            copy_from_Xmipp_real_mask=FALSE;
	    if (dim[2]==1) ctf.Generate_CTF(v,mask);
	    else REPORT_ERROR(1,
	       "Fourier_mask::generate: Cannot apply CTF to a volume\n");
            break;
      }

      // Copy mask from real Xmipp mask
      if (copy_from_Xmipp_real_mask) {
      	 if (dim[2]==1 && dim[1]==1) {
            real_mask.resize(dim[0]);
            real_mask.generate_1Dmask();
	    // CO: for even sizes the CenterFFT does not work too well
	    if (dim[0]%2==0) shift_for_VTK(real_mask.get_cont_mask1D());
	    xmippArray2VTK(real_mask.get_cont_mask1D(),mask,2);
	    CenterFFT(mask);
	 } else if (dim[2]==1) {
            real_mask.resize(dim[1],dim[0]);
            real_mask.generate_2Dmask();
	    // CO: for even sizes the CenterFFT does not work too well
	    if (dim[1]%2==0) shift_for_VTK(real_mask.get_cont_mask2D(),'y');
	    if (dim[0]%2==0) shift_for_VTK(real_mask.get_cont_mask2D(),'x');
	    xmippArray2VTK(real_mask.get_cont_mask2D(),mask,2);
	    CenterFFT(mask);
	 } else {
            real_mask.resize(dim[2],dim[1],dim[0]);
            real_mask.generate_3Dmask();
	    // CO: for even sizes the CenterFFT does not work too well
	    if (dim[2]%2==0) shift_for_VTK(real_mask.get_cont_mask3D(),'z');
	    if (dim[1]%2==0) shift_for_VTK(real_mask.get_cont_mask3D(),'y');
	    if (dim[0]%2==0) shift_for_VTK(real_mask.get_cont_mask3D(),'x');
	    xmippArray2VTK(real_mask.get_cont_mask3D(),mask,2);
	    CenterFFT(mask);
	 }
      }
   } else {
   // Read from file
      if (Is_FourierImageXmipp(fn_mask)) {
         FourierImageXmipp  I; I.read(fn_mask); xmippFFT2VTK(I,mask);
      } else {
         FourierVolumeXmipp V; V.read(fn_mask); xmippFFT2VTK(V,mask);
      }
   }
   mask->UpdateInformation();
}

/* Save -------------------------------------------------------------------- */
void FourierMask::write_amplitude(const FileName &fn, bool do_not_center) {
   if (mask==NULL) return;
   int dim[3]; mask->GetDimensions(dim);
   vtkImageData *aux=NULL; VTK2VTK(mask,aux);
   if (!do_not_center) CenterFFT(aux);
   if (dim[2]==1 && dim[1]==1) {
      matrix1D<double> v; FFT_magnitude(aux,v);
      // CO: for even sizes the CenterFFT does not work too well
      // Even for odd sizes??
      // if (dim[1]%2==0 || dim[0]%2==0) {
      if (TRUE) shift_for_VTK(v);
      FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
         v(i)=log10(1+v(i)*v(i));
      v.write(fn);
   } else if (dim[2]==1) {
      ImageXmipp  I; FFT_magnitude(aux,I());
      // CO: for even sizes the CenterFFT does not work too well
      // Even for odd sizes??
      // if (dim[1]%2==0 || dim[0]%2==0) {
      if (TRUE) {
         shift_for_VTK(I(),'x');
         shift_for_VTK(I(),'y');
      }
      FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
         I(i,j)=log10(1+I(i,j)*I(i,j));
      I.write(fn);
   } else {
      VolumeXmipp V; FFT_magnitude(aux,V());
      // CO: for even sizes the CenterFFT does not work too well
      // Even for odd sizes??
      // if (dim[2]%2==0 || dim[1]%2==0 || dim[0]%2==0) {
      if (TRUE) {
         shift_for_VTK(V(),'x');
         shift_for_VTK(V(),'y');
         shift_for_VTK(V(),'z');
      }
      FOR_ALL_ELEMENTS_IN_MATRIX3D(V())
         V(k,i,j)=log10(1+V(k,i,j)*V(k,i,j));
      V.write(fn);
   }
}

void FourierMask::write_mask(const FileName &fn) {
   if (mask==NULL) return;
   int dim[3]; mask->GetDimensions(dim);
   if (dim[2]==1 && dim[1]==1) {
      matrix1D<double_complex> v; VTK2xmippFFT(mask,v); v.write(fn);
   } else if (dim[2]==1) {
      FourierImageXmipp  I; VTK2xmippFFT(mask,I); I.write(fn);
   } else {
      FourierVolumeXmipp V; VTK2xmippFFT(mask,V); V.write(fn);
   }
}

/* Apply mask -------------------------------------------------------------- */
void FourierMask::apply_mask(vtkImageData *v) _THROW {
   if (!same_shape(v,mask))
      REPORT_ERROR(1,"FourierMask::apply_mask: mask and input FT do not have"
         " the same shape");
   if (v->GetScalarType()!=VTK_FLOAT || v->GetNumberOfScalarComponents()!=2)
      REPORT_ERROR(1,"FourierMask::apply_mask: input array do not seem to be"
         " a Fourier Transform");
   
   int dim[3]; v->GetDimensions(dim);
   float *mi=(float *) mask->GetScalarPointer();
   float *vi=(float *)    v->GetScalarPointer();
   for (int k=0; k<dim[2]; k++)
       for (int i=0; i<dim[1]; i++)
      	   for (int j=0; j<dim[0]; j++) {
	      double result_re=(*mi)*(*vi)-(*(mi+1))*(*(vi+1));
	      double result_im=(*(mi+1))*(*vi)+(*(mi))*(*(vi+1));
	      *(vi)=result_re; *(vi+1)=result_im;
	      vi+=2; mi+=2;
	   }
}

void FourierMask::apply_mask(matrix1D<double> &v) {
   vtkImageData *FFTv=NULL;
   int startingX=STARTINGX(v);
   FFT_VTK(v,FFTv,TRUE);
   if (mask==NULL) generate_mask(FFTv);
   apply_mask(FFTv);
   IFFT_VTK(FFTv,v,TRUE);
   STARTINGX(v)=startingX;
   FFTv->Delete();
}

void FourierMask::apply_mask(matrix2D<double> &v) {
   vtkImageData *FFTv=NULL;
   int startingX=STARTINGX(v);
   int startingY=STARTINGY(v);
   FFT_VTK(v,FFTv,TRUE);
   if (mask==NULL) generate_mask(FFTv);
   apply_mask(FFTv);
   IFFT_VTK(FFTv,v,TRUE);
   STARTINGX(v)=startingX;
   STARTINGY(v)=startingY;
   FFTv->Delete();
}

void FourierMask::apply_mask(matrix3D<double> &v) {
   vtkImageData *FFTv=NULL;
   int startingX=STARTINGX(v);
   int startingY=STARTINGY(v);
   int startingZ=STARTINGZ(v);
   FFT_VTK(v,FFTv,TRUE);
   if (mask==NULL) generate_mask(FFTv);
   apply_mask(FFTv);
   IFFT_VTK(FFTv,v,TRUE);
   STARTINGX(v)=startingX;
   STARTINGY(v)=startingY;
   STARTINGZ(v)=startingZ;
   FFTv->Delete();
}

/* Resize ------------------------------------------------------------------ */
void FourierMask::resize_mask(int Ydim, int Xdim) {
   if (mask==NULL) return;
   FourierImageXmipp I;
   VTK2xmippFFT(mask,I);
   I().scale_to_size(Ydim,Xdim);
   xmippFFT2VTK(I,mask);
   mask->UpdateInformation();
}

/* Get size ---------------------------------------------------------------- */
void FourierMask::mask_size(int &Zdim, int &Ydim, int &Xdim) {
   int dim[3]; mask->GetDimensions(dim);
   Zdim=dim[2];
   Ydim=dim[1];
   Xdim=dim[0];
}

/* Mask power -------------------------------------------------------------- */
double FourierMask::mask_power(double wmin, double wmax) {
   if (mask==NULL) return 0;
   matrix1D<int>    idx(2);
   matrix1D<double> freq(2);
   double retval=0, N=0;
   SPEED_UP_vtk;
   FOR_ALL_ELEMENTS_IN_VTK(mask) {
      XX(idx)=j; YY(idx)=i;
      FFT_idx2digfreq(mask, idx, freq);
      double w=freq.module();
      if (w>wmin && w<wmax) {
         retval += (*vtkPtr)*(*vtkPtr)+(*(vtkPtr+1))*(*(vtkPtr+1));
         N++;
      }
      vtkPtr+=2;
   }

   if (N!=0) return retval/N;
   else      return 0;
}
#endif
