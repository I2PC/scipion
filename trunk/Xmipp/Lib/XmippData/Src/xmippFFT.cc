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

#include "../xmippFFT.hh"
#include "../xmippArgs.hh"
#include "../Bilib/headers/dft.h"

/* Format conversions ------------------------------------------------------ */
template <class T>
void RealImag2Complex(const T *_real, const T *_imag,
   complex<double> *_complex, int length) {
   T *aux_real   =(T *)_real;
   T *aux_imag   =(T *)_imag;
   double *aux_complex=(double *)_complex;
   for (int i=0; i<length; i++) {
      *aux_complex++=(double)(*aux_real++);
      *aux_complex++=(double)(*aux_imag++);
   }
}

template <class T>
void AmplPhase2Complex(const T *_ampl, const T *_phase,
   complex<double> *_complex, int length) {
   T *aux_ampl   =(T *)_ampl;
   T *aux_phase  =(T *)_phase;
   double *aux_complex=(double *)_complex;
   for (int i=0; i<length; i++) {
      double ampl =(double)(*aux_ampl++);
      double phase=(double)(*aux_phase++);
      *aux_complex++=ampl*cos(phase);
      *aux_complex++=ampl*sin(phase);
   }
}

template <class T>
void Complex2RealImag(const complex<double> *_complex,
   T *_real, T *_imag, int length) {
   T *aux_real   =(T *)_real;
   T *aux_imag   =(T *)_imag;
   double *aux_complex=(double *)_complex;
   for (int i=0; i<length; i++) {
      *aux_real++=(T)(*aux_complex++);
      *aux_imag++=(T)(*aux_complex++);
   }
}

template <class T>
void Complex2AmplPhase(const complex<double> *_complex,
   T *_ampl, T *_phase, int length) {
   T *aux_ampl   =(T *)_ampl;
   T *aux_phase  =(T *)_phase;
   double *aux_complex=(double *)_complex;
   for (int i=0; i<length; i++) {
      double re=*aux_complex++;
      double im=*aux_complex++;
      *aux_ampl++=sqrt(re*re+im*im);
      *aux_phase++=atan2(im,re);
   }
}


/** Convert whole -> half of (centro-symmetric) Fourier transforms 2D. -- */
void Whole2Half(const matrix2D<complex<double> > &in, matrix2D<complex<double> > &out) {

  int ydim=(int)(YSIZE(in)/2)+1;

  out.resize(ydim,XSIZE(in));
  for (int i=0; i<ydim; i++) 
    for (int j=0; j<XSIZE(in); j++)
      dMij(out,i,j)=dMij(in,i,j);

}

/** Convert half -> whole of (centro-symmetric) Fourier transforms 2D. -- */
void Half2Whole(const matrix2D<complex<double> > &in, matrix2D<complex<double> > &out, int oriydim) {

  int yshift=2*(YSIZE(in)-1);
  out.resize(oriydim,XSIZE(in));

  // Old part
  for (int i=0; i<YSIZE(in); i++) 
    for (int j=0; j<XSIZE(in); j++)
      dMij(out,i,j)=dMij(in,i,j);
  // New part
  for (int i=YSIZE(in); i<oriydim; i++) {
    dMij(out,i,0)=conj(dMij(in,oriydim-i,0));
    for (int j=1; j<XSIZE(in); j++)
      dMij(out,i,j)=conj(dMij(in,oriydim-i,XSIZE(in)-j));
  }    
}

/** Direct Fourier Transform 1D ------------------------------------------- */
void FourierTransform(const matrix1D<double> &in,
   matrix1D< complex<double> > &out) {
   int N=XSIZE(in);
   matrix1D<double> re(in), tmp(N), im(N), cas(N);
   out.resize(N);
   
   GetCaS(MULTIDIM_ARRAY(cas),N);
   DftRealToRealImaginary (MULTIDIM_ARRAY(re),MULTIDIM_ARRAY(im),
      MULTIDIM_ARRAY(tmp), MULTIDIM_ARRAY(cas), N);
   RealImag2Complex(MULTIDIM_ARRAY(re),MULTIDIM_ARRAY(im),
      MULTIDIM_ARRAY(out),N);
}

/** Direct Fourier Transform 2D. ------------------------------------------ */
void FourierTransform(const matrix2D<double> &in,
   matrix2D< complex<double> > &out) {
   int Status;
   matrix2D<double> re(in), im;
   im.resize(in);
   out.resize(in);
   VolumeDftRealToRealImaginary(MULTIDIM_ARRAY(re),
      MULTIDIM_ARRAY(im), XSIZE(in), YSIZE(in), 1, &Status); 
   RealImag2Complex(MULTIDIM_ARRAY(re),MULTIDIM_ARRAY(im),
      MULTIDIM_ARRAY(out),XSIZE(in)*YSIZE(in));
}

/** Direct Fourier Transform 3D. ------------------------------------------ */
void FourierTransform(const matrix3D<double> &in,
   matrix3D< complex<double> > &out) {
   int Status;
   matrix3D<double> re(in), im;
   im.resize(in);
   out.resize(in);
   VolumeDftRealToRealImaginary(MULTIDIM_ARRAY(re),
      MULTIDIM_ARRAY(im), XSIZE(in), YSIZE(in), ZSIZE(in), &Status); 
   RealImag2Complex(MULTIDIM_ARRAY(re),MULTIDIM_ARRAY(im),
      MULTIDIM_ARRAY(out),XSIZE(in)*YSIZE(in)*ZSIZE(in));
}

/** Inverse Fourier Transform 1D. ----------------------------------------- */
void InverseFourierTransform(const matrix1D< complex<double> > &in,
   matrix1D<double> &out) {
   int N=XSIZE(in);
   matrix1D<double> tmp(N), im(N), cas(N);
   out.resize(N);
   
   GetCaS(MULTIDIM_ARRAY(cas),N);
   Complex2RealImag(MULTIDIM_ARRAY(in),MULTIDIM_ARRAY(out),
      MULTIDIM_ARRAY(im),N);
   InvDftRealImaginaryToReal(MULTIDIM_ARRAY(out),MULTIDIM_ARRAY(im),
      MULTIDIM_ARRAY(tmp), MULTIDIM_ARRAY(cas), N);
}

/** Inverse Fourier Transform 2D. ----------------------------------------- */
void InverseFourierTransform(const matrix2D< complex<double> > &in,
   matrix2D<double> &out) {
   int Status;
   matrix2D<double> im;
   out.resize(in);
   im.resize(in);
   Complex2RealImag(MULTIDIM_ARRAY(in),MULTIDIM_ARRAY(out),
      MULTIDIM_ARRAY(im),XSIZE(in)*YSIZE(in));
   VolumeInvDftRealImaginaryToReal(MULTIDIM_ARRAY(out),
      MULTIDIM_ARRAY(im), XSIZE(in), YSIZE(in), 1, &Status);
}

/** Inverse Fourier Transform 3D. ----------------------------------------- */
void InverseFourierTransform(const matrix3D< complex<double> > &in,
   matrix3D<double> &out) {
   int Status;
   matrix3D<double> im;
   out.resize(in);
   im.resize(in);
   Complex2RealImag(MULTIDIM_ARRAY(in),MULTIDIM_ARRAY(out),
      MULTIDIM_ARRAY(im),XSIZE(in)*YSIZE(in)*ZSIZE(in));
   VolumeInvDftRealImaginaryToReal(MULTIDIM_ARRAY(out),
      MULTIDIM_ARRAY(im), XSIZE(in), YSIZE(in), ZSIZE(in), &Status); 
}


/** Direct Fourier Transform 2D, output half of (centro-symmetric) transform ---- */
void FourierTransformHalf(const matrix2D<double> &in,
   matrix2D< complex<double> > &out) {

  matrix2D<complex <double> > aux;
  FourierTransform(in,aux);
  Whole2Half(aux,out);
}

/** Inverse Fourier Transform 2D, input half of (centro-symmetric) transform ---- */
void InverseFourierTransformHalf(const matrix2D< complex<double> > &in,
   matrix2D<double> &out, int oriydim) {

  matrix2D< complex<double> > aux;
  Half2Whole(in,aux,oriydim);
  InverseFourierTransform(aux,out);
  out.set_Xmipp_origin();
}

/* CenterFFT 1D. ----------------------------------------------------------- */
template <class T>
void CenterFFT(matrix1D<T> &v, bool forward) {
   matrix1D<T> aux;
   int l, shift;
   
   l=XSIZE(v); aux.resize(l); shift=(int)(l/2);
   if (!forward) shift=-shift;
   // Shift the input in an auxiliar vector
   for (int i=0; i<l; i++) {
      int ip=i+shift; if (ip<0) ip+=l; else if (ip>=l) ip-=l;
      aux(ip)=DIRECT_VEC_ELEM(v,i);
   }
   // Copy the vector
   for (int i=0; i<l; i++)
      DIRECT_VEC_ELEM(v,i)=DIRECT_VEC_ELEM(aux,i);
}

/* CenterFFT 2D. ----------------------------------------------------------- */
template <class T>
void CenterFFT(matrix2D<T> &v, bool forward) {
   matrix1D<T> aux;
   int l, shift;

   // Shift in the X direction
   l=XSIZE(v); aux.resize(l); shift=(int)(l/2);
   if (!forward) shift=-shift;
   for (int i=0; i<YSIZE(v); i++) {
      // Shift the input in an auxiliar vector
      for (int j=0; j<l; j++) {
         int jp=j+shift; if (jp<0) jp+=l; else if (jp>=l) jp-=l;
         aux(jp)=DIRECT_MAT_ELEM(v,i,j);
      }
      // Copy the vector
      for (int j=0; j<l; j++)
         DIRECT_MAT_ELEM(v,i,j)=DIRECT_VEC_ELEM(aux,j);
   }

   // Shift in the Y direction
   l=YSIZE(v); aux.resize(l); shift=(int)(l/2);
   if (!forward) shift=-shift;
   for (int j=0; j<XSIZE(v); j++) {
      // Shift the input in an auxiliar vector
      for (int i=0; i<l; i++) {
         int ip=i+shift; if (ip<0) ip+=l; else if (ip>=l) ip-=l;
         aux(ip)=DIRECT_MAT_ELEM(v,i,j);
      }
      // Copy the vector
      for (int i=0; i<l; i++)
         DIRECT_MAT_ELEM(v,i,j)=DIRECT_VEC_ELEM(aux,i);
   }
}

/* CenterFFT 3D. ----------------------------------------------------------- */
template <class T>
void CenterFFT(matrix3D<T> &v, bool forward) {
   matrix1D<T> aux;
   int l, shift;

   // Shift in the X direction
   l=XSIZE(v); aux.resize(l); shift=(int)(l/2);
   if (!forward) shift=-shift;
   for (int k=0; k<ZSIZE(v); k++)
      for (int i=0; i<YSIZE(v); i++) {
         // Shift the input in an auxiliar vector
         for (int j=0; j<l; j++) {
            int jp=j+shift; if (jp<0) jp+=l; else if (jp>=l) jp-=l;
            aux(jp)=DIRECT_VOL_ELEM(v,k,i,j);
         }
         // Copy the vector
         for (int j=0; j<l; j++)
            DIRECT_VOL_ELEM(v,k,i,j)=DIRECT_VEC_ELEM(aux,j);
      }

   // Shift in the Y direction
   l=YSIZE(v); aux.resize(l); shift=(int)(l/2);
   if (!forward) shift=-shift;
   for (int k=0; k<ZSIZE(v); k++)
      for (int j=0; j<XSIZE(v); j++) {
         // Shift the input in an auxiliar vector
         for (int i=0; i<l; i++) {
            int ip=i+shift; if (ip<0) ip+=l; else if (ip>=l) ip-=l;
            aux(ip)=DIRECT_VOL_ELEM(v,k,i,j);
         }
         // Copy the vector
         for (int i=0; i<l; i++)
            DIRECT_VOL_ELEM(v,k,i,j)=DIRECT_VEC_ELEM(aux,i);
      }

   // Shift in the Z direction
   l=ZSIZE(v); aux.resize(l); shift=(int)(l/2);
   if (!forward) shift=-shift;
   for (int i=0; i<YSIZE(v); i++)
      for (int j=0; j<XSIZE(v); j++) {
         // Shift the input in an auxiliar vector
         for (int k=0; k<l; k++) {
            int kp=k+shift; if (kp<0) kp+=l; else if (kp>=l) kp-=l;
            aux(kp)=DIRECT_VOL_ELEM(v,k,i,j);
         }
         // Copy the vector
         for (int k=0; k<l; k++)
            DIRECT_VOL_ELEM(v,k,i,j)=DIRECT_VEC_ELEM(aux,k);
      }
}

/* FFT Magnitude 1D. ------------------------------------------------------- */
void FFT_magnitude(const matrix1D< complex<double> > &v,
   matrix1D<double> &mag) {
   mag.resize(v);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(v) mag(i)=abs(v(i));
}

void FFT_magnitude(const matrix2D< complex<double> > &v,
   matrix2D<double> &mag) {
   mag.resize(v);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(v) mag(i,j)=abs(v(i,j));
}

void FFT_magnitude(const matrix3D< complex<double> > &v,
   matrix3D<double> &mag) {
   mag.resize(v);
   FOR_ALL_ELEMENTS_IN_MATRIX3D(v) mag(k,i,j)=abs(v(k,i,j));
}

/* FFT Phase 1D. ------------------------------------------------------- */
void FFT_phase(const matrix1D< complex<double> > &v,
   matrix1D<double> &phase) {
   phase.resize(v);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(v) phase(i)=atan2(v(i).imag(),v(i).real());
}

void FFT_phase(const matrix2D< complex<double> > &v,
   matrix2D<double> &phase) {
   phase.resize(v);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(v) phase(i,j)=atan2(v(i,j).imag(),v(i,j).real());
}

void FFT_phase(const matrix3D< complex<double> > &v,
   matrix3D<double> &phase) {
   phase.resize(v);
   FOR_ALL_ELEMENTS_IN_MATRIX3D(v) phase(k,i,j)=atan2(v(k,i,j).imag(),
      v(k,i,j).real());
}

/* Correlation and Autocorrelation ----------------------------------------- */
template <class T>
void auto_correlation_matrix(const matrix2D<T> &Img, matrix2D<double> &R) {
   // Compute the Fourier Transform
   matrix2D< complex<double> > FFT1; FourierTransform(Img, FFT1);

   // Multiply FFT1 * FFT1'
   FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT1)
      FFT1(i,j)*=conj(FFT1(i,j));

   // Invert the product, in order to obtain the correlation image
   InverseFourierTransform(FFT1,R);
	  
   // Center the resulting image to obtain a centered autocorrelation
   CenterFFT(R,true);
}

template <class T>
void correlation_matrix(const matrix2D<T> &m1, const matrix2D<T> &m2,
   matrix2D<double> &R) {
   // Compute the Fourier Transforms
   matrix2D< complex<double> > FFT1; FourierTransform(m1, FFT1);
   matrix2D< complex<double> > FFT2; FourierTransform(m2, FFT2);

   // Multiply FFT1 * FFT2'
   FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT1)
      FFT1(i,j)*=conj(FFT2(i,j));

   // Invert the product, in order to obtain the correlation image
   InverseFourierTransform(FFT1,R);
	  
   // Center the resulting image to obtain a centered autocorrelation
   CenterFFT(R,true);
}

/* 2D Fourier Ring Correlation ------------------------------------------------ */
template <class T>
void fourier_ring_correlation(matrix2D<T> const &m1, matrix2D<T> const &m2, double sam,
		 matrix1D<double> &freq, matrix1D<double> &frc, matrix1D<double> &frc_noise) {

  if (!m1.same_shape(m2)) {
    cerr << "Error: matrices have different shapes!"<< endl;
    exit(0);
  }

  matrix2D<T> aux(m1);
  matrix1D<int> origin(3),radial_count;
  matrix1D<double> tmp1,tmp2;
  matrix1D<complex <double> > tmp3;
  matrix2D<complex <double> > FT1; FourierTransform(m1,FT1); CenterFFT(FT1,true);
  matrix2D<complex <double> > FT2; FourierTransform(m2,FT2); CenterFFT(FT2,true);
  int dim=(int)FT1.RowNo()/2;
  origin.init_zeros();

  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux) {
    dMij(aux,i,j)=abs(dMij(FT1,i,j))*abs(dMij(FT1,i,j));
  }
  tmp1.init_zeros();
  radial_average(aux,origin,tmp1,radial_count,TRUE);
  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux) {
    dMij(aux,i,j)=abs(dMij(FT2,i,j))*abs(dMij(FT2,i,j));
  }
  tmp2.init_zeros();
  radial_average(aux,origin,tmp2,radial_count,TRUE);
  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(FT1) {
    dMij(FT1,i,j)=conj(dMij(FT1,i,j))*dMij(FT2,i,j);
  }
  tmp3.init_zeros();
  radial_average(FT1,origin,tmp3,radial_count,TRUE);
  FFT_magnitude(tmp3,frc);
  frc.resize(dim);
  frc_noise.resize(dim);
  freq.resize(dim);
  FOR_ALL_ELEMENTS_IN_MATRIX1D(freq) {
    int j=i;
    VEC_ELEM(freq,i)=(double)j/(dim*2*sam);
    VEC_ELEM(frc,i)=VEC_ELEM(frc,i)/sqrt(VEC_ELEM(tmp1,i)*VEC_ELEM(tmp2,i));
    VEC_ELEM(frc_noise,i)=2/sqrt((double)VEC_ELEM(radial_count,i));
  }

}

/* 3D Fourier Ring Correlation ------------------------------------------------ */
template <class T>
void fourier_ring_correlation(matrix3D<T> const &m1, matrix3D<T> const &m2, double sam,
		 matrix1D<double> &freq, matrix1D<double> &frc, matrix1D<double> &frc_noise) {

  if (!m1.same_shape(m2)) {
    cerr << "Error: matrices have different shapes!"<< endl;
    exit(0);
  }

  matrix3D<T> aux(m1);
  matrix1D<int> origin(3),radial_count;
  matrix1D<double> tmp1,tmp2;
  matrix1D<complex <double> > tmp3;
  matrix3D<complex <double> > FT1; FourierTransform(m1,FT1); CenterFFT(FT1,true);
  matrix3D<complex <double> > FT2; FourierTransform(m2,FT2); CenterFFT(FT2,true);
  int dim=(int)FT1.RowNo()/2;
  origin.init_zeros();

  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(aux) {
    dVkij(aux,k,i,j)=abs(dVkij(FT1,k,i,j))*abs(dVkij(FT1,k,i,j));
  }
  tmp1.init_zeros();
  radial_average(aux,origin,tmp1,radial_count,TRUE);
  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(aux) {
    dVkij(aux,k,i,j)=abs(dVkij(FT2,k,i,j))*abs(dVkij(FT2,k,i,j));
  }
  tmp2.init_zeros();
  radial_average(aux,origin,tmp2,radial_count,TRUE);
  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FT1) {
    dVkij(FT1,k,i,j)=conj(dVkij(FT1,k,i,j))*dVkij(FT2,k,i,j);
  }
  tmp3.init_zeros();
  radial_average(FT1,origin,tmp3,radial_count,TRUE);
  FFT_magnitude(tmp3,frc);
  frc.resize(dim);
  frc_noise.resize(dim);
  freq.resize(dim);
  FOR_ALL_ELEMENTS_IN_MATRIX1D(freq) {
    int j=i;
    VEC_ELEM(freq,i)=(double)j/(dim*2*sam);
    VEC_ELEM(frc,i)=VEC_ELEM(frc,i)/sqrt(VEC_ELEM(tmp1,i)*VEC_ELEM(tmp2,i));
    VEC_ELEM(frc_noise,i)=2/sqrt((double)VEC_ELEM(radial_count,i));
  }

}

/* 2D Differential Phase residual ------------------------------------------------ */
template <class T>
void differential_phase_residual(matrix2D<T> const &m1, matrix2D<T> const &m2, double sam,
		 matrix1D<double> &freq, matrix1D<double> &dpr) {

  if (!m1.same_shape(m2)) {
    cerr << "Error: matrices have different shapes!"<< endl;
    exit(0);
  }

  matrix2D<T> aux(m1);
  matrix1D<int> origin(3),radial_count;
  matrix1D<double> tmp1,tmp2;
  matrix1D<complex <double> > tmp3;
  matrix2D<complex <double> > FT1; FourierTransform(m1,FT1); CenterFFT(FT1,true);
  matrix2D<complex <double> > FT2; FourierTransform(m2,FT2); CenterFFT(FT2,true);
  int dim=(int)FT1.RowNo()/2;

  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux) {
    dMij(aux,i,j)=abs(dMij(FT1,i,j))+abs(dMij(FT2,i,j));
  }
  tmp1.init_zeros();
  radial_average(aux,origin,tmp1,radial_count,TRUE);
  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux) {
    dMij(aux,i,j)*=realWRAP(RAD2DEG(
                     (atan2(dMij(FT1,i,j).imag(),dMij(FT1,i,j).real())) -
                     (atan2(dMij(FT2,i,j).imag(),dMij(FT2,i,j).real()))) 
				  ,-180,180);
    dMij(aux,i,j)*=realWRAP(RAD2DEG(
                     (atan2(dMij(FT1,i,j).imag(),dMij(FT1,i,j).real())) -
                     (atan2(dMij(FT2,i,j).imag(),dMij(FT2,i,j).real()))) 
				  ,-180,180);
  }
  tmp2.init_zeros();
  radial_average(aux,origin,tmp2,radial_count,TRUE);
  dpr=tmp2/tmp1;
  dpr=SQRTnD(dpr);
  dpr.resize(dim);
  freq.resize(dim);
  FOR_ALL_ELEMENTS_IN_MATRIX1D(freq) {
    int j=i;
    VEC_ELEM(freq,i)=(double)j/(dim*2*sam);
  }

}

/* 3D Differential Phase residual ------------------------------------------------ */
template <class T>
void differential_phase_residual(matrix3D<T> const &m1, matrix3D<T> const &m2, double sam,
		 matrix1D<double> &freq, matrix1D<double> &dpr) {

  if (!m1.same_shape(m2)) {
    cerr << "Error: matrices have different shapes!"<< endl;
    exit(0);
  }
  matrix3D<T> aux(m1);
  matrix1D<int> origin(3),radial_count;
  matrix1D<double> tmp1,tmp2;
  matrix1D<complex <double> > tmp3;
  matrix3D<complex <double> > FT1; FourierTransform(m1,FT1); CenterFFT(FT1,true);
  matrix3D<complex <double> > FT2; FourierTransform(m2,FT2); CenterFFT(FT2,true);
  int dim=(int)FT1.RowNo()/2;

  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(aux) {
    dVkij(aux,k,i,j)=abs(dVkij(FT1,k,i,j))+abs(dVkij(FT2,k,i,j));
  }
  tmp1.init_zeros();
  radial_average(aux,origin,tmp1,radial_count,TRUE);
  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(aux) {
    dVkij(aux,k,i,j)*=realWRAP(RAD2DEG(
                     (atan2(dVkij(FT1,k,i,j).imag(),dVkij(FT1,k,i,j).real())) -
                     (atan2(dVkij(FT2,k,i,j).imag(),dVkij(FT2,k,i,j).real()))) 
				  ,-180,180);
    dVkij(aux,k,i,j)*=realWRAP(RAD2DEG(
                     (atan2(dVkij(FT1,k,i,j).imag(),dVkij(FT1,k,i,j).real())) -
                     (atan2(dVkij(FT2,k,i,j).imag(),dVkij(FT2,k,i,j).real()))) 
				  ,-180,180);
  }
  tmp2.init_zeros();
  radial_average(aux,origin,tmp2,radial_count,TRUE);
  dpr=tmp2/tmp1;
  dpr=SQRTnD(dpr);
  dpr.resize(dim);
  freq.resize(dim);
  FOR_ALL_ELEMENTS_IN_MATRIX1D(freq) {
    int j=i;
    VEC_ELEM(freq,i)=(double)j/(dim*2*sam);
  }

}

/* 2D SSNR ------------------------------------------------ */
#define NGRUPOS 20
template <class T>
void my_ssnr(matrix2D<T> const &AverageImage, SelFile   &SF_sel, double sam,
		 matrix1D<double> &freq, matrix1D<double> &ssnr, 
		 matrix1D<double> &pixel,
		 bool apply_geo) {
  Image       Iaverage_sub_group,Id;
  double      dummy;
  SelFile  SF_tmp, SF_vector[NGRUPOS];
  SelLine  line;
  int N=NGRUPOS;
  
  if(SF_sel.ImgNo(SelLine::ACTIVE) < NGRUPOS)
      REPORT_ERROR(1,(string)"FFT::my_ssnr: my_ssnr I need more than " +
                   ItoA(SF_sel.ImgNo(SelLine::ACTIVE)) +
                    " images in file " + SF_sel.name());
                                                                                
  int contImg=SF_sel.ImgNo();
//  SF_tmp=SF_sel.randomize();
SF_tmp=SF_sel;
  int Nsub=(int)((double)contImg/N);//floor()
  //assing images
  for (int i=0; i<N; i++) {
    SF_vector[i].clear();
    SF_tmp.go_beginning();
    SF_tmp.jump_lines(Nsub*i);
    for (int nn=0; nn<Nsub; nn++) {
      SF_vector[i].insert(SF_tmp.current());
      SF_tmp.NextImg();
    }
  }
  //all the images has been assigned?
  if (contImg%N!=0){
     for(int i=(N-(contImg%N));i<N;i++){
	 SF_vector[i].insert(SF_tmp.current());
	 SF_tmp.NextImg();}

  }
//#define DEBUG//check image name subsets
#ifdef DEBUG
  FileName fn_sel,fn_out;
  fn_sel=SF_sel.name();
  for (int i=0; i<N; i++) {
     string num="_"+ItoA(i+1);
     fn_out=fn_sel.insert_before_extension(num);
     SF_vector[i].write(fn_out);
  }
#endif
#undef DEBUG
  //alloc memory for results      
  int dim=(int)AverageImage.RowNo()/2;
  ssnr.resize(dim);
  pixel.resize(dim);
  freq.resize(dim); 
  //compute average and fft
  //fft radial averaged image
  matrix2D<complex <double> > FTAverageImage,
                              FTaverageSubGroup;
  FourierTransform(AverageImage,FTAverageImage);
  CenterFFT(FTAverageImage,true);
  for (int i=0; i<N; i++) {
     SF_vector[i].get_statistics(Iaverage_sub_group,Id,dummy,dummy,apply_geo);
     FourierTransform(Iaverage_sub_group(),FTaverageSubGroup);
     CenterFFT(FTaverageSubGroup,true);   
     //calculate ssNR
     my_ssnr_step(FTAverageImage,FTaverageSubGroup,
                  ssnr,pixel,SF_vector[i].ImgNo()); 
	  
//#define DEBUG//check average images
#ifdef DEBUG
     {
     ImageXmipp IX; 
     IX=Iaverage_sub_group; 
     FileName fn_sel,fn_out;
     fn_sel=SF_sel.name();
     fn_out.compose(fn_sel.get_root(), i+1, (string) "xmp");
     IX.write(fn_out);
     }
#endif
#undef DEBUG      
   }//end for (int i=0; i<N; i++) {
   int dim2=(int)AverageImage.RowNo();
   FOR_ALL_ELEMENTS_IN_MATRIX1D(ssnr)
       {
       VEC_ELEM(ssnr,i)=(contImg*(NGRUPOS-1))/(VEC_ELEM(ssnr,i)-1);
       VEC_ELEM(freq,i)=((double)i+0.5)/(dim2*sam);
       if(VEC_ELEM(ssnr,i)>NGRUPOS)
          VEC_ELEM(ssnr,i) = NGRUPOS;
       }

}
/****************************************************************************/
/*            Square distance from the point (x,y) to (m/2,m/2)             */
/****************************************************************************/
 double distancia2( int x, int y, int m )
{
	double x1, y1;
	
	x1 = (x-m/2)*(x-m/2);
	y1 = (y-m/2)*(y-m/2);
	
	return( x1 + y1 );
}

/* SSNR process single image */
void my_ssnr_step(const matrix2D< complex<double> > &FTAverageImage,
                  const matrix2D< complex<double> > &FTaverageSubGroup,
		  matrix1D<double> &ssnr, matrix1D<double> &pixel,int z)
{
   //n -> number of images in the subset
   int top, cont;
   double d, w2, w12;
   double w=0.0;
   float preal1, preal2, pimag1, pimag2;
   float resta_real, resta_imag;
   float mod2_diferencia, mod2_media;

//   int dim (int)FTAverageImage.RowNo()/2;
   int n = (int)FTAverageImage.RowNo();
   top = (int) (n/2);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(ssnr){
       mod2_media = 0.0;
       mod2_diferencia = 0.0;
       w=i+0.5;
       w2 = w*w;
       w12 = (w+1)*(w+1);
       cont = 0;
       for( int ii=0; ii<n; ii++ )
	  for( int jj=0; jj<=top; jj++ )
	  {
	     d= distancia2( ii, jj, n);
	     if( d>=w2 && d<w12 ){
                   cont++;
                   mod2_diferencia += norm(dMij(FTAverageImage,ii,jj)-
	                                   dMij(FTaverageSubGroup,ii,jj));
                   mod2_media += norm(dMij(FTAverageImage,ii,jj));						      
	     }
	  }
       mod2_diferencia *= z;
       VEC_ELEM(ssnr,i) += mod2_diferencia/mod2_media;
       VEC_ELEM(pixel,i) = cont;
       }
}		  

/* Convolution of series --------------------------------------------------- */
template <class T>
void series_convolution(matrix1D<T> &series1, matrix1D<T> &series2,
   matrix1D<T> &result, bool FullConvolution) {
   // Store dimension of series
   int dim1=series1.get_dim();
   int dim2=series2.get_dim();

   // Resize series to the size of the resulting series
   // (Zeros are stored in the expanded values)
   series1.resize(dim1+dim2-1);
   series2.resize(series1);
   result .resize(series1);

    // Fourier Transform the two series
    matrix1D<complex <double> > FFT1; FourierTransform(series1, FFT1);
    matrix1D<complex <double> > FFT2; FourierTransform(series2, FFT2);

    // Multiply the vectors element by element to do the 
    // convolution in the Fourier space.
    FFT1*=FFT2;

    // Recover the convolution result by inverse FFT
    InverseFourierTransform(FFT1,result);

    // Restore the dimensions
    series1.resize(dim1);
    series2.resize(dim2);
	
    /* If the full convolution is required, nothing more remains to be done.
       Otherwise, if the valid values are required, return the central ones. */
    if (FullConvolution==FALSE) {
       // First get the maximum dimension of the series,
       // which is the dimension of the result
       int dim=MAX(dim1,dim2);

       // Determine the number of values to discard
       int discard=result.get_dim()-dim;

       // Divide it by two as we have to discard them in both sides of
       // the vector
       discard/=2;  // Integer division is intended

       // Copy required values (simple displacement of values)
       for (int i=STARTINGX(result); i<STARTINGX(result)+dim; i++)
          result(i)=result(i+discard);

       // and finally resize to discard not copied values
       result.resize(dim);
    }
}

/* Numerical derivative of a matrix ----------------------------- */
void numerical_derivative(matrix2D<double> &M, matrix2D<double> &D,
   char direction,int order, int window_size,int polynomial_order) {

   // Set D to be a copy in shape of M
   D.copy_shape(M);

   matrix1D<double> v,rotated; 
   matrix1D<double> ans; // To obtain results 

   // Wrap around version of the Savitzky-Golay coefficients
   int dim=2*window_size+1;
   rotated.resize(dim);

   double *pans=ans.adapt_for_numerical_recipes();
   double *pv=v.adapt_for_numerical_recipes();
   double *protated=rotated.adapt_for_numerical_recipes();

   // Calculate the Savitzky-Golay filter coeficients
   savgol(protated, 2*window_size+1, window_size,
      window_size, order, polynomial_order);

   // Savitzky-Golay filter is returned in wrap-around style, so
   // correct it to use with the convolution routine
   matrix1D<double> coefficients(dim);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(coefficients) {
      int j=i+window_size;
      if(j<dim) coefficients(j)=rotated(i);
      else      coefficients(j-dim)=rotated(i);
   }
   
   // Apply the Savitzky-Golay filter to every row or column
   if (direction=='x') {
      // For every row (values in a row are values of the X direction)
      for (int i=STARTINGY(M);i<=FINISHINGY(M);i++) {
          M.getRow(i,v); 
          series_convolution(v,coefficients,ans,false);
          ans.setRow();
          D.setRow(i,ans);			
      }
   } else if(direction=='y') {
      // For every column (values in a column are values of the Y direction)
      for (int i=STARTINGX(M);i<=FINISHINGX(M);i++) {
          M.getCol(i,v);
          series_convolution(v,coefficients,ans,false);
          ans.setCol();
          D.setCol(i,ans);						
      }
   }	
}

/* Instantiation ---------------------------------------------------------- */
void instantiate_FFT() {
   double          *ptr_ta, *ptr_tb;
   complex<double> *ptr_complex;
   
    RealImag2Complex(ptr_ta, ptr_tb, ptr_complex, 0);
   AmplPhase2Complex(ptr_ta, ptr_tb, ptr_complex, 0);
    Complex2RealImag(ptr_complex, ptr_ta, ptr_tb, 0);
   Complex2AmplPhase(ptr_complex, ptr_ta, ptr_tb, 0);

   matrix2D<double> m1;
   matrix2D<double> R;
   SelFile     SF;
   auto_correlation_matrix(m1,R);
   correlation_matrix(m1,m1,R);
   matrix1D<double> t;
   double s;
   matrix3D<double> m3;
   differential_phase_residual(m1,m1,s,t,t);
   differential_phase_residual(m3,m3,s,t,t);
   fourier_ring_correlation(m1,m1,s,t,t,t);
   fourier_ring_correlation(m3,m3,s,t,t,t);
   //ssnr not defined for volumes
   my_ssnr(m1, SF, s, t, t, t,(bool) 1);

   matrix1D<double> series1;
   series_convolution(series1, series1, series1, false);

   matrix1D<double> v; CenterFFT(v,false);
   matrix2D<double> m; CenterFFT(m,false);
   matrix3D<double> V; CenterFFT(V,false);
   matrix1D< complex<double> > vc; CenterFFT(vc,false);
   matrix2D< complex<double> > mc; CenterFFT(mc,false);
   matrix3D< complex<double> > Vc; CenterFFT(Vc,false);
}
