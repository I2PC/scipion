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
/*****************************************************************************/
/* Fourier Transforms                                                        */
/*****************************************************************************/
#ifndef _XMIPP_FFT_HH
   #define _XMIPP_FFT_HH

#include <complex>
#include "xmippMatrices1D.hh"
#include "xmippMatrices2D.hh"
#include "xmippMatrices3D.hh"
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippFuncs.hh>

/**@name Fourier Transforms */
//@{
/**@name Index <--> Frequency, Continuous <--> Discrete*/
//@{
   /** Index to frequency.
       Given an index and a size of the FFT, this function returns the
       corresponding digital frequency (-1/2 to 1/2). */
   #define FFT_IDX2DIGFREQ(idx,size,freq) \
       freq=((double)((idx)<((size)>>1))?(idx):-(size)+(idx))/(double)(size);

   /** Frequency to index (int).
       Given a frequency and a size of the FFT, this macro returns the
       corresponding integer index. */
   #define DIGFREQ2FFT_IDX(freq,size,idx) {\
       (idx)=(int)(ROUND((size)*(freq))); if ((idx)<0) (idx)+=(int)(size);}

   /** Frequency to index (double).
       Given a frequency and a size of the FFT, this macro returns the
       corresponding double index. */
   #define DIGFREQ2FFT_IDX_DOUBLE(freq,size,idx) {\
       (idx)=((size)*(freq)); if ((idx)<0) (idx)+=(size);}

   /** Index to frequency.
       This function can be used with vectors of any size (1,2,3).
       The Digital spectrum is limited between -1/2 and 1/2.
       If the vector has got more than 3 coordinates, then an exception
       is thrown*/
   template <class T>
   void FFT_idx2digfreq(T &v, const matrix1D<int> &idx,
      matrix1D<double> &freq) {
         if (XSIZE(idx)<1 || XSIZE(idx)>3)
            REPORT_ERROR(1,"FFT_idx2digfreq: Index is not of the correct size");
         freq.resize(XSIZE(idx));

         int size[3]; v.get_size(size);
         FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
            FFT_IDX2DIGFREQ(VEC_ELEM(idx,i),size[i],VEC_ELEM(freq,i));
   }

   /** Frequency to index.
       This function can be used with vectors of any size (1,2,3).
       The Digital spectrum is limited between -1/2 and 1/2.
       If the vector has got more than 3 coordinates, then an exception
       is thrown*/
   template <class T>
   void digfreq2FFT_idx(T &v, const matrix1D<double> &freq,
      matrix1D<int> &idx) {
         if (XSIZE(freq)<1 || XSIZE(freq)>3)
            REPORT_ERROR(1,"digfreq2FFT_idx: freq is not of the correct size");
         idx.resize(XSIZE(freq));

         int size[3]; v.get_size(size);
         FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
            DIGFREQ2FFT_IDX(VEC_ELEM(freq,i),size[i],VEC_ELEM(idx,i));
   }


   /** Digital to Continuous frequency.
       The pixel size must be given in Amstrongs. The digital frequency is
       between [-1/2,1/2].*/
   inline void digfreq2contfreq(const matrix1D<double> &digfreq,
      matrix1D<double> &contfreq, double pixel_size)
      {contfreq=digfreq/pixel_size;}
   
   /** Continuous to Digital frequency.
       The pixel size must be given in Amstrongs. The digital frequency is
       between [-1/2,1/2].*/
   inline void contfreq2digfreq(const matrix1D<double> &contfreq,
      matrix1D<double> &digfreq, double pixel_size)
      {digfreq=contfreq*pixel_size;}
//@}

/**@name Format conversions.*/
//@{
    /** Real/Imaginary --> Complex.
        The output array(s) must be already resized*/
    template <class T>
       void RealImag2Complex(const T *_real, const T *_imag,
          complex<double> *_complex, int length);

    /** Amplitude/Phase --> Complex.
        The output array(s) must be already resized*/
    template <class T>
       void AmplPhase2Complex(const T *_ampl, const T *_phase,
          complex<double> *_complex, int length);

    /** Complex --> Real/Imag.
        The output array(s) must be already resized*/
    template <class T>
       void Complex2RealImag(const complex<double> *_complex,
          T *_real, T *_imag, int length);

    /** Complex --> Amplitude/Phase.
        The output array(s) must be already resized*/
    template <class T>
       void Complex2AmplPhase(const complex<double> *_complex,
          T *_ampl, T *_phase, int length);

    /** Conversion from whole -> half 2D.*/
    void Whole2Half(const matrix2D<complex<double> > &in, 
		    matrix2D<complex<double> > &out );

    /** Conversion from half -> whole 2D.*/
    void Half2Whole(const matrix2D<complex<double> > &in, 
		    matrix2D<complex<double> > &out, int oriydim);

//@}

/**@name Fourier Transforms */
//@{
    /** Direct Fourier Transform 1D.*/
    void FourierTransform(const matrix1D<double> &in,
       matrix1D< complex<double> > &out);

    /** Direct Fourier Transform 2D.*/
    void FourierTransform(const matrix2D<double> &in,
       matrix2D< complex<double> > &out);

    /** Direct Fourier Transform 3D.*/
    void FourierTransform(const matrix3D<double> &in,
       matrix3D< complex<double> > &out);

    /** Inverse Fourier Transform 1D.*/
    void InverseFourierTransform(const matrix1D< complex<double> > &in,
       matrix1D<double> &out);

    /** Inverse Fourier Transform 2D.*/
    void InverseFourierTransform(const matrix2D< complex<double> > &in,
       matrix2D<double> &out);

    /** Inverse Fourier Transform 3D.*/
    void InverseFourierTransform(const matrix3D< complex<double> > &in,
       matrix3D<double> &out);

    /** Direct Fourier Transform 2D, output half of (centro-symmetric) transform */
    void FourierTransformHalf(const matrix2D<double> &in,
       matrix2D< complex<double> > &out);

    /** Inverse Fourier Transform 2D, input half of (centro-symmetric) transform */
    void InverseFourierTransformHalf(const matrix2D< complex<double> > &in,
       matrix2D<double> &out, int oriydim);

//@}

/**@name Operations with the Fourier Transforms */
//@{
    /** CenterFFT 1D. */
    template <class T>
    void CenterFFT(matrix1D<T> &v, bool forward);

    /** CenterFFT 2D. */
    template <class T>
    void CenterFFT(matrix2D<T> &v, bool forward);

    /** CenterFFT 3D. */
    template <class T>
    void CenterFFT(matrix3D<T> &v, bool forward);

    /** FFT Magnitude 1D. */
    void FFT_magnitude(const matrix1D< complex<double> > &v,
       matrix1D<double> &mag);

    /** FFT Magnitude 2D. */
    void FFT_magnitude(const matrix2D< complex<double> > &v,
       matrix2D<double> &mag);

    /** FFT Magnitude 3D. */
    void FFT_magnitude(const matrix3D< complex<double> > &v,
       matrix3D<double> &mag);

    /** FFT Phase 1D. */
    void FFT_phase(const matrix1D< complex<double> > &v,
       matrix1D<double> &phase);

    /** FFT Phase 2D. */
    void FFT_phase(const matrix2D< complex<double> > &v,
       matrix2D<double> &phase);

    /** FFT Phase 3D. */
    void FFT_phase(const matrix3D< complex<double> > &v,
       matrix3D<double> &phase);

    /** Autocorrelation function of an Xmipp matrix.
        Fast calcuation of the autocorrelation matrix of a given one using
        Fast Fourier Transform. (Using the correlation theorem) */
  template <class T>
     void auto_correlation_matrix(matrix2D<T> const &Img,
        matrix2D<double> &R);

   /** Autocorrelation function of an Xmipp matrix.
       Fast calcuation of the correlation matrix on two matrices using
       Fast Fourier Transform. (Using the correlation theorem).
       The output matrix must be already resized*/
  template <class T>
     void correlation_matrix(matrix2D<T> const &m1,matrix2D<T> const &m2,
        matrix2D<double> &R);

   /** Fourier-Ring-Correlation between two 2D-matrices using Fast Fourier Transform.*/
template <class T>
     void fourier_ring_correlation(matrix2D<T> const &m1,matrix2D<T> const &m2, double sampling_rate,
        matrix1D<double> &frequency, matrix1D<double> &frc, matrix1D<double> &frc_noise);

   /** Fourier-Ring-Correlation between two 3D-matrices using Fast Fourier Transform.*/
template <class T>
     void fourier_ring_correlation(matrix3D<T> const &m1,matrix3D<T> const &m2, double sampling_rate,
        matrix1D<double> &frequency, matrix1D<double> &frc, matrix1D<double> &frc_noise);

   /** Differential Phase Residualbetween two 2D-matrices using Fast Fourier Transform.*/
template <class T>
     void differential_phase_residual(matrix2D<T> const &m1,matrix2D<T> const &m2, double sampling_rate,
        matrix1D<double> &frequency, matrix1D<double> &dpr);

   /** Differential Phase Residualbetween two 3D-matrices using Fast Fourier Transform.*/
template <class T>
     void differential_phase_residual(matrix3D<T> const &m1,matrix3D<T> const &m2, double sampling_rate,
        matrix1D<double> &frequency, matrix1D<double> &dpr);
   /** Signal to noise ratio for 2D */
template <class T>
     void my_ssnr(matrix2D<T> const &m1, SelFile  &SF, double sampling_rate,
		 matrix1D<double> &freq, matrix1D<double> &ssnr,
		 matrix1D<double> &pixel,bool apply_geo);

   /** Signal to noise ratio for 2D, process one image. freq and ssnr
       must have correct size */
     void my_ssnr_step(matrix2D< complex<double> > const &AverageImage,
                       matrix2D< complex<double> > const &FTaverageSubGroup,
		       matrix1D<double> &ssnr,matrix1D<double> &pixel, int n);
/** Square distance from the point (x,y) to (m/2,m/2)  (used by my_ssnr_step)
*/
 double distancia2( int x, int y, int m );
		 
   /** Series convolution function. Gives the convolution of two series
       given as Xmipp Vectors. Result is stored in result vector.
       Fast calcuation of the convolution result using
       Fast Fourier Transform. If FullConvolution, set by default to
       FALSE, is TRUE the full convolution series is returned. Otherwise
       the convolution vector refers only to the valid values, whose
       number is the greater dimension of the two series.
       Note: Complex numbers are allowed */
   template <class T>
      void series_convolution(matrix1D<T> &series1,matrix1D<T> &series2,
         matrix1D<T> &result, bool FullConvolution=FALSE);

   /** Numerical_derivative: This function computes the numerical derivative
       of a matrix in Y direction (rows) or X direction (columns) of a given
       matrix, using a Savitzky-Golay filter on every row or column, and then
       convolving.Input matrix is M, result is stored in D, direction can have
       values of 'x' or 'y'. Order is the derivative order.
       Window size and polynomial order are parameters for the
       Savitzky-Golay filter that define the number of points forward
       (+window_size)
       and backward (-window_size) considered to calculate the filter, and the degree of
       the polynomial to interpolate these values, respectively. Default values are
       window_size=2 and polynomial_order=4, which are equivalent to a
       5-point algorithm to calculate the derivate and give good results.
       But they can be changed provided that polynomial_order <= 2*window_size.
       As one can expect, the values of the matrix in a border of size
       window_size are not accurate ones, as there aren't enough points to perform
       the estimation. As a rule of thumb, the greater the window, the more the
       filtering and the less the precission of the derivatives, and the
       greater the order of the polynomial, the greater the precission.  */
    void numerical_derivative(matrix2D<double> &M,matrix2D<double> &D,
       char direction, int order, int window_size=2, int polynomial_order=4);
//@}

//@}
#endif
