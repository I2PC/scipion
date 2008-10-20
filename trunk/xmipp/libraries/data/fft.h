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

#ifndef FFT_H
#define FFT_H

#include <complex>

#include "matrix1d.h"
#include "matrix2d.h"
#include "matrix3d.h"

#include "selfile.h"
#include "funcs.h"

/** @defgroup Fourier Fourier transforms
  * @ingroup DataLibrary
  */

/** @defgroup FourierConverters Index <--> Frequency, Continuous <--> Discrete
 * @ingroup Fourier
 */

/** Index to frequency
 * @ingroup FourierConverters
 *
 * Given an index and a size of the FFT, this function returns the corresponding
 * digital frequency (-1/2 to 1/2)
 */
#define FFT_IDX2DIGFREQ(idx, size, freq) \
    freq = ((double) ((idx) < ((size) >> 1)) ? (idx) : -(size) + (idx)) / \
           (double)(size);

/** Frequency to index (int)
 * @ingroup FourierConverters
 *
 * Given a frequency and a size of the FFT, this macro returns the corresponding
 * integer index
 */
#define DIGFREQ2FFT_IDX(freq, size, idx) { \
        (idx) = (int) (ROUND((size) * (freq))); if ((idx) < 0) (idx) += (int) \
                    (size); }

/** Frequency to index (double)
 * @ingroup FourierConverters
 *
 * Given a frequency and a size of the FFT, this macro returns the corresponding
 * double index
 */
#define DIGFREQ2FFT_IDX_DOUBLE(freq, size, idx) { \
        (idx) = ((size) * (freq)); if ((idx) < 0) (idx) += (size); }

/** Index to frequency
 * @ingroup FourierConverters
 *
 * This function can be used with vectors of any size (1,2,3). The Digital
 * spectrum is limited between -1/2 and 1/2. If the vector has got more than 3
 * coordinates, then an exception is thrown
 */
template <typename T>
void FFT_idx2digfreq(T& v, const Matrix1D< int >& idx, Matrix1D< double >& freq)
{
    if (XSIZE(idx) < 1 || XSIZE(idx) > 3)
        REPORT_ERROR(1, "FFT_idx2digfreq: Index is not of the correct size");

    freq.resize(XSIZE(idx));

    int size[3];
    v.getSize(size);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
        FFT_IDX2DIGFREQ(VEC_ELEM(idx, i), size[i], VEC_ELEM(freq, i));
}

/** Frequency to index
 * @ingroup FourierConverters
 *
 * This function can be used with vectors of any size (1,2,3). The Digital
 * spectrum is lim:confirm bd
 * ited between -1/2 and 1/2. If the vector has got more than 3
 * coordinates, then an exception is thrown
 */
template <typename T>
void digfreq2FFT_idx(T& v, const Matrix1D< double >& freq, Matrix1D< int >& idx)
{
    if (XSIZE(freq) < 1 || XSIZE(freq) > 3)
        REPORT_ERROR(1, "digfreq2FFT_idx: freq is not of the correct size");

    idx.resize(XSIZE(freq));

    int size[3];
    v.getSize(size);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
        DIGFREQ2FFT_IDX(VEC_ELEM(freq, i), size[i], VEC_ELEM(idx, i));
}

/** Digital to Continuous frequency
 * @ingroup FourierConverters
 *
 * The pixel size must be given in Amstrongs. The digital frequency is between
 * [-1/2,1/2]
 */
inline void digfreq2contfreq(const Matrix1D< double >& digfreq,
                             Matrix1D< double >& contfreq,
                             double pixel_size)
{
    contfreq = digfreq / pixel_size;
}

/** Continuous to Digital frequency
 * @ingroup FourierConverters
 *
 * The pixel size must be given in Amstrongs. The digital frequency is between
 * [-1/2,1/2]
 */
inline void contfreq2digfreq(const Matrix1D< double >& contfreq,
                             Matrix1D< double >& digfreq,
                             double pixel_size)
{
    digfreq = contfreq * pixel_size;
}

/** @defgroup FourierFormat Format conversions
 * @ingroup Fourier
 */

/** Conversion from whole -> half 1D
 * @ingroup FourierFormat
 */
void Whole2Half(const Matrix1D< std::complex < double > > & in,
                Matrix1D< std::complex < double > > & out);

/** Conversion from half -> whole 1D
 * @ingroup FourierFormat
 */
void Half2Whole(const Matrix1D< std::complex < double > > & in,
                Matrix1D< std::complex< double > > & out,
                int orixdim);

/** Conversion from whole -> half 2D
 * @ingroup FourierFormat
 */
void Whole2Half(const Matrix2D< std::complex < double > > & in,
                Matrix2D< std::complex < double > > & out);

/** Conversion from half -> whole 2D
 * @ingroup FourierFormat
 */
void Half2Whole(const Matrix2D< std::complex < double > > & in,
                Matrix2D< std::complex< double > > & out,
                int oriydim);

/** Conversion from complex -> real,imag 3D
 * @ingroup FourierFormat
 */
void Complex2RealImag(const Matrix3D< std::complex < double > > & in,
                      Matrix3D< double > & real,
                      Matrix3D< double > & imag);

/** Conversion from real,imag -> complex 3D
 * @ingroup FourierFormat
 */
void RealImag2Complex(const Matrix3D< double > & real,
                      const Matrix3D< double > & imag,
                      Matrix3D< std::complex < double > > & out);

/** @defgroup FourierTransforms Fourier Transforms
 * @ingroup Fourier
 *
 *  The theoretical relationship between the Fourier transform of a discrete
 *  signal and the Fourier transform of the continuous signal is
 *  
 *  X(e^jw)=1/T*X_c(jw/T)
 *  
 *  Xmipp is not computing X(e^jw) but samples from it so that
 *  
 *  X(e^jw)=N*X_XMIPP[k]
 *  
 *  where N is the length of the signal being transformed and X_XMIPP[k]
 *  is the k-th sample.
 *
 *  The following program illustrates how the continuous, discrete and
 *  Xmipp Fourier transform relate
 *
 * @code
 * #include <data/matrix1d.h>
 * #include <data/fft.h>
 * 
 * double discreteTransform(double w, int N1) {
 *    if (w==0) return 2*N1+1;
 *    else return sin(w*(N1+0.5))/sin(0.5*w);
 * }
 * 
 * double continuousTransform(double W, double T1) {
 *    if (W==0) return 2*T1;
 *    else return 2*sin(W*T1)/W;
 * }
 * 
 * int main() {
 *     try {
 *     	 Matrix1D<double> x(65);
 * 	 x.setXmippOrigin();
 * 	 double T=0.5;
 * 	 double T1=6;
 * 	 int N1=(int)CEIL(T1/T);
 * 
 * 	 // Fill x with a pulse from -N1 to N1 (-T1 to T1 in continuous)
 * 	 FOR_ALL_ELEMENTS_IN_MATRIX1D(x)
 * 	    if (ABS(i)<=N1) x(i)=1;
 * 
 * 	 // Compute the Fourier transform
 * 	 Matrix1D< std::complex<double> > X;
 * 	 Matrix1D<double> Xmag;
 * 	 FourierTransform(x,X);
 * 	 FFT_magnitude(X,Xmag);
 * 
 * 	 // Compute the frequency axes
 * 	 Matrix1D<double> contfreq(XSIZE(X)), digfreq(XSIZE(X));
 *          FOR_ALL_ELEMENTS_IN_MATRIX1D(X)
 *              FFT_IDX2DIGFREQ(i,XSIZE(X),digfreq(i));
 * 	 digfreq*=2*PI;
 * 	 contfreq=digfreq/T;
 * 
 * 	 // Show all Fourier transforms
 * 	 FOR_ALL_ELEMENTS_IN_MATRIX1D(X) {
 * 	     if (digfreq(i)>=0)
 *                 std::cout << digfreq(i) << " " << contfreq(i) << " "
 * 		          << XSIZE(X)*Xmag(i) << " "
 * 			  << ABS(discreteTransform(digfreq(i),N1)) << " "
 * 			  << ABS(continuousTransform(contfreq(i),T1)/T)
 * 			  << std::endl;
 * 	 }
 *     } catch (Xmipp_error XE) {
 *     	   std::cout << XE << std::endl;
 *     }
 *     return 0;
 * }
 * @endcode
 */

/** Direct Fourier Transform 1D
 * @ingroup FourierTransforms
 */
void FourierTransform(const Matrix1D< double >& in,
                      Matrix1D< std::complex< double > > & out);

/** Direct Fourier Transform 2D
 * @ingroup FourierTransforms
 */
void FourierTransform(const Matrix2D< double >& in,
                      Matrix2D< std::complex< double > > & out);

/** Direct Fourier Transform 3D
 * @ingroup FourierTransforms
 */
void FourierTransform(const Matrix3D< double >& in,
                      Matrix3D< std::complex< double > > & out);

/** Inverse Fourier Transform 1D
 * @ingroup FourierTransforms
 */
void InverseFourierTransform(const Matrix1D< std::complex< double > > & in,
                             Matrix1D< double >& out);

/** Inverse Fourier Transform 2D
 * @ingroup FourierTransforms
 */
void InverseFourierTransform(const Matrix2D< std::complex< double > > & in,
                             Matrix2D< double >& out);

/** Inverse Fourier Transform 3D
 * @ingroup FourierTransforms
 */
void InverseFourierTransform(const Matrix3D< std::complex< double > > & in,
                             Matrix3D< double >& out);

/** Direct Fourier Transform 1D, output half of (centro-symmetric) transform
 * @ingroup FourierTransforms
 */
void FourierTransformHalf(const Matrix1D< double >& in,
                          Matrix1D< std::complex< double > > & out);

/** Direct Fourier Transform 2D, output half of (centro-symmetric) transform
 * @ingroup FourierTransforms
 */
void FourierTransformHalf(const Matrix2D< double >& in,
                          Matrix2D< std::complex< double > > & out);

/** Inverse Fourier Transform 1D, input half of (centro-symmetric) transform
 * @ingroup FourierTransforms
 */
void InverseFourierTransformHalf(const Matrix1D< std::complex< double > > & in,
                                 Matrix1D< double >& out,
                                 int orixdim);


/** Inverse Fourier Transform 2D, input half of (centro-symmetric) transform
 * @ingroup FourierTransforms
 */
void InverseFourierTransformHalf(const Matrix2D< std::complex< double > > & in,
                                 Matrix2D< double >& out,
                                 int oriydim);

/** Complex Direct Fourier Transform 1D
 * @ingroup FourierTransforms
 */
void FourierTransform(const Matrix1D< std::complex< double > >& in,
                      Matrix1D< std::complex< double > > & out);

/** Complex Inverse Fourier Transform 1D
 * @ingroup FourierTransforms
 */
void InverseFourierTransform(const Matrix1D< std::complex< double > > & in,
                             Matrix1D< std::complex< double > >& out);

/** Complex Direct Fourier Transform 1D, output half of (centro-symmetric) transform
 * @ingroup FourierTransforms
 */
void FourierTransformHalf(const Matrix1D<std::complex<double> > &in,
                          Matrix1D< std::complex<double> > &out);

/** Complex Inverse Fourier Transform 1D, input half of (centro-symmetric) transform
 * @ingroup FourierTransforms
 */
void InverseFourierTransformHalf(const Matrix1D< std::complex<double> > &in,
                                 Matrix1D<std::complex<double> > &out, 
				 int orixdim);


/** @defgroup FourierOperations Operations with the Fourier Transforms
 * @ingroup Fourier
 */

/** CenterFFT 1D
 * @ingroup FourierOperations
 */
template <typename T>
void CenterFFT(Matrix1D< T >& v, bool forward)
{
    Matrix1D< T > aux;
    int l, shift;

    l = XSIZE(v);
    aux.resize(l);
    shift = (int)(l / 2);

    if (!forward)
        shift = -shift;

    // Shift the input in an auxiliar vector
    for (int i = 0; i < l; i++)
    {
        int ip = i + shift;

        if (ip < 0)
            ip += l;
        else if (ip >= l)
            ip -= l;

        aux(ip) = DIRECT_VEC_ELEM(v, i);
    }

    // Copy the vector
    for (int i = 0; i < l; i++)
        DIRECT_VEC_ELEM(v, i) = DIRECT_VEC_ELEM(aux, i);
}

/** CenterFFT 2D
 * @ingroup FourierOperations
 */
template <typename T>
void CenterFFT(Matrix2D< T >& v, bool forward)
{
    Matrix1D< T > aux;
    int l, shift;

    // Shift in the X direction
    l = XSIZE(v);
    aux.resize(l);
    shift = (int)(l / 2);

    if (!forward)
        shift = -shift;

    for (int i = 0; i < YSIZE(v); i++)
    {
        // Shift the input in an auxiliar vector
        for (int j = 0; j < l; j++)
        {
            int jp = j + shift;

            if (jp < 0)
                jp += l;
            else if (jp >= l)
                jp -= l;

            aux(jp) = DIRECT_MAT_ELEM(v, i, j);
        }

        // Copy the vector
        for (int j = 0; j < l; j++)
            DIRECT_MAT_ELEM(v, i, j) = DIRECT_VEC_ELEM(aux, j);
    }

    // Shift in the Y direction
    l = YSIZE(v);
    aux.resize(l);
    shift = (int)(l / 2);

    if (!forward)
        shift = -shift;

    for (int j = 0; j < XSIZE(v); j++)
    {
        // Shift the input in an auxiliar vector
        for (int i = 0; i < l; i++)
        {
            int ip = i + shift;

            if (ip < 0)
                ip += l;
            else if (ip >= l)
                ip -= l;

            aux(ip) = DIRECT_MAT_ELEM(v, i, j);
        }

        // Copy the vector
        for (int i = 0; i < l; i++)
            DIRECT_MAT_ELEM(v, i, j) = DIRECT_VEC_ELEM(aux, i);
    }
}

/** CenterFFT 3D
 * @ingroup FourierOperations
 */
template <class T>
void CenterFFT(Matrix3D< T >& v, bool forward)
{
    Matrix1D< T > aux;
    int l, shift;

    // Shift in the X direction
    l = XSIZE(v);
    aux.resize(l);
    shift = (int)(l / 2);

    if (!forward)
        shift = -shift;

    for (int k = 0; k < ZSIZE(v); k++)
        for (int i = 0; i < YSIZE(v); i++)
        {
            // Shift the input in an auxiliar vector
            for (int j = 0; j < l; j++)
            {
                int jp = j + shift;

                if (jp < 0)
                    jp += l;
                else if (jp >= l)
                    jp -= l;

                aux(jp) = DIRECT_VOL_ELEM(v, k, i, j);
            }

            // Copy the vector
            for (int j = 0; j < l; j++)
                DIRECT_VOL_ELEM(v, k, i, j) = DIRECT_VEC_ELEM(aux, j);
        }

    // Shift in the Y direction
    l = YSIZE(v);
    aux.resize(l);
    shift = (int)(l / 2);

    if (!forward)
        shift = -shift;

    for (int k = 0; k < ZSIZE(v); k++)
        for (int j = 0; j < XSIZE(v); j++)
        {
            // Shift the input in an auxiliar vector
            for (int i = 0; i < l; i++)
            {
                int ip = i + shift;

                if (ip < 0)
                    ip += l;
                else if (ip >= l)
                    ip -= l;

                aux(ip) = DIRECT_VOL_ELEM(v, k, i, j);
            }

            // Copy the vector
            for (int i = 0; i < l; i++)
                DIRECT_VOL_ELEM(v, k, i, j) = DIRECT_VEC_ELEM(aux, i);
        }

    // Shift in the Z direction
    l = ZSIZE(v);
    aux.resize(l);
    shift = (int)(l / 2);

    if (!forward)
        shift = -shift;

    for (int i = 0; i < YSIZE(v); i++)
        for (int j = 0; j < XSIZE(v); j++)
        {
            // Shift the input in an auxiliar vector
            for (int k = 0; k < l; k++)
            {
                int kp = k + shift;
                if (kp < 0)
                    kp += l;
                else if (kp >= l)
                    kp -= l;

                aux(kp) = DIRECT_VOL_ELEM(v, k, i, j);
            }

            // Copy the vector
            for (int k = 0; k < l; k++)
                DIRECT_VOL_ELEM(v, k, i, j) = DIRECT_VEC_ELEM(aux, k);
        }
}

/** FFT shift 1D
 * @ingroup FourierOperations
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(Matrix1D< std::complex< double > > & v, double xshift);

/** FFT shift 2D
 * @ingroup FourierOperations
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(Matrix2D< std::complex< double > > & v, double xshift, double yshift);

/** FFT shift 3D
 * @ingroup FourierOperations
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(Matrix3D< std::complex< double > > & v,
              double xshift,
              double yshift,
              double zshift);

/** Place the origin of the 1D FFT at the center of the vector and back
 * @ingroup FourierOperations
 *
 * Changes the real and the fourier space origin
 */
void CenterOriginFFT(Matrix1D< std::complex< double > > & v, bool forward);

/** Place the origin of the 2D FFT at the center of the image and back
 * @ingroup FourierOperations
 *
 * Changes the real and the fourier space origin
 */
void CenterOriginFFT(Matrix2D< std::complex< double > > & v, bool forward);

/** Place the origin of the 3D FFT at the center of the volume and back
 * @ingroup FourierOperations
 *
 * Changes the real and the fourier space origin
 */
void CenterOriginFFT(Matrix3D< std::complex< double > > & v, bool forward);

/** FFT Magnitude 1D
 * @ingroup FourierOperations
 */
void FFT_magnitude(const Matrix1D< std::complex< double > > & v,
                   Matrix1D< double >& mag);

/** FFT Magnitude 2D
 * @ingroup FourierOperations
 */
void FFT_magnitude(const Matrix2D< std::complex< double > > & v,
                   Matrix2D< double >& mag);

/** FFT Magnitude 3D
 * @ingroup FourierOperations
 */
void FFT_magnitude(const Matrix3D< std::complex< double > > & v,
                   Matrix3D< double >& mag);

/** FFT Phase 1D
 * @ingroup FourierOperations
 */
void FFT_phase(const Matrix1D< std::complex< double > > & v,
               Matrix1D< double >& phase);

/** FFT Phase 2D
 * @ingroup FourierOperations
 */
void FFT_phase(const Matrix2D< std::complex< double > > & v,
               Matrix2D< double >& phase);

/** FFT Phase 3D
 * @ingroup FourierOperations
 */
void FFT_phase(const Matrix3D< std::complex< double > > & v,
               Matrix3D< double >& phase);

/** Autocorrelation function of a Xmipp vector
 * @ingroup FourierOperations
 *
 * Fast calcuation of the autocorrelation vector of a given one using Fast
 * Fourier Transform. (Using the correlation theorem)
 */
template <typename T>
void auto_correlation_vector(Matrix1D< T > const & Img, Matrix1D< double >& R)
{
    // Compute the Fourier Transform
    Matrix1D< std::complex< double > > FFT1;
    FourierTransform(Img, FFT1);

    // Multiply FFT1 * FFT1'
    FOR_ALL_ELEMENTS_IN_MATRIX1D(FFT1)
        FFT1(i) *= conj(FFT1(i));

    // Invert the product, in order to obtain the correlation image
    InverseFourierTransform(FFT1, R);

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Autocorrelation function of a Xmipp vector
 * @ingroup FourierOperations
 *
 * Fast calcuation of the correlation matrix on two matrices using Fast Fourier
 * Transform. (Using the correlation theorem). The output matrix must be already
 * resized
 */
template <typename T>
void correlation_vector(Matrix1D< T > const & m1,
                        Matrix1D< T > const & m2,
                        Matrix1D< double >& R)
{
    // Compute the Fourier Transforms
    Matrix1D< std::complex< double > > FFT1;
    FourierTransform(m1, FFT1);
    Matrix1D< std::complex< double > > FFT2;
    FourierTransform(m2, FFT2);

    // Multiply FFT1 * FFT2'
    FOR_ALL_ELEMENTS_IN_MATRIX1D(FFT1)
        FFT1(i) *= conj(FFT2(i));

    // Invert the product, in order to obtain the correlation image
    InverseFourierTransform(FFT1, R);

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Autocorrelation function of a Xmipp matrix
 * @ingroup FourierOperations
 *
 * Fast calcuation of the autocorrelation matrix of a given one using Fast
 * Fourier Transform. (Using the correlation theorem)
 */
template <typename T>
void auto_correlation_matrix(Matrix2D< T > const & Img, Matrix2D< double >& R)
{
    // Compute the Fourier Transform
    Matrix2D< std::complex< double > > FFT1;
    FourierTransform(Img, FFT1);

    // Multiply FFT1 * FFT1'
    FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT1)
        FFT1(i, j) *= conj(FFT1(i, j));

    // Invert the product, in order to obtain the correlation image
    InverseFourierTransform(FFT1, R);

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Autocorrelation function of a Xmipp matrix
 * @ingroup FourierOperations
 *
 * Fast calcuation of the correlation matrix on two matrices using Fast Fourier
 * Transform. (Using the correlation theorem). The output matrix must be already
 * resized
 */
template <typename T>
void correlation_matrix(Matrix2D< T > const & m1,
                        Matrix2D< T > const & m2,
                        Matrix2D< double >& R)
{
    // Compute the Fourier Transforms
    Matrix2D< std::complex< double > > FFT1;
    FourierTransform(m1, FFT1);
    Matrix2D< std::complex< double > > FFT2;
    FourierTransform(m2, FFT2);

    // Multiply FFT1 * FFT2'
    FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT1)
        FFT1(i, j) *= conj(FFT2(i, j));

    // Invert the product, in order to obtain the correlation image
    InverseFourierTransform(FFT1, R);

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Series convolution function.
 * @ingroup FourierOperations
 *
 * Gives the convolution of two series given as Xmipp Vectors. Result is stored
 * in result vector. Fast calcuation of the convolution result using Fast
 * Fourier Transform. If FullConvolution -set by default to FALSE- is TRUE the
 * full convolution series is returned. Otherwise the convolution vector refers
 * only to the valid values, whose number is the greater dimension of the two
 * series.
 *
 * Note: Complex numbers are allowed
 */
template <typename T>
void series_convolution(Matrix1D< T >& series1,
                        Matrix1D< T >& series2,
                        Matrix1D< T >& result,
                        bool FullConvolution = false)
{
    // Store dimension of series
    int dim1 = XSIZE(series1);
    int dim2 = XSIZE(series2);

    // Resize series to the size of the resulting series
    // (Zeros are stored in the expanded values)
    series1.resize(dim1 + dim2 - 1);
    series2.resize(series1);
    result.resize(series1);

    // Fourier Transform the two series
    Matrix1D< std::complex< double> > FFT1;
    FourierTransform(series1, FFT1);

    Matrix1D< std::complex< double > > FFT2;
    FourierTransform(series2, FFT2);

    // Multiply the vectors element by element to do the convolution in the
    // Fourier space
    FFT1 *= FFT2;

    // Recover the convolution result by inverse FFT
    InverseFourierTransform(FFT1, result);

    // Restore the dimensions
    series1.resize(dim1);
    series2.resize(dim2);

    // If the full convolution is required, nothing more remains to be done.
    // Otherwise, if the valid values are required, return the central ones
    if (!FullConvolution)
    {
        // First get the maximum dimension of the series, which is the dimension
        // of the result
        int dim = XMIPP_MAX(dim1, dim2);

        // Determine the number of values to discard
        int discard = XSIZE(result) - dim;

        // Divide it by two as we have to discard them in both sides of the
        // vector
        discard /= 2;  // Integer division is intended

        // Copy required values (simple displacement of values)
        for (int i = STARTINGX(result); i < STARTINGX(result) + dim; i++)
            result(i) = result(i + discard);

        // And finally resize to discard not copied values
        result.resize(dim);
    }
}

/** Numerical_derivative
 * @ingroup FourierOperations
 *
 * This function computes the numerical derivative of a matrix in Y direction
 * (rows) or X direction (columns) of a given matrix, using a Savitzky-Golay
 * filter on every row or column, and then convolving.
 *
 * Input matrix is M, result is stored in D, direction can have values of 'x' or
 *  'y'. Order is the derivative order. Window size and polynomial order are
 * parameters for the Savitzky-Golay filter that define the number of points
 * forward (+window_size) and backward (-window_size) considered to calculate
 * the filter, and the degree of the polynomial to interpolate these values,
 * respectively.
 *
 * Default values are window_size=2 and polynomial_order=4, which are equivalent
 * to a 5-point algorithm to calculate the derivate and give good results. But
 * they can be changed provided that polynomial_order <= 2*window_size.
 *
 * As one can expect, the values of the matrix in a border of size window_size
 * are not accurate ones, as there aren't enough points to perform the
 * estimation. As a rule of thumb, the greater the window, the more the
 * filtering and the less the precission of the derivatives, and the greater the
 * order of the polynomial, the greater the precision.
 */
void numerical_derivative(Matrix2D< double >& M,
                          Matrix2D< double >& D,
                          char direction,
                          int order,
                          int window_size = 2,
                          int polynomial_order = 4);

#endif
