/***************************************************************************
*
* Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
*  e-mail address 'xmipp@cnb.csic.es'
***************************************************************************/

#ifndef FFT_H
#define FFT_H

#include <complex>

#include "multidim_array.h"
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
    if (idx.size() < 1 || idx.size() > 3)
        REPORT_ERROR(1, "FFT_idx2digfreq: Index is not of the correct size");

    freq.resize(idx.size());

    int size[3];
    v.getSize(size);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
        FFT_IDX2DIGFREQ(idx(i), size[i], freq(i));
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
    if (freq.size() < 1 || freq.size() > 3)
        REPORT_ERROR(1, "digfreq2FFT_idx: freq is not of the correct size");

    idx.resize(freq.size());

    int size[3];
    v.getSize(size);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
        DIGFREQ2FFT_IDX(freq(i), size[i], idx(i));
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
    contfreq.resize(digfreq);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(digfreq)
		contfreq(i) = digfreq(i) / pixel_size;
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
    digfreq.resize(contfreq);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(contfreq)
		digfreq(i) = contfreq(i) * pixel_size;
}

/** @defgroup FourierFormat Format conversions
 * @ingroup Fourier
 */

/** Conversion from whole -> half 
 * @ingroup FourierFormat
 */
void Whole2Half(const MultidimArray< std::complex < double > > & in,
                MultidimArray< std::complex < double > > & out);

/** Conversion from half -> whole 
 * @ingroup FourierFormat
 */
void Half2Whole(const MultidimArray< std::complex < double > > & in,
                MultidimArray< std::complex< double > > & out,
                int oridim);

/** Conversion from complex -> real,imag
 * @ingroup FourierFormat
 */
void Complex2RealImag(const MultidimArray< std::complex < double > > & in,
                      MultidimArray< double > & real,
                      MultidimArray< double > & imag);

/** Conversion from real,imag -> complex
 * @ingroup FourierFormat
 */
void RealImag2Complex(const MultidimArray< double > & real,
                      const MultidimArray< double > & imag,
                      MultidimArray< std::complex < double > > & out);

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
 *     	 MultidimArray<double> x(65);
 * 	 x.setXmippOrigin();
 * 	 double T=0.5;
 * 	 double T1=6;
 * 	 int N1=(int)CEIL(T1/T);
 * 
 * 	 // Fill x with a pulse from -N1 to N1 (-T1 to T1 in continuous)
 * 	 FOR_ALL_ELEMENTS_IN_ARRAY1D(x)
 * 	    if (ABS(i)<=N1) x(i)=1;
 * 
 * 	 // Compute the Fourier transform
 * 	 MultidimArray< std::complex<double> > X;
 * 	 MultidimArray<double> Xmag;
 * 	 FourierTransform(x,X);
 * 	 FFT_magnitude(X,Xmag);
 * 
 * 	 // Compute the frequency axes
 * 	 MultidimArray<double> contfreq(XSIZE(X)), digfreq(XSIZE(X));
 *          FOR_ALL_ELEMENTS_IN_ARRAY1D(X)
 *              FFT_IDX2DIGFREQ(i,XSIZE(X),digfreq(i));
 * 	 digfreq*=2*PI;
 * 	 contfreq=digfreq/T;
 * 
 * 	 // Show all Fourier transforms
 * 	 FOR_ALL_ELEMENTS_IN_ARRAY1D(X) {
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

/** Direct Fourier Transform
 * @ingroup FourierTransforms
 */
void FourierTransform(const MultidimArray< double >& in,
                      MultidimArray< std::complex< double > > & out);

/** Inverse Fourier Transform
 * @ingroup FourierTransforms
 */
void InverseFourierTransform(const MultidimArray< std::complex< double > > & in,
                             MultidimArray< double >& out);

/** Direct Fourier Transform, output half of (centro-symmetric) transform
 * @ingroup FourierTransforms
 */
void FourierTransformHalf(const MultidimArray< double >& in,
                          MultidimArray< std::complex< double > > & out);

/** Inverse Fourier Transform 1D, input half of (centro-symmetric) transform
 * @ingroup FourierTransforms
 */
void InverseFourierTransformHalf(const MultidimArray< std::complex< double > > & in,
                                 MultidimArray< double >& out,
                                 int oridim);

/** @defgroup FourierOperations Operations with the Fourier Transforms
 * @ingroup Fourier
 */

/** CenterFFT
 * @ingroup FourierOperations
 */
template <typename T>
void CenterFFT(MultidimArray< T >& v, bool forward)
{
    if ( v.getDim() == 1 )
    {
        // 1D
        MultidimArray< T > aux;
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

            aux(ip) = DIRECT_A1D_ELEM(v, i);
        }

        // Copy the vector
        for (int i = 0; i < l; i++)
            DIRECT_A1D_ELEM(v, i) = DIRECT_A1D_ELEM(aux, i);
    }
    else if ( v.getDim() == 2 )
    {
        // 2D
        MultidimArray< T > aux;
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
                
                aux(jp) = DIRECT_A2D_ELEM(v, i, j);
            }
            
            // Copy the vector
            for (int j = 0; j < l; j++)
                DIRECT_A2D_ELEM(v, i, j) = DIRECT_A1D_ELEM(aux, j);
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
                
                aux(ip) = DIRECT_A2D_ELEM(v, i, j);
            }
            
            // Copy the vector
            for (int i = 0; i < l; i++)
                DIRECT_A2D_ELEM(v, i, j) = DIRECT_A1D_ELEM(aux, i);
        }
    }
    else if ( v.getDim() == 3 )
    {
        // 3D
        MultidimArray< T > aux;
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

                    aux(jp) = DIRECT_A3D_ELEM(v, k, i, j);
                }
                
                // Copy the vector
                for (int j = 0; j < l; j++)
                    DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, j);
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
                    
                    aux(ip) = DIRECT_A3D_ELEM(v, k, i, j);
                }
                
                // Copy the vector
                for (int i = 0; i < l; i++)
                    DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, i);
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
                    
                    aux(kp) = DIRECT_A3D_ELEM(v, k, i, j);
                }
                
                // Copy the vector
                for (int k = 0; k < l; k++)
                    DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, k);
            }
    }
    else
        REPORT_ERROR(1,"CenterFFT ERROR: Dimension should be 1, 2 or 3");
}

/** FFT shift 1D
 * @ingroup FourierOperations
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(MultidimArray< std::complex< double > > & v, double xshift);

/** FFT shift 2D
 * @ingroup FourierOperations
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(MultidimArray< std::complex< double > > & v, double xshift, double yshift);

/** FFT shift 3D
 * @ingroup FourierOperations
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(MultidimArray< std::complex< double > > & v,
              double xshift,
              double yshift,
              double zshift);

/** Place the origin of the FFT at the center of the vector and back
 * @ingroup FourierOperations
 *
 * Changes the real and the fourier space origin
 */
void CenterOriginFFT(MultidimArray< std::complex< double > > & v, bool forward);


#endif
