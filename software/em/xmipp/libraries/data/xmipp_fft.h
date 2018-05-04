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
#include "xmipp_funcs.h"

/** @defgroup Fourier Fourier transforms
  * @ingroup DataLibrary
  */
//@{
/** @name Index <--> Frequency, Continuous <--> Discrete
 */
//@{
/** Index to frequency
 *
 * Given an index and a size of the FFT, this function returns the corresponding
 * digital frequency (-1/2 to 1/2)
 */
#define FFT_IDX2DIGFREQ(idx, size, freq) \
    freq = (size<=1)? 0:(( (((int)idx) <= (((int)(size)) >> 1)) ? ((int)(idx)) : -((int)(size)) + ((int)(idx))) / \
           (double)(size));

#define FFT_IDX2DIGFREQ_DOUBLE(idx, size, freq) \
    freq = (size<=1)? 0:(( (((double)idx) <= (((double)(size)) / 2.0)) ? ((double)(idx)) : -((double)(size)) + ((double)(idx))) / \
           (double)(size));

#define FFT_IDX2DIGFREQ_FAST(idx, size, size_2, isize, freq) \
    freq = ( ((idx) <= (size_2)) ? (idx) : -(size) + (idx) ) * (isize);
/** Frequency to index (int)
 *
 * Given a frequency and a size of the FFT, this macro returns the corresponding
 * integer index
 */
#define DIGFREQ2FFT_IDX(freq, size, idx) { \
        (idx) = (int) (round((size) * (freq))); if ((idx) < 0) (idx) += (int) \
                    (size); }

/** Frequency to index (double)
 *
 * Given a frequency and a size of the FFT, this macro returns the corresponding
 * double index
 */
#define DIGFREQ2FFT_IDX_DOUBLE(freq, size, idx) { \
        (idx) = ((size) * (freq)); if ((idx) < 0) (idx) += (size); }

/** Index to frequency
 *
 * This function can be used with vectors of any size (1,2,3). The Digital
 * spectrum is limited between -1/2 and 1/2. If the vector has got more than 3
 * coordinates, then an exception is thrown
 */
template <typename T>
void FFT_idx2digfreq(T& v, const Matrix1D< int >& idx, Matrix1D< double >& freq)
{
    if (VEC_XSIZE(idx) < 1 || VEC_XSIZE(idx) > 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "FFT_idx2digfreq: Index is not of the correct size");

    freq.resizeNoCopy(VEC_XSIZE(idx));

    switch (VEC_XSIZE(idx))
    {
    case 3:
        FFT_IDX2DIGFREQ(VEC_ELEM(idx,2), ZSIZE(v), VEC_ELEM(freq,2));
    case 2:
        FFT_IDX2DIGFREQ(VEC_ELEM(idx,1), YSIZE(v), VEC_ELEM(freq,1));
    case 1:
        FFT_IDX2DIGFREQ(VEC_ELEM(idx,0), XSIZE(v), VEC_ELEM(freq,0));
    }
}

/** Frequency to index
 *
 * This function can be used with vectors of any size (1,2,3). The Digital
 * spectrum is lim:confirm bd
 * ited between -1/2 and 1/2. If the vector has got more than 3
 * coordinates, then an exception is thrown
 */
template <typename T>
void digfreq2FFT_idx(T& v, const Matrix1D< double >& freq, Matrix1D< int >& idx)
{
    if (VEC_XSIZE(freq) < 1 || VEC_XSIZE(freq) > 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "digfreq2FFT_idx: freq is not of the correct size");

    idx.resizeNoCopy(VEC_XSIZE(freq));

    int size[3];
    v.getSize(size);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
    DIGFREQ2FFT_IDX(VEC_ELEM(freq,i), size[i], VEC_ELEM(idx,i));
}

/** Digital to Continuous frequency
 *
 * The pixel size must be given in Amstrongs. The digital frequency is between
 * [-1/2,1/2]
 */
inline void digfreq2contfreq(const Matrix1D< double >& digfreq,
                             Matrix1D< double >& contfreq,
                             double pixel_size)
{
    contfreq.resizeNoCopy(digfreq);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(digfreq)
    VEC_ELEM(contfreq,i) = VEC_ELEM(digfreq,i) / pixel_size;
}

/** Continuous to Digital frequency
 *
 * The pixel size must be given in Amstrongs. The digital frequency is between
 * [-1/2,1/2]
 */
inline void contfreq2digfreq(const Matrix1D< double >& contfreq,
                             Matrix1D< double >& digfreq,
                             double pixel_size)
{
    digfreq.resizeNoCopy(contfreq);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(contfreq)
    VEC_ELEM(digfreq,i) = VEC_ELEM(contfreq,i) * pixel_size;
}
//@}

/** @name Format conversions
 */
//@{
/** Conversion from whole -> half
 */
void Whole2Half(const MultidimArray< std::complex < double > > & in,
                MultidimArray< std::complex < double > > & out);

/** Conversion from half -> whole
 */
void Half2Whole(const MultidimArray< std::complex < double > > & in,
                MultidimArray< std::complex< double > > & out,
                size_t oridim);

/** Conversion from complex -> real,imag
 */
void Complex2RealImag(const MultidimArray< std::complex < double > > & in,
                      MultidimArray< double > & real,
                      MultidimArray< double > & imag);

/** Conversion from real,imag -> complex
 */
void RealImag2Complex(const MultidimArray< double > & real,
                      const MultidimArray< double > & imag,
                      MultidimArray< std::complex < double > > & out);
//@}

/** @name Fourier Transforms
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
 * #include <data/xmipp_fft.h>
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
 *       MultidimArray<double> x(65);
 *   x.setXmippOrigin();
 *   double T=0.5;
 *   double T1=6;
 *   int N1=(int)CEIL(T1/T);
 *
 *   // Fill x with a pulse from -N1 to N1 (-T1 to T1 in continuous)
 *   FOR_ALL_ELEMENTS_IN_ARRAY1D(x)
 *      if (ABS(i)<=N1) x(i)=1;
 *
 *   // Compute the Fourier transform
 *   MultidimArray< std::complex<double> > X;
 *   MultidimArray<double> Xmag;
 *   FourierTransform(x,X);
 *   FFT_magnitude(X,Xmag);
 *
 *   // Compute the frequency axes
 *   MultidimArray<double> contfreq(XSIZE(X)), digfreq(XSIZE(X));
 *          FOR_ALL_ELEMENTS_IN_ARRAY1D(X)
 *              FFT_IDX2DIGFREQ(i,XSIZE(X),digfreq(i));
 *   digfreq*=2*PI;
 *   contfreq=digfreq/T;
 *
 *   // Show all Fourier transforms
 *   FOR_ALL_ELEMENTS_IN_ARRAY1D(X) {
 *       if (digfreq(i)>=0)
 *                 std::cout << digfreq(i) << " " << contfreq(i) << " "
 *             << XSIZE(X)*Xmag(i) << " "
 *      << ABS(discreteTransform(digfreq(i),N1)) << " "
 *      << ABS(continuousTransform(contfreq(i),T1)/T)
 *      << std::endl;
 *   }
 *     } catch (XmippError XE) {
 *         std::cout << XE << std::endl;
 *     }
 *     return 0;
 * }
 * @endcode
 */
//@{
/** Direct Fourier Transform
 */
void FourierTransform(const MultidimArray< double >& in,
                      MultidimArray< std::complex< double > > & out);

/** Inverse Fourier Transform
 */
void InverseFourierTransform(const MultidimArray< std::complex< double > > & in,
                             MultidimArray< double >& out);

/** Direct Fourier Transform, output half of (centro-symmetric) transform
 */
void FourierTransformHalf(const MultidimArray< double >& in,
                          MultidimArray< std::complex< double > > & out);

/** Inverse Fourier Transform 1D, input half of (centro-symmetric) transform
 */
void InverseFourierTransformHalf(const MultidimArray< std::complex< double > > & in,
                                 MultidimArray< double >& out,
                                 int oridim);
//@}

/** @name Operations with the Fourier Transforms
 */
//@{

/** Faster version of CenterFFT (now just for even images)
 */

#define SWAP_ARRAY(a, b, n)  memcpy(buffer, a, n); memcpy(a, b, n); memcpy(b, buffer, n);

/** Center FFT for 2D arrays.
 * The function is optimized for the particular case of 2D.
 */
void centerFFT2(MultidimArray<double> &v);


/** CenterFFT
 * Relation with Matlab fftshift: forward true is equals to fftshift and forward false
 * equals to ifftshift
 */
template <typename T>
void CenterFFT(MultidimArray< T >& v, bool forward)
{
	bool firstTime=true;						// First time executing inner loops.

	// Check dimension is between 1 and 3 inclusive.
    if ( v.getDim() > 0 && v.getDim() <= 3)
    {
        // Shift in the X direction
    	size_t l=0;

        // Shift in the X direction
        if (((l = XSIZE(v)) > 1) && (YSIZE(v) == 1) && (ZSIZE(v) == 1))
        {
        	size_t firstHalfSize=0;
        	size_t secondHalfSize=0;
        	MultidimArray< T > aux;

        	if (forward)
        	{
        		firstHalfSize  = (l + 1) / 2;
        		secondHalfSize = l - firstHalfSize;
        		aux.resizeNoCopy( firstHalfSize);

                for (size_t k = 0; k < ZSIZE(v); k++)
                {
                    for (size_t i = 0; i < YSIZE(v); i++)
                    {
                    	memcpy( &dAi(aux, 0), &dAkij(v, k, i, 0), sizeof(T)*firstHalfSize);
                        memcpy( &dAkij(v, k, i, 0), &dAkij(v, k, i, firstHalfSize), sizeof(T)*secondHalfSize);
                        memcpy( &dAkij(v, k, i, secondHalfSize), &dAi(aux, 0), sizeof(T)*firstHalfSize);
                    }
                }
        	}
        	else
        	{
        		secondHalfSize = (l + 1) / 2;
        		firstHalfSize  = l - secondHalfSize;
        		aux.resizeNoCopy( secondHalfSize);

                for (size_t k = 0; k < ZSIZE(v); k++)
                {
                    for (size_t i = 0; i < YSIZE(v); i++)
                    {
                        memcpy( &dAi(aux, 0), &dAkij(v, k, i, firstHalfSize), sizeof(T)*secondHalfSize);
                        memcpy( &dAkij(v, k, i, secondHalfSize), &dAkij(v, k, i, 0), sizeof(T)*firstHalfSize);
                        memcpy( &dAkij(v, k, i, 0), &dAi(aux, 0), sizeof(T)*secondHalfSize);
                    }
                }
        	}
        }
        else
        {
			// 3D
			MultidimArray< T > aux;
			size_t 	l=0;
			long int shift;

			int		i=0;							// Loop counter.
			int		halfRows=0;						// Half rows of the matrix.
			int		rowSize=0;						// Size in bytes of the row.

			// Size in bytes of first and second half row.
			int		firstHalfRowSize=0, secondHalfRowSize=0;

			// # elements in first and second half row.
			int		nElemsFirstHalf_X=0, nElemsSecondHalf_X=0;

			MultidimArray< T >	tempVector;			// Temporary vector.

			bool 	isOdd=false;					// Odd # rows in matrix.
			MultidimArray< T >	savedRow;			// Temporary vector in odd matrix.

			// Execute FFT in 2 dimensions before 3D.
			if (forward)
			{
				for (size_t k=0; k<ZSIZE(v); k++)
				{
					// Pointers to upper and lower half matrix.
					T		*upperHalfPtr, *lowerHalfSrcPtr, *lowerHalfDestPtr;

					// This is computed first time only.
					if (firstTime)
					{
						firstTime = false;

						// Compute # elements in each half vector in X-dimension.
						nElemsFirstHalf_X = (XSIZE(v) + 1) / 2;
						nElemsSecondHalf_X = XSIZE(v) - nElemsFirstHalf_X;

						// Compute X dimension size in bytes.
						rowSize = XSIZE(v)*sizeof(T);

						// Allocate temporary vector.
						firstHalfRowSize = nElemsFirstHalf_X*sizeof(T);
						secondHalfRowSize = rowSize - firstHalfRowSize;
						tempVector.resizeNoCopy( XSIZE(v));

						// Compute # iterations.
						halfRows = YSIZE(v) / 2;

						// If even # rows then save row located in the middle of the matrix.
						if ((YSIZE(v) % 2) != 0)
						{
							isOdd = true;
							savedRow.resizeNoCopy( XSIZE(v));
						}
					}

					// Initialize pointers to upper and lower matrix rows.
					upperHalfPtr     = &dAkij(v, k, 0, 0);
					lowerHalfSrcPtr  = upperHalfPtr + (halfRows*XSIZE(v));
					lowerHalfDestPtr = lowerHalfSrcPtr;

					// If even # rows then save row located in the middle of the matrix.
					if (isOdd)
					{
						memcpy( &dAi(savedRow,0), lowerHalfSrcPtr, rowSize);

						// Source and destination rows are not the same in odd matrix.
						lowerHalfSrcPtr = lowerHalfSrcPtr + XSIZE(v);
					}

					// Half of the rows must be exchanged.
					for (i=0; i<halfRows ;i++)
					{
						memcpy( &dAi(tempVector,0), upperHalfPtr, rowSize);
						memcpy( upperHalfPtr, lowerHalfSrcPtr + nElemsFirstHalf_X, secondHalfRowSize);
						memcpy( upperHalfPtr + nElemsSecondHalf_X, lowerHalfSrcPtr, firstHalfRowSize);
						memcpy( lowerHalfDestPtr, &dAi(tempVector, nElemsFirstHalf_X), secondHalfRowSize);
						memcpy( lowerHalfDestPtr + nElemsSecondHalf_X, &dAi(tempVector,0), firstHalfRowSize);

						upperHalfPtr 	 = upperHalfPtr + XSIZE(v);
						lowerHalfSrcPtr  = lowerHalfSrcPtr + XSIZE(v);
						lowerHalfDestPtr = lowerHalfDestPtr + XSIZE(v);
					}

					// If # rows is odd then restore last row.
					if (isOdd)
					{
						memcpy( lowerHalfDestPtr + nElemsSecondHalf_X, &dAi(savedRow,0), firstHalfRowSize);
						memcpy( lowerHalfDestPtr, &dAi(savedRow, nElemsFirstHalf_X), secondHalfRowSize);
					}
				}
			}
			else
			{
				for (size_t k = 0; k<ZSIZE(v); k++)
				{
					// Pointers to upper and lower half matrix.
					T		*upperHalfSrcPtr, *upperHalfDestPtr, *lowerHalfPtr;

					// This is computed first time only.
					if (firstTime)
					{
						firstTime = false;

						// Compute # elements in each half vector in X-dimension.
						nElemsSecondHalf_X = (XSIZE(v) + 1) / 2;
						nElemsFirstHalf_X  = XSIZE(v) - nElemsSecondHalf_X;

						// Allocate temporary vector.
						rowSize = XSIZE(v)*sizeof(T);
						firstHalfRowSize = nElemsFirstHalf_X*sizeof(T);
						secondHalfRowSize = rowSize - firstHalfRowSize;
						tempVector.resizeNoCopy(rowSize);

						// Compute # iterations.
						halfRows = YSIZE(v) / 2;

						// If odd # rows then save last row.
						if ((YSIZE(v) % 2) != 0)
						{
							isOdd = true;
							savedRow.resizeNoCopy(XSIZE(v));
						}
					}

					// Initialize pointers to quadrants.
					upperHalfSrcPtr  = &dAkij(v, k, halfRows-1, 0);
					lowerHalfPtr  	 = &dAkij(v, k, YSIZE(v)-1, 0);
					upperHalfDestPtr = upperHalfSrcPtr;

					// If odd # rows then save last row.
					if (isOdd)
					{
						// Source and destination rows are not the same in odd matrix.
						upperHalfDestPtr = upperHalfDestPtr + XSIZE(v);

						memcpy( &dAi(savedRow,0), upperHalfDestPtr, rowSize);
					}

					// Half of the rows must be exchanged.
					for (i=0; i<halfRows ;i++)
					{
						memcpy( &dAi(tempVector, 0), upperHalfSrcPtr, rowSize);
						memcpy( upperHalfDestPtr, lowerHalfPtr + nElemsFirstHalf_X, secondHalfRowSize);
						memcpy( upperHalfDestPtr + nElemsSecondHalf_X, lowerHalfPtr, firstHalfRowSize);
						memcpy( lowerHalfPtr, &dAi(tempVector, nElemsFirstHalf_X), secondHalfRowSize);
						memcpy( lowerHalfPtr + nElemsSecondHalf_X, &dAi(tempVector,0), firstHalfRowSize);

						upperHalfSrcPtr  = upperHalfSrcPtr - XSIZE(v);
						lowerHalfPtr     = lowerHalfPtr - XSIZE(v);
						upperHalfDestPtr = upperHalfDestPtr - XSIZE(v);
					}

					// If # rows is odd then restore row at medium position.
					if (isOdd)
					{
						memcpy( upperHalfDestPtr + nElemsSecondHalf_X, &dAi(savedRow, 0), firstHalfRowSize);
						memcpy( upperHalfDestPtr, &dAi(savedRow, nElemsFirstHalf_X), secondHalfRowSize);
					}
				}
			}

			// Shift in the Z direction
			if ((l = ZSIZE(v)) > 1)
			{
				aux.resizeNoCopy(l);
				shift = (long int)(l / 2);

				if (!forward)
					shift = -shift;
				size_t lmax=(l/4)*4;
				for (size_t i = 0; i < YSIZE(v); i++)
				{
					for (size_t j = 0; j < XSIZE(v); j++)
					{
						// Shift the input in an auxiliary vector
						for (size_t k = 0; k < l; k++)
						{
							size_t kp = k + shift;
							if (-shift > (long int)k)
								kp += l;
							else if (kp >= l)
								kp -= l;

							dAi(aux,kp) = dAkij(v, k, i, j);
						}

						// Copy the vector
						const T* ptrAux=&dAi(aux,0);
						for (size_t k = 0; k < lmax; k+=4,ptrAux+=4)
						{
							dAkij(v, k,   i, j) = *ptrAux;
							dAkij(v, k+1, i, j) = *(ptrAux+1);
							dAkij(v, k+2, i, j) = *(ptrAux+2);
							dAkij(v, k+3, i, j) = *(ptrAux+3);
						}
						for (size_t k = lmax; k < l; ++k, ++ptrAux)
						{
							dAkij(v, k, i, j) = *ptrAux;
						}
					}
				}
			}
        }
    }
    else
    {
        REPORT_ERROR(ERR_MULTIDIM_DIM,"CenterFFT ERROR: Dimension should be 1, 2 or 3");
    }
}

/** FFT shift 1D
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(MultidimArray< std::complex< double > > & v, double xshift);

/** FFT shift 2D
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(MultidimArray< std::complex< double > > & v, double xshift, double yshift);

/** FFT shift 3D
 *
 * Calculates the Fourier Transform of the shifted real-space vector
 * by phase shifts in Fourier space
 */
void ShiftFFT(MultidimArray< std::complex< double > > & v,
              double xshift,
              double yshift,
              double zshift);

/** Place the origin of the FFT at the center of the vector and back
 *
 * Changes the real and the fourier space origin
 */
void CenterOriginFFT(MultidimArray< std::complex< double > > & v, bool forward);


/** Xmipp image -> Xmipp PSD.
    The log10 is taken, outliers rejected and the image is reorganized. */
void xmipp2PSD(const MultidimArray<double> &input, MultidimArray<double> &output,
               bool takeLog=true);

/** Xmipp image -> Xmipp CTF.
    The log10 is taken, outliers rejected and the image is reorganized. */
void xmipp2CTF(const MultidimArray<double> &input, MultidimArray<double> &output);
//@}
//@}
#endif
