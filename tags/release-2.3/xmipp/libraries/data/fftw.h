/***************************************************************************
 *
 * Authors:    Roberto Marabini                 (roberto@cnb.uam.es)
 *             Carlos Oscar S. Sorzano          (coss@cnb.csic.es)
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

#ifndef __XmippFFTW_H
#define __XmippFFTW_H

#include <complex>
#include "../../external/fftw-3.1.2/api/fftw3.h"
#include "matrix3d.h"
#include "fft.h"

/** @defgroup FourierW FFTW Fourier transforms
  * @ingroup DataLibrary
  */
 
/** Fourier Transformer class.
 * @ingroup FourierW
 *
 * The memory for the Fourier transform is handled by this object.
 * However, the memory for the real space image is handled externally
 * and this object only has a pointer to it.
 */
class XmippFftw
{
public:
    /** Real array, in fact a pointer to the user array is stored. */
    MultidimArray<double> *fReal;

    /** Fourier array  */
    MultidimArray< std::complex<double> > fFourier;

    /* fftw Forawrd plan */
    fftw_plan fPlanForward;

    /* fftw Backward plan */
    fftw_plan fPlanBackward;
// Public methods
public:
    /** Default constructor */
    XmippFftw();

    /** Destructor */
    ~XmippFftw();

    /** Compute the Fourier transform of a Matrix1D, 2D and 3D.
        If getCopy is false, an alias to the transformed data is returned.
        This is a faster option since a copy of all the data is avoided,
        but you need to be careful that an inverse Fourier transform may
        change the data.
        */
    template <typename T, typename T1>
        void FourierTransform(T& v, T1& V, bool getCopy=true)
        {
            setReal(v);
            Transform(FFTW_FORWARD);
            if (getCopy) getFourierCopy(V);
            else         getFourierAlias(V);
        }

    /** Compute the Fourier transform.
        The data is taken from the matrix with which the object was
        created. */
    void FourierTransform();

    /** Inforce Hermitian symmetry.
        If the Fourier transform risks of losing Hermitian symmetry,
        use this function to renforce it. */
    void enforceHermitianSymmetry();

    /** Compute the inverse Fourier transform.
        The result is stored in the same real data that was passed for
        the forward transform. The Fourier coefficients are taken from
        the internal Fourier coefficients */
    void inverseFourierTransform();

    /** Compute the inverse Fourier transform.
        New data is provided for the Fourier coefficients and the output
        can be any matrix1D, 2D or 3D. It is important that the output
        matrix is already resized to the right size before entering
        in this function. */
    template <typename T, typename T1>
        void inverseFourierTransform(T& V, T1& v)
        {
            setReal(v);
            setFourier(V);
            Transform(FFTW_BACKWARD);
        }

    /** Get Fourier coefficients. */
    template <typename T>
        void getFourierAlias(T& V) {V.alias(fFourier); return;}

    /** Get Fourier coefficients. */
    template <typename T>
        void getFourierCopy(T& V) {
            V.resize(fFourier);
            memcpy(MULTIDIM_ARRAY(V),MULTIDIM_ARRAY(fFourier),
                MULTIDIM_SIZE(fFourier)*2*sizeof(double));
        }

    /** Return a complete Fourier transform (two halves).
    */
    template <typename T>
        void getCompleteFourier(T& V) {
            V.resize(*fReal);
            int ndim=3;
            if (ZSIZE(*fReal)==1)
            {
                ndim=2;
                if (YSIZE(*fReal)==1)
                    ndim=1;
            }
            switch (ndim)
            {
                case 1:
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(V)
                        if (i<XSIZE(fFourier))
                            DIRECT_VEC_ELEM(V,i)=DIRECT_VEC_ELEM(fFourier,i);
                        else
                            DIRECT_VEC_ELEM(V,i)=
                                conj(DIRECT_VEC_ELEM(fFourier,
                                    XSIZE(*fReal)-i));
                    break;
                case 2:
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(V)
                        if (j<XSIZE(fFourier))
                            DIRECT_MAT_ELEM(V,i,j)=
                                DIRECT_MAT_ELEM(fFourier,i,j);
                        else
                            DIRECT_MAT_ELEM(V,i,j)=
                                conj(DIRECT_MAT_ELEM(fFourier,
                                    (YSIZE(*fReal)-i)%YSIZE(*fReal),
                                     XSIZE(*fReal)-j));
                    break;
                case 3:
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(V)
                        if (j<XSIZE(fFourier))
                            DIRECT_VOL_ELEM(V,k,i,j)=
                                DIRECT_VOL_ELEM(fFourier,k,i,j);
                        else
                            DIRECT_VOL_ELEM(V,k,i,j)=
                                conj(DIRECT_VOL_ELEM(fFourier,
                                    (ZSIZE(*fReal)-k)%ZSIZE(*fReal),
                                    (YSIZE(*fReal)-i)%YSIZE(*fReal),
                                     XSIZE(*fReal)-j));
                    break;
            }
        }

// Internal methods
public:
    /* Pointer to the array of doubles with which the plan was computed */
    double * dataPtr;

    /** Clear object */
    void clear();

    /** Computes the transform, specified in Init() function
        If normalization=true the forward transform is normalized
        (no normalization is made in the inverse transform)
        If normalize=false no normalization is performed and therefore
        the image is scaled by the number of pixels.
    */
    void Transform(int sign);

    /** Set a Multidimarray for input.
        The data of img will be the one of fReal. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<double> &img);

    /** Set a Multidimarray for the Fourier transform.
        The values of the input array are copied in the internal array.
        It is assumed that the container for the real image as well as
        the one for the Fourier array are already resized.
        No plan is updated. */
    void setFourier(MultidimArray<std::complex<double> > &imgFourier);
};

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
void auto_correlation_vector(const Matrix1D< T > & Img, Matrix1D< double >& R)
{
    // Compute the Fourier Transform
    Matrix1D< std::complex< double > > FFT1;
    XmippFftw transformer1;
    R=Img;
    transformer1.FourierTransform(R, FFT1, false);

    // Multiply FFT1 * FFT1'
    double dSize=XSIZE(Img);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(FFT1)
        FFT1(i) *= dSize * conj(FFT1(i));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

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
void correlation_vector(const Matrix1D< T > & m1,
                        const Matrix1D< T > & m2,
                        Matrix1D< double >& R)
{
    // Compute the Fourier Transforms
    Matrix1D< std::complex< double > > FFT1, FFT2;
    XmippFftw transformer1, transformer2;
    R=m1;
    transformer1.FourierTransform(R, FFT1, false);
    transformer2.FourierTransform((Matrix1D<T> &)m2, FFT2, false);

    // Multiply FFT1 * FFT2'
    double dSize=XSIZE(m1);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(FFT1)
        FFT1(i) *= dSize * conj(FFT2(i));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Compute the correlation vector without using Fourier.
 * @ingroup FourierOperations
 *
 *  The results are the same as the previous ones but this function
 *  is threadsafe while the previous one is not.
 *  
 *  It is assumed that the two vectors v1, and v2 are of the same size. */
template <class T>
void correlation_vector_no_Fourier(const Matrix1D<T> &v1, const Matrix1D<T> &v2,
    Matrix1D<T> &result)
{
    result.initZeros(v1);
    result.setXmippOrigin();
    int N=XSIZE(v1)-1;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(result)
        for (int k=0; k<XSIZE(v1); k++)
            result(i)+=DIRECT_VEC_ELEM(v1,intWRAP(k+i,0,N))*
                       DIRECT_VEC_ELEM(v2,k);
    STARTINGX(result)=0;
}


/** Correlation of two matrices
 * @ingroup FourierOperations
 *
 * Fast calcuation of the correlation matrix on two matrices using Fast Fourier
 * Transform. (Using the correlation theorem). The output matrix must be already
 * resized
 */
template <typename T>
void correlation_matrix(const Matrix2D< T > & m1,
                        const Matrix2D< T > & m2,
                        Matrix2D< double >& R)
{
    // Compute the Fourier Transforms
    Matrix2D< std::complex< double > > FFT1, FFT2;
    XmippFftw transformer1, transformer2;
    R=m1;
    transformer1.FourierTransform(R, FFT1, false);
    transformer2.FourierTransform((Matrix2D<T> &)m2, FFT2, false);

    // Multiply FFT1 * FFT2'
    double dSize=MULTIDIM_SIZE(R);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT1)
        FFT1(i, j) *= dSize * conj(FFT2(i, j));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

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
void auto_correlation_matrix(const Matrix2D< T > & Img, Matrix2D< double >& R)
{
    // Compute the Fourier Transform
    Matrix2D< std::complex< double > > FFT1;
    XmippFftw transformer1;
    R=Img;
    transformer1.FourierTransform(R, FFT1);

    // Multiply FFT1 * FFT1'
    double dSize=MULTIDIM_SIZE(Img);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT1)
        FFT1(i, j) *= dSize * conj(FFT1(i, j));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Fourier-Ring-Correlation between two 2D-matrices using FFT
 * @ingroup FourierOperations
 */
void frc_dpr(Matrix2D< double > & m1,
             Matrix2D< double > & m2,
             double sampling_rate,
             Matrix1D< double >& freq,
             Matrix1D< double >& frc,
             Matrix1D< double >& frc_noise,
             Matrix1D< double >& dpr,
             bool skipdpr=false);

/** Fourier-Ring-Correlation between two 3D-matrices using FFT
 * @ingroup FourierOperations
 */
void frc_dpr(Matrix3D< double > & m1,
             Matrix3D< double > & m2,
             double sampling_rate,
             Matrix1D< double >& freq,
             Matrix1D< double >& frc,
             Matrix1D< double >& frc_noise,
             Matrix1D< double >& dpr,
             bool skipdpr=false);
#endif
