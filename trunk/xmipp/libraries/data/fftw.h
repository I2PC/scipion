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
                    FOR_ALL_ELEMENTS_IN_MATRIX1D(V)
                        if (i<XSIZE(fFourier))
                            DIRECT_VEC_ELEM(V,i)=DIRECT_VEC_ELEM(fFourier,i);
                        else
                            DIRECT_VEC_ELEM(V,i)=
                                conj(DIRECT_VEC_ELEM(fFourier,
                                    XSIZE(*fReal)-i));
                    break;
                case 2:
                    FOR_ALL_ELEMENTS_IN_MATRIX2D(V)
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
                    FOR_ALL_ELEMENTS_IN_MATRIX3D(V)
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

/** Fourier-Ring-Correlation between two 2D-matrices using FFT
 * @ingroup FourierOperations
 */
void frc_dpr(Matrix2D< double > & m1,
             Matrix2D< double > & m2,
             double sampling_rate,
             Matrix1D< double >& freq,
             Matrix1D< double >& frc,
             Matrix1D< double >& frc_noise,
             Matrix1D< double >& dpr);

/** Fourier-Ring-Correlation between two 3D-matrices using FFT
 * @ingroup FourierOperations
 */
void frc_dpr(Matrix3D< double > & m1,
             Matrix3D< double > & m2,
             double sampling_rate,
             Matrix1D< double >& freq,
             Matrix1D< double >& frc,
             Matrix1D< double >& frc_noise,
             Matrix1D< double >& dpr);
#endif
