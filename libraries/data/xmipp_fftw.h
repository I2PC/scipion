/***************************************************************************
 *
 * Authors:    Roberto Marabini                 (roberto@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef __XmippFFTW_H
#define __XmippFFTW_H

#include <complex>
#include "../../external/fftw-3.3.1/api/fftw3.h"
#include "multidim_array.h"
#include "multidim_array_generic.h"
#include "xmipp_fft.h"


/** @defgroup FourierW FFTW Fourier transforms
  * @ingroup DataLibrary
  *@{
  */

/** Fourier Transformer class.
 * @ingroup FourierW
 *
 * The memory for the Fourier transform is handled by this object.
 * However, the memory for the real space image is handled externally
 * and this object only has a pointer to it.
 *
 * Here you have an example of use
 * @code
 * FourierTransformer transformer;
 * MultidimArray< std::complex<double> > Vfft;
 * transformer.FourierTransform(V(),Vfft,false);
 * MultidimArray<double> Vmag;
 * Vmag.resize(Vfft);
 * FOR_ALL_ELEMENTS_IN_ARRAY3D(Vmag)
 *     Vmag(k,i,j)=20*log10(abs(Vfft(k,i,j)));
 * @endcode
 */
class FourierTransformer
{
public:
    /** Real array, in fact a pointer to the user array is stored. */
    MultidimArray<double> *fReal;

    /** Complex array, in fact a pointer to the user array is stored. */
    MultidimArray<std::complex<double> > *fComplex;

    /** Fourier array  */
    MultidimArray< std::complex<double> > fFourier;

    /* fftw Forawrd plan */
    fftw_plan fPlanForward;

    /* fftw Backward plan */
    fftw_plan fPlanBackward;

    /* number of threads*/
    int nthreads;

    /* Threads has been used in this program*/
    bool threadsSetOn;

    /* Sign where the normalization is applied */
    int normSign;

    // Public methods
public:
    /** Default constructor */
    FourierTransformer();

    /* Constructor setting the sign of normalization application*/
    FourierTransformer(int _normSign);
    /** Destructor */
    ~FourierTransformer();

    /** Set Number of threads
     * This function, which should be called once, performs any
     * one-time initialization required to use threads on your
     * system.
     *
     *  The nthreads argument indicates the number of threads you
     *  want FFTW to use (or actually, the maximum number). All
     *  plans subsequently created with any planner routine will use
     *  that many threads. You can call fftw_plan_with_nthreads,
     *  create some plans, call fftw_plan_with_nthreads again with a
     *  different argument, and create some more plans for a new
     *  number of threads. Plans already created before a call to
     *  fftw_plan_with_nthreads are unaffected. If you pass an
     *  nthreads argument of 1 (the default), threads are
     *  disabled for subsequent plans. */
    void setThreadsNumber(int tNumber)
    {
        if (tNumber!=1)
        {
            threadsSetOn=true;
            nthreads = tNumber;
            if(fftw_init_threads()==0)
                REPORT_ERROR(ERR_THREADS_NOTINIT, (std::string)"FFTW cannot init threads (setThreadsNumber)");
            fftw_plan_with_nthreads(nthreads);
        }
    }
    /** Change Number of threads.
     *
     *  The nthreads argument indicates the number of threads you want FFTW to use
     *  (or actually, the maximum number). All plans subsequently
     *  created with any planner routine will use that many
     *  threads. You can call fftw_plan_with_nthreads, create
     *  some plans, call fftw_plan_with_nthreads again with a
     *  different argument, and create some more plans for a new
     *  number of threads. Plans already created before a call to
     *  fftw_plan_with_nthreads are unaffected. If you pass an
     *  nthreads argument of 1 (the default), threads are
     *  disabled for subsequent plans. */
    void changeThreadsNumber(int tNumber)
    {
        nthreads = tNumber;
        fftw_plan_with_nthreads(nthreads);
    }

    /** Destroy Threads. Do not execute any previously created
     *  plans after calling this function   */
    void destroyThreads(void )
    {
        nthreads = 1;
        if(threadsSetOn)
            fftw_cleanup_threads();

        threadsSetOn=false;
    }

    /** Compute the Fourier transform of a MultidimArray, 2D and 3D.
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
        if (getCopy)
            getFourierCopy(V);
        else
            getFourierAlias(V);
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
    void getFourierAlias(T& V)
    {
        V.alias(fFourier);
        return;
    }

    /** Get Fourier coefficients. */
    template <typename T>
    void getFourierCopy(T& V)
    {
        V.resizeNoCopy(fFourier);
        memcpy(MULTIDIM_ARRAY(V),MULTIDIM_ARRAY(fFourier),
               MULTIDIM_SIZE(fFourier)*2*sizeof(double));
    }

    /** Return a complete Fourier transform (two halves).
    */
    template <typename T>
    void getCompleteFourier(T& V)
    {
        V.resizeNoCopy(*fReal);
        int ndim=3;
        if (ZSIZE(*fReal)==1)
        {
            ndim=2;
            if (YSIZE(*fReal)==1)
                ndim=1;
        }
        double *ptrSource=NULL;
        double *ptrDest=NULL;
        switch (ndim)
        {
        case 1:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(V)
            {
                ptrDest=(double*)&DIRECT_A1D_ELEM(V,i);
                if (i<XSIZE(fFourier))
                {
                    ptrSource=(double*)&DIRECT_A1D_ELEM(fFourier,i);
                    *ptrDest=*ptrSource;
                    *(ptrDest+1)=*(ptrSource+1);
                }
                else
                {
                    ptrSource=(double*)&DIRECT_A1D_ELEM(fFourier,XSIZE(*fReal)-i);
                    *ptrDest=*ptrSource;
                    *(ptrDest+1)=-(*(ptrSource+1));
                }
            }
            break;
        case 2:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(V)
            {
                ptrDest=(double*)&DIRECT_A2D_ELEM(V,i,j);
                if (j<XSIZE(fFourier))
                {
                    ptrSource=(double*)&DIRECT_A2D_ELEM(fFourier,i,j);
                    *ptrDest=*ptrSource;
                    *(ptrDest+1)=*(ptrSource+1);
                }
                else
                {
                    ptrSource=(double*)&DIRECT_A2D_ELEM(fFourier,
                                                        (YSIZE(*fReal)-i)%YSIZE(*fReal),
                                                        XSIZE(*fReal)-j);
                    *ptrDest=*ptrSource;
                    *(ptrDest+1)=-(*(ptrSource+1));
                }
            }
            break;
        case 3:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V)
            {
                ptrDest=(double*)&DIRECT_A3D_ELEM(V,k,i,j);
                if (j<XSIZE(fFourier))
                {
                    ptrSource=(double*)&DIRECT_A3D_ELEM(fFourier,k,i,j);
                    *ptrDest=*ptrSource;
                    *(ptrDest+1)=*(ptrSource+1);
                }
                else
                {
                    ptrSource=(double*)&DIRECT_A3D_ELEM(fFourier,
                                                        (ZSIZE(*fReal)-k)%ZSIZE(*fReal),
                                                        (YSIZE(*fReal)-i)%YSIZE(*fReal),
                                                        XSIZE(*fReal)-j);
                    *ptrDest=*ptrSource;
                    *(ptrDest+1)=-(*(ptrSource+1));
                }
            }
            break;
        }
    }


    /** Compute the Fourier transform of a MultidimArray, 2D and 3D and
     * returns a complete copy.
        */
    template <typename T, typename T1>
    void completeFourierTransform(T& v, T1& V)
    {
        setReal(v);
        Transform(FFTW_FORWARD);
        getCompleteFourier(V);
    }

    /** Set one half of the FT in fFourier from the input complete Fourier transform (two halves).
        The fReal and fFourier already should have the right sizes
    */
    template <typename T>
    void setFromCompleteFourier(T& V)
    {
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
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fFourier)
            DIRECT_A1D_ELEM(fFourier,i)=DIRECT_A1D_ELEM(V,i);
            break;
        case 2:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(fFourier)
            DIRECT_A2D_ELEM(fFourier,i,j) = DIRECT_A2D_ELEM(V,i,j);
            break;
        case 3:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(fFourier)
            DIRECT_A3D_ELEM(fFourier,k,i,j) = DIRECT_A3D_ELEM(V,k,i,j);
            break;
        }
    }

    // Internal methods
public:
    /* Pointer to the array of doubles with which the plan was computed */
    double * dataPtr;

    /* Pointer to the array of complex<double> with which the plan was computed */
    std::complex<double> * complexDataPtr;

    /* Init object*/
    void init();
    /** Clear object */
    void clear();
    /**     FFTW's planner saves some other persistent data,
     * such as the accumulated wisdom and a list of algorithms available
     * in the current configuration. If you want to deallocate all of that
     * and reset FFTW to the pristine state it was in when
     * you started your program, you can call:
     */
    void cleanup(void)
    {
        fftw_cleanup();
    }
    /** Computes the transform, specified in Init() function
        If normalization=true the forward transform is normalized
        (no normalization is made in the inverse transform)
        If normalize=false no normalization is performed and therefore
        the image is scaled by the number of pixels.
    */
    void Transform(int sign);

    /** Get the Multidimarray that is being used as input. */
    const MultidimArray<double> &getReal() const;
    const MultidimArray<std::complex<double> > &getComplex() const;

    /** Set a Multidimarray for input.
        The data of img will be the one of fReal. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<double> &img);

    /** Set a Multidimarray for input.
        The data of img will be the one of fComplex. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<std::complex<double> > &img);

    /** Set a Multidimarray for the Fourier transform.
        The values of the input array are copied in the internal array.
        It is assumed that the container for the real image as well as
        the one for the Fourier array are already resized.
        No plan is updated. */
    void setFourier(const MultidimArray<std::complex<double> > &imgFourier);

    /* Set normalization sign.
     * It defines when the normalization must be applied, when doing
     * FFTW_FORWARD OR FFTW_BACKWARD. By default, FFTW_FORWARD.*/
    void setNormalizationSign(int _normSign)
    {
        normSign = _normSign;
    }

};

/** FFT Magnitude 1D
 * @ingroup FourierOperations
 */
void FFT_magnitude(const MultidimArray< std::complex< double > > & v,
                   MultidimArray< double >& mag);

/** FFT Phase 1D
 * @ingroup FourierOperations
 */
void FFT_phase(const MultidimArray< std::complex< double > > & v,
               MultidimArray< double >& phase);

/** Autocorrelation function of a Xmipp vector
 * @ingroup FourierOperations
 *
 * Fast calculation of the autocorrelation vector of a given one using Fast
 * Fourier Transform. (Using the correlation theorem)
 */
template <typename T>
void auto_correlation_vector(const MultidimArray< T > & Img, MultidimArray< double >& R)
{
    Img.checkDimension(1);

    // Compute the Fourier Transform
    MultidimArray< std::complex< double > > FFT1;
    FourierTransformer transformer1;
    R=Img;
    transformer1.FourierTransform(R, FFT1, false);

    // Multiply FFT1 * FFT1'
    double dSize=XSIZE(Img);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FFT1)
    {
        double *ptr=(double*)&DIRECT_MULTIDIM_ELEM(FFT1,n);
        double &realPart=*ptr;
        double &imagPart=*(ptr+1);
        realPart=dSize*(realPart*realPart+imagPart*imagPart);
        imagPart=0;
    }

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
void correlation_vector(const MultidimArray< T > & m1,
                        const MultidimArray< T > & m2,
                        MultidimArray< double >& R)
{
    m1.checkDimension(2);
    m2.checkDimension(2);

    // Compute the Fourier Transforms
    MultidimArray< std::complex< double > > FFT1, FFT2;
    FourierTransformer transformer1, transformer2;
    R=m1;
    transformer1.FourierTransform(R, FFT1, false);
    transformer2.FourierTransform((MultidimArray<T> &)m2, FFT2, false);

    // Multiply FFT1 * FFT2'
    double dSize=XSIZE(m1);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(FFT1)
    FFT1(i) *= dSize * conj(FFT2(i));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Autocorrelation function of a Xmipp vector
 * @ingroup FourierOperations
 *
 * Same as correlation_vector, but the Fourier transforms have already been computed
 * and the transformer is reused. The transformer internal variables must have already
 * been resized.
 */
void fast_correlation_vector(const MultidimArray< std::complex<double> > & FFT1,
                        const MultidimArray< std::complex<double> > & FFT2,
                        MultidimArray< double >& R,
                        FourierTransformer &transformer);

/** Compute the correlation vector without using Fourier.
 * @ingroup FourierOperations
 *
 *  The results are the same as the previous ones but this function
 *  is threadsafe while the previous one is not.
 *
 *  It is assumed that the two vectors v1, and v2 are of the same size. */
template <class T>
void correlation_vector_no_Fourier(const MultidimArray<T> &v1, const MultidimArray<T> &v2,
                                   MultidimArray<T> &result)
{
    v1.checkDimension(1);
    v2.checkDimension(1);

    result.initZeros(v1);
    result.setXmippOrigin();
    int N=XSIZE(v1)-1;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(result)
    for (int k=0; k<XSIZE(v1); ++k)
        A1D_ELEM(result,i)+=DIRECT_A1D_ELEM(v1,intWRAP(k+i,0,N))*
                            DIRECT_A1D_ELEM(v2,k);
    STARTINGX(result)=0;
}

/** Correlation auxiliary. */
class CorrelationAux
{
public:
    MultidimArray< std::complex< double > > FFT1, FFT2;
    FourierTransformer transformer1, transformer2;
};

/** Correlation of two nD images
 * @ingroup FourierOperations
 *
 * Fast calcuation of the correlation matrix on two matrices using Fast Fourier
 * Transform. (Using the correlation theorem). The output matrix must be already
 * resized
 */

void correlation_matrix(const MultidimArray<double> & m1,
                        const MultidimArray<double> & m2,
                        MultidimArray< double >& R,
                        CorrelationAux &aux,
                        bool center=true);

void correlation_matrix(const MultidimArray< std::complex< double > > & FFT1,
                        const MultidimArray<double> & m2,
                        MultidimArray<double>& R,
                        CorrelationAux &aux,
                        bool center=true);

/** Autocorrelation function of an image
 * @ingroup FourierOperations
 *
 * Fast calcuation of the autocorrelation matrix of a given one using Fast
 * Fourier Transform. (Using the correlation theorem)
 */
template <typename T>
void auto_correlation_matrix(const MultidimArray< T > & Img, MultidimArray< double >& R)
{
    // Compute the Fourier Transform
    MultidimArray< std::complex< double > > FFT1;
    FourierTransformer transformer1;
    R=Img;
    transformer1.FourierTransform(R, FFT1, false);

    // Multiply FFT1 * FFT1'
    double dSize=MULTIDIM_SIZE(Img);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FFT1)
    {
        double *ptr=(double*)&DIRECT_MULTIDIM_ELEM(FFT1,n);
        double &realPart=*ptr;
        double &imagPart=*(ptr+1);
        realPart=dSize*(realPart*realPart+imagPart*imagPart);
        imagPart=0;
    }

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}


void convolutionFFT(const MultidimArray<double> &img,
                    const MultidimArray<double> &kernel,
                    MultidimArray<double> &result);

/** Fourier-Ring-Correlation between two multidimArrays using FFT
 * @ingroup FourierOperations
 */
void frc_dpr(MultidimArray< double > & m1,
             MultidimArray< double > & m2,
             double sampling_rate,
             MultidimArray< double >& freq,
             MultidimArray< double >& frc,
             MultidimArray< double >& frc_noise,
             MultidimArray< double >& dpr,
             MultidimArray< double >& error_l2,
             bool skipdpr=false);

/** Scale matrix using Fourier transform
 * @ingroup FourierOperations
 * Ydim and Xdim define the output size, Mpmem is the matrix to scale
 */
void selfScaleToSizeFourier(int Zdim, int Ydim, int Xdim, MultidimArray<double> &Mpmem, int nthreads=1);

void selfScaleToSizeFourier(int Ydim, int Xdim, MultidimArray<double> &Mpmem, int nthreads=1);
/** MultidimArrayGeneric version */
void selfScaleToSizeFourier(int Zdim, int Ydim, int Xdim, MultidimArrayGeneric &Mpmem, int nthreads=1);
void selfScaleToSizeFourier(int Ydim, int Xdim, MultidimArrayGeneric &Mpmem, int nthreads=1);

#define POWER_SPECTRUM 0
#define AMPLITUDE_SPECTRUM 1

/** Get the amplitude or power spectrum of the map in Fourier space.
 * @ingroup FourierOperations
    i.e. the radial average of the (squared) amplitudes of all Fourier components
*/
void getSpectrum(MultidimArray<double> &Min,
                 MultidimArray<double> &spectrum,
                 int spectrum_type=AMPLITUDE_SPECTRUM);

/** Divide the input map in Fourier-space by the spectrum provided.
 * @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void divideBySpectrum(MultidimArray<double> &Min,
                      MultidimArray<double> &spectrum,
                      bool leave_origin_intact=false);

/** Multiply the input map in Fourier-space by the spectrum provided.
 * @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void multiplyBySpectrum(MultidimArray<double> &Min,
                        MultidimArray<double> &spectrum,
                        bool leave_origin_intact=false);

/** Perform a whitening of the amplitude/power spectrum of a 3D map
 * @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void whitenSpectrum(MultidimArray<double> &Min,
                    MultidimArray<double> &Mout,
                    int spectrum_type=AMPLITUDE_SPECTRUM,
                    bool leave_origin_intact=false);

/** Adapts Min to have the same spectrum as spectrum_ref
 * @ingroup FourierOperations
    If only_amplitudes==true, the amplitude rather than the power spectrum will be equalized
*/
void adaptSpectrum(MultidimArray<double> &Min,
                   MultidimArray<double> &Mout,
                   const MultidimArray<double> &spectrum_ref,
                   int spectrum_type=AMPLITUDE_SPECTRUM,
                   bool leave_origin_intact=false);

/** Randomize phases beyond a certain frequency
 * @ingroup FourierOperations
*/
void randomizePhases(MultidimArray<double> &Min, double w);
#endif
/** @} */
