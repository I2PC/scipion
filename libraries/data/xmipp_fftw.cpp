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

#include "xmipp_fftw.h"
#include "args.h"
#include <string.h>
#include <pthread.h>

static pthread_mutex_t fftw_plan_mutex = PTHREAD_MUTEX_INITIALIZER;

// Constructors and destructors --------------------------------------------
FourierTransformer::FourierTransformer()
{
    init();
    nthreads=1;
    threadsSetOn=false;
    normSign = FFTW_FORWARD;
}

FourierTransformer::FourierTransformer(int _normSign)
{
    init();
    nthreads=1;
    threadsSetOn=false;
    normSign = _normSign;
}

void FourierTransformer::init()
{
    fReal=NULL;
    fComplex=NULL;
    fPlanForward     = NULL;
    fPlanBackward    = NULL;
    dataPtr          = NULL;
    complexDataPtr   = NULL;
}

void FourierTransformer::clear()
{
    fFourier.clear();
    // Anything to do with plans has to be protected for threads!
    pthread_mutex_lock(&fftw_plan_mutex);
    if (fPlanForward !=NULL)
        fftw_destroy_plan(fPlanForward);
    if (fPlanBackward!=NULL)
        fftw_destroy_plan(fPlanBackward);
    pthread_mutex_unlock(&fftw_plan_mutex);
    init();
}

FourierTransformer::~FourierTransformer()
{
    clear();
}

// Initialization ----------------------------------------------------------
const MultidimArray<double> &FourierTransformer::getReal() const
{
    return (*fReal);
}

const MultidimArray<std::complex<double> > &FourierTransformer::getComplex() const
{
    return (*fComplex);
}


void FourierTransformer::setReal(MultidimArray<double> &input)
{
    bool recomputePlan=false;
    if (fReal==NULL)
        recomputePlan=true;
    else if (dataPtr!=MULTIDIM_ARRAY(input))
        recomputePlan=true;
    else
        recomputePlan=!(fReal->sameShape(input));
    fFourier.resizeNoCopy(ZSIZE(input),YSIZE(input),XSIZE(input)/2+1);
    fReal=&input;

    if (recomputePlan)
    {
        int ndim=3;
        if (ZSIZE(input)==1)
        {
            ndim=2;
            if (YSIZE(input)==1)
                ndim=1;
        }
        int N[3];
        switch (ndim)
        {
        case 1:
            N[0]=XSIZE(input);
            break;
        case 2:
            N[0]=YSIZE(input);
            N[1]=XSIZE(input);
            break;
        case 3:
            N[0]=ZSIZE(input);
            N[1]=YSIZE(input);
            N[2]=XSIZE(input);
            break;
        }

        pthread_mutex_lock(&fftw_plan_mutex);
        if (fPlanForward!=NULL)
            fftw_destroy_plan(fPlanForward);
        fPlanForward=NULL;
        fPlanForward = fftw_plan_dft_r2c(ndim, N, MULTIDIM_ARRAY(*fReal),
                                         (fftw_complex*) MULTIDIM_ARRAY(fFourier), FFTW_ESTIMATE);
        if (fPlanBackward!=NULL)
            fftw_destroy_plan(fPlanBackward);
        fPlanBackward=NULL;
        fPlanBackward = fftw_plan_dft_c2r(ndim, N,
                                          (fftw_complex*) MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(*fReal),
                                          FFTW_ESTIMATE);
        if (fPlanForward == NULL || fPlanBackward == NULL)
            REPORT_ERROR(ERR_PLANS_NOCREATE, "FFTW plans cannot be created");
        dataPtr=MULTIDIM_ARRAY(*fReal);
        pthread_mutex_unlock(&fftw_plan_mutex);
    }
}

void FourierTransformer::setReal(MultidimArray<std::complex<double> > &input)
{
    bool recomputePlan=false;
    if (fComplex==NULL)
        recomputePlan=true;
    else if (complexDataPtr!=MULTIDIM_ARRAY(input))
        recomputePlan=true;
    else
        recomputePlan=!(fComplex->sameShape(input));
    fFourier.resizeNoCopy(input);
    fComplex=&input;

    if (recomputePlan)
    {
        int ndim=3;
        if (ZSIZE(input)==1)
        {
            ndim=2;
            if (YSIZE(input)==1)
                ndim=1;
        }
        int *N = new int[ndim];
        switch (ndim)
        {
        case 1:
            N[0]=XSIZE(input);
            break;
        case 2:
            N[0]=YSIZE(input);
            N[1]=XSIZE(input);
            break;
        case 3:
            N[0]=ZSIZE(input);
            N[1]=YSIZE(input);
            N[2]=XSIZE(input);
            break;
        }

        pthread_mutex_lock(&fftw_plan_mutex);
        if (fPlanForward!=NULL)
            fftw_destroy_plan(fPlanForward);
        fPlanForward=NULL;
        fPlanForward = fftw_plan_dft(ndim, N, (fftw_complex*) MULTIDIM_ARRAY(*fComplex),
                                     (fftw_complex*) MULTIDIM_ARRAY(fFourier), FFTW_FORWARD, FFTW_ESTIMATE);
        if (fPlanBackward!=NULL)
            fftw_destroy_plan(fPlanBackward);
        fPlanBackward=NULL;
        fPlanBackward = fftw_plan_dft(ndim, N, (fftw_complex*) MULTIDIM_ARRAY(fFourier),
                                      (fftw_complex*) MULTIDIM_ARRAY(*fComplex), FFTW_BACKWARD, FFTW_ESTIMATE);
        if (fPlanForward == NULL || fPlanBackward == NULL)
            REPORT_ERROR(ERR_PLANS_NOCREATE, "FFTW plans cannot be created");
        delete [] N;
        complexDataPtr=MULTIDIM_ARRAY(*fComplex);
        pthread_mutex_unlock(&fftw_plan_mutex);
    }
}

void FourierTransformer::setFourier(const MultidimArray<std::complex<double> > &inputFourier)
{
    memcpy(MULTIDIM_ARRAY(fFourier),MULTIDIM_ARRAY(inputFourier),
           MULTIDIM_SIZE(inputFourier)*2*sizeof(double));
}

// Transform ---------------------------------------------------------------
void FourierTransformer::Transform(int sign)
{
    if (sign == FFTW_FORWARD)
    {
        fftw_execute(fPlanForward);

        if (sign == normSign)
        {
            unsigned long int size=0;
            if(fReal!=NULL)
                size = MULTIDIM_SIZE(*fReal);
            else if (fComplex!= NULL)
                size = MULTIDIM_SIZE(*fComplex);
            else
                REPORT_ERROR(ERR_UNCLASSIFIED,"No complex nor real data defined");

            double isize=1.0/size;
            double *ptr=(double*)MULTIDIM_ARRAY(fFourier);
            size_t nmax=(fFourier.nzyxdim/4)*4;
            for (size_t n=0; n<nmax; n+=4)
            {
                *ptr++ *= isize;
                *ptr++ *= isize;
                *ptr++ *= isize;
                *ptr++ *= isize;
                *ptr++ *= isize;
                *ptr++ *= isize;
                *ptr++ *= isize;
                *ptr++ *= isize;
            }
            for (size_t n=nmax; n<fFourier.nzyxdim; ++n)
            {
                *ptr++ *= isize;
                *ptr++ *= isize;
            }
        }
    }
    else if (sign == FFTW_BACKWARD)
    {
        fftw_execute(fPlanBackward);

        if (sign == normSign)
        {
            unsigned long int size=0;
            if(fReal!=NULL)
            {
                size = MULTIDIM_SIZE(*fReal);
                double isize=1.0/size;
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*fReal)
                DIRECT_MULTIDIM_ELEM(*fReal,n) *= isize;
            }
            else if (fComplex!= NULL)
            {
                size = MULTIDIM_SIZE(*fComplex);
                double isize=1.0/size;
                double *ptr=(double*)MULTIDIM_ARRAY(*fComplex);
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*fComplex)
                {
                    *ptr++ *= isize;
                    *ptr++ *= isize;
                }
            }
            else
                REPORT_ERROR(ERR_UNCLASSIFIED,"No complex nor real data defined");
        }
    }
}

void FourierTransformer::FourierTransform()
{
    Transform(FFTW_FORWARD);
}

void FourierTransformer::inverseFourierTransform()
{
    Transform(FFTW_BACKWARD);
}

// Inforce Hermitian symmetry ---------------------------------------------
void FourierTransformer::enforceHermitianSymmetry()
{
    int ndim=3;
    int Zdim=ZSIZE(*fReal);
    int Ydim=YSIZE(*fReal);
    if (Zdim==1)
    {
        ndim=2;
        if (Ydim==1)
            ndim=1;
    }
    int yHalf=Ydim/2;
    if (Ydim%2==0)
        yHalf--;
    int zHalf=Zdim/2;
    if (Zdim%2==0)
        zHalf--;
    switch (ndim)
    {
    case 2:
        for (int i=1; i<=yHalf; i++)
        {
            int isym=intWRAP(-i,0,Ydim-1);
            std::complex<double> mean=0.5*(
                                          DIRECT_A2D_ELEM(fFourier,i,0)+
                                          conj(DIRECT_A2D_ELEM(fFourier,isym,0)));
            DIRECT_A2D_ELEM(fFourier,i,0)=mean;
            DIRECT_A2D_ELEM(fFourier,isym,0)=conj(mean);
        }
        break;
    case 3:
        for (int k=0; k<Zdim; k++)
        {
            int ksym=intWRAP(-k,0,Zdim-1);
            for (int i=1; i<=yHalf; i++)
            {
                int isym=intWRAP(-i,0,Ydim-1);
                std::complex<double> mean=0.5*(
                                              DIRECT_A3D_ELEM(fFourier,k,i,0)+
                                              conj(DIRECT_A3D_ELEM(fFourier,ksym,isym,0)));
                DIRECT_A3D_ELEM(fFourier,k,i,0)=mean;
                DIRECT_A3D_ELEM(fFourier,ksym,isym,0)=conj(mean);
            }
        }
        for (int k=1; k<=zHalf; k++)
        {
            int ksym=intWRAP(-k,0,Zdim-1);
            std::complex<double> mean=0.5*(
                                          DIRECT_A3D_ELEM(fFourier,k,0,0)+
                                          conj(DIRECT_A3D_ELEM(fFourier,ksym,0,0)));
            DIRECT_A3D_ELEM(fFourier,k,0,0)=mean;
            DIRECT_A3D_ELEM(fFourier,ksym,0,0)=conj(mean);
        }
        break;
    }
}

/* FFT Magnitude  ------------------------------------------------------- */
void FFT_magnitude(const MultidimArray< std::complex<double> > &v,
                   MultidimArray<double> &mag)
{
    mag.resizeNoCopy(v);
    double * ptrv=(double *)MULTIDIM_ARRAY(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v)
    {
        double re=*ptrv;
        double im=*(ptrv+1);
        DIRECT_MULTIDIM_ELEM(mag, n) = sqrt(re*re+im*im);
        ptrv+=2;
    }
}

/* FFT Phase ------------------------------------------------------- */
void FFT_phase(const MultidimArray< std::complex<double> > &v,
               MultidimArray<double> &phase)
{
    phase.resizeNoCopy(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v)
    DIRECT_MULTIDIM_ELEM(phase, n) = atan2(DIRECT_MULTIDIM_ELEM(v,n).imag(), DIRECT_MULTIDIM_ELEM(v, n).real());
}


void convolutionFFT(const MultidimArray<double> &img,
                    const MultidimArray<double> &kernel,
                    MultidimArray<double> &result)
{
    if (&result != &img)
        result = img;

    MultidimArray<double> imgTemp;
    MultidimArray< std::complex< double> > FFTIm, FFTK;
    FourierTransformer transformer1(FFTW_BACKWARD), transformer2(FFTW_BACKWARD);

    transformer2.FourierTransform((MultidimArray<double> &)kernel, FFTK, false);

    for (size_t n = 0; n < ZSIZE(result); n++)
    {
        imgTemp.aliasSlice(result, n);
        transformer1.FourierTransform(imgTemp, FFTIm, false);

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FFTIm)
        DIRECT_MULTIDIM_ELEM(FFTIm,n) *= DIRECT_MULTIDIM_ELEM(FFTK,n);

        transformer1.inverseFourierTransform();

        CenterFFT(imgTemp, false);
    }

}


// Fourier ring correlation -----------------------------------------------
//#define SAVE_REAL_PART
void frc_dpr(MultidimArray< double > & m1,
             MultidimArray< double > & m2,
             double sampling_rate,
             MultidimArray< double >& freq,
             MultidimArray< double >& frc,
             MultidimArray< double >& frc_noise,
             MultidimArray< double >& dpr,
             MultidimArray< double >& error_l2,
             bool dodpr)
{
    if (!m1.sameShape(m2))
        REPORT_ERROR(ERR_MULTIDIM_SIZE,"MultidimArrays have different shapes!");

    int m1sizeX = XSIZE(m1), m1sizeY = YSIZE(m1), m1sizeZ = ZSIZE(m1);

    MultidimArray< std::complex< double > > FT1;
    FourierTransformer transformer1(FFTW_BACKWARD);
    transformer1.FourierTransform(m1, FT1, false);
    m1.clear(); // Free memory

    MultidimArray< std::complex< double > > FT2;
    FourierTransformer transformer2(FFTW_BACKWARD);
    transformer2.FourierTransform(m2, FT2, false);
    m2.clear(); // Free memory

    MultidimArray< int > radial_count(m1sizeX/2+1);
    MultidimArray<double> num, den1, den2, den_dpr;
    Matrix1D<double> f(3);
    num.initZeros(radial_count);
    den1.initZeros(radial_count);
    den2.initZeros(radial_count);

    //dpr calculation takes for ever in large volumes
    //since atan2 is called many times
    //untill atan2 is changed by a table let us make dpr an option
    if (dodpr)
    {
        dpr.initZeros(radial_count);
        den_dpr.initZeros(radial_count);
    }
    freq.initZeros(radial_count);
    frc.initZeros(radial_count);
    frc_noise.initZeros(radial_count);
    error_l2.initZeros(radial_count);
#ifdef SAVE_REAL_PART

    std::vector<double> *realPart=
        new std::vector<double>[XSIZE(radial_count)];
#endif

    int sizeZ_2 = m1sizeZ/2;
    if (sizeZ_2==0)
        sizeZ_2=1;
    double isizeZ = 1.0/m1sizeZ;
    int sizeY_2 = m1sizeY/2;
    double iysize = 1.0/m1sizeY;
    int sizeX_2 = m1sizeX/2;
    double ixsize = 1.0/m1sizeX;

    int ZdimFT1=(int)ZSIZE(FT1);
    int YdimFT1=(int)YSIZE(FT1);
    int XdimFT1=(int)XSIZE(FT1);
    for (int k=0; k<ZdimFT1; k++)
    {
        FFT_IDX2DIGFREQ_FAST(k,m1sizeZ,sizeZ_2,isizeZ,ZZ(f));
        double fz2=ZZ(f)*ZZ(f);
        for (int i=0; i<YdimFT1; i++)
        {
            FFT_IDX2DIGFREQ_FAST(i,YSIZE(m1),sizeY_2, iysize, YY(f));
            double fz2_fy2=fz2 + YY(f)*YY(f);
            for (int j=0; j<XdimFT1; j++)
            {
                FFT_IDX2DIGFREQ_FAST(j,m1sizeX, sizeX_2, ixsize, XX(f));

                double R2 =fz2_fy2 + XX(f)*XX(f);

                if (R2>0.25)
                    continue;

                double R = sqrt(R2);
                int idx = (int)round(R * m1sizeX);
                std::complex<double> &z1 = dAkij(FT1, k, i, j);
                std::complex<double> &z2 = dAkij(FT2, k, i, j);
                double absz1 = abs(z1);
                double absz2 = abs(z2);
                dAi(num,idx) += real(conj(z1) * z2);
                dAi(den1,idx) += absz1*absz1;
                dAi(den2,idx) += absz2*absz2;
                dAi(error_l2,idx) += abs(z1-z2);
                if (dodpr) //this takes to long for a huge volume
                {
                    double phaseDiff=atan2(z1.imag(),z1.real()) - atan2(z2.imag(),z2.real());
                    phaseDiff = RAD2DEG(phaseDiff);
                    phaseDiff = realWRAP(phaseDiff,-180, 180);
                    dAi(dpr,idx) += ((absz1+absz2)*phaseDiff*phaseDiff);
                    dAi(den_dpr,idx) += (absz1+absz2);
#ifdef SAVE_REAL_PART

                    realPart[idx].push_back(z1.real());
#endif

                }
                dAi(radial_count,idx)++;
            }
        }
    }

    FOR_ALL_ELEMENTS_IN_ARRAY1D(freq)
    {
        dAi(freq,i) = (double) i / (m1sizeX * sampling_rate);
        dAi(frc,i) = dAi(num,i)/sqrt(dAi(den1,i)*dAi(den2,i));
        dAi(frc_noise,i) = 2 / sqrt((double) dAi(radial_count,i));
        dAi(error_l2,i) /= dAi(radial_count,i);

        if (dodpr)
            dAi(dpr,i) = sqrt(dAi(dpr,i) / dAi(den_dpr,i));
#ifdef SAVE_REAL_PART

        std::ofstream fhOut;
        fhOut.open(((std::string)"PPP_RealPart_"+integerToString(i)+".txt").
                   c_str());
        for (int j=0; j<realPart[i].size(); j++)
            fhOut << realPart[i][j] << std::endl;
        fhOut.close();
#endif

    }
}



void scaleToSizeFourier(int Zdim, int Ydim, int Xdim, MultidimArray<double> &mdaIn, MultidimArray<double> &mdaOut, int nThreads)
{
    //Mmem = *this
    //memory for fourier transform output
    MultidimArray<std::complex<double> > MmemFourier;
    // Perform the Fourier transform
    FourierTransformer transformerM;
    transformerM.setThreadsNumber(nThreads);
    transformerM.FourierTransform(mdaIn, MmemFourier, false);

    // Create space for the downsampled image and its Fourier transform
    mdaOut.resizeNoCopy(Zdim, Ydim, Xdim);
    MultidimArray<std::complex<double> > MpmemFourier;
    FourierTransformer transformerMp;
    transformerMp.setReal(mdaOut);
    transformerMp.getFourierAlias(MpmemFourier);

    size_t ihalf = XMIPP_MIN((YSIZE(MpmemFourier)/2+1),(YSIZE(MmemFourier)/2+1));
    size_t zhalf = XMIPP_MIN((ZSIZE(MpmemFourier)/2+1),(ZSIZE(MmemFourier)/2+1));
    size_t xsize = XMIPP_MIN((XSIZE(MmemFourier)),(XSIZE(MpmemFourier)));
    //Init with zero
    MpmemFourier.initZeros();

    //TODO: Check if we can accelerate the last loop by using memcpy ???
    for (size_t k = 0; k < zhalf; ++k)
    {
        for (size_t i=0; i<ihalf; i++)
            for (size_t j=0; j<xsize; j++)
                dAkij(MpmemFourier,k,i,j) = dAkij(MmemFourier,k,i,j);
        for (size_t i=YSIZE(MpmemFourier)-1; i>=ihalf; i--)
        {
            size_t ip = i + YSIZE(MmemFourier)-YSIZE(MpmemFourier) ;
            for (size_t j=0; j<XSIZE(MpmemFourier); j++)
                dAkij(MpmemFourier,k,i,j) = dAkij(MmemFourier,k,ip,j);
        }
    }
    for (size_t k = ZSIZE(MpmemFourier)-1; k >= zhalf; --k)
    {
        size_t kp = k + ZSIZE(MmemFourier)-ZSIZE(MpmemFourier) ;
        for (size_t i=0; i<ihalf; i++)
            for (size_t j=0; j<xsize; j++)
                dAkij(MpmemFourier,k,i,j) = dAkij(MmemFourier,kp,i,j);
        for (size_t i=YSIZE(MpmemFourier)-1; i>=ihalf; i--)
        {
            size_t ip = i + YSIZE(MmemFourier)-YSIZE(MpmemFourier) ;
            for (size_t j=0; j<XSIZE(MpmemFourier); j++)
                dAkij(MpmemFourier,k,i,j) = dAkij(MmemFourier,kp,ip,j);
        }
    }
    // Transform data
    transformerMp.inverseFourierTransform();
}


void selfScaleToSizeFourier(int Zdim, int Ydim, int Xdim, MultidimArray<double> &mda, int nThreads)
{
  scaleToSizeFourier(Zdim, Ydim, Xdim, mda, mda, nThreads);
}


void selfScaleToSizeFourier(int Ydim, int Xdim, MultidimArray<double>& Mpmem,int nThreads)
{
    selfScaleToSizeFourier(1, Ydim, Xdim, Mpmem, nThreads);
}

void selfScaleToSizeFourier(int Zdim, int Ydim, int Xdim, MultidimArrayGeneric &Mpmem, int nThreads)
{
    MultidimArray<double> aux;
    Mpmem.getImage(aux);
    selfScaleToSizeFourier(Zdim, Ydim, Xdim, aux, nThreads);
    Mpmem.setImage(aux);
}

void selfScaleToSizeFourier(int Ydim, int Xdim, MultidimArrayGeneric &Mpmem, int nThreads)
{
    MultidimArray<double> aux;
    Mpmem.getImage(aux);
    selfScaleToSizeFourier(1, Ydim, Xdim, aux, nThreads);
    Mpmem.setImage(aux);
}

void getSpectrum(MultidimArray<double> &Min,
                 MultidimArray<double> &spectrum,
                 int spectrum_type)
{
    Min.checkDimension(3);

    MultidimArray<std::complex<double> > Faux;
    int xsize = XSIZE(Min);
    Matrix1D<double> f(3);
    MultidimArray<double> count(xsize);
    FourierTransformer transformer;

    spectrum.initZeros(xsize);
    count.initZeros();
    transformer.FourierTransform(Min, Faux, false);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Faux)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(Faux),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(Faux),ZZ(f));
        double R=f.module();
        //if (R>0.5) continue;
        int idx = round(R*xsize);
        double F=abs(dAkij(Faux, k, i, j));
        if (spectrum_type == AMPLITUDE_SPECTRUM)
            A1D_ELEM(spectrum,idx) += F;
        else
            A1D_ELEM(spectrum,idx) += F*F;
        A1D_ELEM(count,idx) += 1.;
    }
    for (int i = 0; i < xsize; i++)
        if (A1D_ELEM(count,i) > 0.)
            A1D_ELEM(spectrum,i) /= A1D_ELEM(count,i);
}

void divideBySpectrum(MultidimArray<double> &Min,
                      MultidimArray<double> &spectrum,
                      bool leave_origin_intact)
{

    Min.checkDimension(3);

    MultidimArray<double> div_spec(spectrum);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(spectrum)
    {
        if (ABS(dAi(spectrum,i)) > 0.)
            dAi(div_spec,i) = 1./dAi(spectrum,i);
        else
            dAi(div_spec,i) = 1.;
    }
    multiplyBySpectrum(Min,div_spec,leave_origin_intact);
}

void multiplyBySpectrum(MultidimArray<double> &Min,
                        MultidimArray<double> &spectrum,
                        bool leave_origin_intact)
{
    Min.checkDimension(3);

    MultidimArray<std::complex<double> > Faux;
    Matrix1D<double> f(3);
    MultidimArray<double> lspectrum;
    FourierTransformer transformer;
    double dim3 = XSIZE(Min)*YSIZE(Min)*ZSIZE(Min);

    transformer.FourierTransform(Min, Faux, false);
    lspectrum=spectrum;
    if (leave_origin_intact)
        lspectrum(0)=1.;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Faux)
    {
        FFT_IDX2DIGFREQ(j,XSIZE(Min), XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(Faux),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(Faux),ZZ(f));
        double R=f.module();
        //if (R > 0.5) continue;
        int idx=ROUND(R*XSIZE(Min));
        dAkij(Faux, k, i, j) *=  lspectrum(idx) * dim3;
    }
    transformer.inverseFourierTransform();
}

void whitenSpectrum(MultidimArray<double> &Min,
                    MultidimArray<double> &Mout,
                    int spectrum_type,
                    bool leave_origin_intact)
{
    Min.checkDimension(3);

    MultidimArray<double> spectrum;
    getSpectrum(Min,spectrum,spectrum_type);
    Mout=Min;
    divideBySpectrum(Mout,spectrum,leave_origin_intact);

}

void adaptSpectrum(MultidimArray<double> &Min,
                   MultidimArray<double> &Mout,
                   const MultidimArray<double> &spectrum_ref,
                   int spectrum_type,
                   bool leave_origin_intact)
{

    Min.checkDimension(3);

    MultidimArray<double> spectrum;
    getSpectrum(Min,spectrum,spectrum_type);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(spectrum)
    {
        dAi(spectrum, i) = (dAi(spectrum, i) > 0.) ? dAi(spectrum_ref,i)/ dAi(spectrum, i) : 1.;
    }
    Mout=Min;
    multiplyBySpectrum(Mout,spectrum,leave_origin_intact);

}

void correlation_matrix(const MultidimArray<double> & m1,
                        const MultidimArray<double> & m2,
                        MultidimArray< double >& R,
                        CorrelationAux &aux,
                        bool center)
{
    // Compute the Fourier Transforms
    aux.transformer1.FourierTransform((MultidimArray<double> &)m1, aux.FFT1, false);
    correlation_matrix(aux.FFT1,m2,R,aux,center);
}


void correlation_matrix(const MultidimArray< std::complex< double > > & FF1,
                        const MultidimArray<double> & m2,
                        MultidimArray<double>& R,
                        CorrelationAux &aux,
                        bool center)
{
    R=m2;
    aux.transformer2.FourierTransform(R, aux.FFT2, false);

    // Multiply FFT1 * FFT2'
    double dSize=MULTIDIM_SIZE(R);
    double mdSize=-dSize;
    double a, b, c, d; // a+bi, c+di
    double *ptrFFT2=(double*)MULTIDIM_ARRAY(aux.FFT2);
    double *ptrFFT1=(double*)MULTIDIM_ARRAY(FF1);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FF1)
    {
        a=*ptrFFT1++;
        b=*ptrFFT1++;
        c=(*ptrFFT2)*dSize;
        d=(*(ptrFFT2+1))*mdSize;
        *ptrFFT2++ = a*c-b*d;
        *ptrFFT2++ = b*c+a*d;
    }

    // Invert the product, in order to obtain the correlation image
    aux.transformer2.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    if (center)
        CenterFFT(R, true);
}


void fast_correlation_vector(const MultidimArray< std::complex<double> > & FFT1,
                        const MultidimArray< std::complex<double> > & FFT2,
                        MultidimArray< double >& R,
                        FourierTransformer &transformer)
{
	transformer.setFourier(FFT1);

    // Multiply FFT1 * FFT2'
    double a, b, c, d; // a+bi, c+di
    double *ptrFFT2=(double*)MULTIDIM_ARRAY(FFT2);
    double *ptrFFT1=(double*)&(A1D_ELEM(transformer.fFourier,0));
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FFT1)
    {
        a=*ptrFFT1;
        b=*(ptrFFT1+1);
        c=(*ptrFFT2++);
        d=(*ptrFFT2++)*(-1);
        *ptrFFT1++ = a*c-b*d;
        *ptrFFT1++ = b*c+a*d;
    }

    // Invert the product, in order to obtain the correlation image
    transformer.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    R=*transformer.fReal;
    CenterFFT(R, true);
    R.setXmippOrigin();
}

/** Fast autocorrelation matrix */
void auto_correlation_matrix(const MultidimArray<double> & Img, MultidimArray< double >& R, CorrelationAux &aux)
{
    // Compute the Fourier Transform
    R=Img;
    aux.transformer1.FourierTransform(R, aux.FFT1, false);

    // Multiply FFT1 * FFT1'
    double dSize=MULTIDIM_SIZE(Img);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(aux.FFT1)
    {
        double *ptr=(double*)&DIRECT_MULTIDIM_ELEM(aux.FFT1,n);
        double &realPart=*ptr;
        double &imagPart=*(ptr+1);
        realPart=dSize*(realPart*realPart+imagPart*imagPart);
        imagPart=0;
    }

    // Invert the product, in order to obtain the correlation image
    aux.transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

void randomizePhases(MultidimArray<double> &Min, double wRandom)
{
	FourierTransformer transformer;
	MultidimArray< std::complex<double> > F;
	transformer.FourierTransform(Min,F,false);

	Matrix1D<double> f(3);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(F)
    {
        FFT_IDX2DIGFREQ(j,XSIZE(Min), XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(F),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(F),ZZ(f));
        double w=f.module();
        if (w > wRandom)
        {
        	double alpha=rnd_unif(0,2*PI);
        	double c,s;
        	sincos(alpha,&s,&c);
        	DIRECT_A3D_ELEM(F,k,i,j)*=std::complex<double>(s,c);
        }
    }
	transformer.inverseFourierTransform();
}
