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

#include "fftw.h"
#include <string.h>
#include <pthread.h>

static pthread_mutex_t fftw_plan_mutex = PTHREAD_MUTEX_INITIALIZER; 

// Constructors and destructors --------------------------------------------
XmippFftw::XmippFftw()
{    
    fPlanForward=NULL;
    fPlanBackward=NULL;
    fReal=NULL;
    fComplex=NULL;
    dataPtr=NULL;
    complexDataPtr=NULL;
    nthreads=1;
    threadsSetOn=false;
}

void XmippFftw::clear()
{
    fReal=NULL;
    fComplex=NULL;
    fFourier.clear();
    // Anything to do with plans has to be protected for threads!
    pthread_mutex_lock(&fftw_plan_mutex);
        if (fPlanForward !=NULL) fftw_destroy_plan(fPlanForward);
        if (fPlanBackward!=NULL) fftw_destroy_plan(fPlanBackward);
    pthread_mutex_unlock(&fftw_plan_mutex);
    fPlanForward     = NULL;
    fPlanBackward    = NULL;
    dataPtr          = NULL;
    complexDataPtr   = NULL;
}

XmippFftw::~XmippFftw()
{
    clear();
}

// Initialization ----------------------------------------------------------
const MultidimArray<double> &XmippFftw::getReal() const
{
    return (*fReal);
}

const MultidimArray<std::complex<double> > &XmippFftw::getComplex() const
{
    return (*fComplex);
}


void XmippFftw::setReal(MultidimArray<double> &input)
{
    bool recomputePlan=false;
    if (fReal==NULL) recomputePlan=true;
    else if (dataPtr!=MULTIDIM_ARRAY(input)) recomputePlan=true;
    else recomputePlan=!(fReal->sameShape(input));
    fFourier.resize(ZSIZE(input),YSIZE(input),XSIZE(input)/2+1);
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
            if (fPlanForward!=NULL)  fftw_destroy_plan(fPlanForward);
            fPlanForward=NULL;
            fPlanForward = fftw_plan_dft_r2c(ndim, N, MULTIDIM_ARRAY(*fReal),
                (fftw_complex*) MULTIDIM_ARRAY(fFourier), FFTW_ESTIMATE);
            if (fPlanBackward!=NULL) fftw_destroy_plan(fPlanBackward);
            fPlanBackward=NULL;
            fPlanBackward = fftw_plan_dft_c2r(ndim, N,
                (fftw_complex*) MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(*fReal),
                FFTW_ESTIMATE);
            if (fPlanForward == NULL || fPlanBackward == NULL)
                REPORT_ERROR(1, "FFTW plans cannot be created");
            delete [] N;
            dataPtr=MULTIDIM_ARRAY(*fReal);
        pthread_mutex_unlock(&fftw_plan_mutex);
    }
}

void XmippFftw::setReal(MultidimArray<std::complex<double> > &input)
{
    bool recomputePlan=false;
    if (fComplex==NULL) recomputePlan=true;
    else if (complexDataPtr!=MULTIDIM_ARRAY(input)) recomputePlan=true;
    else recomputePlan=!(fComplex->sameShape(input));
    fFourier.resize(input);
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
            if (fPlanForward!=NULL)  fftw_destroy_plan(fPlanForward);
            fPlanForward=NULL;
            fPlanForward = fftw_plan_dft(ndim, N, (fftw_complex*) MULTIDIM_ARRAY(*fComplex),
                                         (fftw_complex*) MULTIDIM_ARRAY(fFourier), FFTW_FORWARD, FFTW_ESTIMATE);
            if (fPlanBackward!=NULL) fftw_destroy_plan(fPlanBackward);
            fPlanBackward=NULL;
            fPlanBackward = fftw_plan_dft(ndim, N, (fftw_complex*) MULTIDIM_ARRAY(fFourier), 
                                          (fftw_complex*) MULTIDIM_ARRAY(*fComplex), FFTW_BACKWARD, FFTW_ESTIMATE);
            if (fPlanForward == NULL || fPlanBackward == NULL)
                REPORT_ERROR(1, "FFTW plans cannot be created");
            delete [] N;
            complexDataPtr=MULTIDIM_ARRAY(*fComplex);
        pthread_mutex_unlock(&fftw_plan_mutex);
    }
}

void XmippFftw::setFourier(MultidimArray<std::complex<double> > &inputFourier)
{
    memcpy(MULTIDIM_ARRAY(fFourier),MULTIDIM_ARRAY(inputFourier),
        MULTIDIM_SIZE(inputFourier)*2*sizeof(double));
}

// Transform ---------------------------------------------------------------
void XmippFftw::Transform(int sign)
{
    if (sign == FFTW_FORWARD)
    {
        fftw_execute(fPlanForward);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fFourier)
            DIRECT_MULTIDIM_ELEM(fFourier,n) /= MULTIDIM_SIZE(*fReal);
    }
    else if (sign == FFTW_BACKWARD)
        fftw_execute(fPlanBackward);
}

void XmippFftw::FourierTransform()
{
    Transform(FFTW_FORWARD);
}

void XmippFftw::inverseFourierTransform()
{
    Transform(FFTW_BACKWARD);
}

// Inforce Hermitian symmetry ---------------------------------------------
void XmippFftw::enforceHermitianSymmetry()
{
    int ndim=3;
    if (ZSIZE(*fReal)==1)
    {
        ndim=2;
        if (YSIZE(*fReal)==1)
            ndim=1;
    }
    int yHalf=YSIZE(*fReal)/2;
    if (YSIZE(*fReal)%2==0) yHalf--;
    int zHalf=ZSIZE(*fReal)/2;
    if (ZSIZE(*fReal)%2==0) zHalf--;
    switch (ndim)
    {
        case 2:
            for (int i=1; i<=yHalf; i++)
            {
                int isym=intWRAP(-i,0,YSIZE(*fReal)-1);
                std::complex<double> mean=0.5*(
                    DIRECT_MAT_ELEM(fFourier,i,0)+
                    conj(DIRECT_MAT_ELEM(fFourier,isym,0)));
                DIRECT_MAT_ELEM(fFourier,i,0)=mean;
                DIRECT_MAT_ELEM(fFourier,isym,0)=conj(mean);
            }
            break;
        case 3:
            for (int k=0; k<ZSIZE(*fReal); k++)
            {
                int ksym=intWRAP(-k,0,ZSIZE(*fReal)-1);
                for (int i=1; i<=yHalf; i++)
                {
                    int isym=intWRAP(-i,0,YSIZE(*fReal)-1);
                    std::complex<double> mean=0.5*(
                        DIRECT_VOL_ELEM(fFourier,k,i,0)+
                        conj(DIRECT_VOL_ELEM(fFourier,ksym,isym,0)));
                    DIRECT_VOL_ELEM(fFourier,k,i,0)=mean;
                    DIRECT_VOL_ELEM(fFourier,ksym,isym,0)=conj(mean);
                }
            }
            for (int k=1; k<=zHalf; k++)
            {
                int ksym=intWRAP(-k,0,ZSIZE(*fReal)-1);
                std::complex<double> mean=0.5*(
                    DIRECT_VOL_ELEM(fFourier,k,0,0)+
                    conj(DIRECT_VOL_ELEM(fFourier,ksym,0,0)));
                DIRECT_VOL_ELEM(fFourier,k,0,0)=mean;
                DIRECT_VOL_ELEM(fFourier,ksym,0,0)=conj(mean);
            }
            break;
    }
}

/* FFT Magnitude 1D. ------------------------------------------------------- */
void FFT_magnitude(const Matrix1D< std::complex<double> > &v,
                   Matrix1D<double> &mag)
{
    mag.resize(v);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(v) mag(i) = abs(v(i));
}

void FFT_magnitude(const Matrix2D< std::complex<double> > &v,
                   Matrix2D<double> &mag)
{
    mag.resize(v);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(v) mag(i, j) = abs(v(i, j));
}

void FFT_magnitude(const Matrix3D< std::complex<double> > &v,
                   Matrix3D<double> &mag)
{
    mag.resize(v);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(v) mag(k, i, j) = abs(v(k, i, j));
}

/* FFT Phase 1D. ------------------------------------------------------- */
void FFT_phase(const Matrix1D< std::complex<double> > &v,
               Matrix1D<double> &phase)
{
    phase.resize(v);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(v) phase(i) = atan2(v(i).imag(), v(i).real());
}

void FFT_phase(const Matrix2D< std::complex<double> > &v,
               Matrix2D<double> &phase)
{
    phase.resize(v);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(v) phase(i, j) = atan2(v(i, j).imag(), v(i, j).real());
}

void FFT_phase(const Matrix3D< std::complex<double> > &v,
               Matrix3D<double> &phase)
{
    phase.resize(v);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(v) phase(k, i, j) = atan2(v(k, i, j).imag(),
            v(k, i, j).real());
}

// Fourier ring correlation -----------------------------------------------
void frc_dpr(Matrix2D< double > & m1,
             Matrix2D< double > & m2,
             double sampling_rate,
             Matrix1D< double >& freq,
             Matrix1D< double >& frc,
             Matrix1D< double >& frc_noise,
             Matrix1D< double >& dpr,
             bool skipdpr)
{
    if (!m1.sameShape(m2))
        REPORT_ERROR(1,"Matrices have different shapes!");

    Matrix2D< std::complex< double > > FT1;
    XmippFftw transformer1;
    transformer1.FourierTransform(m1, FT1, false);

    Matrix2D< std::complex< double > > FT2;
    XmippFftw transformer2;
    transformer2.FourierTransform(m2, FT2, false);
/*
std::cout << STARTINGY(FT1) << "," << FINISHINGY(FT1) << std::endl;
std::cout << STARTINGX(FT1) << "," << FINISHINGX(FT1) << std::endl;

std::cout << STARTINGY(m1) << "," << FINISHINGY(m1) << std::endl;
std::cout << STARTINGX(m1) << "," << FINISHINGX(m1) << std::endl;

    for (int i=STARTINGY(FT1); i<=FINISHINGY(FT1); i++)
    {
        for (int j=STARTINGX(FT1); j<=FINISHINGX(FT1); j++)
	{
		FT1(i,j)=std::complex< double >( m1(i-128,j-128), 0.0 );
	}
    }

    for (int i=STARTINGY(FT2); i<=FINISHINGY(FT2); i++)
    {
        for (int j=STARTINGX(FT2); j<=FINISHINGX(FT2); j++)
	{
		FT2(i,j)=std::complex< double >( m2(i-128,j-128), 0.0 );
	}
    }*/

    Matrix1D< int > radial_count(XSIZE(m1)/2+1);
    Matrix1D<double> num, den1, den2;
    Matrix1D<double> f(3);
    num.initZeros(radial_count);
    den1.initZeros(radial_count);
    den2.initZeros(radial_count);
     
    //dpr calculation takes for ever in large volumes
    //since atan2 is called many times
    //untill atan2 is changed by a table let us make dpr an option
    if (skipdpr)
        dpr.initZeros(radial_count);
    freq.initZeros(radial_count);
    frc.initZeros(radial_count);
    frc_noise.initZeros(radial_count);

    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(FT1)
    {
        FFT_IDX2DIGFREQ(j,XSIZE(m1),XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(m1),YY(f));
        double R=f.module();
        if (R>0.5) continue;
        int idx=ROUND(R*XSIZE(m1));
        std::complex<double> z1=dMij(FT1, i, j);
        std::complex<double> z2=dMij(FT2, i, j);
        double absz1=abs(z1);
        double absz2=abs(z2);
        num(idx)+=real(conj(z1) * z2);
        den1(idx)+= absz1*absz1;
        den2(idx)+= absz2*absz2;
        if (skipdpr)
        {
            double phaseDiff=realWRAP(RAD2DEG((atan2(z1.imag(), z1.real())) -
                    (atan2(z2.imag(), z2.real()))),-180, 180);
            dpr(idx)+=sqrt((absz1+absz2)*phaseDiff*phaseDiff/(absz1+absz2));
        }
        radial_count(idx)++;
    }

    FOR_ALL_ELEMENTS_IN_MATRIX1D(freq)
    {
        freq(i) = (double) i / (XSIZE(m1) * sampling_rate);
        frc(i) = num(i)/sqrt(den1(i)*den2(i));
        frc_noise(i) = 2 / sqrt((double) radial_count(i));
        if (skipdpr)
            dpr(i)/=radial_count(i);
    }
}

void frc_dpr(Matrix3D< double > & m1,
             Matrix3D< double > & m2, double sampling_rate,
             Matrix1D< double >& freq,
             Matrix1D< double >& frc,
             Matrix1D< double >& frc_noise,
             Matrix1D< double >& dpr,
             bool skipdpr)
{
    if (!m1.sameShape(m2))
        REPORT_ERROR(1,"Volumes have different shapes!");

    Matrix3D< std::complex< double > > FT1;
    XmippFftw transformer1;
    transformer1.FourierTransform(m1, FT1, false);

    Matrix3D< std::complex< double > > FT2;
    XmippFftw transformer2;
    transformer2.FourierTransform(m2, FT2, false);

    Matrix1D< int > radial_count(XSIZE(m1)/2+1);
    Matrix1D<double> num, den1, den2, f(3);
    num.initZeros(radial_count);
    den1.initZeros(radial_count);
    den2.initZeros(radial_count);

    if (skipdpr)
        dpr.initZeros(radial_count);
    freq.initZeros(radial_count);
    frc.initZeros(radial_count);
    frc_noise.initZeros(radial_count);

    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FT1)
    {
        FFT_IDX2DIGFREQ(j,XSIZE(m1),XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(m1),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(m1),ZZ(f));
        double R=f.module();
        if (R>0.5) continue;
        int idx=ROUND(R*XSIZE(m1));
        std::complex<double> z1=dVkij(FT1, k, i, j);
        std::complex<double> z2=dVkij(FT2, k, i, j);
        double absz1=abs(z1);
        double absz2=abs(z2);
        num(idx)+=real(conj(z1) * z2);
        den1(idx)+= absz1*absz1;
        den2(idx)+= absz2*absz2;
        if (skipdpr)
        {    
            double phaseDiff=realWRAP(RAD2DEG((atan2(z1.imag(), z1.real())) -
                    (atan2(z2.imag(), z2.real()))),-180, 180);
            dpr(idx)+=sqrt((absz1+absz2)*phaseDiff*phaseDiff/(absz1+absz2));
        }
        radial_count(idx)++;
    }

    frc.initZeros(radial_count);
    frc_noise.resize(radial_count);
    freq.resize(radial_count);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(freq)
    {
        freq(i) = (double) i / (XSIZE(m1) * sampling_rate);
        if(num(i)!=0)
           frc(i) = num(i)/sqrt(den1(i)*den2(i));
        else
           frc(i) = 0;
        frc_noise(i) = 2 / sqrt((double) radial_count(i));
        if (skipdpr)
            dpr(i)/=radial_count(i);
    }
}
void selfScaleToSizeFourier(int Ydim, int Xdim,Matrix2D<double>& Mpmem,int nThreads) 
 {
     //Mmem = *this
     //memory for fourier transform output
     Matrix2D<std::complex<double> > MmemFourier;
     // Perform the Fourier transform
     XmippFftw transformerM;
     transformerM.setThreadsNumber(nThreads);
     transformerM.FourierTransform(Mpmem, MmemFourier, false);

     // Create space for the downsampled image and its Fourier transform
     Mpmem.resize(Ydim, Xdim);
     Matrix2D<std::complex<double> > MpmemFourier;
     XmippFftw transformerMp;
     transformerMp.setReal(Mpmem);
     transformerMp.getFourierAlias(MpmemFourier);

     int ihalf = XMIPP_MIN((YSIZE(MpmemFourier)/2+1),(YSIZE(MmemFourier)/2+1));
     int xsize = XMIPP_MIN((XSIZE(MmemFourier)),(XSIZE(MpmemFourier)));
     int ysize = XMIPP_MIN((YSIZE(MmemFourier)),(YSIZE(MpmemFourier)));
     //Init with zero
     MpmemFourier.initZeros();
     for (int i=0; i<ihalf; i++)
         for (int j=0; j<xsize; j++)
             MpmemFourier(i,j)=MmemFourier(i,j);
     for (int i=YSIZE(MpmemFourier)-1; i>=ihalf; i--)
     {   
         int ip = i + YSIZE(MmemFourier)-YSIZE(MpmemFourier) ;
         for (int j=0; j<XSIZE(MpmemFourier); j++)
             MpmemFourier(i,j)=MmemFourier(ip,j);
     }

     // Transform data
     transformerMp.inverseFourierTransform();
 }


void getSpectrum(Matrix3D<double> &Min, 
                 Matrix1D<double> &spectrum,
                 int spectrum_type)
{
    Matrix3D<std::complex<double> > Faux;
    int xsize = XSIZE(Min);
    Matrix1D<double> f(3), count(xsize);
    XmippFftw transformer;

    spectrum.initZeros(xsize);
    count.initZeros();
    transformer.FourierTransform(Min, Faux, false);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(Faux),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(Faux),ZZ(f));
        double R=f.module();
        //if (R>0.5) continue;
        int idx=ROUND(R*xsize);
        if (spectrum_type == AMPLITUDE_SPECTRUM)
            spectrum(idx) += abs(dVkij(Faux, k, i, j));
        else
            spectrum(idx) += abs(dVkij(Faux, k, i, j)) * abs(dVkij(Faux, k, i, j));
        count(idx) += 1.;
    }
    for (int i = 0; i < xsize; i++)
        if (count(i) > 0.)
            spectrum(i) /= count(i);
}

void divideBySpectrum(Matrix3D<double> &Min, 
                      Matrix1D<double> &spectrum,
                      bool leave_origin_intact)
{
    
    Matrix1D<double> div_spec(spectrum);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(spectrum)
    {
        if (ABS(dVi(spectrum,i)) > 0.)
            dVi(div_spec,i) = 1./dVi(spectrum,i);
        else
            dVi(div_spec,i) = 1.;
    }
    multiplyBySpectrum(Min,div_spec,leave_origin_intact);
}

void multiplyBySpectrum(Matrix3D<double> &Min, 
                        Matrix1D<double> &spectrum,
                        bool leave_origin_intact)
{
    Matrix3D<std::complex<double> > Faux;
    Matrix1D<double> f(3), lspectrum;
    XmippFftw transformer;
    double dim3 = XSIZE(Min)*YSIZE(Min)*ZSIZE(Min);

    transformer.FourierTransform(Min, Faux, false);
    lspectrum=spectrum;
    if (leave_origin_intact)
        lspectrum(0)=1.;
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
    {
        FFT_IDX2DIGFREQ(j,XSIZE(Min), XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(Faux),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(Faux),ZZ(f));
        double R=f.module();
        //if (R > 0.5) continue;
        int idx=ROUND(R*XSIZE(Min));
        dVkij(Faux, k, i, j) *=  lspectrum(idx) * dim3;
    }
    transformer.inverseFourierTransform();

}


void whitenSpectrum(Matrix3D<double> &Min, 
                    Matrix3D<double> &Mout, 
                    int spectrum_type,
                    bool leave_origin_intact)
{

    Matrix1D<double> spectrum;
    getSpectrum(Min,spectrum,spectrum_type);
    Mout=Min;
    divideBySpectrum(Mout,spectrum,leave_origin_intact);

}

void adaptSpectrum(Matrix3D<double> &Min, 
                   Matrix3D<double> &Mout,
                   const Matrix1D<double> spectrum_ref,
                   int spectrum_type,
                   bool leave_origin_intact)
{
    
    Matrix1D<double> spectrum;
    getSpectrum(Min,spectrum,spectrum_type);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(spectrum)
    {
        if (spectrum(i) > 0.)
            spectrum(i) = spectrum_ref(i)/spectrum(i);
        else
            spectrum(i) = 1.;
    }
    Mout=Min;
    multiplyBySpectrum(Mout,spectrum,leave_origin_intact);

}
