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

#include "fftw.h"
#include "fft.h"
#include <string.h>
#include <pthread.h>

pthread_mutex_t fftw_plan_mutex = PTHREAD_MUTEX_INITIALIZER; 

// Constructors and destructors --------------------------------------------
XmippFftw::XmippFftw()
{
    fPlanForward=NULL;
    fPlanBackward=NULL;
    fReal=NULL;
    dataPtr=NULL;
}

void XmippFftw::clear()
{
    fReal=NULL;
    fFourier.clear();
    if (fPlanForward !=NULL) fftw_destroy_plan(fPlanForward);
    if (fPlanBackward!=NULL) fftw_destroy_plan(fPlanBackward);
    fPlanForward     = NULL;
    fPlanBackward    = NULL;
    dataPtr          = NULL;
}

XmippFftw::~XmippFftw()
{
    clear();
}

// Initialization ----------------------------------------------------------
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
            if (fPlanForward == NULL || fPlanForward == NULL)
                REPORT_ERROR(1, "FFTW plans cannot be created");
            delete [] N;
            dataPtr=MULTIDIM_ARRAY(*fReal);
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
        frc(i) = num(i)/sqrt(den1(i)*den2(i));
        frc_noise(i) = 2 / sqrt((double) radial_count(i));
        if (skipdpr)
            dpr(i)/=radial_count(i);
    }
}
