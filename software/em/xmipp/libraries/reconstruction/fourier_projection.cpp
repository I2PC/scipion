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

#include "fourier_projection.h"
#include <data/xmipp_fft.h>

FourierProjector::FourierProjector(MultidimArray<double> &V, double paddFactor, double maxFreq, int degree)
{
    volume = &V;
    volumeSize=XSIZE(*volume);
    paddingFactor = paddFactor;
    maxFrequency = maxFreq;
    BSplineDeg = degree;
    produceSideInfo();
}

void FourierProjector::project(double rot, double tilt, double psi)
{
    double freqy, freqx;
    std::complex< double > f;
    Euler_angles2matrix(rot,tilt,psi,E);

    projectionFourier.initZeros();
    double shift=-FIRST_XMIPP_INDEX(volumeSize);
    double xxshift = -2 * PI * shift / volumeSize;
    double maxFreq2=maxFrequency*maxFrequency;
    double volumePaddedSize=XSIZE(VfourierRealCoefs);
    for (size_t i=0; i<YSIZE(projectionFourier); ++i)
    {
        FFT_IDX2DIGFREQ(i,volumeSize,freqy);
        double freqy2=freqy*freqy;
        double phasey=(double)(i) * xxshift;

        double freqYvol_X=MAT_ELEM(E,1,0)*freqy;
        double freqYvol_Y=MAT_ELEM(E,1,1)*freqy;
        double freqYvol_Z=MAT_ELEM(E,1,2)*freqy;
        for (size_t j=0; j<XSIZE(projectionFourier); ++j)
        {
            // The frequency of pairs (i,j) in 2D
            FFT_IDX2DIGFREQ(j,volumeSize,freqx);

            // Do not consider pixels with high frequency
            if ((freqy2+freqx*freqx)>maxFreq2)
                continue;

            // Compute corresponding frequency in the volume
            double freqvol_X=freqYvol_X+MAT_ELEM(E,0,0)*freqx;
            double freqvol_Y=freqYvol_Y+MAT_ELEM(E,0,1)*freqx;
            double freqvol_Z=freqYvol_Z+MAT_ELEM(E,0,2)*freqx;

            double c,d;
            if (BSplineDeg==0)
            {
                // 0 order interpolation
                // Compute corresponding index in the volume
                int kVolume=(int)round(freqvol_Z*volumePaddedSize);
                int iVolume=(int)round(freqvol_Y*volumePaddedSize);
                int jVolume=(int)round(freqvol_X*volumePaddedSize);
                c = A3D_ELEM(VfourierRealCoefs,kVolume,iVolume,jVolume);
                d = A3D_ELEM(VfourierImagCoefs,kVolume,iVolume,jVolume);
            }
            else if (BSplineDeg==1)
            {
                // B-spline linear interpolation
                double kVolume=freqvol_Z*volumePaddedSize;
                double iVolume=freqvol_Y*volumePaddedSize;
                double jVolume=freqvol_X*volumePaddedSize;
                c=VfourierRealCoefs.interpolatedElement3D(jVolume,iVolume,kVolume);
                d=VfourierImagCoefs.interpolatedElement3D(jVolume,iVolume,kVolume);
            }
            else
            {
                // B-spline cubic interpolation
                double kVolume=freqvol_Z*volumePaddedSize;
                double iVolume=freqvol_Y*volumePaddedSize;
                double jVolume=freqvol_X*volumePaddedSize;
                c=VfourierRealCoefs.interpolatedElementBSpline3D(jVolume,iVolume,kVolume);
                d=VfourierImagCoefs.interpolatedElementBSpline3D(jVolume,iVolume,kVolume);
            }

            // Phase shift to move the origin of the image to the corner
            double dotp = (double)(j) * xxshift + phasey;
            double a,b;
            sincos(dotp,&b,&a);

            // Multiply Fourier coefficient in volume times phase shift
            double ac = a * c;
            double bd = b * d;
            double ab_cd = (a + b) * (c + d);

            // And store the multiplication
            double *ptrI_ij=(double *)&DIRECT_A2D_ELEM(projectionFourier,i,j);
            *ptrI_ij = ac - bd;
            *(ptrI_ij+1) = ab_cd - ac - bd;
        }
    }
    //VfourierRealCoefs.clear();
    //VfourierImagCoefs.clear();
    transformer2D.inverseFourierTransform();
}

void FourierProjector::produceSideInfo()
{
    // Zero padding
    MultidimArray<double> Vpadded;
    int paddedDim=(int)(paddingFactor*volumeSize);
    volume->window(Vpadded,FIRST_XMIPP_INDEX(paddedDim),FIRST_XMIPP_INDEX(paddedDim),FIRST_XMIPP_INDEX(paddedDim),
                   LAST_XMIPP_INDEX(paddedDim),LAST_XMIPP_INDEX(paddedDim),LAST_XMIPP_INDEX(paddedDim));
    volume->clear();
    // Make Fourier transform, shift the volume origin to the volume center and center it
    MultidimArray< std::complex<double> > Vfourier;
    transformer3D.completeFourierTransform(Vpadded,Vfourier);
    ShiftFFT(Vfourier, FIRST_XMIPP_INDEX(XSIZE(Vpadded)), FIRST_XMIPP_INDEX(YSIZE(Vpadded)), FIRST_XMIPP_INDEX(ZSIZE(Vpadded)));
    CenterFFT(Vfourier,true);
    Vfourier.setXmippOrigin();

    // Compensate for the Fourier normalization factor
    double K=(double)(XSIZE(Vpadded)*XSIZE(Vpadded)*XSIZE(Vpadded))/(double)(volumeSize*volumeSize);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vfourier)
    DIRECT_MULTIDIM_ELEM(Vfourier,n)*=K;
    Vpadded.clear();
    // Compute Bspline coefficients
    if (BSplineDeg==3)
    {
        MultidimArray< double > VfourierRealAux, VfourierImagAux;
        Complex2RealImag(Vfourier, VfourierRealAux, VfourierImagAux);
        Vfourier.clear();
        produceSplineCoefficients(BSPLINE3,VfourierRealCoefs,VfourierRealAux);
        produceSplineCoefficients(BSPLINE3,VfourierImagCoefs,VfourierImagAux);
        //VfourierRealAux.clear();
        //VfourierImagAux.clear();
    }
    else
        Complex2RealImag(Vfourier, VfourierRealCoefs, VfourierImagCoefs);

    // Allocate memory for the 2D Fourier transform
    projection().initZeros(volumeSize,volumeSize);
    projection().setXmippOrigin();
    transformer2D.FourierTransform(projection(),projectionFourier,false);
}

void projectVolume(FourierProjector &projection, Projection &P, int Ydim, int Xdim,
                   double rot, double tilt, double psi)
{
    projection.project(rot,tilt,psi);
    P() = projection.projection();
}
