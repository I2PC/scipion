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

void FourierProjector::project(double rot, double tilt, double psi, const MultidimArray<double> *ctf)
{
    double freqy, freqx;
    std::complex< double > f;
    Euler_angles2matrix(rot,tilt,psi,E);

    projectionFourier.initZeros();
    double maxFreq2=maxFrequency*maxFrequency;
    int Xdim=(int)XSIZE(VfourierRealCoefs);
    int Ydim=(int)YSIZE(VfourierRealCoefs);
    int Zdim=(int)ZSIZE(VfourierRealCoefs);

    for (size_t i=0; i<YSIZE(projectionFourier); ++i)
    {
        FFT_IDX2DIGFREQ(i,volumeSize,freqy);
        double freqy2=freqy*freqy;

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

                // Commented for speed-up, the corresponding code is below
                // c=VfourierRealCoefs.interpolatedElementBSpline3D(jVolume,iVolume,kVolume);
                // d=VfourierImagCoefs.interpolatedElementBSpline3D(jVolume,iVolume,kVolume);

                // The code below is a replicate for speed reasons of interpolatedElementBSpline3D
                double z=kVolume;
                double y=iVolume;
                double x=jVolume;

                // Logical to physical
                z -= STARTINGZ(VfourierRealCoefs);
                y -= STARTINGY(VfourierRealCoefs);
                x -= STARTINGX(VfourierRealCoefs);

                int l1 = (int)ceil(x - 2);
                int l2 = l1 + 3;

                int m1 = (int)ceil(y - 2);
                int m2 = m1 + 3;

                int n1 = (int)ceil(z - 2);
                int n2 = n1 + 3;

                c = d = 0.0;
                double aux;
                for (int nn = n1; nn <= n2; nn++)
                {
                    int equivalent_nn=nn;
                    if      (nn<0)
                        equivalent_nn=-nn-1;
                    else if (nn>=Zdim)
                        equivalent_nn=2*Zdim-nn-1;
                    double yxsumRe = 0.0, yxsumIm = 0.0;
                    for (int m = m1; m <= m2; m++)
                    {
                        int equivalent_m=m;
                        if      (m<0)
                            equivalent_m=-m-1;
                        else if (m>=Ydim)
                            equivalent_m=2*Ydim-m-1;
                        double xsumRe = 0.0, xsumIm = 0.0;
                        for (int l = l1; l <= l2; l++)
                        {
                            double xminusl = x - (double) l;
                            int equivalent_l=l;
                            if      (l<0)
                                equivalent_l=-l-1;
                            else if (l>=Xdim)
                                equivalent_l=2*Xdim-l-1;
                            double CoeffRe = (double) DIRECT_A3D_ELEM(VfourierRealCoefs,equivalent_nn,equivalent_m,equivalent_l);
                            double CoeffIm = (double) DIRECT_A3D_ELEM(VfourierImagCoefs,equivalent_nn,equivalent_m,equivalent_l);
                            BSPLINE03(aux,xminusl);
                            xsumRe += CoeffRe * aux;
                            xsumIm += CoeffIm * aux;
                        }

                        double yminusm = y - (double) m;
                        BSPLINE03(aux,yminusm);
						yxsumRe += xsumRe * aux;
						yxsumIm += xsumIm * aux;
                    }

                    double zminusn = z - (double) nn;
                    BSPLINE03(aux,zminusn);
					c += yxsumRe * aux;
					d += yxsumIm * aux;
                }
            }

            // Phase shift to move the origin of the image to the corner
            double a=DIRECT_A2D_ELEM(phaseShiftImgA,i,j);
            double b=DIRECT_A2D_ELEM(phaseShiftImgB,i,j);
            if (ctf!=NULL)
            {
            	double ctfij=DIRECT_A2D_ELEM(*ctf,i,j);
            	a*=ctfij;
            	b*=ctfij;
            }

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
    transformer2D.inverseFourierTransform();
}

void FourierProjector::produceSideInfo()
{
    // Zero padding
    MultidimArray<double> Vpadded;
    int paddedDim=(int)(paddingFactor*volumeSize);
    // JMRT: TODO: I think it is a very poor design to modify the volume passed
    // in the construct, it will be padded anyway, so new memory should be allocated
    volume->window(Vpadded,FIRST_XMIPP_INDEX(paddedDim),FIRST_XMIPP_INDEX(paddedDim),FIRST_XMIPP_INDEX(paddedDim),
                   LAST_XMIPP_INDEX(paddedDim),LAST_XMIPP_INDEX(paddedDim),LAST_XMIPP_INDEX(paddedDim));
    volume->clear();
    // Make Fourier transform, shift the volume origin to the volume center and center it
    MultidimArray< std::complex<double> > Vfourier;
    FourierTransformer transformer3D;
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

        // Release memory as soon as you can
        VfourierRealAux.clear();

        // Remove all those coefficients we are sure we will not use during the projections
        volumePaddedSize=XSIZE(VfourierRealCoefs);
        int idxMax=maxFrequency*XSIZE(VfourierRealCoefs)+10; // +10 is a safety guard
        idxMax=std::min(FINISHINGX(VfourierRealCoefs),idxMax);
        int idxMin=std::max(-idxMax,STARTINGX(VfourierRealCoefs));
        VfourierRealCoefs.selfWindow(idxMin,idxMin,idxMin,idxMax,idxMax,idxMax);

        produceSplineCoefficients(BSPLINE3,VfourierImagCoefs,VfourierImagAux);
        VfourierImagAux.clear();
        VfourierImagCoefs.selfWindow(idxMin,idxMin,idxMin,idxMax,idxMax,idxMax);
    }
    else
        Complex2RealImag(Vfourier, VfourierRealCoefs, VfourierImagCoefs);

    // Allocate memory for the 2D Fourier transform
    projection().initZeros(volumeSize,volumeSize);
    projection().setXmippOrigin();
    transformer2D.FourierTransform(projection(),projectionFourier,false);

    // Calculate phase shift terms
    phaseShiftImgA.initZeros(projectionFourier);
    phaseShiftImgB.initZeros(projectionFourier);
    double shift=-FIRST_XMIPP_INDEX(volumeSize);
    double xxshift = -2 * PI * shift / volumeSize;
    for (size_t i=0; i<YSIZE(projectionFourier); ++i)
    {
        double phasey=(double)(i) * xxshift;
        for (size_t j=0; j<XSIZE(projectionFourier); ++j)
        {
            // Phase shift to move the origin of the image to the corner
            double dotp = (double)(j) * xxshift + phasey;
            sincos(dotp,&DIRECT_A2D_ELEM(phaseShiftImgB,i,j),&DIRECT_A2D_ELEM(phaseShiftImgA,i,j));
        }
    }
}

void projectVolume(FourierProjector &projector, Projection &P, int Ydim, int Xdim,
                   double rot, double tilt, double psi, const MultidimArray<double> *ctf)
{
	projector.project(rot,tilt,psi,ctf);
    P() = projector.projection();
}

