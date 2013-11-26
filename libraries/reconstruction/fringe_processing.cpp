/***************************************************************************
 * Authors:     Javier Vargas (jvargas@cnb.csic.es)
 *
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

#include "fringe_processing.h"
#include <data/xmipp_polynomials.h>
#include <data/xmipp_image.h>
#include <data/multidim_array.h>
#include <data/xmipp_funcs.h>
#include <data/xmipp_fftw.h>
#include "fourier_filter.h"

#include <queue>

//This member function simulates different types of fringe patterns.
// im is the resultant multidimarray
// type is the type of fringe pattern that we want to generate. Possible values are:
// SIMPLY_OPEN_FRINGES, SIMPLY_CLOSED_FRINGES, COMPLEX_OPEN_FRINGES, COMPLEX_CLOSED_FRINGES
// xdim is the image dimension in x axis
// ydim is the image dimension in y axis
void simulPattern(MultidimArray<double> & im, enum FP_TYPE type, int xdim, int ydim, double noiseLevel, const double freq, Matrix1D<int> coefs)
{

    //FIRST WE TRAVEL ALL ELEMENTS fn im
    int ndim = 1;
    int zdim = 1;
    im.resize(ndim,zdim,ydim,xdim,false);
    im.setXmippOrigin();

    double iMaxDim2 = 2./std::max(xdim,ydim);
    PolyZernikes polynom;

    if ( (type == COMPLEX_OPEN_FRINGES) || (type == COMPLEX_CLOSED_FRINGES) )
    {
        MultidimArray< bool > ROI;
        ROI.resizeNoCopy(im);

        ROI.setXmippOrigin();
        FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
        A2D_ELEM(ROI,i,j) = true;

        polynom.zernikePols(coefs,im,ROI);
        im.setXmippOrigin();
    }

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        switch (type)
        {

        case SIMPLY_OPEN_FRINGES :
            {
                A2D_ELEM(im,i,j) = std::cos(j*iMaxDim2*freq)+rnd_gaus(0,noiseLevel);
                break;
            }
        case SIMPLY_CLOSED_FRINGES :
            {
                A2D_ELEM(im,i,j) = std::cos(50*std::exp(-0.5*(std::pow(i*iMaxDim2*freq,2)+std::pow(j*iMaxDim2*freq,2))))+rnd_gaus(0,noiseLevel);
                break;
            }

        case COMPLEX_OPEN_FRINGES :
            {
                A2D_ELEM(im,i,j) = std::cos(j*iMaxDim2*freq+A2D_ELEM(im,i,j))+rnd_gaus(0,noiseLevel);
                break;
            }
        case COMPLEX_CLOSED_FRINGES :
            {
                A2D_ELEM(im,i,j) = std::cos(50*std::exp(-0.5*(std::pow(i*iMaxDim2,2)+std::pow(j*iMaxDim2,2)))+A2D_ELEM(im,i,j))+rnd_gaus(0,noiseLevel);
                break;
            }

        case SIMPLY_CLOSED_FRINGES_MOD :
            {
                A2D_ELEM(im,i,j) =  (std::exp(-((i*i)+(j*j))/(2e3*freq*freq)))*(std::cos(50*std::exp(-0.5*(std::pow(i*iMaxDim2*freq,2)+std::pow(j*iMaxDim2*freq,2))))+rnd_gaus(0,noiseLevel));
                break;
            }
        }
    }

    STARTINGX(im)=STARTINGY(im)=0;
}

//Function to simulate some test fringe patterns.
//sd=SPHT(c) computes the quadrature term of c still affected by the
//direction phase factor. Therefore for a real c=b*cos(phi)
//sd=SPHT(c)=i*exp(i*dir)*b*sin(phi)
//Ref: Kieran G. Larkin, Donald J. Bone, and Michael A. Oldfield, "Natural
//demodulation of two-dimensional fringe patterns. I. General background of the spiral phase quadrature transform," J. Opt. Soc. Am. A 18, 1862-1870 (2001)
void SPTH(FourierTransformer &ftrans, MultidimArray<double> & im, MultidimArray< std::complex <double> > & imProc)
{
    im.setXmippOrigin();
    MultidimArray<std::complex<double> > H, fftIm, imComplex;
    typeCast(im, imComplex);

    // Fourier Transformer
    ftrans.FourierTransform(imComplex, fftIm, false);

    H.resizeNoCopy(fftIm);
    H.setXmippOrigin();

    //std::complex<double> compTemp(0, 0);
    //i -> for exterior filas o Y
    //j -> for interior columns o X
    for (int i=STARTINGY(H); i<=FINISHINGY(H); i++)
    {
    	double i2=((double)i)*i;
    	for (int j=STARTINGX(H); j<=FINISHINGX(H); j++)
    	{
			double j2=((double)j)*j;
			double iR;
			if ( (i!=0) || (j!=0))
				iR=1.0/std::sqrt(i2+j2);
			else
				iR=0.;
			double *ptr=(double*)&A2D_ELEM(H,i,j);
			*ptr=j*iR;
			*(ptr+1)=i*iR;
    	}
    }

    fftIm *= H;
    ftrans.inverseFourierTransform();
    //Here in the Matlab code there is a complex conjugate s
    imProc = imComplex;
}

void orMinDer(const MultidimArray<double> & im, MultidimArray<double > & orMap, MultidimArray<double > & orModMap, int wSize, MultidimArray<bool > & ROI)
{
    size_t NR, NC,NZ, NDim;
    im.getDimensions(NC,NR,NZ,NDim);

    if ( (NZ!=1) || (NDim!=1) )
        REPORT_ERROR(ERR_MULTIDIM_DIM,(std::string)"ZDim and NDim has to be equals to one");

    MultidimArray<double > d0,d45,d90,d135,mask,orn,ornMod;
    d0.resizeNoCopy(im);
    d45.resizeNoCopy(im);
    d90.resizeNoCopy(im);
    d135.resizeNoCopy(im);
    //mask is the region to perform average in d0, d90, d90, d135 maps
    //ROI is the region of interested. If ROI is 1 then process in this px if ROI is 0 not process
    mask.resizeNoCopy(im);
    orn.resizeNoCopy(im);
    ornMod.resizeNoCopy(im);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(im)
    {
        A2D_ELEM(mask,i,j)  = 0;

        if ( (i==0) || (i== (NR-1)) || (j == 0) || (j== (NC-1)) )
        {
            DIRECT_A2D_ELEM(d0,i,j)  = 0;
            DIRECT_A2D_ELEM(d45,i,j) = 0;
            DIRECT_A2D_ELEM(d90,i,j) = 0;
            DIRECT_A2D_ELEM(d135,i,j)= 0;
        }
        else
        {
            DIRECT_A2D_ELEM(d0,i,j)  = std::sqrt(2)*std::abs(A2D_ELEM(im,i+1,j)-A2D_ELEM(im,i-1,j));
            DIRECT_A2D_ELEM(d45,i,j) = std::abs(A2D_ELEM(im,i+1,j-1)-A2D_ELEM(im,i-1,j+1));
            DIRECT_A2D_ELEM(d90,i,j) = std::sqrt(2)*std::abs(A2D_ELEM(im,i,j+1)-A2D_ELEM(im,i,j-1));
            DIRECT_A2D_ELEM(d135,i,j)= std::abs(A2D_ELEM(im,i+1,j+1)-A2D_ELEM(im,i-1,j-1));
        }
    }

    mask.setXmippOrigin();

    for (int i = -std::abs(wSize); i <= std::abs(wSize); i++ )
        for (int j = -std::abs(wSize); j <= std::abs(wSize); j++ )
            A2D_ELEM(mask,i,j) = double(1)/ double(wSize*wSize+1);


    convolutionFFT(d0,mask,d0);
    convolutionFFT(d45,mask,d45);
    convolutionFFT(d135,mask,d135);
    convolutionFFT(d90,mask,d90);

    /*
     *     b=0.5*(D0-D90);
     *     c=-0.5*(D45-D135); %hay que tener en cuenta que la "y" esta downwards
     *     Or=0.5*atan2(c,b);
     *     orn = mod(Or+pi/2,pi);
     *    ornMod=mat2gray(abs(c+1i*b));
     */

    STARTINGX(ROI)=STARTINGY(ROI)=0;
    double tempx = 0;
    double tempy = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        if (A2D_ELEM(ROI,i,j))
        {
            tempx =  0.5*(A2D_ELEM(d0,i,j)-A2D_ELEM(d90,i,j));
            tempy = -0.5*(A2D_ELEM(d45,i,j)-A2D_ELEM(d135,i,j));
            A2D_ELEM(orn,i,j)   =  std::fmod((0.5*std::atan2(tempy,tempx))+double(PI)/2,double(PI));
            A2D_ELEM(ornMod,i,j) = std::sqrt(std::pow(tempx,2)+std::pow(tempy,2));
        }
        else
        {
            A2D_ELEM(orn,i,j) = 0;
            A2D_ELEM(ornMod,i,j) = 0;
        }
    }

    //we set again the ROI in the xmipp origin
    ROI.setXmippOrigin();
    orMap = orn;
    orModMap = ornMod;
}

void normalize(FourierTransformer &ftrans, MultidimArray<double> & im, MultidimArray<double > & imN,  MultidimArray<double > & imModMap,
		double R, double S, MultidimArray<bool> & ROI)
{
    // H is an Annular filter with radius=R and sigma=S and a Gaussian DC filter with sigma=1
    MultidimArray< std::complex<double> > H;
    H.resizeNoCopy(im);

    im.setXmippOrigin();
    H.setXmippOrigin();

    MultidimArray<std::complex<double> > fftIm, imComplex;
    typeCast(im, imComplex);

    //Fourier Transformer
    ftrans.FourierTransform(imComplex, fftIm, false);

    double temp = 0;
    double K=1.0/(2*S*S);
    for (int i=STARTINGY(im); i<=FINISHINGY(im); ++i)
    {
    	double di=i;
    	double i2=di*di;
    	for (int j=STARTINGX(im); j<=FINISHINGX(im); ++j)
		{
        	double dj=j;
        	double j2=dj*dj;
        	double i2_j2=i2+j2;
			// temp= std::exp(-std::pow((std::sqrt(std::pow((double)i,2)+std::pow((double)j,2))-R),2)/(2*std::pow(S,2)))*(1-(std::exp((-1)*(std::pow(double(i),2) + std::pow(double(j),2)) /(2*1))));
			temp= std::exp(-std::pow((std::sqrt(i2_j2)-R),2)*K)*(1-(std::exp(-0.5*(i2_j2))));
			double *ptr=(double*)&A2D_ELEM(H,i,j);
			*ptr=*(ptr+1)=temp;
		}
    }

    CenterFFT(H,false);
    fftIm *= H;
    ftrans.inverseFourierTransform();

    //output of the program
    imN.setXmippOrigin();
    imModMap.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
        A2D_ELEM(imN,i,j) = A2D_ELEM(imComplex,i,j).real();

    SPTH(ftrans,imN,H);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        if (A2D_ELEM(ROI,i,j))
        {
            temp = std::abs(A2D_ELEM(H,i,j));
            A2D_ELEM(imModMap,i,j) = std::sqrt(std::pow(temp,2)+std::pow(A2D_ELEM(imN,i,j),2));
            A2D_ELEM(imN,i,j)      = std::cos(std::atan2(temp, A2D_ELEM(imN,i,j)));
        }
        else
        {
            A2D_ELEM(imModMap,i,j) = 0;
            A2D_ELEM(imN,i,j)      = 0;
        }
    }

    STARTINGX(imN)=STARTINGY(imN)=0;
}

void normalizeWB(MultidimArray<double> & im, MultidimArray<double > & imN,  MultidimArray<double > & imModMap, double rmax, double rmin, MultidimArray<bool> & ROI)
{

    // H is an Annular bandpass filter.
    MultidimArray< std::complex<double> > H;
    H.resizeNoCopy(im);

    im.setXmippOrigin();
    H.setXmippOrigin();

    MultidimArray<std::complex<double> > fftIm, imComplex;
    typeCast(im, imComplex);

    //Fourier Transformer
    FourierTransformer ftrans(FFTW_BACKWARD);
    ftrans.FourierTransform(imComplex, fftIm, false);

    double rang = (rmax-rmin)/2;
    //Inside rang we assume that there will a range of fringes per field from 2 to 10.
    double freq2 = XSIZE(im)/(rang/1);
    double freq1 = XSIZE(im)/(rang/15);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        double r2=i*i+j*j;
        double temp= (1/(1+std::exp(((std::sqrt(r2)-freq1))/(10))))*
                     (1-(std::exp(-r2 /(2*freq2*freq2))));
        A2D_ELEM(H,i,j) = std::complex<double>(temp,temp);
    }

    CenterFFT(H,false);
    fftIm *= H;
    ftrans.inverseFourierTransform();

    //output of the program
    imN.setXmippOrigin();
    imModMap.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        A2D_ELEM(imN,i,j) = A2D_ELEM(imComplex,i,j).real();
    }

    SPTH(ftrans,imN,H);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        if (A2D_ELEM(ROI,i,j))
        {
            double temp = std::abs(A2D_ELEM(H,i,j));
            A2D_ELEM(imModMap,i,j) = std::sqrt(std::pow(temp,2)+std::pow(A2D_ELEM(imN,i,j),2));
            A2D_ELEM(imN,i,j)      = std::cos(std::atan2(temp, A2D_ELEM(imN,i,j)));
        }
        else
        {
            A2D_ELEM(imModMap,i,j) = 0;
            A2D_ELEM(imN,i,j)      = 0;
        }
    }

    STARTINGX(imN)=STARTINGY(imN)=0;
}

void normalizeWB2(MultidimArray<double> & im, MultidimArray<double > & imN,  MultidimArray<double > & imModMap, double rmax, double rmin, MultidimArray<bool> & ROI)
{

    // H is an Annular bandpass filter.
    MultidimArray< std::complex<double> > H;
    H.resizeNoCopy(im);

    im.setXmippOrigin();
    H.setXmippOrigin();

    MultidimArray<std::complex<double> > fftIm, imComplex;
    typeCast(im, imComplex);

    //Fourier Transformer
    FourierTransformer ftrans(FFTW_BACKWARD);
    ftrans.FourierTransform(imComplex, fftIm, false);

    double rang = (rmax-rmin)/2;
    //Inside rang we assume that there will a range of fringes per field from 2 to 10.
    double freq2 = XSIZE(im)/(rang/1);
    double freq1 = XSIZE(im)/(rang/15);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        double r2=i*i+j*j;
        double temp= (1/(1+std::exp(((std::sqrt(r2)-freq1))/(10))))*
                     (1-(std::exp(-r2 /(2*freq2*freq2))));
        A2D_ELEM(H,i,j) = std::complex<double>(temp,temp);
    }

    CenterFFT(H,false);
    fftIm *= H;
    ftrans.inverseFourierTransform();

    //output of the program
    imN.setXmippOrigin();
    imModMap.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        A2D_ELEM(imN,i,j) = A2D_ELEM(imComplex,i,j).real();
    }

    SPTH(ftrans,imN,H);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        if (A2D_ELEM(ROI,i,j))
        {
            double temp = std::abs(A2D_ELEM(H,i,j));
            A2D_ELEM(imModMap,i,j) = std::sqrt(std::pow(temp,2)+std::pow(A2D_ELEM(imN,i,j),2));
            A2D_ELEM(imN,i,j)      = std::cos(std::atan2(temp, A2D_ELEM(imN,i,j)));
        }
        else
        {
            A2D_ELEM(imModMap,i,j) = 0;
            A2D_ELEM(imN,i,j)      = 0;
        }
    }

    STARTINGX(imN)=STARTINGY(imN)=0;
}

class PointQuality
{
public:
    int i;
    int j;
    double quality;
    const bool operator <(const PointQuality& p) const
    {
        return (quality < p.quality);
    }
};

void direction(const MultidimArray<double> & orMap, MultidimArray<double> & qualityMap, double lambda, int size, MultidimArray<double> & dirMap, int x, int y)
{
    //First we perform some setup stuff
    double minQuality = 0;
    int imax, jmax, max_levels=10;


    //We look for the maximun value of qualityMap
    qualityMap.selfABS();
    qualityMap.maxIndex(imax,jmax);

    double maxQualityMapValue = (A2D_ELEM(qualityMap,imax,jmax));

    MultidimArray<double> qualityMapInt = (((qualityMap/maxQualityMapValue)*(max_levels-1)));
    qualityMapInt.selfROUND();

    MultidimArray<double> ds, dc, px, py;
    sincos(orMap,ds,dc);
    px.resizeNoCopy(orMap);
    py.resizeNoCopy(orMap);

    //In hx and hy we store the pixels to process attending to their quality value
    Matrix2D<double> G(2,2), invG(2,2);
    Matrix1D<double> b(2),R(2);

    MultidimArray<bool> processed;
    processed.initZeros(orMap);

    // Process the first point
    size_t i=x;
    size_t j=y;

    A2D_ELEM(dirMap,i,j) = std::atan2(A2D_ELEM(ds,i,j),A2D_ELEM(dc,i,j));
    A2D_ELEM(processed,i,j) = true;
    A2D_ELEM(px,i,j) = -A2D_ELEM(ds,i,j);
    A2D_ELEM(py,i,j) =  A2D_ELEM(dc,i,j);

    PointQuality p;
    p.i=i;
    p.j=j;
    p.quality=A2D_ELEM(qualityMapInt,i,j);
    std::priority_queue<PointQuality> queueToProcess;
    queueToProcess.push(p);

    while (!queueToProcess.empty())
    {
        p=queueToProcess.top();
        i=p.i;
        j=p.j;
        queueToProcess.pop();

        size_t indi[8] = {i-1,i-1,i-1,i,i,i+1,i+1,i+1};
        size_t indj[8] = {j-1,j,j+1,j-1,j+1,j-1,j,j+1};

        for(int k = 0; k< 8 ; k++)
        {

            size_t ni=indi[k];
            size_t nj=indj[k];

            if ( (!A2D_ELEM(processed,ni,nj)) && 
                  (A2D_ELEM(qualityMap,ni,nj) > minQuality) && 
                  ((ni-size)>0) && ((nj-size)>0) && 
                  ((ni+size)<YSIZE(processed)-1) && 
                  ((nj+size)<XSIZE(processed)-1) )
            {
                double G11=0, G12=0, G22=0, b1=0, b2=0;

                for (int li=-size; li <= size; li++)
                {
                    size_t nli=ni+li;
                    for (int lj=-size; lj <= size; lj++)
                    {
                        size_t nlj=nj+lj;
                        if (A2D_ELEM(processed,nli,nlj))
                        {
                            double tempx = A2D_ELEM(ds,nli,nlj);
                            double tempy = A2D_ELEM(dc,nli,nlj);

                            G11 += tempx*tempx+lambda;
                            G12 += tempx*tempy;
                            G22 += tempy*tempy+lambda;

                            b1 += lambda*(A2D_ELEM(px,nli,nlj));
                            b2 += lambda*(A2D_ELEM(py,nli,nlj));
                        }
                    }
                }

                dMij(G,0,0) = G11;
                dMij(G,0,1) = G12;
                dMij(G,1,0) = G12;
                dMij(G,1,1) = G22;

                VEC_ELEM(b,0) = b1;
                VEC_ELEM(b,1) = b2;

                G.inv(invG);
                R = invG*b;

                // Process this point
                A2D_ELEM(dirMap,ni,nj) = std::atan2(-VEC_ELEM(R,1),VEC_ELEM(R,0));
                A2D_ELEM(processed,ni,nj) = true;
                A2D_ELEM(px,ni,nj) = VEC_ELEM(R,0);
                A2D_ELEM(py,ni,nj) = VEC_ELEM(R,1);
                p.i=ni;
                p.j=nj;
                p.quality=A2D_ELEM(qualityMapInt,ni,nj);
                queueToProcess.push(p);
            }

        }
    }
}


void unwrapping(const MultidimArray<double> & wrappedPhase, MultidimArray<double> & qualityMap, double lambda, int size, MultidimArray<double> & unwrappedPhase)
{
    //First we perform some setup stuff
    double minQuality = 0.05;

    int imax, jmax, i, j, max_levels=10;

    //We look for the maximun value of qualityMap
    qualityMap.selfABS();
    qualityMap.maxIndex(imax,jmax);

    double maxQualityMapValue = A2D_ELEM(qualityMap,imax,jmax);
    MultidimArray<double> qualityMapInt = (((qualityMap/maxQualityMapValue)*(max_levels-1)));
    qualityMapInt.selfROUND();

    Matrix2D<double> gaussian;
    gaussian.initGaussian(2*size+1,1.5);

    MultidimArray<bool> processed;
    processed.resizeNoCopy(wrappedPhase);
    processed.initZeros();

    //Here appears the real processing
    A2D_ELEM(unwrappedPhase,imax,jmax) = A2D_ELEM(wrappedPhase,imax,jmax);
    A2D_ELEM(processed,imax,jmax) = true;

    PointQuality p;
    p.i=imax;
    p.j=jmax;
    p.quality=A2D_ELEM(qualityMapInt,imax,jmax);
    std::priority_queue<PointQuality> queueToProcess;
    queueToProcess.push(p);

    //predictor and corrector
    double pred = 0;
    double cor = 0;
    double uw = 0;
    double wp = 0;
    double g = 0;
    int n = 0;
    double t = 0;
    double norm = 0;
    double up = 0;

    while(!queueToProcess.empty())
    {
        p=queueToProcess.top();
        i=p.i;
        j=p.j;
        queueToProcess.pop();

        size_t indi[8] = {i-1,i-1,i-1,i,i,i+1,i+1,i+1};
        size_t indj[8] = {j-1,j,j+1,j-1,j+1,j-1,j,j+1};

        for(int k = 0; k< 8 ; k++)
        {
            size_t ni=indi[k];
            size_t nj=indj[k];
            if ( (!A2D_ELEM(processed,ni,nj)) && (A2D_ELEM(qualityMap,ni,nj) > minQuality) && ((ni-size)>0) && ((nj-size)>0) && ((ni+size)<YSIZE(processed)-1)
                 && ((nj+size)<XSIZE(processed)-1))
            {
                wp =  A2D_ELEM(wrappedPhase,ni,nj);

                for (int li=-size; li <= size; li++)
                {
                    size_t nli=ni+li;
                    for (int lj=-size; lj <= size; lj++)
                    {
                        size_t nlj=nj+lj;
                        if (A2D_ELEM(processed,nli,nlj))
                        {
                            uw =  A2D_ELEM(unwrappedPhase,nli,nlj);
                            //q  =  A2D_ELEM(qualityMap,nli,nlj);
                            g  =  dMij(gaussian,li+size,lj+size);

                            t = (wp - uw);

                            if ( t > 0 )
                                n = (int)(std::floor(t+3.14159265)/(2*3.14159265));
                            else
                                n = (int)(std::ceil(t-3.14159265)/(2*3.14159265));

                            up = t - (2*3.14159265)*n;

                            //pred += (uw*q*g);
                            //cor  += (up*q*g);
                            //norm += (q*g);

                            pred += (uw*g);
                            cor  += (up*g);
                            norm += (g);

                        }
                    }

                }

                A2D_ELEM(unwrappedPhase,ni,nj) = pred/norm + (lambda*cor)/norm;
                A2D_ELEM(processed,ni,nj) = true;

                norm = 0;
                pred = 0;
                cor = 0;

                p.i=ni;
                p.j=nj;
                p.quality=A2D_ELEM(qualityMapInt,ni,nj);
                queueToProcess.push(p);

            }
        }
    }
}

void demodulate(MultidimArray<double> & im, double lambda, int size, int x, int y, int rmin, int rmax,
                                  MultidimArray<double> & phase, MultidimArray<double> & mod, Matrix1D<double> & coeffs, int verbose)
{
    //Initial Setup
    MultidimArray< double > In, orMap, orModMap, dir, wphase;
    MultidimArray< bool > ROI;
    MultidimArray< std::complex<double> > sph;

    In.resizeNoCopy(im);
    orMap.resizeNoCopy(im);
    orModMap.resizeNoCopy(im);
    dir.resizeNoCopy(im);
    wphase.resizeNoCopy(im);
    sph.resizeNoCopy(im);
    ROI.resizeNoCopy(im);
    mod.resizeNoCopy(im);
    phase.resizeNoCopy(im);

    //First we define the Region of Interesting to perform all the operations inside these regions
    ROI.setXmippOrigin();
    dir.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( (temp > rmin) &&  (temp < rmax) )
        {
            A2D_ELEM(ROI,i,j)= true;
            A2D_ELEM(dir,i,j)= -std::atan2(i,j);
        }
        else
            A2D_ELEM(ROI,i,j)= false;
    }

    STARTINGX(dir)=STARTINGY(dir)=0;

    //We obtain previous necessary maps
    //Normalized version of im and modulation map
    normalizeWB(im,In,mod, rmax, rmin, ROI);

    int imax, jmax;
    mod.maxIndex(imax,jmax);
    mod = mod/A2D_ELEM(mod,imax,jmax);

    Image<double> save;
    if (verbose == 1)
    {
        save()=In;
        save.write("PPP1.xmp");
    }

    STARTINGX(mod)=STARTINGY(mod)=0;
    //Orientation map of the fringes

    /*orMinDer(In, orMap, orModMap, size, ROI);
    if (verbose == 2)
{
        save()=orMap;
        save.write("PPP21.xmp");
        save()=orModMap;
        save.write("PPP22.xmp");
}
    */

    //We obtain the direction from the orientation map
    //direction(orMap, orModMap, lambda, size, dir, x, y);
    if (verbose == 2)
    {
        save()=dir;
        save.write("PPP2.xmp");
    }
    //Spiral transform of the normalized image
    FourierTransformer ftrans(FFTW_BACKWARD);
    SPTH(ftrans,In,sph);
    STARTINGX(sph)=STARTINGY(sph)=0;
    STARTINGX(In)=STARTINGY(In)=0;
    STARTINGX(ROI)=STARTINGY(ROI)=0;

    // We obtain the wrapped phase
    std::complex<double> ci = std::complex<double>(0,1.0);
    std::complex<double> temp;

    //tempTheta is the Theta coordinate in polar. This is because we want to subtract in ROI the crux of the CTF
    double tempTheta=0;
    //val is a number that establish the quantity we want to subtract in the phase to avoid the crux
    double val = 0.1;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(im)
    {
        if (A2D_ELEM(ROI,i,j))
        {
            temp = -ci*std::exp(ci*A2D_ELEM(dir,i,j))*A2D_ELEM(sph,i,j);
            A2D_ELEM(wphase,i,j) = std::atan2(temp.real(),A2D_ELEM(In,i,j));
            A2D_ELEM(orModMap,i,j) = std::sqrt( temp.real()*temp.real() + A2D_ELEM(In,i,j)*A2D_ELEM(In,i,j));

            //tempTheta = std::atan2(j-(int)((float) (XSIZE(ROI)) / 2.0), i-(int)((float) (XSIZE(ROI)) / 2.0));
            tempTheta = A2D_ELEM(dir,i,j);
            if (  !(((tempTheta > val) | (tempTheta < -val)) && !((tempTheta > PI-val) || (tempTheta < -PI+val)) &&
                    !( !(tempTheta > PI/2+val) && (tempTheta > PI/2-val)) && !( !(tempTheta > -PI/2+val) && (tempTheta > -PI/2-val))) )
            {
                A2D_ELEM(ROI,i,j) = false;
            }

        }
        else
        {
            A2D_ELEM(orModMap,i,j) = 0;
        }
    }

    if (verbose == 3)
    {
        save()=wphase;
        save.write("PPP3.xmp");
        Image<bool> saveROI;
        saveROI() = ROI;
        saveROI.write("PPP_roi.xmp");
    }

    //unwrapping(wphase, orModMap, lambda, size, phase);
    unwrapping(wphase, mod, lambda, size, phase);

    if (verbose == 4)
    {
        save()=phase;
        save.write("PPP4.xmp");
    }

    ROI.setXmippOrigin();
    mod.setXmippOrigin();
    Matrix1D<int> coefsInit(VEC_XSIZE(coeffs));

    for (size_t i=0; i<VEC_XSIZE(coeffs); i++)
        VEC_ELEM(coefsInit,i) = (int)VEC_ELEM(coeffs,i);

    PolyZernikes polynom;
    MultidimArray< double > mod2;
    mod2.resizeNoCopy(im);
    mod2.initConstant(1.0);
    polynom.fit(coefsInit,phase,mod2,ROI,0);
    polynom.fit(coefsInit,phase,mod2,ROI,verbose);

    int index = 0;
    for (size_t i=0; i<VEC_XSIZE(coeffs); i++)
    {
        if ( VEC_ELEM(coeffs,i) != 0)
        {
            VEC_ELEM(coeffs,i) = VEC_ELEM(polynom.fittedCoeffs,index);
            index++;
        }

        else
            VEC_ELEM(coeffs,i) = 0;
    }

    STARTINGX(im)=STARTINGY(im)=0;
    STARTINGX(phase)=STARTINGY(phase)=0;
    STARTINGX(mod)=STARTINGY(mod)=0;

    Image<bool> save2;
    if (verbose > 5)
    {
        save()=im;
        save.write("PPP0.xmp");
        save()=In;
        save.write("PPP1.xmp");
        save()=orMap;
        save.write("PPP21.xmp");
        save()=orModMap;
        save.write("PPP22.xmp");
        save()=dir;
        save.write("PPP3.xmp");
        save()=wphase;
        save.write("PPP4.xmp");
        save()=phase;
        save.write("PPP5.xmp");
        save()=mod;
        save.write("PPP6.xmp");
        save2()= ROI;
        save2.write("PPP7.xmp");
    }
}

void demodulate2(MultidimArray<double> & im, double lambda, int size, int x, int y, int rmin, int rmax,
                                   Matrix1D<double> & coeffs, int verbose)
{
    //Initial Setup :
    // In : Normalized image after using methods normalizeWB or normalize
    // orModMap : modulation map of im using
    //TODO : subtract orModMap
    // dir : direction map. It can be obtained using direction method but in our case, we impose the direction map shape
    // wphase : wrapped phase
    // phase : absolute phase computed from wphase and method unwrapping
    // mod : modulation map

    MultidimArray< double > In, mod, dir, wphase, phase;
    MultidimArray< bool > ROI;
    MultidimArray< std::complex<double> > sph;

    In.resizeNoCopy(im);
    mod.resizeNoCopy(im);
    dir.resizeNoCopy(im);
    wphase.resizeNoCopy(im);
    sph.resizeNoCopy(im);
    ROI.resizeNoCopy(im);
    phase.resizeNoCopy(im);

    //First we define the Region of Interesting to perform all the operations inside these regions
    ROI.setXmippOrigin();
    dir.setXmippOrigin();

    //We obtain the ROI from rmax and rmin parameters and the direction map
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( (temp > rmin) &&  (temp < rmax) )
        {
            A2D_ELEM(ROI,i,j)= true;
            A2D_ELEM(dir,i,j)= -std::atan2(i,j);
        }
        else
            A2D_ELEM(ROI,i,j)= false;
    }
    STARTINGX(dir)=STARTINGY(dir)=0;

    //We obtain the normalized version of im and the modulation map mod
    normalizeWB(im,In,mod, rmax, rmin, ROI);
    STARTINGX(mod)=STARTINGY(mod)=0;

    //Here we make the max of mod to be 1
    int imax, jmax;
    mod.maxIndex(imax,jmax);
    mod = mod/A2D_ELEM(mod,imax,jmax);

    Image<double> save;
    if (verbose == 1)
    {
        save()=In;
        save.write("PPP1.xmp");
    }

    if (verbose == 2)
    {
        save()=dir;
        save.write("PPP2.xmp");
    }
    //Spiral transform of the normalized image
    FourierTransformer ftrans(FFTW_BACKWARD);
    SPTH(ftrans,In,sph);
    STARTINGX(sph)=STARTINGY(sph)=0;
    STARTINGX(In)=STARTINGY(In)=0;
    STARTINGX(ROI)=STARTINGY(ROI)=0;

    // We obtain the wrapped phase
    std::complex<double> ci = std::complex<double>(0,1.0);
    std::complex<double> temp;

    //tempTheta is the Theta magnitude in polar coordinates. This is because we want to subtract in ROI the crux of the CTF
    double tempTheta=0;
    //val is a number that establish the quantity we want to subtract in the phase to avoid the crux
    double val = 0.05;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(im)
    {
        if (A2D_ELEM(ROI,i,j))
        {
            temp = -ci*std::exp(ci*A2D_ELEM(dir,i,j))*A2D_ELEM(sph,i,j);
            A2D_ELEM(wphase,i,j) = std::atan2(temp.real(),A2D_ELEM(In,i,j));
            //A2D_ELEM(mod,i,j) = std::sqrt( temp.real()*temp.real() + A2D_ELEM(In,i,j)*A2D_ELEM(In,i,j));

            tempTheta = std::atan2(j-(int)((float) (XSIZE(ROI)) / 2.0), i-(int)((float) (XSIZE(ROI)) / 2.0));

            if (  !(((tempTheta > val) | (tempTheta < -val)) && !((tempTheta > PI-val) || (tempTheta < -PI+val)) &&
                    !( !(tempTheta > PI/2+val) && (tempTheta > PI/2-val)) && !( !(tempTheta > -PI/2+val) && (tempTheta > -PI/2-val))) )
            {
                A2D_ELEM(ROI,i,j) = false;
            }
        }
        else
        {
            A2D_ELEM(mod,i,j) = 0;
        }
    }

    if (verbose == 3)
    {
        save()=wphase;
        save.write("PPP3.xmp");
    }

    unwrapping(wphase, mod, lambda, size, phase);

    if (verbose == 4)
    {
        save()=phase;
        save.write("PPP4.xmp");
    }

    ROI.setXmippOrigin();
    mod.setXmippOrigin();
    Matrix1D<int> coefsInit(VEC_XSIZE(coeffs));

    for (size_t i=0; i<VEC_XSIZE(coeffs); i++)
        VEC_ELEM(coefsInit,i) = (int)VEC_ELEM(coeffs,i);

    PolyZernikes polynom;
    MultidimArray< double > onesMatrix;
    onesMatrix.resizeNoCopy(im);
    onesMatrix.initConstant(1.0);
    polynom.fit(coefsInit,phase,onesMatrix,ROI,0);
    polynom.fit(coefsInit,phase,onesMatrix,ROI,verbose);

    int index = 0;
    for (size_t i=0; i<VEC_XSIZE(coeffs); i++)
    {
        if ( VEC_ELEM(coeffs,i) != 0)
        {
            VEC_ELEM(coeffs,i) = VEC_ELEM(polynom.fittedCoeffs,index);
            index++;
        }

        else
            VEC_ELEM(coeffs,i) = 0;
    }


    STARTINGX(im)=STARTINGY(im)=0;
    STARTINGX(phase)=STARTINGY(phase)=0;
    STARTINGX(mod)=STARTINGY(mod)=0;

    if (verbose > 5)
    {
        save()=im;
        save.write("PPP0.xmp");
        save()=In;
        save.write("PPP1.xmp");
        //save()=orModMap;
        //save.write("PPP2.xmp");
        save()=dir;
        save.write("PPP3.xmp");
        save()=wphase;
        save.write("PPP4.xmp");
        save()=phase;
        save.write("PPP5.xmp");
        save()=mod;
        save.write("PPP6.xmp");

        Image<bool> save2;
        save2()= ROI;
        save2.write("PPP7.xmp");
    }
}

void firsPSDZero(MultidimArray<double> & enhancedPSD, Matrix1D<double> & xPoints,Matrix1D<double> & yPoints,
                                   double rmin, double rmax, int numAngles, int verbose)
{
    enhancedPSD.setXmippOrigin();
    Histogram1D hist;
    compute_hist(enhancedPSD, hist, 200);

    double eff0 = hist.percentil(0.1);
    double effF = hist.percentil(98);
    double thrs = (eff0 + effF)*0.5;

    int index = 0;
    double nAngles = 0;

    while (nAngles < 2*PI)
    {
        double minX = (rmin/2)*std::cos(nAngles);
        double minY = (rmin/2)*std::sin(nAngles);
        double maxX = (rmax/2)*std::cos(nAngles);
        double maxY = (rmax/2)*std::sin(nAngles);
        int numPts = numAngles;
        double x = minX;
        double y = 0;
        double step = (maxX-minX)/numPts;
        double z;

        for (int i = 0; i < numPts; i++)
        {
            y = ( ((maxY-minY)/(maxX-minX))*(x-minX)+minY );
            z = A2D_ELEM(enhancedPSD,int(y),int(x));

            if (z < thrs)
            {
                VEC_ELEM(xPoints,index)=x;
                VEC_ELEM(yPoints,index)=y;
                continue;
            }

            x = (x+step);
        }

        index++;
        nAngles += 2*PI/numAngles;
    }


    double x0, y0, majorAxis, minorAxis, ellipseAngle;
    fitEllipse(xPoints, yPoints, x0, y0, majorAxis,minorAxis,ellipseAngle);

    //EXTERNAL FORCE FIELD CALCULATION
    //We perform the gradient of enhancedPSD to calculate the external force field (balloon force)
    //Calculation of the external force field normalized to maximum magnitude 1

    size_t NR, NC,NZ, Ndim;
    enhancedPSD.getDimensions(NR,NC,NZ,Ndim);
    MultidimArray<double > fx, fy;
    fx.resizeNoCopy(enhancedPSD);
    fy.resizeNoCopy(enhancedPSD);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(enhancedPSD)
    {
        if ( (i==STARTINGX(enhancedPSD)) || (i== (FINISHINGX(enhancedPSD)-1)) || (j == STARTINGY(enhancedPSD)) || (j== (FINISHINGY(enhancedPSD)-1)) )
        {
            A2D_ELEM(fx,i,j)  = 0;
            A2D_ELEM(fy,i,j) = 0;
        }
        else
        {
            A2D_ELEM(fx,i,j)  = std::sqrt(2)*std::abs(A2D_ELEM(enhancedPSD,i+1,j)-A2D_ELEM(enhancedPSD,i-1,j));
            A2D_ELEM(fy,i,j)  = std::sqrt(2)*std::abs(A2D_ELEM(enhancedPSD,i,j+1)-A2D_ELEM(enhancedPSD,i,j-1));
        }
    }

    //We calculate accumulative derivatives in the area 2*wSize+1
    int wSize = 5;
    MultidimArray<double> mask;
    mask.resizeNoCopy(enhancedPSD);
    for (int i = -std::abs(wSize); i <= std::abs(wSize); i++ )
        for (int j = -std::abs(wSize); j <= std::abs(wSize); j++ )
            A2D_ELEM(mask,i,j) = double(1)/ double(wSize*wSize+1);

    convolutionFFT(fx,mask,fx);
    convolutionFFT(fy,mask,fy);

    //From the derivatives we calculate kappa:
    MultidimArray<double> absFx;
    absFx.resizeNoCopy(fx);
    //absFx= fx.selfABS();
    /*absFx.maxIndex(imax,jmax);
    double maxFx = A2D_ELEM(fx,imax,jmax);

    MultidimArray<double> absFy = (fy.selfABS());
    absFy.maxIndex(imax,jmax);
    double maxFy = A2D_ELEM(fy,imax,jmax);
    */
    ///////////////////////
    //

    STARTINGX(enhancedPSD)=STARTINGY(enhancedPSD)=0;
    //enhancedPSD

    if (verbose != 0)
    {
        Image<double> save;

        save() = enhancedPSD;
        save.write("PPP0.xmp");

        save() = fx;
        save.write("PPP1.xmp");

        save() = fy;
        save.write("PPP2.xmp");
    }
}

//Obtained from: http://www.mathworks.com/matlabcentral/fileexchange/15125-fitellipse-m/content/fitellipse.m
void fitEllipse(Matrix1D<double> & xPts, Matrix1D<double> & yPts, double & x0, double & y0, double & majorAxis,
                                  double & minorAxis, double & ellipseAngle)
{
    Matrix2D<double> B(VEC_XSIZE(xPts),5), A(2,2);
    Matrix1D<double> V(5), squareXPts(VEC_XSIZE(xPts));

    for (size_t nRow = 0; nRow < VEC_XSIZE(xPts); nRow++)
    {
        double tempX = VEC_ELEM(xPts,nRow);
        double tempY = VEC_ELEM(yPts,nRow);

        dMij(B,nRow,0) = 2*tempX*tempY;
        dMij(B,nRow,1) = tempY*tempY-tempX*tempX;
        dMij(B,nRow,2) = tempX;
        dMij(B,nRow,3) = tempY;
        dMij(B,nRow,4) = 1;
        VEC_ELEM(squareXPts,nRow) = -tempX*tempX;

    }

    Matrix1D<double> bv(2);
    Matrix2D<double> C =(B.inv());
    V = C*(squareXPts);

    dMij(A,0,0) = 1-VEC_ELEM(V,1);
    dMij(A,0,1) = VEC_ELEM(V,0);
    dMij(A,1,0) = VEC_ELEM(V,0);
    dMij(A,1,1) = VEC_ELEM(V,1);
    VEC_ELEM(bv,0) = VEC_ELEM(V,2);
    VEC_ELEM(bv,1) = VEC_ELEM(V,3);

    double c = VEC_ELEM(V,4);

    Matrix2D<double> Ainv(2,2), u, v;
    Matrix1D<double> w;
    SPEED_UP_temps0;
    M2x2_INV(Ainv,A);
    svdcmp(A,u,w,v);

    Matrix1D<double> t = -0.5*Ainv*bv;

    Matrix1D<double> temp1 = (t*(A * t));//c_h = t.transpose() * A * t + bv.transpose() * t + c;
    Matrix1D<double> temp2 = bv.transpose() * t;
    double c_h = temp1.sum()+temp2.sum()+c;
    x0 = VEC_ELEM(t,0);
    y0 = VEC_ELEM(t,1);
    majorAxis = std::sqrt(-c_h / VEC_ELEM(w,0));
    minorAxis = std::sqrt(-c_h / VEC_ELEM(w,1));
    ellipseAngle = std::atan2(dMij(u,0,0),-dMij(u,1,0));

    if ( (VEC_ELEM(w,0)==0) || (VEC_ELEM(w,1)==0) || (VEC_ELEM(w,0)*c_h>0) || ((VEC_ELEM(w,1)*c_h>0)) )
    {
    	majorAxis=0;
    	minorAxis=0;
    }


    double sEllipseAngle, cEllipseAngle;
    sincos(-ellipseAngle,&sEllipseAngle,&cEllipseAngle);

    double angle = 0, deltaAngle=(2*PI)/VEC_XSIZE(xPts);
    for (size_t nPoint = 0; nPoint < VEC_XSIZE(xPts); nPoint++)
    {
    	double sAngle, cAngle;
        sincos(angle,&sAngle,&cAngle);
        //We impose that the origin always is zero and because this whe do not sum it
        double K1=majorAxis*cAngle;
        double K2=minorAxis*sAngle;
        VEC_ELEM(xPts,nPoint) = K1*cEllipseAngle - K2*sEllipseAngle+x0;
        VEC_ELEM(yPts,nPoint) = K1*sEllipseAngle + K2*cEllipseAngle+y0;
        angle += deltaAngle;
    }
}

void fitEllipse(MultidimArray<double> & normImag, double & x0, double & y0, double & majorAxis,
                                  double & minorAxis, double & ellipseAngle, size_t & area)
{
    area = (size_t)normImag.sum();
    Matrix1D<double> xPts(area);
    Matrix1D<double> yPts(area);
    int index = 0;

    FOR_ALL_ELEMENTS_IN_ARRAY2D(normImag)
    {
        if ( A2D_ELEM(normImag,i,j)==1 )
        {
        	VEC_ELEM(xPts,index)=j;
        	VEC_ELEM(yPts,index)=i;
        	index++;
        }
    }

    fitEllipse(xPts,yPts,x0,y0,majorAxis,minorAxis,ellipseAngle);
}

//Equivalent function in JAVA in : xmipp.svn/java/src/xmipp/viewer/ctf/EllipseCTF.java
void calculateDefocus(double & defocusU,double & defocusV, double majorAxis, double minorAxis, double Q0, double lambda,
                                        double Cs, double imgSize, double Ts)
{
    double R = majorAxis / (imgSize * Ts);
    double a = lambda * R * R;
    defocusU = -std::atan(-Q0) / (PI * a) - (1 / a) - 0.5 * (Cs * a * lambda);

    R = minorAxis / (imgSize * Ts);
    a = lambda * R * R;
    defocusV = -std::atan(-Q0) / (PI * a) - (1 / a) - 0.5 * (Cs * a * lambda);

}



