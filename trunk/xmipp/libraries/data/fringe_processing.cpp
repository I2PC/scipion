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
#include "xmipp_polynomials.h"
#include "xmipp_image.h"
#include "multidim_array.h"
#include "xmipp_funcs.h"
#include "xmipp_fftw.h"


//This member function simulates different types of fringe patterns.
// im is the resultant multidimarray
// type is the type of fringe pattern that we want to generate. Possible values are:
// SIMPLY_OPEN_FRINGES, SIMPLY_CLOSED_FRINGES, COMPLEX_OPEN_FRINGES, COMPLEX_CLOSED_FRINGES
// xdim is the image dimension in x axis
// ydim is the image dimension in y axis
void FringeProcessing::simulPattern(MultidimArray<double> & im, enum FP_TYPE type, int xdim, int ydim, double noiseLevel, const double freq, Matrix1D<int> coefs)
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
        polynom.zernikePols(coefs,im);
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
void FringeProcessing::SPTH(MultidimArray<double> & im, MultidimArray< std::complex <double> > & imProc)
{
    im.setXmippOrigin();
    MultidimArray<std::complex<double> > H, fftIm, imComplex;
    typeCast(im, imComplex);

    // Fourier Transformer
    FourierTransformer ftrans(FFTW_BACKWARD);
    ftrans.FourierTransform(imComplex, fftIm, false);

    H.resizeNoCopy(fftIm);
    H.setXmippOrigin();

    //std::complex<double> compTemp(0, 0);
    //i -> for exterior filas o Y
    //j -> for interior columns o X
    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        A2D_ELEM(H,i,j) = (std::complex<double>(j,i))/std::sqrt(std::pow((double)i,2)+(std::pow((double)j,2)));
        //We subtract the singular value in i=0, j=0
        if ( (i==0) && (j==0) )
            A2D_ELEM(H,i,j) = 0;
    }

    CenterFFT(H,false);
    fftIm *= H;
    ftrans.inverseFourierTransform();
    //Here in the Matlab code there is a complex conjugate s
    imProc = imComplex;

}

void FringeProcessing::orMinDer(const MultidimArray<double> & im, MultidimArray<double > & orMap, MultidimArray<double > & orModMap, int wSize)
{
    int NR, NC,NZ;
    size_t NDim;
    im.getDimensions(NC,NR,NZ,NDim);

    if ( (NZ!=1) || (NDim!=1) )
        REPORT_ERROR(ERR_MULTIDIM_DIM,(std::string)"ZDim and NDim has to be equals to one");

    MultidimArray<double > d0,d45,d90,d135,mask,orn,ornMod;
    d0.resizeNoCopy(im);
    d45.resizeNoCopy(im);
    d90.resizeNoCopy(im);
    d135.resizeNoCopy(im);
    mask.resizeNoCopy(im);
    orn.resizeNoCopy(im);
    ornMod.resizeNoCopy(im);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        A2D_ELEM(mask,i,j)  = 0;

        if ( (i==0) || (i== (NR-1)) || (j == 0) || (j== (NC-1)) )
        {
            A2D_ELEM(d0,i,j)  = 0;
            A2D_ELEM(d45,i,j) = 0;
            A2D_ELEM(d90,i,j) = 0;
            A2D_ELEM(d135,i,j)= 0;
        }
        else
        {
            A2D_ELEM(d0,i,j)  = std::sqrt(2)*std::abs(A2D_ELEM(im,i+1,j)-A2D_ELEM(im,i-1,j));
            A2D_ELEM(d45,i,j) = std::abs(A2D_ELEM(im,i+1,j-1)-A2D_ELEM(im,i-1,j+1));
            A2D_ELEM(d90,i,j) = std::sqrt(2)*std::abs(A2D_ELEM(im,i,j+1)-A2D_ELEM(im,i,j-1));
            A2D_ELEM(d135,i,j)= std::abs(A2D_ELEM(im,i+1,j+1)-A2D_ELEM(im,i-1,j-1));
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
    double tempx = 0;
    double tempy = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        tempx =  0.5*(A2D_ELEM(d0,i,j)-A2D_ELEM(d90,i,j));
        tempy = -0.5*(A2D_ELEM(d45,i,j)-A2D_ELEM(d135,i,j));
        A2D_ELEM(orn,i,j)   =  std::fmod((0.5*std::atan2(tempy,tempx))+double(PI)/2,double(PI));
        A2D_ELEM(ornMod,i,j) = std::sqrt(std::pow(tempx,2)+std::pow(tempy,2));
    }

    orMap = orn;
    orModMap = ornMod;
}

void FringeProcessing::normalize(MultidimArray<double> & im, MultidimArray<double > & imN,  MultidimArray<double > & imModMap, double R, double S)
{
    // H is an Annular filter with radius=R and sigma=S and a Gaussian DC filter with sigma=1
    MultidimArray< std::complex<double> > H;
    H.resizeNoCopy(im);

    im.setXmippOrigin();
    H.setXmippOrigin();

    MultidimArray<std::complex<double> > fftIm, imComplex;
    typeCast(im, imComplex);

    //Fourier Transformer
    FourierTransformer ftrans(FFTW_BACKWARD);
    ftrans.FourierTransform(imComplex, fftIm, false);

    double temp = 0;
    std::complex<double> tempCpx;


    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        temp= std::exp(-std::pow((std::sqrt(std::pow((double)i,2)+std::pow((double)j,2))-R),2)/(2*std::pow(S,2)))*(1-(std::exp((-1)*(std::pow(double(i),2) + std::pow(double(j),2)) /(2*1))));
        tempCpx = std::complex<double>(temp,temp);
        A2D_ELEM(H,i,j) = tempCpx;
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

    SPTH(imN,H);
    sph = H;

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        temp = std::abs(A2D_ELEM(H,i,j));
        A2D_ELEM(imModMap,i,j) = std::sqrt(std::pow(temp,2)+std::pow(A2D_ELEM(imN,i,j),2));
        A2D_ELEM(imN,i,j)      = std::cos(std::atan2(temp, A2D_ELEM(imN,i,j)));
    }

    STARTINGX(imN)=STARTINGY(imN)=0;
}

void FringeProcessing::direction(const MultidimArray<double> & orMap, MultidimArray<double> & qualityMap, double lambda, int size, MultidimArray<double> & dirMap)
{
    //First we perform some setup stuff
    int nx = XSIZE(orMap);
    int ny = YSIZE(orMap);
    double minQuality = 0.0;

    Histogram1D hist;
    int no_steps = 10, maxElem = 0, k;
    int imax, jmax, i, j;

    //We look for the maximun value of qualityMap
    qualityMap.selfABS();
    qualityMap.maxIndex(imax,jmax);

    double maxQualityMapValue = A2D_ELEM(qualityMap,imax,jmax);

    //Here we transform the quality map to the no_steps possible values
    MultidimArray<double> qualityMapInt = (((qualityMap/maxQualityMapValue)*(no_steps-1)));
    qualityMapInt.selfROUND();

    //imax = 100;
    //jmax = 100;

    compute_hist(qualityMapInt, hist, no_steps);
    int tempValue = 0;

    //We look for the max element in hist
    FOR_ALL_ELEMENTS_IN_ARRAY1D(hist)
    {
        tempValue = A1D_ELEM(hist,i);
        if (tempValue > maxElem)
            maxElem =  tempValue;
    }

    MultidimArray<double> ds, dc, px, py;
    sincos(orMap,ds,dc);
    px.resizeNoCopy(orMap);
    py.resizeNoCopy(orMap);

    //In hx and hy we store the pixels to process attending to their quality value
    Matrix2D<int> hi(no_steps,maxElem), hj(no_steps,maxElem);
    Matrix1D<int> front(maxElem), final(maxElem);
    Matrix2D<double> G(2,2), invG(2,2);
    Matrix1D<double> b(2),R(2);

    hi.initZeros();
    hj.initZeros();
    front.initZeros();
    final.initZeros();

    MultidimArray<bool> processed;
    processed.resizeNoCopy(orMap);
    processed.initZeros();

    //Here appears the real processing
    A2D_ELEM(dirMap,imax,jmax) = std::atan2(A2D_ELEM(ds,imax,jmax),A2D_ELEM(dc,imax,jmax));
    A2D_ELEM(processed,imax,jmax) = true;
    A2D_ELEM(px,imax,jmax) = -A2D_ELEM(ds,imax,jmax);
    A2D_ELEM(py,imax,jmax) =  A2D_ELEM(dc,imax,jmax);

    int ind = 1;
    double G11 = 0;
    double G12 = 0;
    double G21 = 0;
    double G22 = 0;
    double b1 = 0;
    double b2 = 0;
    int H = 0;
    double tempx,tempy;
    bool tempb;

    i = imax;
    j = jmax;


    while( ind > 0)
    {
        int indi[8] = {i-1,i-1,i-1,i,i,i+1,i+1,i+1};
        int indj[8] = {j-1,j,j+1,j-1,j+1,j-1,j,j+1};

        for(int k = 0; k< 8 ; k++)
            if ( (A2D_ELEM(processed,indi[k],indj[k]) == 0) && (indi[k] > size-1 ) && (indi[k] < nx-size )
                 && (indj[k] > size-1 ) && (indj[k] < ny-size ) &&  (A2D_ELEM(qualityMap,indi[k],indj[k]) > minQuality) )
            {
                for (int li=-size; li <= size; li++)
                    for (int lj=-size; lj <= size; lj++)
                    {
                        tempx = A2D_ELEM(ds,indi[k]+li,indj[k]+lj);
                        tempy = A2D_ELEM(dc,indi[k]+li,indj[k]+lj);
                        tempb = A2D_ELEM(processed,indi[k]+li,indj[k]+lj);

                        G11 += tempx*tempx+lambda*tempb;
                        G12 += tempx*tempy;
                        G22 += tempy*tempy+lambda*tempb;

                        b1 += lambda*(A2D_ELEM(px,indi[k]+li,indj[k]+lj))*tempb;
                        b2 += lambda*(A2D_ELEM(py,indi[k]+li,indj[k]+lj))*tempb;

                    }

                dMij(G,0,0) = G11;
                dMij(G,0,1) = G12;
                dMij(G,1,0) = G12;
                dMij(G,1,1) = G22;

                VEC_ELEM(b,0) = b1;
                VEC_ELEM(b,1) = b2;

                G.inv(invG);
                R = invG*b;

                A2D_ELEM(dirMap,indi[k],indj[k]) = std::atan2(-VEC_ELEM(R,1),VEC_ELEM(R,0));
                A2D_ELEM(processed,indi[k],indj[k]) = true;
                A2D_ELEM(px,indi[k],indj[k]) = VEC_ELEM(R,0);
                A2D_ELEM(py,indi[k],indj[k]) = VEC_ELEM(R,1);


                H = A2D_ELEM(qualityMapInt,indi[k],indj[k]);
                if ( H == 0)
                    H = 1;


                if ( VEC_ELEM(final,H) ==  maxElem )
                    VEC_ELEM(final,H) = 1;
                else
                    VEC_ELEM(final,H) += 1;

                dMij(hi,H,VEC_ELEM(final,H)) = indi[k];
                dMij(hj,H,VEC_ELEM(final,H)) = indj[k];

                if ( VEC_ELEM(front,H) ==  0 )
                    VEC_ELEM(front,H) =  1;


                G11 = 0;
                G12 = 0;
                G21 = 0;
                G22 = 0;
                b1 = 0;
                b2 = 0;

            }

        k = maxElem-1;
        ind = 0;

        //We obtain the next coordinates to process from the histogram looking for the values in hi, hj
        //with highest quality value
        while( (ind == 0) && (k > 0) )
        {
            int temp = VEC_ELEM(front,k);
            ind = temp;
            if (ind > 0)
            {
                i = dMij(hi,k,temp);
                j = dMij(hj,k,temp);

                if ( temp == VEC_ELEM(final,k))
                {
                    VEC_ELEM(front,k) = 0;
                    VEC_ELEM(final,k) = 0;
                }
                else
                    if ( temp == maxElem)
                        VEC_ELEM(front,k) = 1;
                    else
                        VEC_ELEM(front,k) += 1;
            }

            --k;
        }
    }
}


void FringeProcessing::unwrapping(const MultidimArray<double> & wrappedPhase, MultidimArray<double> & qualityMap, double lambda, int size, MultidimArray<double> & unwrappedPhase)
{
    //First we perform some setup stuff
    int nx = XSIZE(wrappedPhase);
    int ny = YSIZE(wrappedPhase);
    double minQuality = 0.01;

    Histogram1D hist;
    int no_steps = 10, maxElem = 0, k;
    int imax, jmax, i, j;

    //We look for the maximun value of qualityMap
    qualityMap.selfABS();
    qualityMap.maxIndex(imax,jmax);

    double maxQualityMapValue = A2D_ELEM(qualityMap,imax,jmax);

    //Here we transform the quality map to the no_steps possible values
    MultidimArray<double> qualityMapInt = (((qualityMap/maxQualityMapValue)*(no_steps-1)));
    qualityMapInt.selfROUND();

    compute_hist(qualityMapInt, hist, no_steps);
    int tempValue = 0;

    //We look for the max element in hist
    FOR_ALL_ELEMENTS_IN_ARRAY1D(hist)
    {
        tempValue = A1D_ELEM(hist,i);
        if (tempValue > maxElem)
            maxElem =  tempValue;
    }

    //In hx and hy we store the pixels to process attending to their quality value
    Matrix2D<int> hi(no_steps,maxElem), hj(no_steps,maxElem);
    Matrix1D<int> front(maxElem), final(maxElem);
    Matrix2D<double> gaussian;
    gaussian.initGaussian(2*size+1,0.6);
    hi.initZeros();
    hj.initZeros();
    front.initZeros();
    final.initZeros();

    MultidimArray<bool> processed;
    processed.resizeNoCopy(wrappedPhase);
    processed.initZeros();

    //Here appears the real processing
    A2D_ELEM(unwrappedPhase,imax,jmax) = A2D_ELEM(wrappedPhase,imax,jmax);
    A2D_ELEM(processed,imax,jmax) = true;

    int ind = 1;
    int H = 0;
    i = imax;
    j = jmax;

    //predictor and corrector
    double pred = 0;
    double cor = 0;
    double uw = 0;
    double wp = 0;
    double q = 0;
    double g = 0;
    int n = 0;
    double t = 0;
    double norm = 0;
    double up = 0;


    while( ind > 0)
    {

        int indi[8] = {i-1,i-1,i-1,i,i,i+1,i+1,i+1};
        int indj[8] = {j-1,j,j+1,j-1,j+1,j-1,j,j+1};

        for(int k = 0; k< 8 ; k++)
            if ( (A2D_ELEM(processed,indi[k],indj[k]) == 0) && (indi[k] > size-1 ) && (indi[k] < nx-size )
                 && (indj[k] > size-1 ) && (indj[k] < ny-size ) &&  (A2D_ELEM(qualityMap,indi[k],indj[k]) > minQuality) )
            {
                wp =  A2D_ELEM(wrappedPhase,indi[k],indj[k]);

                for (int li=-size; li <= size; li++)
                    for (int lj=-size; lj <= size; lj++)
                    {
                        if (A2D_ELEM(processed,indi[k]+li,indj[k]+lj)>0)
                        {
                            uw =  A2D_ELEM(unwrappedPhase,indi[k]+li,indj[k]+lj);
                            q  =  A2D_ELEM(qualityMap,indi[k]+li,indj[k]+lj);
                            g  =  dMij(gaussian,li+size,lj+size);

                            t = (wp - uw);

                            if ( t > 0 )
                                n = std::floor(t+3.14159265)/(2*3.14159265);
                            else
                                n = std::ceil(t-3.14159265)/(2*3.14159265);

                            up = t - (2*3.14159265)*n;

                            pred += (uw*q*g);
                            cor  += (up*q*g);
                            norm += (q*g);
                        }

                    }

                A2D_ELEM(unwrappedPhase,indi[k],indj[k]) = pred/norm + (lambda*cor)/norm;
                A2D_ELEM(processed,indi[k],indj[k]) = true;

                norm = 0;
                pred = 0;
                cor = 0;

                H = A2D_ELEM(qualityMapInt,indi[k],indj[k]);
                if ( H == 0)
                    H = 1;

                if ( VEC_ELEM(final,H) ==  maxElem )
                    VEC_ELEM(final,H) = 1;
                else
                    VEC_ELEM(final,H) += 1;

                dMij(hi,H,VEC_ELEM(final,H)) = indi[k];
                dMij(hj,H,VEC_ELEM(final,H)) = indj[k];

                if ( VEC_ELEM(front,H) ==  0 )
                    VEC_ELEM(front,H) =  1;

            }

        k = maxElem-1;
        ind = 0;

        //We obtain the next coordinates to process from the histogram looking for the values in hi, hj
        //with highest quality value
        while( (ind == 0) && (k > 0) )
        {
            int temp = VEC_ELEM(front,k);
            ind = temp;
            if (ind > 0)
            {
                i = dMij(hi,k,temp);
                j = dMij(hj,k,temp);

                if ( temp == VEC_ELEM(final,k))
                {
                    VEC_ELEM(front,k) = 0;
                    VEC_ELEM(final,k) = 0;
                }
                else
                    if ( temp == maxElem)
                        VEC_ELEM(front,k) = 1;
                    else
                        VEC_ELEM(front,k) += 1;
            }

            k -= 1;
        }

    }

}

void FringeProcessing::demodulate(MultidimArray<double> & im, double R, double S, double lambda, int size, MultidimArray<double> & phase, MultidimArray<double> & mod, int verbose)
{
    //Initial Setup
    MultidimArray< double > In, orMap, orModMap, dir, wphase;
    MultidimArray< std::complex<double> > sph;

    In.resizeNoCopy(im);
    orMap.resizeNoCopy(im);
    orModMap.resizeNoCopy(im);
    dir.resizeNoCopy(im);
    wphase.resizeNoCopy(im);
    sph.resizeNoCopy(im);

    //We obtain previous necessary maps
    //Normalized version of im and modulation map
    normalize(im,In,mod,R,S);

    Image<double> save;
    if (verbose == 1)
    {
	   save()=In;
	   save.write("PPP1.xmp");
    }

    STARTINGX(mod)=STARTINGY(mod)=0;
    //Orientation map of the fringes
    orMinDer(In, orMap, orModMap, size);
    if (verbose == 2)
    {
	   save()=orMap;
	   save.write("PPP2.xmp");
    }

    //We obtain the direction from the orientation map
    direction(orMap, orModMap, lambda, size, dir);
    if (verbose == 3)
    {
	   save()=dir;
	   save.write("PPP3.xmp");
    }
    //Spiral transform of the normalized image
    SPTH(In,sph);
    STARTINGX(sph)=STARTINGY(sph)=0;
    STARTINGX(In)=STARTINGY(In)=0;

    // We obtain the wrapped phase
    std::complex<double> ci = std::complex<double>(0,1.0);
    std::complex<double> temp;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(im)
    {
        temp = -ci*std::exp(ci*A2D_ELEM(dir,i,j))*A2D_ELEM(sph,i,j);
        A2D_ELEM(wphase,i,j) = std::atan2(temp.real(),A2D_ELEM(In,i,j));
    }

    if (verbose == 4)
    {
	   save()=wphase;
	   save.write("PPP3.xmp");
    }
    //Finally we obtain the unwrapped phase
    unwrapping(wphase, mod, lambda, size, phase);

    if (verbose == 5)
    {
	   save()=phase;
	   save.write("PPP4.xmp");
    }

    if (verbose > 5)
    {
	   save()=In;
	   save.write("PPP1.xmp");
	   save()=orMap;
	   save.write("PPP2.xmp");
	   save()=dir;
	   save.write("PPP3.xmp");
	   save()=wphase;
	   save.write("PPP4.xmp");
	   save()=phase;
	   save.write("PPP5.xmp");
	   save()=mod;
	   save.write("PPP6.xmp");

    }

}



