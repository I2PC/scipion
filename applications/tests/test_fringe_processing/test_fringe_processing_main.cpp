/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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

#include <data/multidim_array.h>
#include <reconstruction/fringe_processing.h>
#include <data/xmipp_image.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"

// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class FringeProcessingTests : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
#define len 128

        //get example down1_42_Periodogramavg.psd
        if (!chdir(((String)(getXmippPath() + (String)"/resources/test/fringe")).c_str()))
        	REPORT_ERROR(ERR_UNCLASSIFIED,"Could not change directory");
    }

    //Image to be processed:
    MultidimArray<double> im;

};

TEST_F( FringeProcessingTests, simulPattern)
{
    size_t nx = 311;
    size_t ny = 312;
    double noiseLevel = 0.0;
    double freq = 20;
    Matrix1D<int> coefs(10);

    coefs.initConstant(0);
    VEC_ELEM(coefs,1)=0;
    VEC_ELEM(coefs,6)=5;

    simulPattern(im,SIMPLY_OPEN_FRINGES,nx,ny, noiseLevel,freq, coefs);

    ASSERT_TRUE(XSIZE(im) == nx);
    ASSERT_TRUE(YSIZE(im) == ny);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,0) - 0.521457)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,1) - 0.626272)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,1,0) - 0.521457)<0.01);

    //im.write(imageName);
    freq = 1;
    simulPattern(im,SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    ASSERT_TRUE(XSIZE(im) == nx);
    ASSERT_TRUE(YSIZE(im) == ny);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,0) - 0.943527)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,1) - 0.975946)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,1,0) - 0.976113)<0.01);

}

TEST_F( FringeProcessingTests, SPTH)
{
    //FileName fpName, imProcRealName, imProcImagName;
    //fpName = "txt";
    //imProcRealName  = "IpReal.txt";
    //imProcImagName  = "IpImag.txt";

    MultidimArray< std::complex <double> > imProc;
    MultidimArray< double > imProcReal;
    MultidimArray< double > imProcImag;
    MultidimArray<double> im;

    int nx = 311;
    int ny = 312;
    double noiseLevel = 0.0;
    double freq = 20;
    Matrix1D<int> coefs(10);

    simulPattern(im,SIMPLY_OPEN_FRINGES,nx,ny, noiseLevel,freq, coefs);
    imProc.resizeNoCopy(im);
    FourierTransformer ftrans(FFTW_BACKWARD);
    SPTH(ftrans, im, imProc);

    imProc.getReal(imProcReal);
    imProc.getImag(imProcImag);

    //im.write(fpName);
    //imProcReal.write(imProcRealName);
    //imProcImag.write(imProcImagName);

    ASSERT_TRUE( (A2D_ELEM(imProcReal,10,10)  -  0)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcReal,10,20)  -  0)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcReal,20,10)  -  0)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcReal,20,20)  -  0)  < 1e-3);

    ASSERT_TRUE( (A2D_ELEM(imProcImag,10,10)  -  0.954154)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcImag,10,20)  -  0.536937)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcImag,20,10)  -  0.954154)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcImag,20,20)  -  0.536937)  < 1e-3);

}


TEST_F( FringeProcessingTests, orMinDer)
{

#ifdef DEBUG
    FileName fpName, orName, orMapName;
    fpName   = "txt";
    orName   = "or.txt";
    orMapName= "orMap.txt";
#endif

    MultidimArray<double> im, orMap, orModMap;
    MultidimArray<bool> ROI;

    int nx = 311;
    int ny = 311;
    double noiseLevel = 0.0;
    double freq = 1;
    Matrix1D<int> coefs(10);

    simulPattern(im,SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);
    orMap.resizeNoCopy(im);
    ROI.resizeNoCopy(im);

    int rmin = 20;
    int rmax = 300;

    ROI.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( (temp > rmin) &&  (temp < rmax) )
            A2D_ELEM(ROI,i,j)= true;
        else
            A2D_ELEM(ROI,i,j)= false;
    }

    int wSize = 2;
    orMinDer(im,orMap,orModMap, wSize, ROI);

    ASSERT_TRUE(XSIZE(im) == XSIZE(orMap));
    ASSERT_TRUE(YSIZE(im) == YSIZE(orMap));
    ASSERT_TRUE(XSIZE(im) == XSIZE(orModMap));
    ASSERT_TRUE(YSIZE(im) == YSIZE(orModMap));
    ASSERT_TRUE( std::abs(A2D_ELEM(orMap,1,1) -  2.3562)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(orMap,1,5) -  2.3483)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(orModMap,1,1) -  0.0690)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(orModMap,1,5) -  0.3364)<0.01);


#ifdef DEBUG

    im.write(fpName);
    orMap.write(orName);
    orModMap.write(orMapName);

#endif

}

TEST_F( FringeProcessingTests, normalize)
{

#ifdef DEBUG
    FileName fpName, Iname, ModName;
    fpName = "txt";
    Iname  = "IN.txt";
    ModName= "Mod.txt";
#endif

    MultidimArray<double> im, In, Mod;
    MultidimArray<bool> ROI;

    int nx = 311;
    int ny = 311;
    double noiseLevel = 0.1;
    double freq = 1;
    Matrix1D<int> coefs(10);

    simulPattern(im,SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    In.resizeNoCopy(im);
    Mod.resizeNoCopy(im);
    ROI.resizeNoCopy(im);

    int rmin = 20;
    int rmax = 150;
    ROI.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( (temp > rmin) &&  (temp < rmax) )
            A2D_ELEM(ROI,i,j)= true;
        else
            A2D_ELEM(ROI,i,j)= false;
    }


    //aprox there are 5 fringe in the image
    double R = 5;
    double S = 10;
    FourierTransformer ftrans(FFTW_BACKWARD);
    normalize(ftrans,im,In,Mod, R, S, ROI);

    //We test some values comparing with the values recovered with Matlab
    ASSERT_TRUE( (A2D_ELEM(In,10,10)  -  0.924569)  < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,10,20)  -  0.924981)  < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,20,10)  -  0.924981)  < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,20,20)  -  0.924981)  < 1e-1);

#ifdef DEBUG

    im.write(fpName);
    In.write(Iname);
    Mod.write(ModName);
#endif

}
#undef DEBUG


TEST_F( FringeProcessingTests, normalizeWB)
{

#ifdef DEBUG
    FileName fpName, Iname, ModName;
    fpName = "I.txt";
    Iname  = "IN.txt";
    ModName= "Mod.txt";
#endif

    MultidimArray<double> im, In, Mod;
    MultidimArray<bool> ROI;

    int nx = 311;
    int ny = 311;
    double noiseLevel = 0.1;
    double freq = 1;
    Matrix1D<int> coefs(10);

    simulPattern(im,SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    In.resizeNoCopy(im);
    Mod.resizeNoCopy(im);
    ROI.resizeNoCopy(im);

    int rmin = 20;
    int rmax = 150;
    ROI.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( (temp > rmin) &&  (temp < rmax) )
            A2D_ELEM(ROI,i,j)= true;
        else
            A2D_ELEM(ROI,i,j)= false;
    }

    normalizeWB(im,In,Mod, rmin, rmax, ROI);

    //We test some values comparing with the values recovered with Matlab
    ASSERT_TRUE( (A2D_ELEM(In,100,100)  -  0.99979)    < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,100,200)  -  0.553923)   < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,200,100)  -  0.55329)    < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,200,200)  - -0.368664)   < 1e-1);

#ifdef DEBUG

    im.write(fpName);
    In.write(Iname);
    Mod.write(ModName);
#endif

}
#undef DEBUG


TEST_F( FringeProcessingTests, unwrapping)
{

#ifdef DEBUG
    FileName PName = "P.txt";
    FileName uPName = "uP.txt";
#endif

    int nx = 311;
    int ny = 311;
    double noiseLevel = 0;

    MultidimArray<double> refPhase(nx,ny), wphase, comPhase, im, orMap, orModMap;
    MultidimArray<bool> ROI;

    wphase.resizeNoCopy(refPhase);
    comPhase.resizeNoCopy(refPhase);
    ROI.resizeNoCopy(refPhase);
    im.resizeNoCopy(refPhase);
    orMap.resizeNoCopy(refPhase);
    orModMap.resizeNoCopy(refPhase);

    double iMaxDim2 = 2./std::max(nx,ny);
    double value = 2*3.14159265;

    //we generate a gaussian phase map
    refPhase.setXmippOrigin();
    im.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(refPhase)
    {
        A2D_ELEM(refPhase,i,j) = (50*std::exp(-0.5*(std::pow(i*iMaxDim2,2)+std::pow(j*iMaxDim2,2))))+rnd_gaus(0,noiseLevel);
        A2D_ELEM(im,i,j) = std::cos(A2D_ELEM(refPhase,i,j));
    }
    //We wrap the phase inside [0 2pi]
    mod(refPhase,wphase,value);

    STARTINGX(wphase)=STARTINGY(wphase)=0;
    STARTINGX(im)=STARTINGY(im)=0;

    int rmin = 2;
    int rmax = 300;

    ROI.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( (temp > rmin) &&  (temp < rmax) )
            A2D_ELEM(ROI,i,j)= true;
        else
            A2D_ELEM(ROI,i,j)= false;
    }

    int wSize = 2;
    orMinDer(im, orMap,orModMap, wSize, ROI);

    double lambda = 0.4;
    int size = 4;
    unwrapping(wphase, orModMap, lambda, size, comPhase);

    //Comparing with Matlab results
    ASSERT_TRUE( (A2D_ELEM(comPhase,9,19)  -  -2.87779)   < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(comPhase,19,9)  -  -2.90942)   < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(comPhase,19,19) -  -1.66237)   < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(comPhase,99,49) -   11.7772)   < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(comPhase,49,99) -   11.7503)   < 1e-1);


#ifdef DEBUG

    comPhase.write(uPName);
    wphase.write(PName);
#endif

}
#undef DEBUG

TEST_F( FringeProcessingTests, demodulate)
{

#ifdef DEBUG
    FileName ModName = "Mod.txt";
    FileName fpName = "txt";
    FileName phaseName = "Phase.txt";
#endif

    MultidimArray<double> im, mod, phase;

    int nx = 311;
    int ny = 311;
    int x = 150;
    int y = 65;
    double noiseLevel = 0;
    double freq = 2;
    Matrix1D<double> coefs(13);

    simulPattern(im,SIMPLY_CLOSED_FRINGES_MOD,nx,ny, noiseLevel,freq, coefs);
    mod.resizeNoCopy(im);
    phase.resizeNoCopy(im);

    double lambda = 1;
    int size = 3,rmin=40, rmax=150;
    int verbose=0;

    coefs.initConstant(0);
    VEC_ELEM(coefs,0) = 1;
    VEC_ELEM(coefs,1) = 1;
    VEC_ELEM(coefs,2) = 1;
    VEC_ELEM(coefs,3) = 1;
    VEC_ELEM(coefs,4) = 1;
    VEC_ELEM(coefs,5) = 1;
    VEC_ELEM(coefs,7) = 1;
    VEC_ELEM(coefs,8) = 1;
    VEC_ELEM(coefs,12) = 1;

    demodulate(im, lambda,size, x, y, rmin, rmax,phase,mod, coefs, verbose);

    //Comparing with Matlab results
    ASSERT_TRUE( (A2D_ELEM(phase,30,30)  - 10.5929)  < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(phase,30,50)  - 7.22597)  < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(phase,50,30)  - 7.22597)  < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(phase,100,100)- 34.3254)  < 1e-2);


#ifdef DEBUG

    im.write(fpName);
    mod.write(ModName);
    phase.write(phaseName);

#endif
}

TEST_F( FringeProcessingTests, firsPSDZero)
{
#ifdef DEBUG
    FileName PName = "Psd.txt";
    FileName PName2 = "pX.txt";
    FileName PName3 = "pY.txt";
#endif

    int nx = 311;
    int ny = 311;

    MultidimArray<double> im(nx,ny);
    Matrix2D<double> envM;
    envM.initGaussian(nx,100);

    double noiseLevel = 0;
    double freq = 1;
    Matrix1D<double> coefs(13);

    simulPattern(im,SIMPLY_CLOSED_FRINGES_MOD,nx,ny, noiseLevel,freq, coefs);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        A2D_ELEM(im,i,j) = A2D_ELEM(im,i,j)* A2D_ELEM(im,i,j)*dMij(envM,i,j);
    }

    int numPts = 100;
    Matrix1D<double> ptsX(numPts), ptsY(numPts);
    //firsPSDZero(im,ptsX,ptsY,0.05*XSIZE(im),0.8*XSIZE(im),numPts,1);

    //////////////////////
    //Comparing the results
    //////////////////////
    //In development not working yet!!!!!!!!!!
    /*ASSERT_TRUE( (VEC_ELEM(ptsX,0)  - 29.93375)   < 1e-2);
    ASSERT_TRUE( (VEC_ELEM(ptsX,1)  - 29.874683)  < 1e-2);
    ASSERT_TRUE( (VEC_ELEM(ptsX,2)  - 28.540659)  < 1e-2);
    ASSERT_TRUE( (VEC_ELEM(ptsX,3)  - 28.257948)  < 1e-2);

    ASSERT_TRUE( (VEC_ELEM(ptsY,0)  - 0)          < 1e-2);
    ASSERT_TRUE( (VEC_ELEM(ptsY,1)  - 1.8795557)  < 1e-2);
    ASSERT_TRUE( (VEC_ELEM(ptsY,2)  - 3.6055238)  < 1e-2);
    ASSERT_TRUE( (VEC_ELEM(ptsY,3)  - 5.390492)   < 1e-2);
    */


#ifdef DEBUG

    im.write(PName);
    ptsX.write(PName2);
    ptsY.write(PName3);

#endif
}

/*TEST_F( FringeProcessingTests, fitEllipse)
{

#ifdef DEBUG
    FileName PName2 = "pX.txt";
    FileName PName3 = "pY.txt";

    FileName PName4 = "pXA.txt";
    FileName PName5 = "pYA.txt";

#endif

    int numPoints = 100;
    double angle = 0;
    double angleStep = 2*PI/numPoints;
    double radiusX = 20;
    double radiusY = 50;

    Matrix1D<double> ptsX(numPoints), ptsY(numPoints);
    for (int var = 0; var < numPoints; var++)
    {
        VEC_ELEM(ptsX,var) = radiusX*cos(angle)+radiusY*sin(angle)+10;
        VEC_ELEM(ptsY,var) = -radiusY*sin(angle)+radiusX*cos(angle)+25;
        angle += angleStep;

    }

    double x0,y0,majorAxis, minorAxis, ellipseAngle;
    Matrix1D<double> ptsX2 = ptsX;
    Matrix1D<double> ptsY2 = ptsY;
    fitEllipse(ptsX2, ptsY2, x0, y0, majorAxis, minorAxis, ellipseAngle);

    ASSERT_TRUE( (x0  - 20)   < 1e-2);
    ASSERT_TRUE( (y0  - 50)  < 1e-2);
    ASSERT_TRUE( (majorAxis - 28.2843)  < 1e-2);
    ASSERT_TRUE( (minorAxis  - 70.717)  < 1e-2);
    ASSERT_TRUE( (ellipseAngle  - -0.785398)  < 1e-2);


#ifdef DEBUG

    ptsX.write(PName2);
    ptsY.write(PName3);

    ptsX2.write(PName4);
    ptsY2.write(PName5);


#endif

}*/

/*
TEST_F( FringeProcessingTests, testVahid)
{

    mkdir( ((String)(getXmippPath() + (String)"/resources/test/fringe")).c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
    chdir(((String)(getXmippPath() + (String)"/resources/test/fringe")).c_str());

    FileName fpName, Iname, ModName;
    fpName = "I.txt";
    Iname  = "IN.txt";
    ModName= "Mod.txt";

    MultidimArray<double> im, In, Mod;
    MultidimArray<bool> ROI;

    Image<double> img;
    img.read("KLH_Dataset_I_Test_0001.mrc");

    im = img();
    In.resizeNoCopy(im);
    Mod.resizeNoCopy(im);
    ROI.resizeNoCopy(im);


    int rmin = 30;
    int rmax = 60;

    ROI.setXmippOrigin();
    ROI.initConstant(true);

    double R = 50;
    double S = 15;
    normalize(im,In,Mod, R, S, ROI);

    //aprox there are 5 fringe in the image
    //double R = 1;
    //double S = 20;
    //normalize(im,In,Mod, R, S, ROI);

    In.write(fpName);

}
*/

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}




