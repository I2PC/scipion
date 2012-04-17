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
#include <data/fringe_processing.h>
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
        chdir(((String)(getXmippPath() + (String)"/resources/test/fringe")).c_str());
    }

    //Image to be fitted:
    MultidimArray<double> im;
    //File name of the image to process

};

TEST_F( FringeProcessingTests, simulPattern)
{
    int nx = 311;
    int ny = 312;
    double noiseLevel = 0.0;
    double freq = 20;
    Matrix1D<int> coefs(10);

    coefs.initConstant(0);
    VEC_ELEM(coefs,1)=0;
    VEC_ELEM(coefs,6)=5;

    FringeProcessing fp;
    fp.simulPattern(im,fp.SIMPLY_OPEN_FRINGES,nx,ny, noiseLevel,freq, coefs);

    ASSERT_TRUE(XSIZE(im) == nx);
    ASSERT_TRUE(YSIZE(im) == ny);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,0) - 0.521457)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,1) - 0.626272)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,1,0) - 0.521457)<0.01);

    //im.write(imageName);
    freq = 1;
    fp.simulPattern(im,fp.SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    ASSERT_TRUE(XSIZE(im) == nx);
    ASSERT_TRUE(YSIZE(im) == ny);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,0) - 0.943527)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,1) - 0.975946)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,1,0) - 0.976113)<0.01);

    freq = 20;
    fp.simulPattern(im,fp.COMPLEX_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    ASSERT_TRUE(XSIZE(im) == nx);
    ASSERT_TRUE(YSIZE(im) == ny);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,0) + 0.972063)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,1) + 0.938343)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,1,0) + 0.986396)<0.01);

    //im.write(imageName);

}

TEST_F( FringeProcessingTests, SPTH)
{
    //FileName fpName, imProcRealName, imProcImagName;
    //fpName = "fp.txt";
    //imProcRealName  = "IpReal.txt";
    //imProcImagName  = "IpImag.txt";

    FringeProcessing fp;
    MultidimArray< std::complex <double> > imProc;
    MultidimArray< double > imProcReal;
    MultidimArray< double > imProcImag;
    MultidimArray<double> im;

    int nx = 311;
    int ny = 312;
    double noiseLevel = 0.0;
    double freq = 20;
    Matrix1D<int> coefs(10);

    fp.simulPattern(im,fp.SIMPLY_OPEN_FRINGES,nx,ny, noiseLevel,freq, coefs);
    imProc.resizeNoCopy(im);
    fp.SPTH(im, imProc);

    imProc.getReal(imProcReal);
    imProc.getImag(imProcImag);

    ASSERT_TRUE( (A2D_ELEM(imProcReal,10,10)  -  0)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcReal,10,20)  -  0)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcReal,20,10)  -  0)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcReal,20,20)  -  0)  < 1e-3);

    ASSERT_TRUE( (A2D_ELEM(imProcImag,10,10)  -  0.954154)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcImag,10,20)  -  0.536937)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcImag,20,10)  -  0.954154)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(imProcImag,20,20)  -  0.536937)  < 1e-3);

    //im.write(fpName);
    //imProcReal.write(imProcRealName);
    //imProcImag.write(imProcImagName);
}

TEST_F( FringeProcessingTests, orMinDer)
{

    //FileName fpName;
    //fpName = baseName + "/applications/tests/test_fringe_processing/fp.txt";

    FringeProcessing fp;
    MultidimArray<double> im, orMap, orModMap;

    int nx = 311;
    int ny = 311;
    double noiseLevel = 0.0;
    double freq = 1;
    Matrix1D<int> coefs(10);

    fp.simulPattern(im,fp.SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);
    orMap.resizeNoCopy(im);

    //im.write(fpName);

    int wSize = 2;
    fp.orMinDer(im,orMap,orModMap, wSize);

    ASSERT_TRUE(XSIZE(im) == XSIZE(orMap));
    ASSERT_TRUE(YSIZE(im) == YSIZE(orMap));
    ASSERT_TRUE(XSIZE(im) == XSIZE(orModMap));
    ASSERT_TRUE(YSIZE(im) == YSIZE(orModMap));

    ASSERT_TRUE( std::abs(A2D_ELEM(orMap,1,1) -  2.3562)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(orMap,1,5) -  2.3483)<0.01);

    ASSERT_TRUE( std::abs(A2D_ELEM(orModMap,1,1) -  0.0690)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(orModMap,1,5) -  0.3364)<0.01);

}

TEST_F( FringeProcessingTests, normalize)
{
    //FileName fpName, Iname, ModName;
    //fpName = "fp.txt";
    //Iname  = "IN.txt";
    //ModName= "Mod.txt";
    FringeProcessing fp;
    MultidimArray<double> im, In, Mod;

    int nx = 311;
    int ny = 311;
    double noiseLevel = 0.1;
    double freq = 1;
    Matrix1D<int> coefs(10);

    fp.simulPattern(im,fp.SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    //im.write(fpName);

    In.resizeNoCopy(im);
    Mod.resizeNoCopy(im);

    //aprox there are 5 fringe in the image
    double R = 5;
    double S = 10;

    fp.normalize(im,In,Mod, R, S);

    //We test some values comparing with the values recovered with Matlab
    ASSERT_TRUE( (A2D_ELEM(In,10,10)  -  0.924569)  < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,10,20)  -  0.924981)  < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,20,10)  -  0.924981)  < 1e-1);
    ASSERT_TRUE( (A2D_ELEM(In,20,20)  -  0.924981)  < 1e-1);

    //In.write(Iname);
    //Mod.write(ModName);
}

#define DEBUG
TEST_F( FringeProcessingTests, direction)
{
#ifdef DEBUG
    FileName ModName = "Mod.txt";
    FileName DirName = "Dir.txt";
    FileName OrName = "Or.txt";
#endif
    FringeProcessing fp;
    MultidimArray<double> im, orMap, orModMap, dirMap;

    int nx = 311;
    int ny = 311;
    int x = 150;
    int y = 150;
    double noiseLevel = 0.0;
    double freq = 1;
    Matrix1D<int> coefs(10);

    fp.simulPattern(im,fp.SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);
    orMap.resizeNoCopy(im);
    dirMap.resizeNoCopy(im);

    int wSize = 2;
    //We get the orientation map
    fp.orMinDer(im,orMap,orModMap, wSize);

    //we obtain the phase direction map from the orientation map
    double lambda = 1;
    int size = 3;
    fp.direction(orMap, orModMap, lambda, size, dirMap, x, y);

    //Comparing with Matlab results

    ASSERT_TRUE( (A2D_ELEM(dirMap,10,10)  - 2.35619)  < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(dirMap,10,20)  - 2.33116)  < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(dirMap,20,10)  - 2.38123)  < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(dirMap,20,20)  - 2.35619)  < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(dirMap,100,50) - 2.6420)  < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(dirMap,100,100) -  2.3536) < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(dirMap,200,100) -  -2.4388)< 1e-2);

#ifdef DEBUG
    orModMap.write(ModName);
    orMap.write(OrName);
    dirMap.write(DirName);
#endif
}
#undef DEBUG

#define DEBUG
TEST_F( FringeProcessingTests, unwrapping)
{
#ifdef DEBUG
    FileName PName = "P.txt";
    FileName uPName = "uP.txt";
#endif

    int nx = 311;
    int ny = 311;
    double noiseLevel = 0.0;

    FringeProcessing fp;
    MultidimArray<double> refPhase(nx,ny), wphase, comPhase, Q;
    wphase.resizeNoCopy(refPhase);
    Q.resizeNoCopy(refPhase);
    comPhase.resizeNoCopy(refPhase);

    double iMaxDim2 = 2./std::max(nx,ny);
    double value = 2*3.14159265;

    //we generate a gaussian phase map
    refPhase.setXmippOrigin();
    Q.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(refPhase)
    {
        A2D_ELEM(refPhase,i,j) = (50*std::exp(-0.5*(std::pow(i*iMaxDim2,2)+std::pow(j*iMaxDim2,2))))+rnd_gaus(0,noiseLevel);
    }
    //We wrap the phase inside [0 2pi]
    mod(refPhase,wphase,value);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(refPhase)
    {
        A2D_ELEM(Q,i,j) = (255*std::exp(-0.05*(std::pow(A2D_ELEM(wphase,i,j)-value/2,2)+std::pow(A2D_ELEM(wphase,i,j)-value/2,2))));
    }

    STARTINGX(Q)=STARTINGY(Q)=0;
    STARTINGX(wphase)=STARTINGY(wphase)=0;

    double lambda = 0.5;
    int size = 3;

    fp.unwrapping(wphase, Q, lambda, size, comPhase);

    //Comparing with Matlab results
    /*ASSERT_TRUE( (A2D_ELEM(comPhase,9,19)  -  -9.4841)   < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(comPhase,19,9)  -  -9.4679)   < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(comPhase,19,19) -  -8.2268)   < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(comPhase,99,49) -   5.6156)   < 1e-2);
    ASSERT_TRUE( (A2D_ELEM(comPhase,49,99) -   5.6156)   < 1e-2);
    */

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
	    FileName fpName = "fp.txt";
	    FileName phaseName = "Phase.txt";
#endif

	    FringeProcessing fp;
	    MultidimArray<double> im, mod, phase;

	    int nx = 311;
	    int ny = 311;
	    int x = 150;
	    int y = 150;
	    double noiseLevel = 0;
	    double freq = 2;
	    Matrix1D<int> coefs(10);

	    fp.simulPattern(im,fp.SIMPLY_CLOSED_FRINGES_MOD,nx,ny, noiseLevel,freq, coefs);
	    mod.resizeNoCopy(im);
	    phase.resizeNoCopy(im);

	    double lambda = 1;
	    int size = 3;
	    double R = 10;
	    double S = 5;
	    int verbose=6;

	    fp.demodulate(im,R,S,lambda,size, x, y, phase,mod,verbose);

	    //Comparing with Matlab results
	    /*ASSERT_TRUE( (A2D_ELEM(phase,30,30)  - 10.5929)  < 1e-2);
	    ASSERT_TRUE( (A2D_ELEM(phase,30,50)  - 7.22597)  < 1e-2);
	    ASSERT_TRUE( (A2D_ELEM(phase,50,30)  - 7.22597)  < 1e-2);
	    ASSERT_TRUE( (A2D_ELEM(phase,100,100)- 4.38535)  < 1e-2);
		*/

#ifdef DEBUG
	    im.write(fpName);
	    mod.write(ModName);
	    phase.write(phaseName);
#endif
}


GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}




