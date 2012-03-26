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
        chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
    }

    //Image to be fitted:
    MultidimArray<double> im;
    //File name of the image to process
    FileName tempImageName;

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

    fp.simulPattern(im,fp.SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    ASSERT_TRUE(XSIZE(im) == nx);
    ASSERT_TRUE(YSIZE(im) == ny);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,0) - 0.943527)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,1) - 0.975946)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,1,0) - 0.976113)<0.01);

    fp.simulPattern(im,fp.COMPLEX_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    ASSERT_TRUE(XSIZE(im) == nx);
    ASSERT_TRUE(YSIZE(im) == ny);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,0) + 0.972063)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,0,1) + 0.938343)<0.01);
    ASSERT_TRUE( std::abs(A2D_ELEM(im,1,0) + 0.986396)<0.01);

    //im.write(imageName);
}

//TEST_F( FringeProcessingTests, SPTH)
//{
//    tempImageName.initUniqueName("/tmp/temp_XXXXXX");
//
//	FringeProcessing fp;
//	FileName imageName2, imageName3;
//	imageName2 = "fringe/test2.txt";
//	imageName3 = "fringe/originalImag.txt";
//
//	MultidimArray< std::complex <double> > imProc;
//
//	MultidimArray< double > imProcReal;
//	MultidimArray< double > imProcImag;
//
//	MultidimArray<double> im;
//
//	int nx = 311;
//    int ny = 312;
//    double noiseLevel = 0.0;
//    double freq = 20;
//    Matrix1D<int> coefs(10);
//
//    fp.simulPattern(im,fp.SIMPLY_OPEN_FRINGES,nx,ny, noiseLevel,freq, coefs);
//    imProc.resizeNoCopy(im);
//
//    fp.SPTH(im, imProc);
//
//    imProc.getReal(imProcReal);
//    imProc.getImag(imProcImag);
//
//    imProcReal.write(tempImageName);
//    imProcImag.write(imageName2);
//    im.write(imageName3);
//
//    tempImageName.deleteFile();
//}

TEST_F( FringeProcessingTests, orMinDer)
{

	//FileName fpName;
	//fpName = baseName + "/applications/tests/test_fringe_processing/fp.txt";

	FringeProcessing fp;
	MultidimArray<double> im, orMap, orModMap;

	int nx = 311;
    int ny = 311;
    double noiseLevel = 0.0;
    double freq = 20;
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
	FileName fpName, Iname;
	fpName = "fringe/fp.txt";
	Iname = "fringe/IN.txt";

	FringeProcessing fp;
	MultidimArray<double> im, imN, imNMod;

	int nx = 311;
    int ny = 311;
    double noiseLevel = 0.0;
    double freq = 20;
    Matrix1D<int> coefs(10);

    fp.simulPattern(im,fp.SIMPLY_CLOSED_FRINGES,nx,ny, noiseLevel,freq, coefs);

    im.write(fpName);

    imN.resizeNoCopy(im);
    imNMod.resizeNoCopy(im);

    int fmin,fmax,num;
    fp.normalize(im,imN,imNMod, fmin, fmax, num);

    //imN.write(Iname);

	ASSERT_TRUE(true);

}


GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}




