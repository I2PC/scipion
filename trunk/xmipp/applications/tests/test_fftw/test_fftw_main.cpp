#include <data/image.h>
#include <data/filters.h>
#include <data/fftw.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class FftwTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
#define len 128
        //find binaries directory
        //        char szTmp[len];
        //        char pBuf[len];
        //        sprintf(szTmp, "/proc/%d/exe", getpid());
        //        int bytes = std::min(readlink(szTmp, pBuf, len), (ssize_t)len - 1);
        //        if(bytes >= 0)
        //            pBuf[bytes] = '\0';
        //imageName2 = filename + "/../applications/tests/test_fftw/singleImage2.spi";
        //myImage2.read(imageName2);123321445.xmp
    	mulDouble.resize(3,3);
    	DIRECT_A2D_ELEM(mulDouble,0,0) = 1;
        DIRECT_A2D_ELEM(mulDouble,0,1) = 2;
        DIRECT_A2D_ELEM(mulDouble,0,2) = 3;

        DIRECT_A2D_ELEM(mulDouble,1,0) = 3;
        DIRECT_A2D_ELEM(mulDouble,1,1) = 2;
        DIRECT_A2D_ELEM(mulDouble,1,2) = 1;

        DIRECT_A2D_ELEM(mulDouble,2,0) = 4;
        DIRECT_A2D_ELEM(mulDouble,2,1) = 4;
        DIRECT_A2D_ELEM(mulDouble,2,2) = 5;
}
    MultidimArray<  double  > mulDouble;

    // virtual void TearDown() {}//Destructor


};

TEST_F( FftwTest, directFourierTransform)
{

    MultidimArray< std::complex< double > > FFT1;
    FourierTransformer transformer1;
    transformer1.FourierTransform(mulDouble, FFT1, false);
    MultidimArray<std::complex<double> > auxFFT;
    auxFFT.resize(3,2);
    DIRECT_A2D_ELEM(auxFFT,0,0) = std::complex<double>(2.77778,0);
    DIRECT_A2D_ELEM(auxFFT,0,1) = std::complex<double>(-0.0555556,0.096225);
    DIRECT_A2D_ELEM(auxFFT,1,0) = std::complex<double>(-0.388889,0.673575) ;
    DIRECT_A2D_ELEM(auxFFT,1,1) = std::complex<double>(-0.388889,-0.096225);
    DIRECT_A2D_ELEM(auxFFT,2,0) = std::complex<double>(-0.388889,-0.673575) ;
    DIRECT_A2D_ELEM(auxFFT,2,1) = std::complex<double>(-0.0555556,0.288675) ;

    EXPECT_EQ(FFT1,auxFFT);
}


//correlation_matrix(I1, I2, Mcorr);
//
//TEST_F( FiltersTest, bestShift)
//{
// double x,y;
// bestShift(myImage1(),myImage1(),x,y);
// std::cerr << "x y " << x << " " << y << std::endl;
//
// EXPECT_DOUBLE_EQ(x,0.);
// EXPECT_DOUBLE_EQ(y,0.);
// bestShift(myImage1(),myImage2(),x,y);
// std::cerr << "x y " << x << " " << y << std::endl;
//}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
