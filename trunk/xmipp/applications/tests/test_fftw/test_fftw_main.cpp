#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/xmipp_fftw.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class FftwTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
#define len 128

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

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
