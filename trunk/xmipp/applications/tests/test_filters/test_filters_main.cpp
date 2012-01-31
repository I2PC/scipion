#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/xmipp_fftw.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class FiltersTest : public ::testing::Test
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
        mulDouble1.resize(3,3);
        DIRECT_A2D_ELEM(mulDouble1,0,0) = 1;
        DIRECT_A2D_ELEM(mulDouble1,0,1) = 2;
        DIRECT_A2D_ELEM(mulDouble1,0,2) = 3;

        DIRECT_A2D_ELEM(mulDouble1,1,0) = 3;
        DIRECT_A2D_ELEM(mulDouble1,1,1) = 2;
        DIRECT_A2D_ELEM(mulDouble1,1,2) = 1;

        DIRECT_A2D_ELEM(mulDouble1,2,0) = 4;
        DIRECT_A2D_ELEM(mulDouble1,2,1) = 4;
        DIRECT_A2D_ELEM(mulDouble1,2,2) = 5;

        mulDouble2.resize(4,4);
        DIRECT_A2D_ELEM(mulDouble2,0,0) = 1;
        DIRECT_A2D_ELEM(mulDouble2,0,1) = 2;
        DIRECT_A2D_ELEM(mulDouble2,0,2) = 3;
        DIRECT_A2D_ELEM(mulDouble2,0,3) = 3;

        DIRECT_A2D_ELEM(mulDouble2,0,0) = 1;
        DIRECT_A2D_ELEM(mulDouble2,1,1) = 2;
        DIRECT_A2D_ELEM(mulDouble2,2,2) = 3;
        DIRECT_A2D_ELEM(mulDouble2,3,3) = 3;

        DIRECT_A2D_ELEM(mulDouble2,0,0) = 1;
        DIRECT_A2D_ELEM(mulDouble2,1,1) = 2;
        DIRECT_A2D_ELEM(mulDouble2,2,2) = 3;
        DIRECT_A2D_ELEM(mulDouble2,3,3) = 3;

        DIRECT_A2D_ELEM(mulDouble2,0,0) = 1;
        DIRECT_A2D_ELEM(mulDouble2,1,1) = 2;
        DIRECT_A2D_ELEM(mulDouble2,2,2) = 3;
        DIRECT_A2D_ELEM(mulDouble2,3,3) = 3;

    }
    MultidimArray<  double  > mulDouble1;
    MultidimArray<  double  > mulDouble2;

    // virtual void TearDown() {}//Destructor

};

TEST_F( FiltersTest, bestShift)
{
 double x,y;
 MultidimArray<  double  > auxMul,auxMul2;
 CorrelationAux aux;
 auxMul = mulDouble1;
 auxMul.setXmippOrigin();
 bestShift(auxMul,auxMul,x,y,aux);
 EXPECT_DOUBLE_EQ(x,0.);
 EXPECT_DOUBLE_EQ(y,0.);

}

TEST_F( FiltersTest, correlation_matrix)
{
    MultidimArray<double> Mcorr;
    CorrelationAux aux;
    correlation_matrix(mulDouble1,mulDouble1,Mcorr,aux);
    MultidimArray<  double  > auxMul;
    auxMul.resize(3,3);
    DIRECT_A2D_ELEM(auxMul,0,0) = 64;
    DIRECT_A2D_ELEM(auxMul,0,1) = 62;
    DIRECT_A2D_ELEM(auxMul,0,2) = 66;

    DIRECT_A2D_ELEM(auxMul,1,0) = 78;
    DIRECT_A2D_ELEM(auxMul,1,1) = 85;
    DIRECT_A2D_ELEM(auxMul,1,2) = 78;

    DIRECT_A2D_ELEM(auxMul,2,0) = 66;
    DIRECT_A2D_ELEM(auxMul,2,1) = 62;
    DIRECT_A2D_ELEM(auxMul,2,2) = 64;

    EXPECT_EQ(auxMul,Mcorr);

}

TEST_F( FiltersTest, correlation)
{
 MultidimArray<  double  > auxMul,auxMul2;
 auxMul = mulDouble1;
 double result;
 auxMul.setXmippOrigin();
 result = correlationIndex(auxMul,auxMul);
 EXPECT_DOUBLE_EQ(result,1.);

}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
