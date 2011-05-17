#include <data/image.h>
#include <data/filters.h>
#include <data/fftw.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class PolarTest : public ::testing::Test
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

TEST_F(PolarTest, rotate)
{
    MultidimArray<  double  > auxMul;
    MultidimArray<  double  > aux2Mul;
    aux2Mul.resize(3,3);
    DIRECT_A2D_ELEM(aux2Mul,0,0) = 0;
    DIRECT_A2D_ELEM(aux2Mul,0,1) = 2.1950049;
    DIRECT_A2D_ELEM(aux2Mul,0,2) = 0;

    DIRECT_A2D_ELEM(aux2Mul,1,0) = 2.6541736;
    DIRECT_A2D_ELEM(aux2Mul,1,1) = 2;
    DIRECT_A2D_ELEM(aux2Mul,1,2) = 1.3803737;

    DIRECT_A2D_ELEM(aux2Mul,2,0) = 0;
    DIRECT_A2D_ELEM(aux2Mul,2,1) = 3.9039731;
    DIRECT_A2D_ELEM(aux2Mul,2,2) = 0;

    rotate(BSPLINE3,auxMul, mulDouble,10,DONT_WRAP);
    EXPECT_EQ(aux2Mul,auxMul);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
