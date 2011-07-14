#include <data/matrix2d.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class MatrixTest : public ::testing::Test
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
        A.resize(3,3);
        A(0,0) =-0.9234482 ;
        A(0,1) =  -0.38372311   ;
        A(0,2) =0 ;

        A(1,0) =0.38372311 ;
        A(1,1) =-0.9234482;
        A(1,2) =-0 ;

        A(2,0) =0;
        A(2,1) =0;
        A(2,2) =1 ;
    }
    Matrix2D<double> A;

};



TEST_F( MatrixTest, inverse)
{
    Matrix2D<double> B(3,3);
    Matrix2D<double> auxA(3,3);
    B.initIdentity();
    EXPECT_EQ(B.inv(),B) << "MatrixTest_inverse: identity matrixix failed";

    auxA(0,0) =-0.9234482 ;
    auxA(0,1) =  0.38372311   ;
    auxA(0,2) =0 ;

    auxA(1,0) =-0.38372311 ;
    auxA(1,1) =-0.9234482;
    auxA(1,2) =-0 ;

    auxA(2,0) =0;
    auxA(2,1) =0;
    auxA(2,2) =1 ;
    EXPECT_EQ(auxA,A.inv()) << "MatrixTest_inverse: rotation  matrix failed";


}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
