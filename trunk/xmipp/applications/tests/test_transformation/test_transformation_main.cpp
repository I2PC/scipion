#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/xmipp_fftw.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class TransformationTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
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

        vol.resize(3,3,3);

        DIRECT_A3D_ELEM(vol,0,1,0) = 1;
        DIRECT_A3D_ELEM(vol,0,1,1) = 1;
        DIRECT_A3D_ELEM(vol,0,1,2) = 1;

        DIRECT_A3D_ELEM(vol,0,2,0) = 2;
        DIRECT_A3D_ELEM(vol,0,2,1) = 2;
        DIRECT_A3D_ELEM(vol,0,2,2) = 2;

        DIRECT_A3D_ELEM(vol,0,0,0) = 3;
        DIRECT_A3D_ELEM(vol,0,0,1) = 3;
        DIRECT_A3D_ELEM(vol,0,0,2) = 3;
        //
        DIRECT_A3D_ELEM(vol,1,1,0) = 1;
        DIRECT_A3D_ELEM(vol,1,1,1) = 1;
        DIRECT_A3D_ELEM(vol,1,1,2) = 1;

        DIRECT_A3D_ELEM(vol,1,2,0) = 2;
        DIRECT_A3D_ELEM(vol,1,2,1) = 2;
        DIRECT_A3D_ELEM(vol,1,2,2) = 2;

        DIRECT_A3D_ELEM(vol,1,0,0) = 3;
        DIRECT_A3D_ELEM(vol,1,0,1) = 3;
        DIRECT_A3D_ELEM(vol,1,0,2) = 3;
        //
        DIRECT_A3D_ELEM(vol,2,1,0) = 1;
        DIRECT_A3D_ELEM(vol,2,1,1) = 1;
        DIRECT_A3D_ELEM(vol,2,1,2) = 1;

        DIRECT_A3D_ELEM(vol,2,2,0) = 2;
        DIRECT_A3D_ELEM(vol,2,2,1) = 2;
        DIRECT_A3D_ELEM(vol,2,2,2) = 2;

        DIRECT_A3D_ELEM(vol,2,0,0) = 3;
        DIRECT_A3D_ELEM(vol,2,0,1) = 3;
        DIRECT_A3D_ELEM(vol,2,0,2) = 3;
    }
    MultidimArray<  double  > mulDouble;
    MultidimArray<double> vol;


    // virtual void TearDown() {}//Destructor

};

TEST_F(TransformationTest, rotate)
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

TEST_F(TransformationTest, selfApplyGeometry)
{
    MultidimArray<double> volref;
    volref.resize(3,3,3);
    MultidimArray<double> volout;
    volout.resize(3,3,3);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(vol)
    DIRECT_A3D_ELEM(volref,k,intWRAP(i+1,0,2),j) = DIRECT_A3D_ELEM(vol,k,i,j);
    Matrix1D<double> center(3);
    XX(center)=0.;
    YY(center)=1.;
    ZZ(center)=0.;

    translate(BSPLINE3,volout,vol,center);
    EXPECT_EQ(volref,volout);
}

TEST_F(TransformationTest, scaleToSizeNearest)
{
    MultidimArray<double> imOut, auxMul;
    auxMul.resize(2,2);

    DIRECT_A2D_ELEM(auxMul,0,0) = DIRECT_A2D_ELEM(mulDouble,0,0);
    DIRECT_A2D_ELEM(auxMul,0,1) = DIRECT_A2D_ELEM(mulDouble,0,1);
    DIRECT_A2D_ELEM(auxMul,1,0) = DIRECT_A2D_ELEM(mulDouble,1,0);
    DIRECT_A2D_ELEM(auxMul,1,1) = DIRECT_A2D_ELEM(mulDouble,1,1);

    mulDouble.setXmippOrigin();
    scaleToSize(NEAREST, imOut, mulDouble, 2, 2);

    EXPECT_EQ(auxMul, imOut);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
