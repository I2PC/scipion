#include <iostream>
#include <data/xmipp_fftw.h>
#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/transformations.h>

#include <gtest/gtest.h>
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


    FOR_ALL_ELEMENTS_IN_ARRAY2D(auxMul)
    {
        EXPECT_TRUE( fabs(DIRECT_A2D_ELEM(auxMul,0,0)-DIRECT_A2D_ELEM(mulDouble,0,0))< XMIPP_EQUAL_ACCURACY);
        EXPECT_TRUE( fabs(DIRECT_A2D_ELEM(auxMul,0,1)-DIRECT_A2D_ELEM(mulDouble,0,1))< XMIPP_EQUAL_ACCURACY);
        EXPECT_TRUE( fabs(DIRECT_A2D_ELEM(auxMul,1,0)-DIRECT_A2D_ELEM(mulDouble,1,0))< XMIPP_EQUAL_ACCURACY);
        EXPECT_TRUE( fabs(DIRECT_A2D_ELEM(auxMul,1,1)-DIRECT_A2D_ELEM(mulDouble,1,1))< XMIPP_EQUAL_ACCURACY);
    }
}

TEST_F(TransformationTest, geo2TransformationMatrix)
{
    //It would be nice to have a operator == for mdrow
    double scale = 2.;
    double rot = 40.9090909091;
    double tilt=81.8181818182;
    double psi=54.5454545455;
    double x =1.;
    double y=2.;
    double z=3;
    bool flip = true;
    double scale2 = -1000;
    double rot2 = -1000;
    double tilt2=-1000;
    double psi2=-1000;
    double x2 =-1000;
    double y2=-1000;
    double z2=-1000;
    bool flip2 = false;

    MDRow rowIn,rowOut;
    rowIn.setValue(MDL_SCALE,scale);
    rowIn.setValue(MDL_ANGLE_PSI,psi);
    rowIn.setValue(MDL_SHIFT_X,x);
    rowIn.setValue(MDL_SHIFT_Y,y);
    rowIn.setValue(MDL_FLIP, flip);


    Matrix2D<double> A(3, 3);

    geo2TransformationMatrix(rowIn,A,false);
    transformationMatrix2Geo(A,rowOut);
    rowOut.getValue(MDL_SCALE,scale2);
    rowOut.getValue(MDL_ANGLE_PSI,psi2);
    rowOut.getValue(MDL_SHIFT_X,x2);
    rowOut.getValue(MDL_SHIFT_Y,y2);
    rowOut.getValue(MDL_FLIP, flip2);

    EXPECT_DOUBLE_EQ(scale, scale2);
    EXPECT_DOUBLE_EQ(psi, psi2);
    EXPECT_DOUBLE_EQ(x, x2);
    EXPECT_DOUBLE_EQ(y, y2);
    EXPECT_EQ(flip, flip2);

    rowIn.setValue(MDL_ANGLE_ROT,rot);
    rowIn.setValue(MDL_ANGLE_TILT,tilt);
    rowIn.setValue(MDL_SHIFT_Z,z);

    A.initZeros(4, 4);

    geo2TransformationMatrix(rowIn,A,false);
    transformationMatrix2Geo(A,rowOut);
    rowOut.getValue(MDL_SCALE,scale2);
    rowOut.getValue(MDL_ANGLE_ROT,rot2);
    rowOut.getValue(MDL_ANGLE_TILT,tilt2);
    rowOut.getValue(MDL_ANGLE_PSI,psi2);
    rowOut.getValue(MDL_SHIFT_X,x2);
    rowOut.getValue(MDL_SHIFT_Y,y2);
    rowOut.getValue(MDL_SHIFT_Z,z2);
    rowOut.getValue(MDL_FLIP, flip2);

    EXPECT_DOUBLE_EQ(scale, scale2);
    EXPECT_DOUBLE_EQ(rot, rot2);
    EXPECT_DOUBLE_EQ(tilt, tilt2);
    EXPECT_DOUBLE_EQ(psi, psi2);
    EXPECT_DOUBLE_EQ(x, x2);
    EXPECT_DOUBLE_EQ(y, y2);
    EXPECT_DOUBLE_EQ(z, z2);
    EXPECT_EQ(flip, flip2);
}

TEST_F(TransformationTest, str2TransformationMatrix)
{
    //It would be nice to have a operator == for mdrow

    Matrix2D<double> M(4, 4);
    M.initIdentity();
    dMij(M, 0, 0) = -1.1601138;
    dMij(M, 0, 1) = -1.6291519;
    dMij(M, 0, 2) = 2;
    dMij(M, 1, 0) = -1.6291519;
    dMij(M, 1, 1) = 1.1601138;
    dMij(M, 1, 2) = 4;

    const char * matrixStr1 = " -1.1601138 -1.6291519 2 0 "
                              " -1.6291519  1.1601138 4 0 "
                              "  0          0         1 0 "
                              "  0          0         0 1 ";
    Matrix2D<double> matrix1;
    string2TransformationMatrix(matrixStr1, matrix1);
    EXPECT_EQ(M, matrix1);

    const char * matrixStr2 = " [[-1.1601138 -1.6291519 2 0], "
                              "  [-1.6291519  1.1601138 4 0], "
                              "  [ 0          0         1 0], "
                              "  [ 0          0         0 1]] ";
    Matrix2D<double> matrix2;
    string2TransformationMatrix(matrixStr2, matrix2);
    EXPECT_EQ(M, matrix2);

    Matrix2D<double> M3(3, 3);
    M.initIdentity();
    dMij(M, 0, 0) = -1.1601138;
    dMij(M, 0, 1) = -1.6291519;
    dMij(M, 0, 2) = 2;
    dMij(M, 1, 0) = -1.6291519;
    dMij(M, 1, 1) = 1.1601138;
    dMij(M, 1, 2) = 4;
    Matrix2D<double> matrix3;
    string2TransformationMatrix(matrixStr2, matrix3, 3);
    EXPECT_EQ(M, matrix2);


}

/*




DEBUG_JM: A:
 0.9420923 -1.3393493  1.1483055          2
-0.93493539 0.72492445  1.6125695          4
 1.4961143  1.2963904 0.28462967          6
         0          0          0          1
 */

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
