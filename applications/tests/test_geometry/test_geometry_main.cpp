#include <data/xmipp_image.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class GeometryTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        eulerMatrix.resize(3,3);
    }
    Matrix2D<double> eulerMatrix;
    double rot, tilt, psi;
    MultidimArray<double> vol;


    // virtual void TearDown() {}//Destructor

};

TEST_F(GeometryTest, angles2Matrix2Angles)
{
    double xrad, yrad, zrad;


    for (double _z = 0; _z < 360; _z += 30)
        for (double _y = 0; _y < 360; _y += 30)
            for (double _x = 0; _x < 360; _x += 30)
            {
                Euler_angles2matrix(_x, _y, _z, eulerMatrix);
                Euler_matrix2angles(eulerMatrix, rot, tilt, psi);

                rot = DEG2RAD(rot);
                tilt = DEG2RAD(tilt);
                psi = DEG2RAD(psi);
                xrad = DEG2RAD(_x);
                yrad = DEG2RAD(_y);
                zrad = DEG2RAD(_z);

                double r11 = cos(psi)*cos(tilt)*cos(rot)-sin(psi)*sin(rot);
                double r12 = cos(psi)*cos(tilt)*sin(rot)+sin(psi)*cos(rot)     ;
                double r13 = -cos(psi)*sin(tilt)                               ;
                double r22 = -sin(psi)*cos(tilt)*sin(rot)+cos(psi)*cos(rot)    ;
                double r23 = sin(psi)*sin(tilt)                                ;
                double r33 = cos(tilt)                                         ;

                double r011 = cos(zrad)*cos(yrad)*cos(xrad)-sin(zrad)*sin(xrad);
                double r012 = cos(zrad)*cos(yrad)*sin(xrad)+sin(zrad)*cos(xrad)     ;
                double r013 = -cos(zrad)*sin(yrad)                               ;
                double r022 = -sin(zrad)*cos(yrad)*sin(xrad)+cos(zrad)*cos(xrad)    ;
                double r023 = sin(zrad)*sin(yrad)                                ;
                double r033 = cos(yrad)                                         ;

                EXPECT_NEAR(r011, r11, XMIPP_EQUAL_ACCURACY);
                EXPECT_NEAR(r012, r12, XMIPP_EQUAL_ACCURACY);
                EXPECT_NEAR(r013, r13, XMIPP_EQUAL_ACCURACY);
                EXPECT_NEAR(r022, r22, XMIPP_EQUAL_ACCURACY);
                EXPECT_NEAR(r023, r23, XMIPP_EQUAL_ACCURACY);
                EXPECT_NEAR(r033, r33, XMIPP_EQUAL_ACCURACY);

            }
}
TEST_F(GeometryTest, rotateAngleAroundAxis)
{
    Matrix1D<double> axis(3);
    Matrix2D<double> matrix, rMatrix;


    double ang = 90;
    // Around X
    XX(axis) = 1;
    YY(axis) = 0;
    ZZ(axis) = 0;

    for (int sign = -1; sign < 2; sign+=2)
    {
        rotation3DMatrix(ang*sign, axis, rMatrix, false);

        EXPECT_NEAR(sign, MAT_ELEM(rMatrix,1,2), XMIPP_EQUAL_ACCURACY);
        EXPECT_NEAR(-sign, MAT_ELEM(rMatrix,2,1), XMIPP_EQUAL_ACCURACY);
    }

    // Around Y
    XX(axis) = 0;
    YY(axis) = 1;
    ZZ(axis) = 0;

    for (int sign = -1; sign < 2; sign+=2)
    {
        rotation3DMatrix(ang*sign, axis, rMatrix, false);

        EXPECT_NEAR(-sign, MAT_ELEM(rMatrix,0,2), XMIPP_EQUAL_ACCURACY);
        EXPECT_NEAR(sign, MAT_ELEM(rMatrix,2,0), XMIPP_EQUAL_ACCURACY);
    }

    // Around Z
    XX(axis) = 0;
    YY(axis) = 0;
    ZZ(axis) = 1;

    for (int sign = -1; sign < 2; sign+=2)
    {
        rotation3DMatrix(ang*sign, axis, rMatrix, false);

        EXPECT_NEAR(sign, MAT_ELEM(rMatrix,0,1), XMIPP_EQUAL_ACCURACY);
        EXPECT_NEAR(-sign, MAT_ELEM(rMatrix,1,0), XMIPP_EQUAL_ACCURACY);
    }


    //    EXPECT_NEAR(ang, tilt, XMIPP_EQUAL_ACCURACY);
    //    EXPECT_NEAR(0, rot, XMIPP_EQUAL_ACCURACY);
    //    EXPECT_NEAR(0, psi, XMIPP_EQUAL_ACCURACY);

    // Around Y negative
    //    ang = -90;
    //    rotation3DMatrix(ang, axis, rMatrix, false);
    //    Euler_matrix2angles(rMatrix, rot, tilt, psi);
    //
    //    rot = realWRAP(rot,-180,180);
    //    tilt = realWRAP(tilt,-180,180);
    //    psi = realWRAP(psi,-180,180);
    //
    //    if (ang * tilt < 0)
    //        Euler_another_set(rot, tilt, psi, rot, tilt, psi);
    //
    //    EXPECT_NEAR(ang, tilt, XMIPP_EQUAL_ACCURACY);
    //    EXPECT_NEAR(0, rot, XMIPP_EQUAL_ACCURACY);
    //    EXPECT_NEAR(0, psi, XMIPP_EQUAL_ACCURACY);

}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
