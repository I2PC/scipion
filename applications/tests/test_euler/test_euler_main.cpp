#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include <data/euler.h>
#include <data/geometry.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class EulerTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        if (chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str())==-1)
        	REPORT_ERROR(ERR_UNCLASSIFIED,"Could not change directory");
    }

    Matrix1D<double> origin,xaxis,yaxis,zaxis;
    double alpha, beta, gamma;
    Euler _euler;

};

TEST_F( EulerTest, eulerXmippMatch)
{
    double xp,yp,zp;
    _euler.setOrder(Euler::ZYZ);
    for (int _z = 0; _z < 360; _z += 30)
        for (int _y = 0; _y < 360; _y += 30)
            for (int _x = 0; _x < 360; _x += 30)
            {
                //note euler and xmipp are inverse
                Euler angles(DEG2RAD(_z),DEG2RAD( _y), DEG2RAD(_x), Euler::ZYZ);
                Matrix2D<double> m(4,4);
                Matrix2D<double> M(4,4);
                m.initIdentity();
                angles.toMatrix(m);
                Euler_angles2matrix(_x,_y,_z,M);

                FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
                {
                    EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
                }
                _euler.extract(m);
                _euler.toMatrix(m);
                Euler_matrix2angles(M,xp,yp,zp);
                Euler_angles2matrix(xp,yp,zp,M);
                FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
                {
                    EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
                }
            }
}

////////////////Ruler Rotate
TEST_F( EulerTest, eulerRotateX)
{
    Matrix2D<double> M(4,4);
    Matrix2D<double> m(4,4);
    M.initIdentity();
    m(0,0)=1;
    m(0,1)=0;
    m(0,2)=0;
    m(0,3)=0 ;
    m(1,0)=0;
    m(1,1)=0.87758255;
    m(1,2)=0.47942555;
    m(1,3)=0;
    m(2,0)=0;
    m(2,1)=-0.47942555;
    m(2,2)=0.87758255;
    m(2,3)=0 ;
    m(3,0)=0;
    m(3,1)=0;
    m(3,2)=0;
    m(3,3)=1 ;

    _euler.eulerRotate(M, vectorR3(.5,0.,0.));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
    {
        EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
    }
}

TEST_F( EulerTest, eulerRotateY)
{
    Matrix2D<double> M(4,4);
    Matrix2D<double> m(4,4);
    M.initIdentity();
    m(0,0)=0.96891242;
    m(0,1)=0;
    m(0,2)=-0.24740396;
    m(0,3)=0 ;
    m(1,0)=0;
    m(1,1)=1;
    m(1,2)=0.;
    m(1,3)=0;
    m(2,0)=0.24740396;
    m(2,1)=0.0;
    m(2,2)=0.96891242;
    m(2,3)=0 ;
    m(3,0)=0;
    m(3,1)=0;
    m(3,2)=0;
    m(3,3)=1. ;

    _euler.eulerRotate(M, vectorR3(0.,0.25,0.));

    FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
    {
        EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
    }
}

TEST_F( EulerTest, eulerRotateZ)
{
    Matrix2D<double> M(4,4);
    Matrix2D<double> m(4,4);
    M.initIdentity();
    m(0,0)=0.73168886;
    m(0,1)=0.68163878;
    m(0,2)=0;
    m(0,3)=0 ;
    m(1,0)=-0.68163878;
    m(1,1)=0.73168886;
    m(1,2)=0;
    m(1,3)=0;
    m(2,0)=0;
    m(2,1)=0;
    m(2,2)=1.;
    m(2,3)=0 ;
    m(3,0)=0;
    m(3,1)=0;
    m(3,2)=0;
    m(3,3)=1 ;

    _euler.eulerRotate(M, vectorR3(0.,0.0,0.75));

    FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
    {
        EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
    }
}

TEST_F( EulerTest, eulerRotateXYZ)
{
    Matrix2D<double> M(4,4);
    Matrix2D<double> m(4,4);
    M.initIdentity();
    m(0,0)=-2.310437e-01;
    m(0,1)=2.362753e-01;
    m(0,2)=9.438182e-01;
    m(0,3)=0 ;
    m(1,0)=-6.286172e-01;
    m(1,1)=-7.766573e-01;
    m(1,2)= 4.054479e-02;
    m(1,3)=0;
    m(2,0)=7.426031e-01;
    m(2,1)=-5.839327e-01;
    m(2,2)=3.279685e-01;
    m(2,3)=0 ;
    m(3,0)=0;
    m(3,1)=0;
    m(3,2)=0;
    m(3,3)=1. ;

    _euler.eulerRotate(M, vectorR3(0.123, -1.234, 2.345));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
    {
        EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
    }
}


/////////////////euler Angles

TEST_F( EulerTest, eulerAnglesXYZ)
{
    Matrix2D<double> M(4,4);
    M.initIdentity();
    Matrix2D<double> m(4,4);
    double _z = -3.05844 ;
    double _y =  -0.233197;
    double _x =  0.369401;

    Euler angles(_z, _y, _x, Euler::XYZ);
    angles.toMatrix(M);

    m(0,0)= 9.073022e-01;
    m(0,1)= 3.512840e-01;
    m(0,2)= 2.310892e-01;
    m(0,3)=0 ;
    m(1,0)= 3.777082e-01;
    m(1,1)=-9.223917e-01;
    m(1,2)=-8.080873e-02;
    m(1,3)=0;
    m(2,0)= 1.847679e-01;
    m(2,1)= 1.606022e-01;
    m(2,2)=-9.695709e-01;
    m(2,3)=0 ;
    m(3,0)= 0.000000e+00;
    m(3,1)= 0.000000e+00;
    m(3,2)= 0.000000e+00;
    m(3,3)=1 ;

    FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
    {
        EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
    }

}

TEST_F( EulerTest, eulerAnglesXZY)
{
    Matrix2D<double> M(4,4);
    M.initIdentity();
    Matrix2D<double> m(4,4);
    double _z = -3.05844 ;
    double _y =  -0.233197;
    double _x =  0.369401;

    Euler angles(_x, _y, _z, Euler::XYZ);
    angles.toMatrix(M);

    m(0,0)=-0.96957093;
    m(0,1)=-0.080808729;
    m(0,2)=0.23108916;
    m(0,3)=0.000000e+00;

    m(1,0)=0.16060221;
    m(1,1)=-0.92239171;
    m(1,2)=0.35128403;
    m(1,3)=0.000000e+00;

    m(2,0)=0.18476792;
    m(2,1)=0.3777082;
    m(2,2)=0.9073022;
    m(2,3)=0.000000e+00;

    m(3,0)=0.000000e+00;
    m(3,1)=0.000000e+00;
    m(3,2)=0.000000e+00;
    m(3,3)=1.000000e+00;

    FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
    {
        EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
    }
}


TEST_F( EulerTest, extract)
{
    Matrix2D<double> m(4,4);
    Matrix2D<double> M(4,4);
    Euler _euler;

    for (int _e = 0; _e < eulerOrderNumber; _e++)
    {
        Euler::eulerOrder order = eulerOrderList[_e];//Euler::XYX;
        //Euler::eulerOrder order = Euler::XYZ;
        _euler.init();
        _euler.setOrder(order);
        for (int _z = -360; _z < 360; _z += 30)
            for (int _y = -360; _y < 360; _y += 30)
                for (int _x = -360; _x < 360; _x += 30)
                {
                    //NOTE that x,y and z order should match the order "order" but since
                    //_z,_y and _x are never used but here I do not bother to order them
                    Euler angles(DEG2RAD(_x),DEG2RAD( _y), DEG2RAD(_z), order);
                    angles.toMatrix(m);
                    _euler.extract(m);
                    _euler.toMatrix(M);
                    //std::cerr << "DEBUG_ROB, m:" << m << std::endl;
                    //std::cerr << "DEBUG_ROB, M:" << M << std::endl;
                    FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
                    {
                        EXPECT_TRUE( fabs(M(i,j)-m(i,j))< XMIPP_EQUAL_ACCURACY);
                    }
                }
    }
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
