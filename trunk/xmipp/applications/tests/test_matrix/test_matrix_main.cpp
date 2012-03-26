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
        chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
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
    EXPECT_EQ(B.inv(),B) << "MatrixTest_inverse: identity matrix failed";

    auxA(0,0) =-0.9234482 ;
    auxA(0,1) =  0.38372311   ;
    auxA(0,2) =0 ;

    auxA(1,0) =-0.38372311 ;
    auxA(1,1) =-0.9234482;
    auxA(1,2) =-0 ;

    auxA(2,0) =0;
    auxA(2,1) =0;
    auxA(2,2) =1 ;
    EXPECT_EQ(auxA,A.inv()) << "MatrixTest_inverse: random rotation  matrix failed";

    auxA(0,0) =  1.;
    auxA(0,1) =  0.;
    auxA(0,2) = -0.;

    auxA(1,0) = -0.;
    auxA(1,1) =  1.;
    auxA(1,2) =  0 ;

    auxA(2,0) = 0;
    auxA(2,1) = 0;
    auxA(2,2) = 1;
    EXPECT_EQ(auxA,auxA.inv()) << "MatrixTest_inverse: identitity with negative 0 matrix failed";

}

TEST_F( MatrixTest, solveLinearSystem)
{
	 PseudoInverseHelper pseudoInverter;
	 Matrix2D<double> &A = pseudoInverter.A;
	 Matrix1D<double> &b = pseudoInverter.b;
	 A.resizeNoCopy(4,3);
	 b.resizeNoCopy(4);

	 A(0,0) = 1;
	 A(0,1) = 2;
	 A(0,2) = 3;
	 A(1,0) = 4;
	 A(1,1) = 5;
	 A(1,2) = 6;
	 A(2,0) = 7;
	 A(2,1) = 8;
	 A(2,2) = 9;
	 A(3,0) = 10;
	 A(3,1) = 11;
	 A(3,2) = 12;

	 b(0) = 14;
	 b(1) = 32;
	 b(2) = 50;
	 b(3) = 68;

    Matrix1D<double> auxX(3), X(3);

    auxX(0) = 1;
    auxX(1) = 2;
    auxX(2) = 3;

    solveLinearSystem(pseudoInverter, X);

    EXPECT_EQ(auxX,X) << "MatrixTest_solveLinearSystem failed";

}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
