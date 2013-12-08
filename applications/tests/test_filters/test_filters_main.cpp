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
        if (chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str())==-1)
        	REPORT_ERROR(ERR_UNCLASSIFIED,"Could not change directory");

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
TEST_F( FiltersTest, alignImages)
{
    FileName fileTemp;
    fileTemp.initUniqueName("/tmp/temp_XXXXXX");
    fileTemp = fileTemp + ":spi";
    Image<double> I;
    I.read("filters/test2.spi");
    I().setXmippOrigin();

    // Transform the image: Mirror, shift and rotation
    Image<double> Itransformed, ItransformedMirror;
    Itransformed()=I();
    Matrix2D<double> A;
    rotation2DMatrix(15,A,true);
    MAT_ELEM(A,0,2)=-4;
    MAT_ELEM(A,1,2)= 6;
    selfApplyGeometry(BSPLINE3,Itransformed(),A,IS_NOT_INV,DONT_WRAP);
    Itransformed.write(fileTemp);

    ItransformedMirror()=Itransformed();
    ItransformedMirror().selfReverseX();
    ItransformedMirror().setXmippOrigin();

    // Align the images in 4 different ways
    Matrix2D<double> M1;
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    MultidimArray<double> Ialigned1=ItransformedMirror();
    alignImages(I(),Ialigned1,M1,DONT_WRAP,aux,aux2,aux3);

    MultidimArray<double> Ialigned2;
    applyGeometry(BSPLINE3, Ialigned2, ItransformedMirror(), M1, IS_NOT_INV, DONT_WRAP);

    Matrix2D<double> M2=M1;
    MAT_ELEM(M2,0,0)*=-1;
    MAT_ELEM(M2,1,0)*=-1;
    MultidimArray<double> Ialigned3;
    applyGeometry(BSPLINE3, Ialigned3, Itransformed(), M2, IS_NOT_INV, DONT_WRAP);

    bool flip;
    double scale, shiftX, shiftY, psi;
    transformationMatrix2Parameters2D(M2, flip, scale, shiftX, shiftY, psi);
    MDRow row;
    row.setValue(MDL_FLIP,flip);
    row.setValue(MDL_SCALE,scale);
    row.setValue(MDL_SHIFT_X,shiftX);
    row.setValue(MDL_SHIFT_Y,shiftY);
    row.setValue(MDL_ANGLE_PSI,psi);
    Image<double> Ialigned4;
    Ialigned4.readApplyGeo(fileTemp, row);

    MultidimArray<double> diff;

    diff=Ialigned1-Ialigned2;
    double mean=diff.computeAvg();
    EXPECT_NEAR(mean,0.0,1e-2);

    diff=Ialigned1-Ialigned3;
    mean=diff.computeAvg();
    EXPECT_NEAR(mean,0.0,1e-2);

    Ialigned4().setXmippOrigin();
    diff=Ialigned1-Ialigned4();
    mean=diff.computeAvg();
    EXPECT_NEAR(mean,0.0,1e-2);

    fileTemp.deleteFile();
}
GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
