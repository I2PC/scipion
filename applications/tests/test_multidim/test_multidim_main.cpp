#include <data/multidim_array.h>
#include <data/matrix2d.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// Modify this test so it uses Fixures as test_image and test_metadata
TEST( MultidimTest, Size)
{
    MultidimArray<int> md;
    md.resize(2,3);
    size_t x,y,z, n;
    md.getDimensions(x,y,z,n);
    EXPECT_EQ((size_t)1, n) << "MultidimArray: wrong n size";
    EXPECT_EQ((size_t)1, z) << "MultidimArray: wrong y size";
    EXPECT_EQ((size_t)2, y) << "MultidimArray: wrong z size";
    EXPECT_EQ((size_t)3, x) << "MultidimArray: wrong x size";
}

TEST( MultidimTest, Assign)
{
    MultidimArray<double> * array = new MultidimArray<double>[5];
    MultidimArray<double> source(10);
    MultidimArray<double> source5(5);
    source5.initConstant(5.);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(source)
    source(i)=(double)i;
    array[0] = source;
    EXPECT_EQ(array[0], source);
    array[0] = source5;
    EXPECT_EQ(array[0], source5);
    delete[] array;
}

TEST( MultidimTest, Copy)
{
    MultidimArray<int> mdtarget,mdsource;
    mdsource.resize(2,3);
    DIRECT_MULTIDIM_ELEM(mdsource,0)=1;
    mdtarget=mdsource;
    DIRECT_MULTIDIM_ELEM(mdsource,0)=1;
    EXPECT_EQ(mdtarget,mdsource) << "MultidimArray: copy operator failed";
}

TEST( MultidimTest, CopyFromMatrix2D)
{
    MultidimArray<double> mdTarget;
    Matrix2D<double> mSource(2,2);
    mSource(0,0) = 1;
    mSource(0,1) = 2;
    mSource(1,0) = 3;
    mSource(1,1) = 4;

    mdTarget=mSource;

    EXPECT_EQ(dMn(mSource, 0),DIRECT_MULTIDIM_ELEM(mdTarget,0)) << "MultidimArray: copy from Matrix2D operator failed";
    EXPECT_EQ(dMn(mSource, 1),DIRECT_MULTIDIM_ELEM(mdTarget,1)) << "MultidimArray: copy from Matrix2D operator failed";
    EXPECT_EQ(dMn(mSource, 2),DIRECT_MULTIDIM_ELEM(mdTarget,2)) << "MultidimArray: copy from Matrix2D operator failed";
    EXPECT_EQ(dMn(mSource, 3),DIRECT_MULTIDIM_ELEM(mdTarget,3)) << "MultidimArray: copy from Matrix2D operator failed";
}

TEST( MultidimTest, typeCastComplex)
{
    MultidimArray<std::complex <double> > mdTarget;
    MultidimArray<double> mSource(2,2);

    mSource(0,0) = 1;
    mSource(0,1) = 2;
    mSource(1,0) = 3;
    mSource(1,1) = 4;

    typeCast(mSource,mdTarget);

    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 0),DIRECT_MULTIDIM_ELEM(mdTarget,0).real());
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 1),DIRECT_MULTIDIM_ELEM(mdTarget,1).real());
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 2),DIRECT_MULTIDIM_ELEM(mdTarget,2).real());
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 3),DIRECT_MULTIDIM_ELEM(mdTarget,3).real());

    EXPECT_EQ(0,DIRECT_MULTIDIM_ELEM(mdTarget,0).imag()) ;
    EXPECT_EQ(0,DIRECT_MULTIDIM_ELEM(mdTarget,1).imag()) ;
    EXPECT_EQ(0,DIRECT_MULTIDIM_ELEM(mdTarget,2).imag()) ;
    EXPECT_EQ(0,DIRECT_MULTIDIM_ELEM(mdTarget,3).imag()) ;

}

TEST( MultidimTest, getRealFromComplex)
{
    MultidimArray<std::complex <double> > mSource(2,2);
    MultidimArray<double> mdTarget;

    A2D_ELEM(mSource,0,0) = (std::complex<double>(0,0));
    A2D_ELEM(mSource,1,0) = (std::complex<double>(1,0));
    A2D_ELEM(mSource,0,1) = (std::complex<double>(2,0));
    A2D_ELEM(mSource,1,1) = (std::complex<double>(3,0));

    mSource.getReal(mdTarget);

    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 0).real(),DIRECT_MULTIDIM_ELEM(mdTarget,0));
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 1).real(),DIRECT_MULTIDIM_ELEM(mdTarget,1));
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 2).real(),DIRECT_MULTIDIM_ELEM(mdTarget,2));
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 3).real(),DIRECT_MULTIDIM_ELEM(mdTarget,3));
}

TEST( MultidimTest, getImagFromComplex)
{
    MultidimArray<std::complex <double> > mSource(2,2);
    MultidimArray<double> mdTarget;

    A2D_ELEM(mSource,0,0) = (std::complex<double>(0,0));
    A2D_ELEM(mSource,1,0) = (std::complex<double>(0,1));
    A2D_ELEM(mSource,0,1) = (std::complex<double>(0,2));
    A2D_ELEM(mSource,1,1) = (std::complex<double>(0,3));

    mSource.getImag(mdTarget);

    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 0).imag(),DIRECT_MULTIDIM_ELEM(mdTarget,0));
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 1).imag(),DIRECT_MULTIDIM_ELEM(mdTarget,1));
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 2).imag(),DIRECT_MULTIDIM_ELEM(mdTarget,2));
    EXPECT_EQ(DIRECT_MULTIDIM_ELEM(mSource, 3).imag(),DIRECT_MULTIDIM_ELEM(mdTarget,3));
}

TEST( MultidimTest, sincos)
{
    MultidimArray<double> mSource(2,2);
    MultidimArray<double> mdSIN;
    MultidimArray<double> mdCOS;

    A2D_ELEM(mSource,0,0) = 0;
    A2D_ELEM(mSource,1,0) = 3.14159265/2;
    A2D_ELEM(mSource,0,1) = 3.14159265;
    A2D_ELEM(mSource,1,1) = (3*3.14159265)/2;

    sincos(mSource,mdSIN,mdCOS);

    ASSERT_TRUE( (A2D_ELEM(mdSIN,0,0) -  0)<0.0001);
    ASSERT_TRUE( (A2D_ELEM(mdSIN,1,0) -  1)<0.0001);
    ASSERT_TRUE( (A2D_ELEM(mdSIN,0,1) -  0)<0.0001);
    ASSERT_TRUE( (A2D_ELEM(mdSIN,1,1) +  1)<0.0001);

    ASSERT_TRUE( (A2D_ELEM(mdCOS,0,0) -  1)<0.0001);
    ASSERT_TRUE( (A2D_ELEM(mdCOS,1,0) -  0)<0.0001);
    ASSERT_TRUE( (A2D_ELEM(mdCOS,0,1) +  1)<0.0001);
    ASSERT_TRUE( (A2D_ELEM(mdCOS,1,1) -  0)<0.0001);

}

TEST( MultidimTest, modulo)
{
    MultidimArray<double> mSource(3,3);
    MultidimArray<double> mTarget(3,3);

    A2D_ELEM(mSource,0,0) = 0;
    A2D_ELEM(mSource,0,1) = 10;
    A2D_ELEM(mSource,0,2) = 3.14159265;

    A2D_ELEM(mSource,1,0) = (20*3.14159265);
    A2D_ELEM(mSource,1,1) = (3.14159265/2);
    A2D_ELEM(mSource,1,2) = (3*3.14159265/2);

    A2D_ELEM(mSource,2,0) = (10*3.14159265/2);
    A2D_ELEM(mSource,2,1) = (50*3.14159265);
    A2D_ELEM(mSource,2,2) = (2*3.14159265);

    double value = 2*3.14159265;
    mod(mSource,mTarget,value);

    //We test the obtained values with the results obtained from method "mod" of Matlab
    ASSERT_TRUE( (A2D_ELEM(mTarget,0,0) - 0)       < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(mTarget,0,1) - 3.7168)  < 1e-3 );
    ASSERT_TRUE( (A2D_ELEM(mTarget,0,2) - 3.1416)  < 1e-3);

    ASSERT_TRUE( (A2D_ELEM(mTarget,1,0) - 0)       < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(mTarget,1,1) - 1.5708)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(mTarget,1,2) - 4.7124)  < 1e-3);

    ASSERT_TRUE( (A2D_ELEM(mTarget,2,0) - 3.1416)  < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(mTarget,2,1) - 0)       < 1e-3);
    ASSERT_TRUE( (A2D_ELEM(mTarget,2,2) - 0)       < 1e-3);

}


// check the reslicing is right
TEST( MultidimTest, getImage)
{
    XMIPP_TRY
    MultidimArray<float> imgTgt, imgRef;

    imgRef.resize(3,1,3,3);

    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(imgRef)
    dAi(imgRef, n) = n;


    imgRef.getImage(2, imgTgt);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(imgRef)
        {
            EXPECT_EQ(DIRECT_NZYX_ELEM(imgRef,2,k,i,j), DIRECT_ZYX_ELEM(imgTgt,k,i,j));
        }


    imgTgt.resize(6,1,3,3);

    imgRef.getImage(0, imgTgt, 5);
    imgRef.getImage(1, imgTgt, 3);
    imgRef.getImage(2, imgTgt, 1);


    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(imgRef)
    {
        EXPECT_EQ(DIRECT_NZYX_ELEM(imgRef,0,k,i,j), DIRECT_NZYX_ELEM(imgTgt,5,k,i,j));
        EXPECT_EQ(DIRECT_NZYX_ELEM(imgRef,1,k,i,j), DIRECT_NZYX_ELEM(imgTgt,3,k,i,j));
        EXPECT_EQ(DIRECT_NZYX_ELEM(imgRef,2,k,i,j), DIRECT_NZYX_ELEM(imgTgt,1,k,i,j));
    }

    XMIPP_CATCH
}

// check the reslicing is right
TEST( MultidimTest, reslice)
{
    XMIPP_TRY
    MultidimArray<float> imgSliced, imgRef;

    imgRef.resize(3,3,3);

    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(imgRef)
    dAi(imgRef, n) = n;

    imgRef.reslice(imgSliced, VIEW_Y_NEG);


    //    imgSliced.reslice(ImageGeneric::Y_NEG);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(imgRef)
    {
        EXPECT_EQ(DIRECT_ZYX_ELEM(imgRef,k,i,j), DIRECT_ZYX_ELEM(imgSliced,ZSIZE(imgSliced)-1-i,k,j));
    }

    imgRef.reslice(imgSliced, VIEW_X_NEG);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(imgRef)
    {
        EXPECT_EQ(DIRECT_ZYX_ELEM(imgRef,k,i,j), DIRECT_ZYX_ELEM(imgSliced,ZSIZE(imgSliced)-1-j,i,k));
    }
    XMIPP_CATCH
}



TEST( MultidimTest, mapFile)
{
    XMIPP_TRY
    MultidimArray<double> mda, mdaMap;

    mda.resize(32,32,32);
    mda.initRandom(0,1);

    mdaMap.setMmap(true);
    mdaMap.resize(32,32,32);

    mdaMap = mda;

    ASSERT_EQ(mda, mdaMap);
    XMIPP_CATCH
}

GTEST_API_ int main(int argc, char **argv)
{

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
