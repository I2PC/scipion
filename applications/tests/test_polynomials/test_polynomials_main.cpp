#include <data/xmipp_polynomials.h>
#include <data/xmipp_image.h>
#include <data/matrix1d.h>
#include <data/xmipp_fft.h>
#include <data/matrix2d.h>
#include <data/multidim_array.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class PolynomialsTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        //get example down1_42_Periodogramavg.psd
        chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
        imageName =  "polynomials/down1_42_Periodogramavg.psd";
        try
        {
            im.read(imageName);
        }
        catch (XmippError &xe)
        {
            std::cerr << xe;
            exit(-1);
        }
    }
    //Image to be fitted:
    Image<double> im;
    //File name of the image to process
    FileName imageName;

};

TEST_F( PolynomialsTest, ZernikeFit)
{
    // coefficients obtained using Matlab
    Matrix1D<double> matlabCoeffs;
    // coefficients obtained using Matlab
    Matrix1D<double> xmippCoeffs;

    //xmippCoeffs.resizeNoCopy(10);
    xmippCoeffs.resizeNoCopy(10);
    xmippCoeffs.initConstant(1);
    matlabCoeffs.resizeNoCopy(10);

    VEC_ELEM(matlabCoeffs,0) =   1.24429e-18;
    VEC_ELEM(matlabCoeffs,1) =  -1.36005e-19;
    VEC_ELEM(matlabCoeffs,2) =   9.16062e-19;
    VEC_ELEM(matlabCoeffs,3) =  -5.33021e-19;
    VEC_ELEM(matlabCoeffs,4) =   5.32969e-19;
    VEC_ELEM(matlabCoeffs,5) =  -1.58317e-19;
    VEC_ELEM(matlabCoeffs,6) =   7.30143e-11;
    VEC_ELEM(matlabCoeffs,7) =  -3.39904e-20;
    VEC_ELEM(matlabCoeffs,8) =   2.28636e-19;
    VEC_ELEM(matlabCoeffs,9) =   -8.19429e-10;

    int rmin = -1;
    int rmax = XSIZE(MULTIDIM_ARRAY(im))/2;

    MultidimArray< bool > ROI;
    ROI.resizeNoCopy(MULTIDIM_ARRAY(im));
    ROI.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( (temp > rmin) &&  (temp < rmax) )
            A2D_ELEM(ROI,i,j)= true;
        else
            A2D_ELEM(ROI,i,j)= false;
    }

    CenterFFT(MULTIDIM_ARRAY(im), true);
    PolyZernikes polynom;
    Matrix1D<int> coefs(10);
    coefs.initConstant(1);
    MultidimArray<double> weight;
    weight.resizeNoCopy(im());
    weight.initConstant(1.0);

    polynom.fit(coefs,MULTIDIM_ARRAY(im),weight,ROI,0);
    xmippCoeffs = COEFFICIENTS(polynom);

    Matrix1D<double> error = xmippCoeffs - matlabCoeffs;

    for(int i =0; i<  VEC_XSIZE(error); i++)
    {
        ASSERT_TRUE(std::abs(VEC_ELEM(error,i))<0.01) << "Zernike fit: no correspondence between matlab and xmipp zernike coefficients";
    }
}

TEST_F( PolynomialsTest, ZernikePols)
{
    CenterFFT(MULTIDIM_ARRAY(im), true);
    Matrix1D<int> coefs(8);
    coefs.initConstant(0);

    VEC_ELEM(coefs,3)=1;


    PolyZernikes polynom;

    int rmin = 100;
    int rmax = 1000;

    MultidimArray< bool > ROI;
    ROI.resizeNoCopy(MULTIDIM_ARRAY(im));
    ROI.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( (temp > rmin) &&  (temp < rmax) )
            A2D_ELEM(ROI,i,j)= true;
        else
            A2D_ELEM(ROI,i,j)= false;
    }

    im().initZeros();
    polynom.zernikePols(coefs,MULTIDIM_ARRAY(im),ROI);

    ASSERT_TRUE(std::abs(A2D_ELEM(MULTIDIM_ARRAY(im),0,0)+       0)<0.01) << "Zernike Pols: no correspondence between matlab and xmipp zernike coefficients";
    ASSERT_TRUE(std::abs(A2D_ELEM(MULTIDIM_ARRAY(im),0,1)+0.00779724)<0.01) << "Zernike Pols: no correspondence between matlab and xmipp zernike coefficients";
    ASSERT_TRUE(std::abs(A2D_ELEM(MULTIDIM_ARRAY(im),1,0)-0.00779724)<0.01) << "Zernike Pols: no correspondence between matlab and xmipp zernike coefficients";
    ASSERT_TRUE(std::abs(A2D_ELEM(MULTIDIM_ARRAY(im),250,10)-0.922852)<0.01) << "Zernike Pols: no correspondence between matlab and xmipp zernike coefficients";
    ASSERT_TRUE(std::abs(A2D_ELEM(MULTIDIM_ARRAY(im),10, 250)+0.922852)<0.01) << "Zernike Pols: no correspondence between matlab and xmipp zernike coefficients";

}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
