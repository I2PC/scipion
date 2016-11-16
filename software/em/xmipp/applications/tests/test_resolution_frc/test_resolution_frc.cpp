
#include <stdlib.h>
#include <data/xmipp_image.h>
#include <data/xmipp_image_extension.h>
#include <iostream>
#include <gtest/gtest.h>
#include <data/xmipp_fftw.h>

// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class ResolutionFSCTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        XMIPP_TRY
        ;
        XMIPP_CATCH
    }


};


TEST_F( ResolutionFSCTest, copy)
{

    Image<double> vol1, vol2;
    vol1().initZeros(3, 3, 3);
    vol2().initZeros(3, 3, 3);
    
    DIRECT_A3D_ELEM(vol1(),0,0,0) = 1;
    DIRECT_A3D_ELEM(vol1(),0,0,1) = 2;
    DIRECT_A3D_ELEM(vol1(),0,0,2) = 3;

    DIRECT_A3D_ELEM(vol1(),0,1,0) = 4;
    DIRECT_A3D_ELEM(vol1(),0,1,1) = 5;
    DIRECT_A3D_ELEM(vol1(),0,1,2) = 6;

    DIRECT_A3D_ELEM(vol1(),0,2,0) = 7;
    DIRECT_A3D_ELEM(vol1(),0,2,1) = 8;
    DIRECT_A3D_ELEM(vol1(),0,2,2) = 9;

    DIRECT_A3D_ELEM(vol1(),1,0,0) = 10;
    DIRECT_A3D_ELEM(vol1(),1,0,1) = 11;
    DIRECT_A3D_ELEM(vol1(),1,0,2) = 12;

    DIRECT_A3D_ELEM(vol1(),1,1,0) = 13;
    DIRECT_A3D_ELEM(vol1(),1,1,1) = 14;
    DIRECT_A3D_ELEM(vol1(),1,1,2) = 15;

    DIRECT_A3D_ELEM(vol1(),1,2,0) = 17;
    DIRECT_A3D_ELEM(vol1(),1,2,1) = 18;
    DIRECT_A3D_ELEM(vol1(),1,2,2) = 19;

    DIRECT_A3D_ELEM(vol1(),2,0,0) = 20;
    DIRECT_A3D_ELEM(vol1(),2,0,1) = 21;
    DIRECT_A3D_ELEM(vol1(),2,0,2) = 22;

    DIRECT_A3D_ELEM(vol1(),2,1,0) = 23;
    DIRECT_A3D_ELEM(vol1(),2,1,1) = 24;
    DIRECT_A3D_ELEM(vol1(),2,1,2) = 25;

    DIRECT_A3D_ELEM(vol1(),2,2,0) = 26.4;
    DIRECT_A3D_ELEM(vol1(),2,2,1) = 27.5;
    DIRECT_A3D_ELEM(vol1(),2,2,2) = 28.5;
    
    //vol2
    DIRECT_A3D_ELEM(vol2(),2,0,0) = 1.5;
    DIRECT_A3D_ELEM(vol2(),2,0,1) = 2.4;
    DIRECT_A3D_ELEM(vol2(),2,0,2) = 3.3;

    DIRECT_A3D_ELEM(vol2(),2,1,0) = 4.6;
    DIRECT_A3D_ELEM(vol2(),2,1,1) = 5.7;
    DIRECT_A3D_ELEM(vol2(),2,1,2) = 6.4;

    DIRECT_A3D_ELEM(vol2(),2,2,0) = 7.3;
    DIRECT_A3D_ELEM(vol2(),2,2,1) = 8.2;
    DIRECT_A3D_ELEM(vol2(),2,2,2) = 9.5;

    DIRECT_A3D_ELEM(vol2(),1,0,0) = 10.2;
    DIRECT_A3D_ELEM(vol2(),1,0,1) = 11.4;
    DIRECT_A3D_ELEM(vol2(),1,0,2) = 12.5;

    DIRECT_A3D_ELEM(vol2(),1,1,0) = 13.6;
    DIRECT_A3D_ELEM(vol2(),1,1,1) = 14.5;
    DIRECT_A3D_ELEM(vol2(),1,1,2) = 15.7;

    DIRECT_A3D_ELEM(vol2(),1,2,0) = 17.3;
    DIRECT_A3D_ELEM(vol2(),1,2,1) = 18.2;
    DIRECT_A3D_ELEM(vol2(),1,2,2) = 19.4;

    DIRECT_A3D_ELEM(vol2(),0,0,0) = 20.3;
    DIRECT_A3D_ELEM(vol2(),0,0,1) = 21.4;
    DIRECT_A3D_ELEM(vol2(),0,0,2) = 22.5;

    DIRECT_A3D_ELEM(vol2(),0,1,0) = 23.4;
    DIRECT_A3D_ELEM(vol2(),0,1,1) = 24.5;
    DIRECT_A3D_ELEM(vol2(),0,1,2) = 25.6;

    DIRECT_A3D_ELEM(vol2(),0,2,0) = 26.7;
    DIRECT_A3D_ELEM(vol2(),0,2,1) = 24;
    DIRECT_A3D_ELEM(vol2(),0,2,2) = 23;
/*

    DIRECT_A3D_ELEM(vol1(),0,0,0) = 1;
    DIRECT_A3D_ELEM(vol1(),0,0,1) = 2;
    DIRECT_A3D_ELEM(vol1(),0,1,0) = 3;
    DIRECT_A3D_ELEM(vol1(),0,1,1) = 4;

    DIRECT_A3D_ELEM(vol1(),1,0,0) = 5;
    DIRECT_A3D_ELEM(vol1(),1,0,1) = 6.4;
    DIRECT_A3D_ELEM(vol1(),1,1,0) = 7.5;
    DIRECT_A3D_ELEM(vol1(),1,1,1) = 8.5;
    
    //vol2
    DIRECT_A3D_ELEM(vol2(),0,0,0) = 1.5;
    DIRECT_A3D_ELEM(vol2(),0,0,1) = 2.4;
    DIRECT_A3D_ELEM(vol2(),0,1,0) = 3.3;
    DIRECT_A3D_ELEM(vol2(),0,1,1) = 4.6;

    DIRECT_A3D_ELEM(vol2(),1,0,0) = 5.7;
    DIRECT_A3D_ELEM(vol2(),1,0,1) = 6.4;
    DIRECT_A3D_ELEM(vol2(),1,1,0) = 7.3;
    DIRECT_A3D_ELEM(vol2(),1,1,1) = 8.2;
    
   */ 
    double sam = 2.;
    //std::cerr << "Vol1: " << vol1 << std::endl;
    //std::cerr << vol1() << std::endl;
    //std::cerr << "Vol2: " << vol2 << std::endl;
    //std::cerr << vol2() << std::endl;
    MultidimArray<double> freq, frc, dpr, frc_noise, error_l2;
    bool do_dpr = false;
    bool do_rfactor = true;
    double min_sam =  -1.; //no restrictionn
    double max_sam =  2; //no restriction
    double  rFactor = -1;
    frc_dpr(vol1(), vol2(), sam, freq, frc, frc_noise, dpr, error_l2, do_dpr, do_rfactor, sam/min_sam, sam/max_sam, &rFactor);

    //double Rfactor = 0.0862;
    double Rfactor = 0.134661;
    //std::cerr << "rFactor: " << fabs(Rfactor-rFactor) << " " << Rfactor << " " << rFactor << std::endl;
    ASSERT_TRUE( fabs(Rfactor-rFactor) < 0.00001);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
