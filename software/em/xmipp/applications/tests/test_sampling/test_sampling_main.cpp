#include "data/sampling.h"

#include <iostream>
#include <gtest/gtest.h>
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class SamplingTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
        fn_root = "sampling/";

        //Create the sampling
        double sampling_rate = 3.0;
        double angular_distance = 5.0;
        int  symmetry, sym_order;
        double max_tilt_angle = 180;
        double min_tilt_angle = 0;

        FileName fn_sym("i3h");
        FileName fnExperimentalImages(fn_root + "experimental_images.xmd");

        mysampling.setSampling(sampling_rate);//degrees
        mysampling.computeSamplingPoints(false,max_tilt_angle,min_tilt_angle);
        mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order);
        mysampling.SL.readSymmetryFile(fn_sym);
        mysampling.fillLRRepository();
        mysampling.removeRedundantPoints(symmetry, sym_order);
        mysampling.setNeighborhoodRadius(angular_distance);
        mysampling.fillExpDataProjectionDirectionByLR(fnExperimentalImages);
        mysampling.removePointsFarAwayFromExperimentalData();

    }

    FileName fn_root;
    Sampling mysampling;

    // virtual void TearDown() {}//Destructor

};

TEST_F(SamplingTest, computeSamplingPoints)
{
    XMIPP_TRY
    Sampling s1;
    //by default no sampling points are read
    s1.readSamplingFile(fn_root + "ref");
    Sampling s2;
    s2.setSampling(3.0);
    s2.computeSamplingPoints(false, 180., 0.);
    //s1.saveSamplingFile("/tmp/s1");
    //s2.saveSamplingFile("/tmp/s2");
    EXPECT_EQ(s1, s2);
    XMIPP_CATCH
}

TEST_F(SamplingTest, removeRedundantPointsI3H)
{
    XMIPP_TRY
    int  symmetry, sym_order;
    Sampling s1;//<- check number of sampling points...
    s1.readSamplingFile(fn_root + "ref_i3h",true,false);
    FileName fn_sym("i3h");
    Sampling s2;
    s2.setSampling(3.);
    s2.computeSamplingPoints(false, 180, 0);
    s2.SL.isSymmetryGroup(fn_sym, symmetry, sym_order);
    s2.SL.readSymmetryFile(fn_sym);
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    //s2.saveSamplingFile("/tmp/s2");
    //s1.saveSamplingFile("/tmp/s1");
    EXPECT_EQ(s1, s2);
    XMIPP_CATCH
}

TEST_F(SamplingTest, removeRedundantPointsC1)
{
    XMIPP_TRY
    int  symmetry, sym_order;
    Sampling s1;
    s1.readSamplingFile(fn_root + "ref_c1");
    Sampling s2;
    s2.setSampling(3.);
    s2.computeSamplingPoints(false, 180, 0.);
    s2.SL.isSymmetryGroup("c1", symmetry, sym_order);
    s2.SL.readSymmetryFile("c1");
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    //s2.saveSamplingFile("/tmp/removeRedundantPointsC1");
    EXPECT_EQ(s1, s2);
    XMIPP_CATCH
}

TEST_F(SamplingTest, removePointsFarAwayFromExperimentalDataI3H)
{
    XMIPP_TRY
    int  symmetry, sym_order;
    Sampling s1;
    s1.readSamplingFile(fn_root + "ref_i3h_exp");
    Sampling s2;
    s2.setSampling(3);//degrees
    s2.computeSamplingPoints(false,180.,0.);
    s2.SL.isSymmetryGroup("i3h", symmetry, sym_order);
    s2.SL.readSymmetryFile("i3h");
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    s2.setNeighborhoodRadius(5);
    s2.fillExpDataProjectionDirectionByLR(fn_root + "experimental_images.xmd");
    s2.removePointsFarAwayFromExperimentalData();
    EXPECT_EQ(s1, s2);
    XMIPP_CATCH
}

TEST_F(SamplingTest, removePointsFarAwayFromExperimentalDataC1)
{
    XMIPP_TRY
    int  symmetry, sym_order;
    Sampling s1;
    s1.readSamplingFile(fn_root + "ref_c1_exp");
    Sampling s2;
    s2.setSampling(3);//degrees
    s2.computeSamplingPoints(false,180.,0.);
    s2.SL.isSymmetryGroup("c1", symmetry, sym_order);
    s2.SL.readSymmetryFile("c1");
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    s2.setNeighborhoodRadius(5);
    s2.fillExpDataProjectionDirectionByLR(fn_root + "experimental_images.xmd");
    s2.removePointsFarAwayFromExperimentalData();
    //s2.saveSamplingFile("/tmp/removePointsFarAwayFromExperimentalDataC1");
    EXPECT_EQ(s1, s2);
    XMIPP_CATCH
}

TEST_F(SamplingTest, saveReadSamplingFile)
{
    XMIPP_TRY
    FileName fn;
    fn.initUniqueName("/tmp/temp_XXXXXX");
    mysampling.saveSamplingFile(fn, true, true);
    Sampling s;
    s.readSamplingFile(fn, true, true);
    EXPECT_EQ(mysampling, s);
    fn.deleteFile();
    fn = fn + "_sampling.xmd";
    fn.deleteFile();
    XMIPP_CATCH
}

TEST_F(SamplingTest, computeNeighborsI3H)
{
    XMIPP_TRY
    int  symmetry, sym_order;
    Sampling s1;
    s1.readSamplingFile(fn_root + "neigh_ref_i3h_exp");
    Sampling s2;
    s2.setSampling(3);//degrees
    s2.computeSamplingPoints(false,180.,0.);
    s2.SL.isSymmetryGroup("i3h", symmetry, sym_order);
    s2.SL.readSymmetryFile("i3h");
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    s2.setNeighborhoodRadius(5);
    s2.fillExpDataProjectionDirectionByLR(fn_root + "experimental_images.xmd");
    s2.removePointsFarAwayFromExperimentalData();
    s2.computeNeighbors();
    //s2.saveSamplingFile("/tmp/computeNeighborsI3H");
    EXPECT_EQ(s1, s2);
    XMIPP_CATCH
}

TEST_F(SamplingTest, computeNeighborsC1)
{
    XMIPP_TRY
    int  symmetry, sym_order;
    Sampling s1;
    s1.readSamplingFile(fn_root + "neigh_ref_c1_exp");
    Sampling s2;
    s2.setSampling(3);//degrees
    s2.computeSamplingPoints(false,180.,0.);
    s2.SL.isSymmetryGroup("c1", symmetry, sym_order);
    s2.SL.readSymmetryFile("c1");
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    s2.setNeighborhoodRadius(5);
    s2.fillExpDataProjectionDirectionByLR(fn_root + "experimental_images.xmd");
    s2.removePointsFarAwayFromExperimentalData();
    s2.computeNeighbors();
    //    s2.saveSamplingFile(fn_root + "/tmp/ref_c1_computeNeighborsC1");
    EXPECT_EQ(s1, s2);
    XMIPP_CATCH
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
