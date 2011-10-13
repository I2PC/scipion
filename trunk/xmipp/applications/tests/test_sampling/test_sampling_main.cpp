#include "data/sampling.h"

#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class SamplingTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
#define len 128
        //find binaries directory
        char szTmp[len];
        char pBuf[len];
        sprintf(szTmp, "/proc/%d/exe", getpid());
        int bytes = std::min(readlink(szTmp, pBuf, len), (ssize_t)len - 1);
        if(bytes >= 0)
            pBuf[bytes] = '\0';
        //remove last token
        FileName bin_path = FileName(pBuf).removeFilename();
        fn_root = bin_path + "/../applications/tests/test_sampling/";

        //Create the sampling
        double sampling_rate = 3.0;
        double angular_distance = 5.0;
        int  symmetry, sym_order;
        double max_tilt_angle = 91;
        double min_tilt_angle = -91;
        FileName fn_sym("i3h");
        FileName fnExperimentalImages(fn_root + "experimental_images.xmd");

        mysampling.setSampling(sampling_rate);//degrees
        mysampling.computeSamplingPoints(false,max_tilt_angle,min_tilt_angle);
        mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order);
        mysampling.SL.read_sym_file(fn_sym);
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
    Sampling s1;
    s1.readSamplingFile(fn_root + "ref");
    s1.saveSamplingFile(fn_root + "refkk");
    Sampling s2;
    s2.setSampling(3.);
    s2.computeSamplingPoints(false, 91., -91.);
    EXPECT_EQ(s1, s2);
}

TEST_F(SamplingTest, removeRedundantPoints)
{
    int  symmetry, sym_order;
    Sampling s1;
    s1.readSamplingFile(fn_root + "ref_i3h");
    Sampling s2;
    s2.setSampling(3.);
    s2.computeSamplingPoints(false, 91., -91.);
    s2.SL.isSymmetryGroup("i3h", symmetry, sym_order);
    s2.SL.read_sym_file("i3h");
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    EXPECT_EQ(s1, s2);
}
TEST_F(SamplingTest, removePointsFarAwayFromExperimentalData)
{
    int  symmetry, sym_order;
    Sampling s1;
    s1.readSamplingFile(fn_root + "ref_i3h_exp");
    Sampling s2;
    s2.setSampling(3);//degrees
    s2.computeSamplingPoints(false,91.,-91.);
    s2.SL.isSymmetryGroup("i3h", symmetry, sym_order);
    s2.SL.read_sym_file("i3h");
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    s2.setNeighborhoodRadius(5);
    s2.fillExpDataProjectionDirectionByLR(fn_root + "experimental_images.xmd");
    s2.removePointsFarAwayFromExperimentalData();
    EXPECT_EQ(s1, s2);
}

TEST_F(SamplingTest, saveReadSamplingFile)
{
  FileName fn = fn_root + "test";
    mysampling.saveSamplingFile(fn, true, true);
    Sampling s;
    s.readSamplingFile(fn, true, true);
    EXPECT_EQ(mysampling, s);
}

TEST_F(SamplingTest, computeNeighbors)
{
    int  symmetry, sym_order;
    Sampling s1;
    s1.readSamplingFile(fn_root + "ref_i3h_exp");
    Sampling s2;
    s2.setSampling(3);//degrees
    s2.computeSamplingPoints(false,91.,-91.);
    s2.SL.isSymmetryGroup("i3h", symmetry, sym_order);
    s2.SL.read_sym_file("i3h");
    s2.fillLRRepository();
    s2.removeRedundantPoints(symmetry, sym_order);
    s2.setNeighborhoodRadius(5);
    s2.fillExpDataProjectionDirectionByLR(fn_root + "experimental_images.xmd");
    s2.removePointsFarAwayFromExperimentalData();
    s2.computeNeighbors();
    s2.saveSamplingFile(fn_root + "ref_i3h_expkkkkk");
    EXPECT_EQ(s1, s2);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
