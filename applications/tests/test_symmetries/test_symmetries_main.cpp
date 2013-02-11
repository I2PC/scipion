#include "data/sampling.h"

#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class SamplingTest : public ::testing::Test
{
protected:
    //init symmetries
    virtual void SetUp()
    {
        //there is some overlapping with test_sampling
        chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
        fn_root = "symmetries/";
        FileName fnExperimentalImages(fn_root + "experimental_images.xmd");

        //SL.isSymmetryGroup(fn_sym, symmetry, pg_I3);
        //SL.readSymmetryFile(fn_sym);
    }

    FileName fn_root;
    SymList SL;

    // virtual void TearDown() {}//Destructor

};

TEST_F(SamplingTest, isSymmetryGroup)
{
    XMIPP_TRY
    int symmetry, sym_order;
    FileName fn_sym("i3h");
    SL.isSymmetryGroup(fn_sym, symmetry, sym_order);
    EXPECT_EQ(symmetry, pg_I3H);
    EXPECT_EQ(sym_order, -1);
    fn_sym="c5";
    SL.isSymmetryGroup(fn_sym, symmetry, sym_order);
    EXPECT_EQ(symmetry, pg_CN);
    EXPECT_EQ(sym_order, 5);

    XMIPP_CATCH
}

TEST_F(SamplingTest, readSymmetryFile)
{
    XMIPP_TRY
    int trueSymsNo;
    FileName fn_sym("i3h");
    SL.readSymmetryFile(fn_sym);
    trueSymsNo = SL.trueSymsNo();//true_symNo;
    EXPECT_EQ(trueSymsNo, 119);
    XMIPP_CATCH
}

TEST_F(SamplingTest, computeDistance)
{
    XMIPP_TRY
    FileName fn_sym("i3h");
    SL.readSymmetryFile(fn_sym);
    double rot2=6.;
    double tilt2=5.;
    double psi2=4.;
    double total;
    total = SL.computeDistance(1., 2., 3., rot2, tilt2, psi2,false,false,false);
    EXPECT_NEAR (total, 5.23652,0.00001);
    XMIPP_CATCH
}

TEST_F(SamplingTest, computeDistanceMetadata)
{
    XMIPP_TRY
    FileName fn_sym("i3h");
    SL.readSymmetryFile(fn_sym);
    MetaData md,mdOut;
    MDRow row;
    row.setValue(MDL_ANGLE_ROT,1.);
    row.setValue(MDL_ANGLE_TILT,2.);
    row.setValue(MDL_ANGLE_PSI,3.);
    row.setValue(MDL_ANGLE_ROT2,6.);
    row.setValue(MDL_ANGLE_TILT2,5.);
    row.setValue(MDL_ANGLE_PSI2,4.);
    md.addRow(row);
    row.setValue(MDL_ANGLE_ROT,11.);
    row.setValue(MDL_ANGLE_TILT,12.);
    row.setValue(MDL_ANGLE_PSI,13.);
    row.setValue(MDL_ANGLE_ROT2,16.);
    row.setValue(MDL_ANGLE_TILT2,15.);
    row.setValue(MDL_ANGLE_PSI2,14.);
    md.addRow(row);
    mdOut = md;
    SL.computeDistance(md,false,false,false);
    double total;
    size_t id = md.firstObject();
    md.getValue(MDL_ANGLE_DIFF,total,id);
    EXPECT_NEAR (total, 5.23652,0.00001);
    XMIPP_CATCH
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
