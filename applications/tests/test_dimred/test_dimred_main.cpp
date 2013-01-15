#include <dimred/generate_data.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide

class DimRedTest : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
    }
};

TEST_F( DimRedTest, generate_data)
{
	Matrix2D<double> expectedX;
	GenerateData generator;
	generator.generateNewDataset("swiss",1000,0);

    EXPECT_EQ(expectedX,generator.X) << "DimRedTest: generate_data failed";
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
