#include <data/image.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class FuncTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        a=0;
    }

    // virtual void TearDown() {}//Destructor
    int a;
};


TEST_F( FuncsTest, CompareTwoFiles)
{
    ASSERT_TRUE(myImage==myImage);
    ASSERT_FALSE(myImage==myStack);
}


GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
