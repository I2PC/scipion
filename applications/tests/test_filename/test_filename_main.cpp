#include <iostream>
#include "data/xmipp_filename.h"
#include "data/xmipp_strings.h"

#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class FileNameTest : public ::testing::Test
{
protected:

    // virtual void TearDown() {}//Destructor
    FileName source1;
    FileName source2;
    //init metadatas
    virtual void SetUp()
    {

    }
};


TEST_F( FileNameTest, getBlockName)
{
    String empty = "";
    String a = "a";
    FileName fn;
    fn = "abc.xmd";
    ASSERT_EQ(fn.getBlockName(), empty);
    fn = "@abc.xmd";
    ASSERT_EQ(fn.getBlockName(), empty);
    fn = "1@abc.xmd";
    ASSERT_EQ(fn.getBlockName(), empty);
    fn = "a@abc.xmd";
    ASSERT_EQ(fn.getBlockName(), a);
    fn = "1,a@abc.xmd";
    ASSERT_EQ(fn.getBlockName(), a);
    fn = "1,@abc.xmd";
    ASSERT_EQ(fn.getBlockName(), empty);

}


TEST_F( FileNameTest, removeBlockName)
{
    String empty = "";
    String a = "a";
    FileName fn;
    fn = "abc.xmd";
    ASSERT_EQ(fn.removeBlockName(), fn);
    fn = "@abc.xmd";
    ASSERT_EQ(fn.removeBlockName(), "abc.xmd");
    fn = "1@abc.xmd";
    ASSERT_EQ(fn.removeBlockName(), fn);
    fn = "a@abc.xmd";
    ASSERT_EQ(fn.removeBlockName(), String("abc.xmd"));
    fn = "1,a@abc.xmd";
    ASSERT_EQ(fn.removeBlockName(), "1@abc.xmd");
    fn = "1,@abc.xmd";
    ASSERT_EQ(fn.removeBlockName(), "1@abc.xmd");

}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
