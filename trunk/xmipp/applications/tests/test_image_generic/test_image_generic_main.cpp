#include <data/image.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class ImageGenericTest : public ::testing::Test
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
        FileName filename(pBuf);
        filename = filename.removeFilename();
        //get example images/staks
        imageName = filename + "/../applications/tests/test_image/singleImage.spi";
        stackName = filename + "/../applications/tests/test_image/smallStack.stk";
        myImageGeneric.readMapped(imageName);
        myImageGeneric2.readMapped(stackName,1);
    }

    // virtual void TearDown() {}//Destructor
    ImageGeneric myImageGeneric;
    ImageGeneric myImageGeneric2;
    Image<float> myImageFloat;
    FileName imageName;
    FileName stackName;

};


TEST_F( ImageGenericTest, similarTo)
{
    ImageGeneric auxImageG;
    auxImageG.readMapped(imageName);
    ASSERT_TRUE(myImageGeneric==myImageGeneric);
    ASSERT_TRUE(myImageGeneric==auxImageG);
    ASSERT_FALSE(myImageGeneric==myImageGeneric2);
}
// Check if swapped images are read correctly, mapped and unmapped.
TEST_F( ImageGenericTest, readMapSwapFile)
{
    FileName auxFilename(imageName);
    ImageGeneric auxImageGeneric;
    auxFilename=auxFilename.removeExtension((String)"spi");
    auxFilename=auxFilename + "_swap.spi";
    auxImageGeneric.read(auxFilename);
    EXPECT_EQ(myImageGeneric,auxImageGeneric);

    auxImageGeneric.readMapped(auxFilename);
    EXPECT_EQ(myImageGeneric,auxImageGeneric);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
