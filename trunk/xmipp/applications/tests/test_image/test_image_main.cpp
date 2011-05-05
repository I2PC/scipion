#include <data/image.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class ImageTest : public ::testing::Test
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
        myImage.read(imageName);
        myStack.read(stackName);
    }

    // virtual void TearDown() {}//Destructor
    Image<double> myImage;
    Image<double> myStack;
    FileName imageName;
    FileName stackName;

};


TEST_F( ImageTest, similarTo)
{
    ASSERT_TRUE(myImage==myImage);
    ASSERT_FALSE(myImage==myStack);
}

TEST_F( ImageTest, writeIMAGICimage)
{
    FileName auxFilename(imageName);
    auxFilename=auxFilename.removeExtension((String)"spi");
    auxFilename=auxFilename.addExtension("img");
    myImage.write(auxFilename);
    Image<double> auxImage;
    auxImage.read(auxFilename);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, writeIMAGICstack)
{
    FileName auxFilename(stackName);
    auxFilename=auxFilename.removeExtension((String)"stk");
    auxFilename=auxFilename.addExtension("img");
    myStack.write(auxFilename);
    Image<double> auxStack;
    auxStack.read(auxFilename);
    EXPECT_EQ(myStack,auxStack);
}

TEST_F( ImageTest, writeMRCimage)
{
    FileName auxFilename(imageName);
    auxFilename=auxFilename.removeExtension((String)"spi");
    auxFilename=auxFilename.addExtension("mrc");
    myImage.write(auxFilename);
    Image<double> auxImage;
    auxImage.read(auxFilename);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, writeMRCstack)
{
    FileName auxFilename(stackName);
    auxFilename=auxFilename.removeExtension((String)"stk");
    auxFilename=auxFilename.addExtension("mrcs");
    myStack.write(auxFilename);
    Image<double> auxStack;
    auxStack.read(auxFilename);
    EXPECT_EQ(myStack,auxStack);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
