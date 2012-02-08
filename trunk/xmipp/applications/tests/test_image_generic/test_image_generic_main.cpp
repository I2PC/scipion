#include <data/xmipp_image.h>
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
        filename = getenv("XMIPP_HOME");
        //get example images/staks
        imageName = filename + "/applications/tests/test_image/singleImage.spi";
        stackName = filename + "/applications/tests/test_image/smallStack.stk";
        myImageGeneric.readMapped(imageName);
        myImageGeneric2.readMapped(stackName,1);
    }

    // virtual void TearDown() {}//Destructor
    ImageGeneric myImageGeneric;
    ImageGeneric myImageGeneric2;
    Image<float> myImageFloat;
    FileName imageName;
    FileName stackName;
    FileName filename;

};


TEST_F( ImageGenericTest, equals)
{
    ImageGeneric auxImageG;
    auxImageG.readMapped(imageName);
    ASSERT_TRUE(myImageGeneric==myImageGeneric);
    ASSERT_TRUE(myImageGeneric==auxImageG);
    ASSERT_FALSE(myImageGeneric==myImageGeneric2);
}

TEST_F( ImageGenericTest, copy)
{
    ImageGeneric img1(myImageGeneric);
    ASSERT_TRUE(img1 == myImageGeneric);
    ImageGeneric img2;
    img2 = img1;
    ASSERT_TRUE(img2 == myImageGeneric);
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

TEST_F( ImageGenericTest, add)
{
    FileName auxFilename1((String)"1@"+stackName);
    FileName auxFilename2((String)"2@"+stackName);
    ImageGeneric auxImageGeneric1(auxFilename1);
    ImageGeneric auxImageGeneric2(auxFilename2);
    auxImageGeneric1.add(auxImageGeneric2);
    auxFilename2 = filename + "/applications/tests/test_image_generic/sum.spi";
    auxImageGeneric2.read(auxFilename2);
    EXPECT_TRUE(auxImageGeneric1==auxImageGeneric2);
    auxImageGeneric1.add(auxImageGeneric2);
    EXPECT_FALSE(auxImageGeneric1==auxImageGeneric2);
}

TEST_F( ImageGenericTest, subtract)
{
    FileName sumFn = filename + "/applications/tests/test_image_generic/sum.spi";
    ImageGeneric sumImg(sumFn);
    FileName fn1((String)"1@"+stackName);
    ImageGeneric img1(fn1);
    sumImg.subtract(img1);
    FileName fn2((String)"2@"+stackName);
    ImageGeneric img2(fn2);

    EXPECT_TRUE(sumImg == fn2);
}

// check if an empty file is correctly created
TEST_F( ImageGenericTest, createEmptyFile)
{
    FileName Fn(filename + "/applications/tests/test_image_generic/emptyFile.stk");
    const int size = 16;
    createEmptyFile(Fn,size,size,size,size);
    FileName Fn2;
    size_t dump;
    Fn.decompose(dump, Fn2);
    ImageGeneric Img(Fn2);

    int Xdim, Ydim, Zdim;
    size_t Ndim;
    Img.getDimensions(Xdim, Ydim, Zdim, Ndim);
    EXPECT_TRUE( Xdim == size && Ydim == size && Zdim == size && Ndim == size);
    double std, avg, min, max;
    Img().computeStats(avg, std, min, max);
    EXPECT_DOUBLE_EQ(avg, 0);
    EXPECT_DOUBLE_EQ(std, 0);
    EXPECT_DOUBLE_EQ(min, 0);
    EXPECT_DOUBLE_EQ(max, 0);
}

TEST_F( ImageGenericTest, initConstant)
{
    ImageGeneric img;
    img.setDatatype(Double);
    img.data->im->setDimensions(3,3,1,1);
    img.data->im->coreAllocateReuse();
    img.initConstant(1.);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            double p = img.getPixel(i,j);
            EXPECT_DOUBLE_EQ(p, 1.);
        }
    }
}

TEST_F( ImageGenericTest, initRandom)
{
    ImageGeneric img;
    img.setDatatype(Double);
    img.data->im->setDimensions(1024,1024,1,1);
    img.data->im->coreAllocateReuse();
    img.initRandom(0,1, RND_Gaussian);
    double mean, dev, min, max;
    img.data->computeStats(mean, dev, min, max);
    EXPECT_TRUE(ABS(mean) < 0.001);
    EXPECT_TRUE(ABS(dev-1) < 0.01);
}

// check if a pointer to data array is correctly passed
TEST_F( ImageGenericTest, getArrayPointer)
{
    ImageGeneric img, img2;
    img.read(imageName);
    MultidimArrayGeneric & mag = MULTIDIM_ARRAY_GENERIC(img);

    ImageInfo info;
    img.getInfo(info);
    img2.setDatatype(Float);
    MultidimArrayGeneric & mag2 = MULTIDIM_ARRAY_GENERIC(img2);
    mag2.resize(info.adim, false);
    float * data, *data2;
    mag.getArrayPointer(data);
    mag2.getArrayPointer(data2);

    memcpy(data2, data, info.adim.nzyxdim*sizeof(float));

    EXPECT_TRUE(img == img2);
}

// check if a pointer to MultidimArray is correctly passed
TEST_F( ImageGenericTest, getMultidimArrayPointer)
{
    ImageGeneric img, img2;
    img.read(imageName);
    MultidimArrayGeneric & mag = MULTIDIM_ARRAY_GENERIC(img);

    ImageInfo info;
    img.getInfo(info);
    img2.setDatatype(Float);
    MultidimArrayGeneric & mag2 = MULTIDIM_ARRAY_GENERIC(img2);
    mag2.resize(info.adim, false);
    MultidimArray<float> * data, *data2;
    mag.getMultidimArrayPointer(data);
    mag2.getMultidimArrayPointer(data2);

    typeCast(*data, *data2);

    EXPECT_TRUE(img == img2);
}
GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
