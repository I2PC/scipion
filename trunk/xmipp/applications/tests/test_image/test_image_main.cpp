
#include <stdlib.h>
#include <data/xmipp_image.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
#include <data/metadata.h>
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class ImageTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        //get example images/staks
        baseName = getenv("XMIPP_HOME");
        imageName = baseName + "/applications/tests/test_image/singleImage.spi";
        stackName = baseName + "/applications/tests/test_image/smallStack.stk";
        myImage.read(imageName);
        myStack.read(stackName);
    }

    // virtual void TearDown() {}//Destructor
    Image<double> myImage;
    Image<double> myStack;
    FileName imageName;
    FileName stackName;
    FileName baseName;

};


TEST_F( ImageTest, similarTo)
{
    ASSERT_TRUE(myImage==myImage);
    ASSERT_FALSE(myImage==myStack);
}

TEST_F( ImageTest, getEulerAngles)
{
    Image<double> auxImage;
    auxImage.setEulerAngles(10., 20., 30.);
    double rot, psi, tilt;
    auxImage.getEulerAngles(rot, tilt, psi);
    EXPECT_DOUBLE_EQ(rot, 10.);
    EXPECT_DOUBLE_EQ(tilt, 20.);
    EXPECT_DOUBLE_EQ(psi, 30.);
}

TEST_F( ImageTest, readApplyGeo)
{
    FileName auxFilename(baseName + "/applications/tests/test_image/test2.spi");
    MetaData MD;
    size_t id = MD.addObject();
    MD.setValue(MDL_IMAGE, auxFilename, id);
    MD.setValue(MDL_ANGLEPSI, 45., id);
    Image<double> auxImage, auxImage2;
    auxImage.readApplyGeo(MD,id, false, DATA, ALL_IMAGES, false);
    auxImage2.read(auxFilename.removeAllExtensions()+"_wrap_false.spi");
    EXPECT_TRUE(auxImage == auxImage2);
    auxImage.readApplyGeo(MD,id, false, DATA, ALL_IMAGES, true);
    auxImage2.read(auxFilename.removeAllExtensions()+"_wrap_true.spi");
    EXPECT_TRUE(auxImage == auxImage2);
}

TEST_F( ImageTest, readImageFromStackMetadata)
{
    FileName stackSliceFn, auxFn;
    stackSliceFn.compose(2, stackName);
    Image<double> img1;
    img1.read(stackSliceFn);
    MetaData md(stackSliceFn);
    size_t id = md.firstObject();
    md.getValue(MDL_IMAGE, auxFn, id);
    Image<double> img2;
    img2.read(auxFn);

    EXPECT_TRUE(img1 == img2);
}

//ROB ask kino is angles are saved in header
TEST_F( ImageTest, saveImageinStackwithHeaderAngleRot)
{
    FileName stackSliceFn, auxFn;
    stackSliceFn.compose(2, stackName);
    Image<double> img1;
    img1.read(stackSliceFn);
    stackSliceFn.compose(2, "/tmp/saveImageinStackwithHeaderAngleRot.stk");
    img1.setEulerAngles(10.,20.,30.);
    img1.write(stackSliceFn);
    img1.clear();
    img1.read(stackSliceFn);

    double rot,tilt,psi;
    img1.getEulerAngles(rot,tilt,psi);
    unlink(stackSliceFn.c_str());
    EXPECT_DOUBLE_EQ(10. ,rot);
    EXPECT_DOUBLE_EQ(20.,tilt);
    EXPECT_DOUBLE_EQ(30. ,psi);
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

TEST_F( ImageTest, writeMRCstack)//show -i kk.mrcs for stacks fails for mrc
//ml_tomo anotate bugs
{
    FileName auxFilename(stackName);
    auxFilename=auxFilename.removeExtension((String)"stk");
    auxFilename=auxFilename.addExtension("mrcs");
    myStack.write(auxFilename);
    Image<double> auxStack;
    auxStack.read(auxFilename);
    EXPECT_EQ(myStack,auxStack);
}

TEST_F( ImageTest, writeTIFimage)
{
    FileName auxFilename(imageName);
    auxFilename=auxFilename.removeExtension((String)"spi");
    auxFilename=auxFilename.addExtension("tif");
    myImage.write(auxFilename);
    Image<double> auxImage;
    auxImage.read(auxFilename);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, writeINFimage)
{
    FileName auxFilename(imageName);
    auxFilename=auxFilename.removeExtension((String)"spi");
    auxFilename=auxFilename.addExtension("inf");
    myImage.write(auxFilename);
    Image<double> auxImage;
    auxImage.read(auxFilename);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, writeRAWimage)
{
    FileName auxFilename(imageName);
    auxFilename=auxFilename.removeExtension((String)"spi");
    auxFilename=auxFilename.addExtension("raw#3,3");
    myImage.write(auxFilename);
    Image<double> auxImage;
    auxImage.read(auxFilename);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, readPreview)
{
    FileName auxFilename(baseName + "/applications/tests/test_image/smallVolume.vol");
    Image<double> img1, img2;
    img1.read(auxFilename);

    img1().setXmippOrigin();
    selfScaleToSize(NEAREST, img1(),32,32,4);

    img2.readPreview(auxFilename, 32,32, ALL_SLICES);
    img1().setXmippOrigin();
    img2().setXmippOrigin();

    EXPECT_TRUE(img1 == img2);
}

TEST_F( ImageTest, mapFile2Write)
{
    FileName auxFilename(baseName + "/applications/tests/test_image/smallVolume.vol");
    FileName auxMappedFilename(baseName + "/applications/tests/test_image/mappedFile.vol");
    Image<float> img1, img2;
    img1.read(auxFilename);
    ArrayDim aDim;
    img1().getDimensions(aDim);

    auxMappedFilename.deleteFile();
    img2.mapFile2Write(aDim.xdim, aDim.ydim, aDim.zdim, auxMappedFilename);
    typeCast(img1(), img2());
    img2.write(auxMappedFilename);
    img2.clear();
    img2.read(auxMappedFilename);

    EXPECT_TRUE(img1 == img2);
}
TEST_F( ImageTest, movePointerToSlice)
{
    FileName auxFilename(baseName + "/applications/tests/test_image/smallVolume.vol");
    Image<double> img1, img2;
    img1.read(auxFilename);
    img1().setXmippOrigin();
    selfScaleToSize(NEAREST, img1(),32,32,4);

    ArrayDim aDim;
    img1().getDimensions(aDim);
    img2.readPreview(auxFilename, 32,32, ALL_SLICES);

    for (int k = 1; k <= aDim.zdim; ++k)
    {
        img1.movePointerToSlice(k);
        img2.movePointerToSlice(k);

        img1().setXmippOrigin();
        img2().setXmippOrigin();

        EXPECT_TRUE(img1 == img2);
    }
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
