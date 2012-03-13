
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
        xmippPath = getXmippPath();
       // testBaseName = xmippPath + "/resources/test";
        imageName = TEST_FILENAME("singleImage.spi");
        stackName = TEST_FILENAME("smallStack.stk");
        myImage.read(imageName);
        myStack.read(stackName);
    }

    // virtual void TearDown() {}//Destructor
    Image<double> myImage;
    Image<double> myStack;
    FileName imageName;
    FileName stackName;
    FileName xmippPath;

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
    FileName auxFn = TEST_FILENAME("test2.spi");
    MetaData MD;
    size_t id = MD.addObject();
    MD.setValue(MDL_IMAGE, auxFn, id);
    MD.setValue(MDL_ANGLEPSI, 45., id);
    Image<double> auxImage, auxImage2;
    ApplyGeoParams params;
    params.wrap = false;
    auxImage.readApplyGeo(MD,id, params);
    auxImage2.read(auxFn.insertBeforeExtension("_wrap_false"));
    EXPECT_TRUE(auxImage == auxImage2);
    params.wrap = true;
    auxImage.readApplyGeo(MD,id, params);
    auxImage2.read(auxFn.insertBeforeExtension("_wrap_true"));
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
    img1.setDataMode(_DATA_ALL);
    img1.write(stackSliceFn);
    img1.clear();
    img1.read(stackSliceFn, _DATA_ALL);

    double rot,tilt,psi;
    img1.getEulerAngles(rot,tilt,psi);
    unlink(stackSliceFn.c_str());
    EXPECT_DOUBLE_EQ(10. ,rot);
    EXPECT_DOUBLE_EQ(20.,tilt);
    EXPECT_DOUBLE_EQ(30. ,psi);
}

TEST_F( ImageTest, writeIMAGICimage)
{
    FileName auxFn(imageName);
    auxFn=auxFn.removeExtension((String)"spi");
    auxFn=auxFn.addExtension("img");
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, writeIMAGICstack)
{
    FileName auxFn(stackName);
    auxFn=auxFn.removeExtension((String)"stk");
    auxFn=auxFn.addExtension("img");
    myStack.write(auxFn);
    Image<double> auxStack;
    auxStack.read(auxFn);
    EXPECT_EQ(myStack,auxStack);
}

TEST_F( ImageTest, writeMRCimage)
{
    FileName auxFn(imageName);
    auxFn=auxFn.removeExtension((String)"spi");
    auxFn=auxFn.addExtension("mrc");
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, writeMRCstack)//show -i kk.mrcs for stacks fails for mrc
//ml_tomo anotate bugs
{
    FileName auxFn(stackName);
    auxFn=auxFn.removeExtension((String)"stk");
    auxFn=auxFn.addExtension("mrcs");
    myStack.write(auxFn);
    Image<double> auxStack;
    auxStack.read(auxFn);
    EXPECT_EQ(myStack,auxStack);
}

TEST_F( ImageTest, writeTIFimage)
{
    FileName auxFn(imageName);
    auxFn=auxFn.removeExtension((String)"spi");
    auxFn=auxFn.addExtension("tif");
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, writeINFimage)
{
    FileName auxFn(imageName);
    auxFn=auxFn.removeExtension((String)"spi");
    auxFn=auxFn.addExtension("inf");
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, writeRAWimage)
{
    FileName auxFn(imageName);
    auxFn=auxFn.removeExtension((String)"spi");
    auxFn=auxFn.addExtension("raw#3,3");
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
}

TEST_F( ImageTest, readPreview)
{
    FileName auxFn = TEST_FILENAME("smallVolume.vol");
    Image<double> img1, img2;
    img1.read(auxFn);

    img1().setXmippOrigin();
    selfScaleToSize(NEAREST, img1(),32,32,4);

    img2.readPreview(auxFn, 32,32, ALL_SLICES);
    img1().setXmippOrigin();
    img2().setXmippOrigin();

    EXPECT_TRUE(img1 == img2);
}

TEST_F( ImageTest, getPreview)
{
    FileName auxFn = TEST_FILENAME("smallVolume.vol");
    Image<double> img1, img2;
    img1.read(auxFn);

    img1.getPreview(&img2, 32,32, ALL_SLICES);

    img1().setXmippOrigin();
    selfScaleToSize(NEAREST, img1(),32,32,4);

    img1().setXmippOrigin();
    img2().setXmippOrigin();
    EXPECT_TRUE(img1 == img2);
}

TEST_F( ImageTest, mapFile2Write)
{
     FileName auxFn = TEST_FILENAME("smallVolume.vol");
    FileName auxMappedFilename = TEST_FILENAME("mappedFile.vol");
    Image<float> img1, img2;
    img1.read(auxFn);
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
TEST_F( ImageTest, movePointerTo)
{
    FileName auxFn= TEST_FILENAME("smallVolumeStack.stk");
    Image<double> img1, img2;
    img1.read(auxFn);

    ArrayDim aDim;
    img1().getDimensions(aDim);

    for (size_t n = 1; +n <= aDim.ndim; ++n)
    {
        img2.read(auxFn, DATA, n);

        for (int k = 1; k <= aDim.zdim; ++k)
        {
            img1.movePointerTo(k, n);
            img2.movePointerTo(k);

            img1().setXmippOrigin();
            img2().setXmippOrigin();

            EXPECT_TRUE(img1 == img2);
        }
    }
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
