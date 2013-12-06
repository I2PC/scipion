
#include <stdlib.h>
#include <data/xmipp_image.h>
#include <data/xmipp_image_extension.h>
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
        XMIPP_TRY
        //get example images/staks
        xmippPath = getXmippPath();
        if (chdir(((String)(xmippPath + "/resources/test")).c_str())==-1)
        	REPORT_ERROR(ERR_UNCLASSIFIED,"Could not change directory");
        // testBaseName = xmippPath + "/resources/test";
        imageName = "image/singleImage.spi";
        stackName = "image/smallStack.stk";
        stackVolName= "image/smallVolumeStack.stk";
        myImage.read(imageName);
        myStack.read(stackName);
        myVolStack.read(stackVolName);
        XMIPP_CATCH
    }

    // virtual void TearDown() {}//Destructor
    Image<double> myImage;
    Image<double> myStack;
    Image<double> myVolStack;
    FileName imageName;
    FileName stackName;
    FileName stackVolName;
    FileName xmippPath;

};


TEST_F( ImageTest, similarTo)
{
    ASSERT_TRUE(myImage==myImage);
    ASSERT_FALSE(myImage==myStack);
}

TEST_F( ImageTest, copy)
{
    Image<double> img2(myImage);
    ASSERT_TRUE(img2==myImage);
    ArrayDim aDim1, aDim2;
    myImage.getDimensions(aDim1);
    img2.getDimensions(aDim2);
    ASSERT_TRUE(aDim1 == aDim2);

    Image<double> empty;
    empty().initZeros(32, 32);
    img2 = empty;
    ASSERT_TRUE(empty==img2);

    empty.getDimensions(aDim1);
    img2.getDimensions(aDim2);
    ASSERT_TRUE(aDim1 == aDim2);
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
    XMIPP_TRY
    FileName auxFn = "image/test2.spi";
    MetaData MD;
    size_t id = MD.addObject();
    MD.setValue(MDL_IMAGE, auxFn, id);
    MD.setValue(MDL_ANGLE_PSI, 45., id);
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
    XMIPP_CATCH
}

TEST_F( ImageTest, readImageFromStackMetadata)
{
    XMIPP_TRY
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
    XMIPP_CATCH
}

//ROB ask kino is angles are saved in header
TEST_F( ImageTest, saveImageinStackwithHeaderAngleRot)
{
    XMIPP_TRY
    FileName stackSliceFn, auxFn;
    stackSliceFn.compose(2, stackName);
    Image<double> img1;
    img1.read(stackSliceFn);

    stackSliceFn.initUniqueName("/tmp/temp_stk_XXXXXX");
    stackSliceFn = String("2@") + stackSliceFn  + ":stk";
    img1.setEulerAngles(10.,20.,30.);
    img1.setDataMode(_DATA_ALL);
    img1.write(stackSliceFn);
    img1.clear();
    img1.read(stackSliceFn, _DATA_ALL);

    double rot,tilt,psi;
    img1.getEulerAngles(rot,tilt,psi);
    EXPECT_DOUBLE_EQ(10. ,rot);
    EXPECT_DOUBLE_EQ(20.,tilt);
    EXPECT_DOUBLE_EQ(30. ,psi);
    stackSliceFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, writeIMAGICimage)
{
    XMIPP_TRY
    FileName auxFn;
    auxFn.initUniqueName("/tmp/temp_img_XXXXXX");
    auxFn.deleteFile();
    auxFn = auxFn + ".img";
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
    auxFn.deleteFile();
    auxFn = auxFn.replaceSubstring(".img", ".hed");
    auxFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, mirrorY)
{
    XMIPP_TRY
    Image<double> auxImage(3,3);
    Image<double> auxImageY(3,3);
    int dim ;
    dim=3;//odd
    for (int i=0; i< dim; i++)
        for (int j=0; j< dim ; j++)
            DIRECT_IMGPIXEL(auxImage,i,j)=(double) dim*i+j;
    for (int i=0; i< dim; i++)
        for (int j=0; j< dim ; j++)
            DIRECT_IMGPIXEL(auxImageY,dim-i-1,j)=(double) dim*i+j;
    auxImage.mirrorY();
    EXPECT_EQ(auxImage,auxImageY);
    dim=4;//even
    Image<double> auxImage2(4,4);
    Image<double> auxImage2Y(4,4);
    for (int i=0; i< dim; i++)
        for (int j=0; j< dim ; j++)
            DIRECT_IMGPIXEL(auxImage2,i,j)=(double) dim*i+j;
    for (int i=0; i< dim; i++)
        for (int j=0; j< dim ; j++)
            DIRECT_IMGPIXEL(auxImage2Y,dim-i-1,j)=(double) dim*i+j;
    auxImage2.mirrorY();
    EXPECT_EQ(auxImage2,auxImage2Y);
    XMIPP_CATCH
}

TEST_F( ImageTest, writeIMAGICstack)
{
    XMIPP_TRY
    FileName auxFn;
    auxFn.initUniqueName("/tmp/temp_imgstk_XXXXXX");
    auxFn.deleteFile();
    auxFn = auxFn + ".img";
    myStack.write(auxFn);
    Image<double> auxStack;
    auxStack.read(auxFn);
    EXPECT_EQ(myStack,auxStack);
    auxFn.deleteFile();
    auxFn = auxFn.replaceSubstring(".img", ".hed");
    auxFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, writeMRCimage)
{
    XMIPP_TRY
    FileName auxFn;
    auxFn.initUniqueName("/tmp/temp_mrc_XXXXXX");
    auxFn = auxFn + ":mrc";
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
    auxFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, writeMRCstack)
{
    XMIPP_TRY
    FileName auxFn;
    auxFn.initUniqueName("/tmp/temp_mrcstk_XXXXXX");
    auxFn = auxFn + ":mrcs";
    myStack.write(auxFn);
    Image<double> auxStack;
    auxStack.read(auxFn);
    EXPECT_EQ(myStack,auxStack);
    auxFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, writeMRCVOLstack)
{
    XMIPP_TRY
    FileName auxFn;
    auxFn.initUniqueName("/tmp/temp_mrcvolstk_XXXXXX");
    auxFn = auxFn + ":mrc";
    myVolStack.write(auxFn);
    Image<double> auxStack;
    auxStack.read(auxFn);
    EXPECT_EQ(myVolStack,auxStack);
    ArrayDim auxStackArrayDim;
    ArrayDim StackArrayDim;
    myVolStack.getDimensions(StackArrayDim);
    auxStack.getDimensions(auxStackArrayDim);
    EXPECT_TRUE(StackArrayDim==auxStackArrayDim);
    auxFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, writeTIFimage)
{
    XMIPP_TRY
    FileName auxFn;
    auxFn.initUniqueName("/tmp/temp_tif_XXXXXX");
    auxFn = auxFn + ":tif";
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
    auxFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, writeINFimage)
{
    XMIPP_TRY
    FileName auxFn;
    auxFn.initUniqueName("/tmp/temp_inf_XXXXXX");
    auxFn.deleteFile();
    auxFn = auxFn + ".raw";
    myImage.write(auxFn);
    Image<double> auxImage;
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
    auxFn.deleteFile();
    auxFn = auxFn.addExtension("inf");
    auxFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, readRAWimage)
{
    XMIPP_TRY
    FileName auxFn;
    auxFn.initUniqueName("/tmp/temp_raw_XXXXXX");
    myImage.write(auxFn + ":spi");
    Image<double> auxImage;
    auxFn = auxFn + "#3,3,1032,float";
    auxImage.read(auxFn);
    EXPECT_EQ(myImage,auxImage);
    auxFn.deleteFile();
    XMIPP_CATCH
}

TEST_F( ImageTest, readPreview)
{
    XMIPP_TRY
    FileName auxFn = "image/smallVolume.vol";
    Image<double> img1, img2;
    img1.read(auxFn);

    img1().setXmippOrigin();
    selfScaleToSize(NEAREST, img1(),32,32,4);

    img2.readPreview(auxFn, 32,32, ALL_SLICES);
    img1().setXmippOrigin();
    img2().setXmippOrigin();

    EXPECT_TRUE(img1 == img2);
    XMIPP_CATCH
}

TEST_F( ImageTest, getPreview)
{
    XMIPP_TRY
    FileName auxFn = "image/smallVolume.vol";
    Image<double> img1, img2;
    img1.read(auxFn);

    img1.getPreview(&img2, 32,32, ALL_SLICES);

    img1().setXmippOrigin();
    selfScaleToSize(NEAREST, img1(),32,32,4);

    img1().setXmippOrigin();
    img2().setXmippOrigin();
    EXPECT_TRUE(img1 == img2);
    XMIPP_CATCH
}

TEST_F( ImageTest, mapFile2Write)
{
    XMIPP_TRY
    FileName auxFn = "image/smallVolume.vol";
    FileName auxMappedFilename;
    auxMappedFilename.initUniqueName("/tmp/temp_vol_XXXXXX");
    auxMappedFilename = auxMappedFilename + ":vol";
    Image<float> img1, img2;
    img1.read(auxFn);
    ArrayDim aDim;
    img1().getDimensions(aDim);

    img2.mapFile2Write(aDim.xdim, aDim.ydim, aDim.zdim, auxMappedFilename);
    typeCast(img1(), img2());
    img2.write(auxMappedFilename);
    img2.clear();
    img2.read(auxMappedFilename);

    EXPECT_TRUE(img1 == img2);

    auxMappedFilename.deleteFile();
    XMIPP_CATCH
}
TEST_F( ImageTest, movePointerTo)
{
    XMIPP_TRY
    FileName auxFn= "image/smallVolumeStack.stk";
    Image<double> img1, img2;
    img1.read(auxFn);

    ArrayDim aDim;
    img1().getDimensions(aDim);

    for (size_t n = 1; +n <= aDim.ndim; ++n)
    {
        img2.read(auxFn, DATA, n);

        for (size_t k = 1; k <= aDim.zdim; ++k)
        {
            img1.movePointerTo(k, n);
            img2.movePointerTo(k);

            img1().setXmippOrigin();
            img2().setXmippOrigin();

            EXPECT_TRUE(img1 == img2);
        }
    }
    XMIPP_CATCH
}

TEST_F( ImageTest, checkImageFileSize)
{
    XMIPP_TRY
    EXPECT_TRUE(checkImageFileSize("image/smallVolumeStack.stk"));
    EXPECT_FALSE(checkImageFileSize("image/smallVolumeStackCorrupted.stk"));
    XMIPP_CATCH
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
