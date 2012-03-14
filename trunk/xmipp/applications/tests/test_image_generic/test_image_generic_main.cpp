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
        //filename = getenv("XMIPP_HOME");
        //get example images/staks
        imageName = TEST_FILENAME("singleImage.spi");
        stackName = TEST_FILENAME("smallStack.stk");
        myImageGeneric.readMapped(imageName);
        myImageGeneric2.readMapped(stackName, 1);
    }

    // virtual void TearDown() {}//Destructor
    ImageGeneric myImageGeneric;
    ImageGeneric myImageGeneric2;
    Image<float> myImageFloat;
    FileName imageName;
    FileName stackName;
    // FileName filename;

};


TEST_F( ImageGenericTest, equalsOperator)
{
    ImageGeneric auxImageG;
    auxImageG.readMapped(imageName);
    ASSERT_TRUE(myImageGeneric==myImageGeneric);
    ASSERT_TRUE(myImageGeneric==auxImageG);
    ASSERT_FALSE(myImageGeneric==myImageGeneric2);
}

TEST_F( ImageGenericTest, equalsFunction)
{
    ImageGeneric auxImageG;
    auxImageG.read(imageName);
    ASSERT_TRUE (myImageGeneric.equal(auxImageG) );
    double f = auxImageG.getPixel(1,1);
    auxImageG.setPixel(1,1,f+XMIPP_EQUAL_ACCURACY/2.);
    ASSERT_TRUE (myImageGeneric.equal(auxImageG)      );
    auxImageG.setPixel(1,1,f+(XMIPP_EQUAL_ACCURACY));
    //std::cerr << *((MultidimArray<float>*)myImageGeneric.data->im) <<std::endl;
    //std::cerr << *((MultidimArray<float>*)auxImageG.data->im) <<std::endl;
    ASSERT_FALSE(myImageGeneric.equal(auxImageG));
    ASSERT_TRUE(myImageGeneric.equal(auxImageG,XMIPP_EQUAL_ACCURACY*4.));
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
    FileName auxFn = imageName.insertBeforeExtension("_swap");
    ImageGeneric auxImageGeneric;
    auxImageGeneric.read(auxFn);
    EXPECT_EQ(myImageGeneric,auxImageGeneric);
    auxImageGeneric.readMapped(auxFn);
    EXPECT_EQ(myImageGeneric,auxImageGeneric);
}

TEST_F( ImageGenericTest, add)
{
    FileName auxFilename1((String)"1@"+stackName);
    FileName auxFilename2((String)"2@"+stackName);
    ImageGeneric auxImageGeneric1(auxFilename1);
    ImageGeneric auxImageGeneric2(auxFilename2);
    auxImageGeneric1.add(auxImageGeneric2);
    auxFilename2  = TEST_FILENAME("sum.spi");
    auxImageGeneric2.read(auxFilename2);
    EXPECT_TRUE(auxImageGeneric1==auxImageGeneric2);
    auxImageGeneric1.add(auxImageGeneric2);
    EXPECT_FALSE(auxImageGeneric1==auxImageGeneric2);
}

TEST_F( ImageGenericTest, subtract)
{
    FileName sumFn = TEST_FILENAME("sum.spi");
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
    FileName Fn = TEST_FILENAME("emptyFile.stk");
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

/* check if an image declared as one kind of datatype is correctly
 * converted and/or casted to another datatype
 */
TEST_F( ImageGenericTest, convert2Datatype)
{
    ImageGeneric img;
    img.read(imageName);
    MultidimArray<float> * ma;
    MultidimArray<float> auxMa;
    MultidimArray<unsigned char> *uintMaP;
    MultidimArray<unsigned char> uintMa;
    img().getMultidimArrayPointer(ma);

    // Let's be sure any element is nonzero and higher than 255
    (*ma) += 1024;
    // Now copy it as reference and normalize
    auxMa = *ma;
    auxMa.rangeAdjust(0, 255);

    // Checking of both arrays are different
    EXPECT_NE(auxMa, *ma);

    // Change of datatype and conversion of values
    img.convert2Datatype(UChar, CW_CONVERT);
    typeCast(auxMa, uintMa);
    img().getMultidimArrayPointer(uintMaP);

    // Now both arrays are equal
    EXPECT_EQ(*uintMaP, uintMa);

    // Checking that convert2Datatype works after movePointer2
    ArrayDim befAdim, aftAdim;
    FileName auxFn = TEST_FILENAME("smallVolume.vol");
    myImageGeneric.readMapped(auxFn);
    myImageGeneric.getDimensions(befAdim);
    // Slices number when working with multidimarrays start in zero
    MULTIDIM_ARRAY_GENERIC(myImageGeneric).getSlice(2, auxMa);
    // Slices number when working with Images class start in FIRST_SLICE (1)
    myImageGeneric.movePointerTo(3);
    // Now It should keep the information of the slice only
    myImageGeneric.convert2Datatype(UChar);
    myImageGeneric.movePointerTo(ALL_SLICES);
    myImageGeneric.getDimensions(aftAdim);

    // Dimensions after convert2Datatype never will be the same
    EXPECT_NE(befAdim.zdim, aftAdim.zdim);

    // Now, lets check images are equals
    auxMa.rangeAdjust(0, 255);
    typeCast(auxMa, uintMa);
    MULTIDIM_ARRAY_GENERIC(myImageGeneric).getMultidimArrayPointer(uintMaP);

    EXPECT_EQ(*uintMaP, uintMa);
}

// check the reslicing is right
TEST_F( ImageGenericTest, reslice)
{
    FileName fnVol = TEST_FILENAME("progVol.vol");
    ImageGeneric imgSliced, imgRef;
    imgSliced.read(fnVol);
    imgRef.read(fnVol);

    MultidimArray<float> * dataS, *dataR;
    imgRef().getMultidimArrayPointer(dataR);

    imgSliced.reslice(ImageGeneric::TOP);
    imgSliced().getMultidimArrayPointer(dataS);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(*dataR)
    {
        EXPECT_EQ(DIRECT_ZYX_ELEM(*dataR,k,i,j), DIRECT_ZYX_ELEM(*dataS,ZSIZE(*dataS)-1-i,k,j));
    }

    imgSliced.reslice(ImageGeneric::BOTTOM);
    imgSliced().getMultidimArrayPointer(dataS);

    EXPECT_EQ(*dataR, *dataS);

    imgSliced.reslice(ImageGeneric::LEFT);
    imgSliced().getMultidimArrayPointer(dataS);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(*dataR)
    {
        EXPECT_EQ(DIRECT_ZYX_ELEM(*dataR,k,i,j), DIRECT_ZYX_ELEM(*dataS,ZSIZE(*dataS)-1-j,i,k));
    }
    imgSliced.reslice(ImageGeneric::RIGHT);
    imgSliced().getMultidimArrayPointer(dataS);

    EXPECT_EQ(*dataR, *dataS);
}

TEST_F( ImageGenericTest, getPreview)
{
    FileName auxFn = TEST_FILENAME("smallVolume.vol");
    ImageGeneric img1, img2, imgTemp;
    imgTemp.read(auxFn);

    imgTemp.getPreview(img1, 32,32, ALL_SLICES);
    img2.readPreview(auxFn, 32,32, ALL_SLICES);
    img1().setXmippOrigin();
    img2().setXmippOrigin();
    EXPECT_EQ(img1, img2) << "getPreview the whole volume.";

    for (int k = 1; k <= 4; ++k)
    {
        imgTemp.getPreview(img1, 32,32, k);
        img2.readPreview(auxFn, 32,32, k);
        img1().setXmippOrigin();
        img2().setXmippOrigin();
        EXPECT_EQ(img1, img2) << "getPreview a specific slice";
    }
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
