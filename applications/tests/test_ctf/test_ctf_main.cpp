#include <data/xmipp_image.h>
#include <reconstruction/transform_downsample.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
#include <data/ctf.h>

// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class CtfTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {

        try
        {
            //get example images/staks
            chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
            // testBaseName = xmippPath + "/resources/test";
            imageName = "image/singleImage.spi";

        }
        catch (XmippError &xe)
        {
            std::cerr << xe;
            exit(-1);
        }
    }

    // virtual void TearDown() {}//Destructor
    Image<double> myImageFloat;
    FileName imageName;


};

TEST_F( CtfTest, generateImageWithTwoCTFs)
{
    XMIPP_TRY
    MetaData metadata1;
    long objectId = metadata1.addObject();
    metadata1.setValue(MDL_CTF_SAMPLING_RATE, 1., objectId);
    metadata1.setValue(MDL_CTF_VOLTAGE, 300., objectId);
    metadata1.setValue(MDL_CTF_DEFOCUSU, 26000., objectId);
    metadata1.setValue(MDL_CTF_DEFOCUSV, 2000., objectId);
    metadata1.setValue(MDL_CTF_DEFOCUS_ANGLE, 45., objectId);
    metadata1.setValue(MDL_CTF_CS, 2., objectId);
    metadata1.setValue(MDL_CTF_Q0, 0.1, objectId);
    MultidimArray<double> in;
    generateCTFImageWith2CTFs(metadata1, metadata1, 256, in);
    Image<double> img;
    img().alias(in);
    //img.write("/tmp/kk.mrc");
    EXPECT_TRUE(true);
    XMIPP_CATCH
}

TEST_F( CtfTest, errorBetween2CTFs)
{
    XMIPP_TRY
    MetaData metadata1;
    long objectId = metadata1.addObject();
    metadata1.setValue(MDL_CTF_SAMPLING_RATE, 1., objectId);
    metadata1.setValue(MDL_CTF_VOLTAGE, 300., objectId);
    metadata1.setValue(MDL_CTF_DEFOCUSU, 5000., objectId);
    metadata1.setValue(MDL_CTF_DEFOCUSV, 7500., objectId);
    metadata1.setValue(MDL_CTF_DEFOCUS_ANGLE, -45., objectId);
    metadata1.setValue(MDL_CTF_CS, 2., objectId);
    metadata1.setValue(MDL_CTF_Q0, 0.1, objectId);

    MetaData metadata2;
    objectId = metadata2.addObject();
    metadata2.setValue(MDL_CTF_SAMPLING_RATE, 1., objectId);
    metadata2.setValue(MDL_CTF_VOLTAGE, 300., objectId);
    metadata2.setValue(MDL_CTF_DEFOCUSU, 10000., objectId);
    metadata2.setValue(MDL_CTF_DEFOCUSV, 10000., objectId);
    metadata2.setValue(MDL_CTF_DEFOCUS_ANGLE, 45., objectId);
    metadata2.setValue(MDL_CTF_CS, 2., objectId);
    metadata2.setValue(MDL_CTF_Q0, 0.1, objectId);

    double error = errorBetween2CTFs(metadata1,
                             metadata2,
                             256,
                             0.05,0.25);
    EXPECT_FLOAT_EQ(error,10441.1);
    XMIPP_CATCH
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

