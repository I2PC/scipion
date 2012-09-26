#include <data/xmipp_funcs.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.



/* I UPLOAD THIS SO XMIPP MAY NE COMPILE, KINO UPLOAD THE RIGHT FILE */ 
class Geometry : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
       //get example images/staks
    }

    FileName source1;
};


TEST_F( Geometry, CompareTwoFiles)
{
    ASSERT_TRUE(true);
}


GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
