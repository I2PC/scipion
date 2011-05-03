#include <data/metadata_extension.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"

// This test is named "Size", and belongs to the "MetadataTest"
// test case.
 TEST( MetadataTest, Size)
{
    MetaData md;
    md.addObject();
    EXPECT_EQ(1, md.size()) << "Metadata size provides wrong size";
}

 /*
  * These assertions can work with a user-defined type, but only if you
  * define the corresponding comparison operator (e.g. ==, <, etc).
  *  If the corresponding operator is defined, prefer using the ASSERT_*()
  *  macros because they will print out not only the result of the c
  *  omparison, but the two operands as well.
  */
 TEST( MetadataTest, Copy)
{
    MetaData mDsource,mDtarget;
    size_t id1 = mDsource.addObject();
    mDsource.setValue(MDL_X,1.,id1);
    mDsource.setValue(MDL_Y,2.,id1);
    size_t id2 = mDsource.addObject();
    mDsource.setValue(MDL_X,1.,id2);
    mDsource.setValue(MDL_Y,2.,id2);
    mDtarget=mDsource;
    ///////////////////////////////////DELETE
    //id1 = mDsource.addObject();
    EXPECT_EQ(mDsource,mDtarget);
}

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running main() from gtest_main.cc\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
