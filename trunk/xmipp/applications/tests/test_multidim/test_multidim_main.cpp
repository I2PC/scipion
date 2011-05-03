#include <data/multidim_array.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
TEST( MultidimTest, Size)
{
    MultidimArray<int> md;
    md.resize(2,3);
    int x,y,z;
    size_t n;
    md.getDimensions(x,y,z,n);
    EXPECT_EQ(1, n) << "MultidimArray: wrong n size";
    EXPECT_EQ(1, z) << "MultidimArray: wrong y size";
    EXPECT_EQ(2, y) << "MultidimArray: wrong z size";
    EXPECT_EQ(3, x) << "MultidimArray: wrong x size";
}

TEST( MultidimTest, Copy)
{
    MultidimArray<int> mdtarget,mdsource;
    mdsource.resize(2,3);
    DIRECT_MULTIDIM_ELEM(mdsource,0)=1;
    mdtarget=mdsource;
    DIRECT_MULTIDIM_ELEM(mdsource,0)=1;
    EXPECT_EQ(mdtarget,mdsource) << "MultidimArray: copy operator failed";
}

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running main() from gtest_main.cc\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
