/***************************************************************************
 *
 * Authors:
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

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
    mDsource.setValue(MDL_Y,21122222.,id2);
    EXPECT_EQ(mDsource,mDtarget);
}

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running main() from gtest_main.cc\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
