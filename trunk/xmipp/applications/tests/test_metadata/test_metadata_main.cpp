#include <data/metadata_extension.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"

/*
 * Define a "Fixure so we may reuse the metadatas
 */
class MetadataTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
    	//Md1
        id = mDsource.addObject();
        mDsource.setValue(MDL_X,1.,id);
        mDsource.setValue(MDL_Y,2.,id);
        id = mDsource.addObject();
        mDsource.setValue(MDL_X,3.,id);
        mDsource.setValue(MDL_Y,4.,id);
        //Md2
        id = mDanotherSource.addObject();
        mDanotherSource.setValue(MDL_X,11.,id);
        mDanotherSource.setValue(MDL_Y,22.,id);
        id = mDanotherSource.addObject();
        mDanotherSource.setValue(MDL_X,33.,id);
        mDanotherSource.setValue(MDL_Y,44.,id);

        //Md UnionAll
        mDunion = mDsource;
        id1 = mDunion.addObject();
        mDunion.setValue(MDL_X,11.,id1);
        mDunion.setValue(MDL_Y,22.,id1);
        id2 = mDunion.addObject();
        mDunion.setValue(MDL_X,33.,id2);
        mDunion.setValue(MDL_Y,44.,id2);
    }

    // virtual void TearDown() {}//Destructor

    MetaData mDsource,mDtarget,mDanotherSource;
    MetaData mDunion, auxMetadata;
    size_t id, id1,id2;
};


TEST_F( MetadataTest, Copy)
{
    mDtarget=mDsource;
    EXPECT_EQ(mDsource,mDtarget);
}

TEST_F( MetadataTest, Size)
{
    EXPECT_EQ(2, mDsource.size());
}

TEST_F( MetadataTest, Clear)
{
    auxMetadata = mDsource;
    EXPECT_EQ(2,auxMetadata.size());
    auxMetadata.clear();
    EXPECT_EQ(0,auxMetadata.size());
}

TEST_F( MetadataTest, importObjects)
{
	//FIXME importObjects is overloaded, only one case is tested
    auxMetadata = mDsource;
    auxMetadata.importObject(mDunion,id1,false);
    auxMetadata.importObject(mDunion,id2,false);
    EXPECT_EQ(auxMetadata,mDunion);
}

/*
 * Operations in sets
 */
TEST_F( MetadataTest, union)
{
	//FIXME union all is missing
    auxMetadata = mDsource;
    auxMetadata.unionAll(mDanotherSource);
    EXPECT_EQ(auxMetadata,mDunion);
}

TEST_F( MetadataTest, intersect)
{
    auxMetadata = mDunion;
    auxMetadata.intersection(mDsource,MDL_X);
    EXPECT_EQ(auxMetadata,mDsource);
}

TEST_F( MetadataTest, substraction)
{
    auxMetadata = mDunion;
    auxMetadata.subtraction(mDanotherSource,MDL_X);
    EXPECT_EQ(auxMetadata,mDsource);
}

GTEST_API_ int main(int argc, char **argv)
{
    std::cout << "Running main() from gtest_main.cc\n";

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
