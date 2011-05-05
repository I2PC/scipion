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
    	//Mdjoin
        id = mDjoin.addObject();
        mDjoin.setValue(MDL_X,1.,id);
        mDjoin.setValue(MDL_Z,222.,id);
        id = mDjoin.addObject();
        mDjoin.setValue(MDL_X,3.,id);
        mDjoin.setValue(MDL_Z,444.,id);
        //mDanotherSource
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

    MetaData mDsource,mDanotherSource;
    MetaData mDunion, mDjoin;
    size_t id, id1,id2;
};

TEST_F( MetadataTest, similarToOperator)
{
    EXPECT_TRUE(mDsource==mDsource);
    EXPECT_FALSE(mDsource==mDanotherSource);
}
/** SORT FOR ROUTINE ALPHABETIC ORDER
 *
 */


TEST_F( MetadataTest, Addlabel)
{
	MetaData auxMetadata = mDunion;
    auxMetadata.addLabel(MDL_Z);
    std::vector<MDLabel> v1,v2;
    v1.push_back(MDL_X);
    v1.push_back(MDL_Y);
    v1.push_back(MDL_Z);
    v2 = auxMetadata.getActiveLabels();
    EXPECT_EQ(v2,v1);
}

TEST_F( MetadataTest, Clear)
{
	MetaData auxMetadata = mDsource;
    EXPECT_EQ(2,auxMetadata.size());
    auxMetadata.clear();
    EXPECT_EQ(0,auxMetadata.size());
}

TEST_F( MetadataTest, Copy)
{
	MetaData auxMetadata = mDsource;
    EXPECT_EQ(mDsource,auxMetadata);
}


TEST_F( MetadataTest, ImportObject)
{
	//FIXME importObjects test is in the test named select
	MetaData auxMetadata = mDsource;
    auxMetadata.importObject(mDunion,id1,false);
    auxMetadata.importObject(mDunion,id2,false);
    EXPECT_EQ(auxMetadata,mDunion);
}

TEST_F( MetadataTest, InnerJoin)
{
	MetaData auxMetadata;
	MetaData auxMetadata2 = mDsource;
    auxMetadata2.setValue(MDL_Z,222.,auxMetadata2.firstObject());
    auxMetadata2.setValue(MDL_Z,444.,auxMetadata2.firstObject()+1);//A little bit irregular
	auxMetadata.join(mDsource,mDjoin,MDL_X);
    EXPECT_EQ(auxMetadata,auxMetadata2)<< mDjoin;//print mDjoin if error
}

TEST_F( MetadataTest, Intersect)
{
	MetaData auxMetadata = mDunion;
    auxMetadata.intersection(mDsource,MDL_X);
    EXPECT_EQ(auxMetadata,mDsource);
}

TEST_F( MetadataTest, Merge)
{
	//FIXME is columns not in the same order equal to operator does not return OK
	//should not be like this
	MetaData auxMetadata3, auxMetadata,auxMetadata2;
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_Z,222.,id);
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_Z,444.,id);
	auxMetadata.join(mDsource,mDjoin,MDL_X);
    auxMetadata2 = mDsource;
	auxMetadata2.merge(auxMetadata3);
    EXPECT_EQ(auxMetadata,auxMetadata2);
}

TEST_F( MetadataTest, NaturalJoin)
{
	MetaData auxMetadata;
	MetaData auxMetadata3;
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_X,1.,id);
    auxMetadata3.setValue(MDL_Y,2.,id);
    auxMetadata3.setValue(MDL_Z,222.,id);
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_X,3.,id);
    auxMetadata3.setValue(MDL_Y,4.,id);
    auxMetadata3.setValue(MDL_Z,333.,id);
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_X,5.,id);
    auxMetadata3.setValue(MDL_Y,6.,id);
    auxMetadata3.setValue(MDL_Z,444.,id);

	auxMetadata.join(mDsource,auxMetadata3,MDL_X,NATURAL);
	auxMetadata3.removeObject(id);
    EXPECT_EQ(auxMetadata,auxMetadata3)<< auxMetadata3;//print mDjoin if error
}

TEST_F( MetadataTest, Operate)
{
	MetaData auxMetadata = mDunion;
	MetaData auxMetadata2 = mDunion;
    auxMetadata.operate((String)"X=2*X");
    double x;
    FOR_ALL_OBJECTS_IN_METADATA(auxMetadata2)
    {
    	auxMetadata2.getValue(MDL_X,x,__iter.objId);
    	auxMetadata2.setValue(MDL_X,x*2,__iter.objId);
    }

    EXPECT_EQ(auxMetadata,auxMetadata2);
}

TEST_F( MetadataTest, Randomize)
{
	int different,equal;
	different=-1;
	equal=-2;
	MetaData auxMetadata;
	for (int var = 0; var < 20; var++)
	{
	    auxMetadata.randomize(mDsource);
	    if(mDsource==auxMetadata)
	    	equal=1;
	    else
	    	different=1;
	}
    EXPECT_EQ(different,equal);
}

TEST_F( MetadataTest, Removelabel)
{
	MetaData auxMetadata = mDunion;
    auxMetadata.removeLabel(MDL_X);
    std::vector<MDLabel> v1,v2;
    v1.push_back(MDL_Y);
    v2 = auxMetadata.getActiveLabels();
    EXPECT_EQ(v2,v1);
}

TEST_F( MetadataTest, Size)
{
    EXPECT_EQ(2, mDsource.size());
}

TEST_F( MetadataTest, Sort)
{
	MetaData auxMetadata,auxMetadata2;
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_X,3.,id);
    auxMetadata.setValue(MDL_Y,4.,id);
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_X,1.,id);
    auxMetadata.setValue(MDL_Y,2.,id);
    auxMetadata2.sort(auxMetadata,MDL_X);
    EXPECT_EQ(auxMetadata2,mDsource);
}

TEST_F( MetadataTest, Substraction)
{
	MetaData auxMetadata = mDunion;
    auxMetadata.subtraction(mDanotherSource,MDL_X);
    EXPECT_EQ(auxMetadata,mDsource);
}

TEST_F( MetadataTest, Union)
{
	//FIXME union all is missing
	MetaData auxMetadata = mDsource;
    auxMetadata.unionAll(mDanotherSource);
    EXPECT_EQ(auxMetadata,mDunion);
}


//TEST_F( MetadataTest, select)
//{
//	MetaData auxMetadata;
//    auxMetadata.importObjects(mDsource,MDExpression((String)"x>-2"));
//    EXPECT_EQ(auxMetadata,mDsource);
//}

GTEST_API_ int main(int argc, char **argv)
{
    std::cout << "Running main() from gtest_main.cc\n";

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
