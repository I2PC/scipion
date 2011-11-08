#include <data/metadata_extension.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
#include <stdlib.h>
#include <string.h>
/*
 * Define a "Fixure so we may reuse the metadatas
 */
class MetadataTest : public ::testing::Test
{
protected:
    //init metadatas
#define len 128

    virtual void SetUp()
    {
        //find binaries directory
        char szTmp[len];
        char pBuf[len];
        sprintf(szTmp, "/proc/%d/exe", getpid());
        int bytes = std::min(readlink(szTmp, pBuf, len), (ssize_t)len - 1);
        if(bytes >= 0)
            pBuf[bytes] = '\0';
        //remove last token
        FileName filename(pBuf);
        localDir = filename.removeFilename();

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
    FileName localDir;
};

TEST_F( MetadataTest, SimilarToOperator)
{
    ASSERT_EQ(mDsource,mDsource);
    ASSERT_FALSE(mDsource==mDanotherSource);
    //attribute order should not be important
    MetaData auxMetadata ;
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_Y,2.,id);
    auxMetadata.setValue(MDL_X,1.,id);
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_Y,4.,id);
    auxMetadata.setValue(MDL_X,3.,id);
    ASSERT_EQ(auxMetadata,mDsource);
    //Test form double with a given precission
    auxMetadata.clear();
    auxMetadata.setPrecission(2);
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_Y,2.001,id);
    auxMetadata.setValue(MDL_X,1.,id);
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_Y,4.,id);
    auxMetadata.setValue(MDL_X,3.,id);
    ASSERT_TRUE(auxMetadata==mDsource);
    auxMetadata.setPrecission(4);
    ASSERT_FALSE(auxMetadata==mDsource);


}
/** SORT FOR ROUTINE ALPHABETIC ORDER
 *
 */


TEST_F( MetadataTest, AddLabel)
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

TEST_F( MetadataTest, AddRow)
{
    MetaData md;

    MDRow row;
    row.setValue(MDL_X, 1.);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);
    row.setValue(MDL_X, 3.);
    row.setValue(MDL_Y, 4.);
    md.addRow(row);

    EXPECT_EQ(md, mDsource);
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

TEST_F( MetadataTest, ReadEmptyBlock)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testGetBlocks_XXXXXX", sizeof sfn);
    mkstemp(sfn);
    MetaData md;
    FileName fn = (String)"block_Empty@"+sfn;
    md.write(fn, MD_OVERWRITE);
    md.clear();
    md.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",md.addObject());
    md.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",md.addObject());
    md.write((String)"block_B1@"+sfn,MD_APPEND);

    EXPECT_NO_THROW(MetaData md2(fn););

    unlink(sfn);
}

TEST_F( MetadataTest, GetBlocksInMetadata)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testGetBlocks_XXXXXX", sizeof sfn);
    mkstemp(sfn);

    MetaData auxMetadata;
    auxMetadata.setValue(MDL_IMAGE,(String)"image_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_2.xmp",auxMetadata.addObject());
    auxMetadata.write(sfn,MD_OVERWRITE);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"block_000001@"+sfn,MD_APPEND);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"block_000002@"+sfn,MD_APPEND);
    auxMetadata.clear();

    StringVector compBlockList;
    compBlockList.push_back("");
    compBlockList.push_back("block_000001");
    compBlockList.push_back("block_000002");

    StringVector readBlockList;
    getBlocksInMetaDataFile(sfn,readBlockList);

    EXPECT_EQ(compBlockList,readBlockList);
    unlink(sfn);
}

TEST_F( MetadataTest, CheckRegularExpression)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testGetBlocks_XXXXXX", sizeof sfn);
    mkstemp(sfn);

    MetaData auxMd, auxMd2;
    auxMd.setValue(MDL_IMAGE,(String)"image_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_2.xmp",auxMd.addObject());
    auxMd.write(sfn,MD_OVERWRITE);
    auxMd.clear();
    auxMd.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",auxMd.addObject());
    auxMd.write((String)"block_000001@"+sfn,MD_APPEND);
    auxMd.clear();
    auxMd.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",auxMd.addObject());
    auxMd.write((String)"block_000002@"+sfn,MD_APPEND);
    auxMd.clear();
    auxMd.setValue(MDL_IMAGE,(String)"image_data_A_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_A_2.xmp",auxMd.addObject());
    auxMd.write((String)"block_A@"+sfn,MD_APPEND);
    auxMd.clear();

    auxMd.clear();
    auxMd.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",auxMd.addObject());

    auxMd2.read((String)"block_000[0-9][0-9][12]@" + sfn);
    EXPECT_EQ(auxMd, auxMd2);

    unlink(sfn);
}

TEST_F( MetadataTest, ImportObject)
{
    //FIXME importObjects test is in the test named select
    MetaData auxMetadata = mDsource;
    auxMetadata.importObject(mDunion,id1,false);
    auxMetadata.importObject(mDunion,id2,false);
    EXPECT_EQ(auxMetadata,mDunion);
}

TEST_F( MetadataTest, LeftJoin)
{
    MetaData auxMetadata;
    MetaData auxMetadata2 = mDsource;
    auxMetadata2.setValue(MDL_Z,222.,auxMetadata2.firstObject());
    auxMetadata2.setValue(MDL_Z,444.,auxMetadata2.firstObject()+1);//A little bit irregular
    auxMetadata.join(mDsource,mDjoin,MDL_X);
    EXPECT_EQ(auxMetadata,auxMetadata2)<< mDjoin;//print mDjoin if error
}

TEST_F( MetadataTest, InnerJoin1)
{
    MetaData auxMetadata;
    MetaData auxMetadataResult;
    MetaData auxMetadataLeft = mDsource;
    MetaData auxMetadataRight;

    auxMetadataRight.setValue(MDL_Z,1.,auxMetadataRight.firstObject());
    auxMetadataRight.setValue(MDL_ANGLEPSI,11.,auxMetadataRight.firstObject());

    auxMetadata.join(auxMetadataLeft,auxMetadataRight,MDL_X,MDL_Z,INNER);
    auxMetadataResult.setValue(MDL_X,1.,auxMetadataRight.firstObject());
    auxMetadataResult.setValue(MDL_Y,2.,auxMetadataRight.firstObject());
    auxMetadataResult.setValue(MDL_ANGLEPSI,11.,auxMetadataRight.firstObject());

    EXPECT_EQ(auxMetadata,auxMetadataResult)<< mDjoin;//print mDjoin if error
}

TEST_F( MetadataTest, InnerJoin2)
{
    MetaData auxMetadata;
    MetaData auxMetadataResult;
    MetaData auxMetadataLeft = mDsource;
    MetaData auxMetadataRight;

    auxMetadataRight.setValue(MDL_Z,1.,auxMetadataRight.firstObject());
    auxMetadataRight.setValue(MDL_Y,11.,auxMetadataRight.firstObject());

    auxMetadata.join(auxMetadataLeft,auxMetadataRight,MDL_X,MDL_Z,INNER);
    auxMetadataResult.setValue(MDL_X,1.,auxMetadataRight.firstObject());
    auxMetadataResult.setValue(MDL_Y,2.,auxMetadataRight.firstObject());

    EXPECT_EQ(auxMetadata,auxMetadataResult)<< mDjoin;//print mDjoin if error
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
    EXPECT_EQ(auxMetadata,auxMetadata3);
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

TEST_F( MetadataTest, ReadMultipleBlocks)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testReadMultipleBlocks_XXXXXX", sizeof sfn);
    mkstemp(sfn);

    MetaData auxMetadata;
    auxMetadata.setValue(MDL_IMAGE,(String)"image_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_2.xmp",auxMetadata.addObject());
    auxMetadata.write(sfn,MD_OVERWRITE);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"block_000001@"+sfn,MD_APPEND);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"block_000002@"+sfn,MD_APPEND);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_no_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_no_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"noblock@"+sfn,MD_APPEND);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_3_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_3_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"block_000003@"+sfn,MD_APPEND);
    auxMetadata.clear();

    MetaData compMetadata;
    compMetadata.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",compMetadata.addObject());
    compMetadata.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",compMetadata.addObject());
    compMetadata.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",compMetadata.addObject());
    compMetadata.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",compMetadata.addObject());
    auxMetadata.read((String)"block_00000[12]@"+sfn);
    EXPECT_EQ(compMetadata,auxMetadata);
    compMetadata.clear();

    compMetadata.setValue(MDL_IMAGE,(String)"image_data_3_1.xmp",compMetadata.addObject());
    compMetadata.setValue(MDL_IMAGE,(String)"image_data_3_2.xmp",compMetadata.addObject());
    auxMetadata.read((String)"block_000003@"+sfn);
    EXPECT_EQ(compMetadata,auxMetadata);
    compMetadata.clear();

    compMetadata.setValue(MDL_IMAGE,(String)"image_1.xmp",compMetadata.addObject());
    compMetadata.setValue(MDL_IMAGE,(String)"image_2.xmp",compMetadata.addObject());
    auxMetadata.read(sfn);
    EXPECT_EQ(compMetadata,auxMetadata);
    compMetadata.clear();
    unlink(sfn);
}

TEST_F( MetadataTest, ReadEmptyBlocks)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testReadMultipleBlocks_XXXXXX", sizeof sfn);
    mkstemp(sfn);

    MetaData auxMetadata;
    id=auxMetadata.addObject();
    auxMetadata.setValue(MDL_X,1.,id);
    auxMetadata.setValue(MDL_Y,2.,id);
    auxMetadata.setValue(MDL_Z,222.,id);
    auxMetadata.write((String)"block_000001@"+sfn,MD_APPEND);

    auxMetadata.clear();
    auxMetadata.addLabel(MDL_X);
    auxMetadata.addLabel(MDL_Y);
    auxMetadata.addLabel(MDL_Z);
    auxMetadata.write((String)"block_000002@"+sfn,MD_APPEND);

    auxMetadata.clear();
    id=auxMetadata.addObject();
    auxMetadata.setValue(MDL_X,1.,id);
    auxMetadata.setValue(MDL_Y,2.,id);
    auxMetadata.setValue(MDL_Z,222.,id);
    auxMetadata.write((String)"block_000003@"+sfn,MD_APPEND);

    auxMetadata.clear();
    auxMetadata.addLabel(MDL_X);
    auxMetadata.addLabel(MDL_Y);
    auxMetadata.addLabel(MDL_Z);
    auxMetadata.write((String)"block_000004@"+sfn,MD_APPEND);

    auxMetadata.read((String)"block_000002@"+sfn);
    EXPECT_EQ(auxMetadata.size(),0);

    auxMetadata.read((String)"block_000004@"+sfn);
    EXPECT_EQ(auxMetadata.size(),0);

    unlink(sfn);
}

TEST_F( MetadataTest, ReadEmptyBlocksII)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testReadMultipleBlocks_XXXXXX", sizeof sfn);
    mkstemp(sfn);

    MetaData auxMetadata;

    auxMetadata.addLabel(MDL_X);
    auxMetadata.addLabel(MDL_Y);
    auxMetadata.addLabel(MDL_Z);
    auxMetadata.write((String)"block_000002@"+sfn,MD_APPEND);

    auxMetadata.read((String)"block_000002@"+sfn);
    EXPECT_EQ(auxMetadata.size(),0);
    unlink(sfn);
}

TEST_F( MetadataTest, ReadWrite)
{
    //temp file name
    char sfn[32] = "";
    strncpy(sfn, "/tmp/testWrite_XXXXXX", sizeof sfn);
    mkstemp(sfn);
    mDsource.write(sfn);

    MetaData auxMetadata;
    auxMetadata.read(sfn);

    EXPECT_EQ(mDsource,auxMetadata);
    unlink(sfn);
}

TEST_F( MetadataTest, ExistsBlock)
{
    //temp file name
    char sfn[32] = "";
    strncpy(sfn, "/tmp/testWrite_XXXXXX", sizeof sfn);
    mkstemp(sfn);
    FileName tmpFileName((String) "kk@" + sfn);
    mDsource.write(tmpFileName);
    MetaData auxMetadata;
    bool result1 = auxMetadata.existsBlock(tmpFileName);
    EXPECT_EQ(result1,true);
    tmpFileName =(String) "kk2@" + sfn;
    result1 = auxMetadata.existsBlock(tmpFileName);
    EXPECT_EQ(result1,false);
    unlink(sfn);
}

TEST_F( MetadataTest, ReadWriteAppendBlock)
{
    //temp file name
    char sfn[32] = "";
    strncpy(sfn, "/tmp/testWrite_XXXXXX", sizeof sfn);
    mkstemp(sfn);
    mDsource.write((String)"one@"+sfn);
    mDsource.write((String)"two@"+sfn,MD_APPEND);
    mDsource.write((String)"three@"+sfn,MD_APPEND);
    MetaData auxMetadata;
    FileName sfn2(localDir+ "/../applications/tests/test_metadata/ReadWriteAppendBlock.xmd");
    auxMetadata.read(sfn);
    EXPECT_TRUE(compareTwoFiles(sfn,sfn2,0));
    unlink(sfn);
}

TEST_F( MetadataTest, RemoveDuplicates)
{
    MetaData auxMetadata1,auxMetadata3;
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_X,1.,id);
    auxMetadata3.setValue(MDL_Y,2.,id);
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_X,3.,id);
    auxMetadata3.setValue(MDL_Y,4.,id);
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_X,1.,id);
    auxMetadata3.setValue(MDL_Y,2.,id);
    auxMetadata1.removeDuplicates(auxMetadata3);
    EXPECT_EQ(auxMetadata1,mDsource);//print mDjoin if error
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

TEST_F( MetadataTest, Select)
{
    MetaData auxMetadata;
    MetaData auxMetadata2;
    id = auxMetadata2.addObject();
    auxMetadata2.setValue(MDL_X,3.,id);
    auxMetadata2.setValue(MDL_Y,4.,id);

    auxMetadata.importObjects(mDsource,MDExpression((String)"x>2"));
    EXPECT_EQ(auxMetadata,auxMetadata2);
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

//check if mdl label match its type and
//check if int is different from size_t
TEST_F( MetadataTest, setGetValue)
{
    size_t t;
    int i;
    EXPECT_EQ(MDL::labelType(MDL_ORDER),LABEL_LONG);
    MetaData auxMetadata;
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_ORDER,(size_t)1, id);
    auxMetadata.getValue(MDL_ORDER,t, id);
    EXPECT_EQ((size_t)1,t);
    //We expect that MetaData will throw an exception
    //if you use getValue with a variable of type that
    // doesn't match the label type
    EXPECT_THROW(auxMetadata.getValue(MDL_ORDER, i, id), XmippError);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
