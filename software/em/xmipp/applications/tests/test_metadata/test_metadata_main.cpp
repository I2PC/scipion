#include <unistd.h> 
#include <stdlib.h>
#include <stdio.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image_convert.h>
#include <data/xmipp_funcs.h>
#include <iostream>
#include <gtest/gtest.h>
#include <string.h>
#include <fstream>
#include <sys/time.h>

#define N_ROWS_TEST		2
#define N_ROWS_PERFORMANCE_TEST		8000


/*
 * Define a "Fixture so we may reuse the metadatas
 */
class MetadataTest : public ::testing::Test
{
protected:
    //init metadatas

    virtual void SetUp()
    {
        if (chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str())==-1)
        	REPORT_ERROR(ERR_UNCLASSIFIED,"Could not change directory");
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

TEST_F( MetadataTest, AddIndex)
{
    MetaData auxMetadata = mDunion;
    auxMetadata.addIndex(MDL_X);
    EXPECT_EQ(1,1);
}


TEST_F( MetadataTest, AddRow)
{
	MetaData md, md2;

    MDRow row;
    row.setValue(MDL_X, 1.);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);
    row.setValue(MDL_X, 3.);
    row.setValue(MDL_Y, 4.);
    md.addRow(row);

    row.setValue(MDL_X, 1.);
    row.setValue(MDL_Y, 2.);
    md2.addRow2(row);
    row.setValue(MDL_X, 3.);
    row.setValue(MDL_Y, 4.);
    md2.addRow2(row);

    EXPECT_EQ(md, mDsource);
    EXPECT_EQ(md2, mDsource);
}

TEST_F( MetadataTest, AddRows)
{
    int 	i;					// Loop counter.
    bool	inserted;			// Insertion return value.
    MetaData md;
    MDRow row[N_ROWS_TEST];		// Rows array.

    // Initialize rows.
    row[0].setValue(MDL_X, 1.);
    row[0].setValue(MDL_Y, 2.);
    row[1].setValue(MDL_X, 3.);
    row[1].setValue(MDL_Y, 4.);

    // Initialize insertion.
    if (md.initAddRow( row[0]))
    {
        // Add rows loop.
    	i=0;
    	do
    	{
    		// Insert row and increase number of insertions.
    		inserted = md.execAddRow(row[i]);
        	i++;
    	}
    	while ((i<N_ROWS_TEST) && (inserted));

        // Finalize statement.
        md.finalizeAddRow();
    }

    // Check result.
    EXPECT_EQ(md, mDsource);
}

TEST_F( MetadataTest, AddRowsPerformance)
{
    int 	i,j;				// Loop counters.
    bool	inserted;			// Insertion return value.
    MetaData md, md2, md3;
    MDRow 	row;	// Sample row

    printf("N_ROWS_PERFORMANCE_TEST = %d\n", N_ROWS_PERFORMANCE_TEST);

    // Initialize row.
    row.setValue(MDL_X,1.);
    row.setValue(MDL_Y,2.);
    row.setValue(MDL_ZSCORE,3.);
    row.setValue(MDL_ZSCORE_HISTOGRAM,4.);
    row.setValue(MDL_ZSCORE_RESMEAN,5.);
    row.setValue(MDL_ZSCORE_RESVAR,6.);
    row.setValue(MDL_ZSCORE_RESCOV,7.);
    row.setValue(MDL_ZSCORE_SHAPE1,8.);
    row.setValue(MDL_ZSCORE_SHAPE2,9.);
    row.setValue(MDL_ZSCORE_SNR1,10.);
    row.setValue(MDL_ZSCORE_SNR2,11.);
    row.setValue(MDL_IMAGE, String("particles.stk"));
    row.setValue(MDL_SHIFT_X_DIFF, 1.5);
    row.setValue(MDL_SHIFT_Y_DIFF, 2.5);
    row.setValue(MDL_CONTINUOUS_X, 1.5);
    row.setValue(MDL_CONTINUOUS_Y, 2.5);
    row.setValue(MDL_SHIFT_X, 1.5);
    row.setValue(MDL_SHIFT_Y, 2.5);
    row.setValue(MDL_SHIFT_Z, 3.5);

    Timer t;
    size_t s1, s2, s3;
    float r;


    t.tic();

    for (i=0; i<N_ROWS_PERFORMANCE_TEST; i++)
    {
    	md.addRow(row);
    }
    s1 = t.toc("Time original:", false);

    t.tic();
    for (i=0; i<N_ROWS_PERFORMANCE_TEST; i++)
    {
    	md2.addRow2(row);
    }

    s2 = t.toc("Time by row: ", false);
    printf("    Speed up from original: %f\n", ((float) s1 / (float) s2));

    // Initialize insertion.
    t.tic();
    if (md3.initAddRow(row))
    {
        // Add rows loop.
    	i=0;
    	do
    	{
    		// Insert row and increase number of insertions.
    		md3.execAddRow(row);
        	i++;
    	}
    	while (i<N_ROWS_PERFORMANCE_TEST);

        // Finalize statement.
    	md3.finalizeAddRow();
    }
    s3 = t.toc("Time by set:", false);
    printf("    Speed up from original: %f\n", ((float) s1 / (float) s3));
    printf("    Speed up from row: %f\n", ((float) s2 / (float) s3));
    // Check result.
    EXPECT_EQ(md, md2);
    EXPECT_EQ(md2, md3);
}

TEST_F( MetadataTest, addLabelAlias)
{
    //metada with no xmipp labels
    FileName fnNonXmippSTAR =(String)"metadata/noXmipp.xmd";
    MDL::addLabelAlias(MDL_Y,(String)"noExixtingLabel");
    MetaData md = MetaData(fnNonXmippSTAR);
    EXPECT_EQ(mDsource,md);
}

TEST_F( MetadataTest, Aggregate1)
{
    //simple aggregation
    MetaData md,mdOut;
    size_t count;

    MDRow row;
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_Y, 2.);
    row.setValue(MDL_DEFGROUP, 2);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_Y, 4.);
    row.setValue(MDL_DEFGROUP, 23);
    md.addRow(row);

    mdOut.aggregate(md, AGGR_COUNT, MDL_ORDER, MDL_ORDER, MDL_COUNT);
    mdOut.getValue(MDL_COUNT,count,mdOut.firstObject());
    EXPECT_EQ(count, (size_t)2);
    mdOut.clear();
    mdOut.aggregate(md, AGGR_COUNT, MDL_Y, MDL_Y, MDL_COUNT);
    mdOut.getValue(MDL_COUNT,count,mdOut.firstObject());
    EXPECT_EQ(count,(size_t)1);

    MDObject mdValueOut(MDL_Y);
    double d;
    md.aggregateSingle(mdValueOut, AGGR_MAX ,MDL_Y);
    mdValueOut.getValue(d);
    EXPECT_EQ(d,4);

    MDObject mdValueOut2(MDL_ORDER);
    size_t t;
    md.aggregateSingleSizeT(mdValueOut2, AGGR_MAX ,MDL_ORDER);
    mdValueOut2.getValue(t);
    EXPECT_EQ(t,(size_t)1);

    MDObject mdValueOut3(MDL_DEFGROUP);
    int i;
    md.aggregateSingleInt(mdValueOut3, AGGR_MAX ,MDL_DEFGROUP);
    mdValueOut3.getValue(i);
    EXPECT_EQ(i,(int)23);

}

TEST_F( MetadataTest, Aggregate2)
{
    //multiple aggregarion
    MetaData md,mdOut;

    MDRow row;
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_Y, 4.);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)2);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);

    const AggregateOperation MyaggregateOperations[] =
        {
            AGGR_COUNT, AGGR_SUM, AGGR_MIN, AGGR_MAX, AGGR_AVG
        };
    std::vector<AggregateOperation> aggregateOperations(MyaggregateOperations,MyaggregateOperations+5);

    const MDLabel MyoperateLabels[]       =
        {
            MDL_ORDER,MDL_ORDER, MDL_Y, MDL_Y, MDL_Y
        };
    std::vector<MDLabel> operateLabels(MyoperateLabels,MyoperateLabels+5);

    const MDLabel MyresultLabels[]        =
        {
            MDL_ORDER,MDL_COUNT, MDL_SUM,  MDL_MIN, MDL_MAX, MDL_AVG
        };
    std::vector<MDLabel> resultLabels(MyresultLabels,MyresultLabels+6);

    mdOut.aggregate(md,aggregateOperations,operateLabels,resultLabels);
    md.clear();
    row.clear();
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_COUNT, (size_t)2);
    row.setValue(MDL_SUM, 2.);
    row.setValue(MDL_MIN, 2.);
    row.setValue(MDL_MAX, 4.);
    row.setValue(MDL_AVG, 3.);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)2);
    row.setValue(MDL_COUNT, (size_t)1);
    row.setValue(MDL_SUM, 2.);
    row.setValue(MDL_MIN, 2.);
    row.setValue(MDL_MAX, 2.);
    row.setValue(MDL_AVG, 2.);
    md.addRow(row);


    EXPECT_EQ(mdOut,md);
}

TEST_F( MetadataTest, AggregateGroupBy)
{
    //aggregation simple grouped by several attributes
    MetaData md,mdOut;

    MDRow row;
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_DEFGROUP, 2);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_DEFGROUP, 2);
    row.setValue(MDL_Y, 4.);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)2);
    row.setValue(MDL_DEFGROUP, 2);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);

    const MDLabel myGroupByLabels[] =
        {
            MDL_ORDER, MDL_DEFGROUP
        };
    std::vector<MDLabel> groupbyLabels(myGroupByLabels,myGroupByLabels+2);
    mdOut.aggregateGroupBy(md, AGGR_COUNT, groupbyLabels, MDL_Y, MDL_COUNT);

    md.clear();
    row.clear();
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_DEFGROUP, 2);
    row.setValue(MDL_COUNT, (size_t)2);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)2);
    row.setValue(MDL_DEFGROUP, 2);
    row.setValue(MDL_COUNT, (size_t)1);
    md.addRow(row);

    EXPECT_EQ(mdOut,md);
}

TEST_F( MetadataTest, Clear)
{
    MetaData auxMetadata = mDsource;
    EXPECT_EQ((size_t)2,auxMetadata.size());
    auxMetadata.clear();
    EXPECT_EQ((size_t)0,auxMetadata.size());
}

TEST_F( MetadataTest, Copy)
{
    MetaData auxMetadata = mDsource;
    EXPECT_EQ(mDsource,auxMetadata);
}

TEST_F( MetadataTest, MDInfo)
{
    //char sfnStar[64] = "";
    //char sfnSqlite[64] = "";
    //strncpy(sfnStar, "MDInfo_XXXXXX.xmd", sizeof "MDInfo_XXXXXX.xmd");
    //strncpy(sfnSqlite, "MDInfo_XXXXXX.sqlite", sizeof "MDInfo_XXXXXX.sqlite");
    FileName fn;
    fn.initUniqueName("MDInfo_XXXXXX");
    FileName fnDB;
    fnDB = fn + ".sqlite";
    FileName fnSTAR;
    fnSTAR = fn + ".xmd";

    XMIPP_TRY

    mDsource.write(fnDB);
    mDsource.write(fnSTAR);

    MetaData md;
    //Read from sqlite
    md.read(fnDB);
    MetaData mdOnlyOne;
    mdOnlyOne.setMaxRows(1);
    mdOnlyOne.read(fnDB);
    EXPECT_EQ(md.size(), mdOnlyOne.getParsedLines());

    mdOnlyOne.read(fnSTAR);
    EXPECT_EQ(md.size(), mdOnlyOne.getParsedLines());

    // Read from STAR
    md.read(fnSTAR);
    MDLabelVector labels = md.getActiveLabels();
    // Check containsLabel is true for all md labels
    for (size_t i = 0; i < labels.size(); ++i)
      EXPECT_TRUE(mdOnlyOne.containsLabel(labels[i]));

    XMIPP_CATCH

    unlink(fn.c_str());
    unlink(fnDB.c_str());
    unlink(fnSTAR.c_str());
}

TEST_F( MetadataTest,multiWrite)
{
/*
	char sfnStar[64] = "";
	char sfnSqlite[64] = "";
	char sfnXML[64] = "";
    strncpy(sfnStar, "/tmp/testReadMultipleBlocks_XXXXXX.xmd", sizeof sfnStar);
    if (mkstemps(sfnStar,4)==-1)
	    REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary STAR file");
    strncpy(sfnSqlite, "/tmp/testReadMultipleBlocks_XXXXXX.sqlite", sizeof sfnSqlite);
    if (mkstemps(sfnSqlite,7)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary SQLITE file");
    strncpy(sfnXML, "/tmp/testReadMultipleBlocks_XXXXXX.xml", sizeof sfnXML);
    if (mkstemps(sfnXML,4)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary XML file");
*/

    FileName fn   ;
    FileName fnDB   ;
    FileName fnXML  ;
    FileName fnSTAR ;
    fn.initUniqueName("/tmp/testReadMultipleBlocks_XXXXXX");
    fnDB = fn + ".sqlite";
    fnXML = fn + ".xml";
    fnSTAR = fn + ".xmd";
    
    FileName fnDBref   =(String)"metadata/mDsource.sqlite";
    FileName fnXMLref  =(String)"metadata/mDsource.xml";
    FileName fnSTARref =(String)"metadata/mDsource.xmd";

    XMIPP_TRY
    mDsource.write((String)"myblock@"+fnDB);
    mDsource.write((String)"myblock@"+fnXML);
    mDsource.write((String)"myblock@"+fnSTAR);
    XMIPP_CATCH

    EXPECT_TRUE(compareTwoFiles(fnDB, fnDBref));
    EXPECT_TRUE(compareTwoFiles(fnXML, fnXMLref));
    EXPECT_TRUE(compareTwoFiles(fnSTAR, fnSTARref));

    unlink(fn.c_str());
    unlink(fnDB.c_str());
    unlink(fnXML.c_str());
    unlink(fnSTAR.c_str());
}

TEST_F( MetadataTest,multiWriteSqlite)
{
    //char sfn[64] = "";
    //strncpy(sfn, "/tmp/multiWriteSqlite_XXXXXX.sqlite", sizeof sfn);
    //if (mkstemps(sfn,7)==-1)
    // 	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
    FileName fn   ;
    FileName fnDB   ;
    fn.initUniqueName("/tmp/multiWriteSqlite_XXXXXX");
    fnDB = fn + ".sqlite";
    FileName fnDBref   =(String)"metadata/mDsource.sqlite";

    MetaData md = MetaData();
    MetaData mdRead = MetaData();
    MDRow row;

    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_DEFGROUP, 2);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)1);
    row.setValue(MDL_DEFGROUP, 2);
    row.setValue(MDL_Y, 4.);
    md.addRow(row);
    row.setValue(MDL_ORDER, (size_t)2);
    row.setValue(MDL_DEFGROUP, 2);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);

    FileName fileNameA = FileName();
    XMIPP_TRY
    fileNameA.compose("block001", fnDB);

    md.setValue(MDL_ORDER,(size_t)11, md.firstObject());
    md.write(fileNameA);
    mdRead.read(fileNameA);
    EXPECT_EQ(md,mdRead);

    md.setValue(MDL_ORDER,(size_t)22, md.firstObject());
    fileNameA.compose("block002", fnDB);
    md.write(fileNameA,MD_APPEND);
    mdRead.read(fileNameA);
    EXPECT_EQ(md,mdRead);

    md.setValue(MDL_ORDER,(size_t)33, md.firstObject());
    fileNameA.compose("block003", fnDB);
    md.write(fileNameA,MD_APPEND);
    mdRead.read(fileNameA);
    EXPECT_EQ(md,mdRead);

    md.setValue(MDL_ORDER,(size_t)44, md.firstObject());
    fileNameA.compose("block003", fnDB);
    md.write(fileNameA,MD_APPEND);
    mdRead.read(fileNameA);
    EXPECT_EQ(md,mdRead);

    StringVector readBlockList;
    getBlocksInMetaDataFile(fnDB,readBlockList);
    EXPECT_EQ(readBlockList[0],"block001");
    EXPECT_EQ(readBlockList[1],"block002");
    EXPECT_EQ(readBlockList[2],"block003");

    XMIPP_CATCH
    unlink(fn.c_str());
    unlink(fnDB.c_str());
}

TEST_F( MetadataTest, ReadEmptyBlock)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testGetBlocks_XXXXXX", sizeof sfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
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
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

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
    compBlockList.push_back(DEFAULT_BLOCK_NAME);
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
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

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
    auxMd.setValue(MDL_IMAGE,(String)"image_data_3_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_3_2.xmp",auxMd.addObject());
    auxMd.write((String)"block_000003@"+sfn,MD_APPEND);
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
    auxMd.setValue(MDL_IMAGE,(String)"image_data_3_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_3_2.xmp",auxMd.addObject());

    auxMd2.read((String)"block_000[0-9][0-9][123]@" + sfn);
    EXPECT_EQ(auxMd, auxMd2);

    unlink(sfn);
}

TEST_F( MetadataTest, CheckRegularExpression2)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testGetBlocks_XXXXXX", sizeof sfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

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
    auxMd.write((String)"block_0000023@"+sfn,MD_APPEND);
    auxMd.clear();

    auxMd.clear();
    auxMd.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",auxMd.addObject());

    auxMd2.read((String)"block_000[0-9][0-9][0-9]$@" + sfn);
    EXPECT_EQ(auxMd, auxMd2);

    unlink(sfn);
}

TEST_F( MetadataTest, compareTwoMetadataFiles)
{
    XMIPP_TRY
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testGetBlocks_XXXXXX", sizeof sfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
    char sfn2[64] = "";
    strncpy(sfn2, "/tmp/testGetBlocks_XXXXXX", sizeof sfn2);
    if (mkstemp(sfn2)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
    char sfn3[64] = "";
    strncpy(sfn3, "/tmp/testGetBlocks_XXXXXX", sizeof sfn3);
    if (mkstemp(sfn3)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

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
    auxMd.write(sfn2, MD_OVERWRITE);
    auxMd.clear();
    auxMd.setValue(MDL_IMAGE,(String)"image_data_A_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_A_2.xmp",auxMd.addObject());
    auxMd.write((String)"block_000001@"+sfn2,MD_APPEND);
    auxMd.clear();

    EXPECT_FALSE(compareTwoMetadataFiles(sfn, sfn2));
    EXPECT_TRUE(compareTwoMetadataFiles(sfn, sfn));

    auxMd.setValue(MDL_IMAGE,(String)"image_1.xmpSPACE",auxMd.addObject());//extra space
    auxMd.setValue(MDL_IMAGE,(String)"image_2.xmp",auxMd.addObject());
    auxMd.write(sfn2,MD_OVERWRITE);
    auxMd.clear();
    auxMd.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",auxMd.addObject());
    auxMd.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",auxMd.addObject());
    auxMd.write((String)"block_000001@"+sfn2,MD_APPEND);

    String command=(String)"sed 's/SPACE/ /g' " + sfn2 + (String) ">" + sfn3;
    if (system (command.c_str())==-1)
    	REPORT_ERROR(ERR_UNCLASSIFIED,"Could not open shell");

    EXPECT_TRUE(compareTwoMetadataFiles(sfn, sfn3));

    unlink(sfn);
    unlink(sfn2);
    unlink(sfn3);
    XMIPP_CATCH
}

TEST_F( MetadataTest, FillExpand)
{
    XMIPP_TRY
	//create 2 temporary CTFs plus a metadata with dependences
    char sfn1[64] = "";
    strncpy(sfn1, "/tmp/FillExpandCTF2_XXXXXX", sizeof sfn1);
    if (mkstemp(sfn1)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
    char sfn2[64] = "";
    strncpy(sfn2, "/tmp/FillExpandMD_XXXXXX", sizeof sfn2);
    if (mkstemp(sfn2)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

    //create 2 CTFs
    MetaData ctfMd1, ctfMd2, md;
    ctfMd1.setColumnFormat(false);
    ctfMd2.setColumnFormat(false);
    size_t id1 =ctfMd1.addObject();
    size_t id2 =ctfMd2.addObject();
    ctfMd1.setValue(MDL_CTF_SAMPLING_RATE,1.,id1);
    ctfMd2.setValue(MDL_CTF_SAMPLING_RATE,1.,id2);
    ctfMd1.setValue(MDL_CTF_VOLTAGE,100.,id1);
    ctfMd2.setValue(MDL_CTF_VOLTAGE,100.,id2);
    ctfMd1.setValue(MDL_CTF_DEFOCUSU,1000.,id1);
    ctfMd2.setValue(MDL_CTF_DEFOCUSU,1500.,id2);
    ctfMd1.write(sfn1,MD_OVERWRITE);
    ctfMd2.write(sfn2,MD_OVERWRITE);
    //create 1 md referring the ctf
    size_t id =md.addObject();
    md.setValue(MDL_IMAGE,(String)"image1",id);
    md.setValue(MDL_CTF_MODEL,(String)sfn1,id);
    id =md.addObject();
    md.setValue(MDL_IMAGE,(String)"image2",id);
    md.setValue(MDL_CTF_MODEL,(String)sfn1,id);
    id =md.addObject();
    md.setValue(MDL_IMAGE,(String)"image3",id);
    md.setValue(MDL_CTF_MODEL,(String)sfn2,id);
    // call fillExpand
    md.fillExpand(MDL_CTF_MODEL);
    //create md with results
    MetaData mdResults;
    id =mdResults.addObject();
    mdResults.setValue(MDL_IMAGE,(String)"image1",id);
    mdResults.setValue(MDL_CTF_MODEL,(String)sfn1,id);
    mdResults.setValue(MDL_CTF_SAMPLING_RATE,1.,id);
    mdResults.setValue(MDL_CTF_VOLTAGE,100.,id);
    mdResults.setValue(MDL_CTF_DEFOCUSU,1000.,id);
    id =mdResults.addObject();
    mdResults.setValue(MDL_IMAGE,(String)"image2",id);
    mdResults.setValue(MDL_CTF_MODEL,(String)sfn1,id);
    mdResults.setValue(MDL_CTF_SAMPLING_RATE,1.,id);
    mdResults.setValue(MDL_CTF_VOLTAGE,100.,id);
    mdResults.setValue(MDL_CTF_DEFOCUSU,1000.,id);
    id =mdResults.addObject();
    mdResults.setValue(MDL_IMAGE,(String)"image3",id);
    mdResults.setValue(MDL_CTF_MODEL,(String)sfn2,id);
    mdResults.setValue(MDL_CTF_SAMPLING_RATE,1.,id);
    mdResults.setValue(MDL_CTF_VOLTAGE,100.,id);
    mdResults.setValue(MDL_CTF_DEFOCUSU,1500.,id);
    EXPECT_EQ(md,mdResults);
    //mdResults.setValue(MDL_CTF_DEFOCUSU,15000.,id);
    //EXPECT_NE(md,mdResults);

    unlink(sfn1);
    unlink(sfn2);
    XMIPP_CATCH
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
    auxMetadata.join1(mDsource, mDjoin, MDL_X);
    EXPECT_EQ(auxMetadata,auxMetadata2)<< mDjoin;//print mDjoin if error
}

TEST_F( MetadataTest, InnerJoin1)
{
    MetaData auxMetadata;
    MetaData auxMetadataResult;
    MetaData auxMetadataLeft = mDsource;
    MetaData auxMetadataRight;

    auxMetadataRight.setValue(MDL_Z,1.,auxMetadataRight.firstObject());
    auxMetadataRight.setValue(MDL_ANGLE_PSI,11.,auxMetadataRight.firstObject());

    auxMetadata.join2(auxMetadataLeft,auxMetadataRight,MDL_X,MDL_Z,INNER);
    auxMetadataResult.setValue(MDL_X,1.,auxMetadataRight.firstObject());
    auxMetadataResult.setValue(MDL_Y,2.,auxMetadataRight.firstObject());
    auxMetadataResult.setValue(MDL_ANGLE_PSI,11.,auxMetadataRight.firstObject());

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

    auxMetadata.join2(auxMetadataLeft,auxMetadataRight,MDL_X,MDL_Z,INNER);
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
    auxMetadata.join1(mDsource,mDjoin,MDL_X);
    auxMetadata2 = mDsource;
    auxMetadata2.merge(auxMetadata3);
    EXPECT_EQ(auxMetadata,auxMetadata2);
}


TEST_F( MetadataTest, MultiQuery)
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
    auxMetadata3.setValue(MDL_X,3.,id);
    auxMetadata3.setValue(MDL_Y,4.,id);
    auxMetadata3.setValue(MDL_Z,444.,id);

    MDValueEQ eq1(MDL_X, 3.);
    MDValueEQ eq2(MDL_Y, 4.);
    MDMultiQuery multi;

    //Test empty query
    auxMetadata.importObjects(auxMetadata3, multi);
    EXPECT_EQ(auxMetadata3,auxMetadata);

    multi.addAndQuery(eq1);
    multi.addAndQuery(eq2);

    auxMetadata.importObjects(auxMetadata3, multi);

    MetaData outMetadata;
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,3.,id);
    outMetadata.setValue(MDL_Y,4.,id);
    outMetadata.setValue(MDL_Z,333.,id);
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,3.,id);
    outMetadata.setValue(MDL_Y,4.,id);
    outMetadata.setValue(MDL_Z,444.,id);

    EXPECT_EQ(outMetadata,auxMetadata);
}

TEST_F( MetadataTest, JoinVector)
{
    MetaData auxMetadata, auxMetadata2;
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
    auxMetadata3.setValue(MDL_X,3.,id);
    auxMetadata3.setValue(MDL_Y,4.,id);
    auxMetadata3.setValue(MDL_Z,444.,id);

    id = auxMetadata2.addObject();
    auxMetadata2.setValue(MDL_X,1.,id);
    auxMetadata2.setValue(MDL_Y,2.,id);
    auxMetadata2.setValue(MDL_Z,3.,id);
    auxMetadata2.setValue(MDL_ANGLE_ROT,0.,id);
    id = auxMetadata2.addObject();
    auxMetadata2.setValue(MDL_X,3.,id);
    auxMetadata2.setValue(MDL_Y,4.,id);
    auxMetadata2.setValue(MDL_Z,5.,id);
    auxMetadata2.setValue(MDL_ANGLE_ROT,180.,id);

    std::vector<MDLabel> labels;
    labels.push_back(MDL_X);
    labels.push_back(MDL_Y);
    auxMetadata.join1(auxMetadata2,auxMetadata3,labels,LEFT);

    MetaData outMetadata;
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,1.,id);
    outMetadata.setValue(MDL_Y,2.,id);
    outMetadata.setValue(MDL_Z,3.,id);
    outMetadata.setValue(MDL_ANGLE_ROT,0.,id);
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,3.,id);
    outMetadata.setValue(MDL_Y,4.,id);
    outMetadata.setValue(MDL_Z,5.,id);
    outMetadata.setValue(MDL_ANGLE_ROT,180.,id);
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,3.,id);
    outMetadata.setValue(MDL_Y,4.,id);
    outMetadata.setValue(MDL_Z,5.,id);
    outMetadata.setValue(MDL_ANGLE_ROT,180.,id);

    EXPECT_EQ(outMetadata,auxMetadata);
}

TEST_F( MetadataTest, MDValueEQ)
{
    try
    {
        MetaData md;
        md.setValue(MDL_IMAGE, (String)"a", md.addObject());
        md.setValue(MDL_IMAGE, (String)"b", md.addObject());
        md.setValue(MDL_IMAGE, (String)"c", md.addObject());
        md.setValue(MDL_IMAGE, (String)"a", md.addObject());

        MetaData md2;
        md2.setValue(MDL_IMAGE, (String)"a", md2.addObject());
        md2.setValue(MDL_IMAGE, (String)"a", md2.addObject());

        MDValueEQ eq(MDL_IMAGE,(String)"a");
        //Test empty query
        MetaData md3;
        md3.importObjects(md, eq);
        EXPECT_EQ(md2, md3);
    }
    catch (XmippError &xe)
    {
        std::cerr << "DEBUG_JM: xe: " << xe << std::endl;
    }
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

    auxMetadata.joinNatural(mDsource,auxMetadata3);
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
#include <math.h>
TEST_F( MetadataTest, OperateExt)
{
    MetaData auxMetadata = mDunion;
    MetaData auxMetadata2 = mDunion;
    MDSql::activateMathExtensions();
    auxMetadata.operate((String)"X=sqrt(X)");
    double x;
    FOR_ALL_OBJECTS_IN_METADATA(auxMetadata2)
    {
        auxMetadata2.getValue(MDL_X,x,__iter.objId);
        auxMetadata2.setValue(MDL_X,sqrt(x),__iter.objId);
    }

    EXPECT_EQ(auxMetadata,auxMetadata2);
}

TEST_F( MetadataTest, RegularExp)
{
    XMIPP_TRY

    //create temporal file with three metadas
    //char sfnStar[64] = "";
    //char sfnSqlite[64] = "";
    //strncpy(sfnStar, "/tmp/testReadMultipleBlocks_XXXXXX.xmd", sizeof sfnStar);
    //if (mkstemps(sfnStar,4)==-1)
    //	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary STAR file");
//    strncpy(sfnSqlite, "/tmp/testReadMultipleBlocks_XXXXXX.sqlite", sizeof sfnSqlite);
//    if (mkstemps(sfnSqlite,7)==-1)
//    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary SQLITE file");

    FileName fn;
    fn.initUniqueName("/tmp/testReadMultipleBlocks_XXXXXX");
    FileName sfnStar;
    sfnStar = fn + ".xmd";
    MetaData auxMetadata;
    auxMetadata.setValue(MDL_IMAGE,(String)"image_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_2.xmp",auxMetadata.addObject());
    auxMetadata.write(sfnStar,MD_OVERWRITE);
//    auxMetadata.write(sfnSqlite,MD_OVERWRITE);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"block_000001@"+sfnStar,MD_APPEND);
//    auxMetadata.write((String)"block_000001@"+sfnSqlite,MD_APPEND);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"block_000002@"+sfnStar,MD_APPEND);
//    auxMetadata.write((String)"block_000002@"+sfnSqlite,MD_APPEND);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_no_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_no_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"noblock@"+sfnStar,MD_APPEND);
//    auxMetadata.write((String)"noblock@"+sfnSqlite,MD_APPEND);
    auxMetadata.clear();
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_3_1.xmp",auxMetadata.addObject());
    auxMetadata.setValue(MDL_IMAGE,(String)"image_data_3_2.xmp",auxMetadata.addObject());
    auxMetadata.write((String)"block_000003@"+sfnStar,MD_APPEND);
//    auxMetadata.write((String)"block_000003@"+sfnSqlite,MD_APPEND);
    auxMetadata.clear();
    MetaData auxMetadata2;
    auxMetadata2.setValue(MDL_IMAGE,(String)"image_data_1_1.xmp",auxMetadata2.addObject());
    auxMetadata2.setValue(MDL_IMAGE,(String)"image_data_1_2.xmp",auxMetadata2.addObject());
    auxMetadata2.setValue(MDL_IMAGE,(String)"image_data_2_1.xmp",auxMetadata2.addObject());
    auxMetadata2.setValue(MDL_IMAGE,(String)"image_data_2_2.xmp",auxMetadata2.addObject());
    MDSql::activateRegExtensions();
    auxMetadata2.write("/tmp/kk.sqlite");
    //query file
    MetaData md;
    FileName blockFileName;
    md.read((String)"block_000001@"+sfnStar);
    //compare with reference metada
//    EXPECT_EQ(md,auxMetadata2);
//    md.read((String)"block_000001@"+sfnSqlite);
    //compare with reference metada
//    EXPECT_EQ(md,auxMetadata2);
    unlink(fn.c_str());
    unlink(sfnStar.c_str());
    //unlink(sfnSqlite);

    XMIPP_CATCH
}

TEST_F( MetadataTest, Query)
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
    auxMetadata3.setValue(MDL_X,3.,id);
    auxMetadata3.setValue(MDL_Y,4.,id);
    auxMetadata3.setValue(MDL_Z,444.,id);

    auxMetadata.importObjects(auxMetadata3, MDValueEQ(MDL_X, 3.));

    MetaData outMetadata;
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,3.,id);
    outMetadata.setValue(MDL_Y,4.,id);
    outMetadata.setValue(MDL_Z,333.,id);
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,3.,id);
    outMetadata.setValue(MDL_Y,4.,id);
    outMetadata.setValue(MDL_Z,444.,id);

    EXPECT_EQ(outMetadata,auxMetadata);
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
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

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
#define sizesfn 64
    char sfn[sizesfn] = "";
    strncpy(sfn, "/tmp/testReadMultipleBlocks_XXXXXX", sizesfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

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
    EXPECT_EQ(auxMetadata.size(),(size_t)0);

    auxMetadata.read((String)"block_000004@"+sfn);
    EXPECT_EQ(auxMetadata.size(),(size_t)0);

    unlink(sfn);
}

TEST_F( MetadataTest, ReadEmptyBlocksII)
{
    char sfn[64] = "";
    strncpy(sfn, "/tmp/testReadMultipleBlocks_XXXXXX", sizeof sfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

    MetaData auxMetadata;

    auxMetadata.addLabel(MDL_X);
    auxMetadata.addLabel(MDL_Y);
    auxMetadata.addLabel(MDL_Z);
    auxMetadata.write((String)"block_000002@"+sfn,MD_APPEND);

    auxMetadata.read((String)"block_000002@"+sfn);
    EXPECT_EQ(auxMetadata.size(),(size_t)0);
    unlink(sfn);
}

TEST_F( MetadataTest, ReadWrite)
{
    //temp file name
    char sfn[32] = "";
    strncpy(sfn, "/tmp/testWrite_XXXXXX", sizeof sfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
    mDsource.write(sfn);
    MetaData auxMetadata;
    auxMetadata.read(sfn);

    EXPECT_EQ(mDsource,auxMetadata);
    unlink(sfn);
}

TEST_F( MetadataTest, WriteIntermediateBlock)
{
    //read metadata block between another two
    FileName filename("metadata/WriteIntermediateBlock.xmd");
    FileName blockFileName;
    blockFileName.compose("two", filename);
    MetaData auxMetadata(blockFileName);
    MDRow row;
    row.setValue(MDL_X, 11.);
    row.setValue(MDL_Y, 22.);
    auxMetadata.addRow(row);
    row.setValue(MDL_X, 33.);
    row.setValue(MDL_Y, 44.);
    auxMetadata.addRow(row);
    auxMetadata.setValue(MDL_X,111.,auxMetadata.firstObject());

    //temporal file for modified metadata
    char sfn2[32] = "";
    strncpy(sfn2, "/tmp/testWrite_XXXXXX", sizeof sfn2);
    if (mkstemp(sfn2)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");

    //copy input metadata file
    std::ifstream src; // the source file
    std::ofstream dest; // the destination file
    src.open (filename.c_str(), std::ios::binary); // open in binary to prevent jargon at the end of the buffer
    dest.open (sfn2, std::ios::binary); // same again, binary
    if (!src.is_open())
        std::cerr << "Can not open file: " << filename.c_str() <<std::endl; // could not be copied
    if (!dest.is_open())
        std::cerr << "Can not open file: " <<sfn2 <<std::endl; // could not be copied
    dest << src.rdbuf (); // copy the content
    dest.close (); // close destination file
    src.close (); // close source file

    blockFileName.compose("two", sfn2);
    auxMetadata.write(blockFileName,MD_APPEND);
    //file with correct values
    FileName fn2("metadata/ReadWriteAppendBlock2.xmd");
    EXPECT_TRUE(compareTwoFiles("metadata/WriteIntermediateBlock2.xmd",sfn2,0));
    unlink(sfn2);
}

TEST_F( MetadataTest, ExistsBlock)
{
    //temp file name
    char sfn[32] = "";
    strncpy(sfn, "/tmp/testWrite_XXXXXX", sizeof sfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
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
    XMIPP_TRY
    //temp file name
    char sfn[32] = "";
    strncpy(sfn, "/tmp/testWrite_XXXXXX", sizeof sfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
    mDsource.write((String)"one@"+sfn);
    mDsource.write((String)"two@"+sfn,MD_APPEND);
    mDsource.write((String)"three@"+sfn,MD_APPEND);
    MetaData auxMetadata;
    FileName sfn2 = "metadata/ReadWriteAppendBlock.xmd";
    EXPECT_TRUE(compareTwoFiles(sfn,sfn2,0));
    unlink(sfn);
    XMIPP_CATCH
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

TEST_F( MetadataTest, Distinct)
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
    auxMetadata1.distinct(auxMetadata3,MDL_X);
    auxMetadata3.clear();
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_X,1.,id);
    id = auxMetadata3.addObject();
    auxMetadata3.setValue(MDL_X,3.,id);

    EXPECT_EQ(auxMetadata1,auxMetadata3);
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
    EXPECT_EQ((size_t)2, mDsource.size());
}

TEST_F( MetadataTest, Sort)
{
    MetaData auxMetadata,auxMetadata2,auxMetadata3,outMetadata;
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_X,3.,id);
    auxMetadata.setValue(MDL_Y,4.,id);
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_X,1.,id);
    auxMetadata.setValue(MDL_Y,2.,id);
    auxMetadata2.sort(auxMetadata,MDL_X);
    EXPECT_EQ(auxMetadata2,mDsource);

    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_X,5.,id);
    auxMetadata.setValue(MDL_Y,6.,id);

    auxMetadata2.clear();
    auxMetadata2.sort(auxMetadata,MDL_X,true,1,0);
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,1.,id);
    outMetadata.setValue(MDL_Y,2.,id);
    EXPECT_EQ(auxMetadata2,outMetadata);

    auxMetadata2.clear();
    auxMetadata2.sort(auxMetadata,MDL_X,true,2,1);
    outMetadata.clear();
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,3.,id);
    outMetadata.setValue(MDL_Y,4.,id);
    id = outMetadata.addObject();
    outMetadata.setValue(MDL_X,5.,id);
    outMetadata.setValue(MDL_Y,6.,id);
    EXPECT_EQ(auxMetadata2,outMetadata);
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
    EXPECT_EQ(MDL::labelType(MDL_ORDER),LABEL_SIZET);
    MetaData auxMetadata;
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_ORDER,(size_t)1, id);
    auxMetadata.getValue(MDL_ORDER,t, id);
    EXPECT_EQ((size_t)1,t);
    //We expect that MetaData will throw an exception
    //if you use getValue with a variable of type that
    // doesn't match the label type
    std::cerr << "TEST COMMENT: you should get the ERROR: Mismatch Label (order_) and value type(INT)" <<std::endl;
    EXPECT_THROW(auxMetadata.getValue(MDL_ORDER, i, id), XmippError);
}
TEST_F( MetadataTest, Comment)
{
    XMIPP_TRY
    char sfn[64] = "";
    MetaData md1(mDsource);
    strncpy(sfn, "/tmp/testComment_XXXXXX", sizeof sfn);
    if (mkstemp(sfn)==-1)
    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
    String s1((String)"This is a very long comment that has more than 80 characters"+
              " Therefore should be split in several lines"+
              " Let us see what happened");
    md1.setComment(s1);
    md1.write(sfn, MD_OVERWRITE);
    MetaData md2;
    md2.read(sfn);
    String s2;
    s2=md2.getComment();
    EXPECT_EQ(s1, s2);
    unlink(sfn);

    XMIPP_CATCH
}
//read file with vector
TEST_F( MetadataTest, getValue)
{
    XMIPP_TRY
    std::vector<double> v1(3);
    std::vector<double> v2(3);
    MetaData auxMD1;
    id = auxMD1.addObject();
    v1[0]=1.;
    v1[1]=2.;
    v1[2]=3.;
    auxMD1.setValue(MDL_CLASSIFICATION_DATA,v1,id);
    id = auxMD1.firstObject();
    auxMD1.getValue(MDL_CLASSIFICATION_DATA,v2, id);

    EXPECT_EQ(v1[0],v2[0]);
    EXPECT_EQ(v1[1],v2[1]);
    EXPECT_EQ(v1[2],v2[2]);
    XMIPP_CATCH
}

TEST_F( MetadataTest, getValueDefault)
{
    XMIPP_TRY
    MetaData auxMD1;
    MetaData auxMD2;
    double rot=1., tilt=2., psi=3.;
    double rot2=0., tilt2=0., psi2=0.;
    id = auxMD1.addObject();
    auxMD1.setValue(MDL_ANGLE_ROT,rot,id);
    auxMD1.setValue(MDL_ANGLE_TILT,tilt,id);
    //psi assigned by defaults
    id = auxMD1.firstObject();
    auxMD1.getValueOrDefault(MDL_ANGLE_ROT,rot2, id, 0.);
    auxMD1.getValueOrDefault(MDL_ANGLE_TILT,tilt2, id, 0.);
    auxMD1.getValueOrDefault(MDL_ANGLE_PSI,psi2, id, 3.);

    EXPECT_EQ(rot,rot2);
    EXPECT_EQ(tilt,tilt2);
    EXPECT_EQ(psi,psi2);

    MDRow  rowIn;
    psi2=0;
    auxMD1.getRow(rowIn, id);
    rowIn.getValueOrDefault(MDL_ANGLE_PSI,psi2,3.);
    EXPECT_EQ(psi,psi2);

    auxMD1.getRow2(rowIn, id);
    rowIn.getValueOrDefault(MDL_ANGLE_PSI,psi2,3.);
    EXPECT_EQ(psi,psi2);

    XMIPP_CATCH
}
TEST_F( MetadataTest, getValueAbort)
{
    XMIPP_TRY
    MetaData auxMD1;
    double rot=1.;
    id = auxMD1.addObject();
    auxMD1.setValue(MDL_ANGLE_ROT,rot,id);
    //psi assigned by defaults
    id = auxMD1.firstObject();
    std::cerr << "TEST COMMENT: You should get the error  Cannot find label: order_" <<std::endl;
    EXPECT_THROW(auxMD1.getValueOrAbort(MDL_ORDER, rot, id), XmippError);
    MDRow  rowIn;
    auxMD1.getRow(rowIn, id);
    std::cerr << "TEST COMMENT: You should get the error  Cannot find label: anglePsi" <<std::endl;
    EXPECT_THROW(rowGetValueOrAbort(rowIn,MDL_ANGLE_PSI,rot), XmippError);
    XMIPP_CATCH
}

TEST_F( MetadataTest, CopyColumn)
{
    XMIPP_TRY
    MetaData md1(mDsource), md2(mDsource);
    double value;

    FOR_ALL_OBJECTS_IN_METADATA(md1)
    {
        md1.getValue(MDL_Y, value, __iter.objId);
        md1.setValue(MDL_Z, value, __iter.objId);
    }

    md2.copyColumn(MDL_Z, MDL_Y);

    EXPECT_EQ(md1, md2);
    XMIPP_CATCH
}

TEST_F( MetadataTest, RenameColumn)
{
    XMIPP_TRY
    MetaData md1(mDsource);
    MetaData md2;
    md1.renameColumn(MDL_Y,MDL_Z);
    id = md2.addObject();
    md2.setValue(MDL_X,1.,id);
    md2.setValue(MDL_Z,2.,id);
    id = md2.addObject();
    md2.setValue(MDL_X,3.,id);
    md2.setValue(MDL_Z,4.,id);


    EXPECT_EQ(md1, md2);
    XMIPP_CATCH
}

TEST_F( MetadataTest, BsoftRemoveLoopBlock)
{
//    XMIPP_TRY
//    FileName fnSTAR =(String)"metadata/symop.star";
//    char sfn[64] = "";
//    strncpy(sfn, "/tmp/BsoftRemoveLoopBlock1_XXXXXX.star", sizeof sfn);
//    if (mkstemps(sfn,5)==-1)
//    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
//    bsoftRemoveLoopBlock(fnSTAR,sfn);
//    size_t id;
//    MetaData readMd;
//
//    //md1a
//    MetaData md1a;
//    md1a.setColumnFormat(false);
//    id = md1a.addObject();
//    md1a.setValue(BSOFT_SYMMETRY_INT_TABLES_NUMBER,1,id);
//    md1a.setValue(BSOFT_SYMMETRY_SPACE_GROUP_NAME_H_M,(String)"P1",id);
//    md1a.setValue(BSOFT_SYMMETRY_CELL_SETTING,(String)"TRICLINIC",id);
//    readMd.read((std::string)"A1@"+sfn);
//    EXPECT_EQ(md1a, readMd);
//
//    //md1b
//    MetaData md1b;
//    md1a.setColumnFormat(true);
//    id = md1b.addObject();
//    md1b.setValue(BSOFT_SYMMETRY_EQUIV_ID,1,id);
//    md1b.setValue(BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ,(String)"X,Y,Z",id);
//    readMd.read((std::string)"loop_1@"+sfn);
//    EXPECT_EQ(md1b, readMd);
//
//    //md2a
//    MetaData md2a;
//    md1a.setColumnFormat(false);
//    id = md2a.addObject();
//    md2a.setValue(BSOFT_SYMMETRY_INT_TABLES_NUMBER,2,id);
//    md2a.setValue(BSOFT_SYMMETRY_SPACE_GROUP_NAME_H_M,(String)"P-1",id);
//    md2a.setValue(BSOFT_SYMMETRY_CELL_SETTING,(String)"TRICLINIC",id);
//    readMd.read((std::string)"A2@"+sfn);
//    EXPECT_EQ(md2a, readMd);
//    //md2b
//    MetaData md2b;
//    md1a.setColumnFormat(true);
//    id = md2b.addObject();
//    md2b.setValue(BSOFT_SYMMETRY_EQUIV_ID,1,id);
//    md2b.setValue(BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ,(String)"X,Y,Z",id);
//    id = md2b.addObject();
//    md2b.setValue(BSOFT_SYMMETRY_EQUIV_ID,2,id);
//    md2b.setValue(BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ,(String)"-X,-Y,-Z",id);
//    readMd.read((std::string)"loop_2@"+sfn);
//    EXPECT_EQ(md2b, readMd);
//
//    //md5090row
//    MetaData md5090row;
//    md5090row.setColumnFormat(true);
//    id = md5090row.addObject();
//    md5090row.setValue(BSOFT_SYMMETRY_EQUIV_ID,1,id);
//    md5090row.setValue(BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ,(String)"X,Y,Z",id);
//    id = md5090row.addObject();
//    md5090row.setValue(BSOFT_SYMMETRY_EQUIV_ID,2,id);
//    md5090row.setValue(BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ,(String)"-X,-Y,Z",id);
//    id = md5090row.addObject();
//    md5090row.setValue(BSOFT_SYMMETRY_EQUIV_ID,3,id);
//    md5090row.setValue(BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ,(String)"-Y,X,Z",id);
//    readMd.read((std::string)"A5090row@"+sfn);
//    EXPECT_EQ(md5090row, readMd);
////    //md5090col
//    MetaData md5090col;
//    md5090col.setColumnFormat(false);
//    id = md5090col.addObject();
//    md5090col.setValue(BSOFT_SYMMETRY_INT_TABLES_NUMBER,5090,id);
//    md5090col.setValue(BSOFT_SYMMETRY_SPACE_GROUP_NAME_H_M,(String)"P4212",id);
//    md5090col.setValue(BSOFT_SYMMETRY_CELL_SETTING,(String)"TETRAGONAL_4axis",id);
//    readMd.read((std::string)"A5090col@"+sfn);
//    EXPECT_EQ(md5090col, readMd);
//
//    unlink(sfn);
//    XMIPP_CATCH

}

TEST_F( MetadataTest, bsoftRestoreLoopBlock)
{
//    XMIPP_TRY
//    FileName fnSTAR =(String)"metadata/symop.star";
//    char sfn[64] = "";
//    strncpy(sfn, "/tmp/BsoftRemoveLoopBlock1_XXXXXX.star", sizeof sfn);
//    if (mkstemps(sfn,5)==-1)
//    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
//    char sfn2[64] = "";
//    strncpy(sfn2, "/tmp/BsoftRemoveLoopBlock2_XXXXXX.star", sizeof sfn);
//    if (mkstemps(sfn2,5)==-1)
//    	REPORT_ERROR(ERR_IO_NOTOPEN,"Cannot create temporary file");
//    size_t id;
//    MetaData readMd;
//
//    //md1a
//    MetaData md1a;
//    md1a.setColumnFormat(false);
//    id = md1a.addObject();
//    md1a.setValue(BSOFT_SYMMETRY_INT_TABLES_NUMBER,1,id);
//    md1a.setValue(BSOFT_SYMMETRY_SPACE_GROUP_NAME_H_M,(String)"P1",id);
//    md1a.setValue(BSOFT_SYMMETRY_CELL_SETTING,(String)"TRICLINIC",id);
//    md1a.write((std::string)"A1@"+sfn, MD_OVERWRITE);
//    //md1b
//    MetaData md1b;
//    md1a.setColumnFormat(true);
//    id = md1b.addObject();
//    md1b.setValue(BSOFT_SYMMETRY_EQUIV_ID,1,id);
//    md1b.setValue(BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ,(String)"X,Y,Z",id);
//    md1b.write((std::string)"loop_1@"+sfn, MD_APPEND);
//
//    bsoftRestoreLoopBlock(sfn,sfn2);
//    bsoftRemoveLoopBlock(sfn2,sfn);
//    readMd.read((std::string)"A1@"+sfn);
//    EXPECT_EQ(md1a, readMd);
//    readMd.read((std::string)"loop_1@"+sfn);
//    EXPECT_EQ(md1b, readMd);
//    unlink(sfn);
//    unlink(sfn2);
//
//    XMIPP_CATCH
}

//Copy images on metadata using ImageConvert logic
TEST_F( MetadataTest, copyImages)
{
    XMIPP_TRY
    FileName fn = "metadata/smallStack.stk";
    FileName out;
    FileName oroot;
    oroot.initUniqueName("/tmp/smallImg_XXXXXX");
    oroot.deleteFile();
    out = oroot.addExtension("xmd");
    oroot = oroot + ":mrc";

    FileName fn1, fn2;
    MetaData md(fn);
    ProgConvImg conv;
    conv.verbose = 0;
    conv.setup(&md, out, oroot);
    conv.tryRun();
    MetaData *mdOut = conv.getOutputMd();

    FOR_ALL_OBJECTS_IN_METADATA2(md, *mdOut)
    {
        md.getValue(MDL_IMAGE, fn1, __iter.objId);
        mdOut->getValue(MDL_IMAGE, fn2, __iter2.objId);
        EXPECT_TRUE(compareImage(fn1, fn2));
        fn2.deleteFile();
    }

    out.deleteFile();
    out.initUniqueName("/tmp/smallStack_XXXXXX");
    out = out + ":mrcs";
    conv.setup(&md, out);
    conv.tryRun();

    FOR_ALL_OBJECTS_IN_METADATA2(md, *mdOut)
    {
        md.getValue(MDL_IMAGE, fn1, __iter.objId);
        mdOut->getValue(MDL_IMAGE, fn2, __iter2.objId);
        EXPECT_TRUE(compareImage(fn1, fn2));
    }
    out.deleteFile();

    out.initUniqueName("/tmp/smallStackVol_XXXXXX");
    out = out + ":mrc";
    conv.setup(&md, out);
    conv.setType("vol");
    conv.tryRun();
    Image<float> imgStk, imgVol;
    imgStk.read(fn);
    imgVol.read(out);

    int n = imgStk.getDimensions().ndim;
    for (int i = FIRST_IMAGE; i <= n; ++i)
    {
        imgStk.movePointerTo(1, i);
        imgVol.movePointerTo(i);
        EXPECT_TRUE(imgStk == imgVol);
    }
    out.deleteFile();
    XMIPP_CATCH
}

TEST_F( MetadataTest, updateRow)
{
    ASSERT_EQ(mDsource,mDsource);
    ASSERT_FALSE(mDsource==mDanotherSource);
    //attribute order should not be important
    MetaData auxMetadata ;
    id1 = auxMetadata.addObject();
    auxMetadata.setValue(MDL_Y,0.,id);
    auxMetadata.setValue(MDL_X,0.,id);
    id2 = auxMetadata.addObject();
    auxMetadata.setValue(MDL_Y,0.,id);
    auxMetadata.setValue(MDL_X,0.,id);
    ASSERT_FALSE(auxMetadata==mDsource);

    MDRow row;
    row.setValue(MDL_X, 1.);
    row.setValue(MDL_Y, 2.);
    auxMetadata.setRow( row, id1);
    row.setValue(MDL_X, 3.);
    row.setValue(MDL_Y, 4.);
    auxMetadata.setRow( row, id2);
    ASSERT_EQ(auxMetadata,mDsource);

    row.setValue(MDL_X, 1.);
    row.setValue(MDL_Y, 2.);
    auxMetadata.setRow2( row, id1);
    row.setValue(MDL_X, 3.);
    row.setValue(MDL_Y, 4.);
    auxMetadata.setRow2( row, id2);
    ASSERT_EQ(auxMetadata,mDsource);

    //Test form double with a given precission.
/*    auxMetadata.clear();
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

    auxMetadata.setValue(MDL_Y,2.,id);
    auxMetadata.setValue(MDL_X,1.,id);
    id = auxMetadata.addObject();
    auxMetadata.setValue(MDL_Y,4.,id);
    auxMetadata.setValue(MDL_X,3.,id);

    MDRow row;
    row.setValue(MDL_X, 1.);
    row.setValue(MDL_Y, 2.);
    md.addRow(row);
    */
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
