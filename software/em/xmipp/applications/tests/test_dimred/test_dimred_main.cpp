#include <dimred/lpp.h>
#include <dimred/nca.h>
#include <dimred/npe.h>
#include <dimred/spe.h>
#include <dimred/ltsa.h>
#include <dimred/gplvm.h>
#include <dimred/lltsa.h>
#include <dimred/kernelPCA.h>
#include <dimred/hessianLLE.h>
#include <dimred/diffusionMaps.h>
#include <dimred/probabilisticPCA.h>
#include <dimred/laplacianEigenmaps.h>
#include <iostream>
#include <stdlib.h>     /* getenv */
#include <gtest/gtest.h>
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide

class DimRedTest : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        if (chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str())==-1)
        	REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot change directory");
    }
};

TEST_F( DimRedTest, generate_data)
{
	GenerateData generator;

	// Swiss
	generator.generateNewDataset("swiss",1000,0);
	// generator.X.write("dimred/swiss.txt");
	// MATLAB: load swiss.txt; plot3(swiss(:,1),swiss(:,2),swiss(:,3),'.');
	Matrix2D<double> expectedX;
	expectedX.resizeNoCopy(generator.X);
	expectedX.read("dimred/swiss.txt");
	ASSERT_TRUE(expectedX.equal(generator.X,1e-5));

	// Helix
	generator.generateNewDataset("helix",1000,0);
	// generator.X.write("dimred/helix.txt");
	// MATLAB: load helix.txt; plot3(helix(:,1),helix(:,2),helix(:,3),'.');
	expectedX.resizeNoCopy(generator.X);
	expectedX.read("dimred/helix.txt");
	ASSERT_TRUE(expectedX.equal(generator.X,1e-5));

	// Twin peaks
	generator.generateNewDataset("twinpeaks",1000,0);
	// generator.X.write("dimred/twinpeaks.txt");
	// MATLAB: load twinpeaks.txt; plot3(twinpeaks(:,1),twinpeaks(:,2),twinpeaks(:,3),'.');
	expectedX.resizeNoCopy(generator.X);
	expectedX.read("dimred/twinpeaks.txt");
	ASSERT_TRUE(expectedX.equal(generator.X,1e-5));

	// Clusters
	generator.generateNewDataset("3d_clusters",1000,0);
	// generator.X.write("dimred/clusters.txt");
	// MATLAB: load clusters.txt; plot3(clusters(:,1),clusters(:,2),clusters(:,3),'.');
	expectedX.resizeNoCopy(generator.X);
	expectedX.read("dimred/clusters.txt");
	ASSERT_TRUE(expectedX.equal(generator.X,1e-5));

	// Intersect
	generator.generateNewDataset("intersect",1000,0);
	// generator.X.write("dimred/intersect.txt");
	// MATLAB: load intersect.txt; plot3(intersect(:,1),intersect(:,2),intersect(:,3),'.');
	expectedX.resizeNoCopy(generator.X);
	expectedX.read("dimred/intersect.txt");
	ASSERT_TRUE(expectedX.equal(generator.X,1e-5));
}

TEST_F( DimRedTest, intrinsic_dimensionality)
{
	GenerateData generator;
	generator.generateNewDataset("swiss",1000,0);
	// generator.X.write("dimred/swiss.txt");
	// MATLAB: load swiss.txt;

	double dimMLE=intrinsicDimensionality(generator.X,"MLE");
	// generator.X.write("dimred/swissNormalized.txt");
	// MATLAB: load swissNormalized.txt; mean(swissNormalized); std(swissNormalized); d=intrinsic_dimension(swissNormalized)
	double expectedDim=1.927789055150985;
	EXPECT_LT(fabs(dimMLE-expectedDim),1e-6);

	double dimCorrDim=intrinsicDimensionality(generator.X,"CorrDim",false);
	expectedDim=1.9244901554639233;
	EXPECT_LT(fabs(dimCorrDim-expectedDim),5e-2);
}

#define INCOMPLETE_TEST(method,DimredClass,dataset,Npoints,file) \
	TEST_F( DimRedTest, method) \
{ \
	GenerateData generator; \
	generator.generateNewDataset(dataset,Npoints,0); \
	DimredClass dimred; \
	dimred.setInputData(generator.X); \
	dimred.setOutputDimensionality(2); \
	dimred.setSpecificParameters(); \
	dimred.reduceDimensionality(); \
	const Matrix2D<double> &Y=dimred.getReducedData(); \
	Y.write(file); \
}

#define COMPLETE_TEST(method,DimredClass,dataset,Npoints,file,doTestAlways) \
	TEST_F( DimRedTest, method) \
{ \
	const char * doTest = getenv ("XMIPP_FULLTEST");\
	if (doTest==NULL && !doTestAlways)\
	{\
       std::cout << "Skipping test: DimRedTest using file " << file << ".\n Define environmental variable XMIPP_FULLTEST to activate test"<< std::endl;\
       ASSERT_TRUE(true);\
       return;\
    }\
	GenerateData generator; \
	generator.generateNewDataset(dataset,Npoints,0); \
	DimredClass dimred; \
	dimred.setInputData(generator.X); \
	dimred.setOutputDimensionality(2); \
	dimred.setSpecificParameters(); \
	dimred.reduceDimensionality(); \
	const Matrix2D<double> &Y=dimred.getReducedData();\
	/* Y.write(file); */ \
	Matrix2D<double> expectedY; \
	expectedY.resizeNoCopy(Y); \
	expectedY.read(file); \
	ASSERT_TRUE(expectedY.equalAbs(Y,1e-4));\
}

COMPLETE_TEST(ltsa,               LTSA,             "helix",1000,"dimred/ltsa.txt",true)
COMPLETE_TEST(diffusionMaps,      DiffusionMaps,    "helix",1000,"dimred/diffusionMaps.txt",true)
COMPLETE_TEST(lltsa,              LLTSA,            "helix",1000,"dimred/lltsa.txt",false)
COMPLETE_TEST(lpp,                LPP,              "helix",1000,"dimred/lpp.txt",true)
COMPLETE_TEST(kernelPCA,          KernelPCA,        "helix",1000,"dimred/kernelPCA.txt",false)
COMPLETE_TEST(probabilisticPCA,   ProbabilisticPCA, "helix",1000,"dimred/probabilisticPCA.txt",true)
COMPLETE_TEST(laplacianEigenmap,LaplacianEigenmap,  "helix",1000,"dimred/laplacianEigenmap.txt",true)
COMPLETE_TEST(hessianlle,         HessianLLE,       "helix",1000,"dimred/hessianlle.txt",true)
COMPLETE_TEST(spe,                SPE,              "helix",1000,"dimred/spe.txt",true)
COMPLETE_TEST(npe,                NPE,              "helix",1000,"dimred/npe.txt",false)

TEST_F( DimRedTest, nca)
{
	GenerateData generator;
	generator.generateNewDataset("helix",1000,0);
	NeighbourhoodCA dimred;
	dimred.setInputData(generator.X);
	dimred.setOutputDimensionality(2);
	dimred.setSpecificParameters();
	Matrix1D<unsigned char> labels;
	labels.resizeNoCopy(1000);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(labels)
		VEC_ELEM(labels,i)=(i<500);
	dimred.setLabels(labels);
	dimred.reduceDimensionality();
	const Matrix2D<double> &Y=dimred.getReducedData();
	// Y.write("dimred/nca.txt");
	Matrix2D<double> expectedY;
	expectedY.resizeNoCopy(Y);
	expectedY.read("dimred/nca.txt");
	ASSERT_TRUE(expectedY.equal(Y,1e-4));
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
