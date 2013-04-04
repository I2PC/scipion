#include <dimred/dimred_tools.h>
#include <dimred/kernelPCA.h>
#include <dimred/ltsa.h>
#include <dimred/diffusionMaps.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide

class DimRedTest : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        chdir(((String)(getXmippPath() + (String)"/resources/test")).c_str());
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
	EXPECT_LT(fabs(dimCorrDim-expectedDim),1e-6);
}

TEST_F( DimRedTest, kernelPCA)
{
	GenerateData generator;
	generator.generateNewDataset("swiss",1000,0);
	generator.X.write("dimred/swiss.txt");
	// MATLAB: load swiss.txt;

	KernelPCA kernelPCA;
	kernelPCA.setInputData(generator.X);
	kernelPCA.setOutputDimensionality(2);
	kernelPCA.setSpecificParameters(0.1);
	kernelPCA.reduceDimensionality();
	const Matrix2D<double> &Y=kernelPCA.getReducedData();

	// Compare C output with Matlab output
	// ***
	Y.write("dimred/kernelPCA.txt");
}

TEST_F( DimRedTest, ltsa)
{
	GenerateData generator;
	generator.generateNewDataset("helix",1000,0);
	// generator.X.write("dimred/helix.txt");
	// MATLAB: load swiss.txt;

	LTSA ltsa;
	ltsa.setInputData(generator.X);
	ltsa.setOutputDimensionality(2);
	ltsa.setSpecificParameters();
	ltsa.reduceDimensionality();
	const Matrix2D<double> &Y=ltsa.getReducedData();

	//Y.write("dimred/ltsa.txt");
	Matrix2D<double> expectedY;
	expectedY.resizeNoCopy(Y);
	expectedY.read("dimred/ltsa.txt");
	ASSERT_TRUE(expectedY.equal(Y,1e-5));
}

TEST_F( DimRedTest, diffusionMaps)
{
	GenerateData generator;
	generator.generateNewDataset("helix",1000,0);
	// generator.X.write("dimred/helix.txt");
	// MATLAB: load swiss.txt;

	DiffusionMaps dimred;
	dimred.setInputData(generator.X);
	dimred.setOutputDimensionality(2);
	dimred.setSpecificParameters();
	dimred.reduceDimensionality();
	const Matrix2D<double> &Y=dimred.getReducedData();

	//Y.write("dimred/diffusionMaps.txt");
	Matrix2D<double> expectedY;
	expectedY.resizeNoCopy(Y);
	expectedY.read("dimred/diffusionMaps.txt");
	ASSERT_TRUE(expectedY.equal(Y,1e-5));
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
