#include <data/image.h>
#include <data/filters.h>
#include <data/fftw.h>
#include <data/polar.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class PolarTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
        mulDouble.resize(3,3);
        DIRECT_A2D_ELEM(mulDouble,0,0) = 1;
        DIRECT_A2D_ELEM(mulDouble,0,1) = 2;
        DIRECT_A2D_ELEM(mulDouble,0,2) = 3;

        DIRECT_A2D_ELEM(mulDouble,1,0) = 3;
        DIRECT_A2D_ELEM(mulDouble,1,1) = 2;
        DIRECT_A2D_ELEM(mulDouble,1,2) = 1;

        DIRECT_A2D_ELEM(mulDouble,2,0) = 4;
        DIRECT_A2D_ELEM(mulDouble,2,1) = 4;
        DIRECT_A2D_ELEM(mulDouble,2,2) = 5;
    }
    MultidimArray<  double  > mulDouble;

    // virtual void TearDown() {}//Destructor

};
TEST_F( PolarTest, computeAverageAndStddev)
{
	Polar<double>                 P;
	double mean, stddev;
	P.getPolarFromCartesianBSpline(mulDouble,0,1);
    P.computeAverageAndStddev(mean,stddev,true);
    EXPECT_DOUBLE_EQ(mean,1.886528450043468);
    EXPECT_DOUBLE_EQ(stddev,0.49643800057938808);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
