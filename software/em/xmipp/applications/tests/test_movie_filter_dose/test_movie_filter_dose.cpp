#include <reconstruction/movie_filter_dose.h>
#include <iostream>
#include <gtest/gtest.h>
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class MovieFilterDoseTest : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
    	int a=0;
    }

    // virtual void TearDown() {}//Destructor

};

TEST_F( MovieFilterDoseTest, doseFilter)
{
	ProgMovieFilterDose pmfdProg(200);
	double doseFilterValue;
	double dose_finish =  4.;
	double current_critical_dose =  412084.3;
	double doseFilter =  0.9999952;

	doseFilterValue = pmfdProg.doseFilter(dose_finish,current_critical_dose);
	EXPECT_FLOAT_EQ(doseFilter,doseFilterValue);

	dose_finish =  4.000000;
	current_critical_dose =  12.82717;
	doseFilter =  0.8556285;
	doseFilterValue = pmfdProg.doseFilter(dose_finish,current_critical_dose);
	EXPECT_FLOAT_EQ(doseFilter,doseFilterValue);
}

TEST_F( MovieFilterDoseTest, movieFilterDose)
{
	double acceleration_voltage = 300;
	double voltage_scaling_factor = 1;
	ProgMovieFilterDose pmfdProg(acceleration_voltage);
	EXPECT_FLOAT_EQ(acceleration_voltage,pmfdProg.acceleration_voltage);
	EXPECT_FLOAT_EQ(voltage_scaling_factor,pmfdProg.voltage_scaling_factor);

	acceleration_voltage = 200;
	voltage_scaling_factor = 0.8;
	ProgMovieFilterDose pmfdProg2(acceleration_voltage);
	EXPECT_FLOAT_EQ(acceleration_voltage,pmfdProg2.acceleration_voltage);
	EXPECT_FLOAT_EQ(voltage_scaling_factor,pmfdProg2.voltage_scaling_factor);

}

TEST_F( MovieFilterDoseTest, criticalDose)
{
	double acceleration_voltage = 300;
	double voltage_scaling_factor = 1;
	ProgMovieFilterDose pmfdProg(acceleration_voltage);
    double criticalDose2;
	double criticalDose=  412084.3;
	double spatial_frequency=  1.8219448E-04;
	criticalDose2 =  pmfdProg.criticalDose(spatial_frequency);
	EXPECT_EQ(int(criticalDose2),int(criticalDose));

	criticalDose=  4.163977;
	spatial_frequency=  0.3587903;
	criticalDose2 = pmfdProg.criticalDose(spatial_frequency);
	EXPECT_FLOAT_EQ(criticalDose2,criticalDose);

	criticalDose=  200000;
	spatial_frequency=  0.3587903;
	criticalDose2 = pmfdProg.criticalDose(spatial_frequency);
	EXPECT_NE(criticalDose2,criticalDose);
}

TEST_F( MovieFilterDoseTest, optimalDoseGivenCriticalDose)
{
	double acceleration_voltage = 300;
	double current_critical_dose = 38.49693;
	double current_optimal_dose = 96.73663;
	double current_optimal_dose2 = 0.;
	ProgMovieFilterDose pmfdProg(acceleration_voltage);
	current_optimal_dose2 =  pmfdProg.optimalDoseGivenCriticalDose(current_critical_dose);
	EXPECT_FLOAT_EQ(current_optimal_dose2,current_optimal_dose);

	current_optimal_dose=  200000;
	EXPECT_NE(current_optimal_dose2,current_optimal_dose);

}

/*
TEST_F( FftwTest, directFourierTransformComplex)
{

    MultidimArray< std::complex< double > > FFT1, complxDouble;
    FourierTransformer transformer1;
    typeCast(mulDouble, complxDouble);
    transformer1.FourierTransform(complxDouble, FFT1, false);
    transformer1.inverseFourierTransform();

    transformer1.inverseFourierTransform();
    MultidimArray<std::complex<double> > auxFFT;
    auxFFT.resize(3,3);
    DIRECT_A2D_ELEM(auxFFT,0,0) = std::complex<double>(2.77778,0);
    DIRECT_A2D_ELEM(auxFFT,0,1) = std::complex<double>(-0.0555556,0.096225);
    DIRECT_A2D_ELEM(auxFFT,0,2) = std::complex<double>(-0.0555556,-0.096225);

    DIRECT_A2D_ELEM(auxFFT,1,0) = std::complex<double>(-0.388889,0.673575) ;
    DIRECT_A2D_ELEM(auxFFT,1,1) = std::complex<double>(-0.388889,-0.096225);
    DIRECT_A2D_ELEM(auxFFT,1,2) = std::complex<double>(-0.0555556,-0.288675);

    DIRECT_A2D_ELEM(auxFFT,2,0) = std::complex<double>(-0.388889,-0.673575) ;
    DIRECT_A2D_ELEM(auxFFT,2,1) = std::complex<double>(-0.0555556,0.288675) ;
    DIRECT_A2D_ELEM(auxFFT,2,2) = std::complex<double>(-0.388889,0.096225) ;
    EXPECT_EQ(FFT1,auxFFT);
    transformer1.cleanup();
}

TEST_F( FftwTest, fft_IDX2DIGFREQ)
{
	double w;
    FFT_IDX2DIGFREQ(0,128,w);
    EXPECT_EQ(0,w);
    FFT_IDX2DIGFREQ(1,128,w);
    EXPECT_EQ(1.0/128.0,w);
    FFT_IDX2DIGFREQ(64,128,w);
    EXPECT_EQ(0.5,w);
    FFT_IDX2DIGFREQ(65,128,w);
    EXPECT_EQ(-63.0/128.0,w);
    FFT_IDX2DIGFREQ(127,128,w);
    EXPECT_EQ(-1.0/128.0,w);

    FFT_IDX2DIGFREQ(0,129,w);
    EXPECT_EQ(0,w);
    FFT_IDX2DIGFREQ(1,129,w);
    EXPECT_EQ(1.0/129.0,w);
    FFT_IDX2DIGFREQ(64,129,w);
    EXPECT_EQ(64.0/129.0,w);
    FFT_IDX2DIGFREQ(65,129,w);
    EXPECT_EQ(-64.0/129.0,w);
    FFT_IDX2DIGFREQ(128,129,w);
    EXPECT_EQ(-1.0/129.0,w);

    size_t i=255;
    size_t dim=256;
    FFT_IDX2DIGFREQ(i,dim,w);
    EXPECT_EQ(-1.0/256.0,w);
}
*/
GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
