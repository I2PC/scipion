/***************************************************************************
 *
 * Authors:    Roberto Marabini       roberto@cnb.csic.es   
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

#include "movie_filter_dose.h"
#include <data/metadata_extension.h>
#include <data/xmipp_fftw.h>
#include <data/filters.h>
#include <typeinfo>

#define OUTSIDE_WRAP 0
#define OUTSIDE_AVG 1
#define OUTSIDE_VALUE 2

// Read arguments ==========================================================
void ProgMovieFilterDose::readParams() {
	user_supplied_input_filename = getParam("-i");
	user_supplied_output_filename = getParam("-o");
	user_supplied_first_frame = getIntParam("--frameRange", 0);
	user_supplied_last_frame = getIntParam("--frameRange", 1);
	pixel_size = getDoubleParam("--sampling");
	dose_per_frame = getDoubleParam("--dosePerFrame");
	acceleration_voltage = getDoubleParam("--accVoltage");
	initVoltage(acceleration_voltage);
	pre_exposure_amount = getDoubleParam("--preExposure");
	//restore_power = checkParam("--restoreNoiPow");

}

// Show ====================================================================
void ProgMovieFilterDose::show() {
	if (!verbose)
		return;
	std::cout << "Input movie:           " << user_supplied_input_filename
			<< std::endl << "Output movie:          "
			<< user_supplied_output_filename << std::endl
			<< "Frame range alignment: " << user_supplied_first_frame << " "
			<< user_supplied_last_frame << std::endl
			<< "Sampling:              " << pixel_size << std::endl
			<< "Dose per Frame (e/A^2):      " << dose_per_frame << std::endl
			<< "Acceleration_Voltage (kV):   " << acceleration_voltage
			<< std::endl << "Pre-exposure Amount (e/A^2): "
			<< pre_exposure_amount << std::endl
			<< "Restore Power after Filter:  " << restore_power << std::endl;
}

// usage ===================================================================
void ProgMovieFilterDose::defineParams() {
	addUsageLine("Align a set of frames by cross-correlation of the frames");
	addParamsLine("   -i <movie>                  : input movie");
	addParamsLine("  [-o <fn=\"out.mrcs\">]        : output filtered movie.");
	addParamsLine(
			"  [--frameRange <n0=-1> <nF=-1>]  : First and last frame to align, frame numbers start at 0");
	addParamsLine("  [--sampling <Ts=1>]          : Sampling rate (A/pixel)");
	addParamsLine("  [--dosePerFrame <dose=2>]    : Dose per frame (e/A^2)");
	addParamsLine(
			"  [--accVoltage <voltage=300>] : Acceleration voltage (kV) min_value=200.0e0,max_value=300.0e0)");
	addParamsLine("  [--preExposure <preExp=0>]   : P (e/A^2)");
//	addParamsLine(
//			"  [--restoreNoiPow]            : Restore noise power after filtering? (default=false)");
	addExampleLine("A typical example", false);
	addExampleLine("xmipp_cxzczxczccxztobedone");
}

double ProgMovieFilterDose::doseFilter(double dose_at_end_of_frame,
		double critical_dose) {
	return exp((-0.5 * dose_at_end_of_frame) / critical_dose);
}

double ProgMovieFilterDose::criticalDose(double spatial_frequency) {
	return ((critical_dose_a * pow(spatial_frequency, critical_dose_b))
			+ critical_dose_c) * voltage_scaling_factor;
}
double ProgMovieFilterDose::optimalDoseGivenCriticalDose(double critical_dose) {
	return 2.51284 * critical_dose;
}

void ProgMovieFilterDose::init(void) {
	/**  Set up the critical curve function (summovie) */
	critical_dose_a = 0.24499;
	critical_dose_b = -1.6649;
	critical_dose_c = 2.8141;
	// init with a very large number
	critical_dose_at_dc = std::numeric_limits<double>::max() * 0.001;

}
ProgMovieFilterDose::ProgMovieFilterDose(void) {
	init();
}

void ProgMovieFilterDose::initVoltage(double accelerationVoltage) {
		if (accelerationVoltage < 301 && accelerationVoltage > 299.) {
			acceleration_voltage = 300.0;
			voltage_scaling_factor = 1.0;
		} else if (accelerationVoltage < 201.0 && accelerationVoltage > 199.0) {
			acceleration_voltage = 200.0;
			voltage_scaling_factor = 0.8;
		} else
			REPORT_ERROR(ERR_ARG_INCORRECT,
					"Bad acceleration voltage (must be 200 or 300 kV");

	}

ProgMovieFilterDose::ProgMovieFilterDose(double accelerationVoltage) {
	init();
	initVoltage(accelerationVoltage);
}


void ProgMovieFilterDose::applyDoseFilterToImage(
		int Ydim, int Xdim,
		const MultidimArray<std::complex<double> > &FFT1,
		const double dose_start, const double dose_finish) {
	double * ptrv = (double *) MULTIDIM_ARRAY(FFT1);
	std::complex<double> complex_zero = 0;
	double x, y;
	double xx, yy;
	double current_critical_dose, current_optimal_dose;
	int sizeX_2 = Xdim / 2;
	double ixsize = 1.0 / Xdim;
	int sizeY_2 = Ydim / 2;
	double iysize = 1.0 / Ydim;
	for (long int i = 0; i < YSIZE(FFT1); i++) { //i->y,j->x xmipp convention is oposite summove
		FFT_IDX2DIGFREQ_FAST(i, Ydim, sizeY_2, iysize, y);
		yy = y * y;
		for (long int j = 0; j < XSIZE(FFT1); j++) {
			FFT_IDX2DIGFREQ_FAST(j, Xdim, sizeX_2, ixsize, x);
#ifdef DEBUG_FFT_FREQ
			std::cerr << "(j+1,i+1)(x,y)=" << j +1 << "," << i + 1<< " "
					                       << x    << " " << y << " = "
										   << DIRECT_A2D_ELEM(FFT1, i, j) << std::endl;
#endif
			if (i == 0 && j == 0)
				//ROB: I do not understand this step
				//It forces the origin to 0
				//Why not keep the value?
				current_critical_dose = critical_dose_at_dc;
			else
				current_critical_dose = criticalDose(
						sqrt(x * x + yy) / pixel_size);
		    current_optimal_dose = optimalDoseGivenCriticalDose(
					current_critical_dose);
			if (abs(dose_finish - current_optimal_dose)
					< abs(dose_start - current_optimal_dose)) {
				DIRECT_A2D_ELEM(FFT1, i, j ) *= doseFilter(dose_finish,
						current_critical_dose);
			} else
				DIRECT_A2D_ELEM(FFT1, i, j ) = complex_zero;

		}
	}

}
void ProgMovieFilterDose::run() {
	show();
	//read movie
	MetaData movie;
	size_t Xdim, Ydim, Zdim, Ndim;

	//if input is an stack create a metadata.
	if (user_supplied_input_filename.isMetaData())
		movie.read(user_supplied_input_filename);
	else {
		ImageGeneric movieStack;
		movieStack.read(user_supplied_input_filename, HEADER);
		movieStack.getDimensions(Xdim, Ydim, Zdim, Ndim);
		if (user_supplied_input_filename.getExtension() == "mrc" and Ndim == 1)
			Ndim = Zdim;
		size_t id;
		FileName fn;
		for (size_t i = 0; i < Ndim; i++) {
			id = movie.addObject();
			fn.compose(i + FIRST_IMAGE, user_supplied_input_filename);
			movie.setValue(MDL_IMAGE, fn, id);
		}
	}
	//movie.write("/tmp/movie.md");
	if (user_supplied_first_frame < 0)
		user_supplied_first_frame = 0;
	if (user_supplied_last_frame < 0)
		user_supplied_last_frame = movie.size() - 1;

	//Fourier transform
	//apply filter frame, by frame
	if (verbose) {
		std::cout << "Computing Fourier transform of frames ..." << std::endl;
		init_progress_bar(movie.size());
	}
	size_t n = 0;
	FourierTransformer transformer;
	FileName fnFrame;
	FileName fnOutFrame;
	Image<double> frame;
	int mode = WRITE_OVERWRITE;

	MultidimArray<std::complex<double> > FFT1;
	FOR_ALL_OBJECTS_IN_METADATA(movie)
	{
		if (n >= user_supplied_first_frame && n <= user_supplied_last_frame) {
			movie.getValue(MDL_IMAGE, fnFrame, __iter.objId);
			frame.read(fnFrame);
			// Now do the Fourier transform and filter
			transformer.FourierTransform(frame(), FFT1, false);
#ifdef PRINT_ORIGINAL_IMAGE
			std::cerr << "n= (original image)" << n << std::endl << frame() << std::endl;
#endif
#undef DEBUG
//#define DEBUG
#ifdef DEBUG
			//same order than summovie
			Image<double> VFourierMag;
			//FFT_magnitude(FFT1,VFourierMag());
			VFourierMag().resizeNoCopy(FFT1);
			double * ptrv=(double *)MULTIDIM_ARRAY(FFT1);
			int counter = 0;
			double x,y;
			int m1sizeX = XSIZE(frame()); int m1sizeY = YSIZE(frame());
			int sizeX_2 = m1sizeX/2; double ixsize = 1.0/m1sizeX;
			int sizeY_2 = m1sizeY/2; double iysize = 1.0/m1sizeY;
			for (long int i=0; i<YSIZE(FFT1); i++) {
				FFT_IDX2DIGFREQ_FAST(i, m1sizeY, sizeY_2,iysize , y);
				for (long int j=0; j<XSIZE(FFT1); j++) {    //i=y, j=x
					FFT_IDX2DIGFREQ_FAST(j, m1sizeX, sizeX_2,ixsize , x);
					std::cerr << "(x,y)=" << j +1 << "," << i + 1<< " = " <<
					DIRECT_A2D_ELEM(FFT1, i, j) <<
					x*x << " " << y*y << " " << (x*x+y*y) << std::endl;
				}
			}
			exit(1);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FFT1)
			{
				counter ++;
				double re=*ptrv;
				double im=*(ptrv+1);
				std::cerr << "DEBUG_JM: n re im: " << counter << " " << re << " " << im << std::endl;
				DIRECT_MULTIDIM_ELEM(VFourierMag(), n) =
				log(sqrt(re*re+im*im)+1);
				ptrv+=2;
			}
			VFourierMag.write("/tmp/frame_mag_1.mrc");
			frame.write("/tmp/frame_1.mrc");
#endif
#undef DEBUG
#ifdef PRINTFFT
			std::cerr << "fourier before filter " << n << std::endl
					<< FFT1 << std::endl;
#endif
			applyDoseFilterToImage((int)YSIZE(frame()), (int)XSIZE(frame()), FFT1,
					(n * dose_per_frame) + pre_exposure_amount,    //dose_star
					((n + 1) * dose_per_frame) + pre_exposure_amount //dose_end
			);
#ifdef PRINTFFTAFTERFILTER
			std::cerr << "Fourier  after filter"<< std::endl << FFT1 << std::endl;
#endif
			transformer.inverseFourierTransform();
#ifdef PRINTREALIMAGEAFTERDOEFILTERING
			std::cerr << "fourier matrix after inverse (very likely garbage)" << std::endl
					<< FFT1 << std::endl;
			std::cerr << "image after filter "  << std::endl
					<< frame() << std::endl;
#endif
			frame.write(user_supplied_output_filename,n+FIRST_IMAGE,true,mode);
			mode = WRITE_APPEND;
		}
		++n;
		if (verbose)
			progress_bar(n);
	}
	if (verbose)
		progress_bar(movie.size());

}
