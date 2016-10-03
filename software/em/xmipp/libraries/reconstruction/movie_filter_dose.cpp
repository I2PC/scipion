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
void ProgMovieFilterDose::readParams()
{
	user_supplied_input_filename = getParam("-i");
	user_supplied_output_filename = getParam("-o");
	user_supplied_first_frame = getIntParam("--frameRange",0);
	user_supplied_last_frame = getIntParam("--frameRange",1);
	pixel_size = getDoubleParam("--sampling");
	dose_per_frame = getDoubleParam("--dosePerFrame");
	acceleration_voltage = getDoubleParam("--accVoltage");
	pre_exposure_amount = getDoubleParam("--preExposure");
	restore_power = checkParam("--restoreNoiPow");

}

// Show ====================================================================
void ProgMovieFilterDose::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input movie:           " << user_supplied_input_filename << std::endl
    << "Output movie:          " << user_supplied_output_filename              << std::endl
    << "Frame range alignment: " << user_supplied_first_frame << " "
	                             << user_supplied_last_frame  << std::endl
    << "Sampling:              " << pixel_size                << std::endl
    << "Dose per Frame (e/A^2):      " << dose_per_frame            << std::endl
    << "Acceleration_Voltage (kV):   " << acceleration_voltage                << std::endl
    << "Pre-exposure Amount (e/A^2): " << pre_exposure_amount                << std::endl
    << "Restore Power after Filter:  " << restore_power                << std::endl
    ;
}

// usage ===================================================================
void ProgMovieFilterDose::defineParams()
{
    addUsageLine("Align a set of frames by cross-correlation of the frames");
    addParamsLine("   -i <movie>                  : input movie");
    addParamsLine("  [-o <fn=\"out.mrcs\">]        : output filtered movie.");
    addParamsLine("  [--frameRange <n0=-1> <nF=-1>]  : First and last frame to align, frame numbers start at 0");
    addParamsLine("  [--sampling <Ts=1>]          : Sampling rate (A/pixel)");
    addParamsLine("  [--dosePerFrame <dose=2>]    : Dose per frame (e/A^2)");
    addParamsLine("  [--accVoltage <voltage=300>] : Acceleration voltage (kV) min_value=200.0e0,max_value=300.0e0)");
    addParamsLine("  [--preExposure <preExp=0>]   : P (e/A^2)");
    addParamsLine("  [--restoreNoiPow]            : Restore noise power after filtering? (default=false)");
    addExampleLine("A typical example",false);
    addExampleLine("xmipp_cxzczxczccxztobedone");
}


double ProgMovieFilterDose::doseFilter(double dose_at_end_of_frame,
		                               double critical_dose)
{
	return exp((-0.5 * dose_at_end_of_frame)/critical_dose);
}

double ProgMovieFilterDose::criticalDose(double spatial_frequency)
    {
     return (critical_dose_a * pow(spatial_frequency,critical_dose_b) +
    		        critical_dose_c) * voltage_scaling_factor;
    }
double ProgMovieFilterDose::OptimalDoseGivenCriticalDose(double critical_dose)
{
	 return 2.51284 * critical_dose;
}

void ProgMovieFilterDose::init(void)
{
	/**  Set up the critical curve function (summovie */
    critical_dose_a =  0.24499;
	critical_dose_b =  -1.6649;
	critical_dose_c =  2.8141;
	// init with a very large number
	critical_dose_at_dc = std::numeric_limits<double>::max() * 0.001;

}
ProgMovieFilterDose::ProgMovieFilterDose(void)
{
	init();
}

ProgMovieFilterDose::ProgMovieFilterDose(double accelerationVoltage)
{

    init();
	if (accelerationVoltage < 301 && accelerationVoltage > 299.){
        acceleration_voltage = 300.0;
        voltage_scaling_factor = 1.0e0;
	}
    else if (accelerationVoltage < 201.0 &&  accelerationVoltage > 199.0){
        acceleration_voltage = 200.0;
        voltage_scaling_factor = 0.80;
	}
    else
		REPORT_ERROR(ERR_ARG_INCORRECT,"Bad acceleration voltage (must be 200 or 300 kv");

}
void ProgMovieFilterDose::ApplyDoseFilterToImage(
		const MultidimArray<double> &frame,
		const MultidimArray< std::complex<double> > &FFT1,
		const double dose_start, const double dose_finish
	){
    double * ptrv=(double *)MULTIDIM_ARRAY(FFT1);
    std::complex<double> complex_zero=0;
	double x,y;
	double xx,yy;
	double current_critical_dose, current_optimal_dose;
	int m1sizeX = XSIZE(frame);  int m1sizeY = YSIZE(frame);
	int sizeX_2 = m1sizeX/2; double ixsize = 1.0/m1sizeX;
	int sizeY_2 = m1sizeY/2; double iysize = 1.0/m1sizeY;
	for (long int i=0; i<YSIZE(FFT1); i++){//i->y,j->x xmipp convention is oposite summove
        FFT_IDX2DIGFREQ_FAST(i, m1sizeY, sizeY_2,iysize , y);
        yy = y*y;
		for (long int j=0; j<XSIZE(FFT1); j++){
            FFT_IDX2DIGFREQ_FAST(j, m1sizeX, sizeX_2,ixsize , x);
            if (i == 1 && j == 1)
            	current_critical_dose = critical_dose_at_dc;
            else
            	current_critical_dose = criticalDose(sqrt(x*x+yy)/pixel_size);
			current_optimal_dose = OptimalDoseGivenCriticalDose(current_critical_dose);
            if (abs(dose_finish-current_optimal_dose) <
            	abs(dose_start-current_optimal_dose))
            	DIRECT_A2D_ELEM(FFT1, i, j ) *= doseFilter(dose_finish,current_critical_dose);
            else
            	DIRECT_A2D_ELEM(FFT1, i, j ) = complex_zero;

    	    std::cerr << "(x,y)=" << j +1 << "," << i + 1<< " = " <<
    		DIRECT_A2D_ELEM(FFT1, i, j) <<
    	    x*x << " " << y*y  << " " << (x*x+y*y) << std::endl;
	     }
	}


}
void ProgMovieFilterDose::run()
{
    show();
    //read movie
    MetaData movie;
    size_t Xdim, Ydim, Zdim, Ndim;

    //if input is an stack create a metadata.
    if (user_supplied_input_filename.isMetaData())
        movie.read(user_supplied_input_filename);
    else
    {
        ImageGeneric movieStack;
        movieStack.read(user_supplied_input_filename,HEADER);
        movieStack.getDimensions(Xdim,Ydim,Zdim,Ndim);
        if (user_supplied_input_filename.getExtension()=="mrc" and Ndim ==1)
            Ndim = Zdim;
        size_t id;
        FileName fn;
        for (size_t i=0;i<Ndim;i++)
        {
            id = movie.addObject();
            fn.compose(i+FIRST_IMAGE,user_supplied_input_filename);
            movie.setValue(MDL_IMAGE, fn, id);
        }
    }
    //movie.write("/tmp/movie.md");
    if (user_supplied_first_frame<0)
    	user_supplied_first_frame=0;
    if (user_supplied_last_frame<0)
    	user_supplied_last_frame=movie.size()-1;

    //Fourier transform
    //apply filter frame, by frame
	if (verbose)
	{
		std::cout << "Computing Fourier transform of frames ..." << std::endl;
		init_progress_bar(movie.size());
	}
	int n=0;
	FourierTransformer transformer;
	FileName fnFrame;
	Image<double> frame;

	FOR_ALL_OBJECTS_IN_METADATA(movie)
	{
		if (n>=user_supplied_first_frame && n<=user_supplied_last_frame)
		{
			movie.getValue(MDL_IMAGE,fnFrame,__iter.objId);
			frame.read(fnFrame);
			// Now do the Fourier transform and filter
			MultidimArray< std::complex<double> > FFT1;
			transformer.FourierTransform(frame(),FFT1,true);
			// apply dose
			ApplyDoseFilterToImage(frame(),FFT1,
                    ( n*dose_per_frame   ) + pre_exposure_amount,//dose_star
                    ((n+1)*dose_per_frame) + pre_exposure_amount //dose_end
	        		);


//#define DEBUG
#ifdef DEBUG
			//same order than summovie
			Image<double> VFourierMag;
			//FFT_magnitude(FFT1,VFourierMag());
			VFourierMag().resizeNoCopy(FFT1);
		    double * ptrv=(double *)MULTIDIM_ARRAY(FFT1);
		    int counter = 0;
			double x,y;
			int m1sizeX = XSIZE(frame());  int m1sizeY = YSIZE(frame());
			int sizeX_2 = m1sizeX/2; double ixsize = 1.0/m1sizeX;
			int sizeY_2 = m1sizeY/2; double iysize = 1.0/m1sizeY;
			for (long int i=0; i<YSIZE(FFT1); i++){
	            FFT_IDX2DIGFREQ_FAST(i, m1sizeY, sizeY_2,iysize , y);
				for (long int j=0; j<XSIZE(FFT1); j++){//i=y, j=x
	                FFT_IDX2DIGFREQ_FAST(j, m1sizeX, sizeX_2,ixsize , x);
		    	    std::cerr << "(x,y)=" << j +1 << "," << i + 1<< " = " <<
		    		DIRECT_A2D_ELEM(FFT1, i, j) <<
		    	    x*x << " " << y*y  << " " << (x*x+y*y) << std::endl;
			     }
			}
		    //exit(1);
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
			//frameFourierVec.push_back(FFT1);
			///call filter here
		}
		++n;
		if (verbose)
			progress_bar(n);
	}
	if (verbose)
		progress_bar(movie.size());
/*
		// Free useless memory
		filter.clear();
		reducedFrame.clear();
		croppedFrame.clear();
		frame.clear();

		// Now compute all shifts
		size_t N=frameFourier.size();
		Matrix2D<double> A(N*(N-1)/2,N-1);
		Matrix1D<double> bX(N*(N-1)/2), bY(N*(N-1)/2);
		if (verbose)
			std::cout << "Computing shifts between frames ..." << std::endl;
		int idx=0;
		MultidimArray<double> Mcorr;
		Mcorr.resizeNoCopy(newYdim,newXdim);
		Mcorr.setXmippOrigin();
		CorrelationAux aux;
		for (size_t i=0; i<N-1; ++i)
		{
			for (size_t j=i+1; j<N; ++j)
			{
				bestShift(*frameFourier[i],*frameFourier[j],Mcorr,bX(idx),bY(idx),aux,NULL,maxShift);
				if (verbose)
					std::cerr << "Frame " << i+nfirst << " to Frame " << j+nfirst << " -> (" << bX(idx) << "," << bY(idx) << ")\n";
				for (int ij=i; ij<j; ij++)
					A(idx,ij)=1;

				idx++;
			}
			delete frameFourier[i];
		}

		// Finally solve the equation system
		Matrix1D<double> shiftX, shiftY, ex, ey;
		WeightedLeastSquaresHelper helper;
		helper.A=A;
		helper.w.initZeros(VEC_XSIZE(bX));
		helper.w.initConstant(1);

		int it=0;
		double mean, varbX, varbY;
		bX.computeMeanAndStddev(mean,varbX);
		varbX*=varbX;
		bY.computeMeanAndStddev(mean,varbY);
		varbY*=varbY;
		if (verbose)
			std::cout << "\nSolving for the shifts ...\n";
		do
		{
			// Solve the equation system
			helper.b=bX;
			weightedLeastSquares(helper,shiftX);
			helper.b=bY;
			weightedLeastSquares(helper,shiftY);

			// Compute residuals
			ex=bX-A*shiftX;
			ey=bY-A*shiftY;

			// Compute R2
			double mean, vareX, vareY;
			ex.computeMeanAndStddev(mean,vareX);
			vareX*=vareX;
			ey.computeMeanAndStddev(mean,vareY);
			vareY*=vareY;
			double R2x=1-vareX/varbX;
			double R2y=1-vareY/varbY;
			if (verbose)
				std::cout << "Iteration " << it << " R2x=" << R2x << " R2y=" << R2y << std::endl;

			// Identify outliers
			double oldWeightSum=helper.w.sum();
			double stddeveX=sqrt(vareX);
			double stddeveY=sqrt(vareY);
			FOR_ALL_ELEMENTS_IN_MATRIX1D(ex)
			if (fabs(VEC_ELEM(ex,i))>3*stddeveX || fabs(VEC_ELEM(ey,i))>3*stddeveY)
				VEC_ELEM(helper.w,i)=0.0;
			double newWeightSum=helper.w.sum();
			if (newWeightSum==oldWeightSum)
			{
				std::cout << "No outlier found\n";
				break;
			}
			else
				std::cout << "Found " << (int)(oldWeightSum-newWeightSum) << " outliers\n";

			it++;
		}
		while (it<solverIterations);

		// Choose reference image as the minimax of shifts
		double worstShiftEver=1e38;
		for (int iref=0; iref<N; ++iref)
		{
			double worstShift=-1;
			for (int j=0; j<N; ++j)
			{
				double totalShiftX, totalShiftY;
				computeTotalShift2(iref, j, shiftX, shiftY,totalShiftX, totalShiftY);
				if (fabs(totalShiftX)>worstShift)
					worstShift=fabs(totalShiftX);
				if (fabs(totalShiftX)>worstShift)
					worstShift=fabs(totalShiftX);
			}
			if (worstShift<worstShiftEver)
			{
				worstShiftEver=worstShift;
				bestIref=iref;
			}
		}
		if (verbose)
			std::cout << "Reference frame: " << bestIref+1+nfirst << std::endl;

	    // Compute shifts
	    int j=0;
	    n=0;
	    FOR_ALL_OBJECTS_IN_METADATA(movie)
	    {
	        if (n>=nfirst && n<=nlast)
	        {
	            double totalShiftX, totalShiftY;
	            computeTotalShift2(bestIref, j, shiftX, shiftY,XX(shift), YY(shift));
	            shift/=sizeFactor;
	            shift*=-1;
	            movie.setValue(MDL_SHIFT_X,XX(shift),__iter.objId);
	            movie.setValue(MDL_SHIFT_Y,YY(shift),__iter.objId);
	            j++;
	            movie.setValue(MDL_ENABLED,1,__iter.objId);
	        }
	        else
	        {
	            movie.setValue(MDL_ENABLED,-1,__iter.objId);
	            movie.setValue(MDL_SHIFT_X,0.0,__iter.objId);
	            movie.setValue(MDL_SHIFT_Y,0.0,__iter.objId);
	        }
	        movie.setValue(MDL_WEIGHT,1.0,__iter.objId);
	        n++;
	    }
    }
    else
    {
    	if (!movie.containsLabel(MDL_SHIFT_X))
    	{
    		movie.addLabel(MDL_SHIFT_X);
    		movie.addLabel(MDL_SHIFT_Y);
    		movie.fillConstant(MDL_SHIFT_X,"0.0");
    		movie.fillConstant(MDL_SHIFT_Y,"0.0");
    	}
    }

    // Apply shifts and compute average
	int n=0;
    int j=0;
	size_t N=0;
    FOR_ALL_OBJECTS_IN_METADATA(movie)
    {
        if (n>=nfirstSum && n<=nlastSum)
        {
            movie.getValue(MDL_IMAGE,fnFrame,__iter.objId);
            movie.getValue(MDL_SHIFT_X,XX(shift),__iter.objId);
            movie.getValue(MDL_SHIFT_Y,YY(shift),__iter.objId);
            std::cout << fnFrame << " shiftX=" << XX(shift) << " shiftY=" << YY(shift) << std::endl;

            frame.read(fnFrame);
            if (XSIZE(dark())>0)
            	croppedFrame()-=dark();
            if (XSIZE(gain())>0)
            	croppedFrame()*=gain();
            if (yDRcorner!=-1)
                frame().window(croppedFrame(), yLTcorner, xLTcorner, yDRcorner, xDRcorner);
            else
            	croppedFrame()=frame();
            if (bin>0)
            {
            	scaleToSizeFourier(1,floor(YSIZE(croppedFrame())/bin),floor(XSIZE(croppedFrame())/bin),croppedFrame(),reducedFrame());
            	shift/=bin;
            	croppedFrame()=reducedFrame();
            }

            if (fnAligned!="" || fnAvg!="")
            {
            	if (outsideMode == OUTSIDE_WRAP)
            		translate(BsplineOrder,shiftedFrame(),croppedFrame(),shift,WRAP);
            	else if (outsideMode == OUTSIDE_VALUE)
            		translate(BsplineOrder,shiftedFrame(),croppedFrame(),shift,DONT_WRAP, outsideValue);
            	else
            		translate(BsplineOrder,shiftedFrame(),croppedFrame(),shift,DONT_WRAP, croppedFrame().computeAvg());
                if (fnAligned!="")
                    shiftedFrame.write(fnAligned,j+1,true,WRITE_REPLACE);
                if (fnAvg!="")
                {
                    if (j==0)
                        averageMicrograph()=shiftedFrame();
                    else
                        averageMicrograph()+=shiftedFrame();
                    N++;
                }
            }

            j++;
        }
        n++;
    }
    if (fnAvg!="")
    {
        averageMicrograph()/=N;
        averageMicrograph.write(fnAvg);
    }
    movie.write((FileName)"frameShifts@"+fnOut);
    if (bestIref>=0)
    {
    	MetaData mdIref;
    	mdIref.setValue(MDL_REF,nfirst+bestIref,mdIref.addObject());
    	mdIref.write((FileName)"referenceFrame@"+fnOut,MD_APPEND);
    }
    */
}
