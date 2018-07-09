/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2016)
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

#include "resolution_monogenic_signal.h"
//#define DEBUG
//#define DEBUG_MASK

void ProgMonogenicSignalRes::readParams()
{
	fnVol = getParam("--vol");
	fnVol2 = getParam("--vol2");
	fnMeanVol = getParam("--meanVol");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	fnMaskOut = getParam("--mask_out");
	fnchim = getParam("--chimera_volume");
	sampling = getDoubleParam("--sampling_rate");
	R = getDoubleParam("--volumeRadius");
	minRes = getDoubleParam("--minRes");
	maxRes = getDoubleParam("--maxRes");
	fnSym = getParam("--sym");
	freq_step = getDoubleParam("--step");
	trimBound = getDoubleParam("--trimmed");
	exactres = checkParam("--exact");
	noiseOnlyInHalves = checkParam("--noiseonlyinhalves");
	fnSpatial = getParam("--filtered_volume");
	significance = getDoubleParam("--significance");
	fnMd = getParam("--md_outputdata");
	automaticMode = checkParam("--automatic");
	nthrs = getIntParam("--threads");
}


void ProgMonogenicSignalRes::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file=\"\">   : Input volume");
	addParamsLine("  [--mask <vol_file=\"\">]  : Mask defining the macromolecule");
	addParamsLine("                          :+ If two half volume are given, the noise is estimated from them");
	addParamsLine("                          :+ Otherwise the noise is estimated outside the mask");
	addParamsLine("  [--mask_out <vol_file=\"\">]  : sometimes the provided mask is not perfect, and contains voxels out of the particle");
	addParamsLine("                          :+ Thus the algorithm calculated a tight mask to the volume");
	addParamsLine("  [--vol2 <vol_file=\"\">]: Half volume 2");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--meanVol <vol_file=\"\">]: Mean volume of half1 and half2 (only it is neccesary the two haves are used)");
	addParamsLine("  --sym <symmetry>: Symmetry (c1, c2, c3,..d1, d2, d3,...)");
	addParamsLine("  [--chimera_volume <output=\"Chimera_resolution_volume.vol\">]: Local resolution volume for chimera viewer (in Angstroms)");
	addParamsLine("  [--sampling_rate <s=1>]   : Sampling rate (A/px)");
	addParamsLine("                            : Use -1 to disable this option");
	addParamsLine("  [--volumeRadius <s=100>]   : This parameter determines the radius of a sphere where the volume is");
	addParamsLine("  [--step <s=0.25>]       : The resolution is computed at a number of frequencies between mininum and");
	addParamsLine("                            : maximum resolution px/A. This parameter determines that number");
	addParamsLine("  [--minRes <s=30>]         : Minimum resolution (A)");
	addParamsLine("  [--maxRes <s=1>]          : Maximum resolution (A)");
	addParamsLine("  [--trimmed <s=0.5>]       : Trimming percentile");
	addParamsLine("  [--exact]                 : The search for resolution will be exact (slower) of approximated (fast).");
	addParamsLine("                            : Usually there are no difference between both in the resolution map.");
	addParamsLine("  [--noiseonlyinhalves]     : The noise estimation is only performed inside the mask");
	addParamsLine("  [--filtered_volume <vol_file=\"\">]       : The input volume is locally filtered at local resolutions.");
	addParamsLine("  [--significance <s=0.95>]    : The level of confidence for the hypothesis test.");
	addParamsLine("  [--md_outputdata <file=\".\">]  : It is a control file. The provided mask can contain voxels of noise.");
	addParamsLine("                                  : Moreover, voxels inside the mask cannot be measured due to an unsignificant");
	addParamsLine("                                  : SNR. Thus, a new mask is created. This metadata file, shows, the number of");
	addParamsLine("                                  : voxels of the original mask, and the created mask");
	addParamsLine("  [--automatic]                   : Resolution range is not neccesary provided");
	addParamsLine("  [--threads <s=4>]               : Number of threads");
}


void ProgMonogenicSignalRes::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;
	Image<double> V;
	if ((fnVol !="") && (fnVol2 !=""))
	{
		Image<double> V1, V2;
		V1.read(fnVol);
		V2.read(fnVol2);
		V()=0.5*(V1()+V2());
		V.write(fnMeanVol);
		halfMapsGiven = true;
	}
	else{
	    V.read(fnVol);
	    halfMapsGiven = false;
	}
	V().setXmippOrigin();


	transformer_inv.setThreadsNumber(nthrs);

	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();
	VRiesz.resizeNoCopy(inputVol);

	if (fnSpatial!="")
		VresolutionFiltered().initZeros(V());

	transformer.FourierTransform(inputVol, fftV);
	iu.initZeros(fftV);

	// Calculate u and first component of Riesz vector
	double uz, uy, ux, uz2, u2, uz2y2;
	long n=0;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
		uz2=uz*uz;

		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
			uz2y2=uz2+uy*uy;

			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
				u2=uz2y2+ux*ux;
				if ((k != 0) || (i != 0) || (j != 0))
					DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
				else
					DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
				++n;
			}
		}
	}
	#ifdef DEBUG
	Image<double> saveiu;
	saveiu = 1/iu;
	saveiu.write("iu.vol");
	#endif

	// Prepare low pass filter
	lowPassFilter.FilterShape = RAISED_COSINE;
	lowPassFilter.raised_w = 0.01;
	lowPassFilter.do_generate_3dmask = false;
	lowPassFilter.FilterBand = LOWPASS;

	// Prepare mask
	MultidimArray<int> &pMask=mask();

	if (fnMask != "")
	{
		mask.read(fnMask);
		mask().setXmippOrigin();
	}
	else
	{
		std::cout << "Error: a mask ought to be provided" << std::endl;
		exit(0);
	}

	NVoxelsOriginalMask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
	{
		if (A3D_ELEM(pMask, k, i, j) == 1)
			++NVoxelsOriginalMask;
		if (i*i+j*j+k*k > R*R)
			A3D_ELEM(pMask, k, i, j) = -1;
	}

	#ifdef DEBUG_MASK
	mask.write("mask.vol");
	#endif

	if (halfMapsGiven)
	{
		Image<double> V1, V2;
		V1.read(fnVol);
		V2.read(fnVol2);

		V1()-=V2();
		V1()/=2;
		fftN=new MultidimArray< std::complex<double> >;
		FourierTransformer transformer2;
		#ifdef DEBUG
		  V1.write("diff_volume.vol");
		#endif
		transformer2.FourierTransform(V1(), *fftN);
	}
	else
	{
		fftN=&fftV;
	}
	V.clear();

	double u;
	int size_fourier = ZSIZE(fftV);
	freq_fourier.initZeros(size_fourier);

	int size = ZSIZE(pMask);

	VEC_ELEM(freq_fourier,0) = 1e-38;
	for(size_t k=0; k<size_fourier; ++k)
	{
		FFT_IDX2DIGFREQ(k,size, u);
		VEC_ELEM(freq_fourier,k) = u;
	}
}


void ProgMonogenicSignalRes::amplitudeMonogenicSignal3D(MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, double freqL, MultidimArray<double> &amplitude, int count, FileName fnDebug)
{
	fftVRiesz.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(freq-freqH);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
	{
		double iun=DIRECT_MULTIDIM_ELEM(iu,n);
		double un=1.0/iun;
		if (freqH<=un && un<=freq)
		{
			//double H=0.5*(1+cos((un-w1)*ideltal));
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
		} else if (un>freq)
		{
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
		}
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	#ifdef DEBUG
	Image<double> filteredvolume;
	filteredvolume = VRiesz;
	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
	#endif

	if (fnSpatial!="")
		Vfiltered()=VRiesz;

	#ifdef DEBUG
	FileName iternumber;
	iternumber = formatString("_Volume_%i.vol", count);
	Image<double> saveImg2;
	saveImg2() = VRiesz;
	  if (fnDebug.c_str() != "")
	  {
		saveImg2.write(fnDebug+iternumber);
	  }
	saveImg2.clear(); 
	#endif

	amplitude.resizeNoCopy(VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	double uz, uy, ux;
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				ux = VEC_ELEM(freq_fourier,j);
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second and third components of Riesz vector
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		uz = VEC_ELEM(freq_fourier,k);
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			uy = VEC_ELEM(freq_fourier,i);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = uz*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	transformer_inv.inverseFourierTransform(fftVRiesz_aux, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
	}
	#ifdef DEBUG
	if (fnDebug.c_str() != "")
	{
	Image<double> saveImg;
	saveImg = amplitude;
	iternumber = formatString("_Amplitude_%i.vol", count);
	saveImg.write(fnDebug+iternumber);
	saveImg.clear();
	}
	#endif // DEBUG
//
	// Low pass filter the monogenic amplitude
	transformer_inv.FourierTransform(amplitude, fftVRiesz, false);
	double raised_w = PI/(freqL-freq);

	n=0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
	{
		double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
//		std::cout << "un = " << un << "  freqL = " << freqL << " freq = " << freq << std::endl;
		if ((freqL)>=un && un>=freq)
		{
			DIRECT_MULTIDIM_ELEM(fftVRiesz,n) *= 0.5*(1 + cos(raised_w*(un-freq)));
		}
		else
		{
			if (un>freqL)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
			}
		}
	}
	transformer_inv.inverseFourierTransform();


	#ifdef DEBUG
	saveImg2 = amplitude;
	FileName fnSaveImg2;
	if (fnDebug.c_str() != "")
	{
		iternumber = formatString("_Filtered_Amplitude_%i.vol", count);
		saveImg2.write(fnDebug+iternumber);
	}
	saveImg2.clear(); 
	#endif // DEBUG
}


void ProgMonogenicSignalRes::firstMonoResEstimation(MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, double freqL, MultidimArray<double> &amplitude,
		int count, FileName fnDebug, double &mean_Signal, double &mean_noise, double &thresholdFirstEstimation)
{
	amplitudeMonogenicSignal3D(myfftV, freq, freqH, freqL, amplitude, count, fnDebug);

	double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
	MultidimArray<int> &pMask = mask();
	std::vector<double> noiseValues;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitude, n);
		if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
		{
			sumS  += amplitudeValue;
			++NS;
		}
		else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
		{
			noiseValues.push_back(amplitudeValue);
			sumN  += amplitudeValue;
			++NN;
		}
	}
	std::sort(noiseValues.begin(),noiseValues.end());

	mean_Signal = sumS/NS;
	mean_noise = sumN/NN;

	std::cout << "means_signal = " << mean_Signal << "mean_noise = " << mean_noise << std::endl;

	thresholdFirstEstimation = noiseValues[size_t(noiseValues.size()*significance)];


}


void ProgMonogenicSignalRes::postProcessingLocalResolutions(MultidimArray<double> &resolutionVol,
		std::vector<double> &list, MultidimArray<double> &resolutionChimera, double &cut_value, MultidimArray<int> &pMask)
{
	MultidimArray<double> resolutionVol_aux = resolutionVol;
	double last_resolution_2 = list[(list.size()-1)];

	double Nyquist;
	Nyquist = 2*sampling;

	// Count number of voxels with resolution
	size_t N=0;
	if (automaticMode)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
			if (DIRECT_MULTIDIM_ELEM(resolutionVol, n)>=(Nyquist)) //the value 0.001 is a tolerance
				++N;
	}
	else
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
			if (DIRECT_MULTIDIM_ELEM(resolutionVol, n)>=(last_resolution_2-0.001)) //the value 0.001 is a tolerance
				++N;
	}


	// Get all resolution values
	MultidimArray<double> resolutions(N);
	size_t N_iter=0;
	if (automaticMode)
	{
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n)>=Nyquist)
			DIRECT_MULTIDIM_ELEM(resolutions,N_iter++)=DIRECT_MULTIDIM_ELEM(resolutionVol, n);
	}
	else
	{
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n)>(last_resolution_2-0.001))
			DIRECT_MULTIDIM_ELEM(resolutions,N_iter++)=DIRECT_MULTIDIM_ELEM(resolutionVol, n);
	}
	// Sort value and get threshold
	std::sort(&A1D_ELEM(resolutions,0),&A1D_ELEM(resolutions,N));
	double filling_value = A1D_ELEM(resolutions, (int)(0.5*N)); //median value
	double trimming_value = A1D_ELEM(resolutions, (int)((1-cut_value)*N));

	double freq, res, init_res, last_res;

	init_res = list[0];
	last_res = list[(list.size()-1)];
	
	std::cout << "--------------------------" << std::endl;
	std::cout << "last_res = " << last_res << std::endl;

	resolutionChimera = resolutionVol;

	if (automaticMode)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
		{
			if (DIRECT_MULTIDIM_ELEM(pMask,n) == 0)
			{
			  DIRECT_MULTIDIM_ELEM(resolutionChimera, n) = filling_value;
			}
		}
	}
	else
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
		{
			if (DIRECT_MULTIDIM_ELEM(resolutionVol, n) < last_res)
			{
				if (DIRECT_MULTIDIM_ELEM(pMask, n) >=1)
				{
					DIRECT_MULTIDIM_ELEM(resolutionChimera, n) = filling_value;
					DIRECT_MULTIDIM_ELEM(resolutionVol, n) = filling_value;
				}
				else
				{
					DIRECT_MULTIDIM_ELEM(resolutionChimera, n) = filling_value;
					DIRECT_MULTIDIM_ELEM(pMask,n) = 0;
				}
			}
			if (DIRECT_MULTIDIM_ELEM(resolutionVol, n) > trimming_value)
			{
			  DIRECT_MULTIDIM_ELEM(pMask,n) = 0;
			  DIRECT_MULTIDIM_ELEM(resolutionVol, n) = filling_value;
			  DIRECT_MULTIDIM_ELEM(resolutionChimera, n) = filling_value;
			}
		}
	}
	//#ifdef DEBUG_MASK
	Image<int> imgMask;
	imgMask = pMask;
	imgMask.write(fnMaskOut);
	//#endif
}

void ProgMonogenicSignalRes::resolution2eval(int &count_res, double step,
								double &resolution, double &last_resolution,
								double &freq, double &freqL,
								int &last_fourier_idx,
								bool &continueIter,	bool &breakIter,
								bool &doNextIteration)
{
	if (automaticMode)
		resolution = minRes + count_res*step;
	else
		resolution = maxRes - count_res*step;
	freq = sampling/resolution;
	++count_res;

	double Nyquist = 2*sampling;
	double aux_frequency;
	int fourier_idx;

	DIGFREQ2FFT_IDX(freq, ZSIZE(VRiesz), fourier_idx);

//	std::cout << "Resolution = " << resolution << "   iter = " << count_res-1 << std::endl;
//	std::cout << "freq = " << freq << "   Fourier index = " << fourier_idx << std::endl;

	FFT_IDX2DIGFREQ(fourier_idx, ZSIZE(VRiesz), aux_frequency);

	freq = aux_frequency;

	if (fourier_idx == last_fourier_idx)
	{
//		std::cout << "entro en el if"  << std::endl;
		continueIter = true;
		return;
	}

	last_fourier_idx = fourier_idx;
	resolution = sampling/aux_frequency;


	if (count_res == 0)
		last_resolution = resolution;

	if ( ( resolution<Nyquist ))// || (resolution > last_resolution) )
	{
		//std::cout << "Nyquist limit reached" << std::endl;
		breakIter = true;
		return;
	}
//	if ( ( freq>0.495))// || (resolution > last_resolution) )
//	{
//		std::cout << "Nyquist limit reached" << std::endl;
//		doNextIteration = false;
//		return;
//	}

	freqL = sampling/(resolution + step);

	int fourier_idx_2;

	DIGFREQ2FFT_IDX(freqL, ZSIZE(VRiesz), fourier_idx_2);

	if (fourier_idx_2 == fourier_idx)
	{
		if (fourier_idx > 0){
			//std::cout << " index low =  " << (fourier_idx - 1) << std::endl;
			FFT_IDX2DIGFREQ(fourier_idx - 1, ZSIZE(VRiesz), freqL);
		}
		else{
			freqL = sampling/(resolution + step);
		}
	}

}


void ProgMonogenicSignalRes::run()
{
	produceSideInfo();

	Image<double> outputResolution;
	outputResolution().initZeros(VRiesz);

	MultidimArray<int> &pMask = mask();
	MultidimArray<double> &pOutputResolution = outputResolution();
	MultidimArray<double> &pVfiltered = Vfiltered();
	MultidimArray<double> &pVresolutionFiltered = VresolutionFiltered();
	MultidimArray<double> amplitudeMS, amplitudeMN;

	std::cout << "Looking for maximum frequency ..." << std::endl;
	double criticalZ=icdf_gauss(significance);
	double criticalW=-1;
	double resolution, resolution_2, last_resolution = 10000;  //A huge value for achieving
												//last_resolution < resolution
	double freq, freqH, freqL, resVal, counter;
	double max_meanS = -1e38;
	double cut_value = 0.025;


//	Image<int> imgMask2;
//	imgMask2 = mask;
//	imgMask2.write("mascara_original.vol");

	double R_ = freq_step;

	if (R_<0.25)
		R_=0.25;

	double w0, wF;
	double Nyquist = 2*sampling;
	if (minRes<2*sampling)
		minRes = Nyquist;
	if (automaticMode)
		minRes = Nyquist;
	else
	{
		w0 = sampling/maxRes;
		wF = sampling/minRes;
	}
	double w=w0;
	bool doNextIteration=true;
	bool lefttrimming = false;
	int fourier_idx, last_fourier_idx = -1, fourier_idx_2;

	//A first MonoRes estimation to get an accurate mask

	double mean_Signal, mean_noise, thresholdFirstEstimation;

	DIGFREQ2FFT_IDX((maxRes+3)/sampling, ZSIZE(VRiesz), fourier_idx);

	FFT_IDX2DIGFREQ(fourier_idx, ZSIZE(VRiesz), freq);
	FFT_IDX2DIGFREQ(fourier_idx + 2, ZSIZE(VRiesz), freqH);
	FFT_IDX2DIGFREQ(fourier_idx - 2, ZSIZE(VRiesz), freqL);

	//std::cout << " freq = " << freq << " freqH = " << freqH << " freqL= " << freq <<std::endl;

	int count_res = 0;
	FileName fnDebug;

	firstMonoResEstimation(fftV, freq, freqH, freqL, amplitudeMS,
			count_res, fnDebug, mean_Signal, mean_noise, thresholdFirstEstimation);

	//refining the mask
//	std::cout << "mean_Signal = " << mean_Signal << std::endl;
	NVoxelsOriginalMask = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
	{
		if (DIRECT_MULTIDIM_ELEM(pMask, n) >=1)
		{
			if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)<thresholdFirstEstimation)
			{
				DIRECT_MULTIDIM_ELEM(pMask, n) = 0;
			}
			else
			{
				++NVoxelsOriginalMask;
			}
		}
	}

//	imgMask2 = mask;
//	imgMask2.write("mascara_refinada.vol");

	int iter=0;
	std::vector<double> list;

	std::cout << "Analyzing frequencies" << std::endl;
	std::vector<double> noiseValues;

	do
	{

		bool continueIter = false;
		bool breakIter = false;

		resolution2eval(count_res, R_,
						resolution, last_resolution,
						freq, freqH,
						last_fourier_idx, continueIter, breakIter, doNextIteration);

		if (continueIter)
			continue;

		if (breakIter)
			break;

		std::cout << "resolution = " << resolution << std::endl;


		list.push_back(resolution);

		if (iter <2)
			resolution_2 = list[0];
		else
			resolution_2 = list[iter - 2];

		fnDebug = "Signal";

		freqL = freq + 0.01;
//		if (freqL>=0.5)
//			freqL = 0.5;

		amplitudeMonogenicSignal3D(fftV, freq, freqH, freqL, amplitudeMS, iter, fnDebug);
		if (halfMapsGiven)
		{
			fnDebug = "Noise";
			amplitudeMonogenicSignal3D(*fftN, freq, freqH, freqL, amplitudeMN, iter, fnDebug);
		}


		double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
		noiseValues.clear();

		if (exactres)
		{
			if (halfMapsGiven)
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
				{
					if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
					{
						double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
						double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
						sumS  += amplitudeValue;
						sumS2 += amplitudeValue*amplitudeValue;
						noiseValues.push_back(amplitudeValueN);
						sumN  += amplitudeValueN;
						sumN2 += amplitudeValueN*amplitudeValueN;
						++NS;
						++NN;
					}
				}
			}
			else
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
				{
					double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
					if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
					{
						sumS  += amplitudeValue;
						sumS2 += amplitudeValue*amplitudeValue;
						++NS;
					}
					else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
					{
						noiseValues.push_back(amplitudeValue);
						sumN  += amplitudeValue;
						sumN2 += amplitudeValue*amplitudeValue;
						++NN;
					}
				}
			}
		}
		else
		{
			if (halfMapsGiven)
			{
				if (noiseOnlyInHalves)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
					{
						double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
						double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
						{
							sumS  += amplitudeValue;
							sumS2 += amplitudeValue*amplitudeValue;
							++NS;
							sumN  += amplitudeValueN;
							sumN2 += amplitudeValueN*amplitudeValueN;
							++NN;
						}
					}
				}
				else
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
					{
						double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
						double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
						{
							sumS  += amplitudeValue;
							sumS2 += amplitudeValue*amplitudeValue;
							++NS;
						}
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=0)
						{
							sumN  += amplitudeValueN;
							sumN2 += amplitudeValueN*amplitudeValueN;
							++NN;
						}
					}
				}
			}
			else
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
				{
					double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
					if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
					{
						sumS  += amplitudeValue;
						sumS2 += amplitudeValue*amplitudeValue;
						++NS;
					}
					else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
					{
						sumN  += amplitudeValue;
						sumN2 += amplitudeValue*amplitudeValue;
						++NN;
					}
				}
			}
		}
	
		#ifdef DEBUG
		std::cout << "NS" << NS << std::endl;
		std::cout << "NVoxelsOriginalMask" << NVoxelsOriginalMask << std::endl;
		std::cout << "NS/NVoxelsOriginalMask = " << NS/NVoxelsOriginalMask << std::endl;
		#endif
		
//		std::cout << "NS = " << NS << "  NVoxelsOriginalMask" << NVoxelsOriginalMask << std::endl;
		if ( (NS/NVoxelsOriginalMask)<cut_value ) //when the 2.5% is reached then the iterative process stops
		{
			std::cout << "Search of resolutions stopped due to mask has been completed" << std::endl;
			doNextIteration =false;
			Nvoxels = 0;

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
			{
				if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n) == 0)
					DIRECT_MULTIDIM_ELEM(pMask, n) = 0;
				else
				{
					Nvoxels++;
					DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
				}
			}

			#ifdef DEBUG_MASK
			mask.write("partial_mask.vol");
			#endif

			lefttrimming = true;
		}
		else
		{
			if (NS == 0)
			{
				std::cout << "There are no points to compute inside the mask" << std::endl;
				std::cout << "If the number of computed frequencies is low, perhaps the provided"
						"mask is not enough tight to the volume, in that case please try another mask" << std::endl;
				break;
			}

			double meanS=sumS/NS;
			double sigma2S=sumS2/NS-meanS*meanS;
			double meanN=sumN/NN;
			double sigma2N=sumN2/NN-meanN*meanN;

			// Check local resolution
			double thresholdNoise;
			if (exactres)
			{
				std::sort(noiseValues.begin(),noiseValues.end());
				thresholdNoise = noiseValues[size_t(noiseValues.size()*significance)];
			}
			else
				thresholdNoise = meanN+criticalZ*sqrt(sigma2N);

			#ifdef DEBUG
			  std::cout << "Iteration = " << iter << ",   Resolution= " << resolution <<
					  ",   Signal = " << meanS << ",   Noise = " << meanN << ",  Threshold = "
					  << thresholdNoise <<std::endl;
			#endif

			double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);

//			std::cout << "z = " << z << "  zcritical = " << criticalZ << std::endl;

			if (automaticMode)
			{
				if (z>criticalZ)
				{
//					std::cout << "z<criticalZ---->ENTRO" <<std::endl;

					// Check local resolution
					double thresholdNoise;
					if (exactres)
					{
						std::sort(noiseValues.begin(),noiseValues.end());
						thresholdNoise = noiseValues[size_t(noiseValues.size()*significance)];
					}
					else
						thresholdNoise = meanN+criticalZ*sqrt(sigma2N);

					#ifdef DEBUG
					  std::cout << "Iteration = " << iter << ",   Resolution= " << resolution <<
							  ",   Signal = " << meanS << ",   Noise = " << meanN << ",  Threshold = "
							  << thresholdNoise <<std::endl;
					#endif
					double NRES=0;
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
					{
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
							if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)<thresholdNoise)
							{
								DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
								DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution;//sampling/freq;
								if (fnSpatial!="")
									DIRECT_MULTIDIM_ELEM(pVresolutionFiltered,n)=DIRECT_MULTIDIM_ELEM(pVfiltered,n);
							}
							else{
								++NRES;
								DIRECT_MULTIDIM_ELEM(pMask, n) += 1;
								if (DIRECT_MULTIDIM_ELEM(pMask, n) >2)
								{
									DIRECT_MULTIDIM_ELEM(pMask, n) = -1;
									DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution_2;//maxRes - counter*R_;
								}
							}
					}
//					std::cout << " NRES" << NRES << std::endl;
					if ( ( NRES/((double)NVoxelsOriginalMask) ) > 0.8 )
					{
						std::cout << " entroooo" << std::endl;
						mask.read(fnMask);
					}
				}
			}
			else
			{
				if (meanS>max_meanS)
					max_meanS = meanS;

				if (meanS<0.001*max_meanS)
				{
					std::cout << "Search of resolutions stopped due to too low signal" << std::endl;
					break;
				}

				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
				{
					if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
						if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
						{
							DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
							DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution;//sampling/freq;
							if (fnSpatial!="")
								DIRECT_MULTIDIM_ELEM(pVresolutionFiltered,n)=DIRECT_MULTIDIM_ELEM(pVfiltered,n);
						}
						else{
							DIRECT_MULTIDIM_ELEM(pMask, n) += 1;
							if (DIRECT_MULTIDIM_ELEM(pMask, n) >2)
							{
								DIRECT_MULTIDIM_ELEM(pMask, n) = -1;
								DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution_2;//maxRes - counter*R_;
							}
						}
				}
				#ifdef DEBUG_MASK
				FileName fnmask_debug;
				fnmask_debug = formatString("maske_%i.vol", iter);
				mask.write(fnmask_debug);
				#endif

				// Is the mean inside the signal significantly different from the noise?
				z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);
				#ifdef DEBUG
					std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
					std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
					std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
					std::cout << "  z=" << z << " (" << criticalZ << ")" << std::endl;
				#endif
				if (z<criticalZ)
				{
					criticalW = freq;
					std::cout << "Search stopped due to z>Z (hypothesis test)" << std::endl;
					doNextIteration=false;
				}
				if (doNextIteration)
				{
					if (resolution <= (minRes-0.001))
						doNextIteration = false;
				}
			}
		}
		iter++;
		last_resolution = resolution;
	} while (doNextIteration);

	if (lefttrimming == false)
	{
	  Nvoxels = 0;
	  FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
	  {
	    if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n) == 0)
	      DIRECT_MULTIDIM_ELEM(pMask, n) = 0;
	    else
	    {
	      Nvoxels++;
	      DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
//	      if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
//	      {
//	    	  //DIRECT_MULTIDIM_ELEM(pMask, n) = 0;
//	    	  DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
//	      }
//	      else
//	      {
////	    	  DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
//	    	  DIRECT_MULTIDIM_ELEM(pMask, n) = 0;
//	      }
	    }
	  }
	#ifdef DEBUG_MASK
	  //mask.write(fnMaskOut);
	#endif
	}
	amplitudeMN.clear();
	amplitudeMS.clear();

	double last_resolution_2 = resolution;
	if (fnSym!="c1")
	{
		SymList SL;
		SL.readSymmetryFile(fnSym);
		MultidimArray<double> VSimetrized;
		symmetrizeVolume(SL, pOutputResolution, VSimetrized, LINEAR, DONT_WRAP);
		outputResolution() = VSimetrized;
		VSimetrized.clear();
	}

	#ifdef DEBUG
		outputResolution.write("resolution_simple_simmetrized.vol");
	#endif

	MultidimArray<double> resolutionFiltered, resolutionChimera;
	postProcessingLocalResolutions(pOutputResolution, list, resolutionChimera, cut_value, pMask);


	Image<double> outputResolutionImage;
	outputResolutionImage() = pOutputResolution;//resolutionFiltered;
	outputResolutionImage.write(fnOut);
	outputResolutionImage() = resolutionChimera;
	outputResolutionImage.write(fnchim);


	#ifdef DEBUG
		outputResolution.write("resolution_simple.vol");
	#endif

	if (fnSpatial!="")
	{
		mask.read(fnMask);
		mask().setXmippOrigin();
		Vfiltered.read(fnVol);
		pVfiltered=Vfiltered();
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pVfiltered)
		if (DIRECT_MULTIDIM_ELEM(pMask,n)==1)
			DIRECT_MULTIDIM_ELEM(pVfiltered,n)-=DIRECT_MULTIDIM_ELEM(pVresolutionFiltered,n);
		Vfiltered.write(fnSpatial);

		VresolutionFiltered().clear();
		Vfiltered().clear();
	}

	MetaData md;
	size_t objId;
	objId = md.addObject();
	md.setValue(MDL_IMAGE, fnOut, objId);
	md.setValue(MDL_COUNT, (size_t) NVoxelsOriginalMask, objId);
	md.setValue(MDL_COUNT2, (size_t) Nvoxels, objId);

	md.write(fnMd);
}
