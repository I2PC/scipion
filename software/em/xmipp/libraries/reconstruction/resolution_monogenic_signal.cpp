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
#define DEBUG
#define DEBUG_MASK

void ProgMonogenicSignalRes::readParams()
{
	fnVol = getParam("--vol");
	//fnVol1 = getParam("--vol1");
	fnVol2 = getParam("--vol2");
	fnMeanVol = getParam("--meanVol");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	fnchim = getParam("--chimera_volume");
	sampling = getDoubleParam("--sampling_rate");
	R = getDoubleParam("--volumeRadius");
	minRes = getDoubleParam("--minRes");
	maxRes = getDoubleParam("--maxRes");
	fnSym = getParam("--sym");
	N_freq = getDoubleParam("--number_frequencies");
	trimBound = getDoubleParam("--trimmed");
	linearchk = checkParam("--linear");
	exactres = checkParam("--exact");
	fnSpatial = getParam("--filtered_volume");
	significance = getDoubleParam("--significance");
}


void ProgMonogenicSignalRes::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file=\"\">   : Input volume");
	addParamsLine("  [--mask <vol_file=\"\">]  : Mask defining the macromolecule");
	addParamsLine("                          :+ If two half volume are given, the noise is estimated from them");
	addParamsLine("                          :+ Otherwise the noise is estimated outside the mask");
	addParamsLine("  [--vol2 <vol_file=\"\">]: Half volume 2");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--meanVol <vol_file=\"\">]: Mean volume of half1 and half2 (only it is neccesary the two haves are used)");
	addParamsLine("  --sym <symmetry>: Symmetry (c1, c2, c3,..d1, d2, d3,...)");
	addParamsLine("  [--chimera_volume <output=\"Chimera_resolution_volume.vol\">]: Local resolution volume for chimera viewer (in Angstroms)");
	addParamsLine("  [--sampling_rate <s=1>]   : Sampling rate (A/px)");
	addParamsLine("                            : Use -1 to disable this option");
	addParamsLine("  [--volumeRadius <s=100>]   : This parameter determines the radius of a sphere where the volume is");
	addParamsLine("  [--number_frequencies <w=50>]       : The resolution is computed at a number of frequencies between mininum and");
	addParamsLine("                            : maximum resolution px/A. This parameter determines that number");
	addParamsLine("  [--minRes <s=30>]         : Minimum resolution (A)");
	addParamsLine("  [--maxRes <s=1>]          : Maximum resolution (A)");
	addParamsLine("  [--trimmed <s=0.5>]         : Trimming percentile");
	addParamsLine("  [--linear]                : The search for resolution is linear (equidistance between resolutions).");
	addParamsLine("  [--exact]                 : The search for resolution will be exact (slower) of approximated (fast).");
	addParamsLine("                            : Usually there are no difference between both in the resolution map.");
	addParamsLine("  [--filtered_volume <vol_file=\"\">]       : The input volume is locally filtered at local resolutions.");
	addParamsLine("  [--significance <s=0.95>]  : The level of confidence for the hypothesis test.");
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

	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
	{
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
}


void ProgMonogenicSignalRes::amplitudeMonogenicSignal3D(MultidimArray< std::complex<double> > &myfftV,
		double w1, MultidimArray<double> &amplitude, int count, FileName fnDebug)
{
	fftVRiesz.initZeros(myfftV);
	amplitude.resizeNoCopy(VRiesz);


	// Filter the input volume and add it to amplitude
	long n=0;
	double iu1=1.0/w1;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
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


	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	double uz, uy, ux;
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
				{
					FFT_IDX2DIGFREQ(j,XSIZE(amplitude),ux);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (ux*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
				}
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(amplitude),uy);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (uy*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate third component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(amplitude),uz);
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (uz*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
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

	// Low pass filter the monogenic amplitude
	lowPassFilter.w1 = w1;
	amplitude.setXmippOrigin();
	lowPassFilter.applyMaskSpace(amplitude);

	#ifdef DEBUG
	//Image<double> saveImg2;
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
	double resolution, last_resolution = 10000;  //A huge value for achieving last_resolution < resolution
	double freq;
	double max_meanS = -1e38;

	double range = maxRes-minRes;
	double R_ = range/N_freq;

	double w0 = sampling/maxRes;
	double wF = sampling/minRes;
	double w=w0;
	double stepW = (wF-w0)/N_freq;
	bool doNextIteration=true;
	int iter=0;
	int count_res = 0;

	std::cout << "Analyzing frequencies" << std::endl;
	std::vector<double> noiseValues;
	FileName fnDebug;
	do
	{
		if (linearchk ==true)
		{
			resolution = maxRes - count_res*R_;
			freq = sampling/resolution;
			++count_res;
		}
		else
		{
			resolution =  sampling/w;
			freq = sampling/resolution;
		}

		std::cout << "Iteration " << iter << " Freq = " << freq << " Resolution = " << resolution << " (A)" << std::endl;

		fnDebug = "Signal";
		amplitudeMonogenicSignal3D(fftV, freq, amplitudeMS, iter, fnDebug);
		if (halfMapsGiven)
		{
			fnDebug = "Noise";
			amplitudeMonogenicSignal3D(*fftN, freq, amplitudeMN, iter, fnDebug);
		}

		double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
		noiseValues.clear();
		if (exactres)
		{
			if (halfMapsGiven)
					{
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
						{
							double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
							double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
							if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
							{
								sumS  += amplitudeValue;
								sumS2 += amplitudeValue*amplitudeValue;
								noiseValues.push_back(amplitudeValueN);
								sumN  += amplitudeValueN;
								sumN2 += amplitudeValueN*amplitudeValueN;
								++NS;
								++NN;
							}
							else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
							{
								noiseValues.push_back(amplitudeValueN);
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
							if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
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
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
						{
							double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
							double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
							if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
							{
								sumS  += amplitudeValue;
								sumS2 += amplitudeValue*amplitudeValue;
								sumN  += amplitudeValueN;
								sumN2 += amplitudeValueN*amplitudeValueN;
								++NS;
								++NN;
							}
							else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
							{
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
							if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
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

		if (NS == 0)
		{
			std::cout << "There are no points to compute inside the mask" << std::endl;
			break;
		}

		double meanS=sumS/NS;
		double sigma2S=sumS2/NS-meanS*meanS;
		double meanN=sumN/NN;
		double sigma2N=sumN2/NN-meanN*meanN;

		if (meanS>max_meanS)
			max_meanS = meanS;

		if (meanS<0.001*max_meanS)
		{
			std::cout << "Search of resolutions stopped due to too low signal" << std::endl;
			break;
		}

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
		  std::cout << "Iteration = " << iter << ",   Resolution= " << resolution << ",   Signal = " << meanS << ",   Noise = " << meanN << ",  Threshold = " << thresholdNoise <<std::endl;
		#endif

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
				if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
				{
					DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = sampling/freq;
					if (fnSpatial!="")
						DIRECT_MULTIDIM_ELEM(pVresolutionFiltered,n)=DIRECT_MULTIDIM_ELEM(pVfiltered,n);
				}
				else
					DIRECT_MULTIDIM_ELEM(pMask, n) = 0;
		}

		// Is the mean inside the signal significantly different from the noise?
		double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);
		if (verbose >=2)
		{
			std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
			std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
			std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
			std::cout << "  z=" << z << " (" << criticalZ << ")" << std::endl;
		}
		if (z<criticalZ)
		{
			criticalW = freq;
			doNextIteration=false;
		}

		if (doNextIteration)
		{
			if (linearchk == false)
			{
				last_resolution = resolution;
				w+=stepW;
				resolution = sampling/w;


			if (last_resolution-resolution<0.1)
			{
				resolution=last_resolution-0.1;
				w = sampling/resolution;
			}

			if (w > wF)
			{
				doNextIteration = false;
				std::cout << "Search of resolutions stopped due to out of resolution range" << std::endl;
			}
			}
			else
			{
				if (resolution == minRes)
					doNextIteration = false;
			}
		}
		iter++;
	} while (doNextIteration);


	double last_resolution_2 = resolution;
	std::cout << "last computed resolution = " << last_resolution_2 << std::endl;

	#ifdef DEBUG
		outputResolution.write("resolution_simple.vol");
	#endif

	if (fnSym!="c1")
	{
		SymList SL;
		SL.readSymmetryFile(fnSym);
		MultidimArray<double> VSimetrized;
		symmetrizeVolume(SL, pOutputResolution, VSimetrized, LINEAR, DONT_WRAP);
		outputResolution() = VSimetrized;
	}

	#ifdef DEBUG
		outputResolution.write("resolution_simple_simmetrized.vol");
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


		double resValue;
		if (trimBound>0)
		{
			// Count number of voxels with resolution
			size_t N=0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
				if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n)>(last_resolution_2-0.001)) //the value 0.001 is a tolerance
					++N;

			// Get all resolution values
			MultidimArray<double> resolutions(N);
			size_t N_iter=0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
				if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n)>(last_resolution_2-0.001))
					DIRECT_MULTIDIM_ELEM(resolutions,N_iter++)=DIRECT_MULTIDIM_ELEM(pOutputResolution, n);

			// Sort value and get threshold
			std::sort(&A1D_ELEM(resolutions,0),&A1D_ELEM(resolutions,N));
			double threshold=A1D_ELEM(resolutions,(int)(trimBound*N/100));
			double filling_value = A1D_ELEM(resolutions, (int)(0.5*N));
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
			{
				if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n)>threshold){
					DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = filling_value;
				}
//				else{
//					DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = 0;
//				}
			}

			outputResolution.write(fnOut);
		}
		else
		{
			outputResolution.write(fnOut);
		}


		// Write volume for Chimera
		if (fnchim != "")
		{
			//Calculating mean
			double SumRes = 0, Nsum = 0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
			{
				resValue = DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
				if (resValue>(last_resolution_2-0.001))
				{
					SumRes  += resValue;
					Nsum += 1;
				}
			}

			double meanRes = SumRes/Nsum;

			std::cout << "mean = " << meanRes << std::endl;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
			{
				resValue = DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
				if (resValue<=(last_resolution_2-0.001))
					DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = meanRes;
			}

			outputResolution.write(fnchim);
		}
}
