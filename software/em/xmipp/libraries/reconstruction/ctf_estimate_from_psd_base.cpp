/***************************************************************************
 *
 * Authors:     Javier Mota Garcia (jmota@cnb.csic.es)
 *              Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
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

#include "ctf_estimate_from_psd_base.h"
#include "ctf_enhance_psd.h"

#include <data/args.h>
#include <data/histogram.h>
#include <data/filters.h>
#include <data/xmipp_fft.h>

const double ProgCTFBasicParams::penalty = 32.0; // Maximum penalization

/* Read parameters --------------------------------------------------------- */
void ProgCTFBasicParams::readBasicParams(XmippProgram *program)
{
	downsampleFactor = program->getDoubleParam("--downSamplingPerformed");
    show_optimization = program->checkParam("--show_optimization");
    min_freq = program->getDoubleParam("--min_freq");
    max_freq = program->getDoubleParam("--max_freq");
    defocus_range = program->getDoubleParam("--defocus_range");
    modelSimplification = program->getIntParam("--model_simplification");
    bootstrap = program->checkParam("--bootstrapFit");
    fastDefocusEstimate = program->checkParam("--fastDefocus");
    if (fastDefocusEstimate)
    {
        lambdaPhase=program->getDoubleParam("--fastDefocus",0);
        sizeWindowPhase=program->getIntParam("--fastDefocus",1);
    }
    ctfmodelSize = program->getIntParam("--ctfmodelSize");
    enhanced_weight = program->getDoubleParam("--enhance_weight");
    if (!program->checkParam("--enhance_min_freq"))
        f1 = (max_freq > 0.35) ? 0.01 : 0.02;
    else
        f1 = program->getDoubleParam("--enhance_min_freq");
    if (!program->checkParam("--enhance_max_freq"))
        f2 = (max_freq > 0.35) ? 0.08 : 0.15;
    else
        f2 = program->getDoubleParam("--enhance_max_freq");

}

void ProgCTFBasicParams::readParams()
{
	fn_psd = getParam("--psd");
    readBasicParams(this);
}

/* Show -------------------------------------------------------------------- */
void ProgCTFBasicParams::show()
{
    if (verbose==0)
        return;
    std::cout
    << "PSD file:            " << fn_psd << std::endl
    << "Downsampling:        " << downsampleFactor << std::endl
    << "Min Freq.:           " << min_freq << std::endl
    << "Max Freq.:           " << max_freq << std::endl
    << "Sampling:            " << Tm << std::endl
    << "Defocus range:       " << defocus_range << std::endl
    << "ctfmodelSize:        " << ctfmodelSize << std::endl
    << "Enhance min freq:    " << f1 << std::endl
    << "Enhance max freq:    " << f2 << std::endl
    << "Enhance weight:      " << enhanced_weight << std::endl
    << "Model simplification:" << modelSimplification << std::endl
    << "Bootstrap:           " << bootstrap << std::endl
    << "Fast defocus:        " << fastDefocusEstimate << std::endl
    ;
    if (fastDefocusEstimate)
        std::cout
        << "Regularization factor:      " << lambdaPhase << std::endl
        << "Window size:                " << sizeWindowPhase << std::endl;

}

/* Usage ------------------------------------------------------------------- */
void ProgCTFBasicParams::defineBasicParams(XmippProgram * program)
{
    program->addSeeAlsoLine("ctf_enhance_psd");
    program->addParamsLine("== Downsampling");
    program->addParamsLine(
        "  [--downSamplingPerformed <F=1>] : Downsampling performed to produce this PSD");
    program->addParamsLine(
        "                                  : Note that the output CTF model is referred to the original sampling rate");
    program->addParamsLine(
        "                                  : not the one of the downsampled image.");
    program->addParamsLine("== CTF fit: Optimization constraints");
    program->addParamsLine(
        "   [--min_freq <fmin=0.03>]     : Minimum digital frequency (<0.5) to use in adjust. Its value");
    program->addParamsLine(
        "                                : should be a little lower than the dig. freq. of the first ");
    program->addParamsLine("                                : CTF zero.");
    program->addParamsLine(
        "   [--max_freq <fmax=0.35>]     : Maximum digital frequency (<0.5) to use in adjust.");
    program->addParamsLine(
        "                                : It should be higher than the last zero of the CTF.");
    program->addParamsLine(
        "   [--fastDefocus <lambda=2> <size=10>] : Estimate first defocus with Zernike polynomials");
    program->addParamsLine(
        "                                :+Lambda is a regularization factor used during the estimation of the CTF phase");
    program->addParamsLine(
        "                                :+During the estimation, the phase values are averaged within a window of this size");
    program->addParamsLine(
        "   [--defocus_range <D=8000>]   : Defocus range in Angstroms");
    program->addParamsLine(
        "   [--show_optimization+]       : Show optimization process");
    program->addParamsLine(
        "   [--radial_noise++]           : By default, noise is astigmatic");
    program->addParamsLine(
        "   [--enhance_weight++ <w=1>]   : Weight of the enhanced term");
    program->addParamsLine(
        "   [--model_simplification++ <s=0>]: 0 (no simplification)");
    program->addParamsLine(
        "                                : 1 (simplified envelope)");
    program->addParamsLine(
        "                                : 2 (last Gaussian removal)");
    program->addParamsLine(
        "                                : 3 (symmetric intermediate Gaussian)");
    program->addParamsLine(
        "   [--bootstrapFit++ <N=-1>]    : Perform bootstrap fit (Fourier pixels are randomly chosen)");
    program->addParamsLine(
        "                                : This is used to test the variability of the fit");
    program->addParamsLine(
        "                                : N defines the number of times the fit is repeated");
    program->addParamsLine("==+ CTF fit: Output CTF models");
    program->addParamsLine(
        "   [--ctfmodelSize <size=256>]  : Size for the ctfmodel thumbnails");
    program->addParamsLine("==+ PSD enhancement");
    program->addParamsLine(
        "   [--enhance_min_freq <f1>]    : Bandpass cutoff. Normalized to 0.5");
    program->addParamsLine(
        "                                : If fmax>0.35, f1 default=0.01");
    program->addParamsLine(
        "                                : If fmax<0.35, f1 default=0.02");
    program->addParamsLine(
        "   [--enhance_max_freq <f2>]    : Bandpass cutoff. Normalized to 0.5.");
    program->addParamsLine(
        "                                : If fmax>0.35, f2 default=0.08");
    program->addParamsLine(
        "                                : If fmax<0.35, f2 default=0.15");
}

void ProgCTFBasicParams::defineParams()
{
    addUsageLine("Adjust a parametric model to a PSD file.");
    addUsageLine(
        "The PSD is enhanced ([[http://www.ncbi.nlm.nih.gov/pubmed/16987671][See article]]). ");
    addUsageLine(
        "And finally, the CTF is fitted to the PSD, being guided by the enhanced PSD ");
    addUsageLine(
        "([[http://www.ncbi.nlm.nih.gov/pubmed/17911028][See article]]).");
    addParamsLine("   --psd <PSDfile> : PSD file");
    addSeeAlsoLine("ctf_estimate_from_micrograph");
    defineBasicParams(this);
}

/* Produce side information ------------------------------------------------ */
void ProgCTFBasicParams::produceSideInfo()
{
    // Resize the frequency
	x_digfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
	y_digfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    w_digfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    w_digfreq_r.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    x_contfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    y_contfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    w_contfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);

    Matrix1D<int> idx(2); // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane

    FOR_ALL_ELEMENTS_IN_ARRAY2D(x_digfreq)
    {
        XX(idx) = j;
        YY(idx) = i;

        // Digital frequency
        FFT_idx2digfreq(*f, idx, freq);
        x_digfreq(i, j) = XX(freq);
        y_digfreq(i, j) = YY(freq);
        w_digfreq(i, j) = freq.module();
        w_digfreq_r(i, j) = (int)(w_digfreq(i,j) * (double)YSIZE(w_digfreq));

        // Continuous frequency
        digfreq2contfreq(freq, freq, Tm);
        x_contfreq(i, j) = XX(freq);
        y_contfreq(i, j) = YY(freq);
        w_contfreq(i, j) = freq.module();
    }

    // Build frequency mask
    mask.initZeros(w_digfreq);
    w_count.initZeros(XSIZE(w_digfreq));
    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    {
        if (w_digfreq(i, j) >= max_freq
            || w_digfreq(i, j) <= min_freq)
            continue;
        mask(i, j) = 1;
        w_count(w_digfreq_r(i, j))++;
    }

    // Enhance PSD for ctfmodels
	ProgCTFEnhancePSD prm;
	prm.filter_w1 = 0.02;
	prm.filter_w2 = 0.2;
	prm.decay_width = 0.02;
	prm.mask_w1 = 0.01;
	prm.mask_w2 = 0.5;
	enhanced_ctftomodel() = ctftomodel();
	prm.applyFilter(enhanced_ctftomodel());
	if (fn_psd.find('@')==std::string::npos)
		enhanced_ctftomodel.write(fn_psd.withoutExtension() + "_enhanced_psd.xmp");
	else
		enhanced_ctftomodel.write(fn_psd.withoutExtension() + "_enhanced_psd.stk");
	CenterFFT(enhanced_ctftomodel(), false);
	enhanced_ctftomodel_fullsize() = enhanced_ctftomodel();

	// Enhance PSD for optimization
	prm.filter_w1 = f1;
	prm.filter_w2 = f2;
	enhanced_ctftomodel() = ctftomodel();
	prm.applyFilter(enhanced_ctftomodel());
	CenterFFT(enhanced_ctftomodel(), false);
	enhanced_ctftomodel().resize(w_digfreq);

	// Divide by the number of count at each frequency
	// and mask between min_freq and max_freq
	double min_val = enhanced_ctftomodel().computeMin();
	FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
	if (mask(i, j) <= 0)
		enhanced_ctftomodel(i, j) = min_val;

	MultidimArray<double> aux;
	medianFilter3x3(enhanced_ctftomodel(), aux);

	enhanced_ctftomodel() = aux;
	enhanced_ctftomodel().rangeAdjust(0, 1);

	FourierFilter Filter;
	Filter.FilterShape = RAISED_COSINE;
	Filter.FilterBand = HIGHPASS;
	Filter.w1 = 0.01;
	Filter.raised_w = 0.005;
	enhanced_ctftomodel().setXmippOrigin();
	Filter.generateMask(enhanced_ctftomodel());
	Filter.applyMaskSpace(enhanced_ctftomodel());
	STARTINGX(enhanced_ctftomodel()) = STARTINGY(enhanced_ctftomodel()) = 0;
	//global_prm->enhanced_ctftomodel()= enhanced_ctftomodel();

	// Compute now radial average of the enhanced_ctftomodel
	psd_exp_enhanced_radial_derivative.initZeros(XSIZE(enhanced_ctftomodel()));
	psd_theo_radial_derivative.initZeros(psd_exp_enhanced_radial_derivative);
	psd_exp_enhanced_radial.initZeros(psd_exp_enhanced_radial_derivative);
	psd_exp_radial.initZeros(psd_exp_enhanced_radial_derivative);
	psd_theo_radial.initZeros(psd_exp_enhanced_radial_derivative);
	w_digfreq_r_iN.initZeros(psd_exp_enhanced_radial_derivative);
	FOR_ALL_ELEMENTS_IN_ARRAY2D(enhanced_ctftomodel())
	{
		if (w_digfreq(i,j)>min_freq && w_digfreq(i,j)<max_freq)
		{
			int r = w_digfreq_r(i, j);
			w_digfreq_r_iN(r)+=1;
			psd_exp_enhanced_radial(r) += enhanced_ctftomodel(i, j);
			psd_exp_radial(r) += ctftomodel(i, j);
		}
	}
	FOR_ALL_ELEMENTS_IN_ARRAY1D(w_digfreq_r_iN)
	if (w_digfreq_r_iN(i)>0)
	{
		w_digfreq_r_iN(i)=1.0/w_digfreq_r_iN(i);
		psd_exp_enhanced_radial(i)*=w_digfreq_r_iN(i);
		psd_exp_radial(i)*=w_digfreq_r_iN(i);
	}

	// Compute its derivative
	int state=0;
	double maxDiff=0;
	FOR_ALL_ELEMENTS_IN_ARRAY1D(psd_exp_enhanced_radial)
	{
		switch (state)
		{
		case 0:
			if (w_digfreq(i,0)>min_freq)
				state=1;
			break;
		case 1:
			state=2; // Skip first sample
			break;
		case 2:
			if (w_digfreq(i,0)>max_freq)
				state=3;
			else
			{
				double diff=psd_exp_enhanced_radial(i)-psd_exp_enhanced_radial(i-1);
				psd_exp_enhanced_radial_derivative(i)=diff;
				maxDiff=std::max(maxDiff,fabs(diff));
			}
			break;
		}
	}
	psd_exp_enhanced_radial_derivative/=maxDiff;
	if (show_optimization)
	{
		psd_exp_enhanced_radial.write("PPPexpEnhanced_fast.txt");
		psd_exp_radial.write("PPPexp_fast.txt");
	}

}

double evaluateIceness(MultidimArray<double> &psd, double Tm)
{
	double R4_4=floor(XSIZE(psd)*Tm/4.4);
	double R4_0=floor(XSIZE(psd)*Tm/4.0);
	double R3_6=floor(XSIZE(psd)*Tm/3.6);
	R4_4*=R4_4;
	R4_0*=R4_0;
	R3_6*=R3_6;
	double N1=0, N2=0, avg1=0, avg2=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(psd)
	{
		double R2=i*i+j*j;
		if (R2>R4_4)
		{
			if (R2<R4_0)
			{
				avg1+=A2D_ELEM(psd,i,j);
				N1+=1;
			}
			else if (R2<R3_6)
			{
				avg2+=A2D_ELEM(psd,i,j);
				N2+=1;
			}
		}
	}
	if (N1>0 && N2>0)
	{
		avg1/=N1;
		avg2/=N2;
		return avg2/avg1;
	}
	return -1;
}
