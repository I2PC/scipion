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

#include "ctf_estimate_from_psd_fast.h"
#include "ctf_enhance_psd.h"

#include <data/args.h>
#include <data/histogram.h>
#include <data/filters.h>
#include <data/xmipp_fft.h>

/* prototypes */
double CTF_fitness_fast(double *, void *);

/* Number of CTF parameters */
#define ALL_CTF_PARAMETERS           26
#define CTF_PARAMETERS               20
#define PARAMETRIC_CTF_PARAMETERS    16
#define BACKGROUND_CTF_PARAMETERS    9
#define SQRT_CTF_PARAMETERS           6
#define ENVELOPE_PARAMETERS          11
#define DEFOCUS_PARAMETERS            3
#define FIRST_SQRT_PARAMETER         14
#define FIRST_ENVELOPE_PARAMETER      2
#define FIRST_DEFOCUS_PARAMETER       0

//#define DEBUG_WITH_TEXTFILES
#ifdef DEBUG_WITH_TEXTFILES
std::ofstream fhDebug;
#define DEBUG_OPEN_TEXTFILE(fnRoot) fhDebug.open((fnRoot+"_debug.txt").c_str());
#define DEBUG_TEXTFILE(str) fhDebug << time (NULL) << " " << str << std::endl;
#define DEBUG_MODEL_TEXTFILE fhDebug << global_ctfmodel << std::endl;
#else
#define DEBUG_OPEN_TEXTFILE(fnRoot);
#define DEBUG_TEXTFILE(str);
#define DEBUG_MODEL_TEXTFILE
#endif

namespace AdjustCTF1D
{
// Some aliases
ProgCTFEstimateFromPSDFast *global_prm;
MultidimArray<double> *f; // The CTF to model
Matrix1D<double> *global_adjust; // Current theoretical adjustment
double global_corr13;

// Frequency of each point in digital units
MultidimArray<double> global_x_digfreq;
MultidimArray<double> global_y_digfreq;
MultidimArray<double> global_x_contfreq;
MultidimArray<double> global_y_contfreq;
MultidimArray<double> global_w_digfreq;
MultidimArray<double> global_w_digfreq_r;
MultidimArray<double> global_w_contfreq;
MultidimArray<double> global_mask_between_zeroes;
MultidimArray<double> global_mask;
MultidimArray<double> global_w_count;

MultidimArray<double> global_psd_exp_radial;
MultidimArray<double> global_psd_exp_enhanced_radial;
MultidimArray<double> global_psd_exp_enhanced_radial_derivative;
MultidimArray<double> global_psd_theo_radial_derivative;
MultidimArray<double> global_psd_theo_radial;
MultidimArray<double> global_w_digfreq_r_iN;

// Penalization factor for the background
bool global_penalize;

// Speed up factor
int global_evaluation_reduction;

// Penalization for forbidden values of the parameters
double global_heavy_penalization;
double global_current_penalty;
const double global_penalty = 32; // Maximum penalization

// Program status
int global_action;// 0: Computing the background (sqrt)
// 1: Computing the full background
// 2: Computing the envelope
// 3: Computing defoci
// 4: Computing all CTF parameters
// 5: Computing all CTF parameters + Gaussian2
// 6: Computing all CTF parameters + Gaussian2
// 7: Produce output

// CTF model and noise model
CTF1D global_ctfmodel;
CTF1D global_ctfmodel_defoci;

// Maximum of the gaussian
double global_max_gauss_freq;

int global_show; // 0: Do not show
int global_contador = 0;

}

//Global values
// Autofocus
double global_value_th;
double global_min_freq;
double global_max_freq;

using namespace AdjustCTF1D;

#define ASSIGN_CTF_PARAM(index, paramName) if (ia <= index && l > 0) { ctf1Dmodel.paramName = p[index]; --l; }

/* Assign ctf1Dmodel from a vector and viceversa ----------------------------- */
void assignCTFfromParameters(double *p, CTF1D &ctf1Dmodel, int ia,
                             int l, int modelSimplification)
{
    ctf1Dmodel.Tm = global_prm->Tm;

    ASSIGN_CTF_PARAM(0, Defocus);
    ASSIGN_CTF_PARAM(1, kV);
    ASSIGN_CTF_PARAM(2, K);
    ASSIGN_CTF_PARAM(3, Cs);
    ASSIGN_CTF_PARAM(4, Ca);
    ASSIGN_CTF_PARAM(5, espr); //deltaV/V
    ASSIGN_CTF_PARAM(6, ispr); //deltaI/I
    ASSIGN_CTF_PARAM(7, alpha);
    ASSIGN_CTF_PARAM(8, DeltaF);
    ASSIGN_CTF_PARAM(9, DeltaR);
    ASSIGN_CTF_PARAM(10, envR0);
    ASSIGN_CTF_PARAM(11, envR1);
    ASSIGN_CTF_PARAM(12, envR2);
    ASSIGN_CTF_PARAM(13, Q0);
    ASSIGN_CTF_PARAM(14, base_line);
    ASSIGN_CTF_PARAM(15, sqrt_K); //Ks
    ASSIGN_CTF_PARAM(16, sq);
    ASSIGN_CTF_PARAM(17, bgR1);
    ASSIGN_CTF_PARAM(18, bgR2);
    ASSIGN_CTF_PARAM(19, bgR3);
    ASSIGN_CTF_PARAM(20, gaussian_K);
    ASSIGN_CTF_PARAM(21, sigma1);
    ASSIGN_CTF_PARAM(22, Gc1);
    ASSIGN_CTF_PARAM(23, gaussian_K2);
    ASSIGN_CTF_PARAM(24, sigma2);
    ASSIGN_CTF_PARAM(25, Gc2);
}//function assignCTFfromParameters

#define COPY_ctfmodel_TO_CURRENT_GUESS \
    assignParametersFromCTF(global_ctfmodel, \
                               MATRIX1D_ARRAY(*global_adjust),0,ALL_CTF_PARAMETERS, \
                               global_prm->modelSimplification);


#define ASSIGN_PARAM_CTF(index, paramName) if (ia <= index && l > 0) { p[index] = ctfmodel.paramName; --l; }

void assignParametersFromCTF(CTF1D &ctfmodel, double *p, int ia,
                             int l, int modelSimplification)
{
	 	ASSIGN_PARAM_CTF(0, Defocus);
	    ASSIGN_PARAM_CTF(1, kV);
	    ASSIGN_PARAM_CTF(2, K);
	    ASSIGN_PARAM_CTF(3, Cs);
	    ASSIGN_PARAM_CTF(4, Ca);
	    ASSIGN_PARAM_CTF(5, espr); //deltaV/V
	    ASSIGN_PARAM_CTF(6, ispr); //deltaI/I
	    ASSIGN_PARAM_CTF(7, alpha);
	    ASSIGN_PARAM_CTF(8, DeltaF);
	    ASSIGN_PARAM_CTF(9, DeltaR);
	    ASSIGN_PARAM_CTF(10, envR0);
	    ASSIGN_PARAM_CTF(11, envR1);
	    ASSIGN_PARAM_CTF(12, envR2);
	    ASSIGN_PARAM_CTF(13, Q0);
	    ASSIGN_PARAM_CTF(14, base_line);
	    ASSIGN_PARAM_CTF(15, sqrt_K); //Ks
	    ASSIGN_PARAM_CTF(16, sq);
	    ASSIGN_PARAM_CTF(17, bgR1);
	    ASSIGN_PARAM_CTF(18, bgR2);
	    ASSIGN_PARAM_CTF(19, bgR3);
	    ASSIGN_PARAM_CTF(20, gaussian_K);
	    ASSIGN_PARAM_CTF(21, sigma1);
	    ASSIGN_PARAM_CTF(22, Gc1);
	    ASSIGN_PARAM_CTF(23, gaussian_K2);
	    ASSIGN_PARAM_CTF(24, sigma2);
	    ASSIGN_PARAM_CTF(25, Gc2);
}

/* Read parameters --------------------------------------------------------- */
void ProgCTFEstimateFromPSDFast::readBasicParams(XmippProgram *program)
{
	downsampleFactor = program->getDoubleParam("--downSamplingPerformedfast");
	show_optimization = program->checkParam("--show_optimizationfast");
	min_freq = program->getDoubleParam("--min_freqfast");
	max_freq = program->getDoubleParam("--max_freqfast");
	defocus_range = program->getDoubleParam("--defocus_rangefast");
	modelSimplification = program->getIntParam("--model_simplification");
	bootstrap = program->checkParam("--bootstrapFitfast");
	fastDefocusEstimate = program->checkParam("--fastDefocusfast");
	activar1D = program->checkParam("--activar1D");
	if (fastDefocusEstimate)
	{
		lambdaPhase=program->getDoubleParam("--fastDefocusfast",0);
		sizeWindowPhase=program->getIntParam("--fastDefocusfast",1);
	}
	ctfmodelSize = program->getIntParam("--ctfmodelSizefast");
	enhanced_weight = program->getDoubleParam("--enhance_weightfast");
	if (!program->checkParam("--enhance_min_freqfast"))
		f1 = (max_freq > 0.35) ? 0.01 : 0.02;
	else
		f1 = program->getDoubleParam("--enhance_min_freqfast");
	if (!program->checkParam("--enhance_max_freqfast"))
		f2 = (max_freq > 0.35) ? 0.08 : 0.15;
	else
		f2 = program->getDoubleParam("--enhance_max_freqfast");

	initial_ctfmodel.enable_CTF = initial_ctfmodel.enable_CTFnoise = true;
	initial_ctfmodel.readParams(program);
	if (initial_ctfmodel.Defocus>100e3) //Nuevo
		REPORT_ERROR(ERR_ARG_INCORRECT,"Defocus cannot be larger than 10 microns (100,000 Angstroms)");
	Tm = initial_ctfmodel.Tm;
}

void ProgCTFEstimateFromPSDFast::readParams()
{
	fn_psd = getParam("--psd");
	readBasicParams(this);
}

/* Show -------------------------------------------------------------------- */
void ProgCTFEstimateFromPSDFast::show()
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
		std::cout << "Starting CTF:\n" << initial_ctfmodel << std::endl;
}

/* Usage ------------------------------------------------------------------- */
void ProgCTFEstimateFromPSDFast::defineBasicParams(XmippProgram * program)
{
	program->addSeeAlsoLine("ctf_enhance_psd");
	    program->addParamsLine("== Downsampling1D");
	    program->addParamsLine(
	        "  [--downSamplingPerformedfast <F=1>] : Downsampling performed to produce this PSD");
	    program->addParamsLine(
	        "                                  : Note that the output CTF model is referred to the original sampling rate");
	    program->addParamsLine(
	        "                                  : not the one of the downsampled image.");
	    program->addParamsLine("== CTF fit: Optimization constraints");
	    program->addParamsLine(
	        "   [--min_freqfast <fmin=0.03>]     : Minimum digital frequency (<0.5) to use in adjust. Its value");
	    program->addParamsLine(
	        "                                : should be a little lower than the dig. freq. of the first ");
	    program->addParamsLine("                                : CTF zero.");
	    program->addParamsLine(
	        "   [--max_freqfast <fmax=0.35>]     : Maximum digital frequency (<0.5) to use in adjust.");
	    program->addParamsLine(
	        "                                : It should be higher than the last zero of the CTF.");
	    program->addParamsLine(
	        "   [--fastDefocusfast <lambda=2> <size=10>] : Estimate first defocus with Zernike polynomials");
	    program->addParamsLine(
	        "                                :+Lambda is a regularization factor used during the estimation of the CTF phase");
	    program->addParamsLine(
	        "                                :+During the estimation, the phase values are averaged within a window of this size");
	    program->addParamsLine(
	        "   [--defocus_rangefast <D=8000>]   : Defocus range in Angstroms");
	    program->addParamsLine(
	        "   [--show_optimizationfast++]       : Show optimization process");
	    program->addParamsLine(
	        "   [--radial_noisefast++]           : By default, noise is astigmatic");
	    program->addParamsLine(
	        "   [--enhance_weightfast++ <w=1>]   : Weight of the enhanced term");
	    program->addParamsLine(
	        "   [--model_simplificationfast++ <s=0>]: 0 (no simplification)");
	    program->addParamsLine(
	        "                                : 1 (simplified envelope)");
	    program->addParamsLine(
	        "                                : 2 (last Gaussian removal)");
	    program->addParamsLine(
	        "                                : 3 (symmetric intermediate Gaussian)");
	    program->addParamsLine(
	        "   [--bootstrapFitfast++ <N=-1>]    : Perform bootstrap fit (Fourier pixels are randomly chosen)");
	    program->addParamsLine(
	        "                                : This is used to test the variability of the fit");
	    program->addParamsLine(
	        "                                : N defines the number of times the fit is repeated");
	    program->addParamsLine("==+ CTF fit: Output CTF models");
	    program->addParamsLine(
	        "   [--ctfmodelSizefast <size=256>]  : Size for the ctfmodel thumbnails");
	    program->addParamsLine("==+ PSD enhancement");
	    program->addParamsLine(
	        "   [--enhance_min_freqfast <f1>]    : Bandpass cutoff. Normalized to 0.5");
	    program->addParamsLine(
	        "                                : If fmax>0.35, f1 default=0.01");
	    program->addParamsLine(
	        "                                : If fmax<0.35, f1 default=0.02");
	    program->addParamsLine(
	        "   [--enhance_max_freqfast <f2>]    : Bandpass cutoff. Normalized to 0.5.");
	    program->addParamsLine(
	        "                                : If fmax>0.35, f2 default=0.08");
	    program->addParamsLine(
	        "                                : If fmax<0.35, f2 default=0.15");
	    program->addParamsLine(
	    	        "   [--activar1D]       : activar 1D");
	    CTF1D::defineParams(program);
}

void ProgCTFEstimateFromPSDFast::defineParams()
{
	addUsageLine("Adjust a parametric model to a PSD file.");
	addUsageLine("The PSD is enhanced ([[http://www.ncbi.nlm.nih.gov/pubmed/16987671][See article]]). ");
	addUsageLine("And finally, the CTF is fitted to the PSD, being guided by the enhanced PSD ");
	addUsageLine("([[http://www.ncbi.nlm.nih.gov/pubmed/17911028][See article]]).");
	addParamsLine("   --psd <PSDfile> : PSD file");
	addSeeAlsoLine("ctf_estimate_from_micrograph");
	defineBasicParams(this);
}

/* Produce side information ------------------------------------------------ */
void ProgCTFEstimateFromPSDFast::produceSideInfo()
{
    adjust.resize(ALL_CTF_PARAMETERS);
    adjust.initZeros();
    global_ctfmodel.clear();
    global_ctfmodel_defoci.clear();
    assignParametersFromCTF(initial_ctfmodel, MATRIX1D_ARRAY(adjust), 0, ALL_CTF_PARAMETERS, true);

    // Read the CTF file, supposed to be the uncentered squared amplitudes

    if (fn_psd != "")
        ctftomodel.read(fn_psd);
    f = &(ctftomodel());

    // Resize the frequency
    global_x_digfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    global_y_digfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    global_w_digfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    global_w_digfreq_r.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    global_x_contfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    global_y_contfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);
    global_w_contfreq.initZeros(YSIZE(*f), XSIZE(*f) / 2);

    Matrix1D<int> idx(2); // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane

    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_x_digfreq)
    {
        XX(idx) = j;
        YY(idx) = i;

        // Digital frequency
        FFT_idx2digfreq(*f, idx, freq);
        //std::cout<< "freq =" << freq << std::endl;
        global_x_digfreq(i, j) = XX(freq);
        global_y_digfreq(i, j) = YY(freq);
        global_w_digfreq(i, j) = freq.module();
        global_w_digfreq_r(i, j) = (int)(global_w_digfreq(i,j) * (double)YSIZE(global_w_digfreq));

        // Continuous frequency
        digfreq2contfreq(freq, freq, global_prm->Tm);
        global_x_contfreq(i, j) = XX(freq);
        global_y_contfreq(i, j) = YY(freq);
        global_w_contfreq(i, j) = freq.module();
    }

    // Precompute frequency related terms in the CTF
    global_ctfmodel.precomputeValues(global_x_contfreq);

    // Build frequency mask
    global_mask.initZeros(global_w_digfreq);
    global_w_count.initZeros(XSIZE(global_w_digfreq));
    //printf("max_freq = %lf\n",max_freq);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_w_digfreq(i, j) >= max_freq
            || global_w_digfreq(i, j) <= min_freq)
            continue;
        global_mask(i, j) = 1;
        global_w_count(global_w_digfreq_r(i, j))++;
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
	enhanced_ctftomodel().resize(global_w_digfreq);

	// Divide by the number of count at each frequency
	// and mask between min_freq and max_freq
	double min_val = enhanced_ctftomodel().computeMin();
	FOR_ALL_ELEMENTS_IN_ARRAY2D(global_mask)
	if (global_mask(i, j) <= 0)
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

	// Compute now radial average of the enhanced_ctftomodel
	global_psd_exp_enhanced_radial_derivative.initZeros(XSIZE(enhanced_ctftomodel()));
	global_psd_theo_radial_derivative.initZeros(global_psd_exp_enhanced_radial_derivative);
	global_psd_exp_enhanced_radial.initZeros(global_psd_exp_enhanced_radial_derivative);
	global_psd_exp_radial.initZeros(global_psd_exp_enhanced_radial_derivative);
	global_psd_theo_radial.initZeros(global_psd_exp_enhanced_radial_derivative);
	global_w_digfreq_r_iN.initZeros(global_psd_exp_enhanced_radial_derivative);
	FOR_ALL_ELEMENTS_IN_ARRAY2D(enhanced_ctftomodel())
	{
		if (global_w_digfreq(i,j)>min_freq && global_w_digfreq(i,j)<max_freq)
		{
			int r = global_w_digfreq_r(i, j);
			global_w_digfreq_r_iN(r)+=1;
			global_psd_exp_enhanced_radial(r) += enhanced_ctftomodel(i, j);
			global_psd_exp_radial(r) += ctftomodel(i, j);
		}
	}
	FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq_r_iN)
	if (global_w_digfreq_r_iN(i)>0)
	{
		global_w_digfreq_r_iN(i)=1.0/global_w_digfreq_r_iN(i);
		global_psd_exp_enhanced_radial(i)*=global_w_digfreq_r_iN(i);
		global_psd_exp_radial(i)*=global_w_digfreq_r_iN(i);
	}

	// Compute its derivative
	int state=0;
	double maxDiff=0;
	FOR_ALL_ELEMENTS_IN_ARRAY1D(global_psd_exp_enhanced_radial)
	{
		switch (state)
		{
		case 0:
			if (global_w_digfreq(i,0)>min_freq)
				state=1;
			break;
		case 1:
			state=2; // Skip first sample
			break;
		case 2:
			if (global_w_digfreq(i,0)>max_freq)
				state=3;
			else
			{
				double diff=global_psd_exp_enhanced_radial(i)-global_psd_exp_enhanced_radial(i-1);
				global_psd_exp_enhanced_radial_derivative(i)=diff;
				maxDiff=std::max(maxDiff,fabs(diff));
			}
			break;
		}
	}
	//printf("maxdiff = %lf\n",maxDiff);
	global_psd_exp_enhanced_radial_derivative/=maxDiff;
	global_psd_exp_enhanced_radial.write("PPPexpEnhanced_fast.txt");
	global_psd_exp_radial.write("PPPexp_fast.txt");
}

/* Generate model at a given size ------------------------------------------ */
void ProgCTFEstimateFromPSDFast::generate_model_quadrant(int Ydim, int Xdim,
        MultidimArray<double> &model)
{

}

void ProgCTFEstimateFromPSDFast::generate_model_halfplane(int Ydim, int Xdim,
        MultidimArray<double> &model)
{
}

/*void printconsole()
{
	std::cout
	<< "   sq=" << global_ctfmodel.sq << "\n"
	<< "   c=" << global_ctfmodel.Gc1 << "\n"
	<< "   sigma=" << global_ctfmodel.sigma1 << "\n"
	<< "   c2=" << global_ctfmodel.Gc2 << "\n"
	<< "   sigma2=" << global_ctfmodel.sigma2 << "\n"
	<< "   gaussian_K=" << global_ctfmodel.gaussian_K << "\n"
	<< "   sqrt_K=" << global_ctfmodel.sqrt_K << "\n"
	<< "   base_line=" << global_ctfmodel.base_line << "\n";
	std::cout << "   u=" << global_ctfmodel.precomputed.u << "u_sqrt=" << global_ctfmodel.precomputed.u_sqrt <<") CTFnoise="
	<< global_ctfmodel.base_line +
	global_ctfmodel.gaussian_K*exp(-global_ctfmodel.sigma1*(global_ctfmodel.precomputed.u - global_ctfmodel.Gc1)*(global_ctfmodel.precomputed.u - global_ctfmodel.Gc1)) +
	global_ctfmodel.sqrt_K*exp(-global_ctfmodel.sq*sqrt(global_ctfmodel.precomputed.u)) -
	global_ctfmodel.gaussian_K2*exp(-global_ctfmodel.sigma2*(global_ctfmodel.precomputed.u - global_ctfmodel.Gc2)*(global_ctfmodel.precomputed.u - global_ctfmodel.Gc2)) << std::endl;
}*/

void generateModelSoFar_fast(MultidimArray<double> &I, bool apply_log = false)
{
    Matrix1D<int> idx(1); // Indexes for Fourier plane
    Matrix1D<double> freq(1); // Frequencies for Fourier plane

    assignCTFfromParameters(MATRIX1D_ARRAY(*global_adjust), global_ctfmodel,
                            0, ALL_CTF_PARAMETERS, global_prm->modelSimplification);
    global_ctfmodel.produceSideInfo();

    I.initZeros(global_psd_exp_radial);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(I)
    {

        XX(idx) = i;
        FFT_idx2digfreq(*f, idx, freq);
        double w=freq.module();
        if (w>global_max_freq)
        	continue;
        digfreq2contfreq(freq, freq, global_prm->Tm);

        // Decide what to save
        global_ctfmodel.precomputeValues(XX(freq));
        if (global_action <= 1)
            I(i) = global_ctfmodel.getValueNoiseAt();
        else if (global_action == 2)
        {
            double E = global_ctfmodel.getValueDampingAt();
            I(i) = global_ctfmodel.getValueNoiseAt() + E * E;
        }
        else if (global_action >= 3 && global_action <= 6)
        {
            double ctf = global_ctfmodel.getValuePureAt();
            I(i) = global_ctfmodel.getValueNoiseAt() + ctf * ctf;
        }
        else
        {
            double ctf = global_ctfmodel.getValuePureAt();
            I(i) = ctf;
        }
        if (apply_log)
            I(i) = 10 * log10(I(i));

    }
    //std::cout << "***" << I << std::endl;
}

void saveIntermediateResults_fast(const FileName &fn_root, bool generate_profiles =
                                 true)
{

    std::ofstream plot_radial;
    MultidimArray<double> save;
    generateModelSoFar_fast(save, false);

    if (!generate_profiles)
        return;
    plot_radial.open((fn_root + "_radial.txt").c_str());
    if (!plot_radial)
        REPORT_ERROR(
            ERR_IO_NOWRITE,
            "save_intermediate_results::Cannot open plot file for writing\n");
    //plotX << "# freq_dig freq_angstrom model psd enhanced logModel logPsd\n";
    //plotY << "# freq_dig freq_angstrom model psd enhanced logModel logPsd\n";
    //plot_radial << "# freq_dig freq_angstrom model psd enhanced logModel logPsd\n";

    // Generate radial average
    MultidimArray<double> radial_CTFmodel_avg;
    radial_CTFmodel_avg.initZeros(global_psd_exp_enhanced_radial);

    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)
    {
        if (global_mask(i) <= 0)
            continue;
        double model2 = save(i);

        int r = global_w_digfreq_r(i);
        radial_CTFmodel_avg(r) = model2;

        plot_radial << global_w_digfreq(r, 0) << " " << global_w_contfreq(r, 0)
		<< " " << radial_CTFmodel_avg(r) << " "
		<< global_psd_exp_enhanced_radial(r) << " " << global_psd_exp_radial(r)
		<< " " << log10(radial_CTFmodel_avg(r)) << " "
		<< std::endl;
    }

    plot_radial.close();
}

/* CTF fitness ------------------------------------------------------------- */
/* This function measures the distance between the estimated CTF and the
 measured CTF */
double CTF_fitness_fast(double *p, void *)
{

	double retval;

    // Generate CTF model
    switch (global_action)
    {
        // Remind that p is a vector whose first element is at index 1
    case 0:
        assignCTFfromParameters(p - FIRST_SQRT_PARAMETER + 1,
                                global_ctfmodel, FIRST_SQRT_PARAMETER, SQRT_CTF_PARAMETERS,
                                global_prm->modelSimplification);
        if (global_show >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= SQRT_CTF_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 1:
            assignCTFfromParameters(p - FIRST_SQRT_PARAMETER + 1,
                                    global_ctfmodel, FIRST_SQRT_PARAMETER,
                                    BACKGROUND_CTF_PARAMETERS, global_prm->modelSimplification);
        if (global_show >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= BACKGROUND_CTF_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 2:
            assignCTFfromParameters(p - FIRST_ENVELOPE_PARAMETER + 1,
                                    global_ctfmodel, FIRST_ENVELOPE_PARAMETER, ENVELOPE_PARAMETERS,
                                    global_prm->modelSimplification);
        if (global_show >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= ENVELOPE_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 3:
            assignCTFfromParameters(p - FIRST_DEFOCUS_PARAMETER + 1,
                                    global_ctfmodel, FIRST_DEFOCUS_PARAMETER, DEFOCUS_PARAMETERS,
                                    global_prm->modelSimplification);
        global_psd_theo_radial_derivative.initZeros();
        if (global_show >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= DEFOCUS_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 4:
            assignCTFfromParameters(p - 0 + 1, global_ctfmodel, 0,
                                    CTF_PARAMETERS, global_prm->modelSimplification);
        global_psd_theo_radial.initZeros();
        if (global_show >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= CTF_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 5:
    case 6:
    case 7:
		assignCTFfromParameters(p - 0 + 1, global_ctfmodel, 0,
								ALL_CTF_PARAMETERS, global_prm->modelSimplification);
        global_psd_theo_radial.initZeros();
        if (global_show >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= ALL_CTF_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    }
    global_ctfmodel.produceSideInfo();
    if (global_show >= 2)
        std::cout << "Model:\n" << global_ctfmodel << std::endl;
    if (!global_ctfmodel.hasPhysicalMeaning())
    {
        if (global_show >= 2)
            std::cout << "Does not have physical meaning\n";
        return global_heavy_penalization;
    }
    if (global_action > 3
        && (fabs((global_ctfmodel.Defocus - global_ctfmodel_defoci.Defocus)
                / global_ctfmodel_defoci.Defocus) > 0.2))
    {
        if (global_show >= 2)
            std::cout << "Too large defocus\n";
        return global_heavy_penalization;
    }

    if (global_prm->initial_ctfmodel.Defocus != 0 && global_action >= 3)
    {
        // If there is an initial model, the true solution
        // cannot be too far
        if (fabs(global_prm->initial_ctfmodel.Defocus - global_ctfmodel.Defocus) > global_prm->defocus_range)
        {
            if (global_show >= 2)
            {
                std::cout << "Too far from hint: Initial (" << global_prm->initial_ctfmodel.Defocus << ")"
                << " current guess (" << global_ctfmodel.Defocus << ") max allowed difference: "
                << global_prm->defocus_range << std::endl;
            }
            return global_heavy_penalization;

        }
    }

    // Now the 1D error
    double distsum = 0;
    int N = 0, Ncorr = 0;
    double enhanced_avg = 0;
    double model_avg = 0;
    double enhanced_model = 0;
    double enhanced2 = 0;
    double model2 = 0;
    double lowerLimit = 1.1 * global_min_freq;
    double upperLimit = 0.9 * global_max_freq;
    const MultidimArray<double>& local_enhanced_ctf = global_psd_exp_enhanced_radial;
    int XdimW=XSIZE(global_w_digfreq);
    global_corr13=0;
    int cont = 0;
    //printf("base_line = %lf\n",global_ctfmodel.base_line);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)
    //for (int i = 0; i < XdimW; i += global_evaluation_reduction)
    {
		if (DIRECT_A1D_ELEM(global_mask, i) <= 0)
			continue;

		// Compute each component
		global_ctfmodel.precomputeValues(i);
		double bg = global_ctfmodel.getValueNoiseAt();

		//printf("bg = %lf\n", bg);
		double envelope=0, ctf_without_damping, ctf_with_damping=0;
		double ctf2_th=0;
		switch (global_action)
		{
		case 0:
		case 1:
			ctf2_th = bg;
			break;
		case 2:
			envelope = global_ctfmodel.getValueDampingAt();
			ctf2_th = bg + envelope * envelope;
			break;
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:

			envelope = global_ctfmodel.getValueDampingAt();
			ctf_without_damping = global_ctfmodel.getValuePureWithoutDampingAt();
			ctf_with_damping = envelope * ctf_without_damping;
			ctf2_th = bg + ctf_with_damping * ctf_with_damping;
			//printf("envelope = %lf nodamping = %lf damping = %lf ctf2_th = %lf bg = %lf cont =%i\n",envelope,ctf_without_damping,ctf_with_damping,ctf2_th,bg,cont);
			break;

		}
		// Compute distance
		double ctf2 = DIRECT_A1D_ELEM(global_psd_exp_radial, i);
		double dist = 0;
		double ctf_with_damping2;
		switch (global_action)
		{
		case 0:
		case 1:
			dist = fabs(ctf2 - bg);
			if (global_penalize && bg > ctf2
				&& DIRECT_A1D_ELEM(global_w_digfreq, i)
				> global_max_gauss_freq)
				dist *= global_current_penalty;
			break;
		case 2:
			dist = fabs(ctf2 - ctf2_th);
			if (global_penalize && ctf2_th < ctf2
				&& DIRECT_A1D_ELEM(global_w_digfreq, i)
				> global_max_gauss_freq)
				dist *= global_current_penalty;
			break;
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
			if (DIRECT_A1D_ELEM(global_w_digfreq,i) < upperLimit
				&& DIRECT_A1D_ELEM(global_w_digfreq,i) > lowerLimit)
			{
				if  (global_action == 3 ||
					 (global_action == 4 && DIRECT_A1D_ELEM(global_mask_between_zeroes,i) == 1) ||
					 (global_action == 7 && DIRECT_A1D_ELEM(global_mask_between_zeroes,i) == 1))
				{
					double enhanced_ctf = DIRECT_A1D_ELEM(local_enhanced_ctf, i);
					ctf_with_damping2 = ctf_with_damping * ctf_with_damping;
					enhanced_model += enhanced_ctf * ctf_with_damping2;
					enhanced2 += enhanced_ctf * enhanced_ctf;
					model2 += ctf_with_damping2 * ctf_with_damping2;
					enhanced_avg += enhanced_ctf;
					model_avg += ctf_with_damping2;
					Ncorr++;

					if (global_action==3)
					{
						int r = A1D_ELEM(global_w_digfreq_r,i);
						A1D_ELEM(global_psd_theo_radial,r) = ctf2_th;
					}
				}
			}
			if (envelope > 1e-2)
				dist = fabs(ctf2 - ctf2_th) / (envelope * envelope);
			else
				dist = fabs(ctf2 - ctf2_th);
			break;
		}
		distsum += dist * DIRECT_A1D_ELEM(global_mask,i);
		N++;
	}

    if (N > 0)
    { retval = distsum/N ;
    }
    else
        retval = global_heavy_penalization;
    if (global_show >=2)
        std::cout << "Fitness1=" << retval << std::endl;
    if ( (((global_action >= 3) && (global_action <= 4)) || (global_action == 7))
         && (Ncorr > 0) && (global_prm->enhanced_weight != 0) )
    {

    	model_avg /= Ncorr;
    	enhanced_avg /= Ncorr;
        double correlation_coeff = enhanced_model/Ncorr - model_avg * enhanced_avg;
        //printf("correlation_coeff = %lf\n",correlation_coeff);
        double sigma1 = sqrt(fabs(enhanced2/Ncorr - enhanced_avg * enhanced_avg));
        double sigma2 = sqrt(fabs(model2/Ncorr - model_avg * model_avg));
        double maxSigma = std::max(sigma1, sigma2);
        if (sigma1 < XMIPP_EQUAL_ACCURACY || sigma2 < XMIPP_EQUAL_ACCURACY
            || (fabs(sigma1 - sigma2) / maxSigma > 0.9 && global_action>=5))
        {
            retval = global_heavy_penalization;
            if (global_show>=2)
                std::cout << "Fitness2=" << global_heavy_penalization << " sigma1=" << sigma1 << " sigma2=" << sigma2 << std::endl;
        }
        else
        {
            correlation_coeff /= sigma1 * sigma2;
            printf("correlation_coeff = %lf\n",correlation_coeff);
            if (global_action == 7)
                global_corr13 = correlation_coeff;
            else
                retval -= global_prm->enhanced_weight * correlation_coeff;
            if (global_show >= 2)
            {
                std::cout << "model_avg=" << model_avg << std::endl;
                std::cout << "enhanced_avg=" << enhanced_avg << std::endl;
                std::cout << "enhanced_model=" << enhanced_model / Ncorr
                << std::endl;
                std::cout << "sigma1=" << sigma1 << std::endl;
                std::cout << "sigma2=" << sigma2 << std::endl;
                std::cout << "Fitness2="
                << -(global_prm->enhanced_weight * correlation_coeff)
                << " (" << correlation_coeff << ")" << std::endl;
            }
        }

        // Correlation of the derivative of the radial profile
        if (global_action==3 || global_evaluation_reduction==1)
        {
            int state=0;
            double maxDiff=0;
            global_psd_theo_radial_derivative.initZeros();
            double lowerlimt=1.1*global_prm->min_freq;
            double upperlimit=0.9*global_prm->max_freq;
            FOR_ALL_ELEMENTS_IN_ARRAY1D(global_psd_theo_radial)
            if (A1D_ELEM(global_w_digfreq_r_iN,i)>0)
            {
                double freq=A1D_ELEM(global_w_digfreq,i);
                switch (state)
                {
                case 0:
                    if (freq>lowerlimt)
                        state=1;
                    break;
                case 1:
                    if (freq>upperlimit)
                        state=2;
                    else
                    {
                        double diff=A1D_ELEM(global_psd_theo_radial,i)-A1D_ELEM(global_psd_theo_radial,i-1);

                        A1D_ELEM(global_psd_theo_radial_derivative,i)=diff;
                        maxDiff=std::max(maxDiff,fabs(diff));
                        //printf("maxdiff = %lf\n",maxDiff);
                    }
                    break;
                }
            }
            double corrRadialDerivative=0,mux=0, muy=0, Ncorr=0, sigmax=0, sigmay=0;
            double iMaxDiff=1.0/maxDiff;

            FOR_ALL_ELEMENTS_IN_ARRAY1D(global_psd_theo_radial)
            {
                A1D_ELEM(global_psd_theo_radial_derivative,i)*=iMaxDiff;
                double x=A1D_ELEM(global_psd_exp_enhanced_radial_derivative,i);
                double y=A1D_ELEM(global_psd_theo_radial_derivative,i);
                corrRadialDerivative+=x*y;
                mux+=x;
                muy+=y;
                sigmax+=x*x;
                sigmay+=y*y;
                Ncorr++;
            }

            if (Ncorr>0)
            {
                double iNcorr=1.0/Ncorr;
                corrRadialDerivative*=iNcorr;
                mux*=iNcorr;
                muy*=iNcorr;
                sigmax=sqrt(fabs(sigmax*iNcorr-mux*mux));
                sigmay=sqrt(fabs(sigmay*iNcorr-muy*muy));
                corrRadialDerivative=(corrRadialDerivative-mux*muy)/(sigmax*sigmay);
            }
            retval-=corrRadialDerivative;
            if (global_show>=2)
            {
                std::cout << "Fitness3=" << -corrRadialDerivative << std::endl;
                if (global_show==3)
                {
                    global_psd_exp_enhanced_radial.write("PPPexpRadial_fast.txt");
                    global_psd_theo_radial.write("PPPtheoRadial.txt");
                    global_psd_exp_enhanced_radial_derivative.write("PPPexpRadialDerivative_fast.txt");
                    global_psd_theo_radial_derivative.write("PPPtheoRadialDerivative.txt");
                }
           }
        }
    }


    return retval;
}

/* Center focus ----------------------------------------------------------- */
void center_optimization_focus_fast(bool adjust_freq, bool adjust_th, double margin = 1)
{
    if (global_prm->show_optimization)
        std::cout << "Freq frame before focusing=" << global_min_freq << ","
        << global_max_freq << std::endl << "Value_th before focusing="
        << global_value_th << std::endl;

    double w1 = global_min_freq, w2 = global_max_freq;
    if (adjust_freq)
    {
        double w1U, w2U;
        //compute_central_region(w1U, w2U, global_ctfmodel.azimuthal_angle);
        //compute_central_region(w1V, w2V, global_ctfmodel.azimuthal_angle + 90);
        w1 = XMIPP_MIN(w1U, w1U);
        w2 = XMIPP_MAX(w2U, w2U);
        global_min_freq = XMIPP_MAX(global_min_freq, w1 - 0.05);
        global_max_freq = XMIPP_MIN(global_max_freq, w2 + 0.01);
    }

    // Compute maximum value within central region
    if (adjust_th)
    {
        MultidimArray<double> save;
        generateModelSoFar_fast(save);
        double max_val = 0;
        FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)
        {
            double w = global_w_digfreq(i);
            if (w >= w1 && w <= w2)
                max_val = XMIPP_MAX(max_val, save(i));
        }
        if (global_value_th != -1)
            global_value_th = XMIPP_MIN(global_value_th, max_val * margin);
        else
            global_value_th = max_val * margin;
    }


    if (global_prm->show_optimization)
        std::cout << "Freq frame after focusing=" << global_min_freq << ","
        << global_max_freq << std::endl << "Value_th after focusing="
        << global_value_th << std::endl;
}
// Estimate sqrt parameters ------------------------------------------------
// Results are written in global_ctfmodel
void estimate_background_sqrt_parameters_fast()
{
    if (global_prm->show_optimization)
        std::cout << "Computing first sqrt background ...\n";

    // Estimate the base line taking the value of the CTF
    // for the maximum X
    double base_line = 0;
    int N = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)  //global_w_digfreq comprobado y funciona bien
    if (global_w_digfreq(i) > 0.4)
    {
        N++;
        base_line += global_psd_exp_radial(i);	//global_psd_exp_enhanced_radial
    }
    global_ctfmodel.base_line = base_line / N;

    // Find the linear least squares solution for the sqrt part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();

    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)
    {
    	//printf("w = %lf ",global_w_digfreq(i,0));
        if (global_mask(i) <= 0)
            continue;

        // Compute weight for this point
        double weight = 1 + global_max_freq - global_w_digfreq(i);

        // Compute error
        global_ctfmodel.precomputeValues(global_x_contfreq(i));
        double explained = global_ctfmodel.getValueNoiseAt(); //PSDna: ruido despues del proceso.
        double unexplained = global_psd_exp_radial(i) - explained;
        if (unexplained <= 0)
            continue;
        unexplained = log(unexplained);

        double X = -sqrt(global_w_contfreq(i));
        A(0, 0) += weight * X * X;
        A(0, 1) += weight * X;
        A(1, 1) += weight * 1;
        b(0) += X * weight * unexplained;
        b(1) += weight * unexplained;

    }

    A(1, 0) = A(0, 1);
    printf("%lf  %lf\n %lf  %lf\n",A(0,0),A(0,1),A(1,0),A(1,1));
    b = A.inv() * b;
    printf("b(0) = %lf\n, b(1) = %lf \n",b(0),b(1));

    global_ctfmodel.sq = b(0);
    global_ctfmodel.sqrt_K = exp(b(1));
    COPY_ctfmodel_TO_CURRENT_GUESS;

    //if (global_prm->show_optimization)
    //{
        std::cout << "First SQRT Fit:\n" << global_ctfmodel << std::endl;
        //printf("Save results\n");
        saveIntermediateResults_fast("step01a_first_sqrt_fit_fast");
    //}

    // Now optimize .........................................................
	double fitness;
	Matrix1D<double> steps;
	steps.resize(SQRT_CTF_PARAMETERS);
	steps.initConstant(1);

	// Optimize without penalization
	if (global_prm->show_optimization)
		std::cout << "Looking for best fitting sqrt ...\n";
	global_penalize = false;
	int iter;
	powellOptimizer(*global_adjust, FIRST_SQRT_PARAMETER + 1,
					SQRT_CTF_PARAMETERS, &CTF_fitness_fast, NULL, 0.05, fitness, iter, steps,
					global_prm->show_optimization);

	// Optimize with penalization
	if (global_prm->show_optimization)
		std::cout << "Penalizing best fitting sqrt ...\n";
	global_penalize = true;
	global_current_penalty = 2;
	int imax = CEIL(log(global_penalty) / log(2.0));

	for (int i = 1; i <= imax; i++)
	{
		if (global_prm->show_optimization)
			std::cout << "     Iteration " << i << " penalty="
			<< global_current_penalty << std::endl;
		powellOptimizer(*global_adjust, FIRST_SQRT_PARAMETER + 1,
						SQRT_CTF_PARAMETERS, &CTF_fitness_fast, NULL, 0.05, fitness, iter,
						steps, global_prm->show_optimization);
		global_current_penalty *= 2;
		global_current_penalty =
			XMIPP_MIN(global_current_penalty, global_penalty);
	}
	// Keep the result in global_prm->adjust
	global_ctfmodel.forcePhysicalMeaning();
	COPY_ctfmodel_TO_CURRENT_GUESS;

	if (global_prm->show_optimization)
	{
		std::cout << "Best penalized SQRT Fit:\n" << global_ctfmodel
		<< std::endl;
		saveIntermediateResults_fast("step01b_best_penalized_sqrt_fit_fast");
	}

   center_optimization_focus_fast(false, true, 1.5);

}

// Estimate gaussian parameters --------------------------------------------
//#define DEBUG
void estimate_background_gauss_parameters_fast()
{

    if (global_prm->show_optimization)
        std::cout << "Computing first background Gaussian parameters ...\n";

    // Compute radial averages
    MultidimArray<double> radial_CTFmodel_avg(YSIZE(*f) / 2);
    MultidimArray<double> radial_CTFampl_avg(YSIZE(*f) / 2);
    MultidimArray<int> radial_N(YSIZE(*f) / 2);
    double w_max_gauss = 0.25;

    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)
    {
        if (global_mask(i) <= 0)
            continue;
        double w = global_w_digfreq(i);

        if (w > w_max_gauss)
            continue;

        int r = FLOOR(w * (double)YSIZE(*f));

        global_ctfmodel.precomputeValues(global_x_contfreq(i));
        radial_CTFmodel_avg(r) += global_ctfmodel.getValueNoiseAt();
        radial_CTFampl_avg(r) += global_psd_exp_radial(i);
        radial_N(r)++;

    }
    //radial_CTFmodel_avg.write("radial_CTFmodel_avg_fast.txt");
    //radial_CTFampl_avg.write("radial_CTFampl_avg_fast.txt");
    //radial_N.write("radial_N_fast.txt");

    // Compute the average radial error
    double  error2_avg = 0;
    int N_avg = 0;
    MultidimArray<double> error;
    error.initZeros(radial_CTFmodel_avg);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        error(i) = (radial_CTFampl_avg(i) - radial_CTFmodel_avg(i))
                   / radial_N(i);
        error2_avg += error(i) * error(i);
        N_avg++;
    }
    if (N_avg != 0)
        error2_avg /= N_avg;

#ifdef DEBUG

    std::cout << "Error2 avg=" << error2_avg << std::endl;
#endif

    // Compute the minimum radial error
    bool first = true, OK_to_proceed = false;
    double error2_min = 0, wmin=0;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;

        double w = global_w_digfreq(i);

        if (error(i) < 0 && first){
        	continue;
        }

        else if (error(i) < 0)
        {
        	break;
        }

        double error2 = error(i) * error(i);
        // If the two lines cross, do not consider any error until
        // the cross is "old" enough
        if (first && error2 > 0.15 * error2_avg)
            OK_to_proceed = true;
        if (first && i > 0)
            OK_to_proceed &= (error(i) < error(i - 1));

        // If the error now is bigger than a 30% (1.69=1.3*1.3) of the error min
        // this must be a rebound. Stop here
        if (!first && error2 > 1.69 * error2_min)
            break;
        if (first && OK_to_proceed)
        {
            wmin = w;
            error2_min = error2;
            first = false;
        }
        if (!first && error2 < error2_min)
        {
            wmin = w;
            error2_min = error2;
        }
#ifdef DEBUG
        std::cout << w << " " << error2 << " " << wmin << " " << std::endl;
#endif

    }

    // Compute the frequency of the minimum error
    global_max_gauss_freq = wmin;
#ifdef DEBUG

    std::cout << "Freq of the minimum error: " << wmin << " " << fmin << std::endl;
#endif

    // Compute the maximum radial error
    first = true;
    double error2_max = 0, wmax=0, fmax;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        double w = global_w_digfreq(i);
        if (w > wmin)
            continue;

        if (error(i) < 0 && first)
            continue;
        else if (error(i) < 0)
            break;
        double error2 = error(i) * error(i);
        if (first)
        {
            wmax = w;
            error2_max = error2;
            first = false;
        }
        if (error2 > error2_max)
        {
            wmax = w;
            error2_max = error2;
        }
#ifdef DEBUG
        std::cout << w << " " << error2 << " " << wmax << std::endl;
#endif

    }
    fmax = global_ctfmodel.Gc1 = wmax / global_prm->Tm;
#ifdef DEBUG

    std::cout << "Freq of the maximum error: " << wmax << " " << fmax << std::endl;
#endif

    // Find the linear least squares solution for the gauss part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();
    double totalunexplained;
    double totalweight;

    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)
    {

        if (global_mask(i) <= 0)
            continue;

        if (global_w_digfreq(i) > wmin)
        	continue;

        double fmod = global_w_contfreq(i);

        // Compute weight for this point
        double weight = 1 + global_max_freq - global_w_digfreq(i);

        // Compute error
        global_ctfmodel.precomputeValues(global_x_contfreq(i));
        double explained = global_ctfmodel.getValueNoiseAt();
        double unexplained = global_psd_exp_radial(i) - explained;
        if (unexplained <= 0)
            continue;
        unexplained = log(unexplained);
        double F = -(fmod - fmax) * (fmod - fmax);

        A(0, 0) += weight * 1;
        A(0, 1) += weight * F;
        A(1, 1) += weight * F * F;
        b(0) += weight * unexplained;
        b(1) += F * weight * unexplained;

        totalunexplained += unexplained;
        totalweight += weight;

    }
    A(1, 0) = A(0, 1);
    printf("A=\n %lf %lf\n %lf %lf\n",A(0,0),A(0,1),A(1,0),A(1,1));
    printf("b(0) = %lf\n, b(1) = %lf \n",b(0),b(1));

    if ( (A(0, 0)== 0) && (A(1, 0)== 0) && (A(1, 1)== 0))
    {
        std::cout << "Matriz A es zero" << std::endl;
    }
    else
    {
        b = A.inv() * b;
        printf("b(0) = %lf\n, b(1) = %lf \n",b(0),b(1));

        global_ctfmodel.sigma1 = XMIPP_MIN(fabs(b(1)), 95e3); // This value should beconformant with the physical

        // meaning routine in CTF.cc
        global_ctfmodel.gaussian_K = exp(b(0));
        // Store the CTF values in global_prm->adjust
        global_ctfmodel.forcePhysicalMeaning();
        COPY_ctfmodel_TO_CURRENT_GUESS;

        //if (global_prm->show_optimization)
        //{
            std::cout << "First Background Fit:\n" << global_ctfmodel << std::endl;
            saveIntermediateResults_fast("step01c_first_background_fit_fast");
        //}

        center_optimization_focus_fast(false, true, 1.5);

    }
}

// Estimate envelope parameters --------------------------------------------
//#define DEBUG
void estimate_envelope_parameters_fast()
{
    if (global_prm->show_optimization)
        std::cout << "Looking for best fitting envelope ...\n";

    // Set the envelope
    global_ctfmodel.Ca = global_prm->initial_ctfmodel.Ca;
    global_ctfmodel.K = 1.0;
    global_ctfmodel.espr = 0.0;
    global_ctfmodel.ispr = 0.0;
    global_ctfmodel.alpha = 0.0;
    global_ctfmodel.DeltaF = 0.0;
    global_ctfmodel.DeltaR = 0.0;
    global_ctfmodel.Q0 = global_prm->initial_ctfmodel.Q0;
    global_ctfmodel.envR0 = 0.0;
    global_ctfmodel.envR1 = 0.0;
    COPY_ctfmodel_TO_CURRENT_GUESS;

    // Now optimize the envelope
    global_penalize = false;
    int iter;
    double fitness;
    Matrix1D<double> steps;
    steps.resize(ENVELOPE_PARAMETERS);
    steps.initConstant(1);
    //steps correctos con mis parametros
    steps(1) = 0; // Do not optimize Cs
    steps(5) = 0; // Do not optimize for alpha, since Ealpha depends on the
    // defocus
    if (global_prm->modelSimplification >= 1)
        steps(6) = steps(7) = 0; // Do not optimize DeltaF and DeltaR
    powellOptimizer(*global_adjust, FIRST_ENVELOPE_PARAMETER + 1,
                    ENVELOPE_PARAMETERS, &CTF_fitness_fast, NULL, 0.05, fitness, iter, steps,
                    global_prm->show_optimization);

    // Keep the result in global_prm->adjust
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    //if (global_prm->show_optimization)
    //{
        std::cout << "Best envelope Fit:\n" << global_ctfmodel << std::endl;
        saveIntermediateResults_fast("step02a_best_envelope_fit_fast");
    //}

    // Optimize with penalization
    if (global_prm->show_optimization)
        std::cout << "Penalizing best fitting envelope ...\n";
    global_penalize = true;
    global_current_penalty = 2;
    int imax = CEIL(log(global_penalty) / log(2.0));
    for (int i = 1; i <= imax; i++)
    {
        if (global_prm->show_optimization)
            std::cout << "     Iteration " << i << " penalty="
            << global_current_penalty << std::endl;
        powellOptimizer(*global_adjust, FIRST_ENVELOPE_PARAMETER + 1,
                        ENVELOPE_PARAMETERS, &CTF_fitness_fast, NULL, 0.05, fitness, iter,
                        steps, global_prm->show_optimization);
        global_current_penalty *= 2;
        global_current_penalty =
            XMIPP_MIN(global_current_penalty, global_penalty);
    }
    // Keep the result in global_prm->adjust
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    //if (global_prm->show_optimization)
    //{
        std::cout << "Best envelope Fit:\n" << global_ctfmodel << std::endl;
        saveIntermediateResults_fast("step02b_best_penalized_envelope_fit_fast");
    //}
}

// Estimate defoci ---------------------------------------------------------
void estimate_defoci_fast()
{
	double fitness;
    int iter;
	Matrix1D<double> steps(DEFOCUS_PARAMETERS);
	steps.initConstant(1);
	steps(1) = 0; // Do not optimize kV
	steps(2) = 0; // Do not optimize K
	(*global_adjust)(0) = global_prm->initial_ctfmodel.Defocus;
	(*global_adjust)(2) = global_ctfmodel.K;
	powellOptimizer(*global_adjust, FIRST_DEFOCUS_PARAMETER + 1,
								DEFOCUS_PARAMETERS, &CTF_fitness_fast, NULL, 0.05,
								fitness, iter, steps, false);

	// Keep the result in global_prm->adjust
	global_ctfmodel.forcePhysicalMeaning();
	COPY_ctfmodel_TO_CURRENT_GUESS;
	global_ctfmodel_defoci = global_ctfmodel;

	std::cout << "First defocus Fit:\n" << global_ctfmodel_defoci << std::endl;
	saveIntermediateResults_fast("step03a_first_defocus_fit_fast");

}

// Estimate second gaussian parameters -------------------------------------
//#define DEBUG
void estimate_background_gauss_parameters2_fast()
{
    if (global_prm->show_optimization)
        std::cout << "Computing first background Gaussian2 parameters ...\n";

    // Compute radial averages
    MultidimArray<double> radial_CTFmodel_avg(YSIZE(*f) / 2);
    MultidimArray<double> radial_CTFampl_avg(YSIZE(*f) / 2);
    MultidimArray<int> radial_N(YSIZE(*f) / 2);
    double w_max_gauss = 0.25;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)
    {
        if (global_mask(i) <= 0)
            continue;
        double w = global_w_digfreq(i);
        if (w > w_max_gauss)
            continue;

        int r = FLOOR(w * (double)YSIZE(*f));
        double f_x = DIRECT_A1D_ELEM(global_x_contfreq, i);
        global_ctfmodel.precomputeValues(f_x);
        double bg = global_ctfmodel.getValueNoiseAt();
        double envelope = global_ctfmodel.getValueDampingAt();
        double ctf_without_damping =
            global_ctfmodel.getValuePureWithoutDampingAt();
        double ctf_with_damping = envelope * ctf_without_damping;
        double ctf2_th = bg + ctf_with_damping * ctf_with_damping;
        radial_CTFmodel_avg(r) += ctf2_th;
        radial_CTFampl_avg(r) += global_psd_exp_radial(i);
        radial_N(r)++;
    }

    // Compute the average radial error
    MultidimArray <double> error;
    error.initZeros(radial_CTFmodel_avg);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        error(i) = (radial_CTFampl_avg(i) - radial_CTFmodel_avg(i))/ radial_N(i);
    }
#ifdef DEBUG
    std::cout << "Error:\n" << error << std::endl;
#endif

    // Compute the frequency of the minimum error
    double wmin = 0.15;

    // Compute the maximum (negative) radial error
    double error_max = 0, wmax=0, fmax;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        double w = global_w_digfreq(i);
        if (w > wmin)
            break;
        if (error(i) < error_max)
        {
            wmax = w;
            error_max = error(i);
        }
    }
    fmax = global_ctfmodel.Gc2 = wmax / global_prm->Tm;
#ifdef DEBUG

    std::cout << "Freq of the maximum error: " << wmax << " " << fmax << std::endl;
#endif

    // Find the linear least squares solution for the gauss part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();
    int N = 0;
    /*MultidimArray<double> w_contfreq;
    MultidimArray<double> x_contfreq;
    w_contfreq.initZeros(YSIZE(global_w_contfreq));
    x_contfreq.initZeros(YSIZE(global_x_contfreq));
    int cont=0;*/

    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq)
    {
        if (global_mask(i) <= 0)
            continue;
        if (global_w_digfreq(i) > wmin)
            continue;
        double fmod = global_w_contfreq(i);
        printf("fmod = %lf\n",fmod);

        // Compute the zero on the direction of this point
        Matrix1D<double> u(1), fzero(1);
        u = global_x_contfreq(i) / fmod;
        printf("XX(u) = %lf u = %lf\n",XX(u),u);
        //YY(u) = global_y_contfreq(i) / fmod;
        global_ctfmodel.lookFor(1, u, fzero, 0);
        printf("fzero = %lf\n",fzero.module());
        if (fmod > fzero.module())
            continue;

        // Compute weight for this point
        double weight = 1 + global_max_freq - global_w_digfreq(i);
        printf("w = %lf\n",weight);

        // Compute error
        double f_x = DIRECT_A1D_ELEM(global_x_contfreq, i);
        global_ctfmodel.precomputeValues(f_x);
        double bg = global_ctfmodel.getValueNoiseAt();
        double envelope = global_ctfmodel.getValueDampingAt();
        double ctf_without_damping =
            global_ctfmodel.getValuePureWithoutDampingAt();
        double ctf_with_damping = envelope * ctf_without_damping;
        double ctf2_th = bg + ctf_with_damping * ctf_with_damping;
        double explained = ctf2_th;
        double unexplained = explained - global_psd_exp_radial(i);
        printf("uneplained = %lf\n",unexplained);
        if (unexplained <= 0)
            continue;
        unexplained = log(unexplained);
        double F = -(fmod - fmax) * (fmod - fmax);
        A(0, 0) += weight * 1;
        A(0, 1) += weight * F;
        A(1, 1) += weight * F * F;
        b(0) += weight * unexplained;
        b(1) += F * weight * unexplained;
        N++;
    }

    printf("A=\n %lf %lf\n %lf %lf\n",A(0,0),A(0,1),A(1,0),A(1,1));
    printf("b(0) = %lf\n, b(1) = %lf \n",b(0),b(1));

    if (N != 0)
    {
        A(1, 0) = A(0, 1);

        double det=A.det();
        if (fabs(det)>1e-9)
        {
            b = A.inv() * b;
            global_ctfmodel.sigma2 = XMIPP_MIN(fabs(b(1)), 95e3); // This value should be conformant with the physical
            // meaning routine in CTF.cc
            global_ctfmodel.gaussian_K2 = exp(b(0));
        }
        else
        {
            global_ctfmodel.sigma2 = 0;
            global_ctfmodel.gaussian_K2 = 0;
        }
    }
    else
    {
        global_ctfmodel.sigma2 = 0;
        global_ctfmodel.gaussian_K2 = 0;
    }

    // Store the CTF values in global_prm->adjust
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (global_prm->show_optimization)
    {
        std::cout << "First Background Gaussian 2 Fit:\n" << global_ctfmodel
        << std::endl;
        saveIntermediateResults_fast("step04a_first_background2_fit_fast");
    }
}
#undef DEBUG


/* Main routine ------------------------------------------------------------ */
//#define DEBUG
double ROUT_Adjust_CTFFast(ProgCTFEstimateFromPSDFast &prm, CTF1D &output_ctfmodel, bool standalone)
{
	DEBUG_OPEN_TEXTFILE(prm.fn_psd.removeLastExtension());
	global_prm = &prm;
	if (standalone || prm.show_optimization)
	   prm.show();
	prm.produceSideInfo();
	DEBUG_TEXTFILE(formatString("After producing side info: Avg=%f",prm.ctftomodel().computeAvg()));
	DEBUG_MODEL_TEXTFILE;
	// Build initial frequency mask
	global_value_th = -1;
	global_min_freq = prm.min_freq;
	global_max_freq = prm.max_freq;

	// Set some global variables
	global_adjust = &prm.adjust;
	global_penalize = false;
	global_max_gauss_freq = 0;
	global_heavy_penalization = f->computeMax() * YSIZE(*f); //xsize = 512
	global_show = 0;

	// Some variables needed by all steps
	int iter;
	double fitness;
	Matrix1D<double> steps;
	/************************************************************************/
	/* STEPs 1, 2, 3 and 4:  Find background which best fits the CTF        */
	/************************************************************************/

	global_ctfmodel.enable_CTFnoise = true;
	global_ctfmodel.enable_CTF = false;
	global_evaluation_reduction = 4;

	// If initial parameters were not supplied for the gaussian curve,
	// estimate them from the CTF file
	global_action = 0;
	if (prm.adjust(FIRST_SQRT_PARAMETER) == 0)
	{
		estimate_background_sqrt_parameters_fast();
		estimate_background_gauss_parameters_fast();
	}

	std::cout << "Best background Fit:\n" << global_ctfmodel << std::endl;
	saveIntermediateResults_fast("step01c_best_background_fit_fast");

	 // Optimize the current background
	global_action = 1;
	global_penalize = true;
	global_current_penalty = global_penalty;
	steps.resize(BACKGROUND_CTF_PARAMETERS);
	steps.initConstant(1);

	if (!global_prm->modelSimplification >= 3)
		steps(1) = steps(2) = steps(4) = 0;
	powellOptimizer(*global_adjust, FIRST_SQRT_PARAMETER + 1,
					BACKGROUND_CTF_PARAMETERS, &CTF_fitness_fast, NULL, 0.01, fitness, iter,
					steps, global_prm->show_optimization);

	// Make sure that the model has physical meaning
	// (In some machines due to numerical imprecission this check is necessary
	// at the end)
	global_ctfmodel.forcePhysicalMeaning();
	COPY_ctfmodel_TO_CURRENT_GUESS;

	//if (global_prm->show_optimization)
	//{
		std::cout << "Best background Fit:\n" << global_ctfmodel << std::endl;
		saveIntermediateResults_fast("step01d_best_background_fit_fast");
	//}
	DEBUG_TEXTFILE(formatString("Step 4: CTF_fitness=%f",CTF_fitness_fast));
	DEBUG_MODEL_TEXTFILE;

	/************************************************************************/
	/* STEPs 5 and 6:  Find envelope which best fits the CTF                */
	/************************************************************************/
	//printf("modelSimplification2 = %i \n", prm.modelSimplification);
	global_action = 2;
	global_ctfmodel.enable_CTF = true;

	global_ctfmodel.kV = prm.initial_ctfmodel.kV;
	global_ctfmodel.Cs = prm.initial_ctfmodel.Cs;
	if (prm.initial_ctfmodel.Q0 != 0)
		global_ctfmodel.Q0 = prm.initial_ctfmodel.Q0;
	estimate_envelope_parameters_fast();

	std::cout << "Best envelope Fit:\n" << global_ctfmodel << std::endl;
	saveIntermediateResults_fast("step02b_best_penalized_envelope_fit_fast");
	DEBUG_TEXTFILE(formatString("Step 6: espr=%f",global_ctfmodel.espr));
	DEBUG_MODEL_TEXTFILE;

	/************************************************************************/
	/* STEP 7:  the defocus and angular parameters                          */
	/************************************************************************/
	global_action = 3;
	estimate_defoci_fast();

	DEBUG_TEXTFILE(formatString("Step 7: Defocus=%f",global_ctfmodel.Defocus));
	DEBUG_MODEL_TEXTFILE;
	/************************************************************************/
	/* STEPs 9, 10 and 11: all parameters included second Gaussian          */
	/************************************************************************/
	global_action = 5;
	if (prm.modelSimplification < 2)
		estimate_background_gauss_parameters2_fast();

	steps.resize(ALL_CTF_PARAMETERS);
	steps.initConstant(1);
	steps(1) = 0; // kV
	steps(3) = 0; // The spherical aberration (Cs) is not optimized
	if (prm.initial_ctfmodel.Q0 != 0)
		steps(13) = 0; // Q0
	if (prm.modelSimplification >= 1)
		steps(8) = steps(9) = 0;

	powellOptimizer(*global_adjust, 0 + 1, ALL_CTF_PARAMETERS, &CTF_fitness_fast,
					NULL, 0.01, fitness, iter, steps, global_prm->show_optimization);

	global_ctfmodel.forcePhysicalMeaning();
	COPY_ctfmodel_TO_CURRENT_GUESS;

	if (global_prm->show_optimization)
	{
		std::cout << "Best fit with Gaussian2:\n" << global_ctfmodel
		<< std::endl;
		saveIntermediateResults_fast("step04b_best_fit_with_gaussian2_fast");
	}
	return fitness;
}

void ProgCTFEstimateFromPSDFast::run()
{
	 	CTF1D ctf1Dmodel;
	    ROUT_Adjust_CTFFast(*this, ctf1Dmodel);
}
