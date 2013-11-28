/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.csic.es)
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

#include "ctf_estimate_psd_with_arma.h"
#include "ctf_estimate_from_psd.h"
#include "ctf_enhance_psd.h"
#include "fourier_filter.h"
#include "fringe_processing.h"

#include <data/args.h>
#include <data/histogram.h>
#include <data/filters.h>
#include <data/xmipp_fft.h>

/* prototypes */
double CTF_fitness(double *, void *);

/* Number of CTF parameters */
#define ALL_CTF_PARAMETERS         30
#define CTF_PARAMETERS             24
#define PARAMETRIC_CTF_PARAMETERS  13
#define BACKGROUND_CTF_PARAMETERS  11
#define SQRT_CTF_PARAMETERS         5
#define ENVELOPE_PARAMETERS         8
#define DEFOCUS_PARAMETERS          5
#define FIRST_SQRT_PARAMETER       13
#define FIRST_ENVELOPE_PARAMETER    4
#define FIRST_DEFOCUS_PARAMETER     0

//#define DEBUG_WITH_TEXTFILES
#ifdef DEBUG_WITH_TEXTFILES
std::ofstream fhDebug;
#define DEBUG_OPEN_TEXTFILE(fnRoot) fhDebug.open((fnRoot+"_debug.txt").c_str());
#define DEBUG_CLOSE_TEXTFILE fhDebug.close();
#define DEBUG_TEXTFILE(str) fhDebug << time (NULL) << " " << str << std::endl;
#define DEBUG_MODEL_TEXTFILE fhDebug << global_ctfmodel << std::endl;
#else
#define DEBUG_OPEN_TEXTFILE(fnRoot);
#define DEBUG_CLOSE_TEXTFILE ;
#define DEBUG_TEXTFILE(str);
#define DEBUG_MODEL_TEXTFILE
#endif

/* Global variables -------------------------------------------------------- */
namespace AdjustCTF
{
// Some aliases
ProgCTFEstimateFromPSD *global_prm;
MultidimArray<double> *f; // The CTF to model
Matrix1D<double> *global_adjust; // Current theoretical adjustment
double global_corr13; // Correlation with enhanced PSD between rings 1 and 3

// Frequency of each point in digital units
MultidimArray<double> global_x_digfreq;
MultidimArray<double> global_y_digfreq;
MultidimArray<double> global_w_digfreq;
MultidimArray<int> global_w_digfreq_r;
MultidimArray<double> global_w_digfreq_r_iN;
MultidimArray<double> global_x_contfreq;
MultidimArray<double> global_y_contfreq;
MultidimArray<double> global_w_contfreq;
MultidimArray<double> global_mask;
MultidimArray<double> global_mask_between_zeroes;
MultidimArray<double> global_w_count;

MultidimArray<double> global_psd_exp_radial_derivative;
MultidimArray<double> global_psd_theo_radial_derivative;
MultidimArray<double> global_psd_exp_radial;
MultidimArray<double> global_psd_theo_radial;

// Penalization for forbidden values of the parameters
double global_heavy_penalization;

// Penalization factor for the background
bool global_penalize;
double global_current_penalty;
const double global_penalty = 32; // Maximum penalization

// Speed up factor
int global_evaluation_reduction;

// CTF model and noise model
CTFDescription global_ctfmodel;
CTFDescription global_ctfmodel_defoci;

// Maximum of the gaussian
double global_max_gauss_freq;

// Autofocus
double global_value_th;
double global_min_freq;
double global_max_freq;

// Program status
int global_action; // 0: Computing the background (sqrt)
// 1: Computing the full background
// 2: Computing the envelope
// 3: Computing defoci
// 4: Computing all CTF parameters
// 5: Computing all CTF parameters + Gaussian2
// 6: Produce output
int global_show; // 0: Do not show
// 1: Partially detailed
// 2: Very detailed
}

using namespace AdjustCTF;

#define ASSIGN_CTF_PARAM(index, paramName) if (ia <= index && l > 0) { ctfmodel.paramName = p[index]; --l; }

/* Assign ctfmodel from a vector and viceversa ----------------------------- */
void assignCTFfromParameters(double *p, CTFDescription &ctfmodel, int ia,
                             int l, int modelSimplification)
{
    ctfmodel.Tm = global_prm->Tm;

    ASSIGN_CTF_PARAM(0, DeltafU);
    ASSIGN_CTF_PARAM(1, DeltafV);
    ASSIGN_CTF_PARAM(2, azimuthal_angle);
    ASSIGN_CTF_PARAM(3, kV);
    ASSIGN_CTF_PARAM(4, K);
    ASSIGN_CTF_PARAM(5, Cs);
    ASSIGN_CTF_PARAM(6, Ca);
    ASSIGN_CTF_PARAM(7, espr);
    ASSIGN_CTF_PARAM(8, ispr);
    ASSIGN_CTF_PARAM(9, alpha);
    ASSIGN_CTF_PARAM(10, DeltaF);
    ASSIGN_CTF_PARAM(11, DeltaR);
    ASSIGN_CTF_PARAM(12, Q0);
    ASSIGN_CTF_PARAM(13, base_line);
    ASSIGN_CTF_PARAM(14, sqrt_K);
    ASSIGN_CTF_PARAM(15, sqU);
    ASSIGN_CTF_PARAM(16, sqV);
    ASSIGN_CTF_PARAM(17, sqrt_angle);
    ASSIGN_CTF_PARAM(18, gaussian_K);
    ASSIGN_CTF_PARAM(19, sigmaU);

    if (ia <= 20 && l > 0)
    {
        if (modelSimplification < 3)
        {
            ctfmodel.sigmaV = p[20];
            l--;
        } //     7 *
        else
        {
            ctfmodel.sigmaV = p[19];
            l--;
        }
    }
    if (ia <= 21 && l > 0)
    {
        if (modelSimplification < 3)
        {
            ctfmodel.gaussian_angle = p[21];
            l--;
        } //     8 *
        else
        {
            ctfmodel.gaussian_angle = 0;
            l--;
        }
    }

    ASSIGN_CTF_PARAM(22, cU);

    if (ia <= 23 && l > 0)
    {
        if (modelSimplification < 3)
        {
            ctfmodel.cV = p[23];
            l--;
        } //    10 *
        else
        {
            ctfmodel.cV = p[22];
            l--;
        }
    }

    ASSIGN_CTF_PARAM(24, gaussian_K2);
    ASSIGN_CTF_PARAM(25, sigmaU2);
    ASSIGN_CTF_PARAM(26, sigmaV2);
    ASSIGN_CTF_PARAM(27, gaussian_angle2);
    ASSIGN_CTF_PARAM(28, cU2);
    ASSIGN_CTF_PARAM(29, cV2);
}//function assignCTFfromParameters


#define ASSIGN_PARAM_CTF(index, paramName) if (ia <= index && l > 0) { p[index] = ctfmodel.paramName; --l; }

void assignParametersFromCTF(CTFDescription &ctfmodel, double *p, int ia,
                             int l, int modelSimplification)
{
    ASSIGN_PARAM_CTF(0, DeltafU);
    ASSIGN_PARAM_CTF(1, DeltafV);
    ASSIGN_PARAM_CTF(2, azimuthal_angle);
    ASSIGN_PARAM_CTF(3, kV);
    ASSIGN_PARAM_CTF(4, K);
    ASSIGN_PARAM_CTF(5, Cs);
    ASSIGN_PARAM_CTF(6, Ca);
    ASSIGN_PARAM_CTF(7, espr);
    ASSIGN_PARAM_CTF(8, ispr);
    ASSIGN_PARAM_CTF(9, alpha);
    ASSIGN_PARAM_CTF(10, DeltaF);
    ASSIGN_PARAM_CTF(11, DeltaR);
    ASSIGN_PARAM_CTF(12, Q0);
    ASSIGN_PARAM_CTF(13, base_line);
    ASSIGN_PARAM_CTF(14, sqrt_K);
    ASSIGN_PARAM_CTF(15, sqU);
    ASSIGN_PARAM_CTF(16, sqV);
    ASSIGN_PARAM_CTF(17, sqrt_angle);
    ASSIGN_PARAM_CTF(18, gaussian_K);
    ASSIGN_PARAM_CTF(19, sigmaU);

    if (ia <= 20 && l > 0)
    {
        if (modelSimplification < 3)
        {
            p[20] = ctfmodel.sigmaV;
            l--;
        }
        else
        {
            p[20] = 0;
            l--;
        }
    }
    if (ia <= 21 && l > 0)
    {
        if (modelSimplification < 3)
        {
            p[21] = ctfmodel.gaussian_angle;
            l--;
        }
        else
        {
            p[21] = 0;
            l--;
        }
    }

    ASSIGN_PARAM_CTF(22, cU);

    if (ia <= 23 && l > 0)
    {
        if (modelSimplification < 3)
        {
            p[23] = ctfmodel.cV;
            l--;
        }
        else
        {
            p[23] = 0;
            l--;
        }
    }

    ASSIGN_PARAM_CTF(24, gaussian_K2);
    ASSIGN_PARAM_CTF(25, sigmaU2);
    ASSIGN_PARAM_CTF(26, sigmaV2);
    ASSIGN_PARAM_CTF(27, gaussian_angle2);
    ASSIGN_PARAM_CTF(28, cU2);
    ASSIGN_PARAM_CTF(29, cV2);
}

#define COPY_ctfmodel_TO_CURRENT_GUESS \
    assignParametersFromCTF(global_ctfmodel, \
                               MATRIX1D_ARRAY(*global_adjust),0,ALL_CTF_PARAMETERS, \
                               global_prm->modelSimplification);

/* Read parameters --------------------------------------------------------- */
void ProgCTFEstimateFromPSD::readBasicParams(XmippProgram *program)
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

    initial_ctfmodel.enable_CTF = initial_ctfmodel.enable_CTFnoise = true;
    initial_ctfmodel.readParams(program);
    if (initial_ctfmodel.DeltafU>100e3 || initial_ctfmodel.DeltafV>100e3)
    	REPORT_ERROR(ERR_ARG_INCORRECT,"Defocus cannot be larger than 10 microns (100,000 Angstroms)");
    Tm = initial_ctfmodel.Tm;
}

void ProgCTFEstimateFromPSD::readParams()
{
    fn_psd = getParam("--psd");
    readBasicParams(this);
}

/* Show -------------------------------------------------------------------- */
void ProgCTFEstimateFromPSD::show()
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
void ProgCTFEstimateFromPSD::defineBasicParams(XmippProgram * program)
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
    CTFDescription::defineParams(program);
}

void ProgCTFEstimateFromPSD::defineParams()
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
void ProgCTFEstimateFromPSD::produceSideInfo()
{
    adjust.resize(ALL_CTF_PARAMETERS);
    adjust.initZeros();
    global_ctfmodel.clear();
    global_ctfmodel_defoci.clear();
    assignParametersFromCTF(initial_ctfmodel, MATRIX1D_ARRAY(adjust), 0,
                            ALL_CTF_PARAMETERS, true);

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
    global_ctfmodel.precomputeValues(global_x_contfreq, global_y_contfreq);

    // Build frequency mask
    global_mask.initZeros(global_w_digfreq);
    global_w_count.initZeros(XSIZE(global_w_digfreq));
    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_w_digfreq(i, j) >= max_freq
            || global_w_digfreq(i, j) <= min_freq)
            continue;
        global_mask(i, j) = 1;
        global_w_count(global_w_digfreq_r(i, j))++;
    }

    // Bootstrap
    if (bootstrap)
    {
        double N=global_mask.sum();
        std::vector< Matrix1D<int> > positions;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(global_mask)
        if (global_mask(i,j))
        {
            Matrix1D<int> r(2);
            XX(r)=j;
            YY(r)=i;
            positions.push_back(r);
        }
        global_mask.initZeros();
        for (int n=0; n<N; n++)
        {
            int idx=ROUND(rnd_unif(0,N-1));
            global_mask(YY(positions[idx]),XX(positions[idx]))+=1;
        }
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
    global_psd_exp_radial_derivative.initZeros(XSIZE(enhanced_ctftomodel()));
    global_psd_theo_radial_derivative.initZeros(global_psd_exp_radial_derivative);
    global_psd_exp_radial.initZeros(global_psd_exp_radial_derivative);
    global_psd_theo_radial.initZeros(global_psd_exp_radial_derivative);
    global_w_digfreq_r_iN.initZeros(global_psd_exp_radial_derivative);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(enhanced_ctftomodel())
    {
        if (global_w_digfreq(i,j)>min_freq && global_w_digfreq(i,j)<max_freq)
        {
            int r = global_w_digfreq_r(i, j);
            global_w_digfreq_r_iN(r)+=1;
            global_psd_exp_radial(r) += enhanced_ctftomodel(i, j);
        }
    }
    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_w_digfreq_r_iN)
    if (global_w_digfreq_r_iN(i)>0)
    {
        global_w_digfreq_r_iN(i)=1.0/global_w_digfreq_r_iN(i);
        global_psd_exp_radial(i)*=global_w_digfreq_r_iN(i);
    }

    // Compute its derivative
    int state=0;
    double maxDiff=0;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(global_psd_exp_radial)
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
                double diff=global_psd_exp_radial(i)-global_psd_exp_radial(i-1);
                global_psd_exp_radial_derivative(i)=diff;
                maxDiff=std::max(maxDiff,fabs(diff));
            }
            break;
        }
    }
    global_psd_exp_radial_derivative/=maxDiff;
}

/* Generate model so far ---------------------------------------------------- */
/* The model is taken from global_adjust and global_ctfmodel is modified */
void generateModelSoFar(Image<double> &I, bool apply_log = false)
{
    Matrix1D<int> idx(2); // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane

    assignCTFfromParameters(MATRIX1D_ARRAY(*global_adjust), global_ctfmodel,
                            0, ALL_CTF_PARAMETERS, global_prm->modelSimplification);
    global_ctfmodel.produceSideInfo();

    I().resize(*f);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I())
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(*f, idx, freq);
        digfreq2contfreq(freq, freq, global_prm->Tm);

        // Decide what to save
        global_ctfmodel.precomputeValues(XX(freq), YY(freq));
        if (global_action <= 1)
            I()(i, j) = global_ctfmodel.getValueNoiseAt();
        else if (global_action == 2)
        {
            double E = global_ctfmodel.getValueDampingAt();
            I()(i, j) = global_ctfmodel.getValueNoiseAt() + E * E;
        }
        else if (global_action >= 3 && global_action <= 5)
        {
            double ctf = global_ctfmodel.getValuePureAt();
            I()(i, j) = global_ctfmodel.getValueNoiseAt() + ctf * ctf;
        }
        else
        {
            double ctf = global_ctfmodel.getValuePureAt();
            I()(i, j) = ctf;
        }
        if (apply_log)
            I()(i, j) = 10 * log10(I()(i, j));
    }
}

/* Save intermediate results ----------------------------------------------- */
/* First call to generate model so far and then save the image, and a couple
 of cuts along X and Y.

 This function returns the fitting error.*/
void saveIntermediateResults(const FileName &fn_root, bool generate_profiles =
                                 true)
{
    std::ofstream plotX, plotY, plot_radial;
    Image<double> save;
    generateModelSoFar(save, false);

    Image<double> save_ctf;
    global_prm->generate_model_halfplane(global_prm->ctfmodelSize,
                                         global_prm->ctfmodelSize, save_ctf());
    if (fn_root.find("@")==std::string::npos)
        save_ctf.write(fn_root + "_ctfmodel_halfplane.xmp");
    else
        save_ctf.write(fn_root + "_ctfmodel_halfplane.stk");
    global_prm->generate_model_quadrant(global_prm->ctfmodelSize,
                                        global_prm->ctfmodelSize, save_ctf());
    if (fn_root.find("@")==std::string::npos)
        save_ctf.write(fn_root + "_ctfmodel_quadrant.xmp");
    else
        save_ctf.write(fn_root + "_ctfmodel_quadrant.stk");

    if (!generate_profiles)
        return;
    plotX.open((fn_root + "X.txt").c_str());
    plotY.open((fn_root + "Y.txt").c_str());
    plot_radial.open((fn_root + "_radial.txt").c_str());
    if (!plotX)
        REPORT_ERROR(
            ERR_IO_NOWRITE,
            "save_intermediate_results::Cannot open plot file for writing\n");
    if (!plotY)
        REPORT_ERROR(
            ERR_IO_NOWRITE,
            "save_intermediate_results::Cannot open plot file for writing\n");
    if (!plot_radial)
        REPORT_ERROR(
            ERR_IO_NOWRITE,
            "save_intermediate_results::Cannot open plot file for writing\n");
    plotX << "# freq_dig freq_angstrom model psd enhanced\n";
    plotY << "# freq_dig freq_angstrom model psd enhanced\n";
    plot_radial << "# freq_dig freq_angstrom model psd enhanced\n";

    // Generate cut along X
    for (int i = STARTINGY(save()); i <= FINISHINGY(save()) / 2; i++)
    {
        if (global_mask(i, 0) <= 0)
            continue;
        plotY << global_w_digfreq(i, 0) << " " << global_w_contfreq(i, 0) << " "
        << save()(i, 0) << " " << (*f)(i, 0) << " "
        << global_prm->enhanced_ctftomodel()(i, 0) << std::endl;
    }

    // Generate cut along Y
    for (int j = STARTINGX(save()); j <= FINISHINGX(save()) / 2; j++)
    {
        if (global_mask(0, j) <= 0)
            continue;
        plotX << global_w_digfreq(0, j) << " " << global_w_contfreq(0, j) << " "
        << save()(0, j) << " " << (*f)(0, j) << " "
        << global_prm->enhanced_ctftomodel()(0, j) << std::endl;
    }

    // Generate radial average
    MultidimArray<double> radial_CTFmodel_avg, radial_CTFampl_avg,
    radial_enhanced_avg;
    MultidimArray<int> radial_N;
    radial_CTFmodel_avg.initZeros(YSIZE(save()) / 2);
    radial_CTFampl_avg.initZeros(YSIZE(save()) / 2);
    radial_enhanced_avg.initZeros(YSIZE(save()) / 2);
    radial_N.initZeros(YSIZE(save()) / 2);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_mask(i, j) <= 0)
            continue;
        double model2 = save()(i, j);

        int r = global_w_digfreq_r(i, j);
        radial_CTFmodel_avg(r) += model2;
        radial_CTFampl_avg(r) += (*f)(i, j);
        radial_enhanced_avg(r) += global_prm->enhanced_ctftomodel()(i, j);
        radial_N(r)++;
    }

    FOR_ALL_ELEMENTS_IN_ARRAY1D (radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        plot_radial << global_w_digfreq(i, 0) << " " << global_w_contfreq(i, 0)
        << " " << radial_CTFmodel_avg(i) / radial_N(i) << " "
        << radial_CTFampl_avg(i) / radial_N(i) << " "
        << radial_enhanced_avg(i) / radial_N(i) << std::endl;
    }

    plotX.close();
    plotY.close();
    plot_radial.close();
}

/* Generate model at a given size ------------------------------------------ */
void ProgCTFEstimateFromPSD::generate_model_quadrant(int Ydim, int Xdim,
        MultidimArray<double> &model)
{
    Matrix1D<int> idx(2); // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane

    // Compute the scaled PSD
    MultidimArray<double> enhancedPSD;
    enhancedPSD = global_prm->enhanced_ctftomodel_fullsize();
    CenterFFT(enhancedPSD, false);
    selfScaleToSize(BSPLINE3, enhancedPSD, Ydim, Xdim);
    CenterFFT(enhancedPSD, true);

    // Generate the CTF model
    assignCTFfromParameters(MATRIX1D_ARRAY(*global_adjust), global_ctfmodel,
                            0, ALL_CTF_PARAMETERS, global_prm->modelSimplification);
    global_ctfmodel.produceSideInfo();

    // Write the two model quadrants
    MultidimArray<int> mask;
    mask.initZeros(enhancedPSD);
    model.resizeNoCopy(enhancedPSD);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(model)
    {
        if ((j >= Xdim / 2 && i >= Ydim / 2)
            || (j < Xdim / 2 && i < Ydim / 2))
        {
            XX(idx) = j;
            YY(idx) = i;
            FFT_idx2digfreq(model, idx, freq);
            if (fabs(XX(freq))>0.03 && fabs(YY(freq))>0.03)
                mask(i,j)=(int)global_mask(i,j);
            digfreq2contfreq(freq, freq, global_prm->Tm);

            global_ctfmodel.precomputeValues(XX(freq), YY(freq));
            model(i, j) = global_ctfmodel.getValuePureAt();
            model(i, j) *= model(i, j);
        }
    }

    // Normalize the left part so that it has similar values to
    // the enhanced PSD
    model.rangeAdjust(enhancedPSD, &mask);

    // Copy the part of the enhancedPSD
    FOR_ALL_ELEMENTS_IN_ARRAY2D(model)
    if (!((j >= Xdim / 2 && i >= Ydim / 2) || (j < Xdim / 2 && i < Ydim / 2)))
        model(i, j) = enhancedPSD(i, j);

    // Produce a centered image
    CenterFFT(model, true);
}

void ProgCTFEstimateFromPSD::generate_model_halfplane(int Ydim, int Xdim,
        MultidimArray<double> &model)
{
    Matrix1D<int> idx(2); // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane

    // Compute the scaled PSD
    MultidimArray<double> enhancedPSD;
    enhancedPSD = global_prm->enhanced_ctftomodel_fullsize();
    CenterFFT(enhancedPSD, false);
    selfScaleToSize(BSPLINE3, enhancedPSD, Ydim, Xdim);
    CenterFFT(enhancedPSD, true);

    // The left part is the CTF model
    assignCTFfromParameters(MATRIX1D_ARRAY(*global_adjust), global_ctfmodel,
                            0, CTF_PARAMETERS, global_prm->modelSimplification);
    global_ctfmodel.produceSideInfo();

    MultidimArray<int> mask;
    mask.initZeros(enhancedPSD);
    model.resizeNoCopy(enhancedPSD);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(model)
    {
        if (j >= Xdim / 2)
            continue;

        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(model, idx, freq);
        if (fabs(XX(freq))>0.03 && fabs(YY(freq))>0.03)
            mask(i,j)=(int)global_mask(i,j);
        digfreq2contfreq(freq, freq, global_prm->Tm);

        global_ctfmodel.precomputeValues(XX(freq), YY(freq));
        model(i, j) = global_ctfmodel.getValuePureAt();
        model(i, j) *= model(i, j);
    }

    // Normalize the left part so that it has similar values to
    // the enhanced PSD
    model.rangeAdjust(enhancedPSD, &mask);

    // Copy the part of the enhancedPSD
    FOR_ALL_ELEMENTS_IN_ARRAY2D(model)
    if (j >= Xdim / 2)
        model(i, j) = enhancedPSD(i, j);

    // Produce a centered image
    CenterFFT(model, true);
}

/* CTF fitness ------------------------------------------------------------- */
/* This function measures the distance between the estimated CTF and the
 measured CTF */
double CTF_fitness(double *p, void *)
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
        && (fabs(
                (global_ctfmodel.DeltafU - global_ctfmodel_defoci.DeltafU)
                / global_ctfmodel_defoci.DeltafU) > 0.2
            || fabs(
                (global_ctfmodel.DeltafV
                 - global_ctfmodel_defoci.DeltafV)
                / global_ctfmodel_defoci.DeltafU) > 0.2))
    {
        if (global_show >= 2)
            std::cout << "Too large defocus\n";
        return global_heavy_penalization;
    }
    if (global_prm->initial_ctfmodel.DeltafU != 0 && global_action >= 3)
    {
        // If there is an initial model, the true solution
        // cannot be too far
        if (fabs(global_prm->initial_ctfmodel.DeltafU - global_ctfmodel.DeltafU) > global_prm->defocus_range ||
            fabs(global_prm->initial_ctfmodel.DeltafV - global_ctfmodel.DeltafV) > global_prm->defocus_range)
        {
            if (global_show >= 2)
            {
                std::cout << "Too far from hint: Initial (" << global_prm->initial_ctfmodel.DeltafU << "," << global_prm->initial_ctfmodel.DeltafV << ")"
                << " current guess (" << global_ctfmodel.DeltafU << "," << global_ctfmodel.DeltafV << ") max allowed difference: "
                << global_prm->defocus_range << std::endl;
            }
            return global_heavy_penalization;
        }
    }

    // Now the 2D error
    double distsum = 0;
    int N = 0, Ncorr = 0;
    double enhanced_avg = 0;
    double model_avg = 0;
    double enhanced_model = 0;
    double enhanced2 = 0;
    double model2 = 0;
    double lowerLimit = 1.1 * global_min_freq;
    double upperLimit = 0.9 * global_max_freq;
    const MultidimArray<double>& local_enhanced_ctf =
        global_prm->enhanced_ctftomodel();
    global_psd_exp_radial.initZeros();
    int XdimW=XSIZE(global_w_digfreq);
    int YdimW=YSIZE(global_w_digfreq);
    for (int i = 0; i < YdimW; i +=
             global_evaluation_reduction)
        for (int j = 0; j < XdimW; j +=
                 global_evaluation_reduction)
        {
            if (DIRECT_A2D_ELEM(global_mask, i, j) <= 0)
                continue;

            // Compute each component
            global_ctfmodel.precomputeValues(i, j);
            double bg = global_ctfmodel.getValueNoiseAt();
            double envelope, ctf_without_damping, ctf_with_damping;
            double ctf2_th;
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
                envelope = global_ctfmodel.getValueDampingAt();
                ctf_without_damping =
                    global_ctfmodel.getValuePureWithoutDampingAt();
                ctf_with_damping = envelope * ctf_without_damping;
                ctf2_th = bg + ctf_with_damping * ctf_with_damping;
                break;
            }

            // Compute distance
            double ctf2 = DIRECT_A2D_ELEM(*f, i, j);
            double dist = 0;
            double ctf_with_damping2;
            switch (global_action)
            {
            case 0:
            case 1:
                dist = fabs(ctf2 - bg);
                if (global_penalize && bg > ctf2
                    && DIRECT_A2D_ELEM(global_w_digfreq, i, j)
                    > global_max_gauss_freq)
                    dist *= global_current_penalty;
                break;
            case 2:
                dist = fabs(ctf2 - ctf2_th);
                if (global_penalize && ctf2_th < ctf2
                    && DIRECT_A2D_ELEM(global_w_digfreq, i, j)
                    > global_max_gauss_freq)
                    dist *= global_current_penalty;
                break;
            case 4:
            case 5:
            case 6:
            case 3:
                if (DIRECT_A2D_ELEM(global_w_digfreq,i, j) < upperLimit
                    && DIRECT_A2D_ELEM(global_w_digfreq,i, j) > lowerLimit)
                {
                    if  (global_action == 3 ||
                         (global_action == 4 && DIRECT_A2D_ELEM(global_mask_between_zeroes,i,j) == 1) ||
                         (global_action == 6 && DIRECT_A2D_ELEM(global_mask_between_zeroes,i,j) == 1))
                    {
                        double enhanced_ctf =
                            DIRECT_A2D_ELEM(local_enhanced_ctf, i, j);
                        ctf_with_damping2 = ctf_with_damping * ctf_with_damping;
                        enhanced_model += enhanced_ctf * ctf_with_damping2;
                        enhanced2 += enhanced_ctf * enhanced_ctf;
                        model2 += ctf_with_damping2 * ctf_with_damping2;
                        enhanced_avg += enhanced_ctf;
                        model_avg += ctf_with_damping2;
                        Ncorr++;
                        if (global_action==3)
                        {
                            int r = A2D_ELEM(global_w_digfreq_r,i, j);
                            A1D_ELEM(global_psd_theo_radial,r) += ctf2_th;
                        }
                    }
                }
                if (envelope > 1e-2)
                    dist = fabs(ctf2 - ctf2_th) / (envelope * envelope);
                else
                    dist = fabs(ctf2 - ctf2_th);
                // This expression comes from mapping any value so that
                // bg becomes 0, and bg+envelope^2 becomes 1
                // This is the transformation
                //        (x-bg)      x-bg
                //    -------------=-------
                //    (bg+env^2-bg)  env^2
                // If we substract two of this scaled values
                //    x-bg      y-bg       x-y
                //   ------- - ------- = -------
                //    env^2     env^2     env^2
                break;
            }
            distsum += dist * DIRECT_A2D_ELEM(global_mask,i,j);
            N++;
        }
    if (N > 0)
        retval = distsum / N;
    else
        retval = global_heavy_penalization;
    if (global_show >=2)
        std::cout << "Fitness1=" << retval << std::endl;
    if ( (((global_action >= 3) && (global_action <= 4)) || (global_action == 6))
         && (Ncorr > 0) && (global_prm->enhanced_weight != 0) )
    {
        model_avg /= Ncorr;
        enhanced_avg /= Ncorr;
        double correlation_coeff = enhanced_model / Ncorr
                                   - model_avg * enhanced_avg;
        double sigma1 = sqrt(
                            fabs(enhanced2 / Ncorr - enhanced_avg * enhanced_avg));
        double sigma2 = sqrt(fabs(model2 / Ncorr - model_avg * model_avg));
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
            if (global_action == 6)
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
                A1D_ELEM(global_psd_theo_radial,i)*=A1D_ELEM(global_w_digfreq_r_iN,i);
                double freq=A2D_ELEM(global_w_digfreq,i,0);
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
                    }
                    break;
                }
            }
            double corrRadialDerivative=0,mux=0, muy=0, Ncorr=0, sigmax=0, sigmay=0;
            double iMaxDiff=1.0/maxDiff;
            FOR_ALL_ELEMENTS_IN_ARRAY1D(global_psd_theo_radial)
            {
                A1D_ELEM(global_psd_theo_radial_derivative,i)*=iMaxDiff;
                double x=A1D_ELEM(global_psd_exp_radial_derivative,i);
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
                    global_psd_exp_radial.write("PPPexpRadial.txt");
                    global_psd_theo_radial.write("PPPtheoRadial.txt");
                    global_psd_exp_radial_derivative.write("PPPexpRadialDerivative.txt");
                    global_psd_theo_radial_derivative.write("PPPtheoRadialDerivative.txt");
                }
            }
        }
    }

    // Show some debugging information
    if (global_show >= 2)
    {
        std::cout << "Fitness=" << retval << std::endl;
        if (global_show == 3)
        {
            saveIntermediateResults("PPP");
            std::cout << "Press any key\n";
            char c;
            std::cin >> c;
        }
    }

    return retval;
}

/* Compute central region -------------------------------------------------- */
void compute_central_region(double &w1, double &w2, double ang)
{
    w1 = global_max_freq;
    w2 = global_min_freq;
    Matrix1D<double> freq(2), dir(2);

    // Compute first and third zero in the given direction
    VECTOR_R2(dir, COSD(ang), SIND(ang));

    // Detect first zero
    global_ctfmodel.lookFor(1, dir, freq, 0);
    if (XX(freq) == -1 && YY(freq) == -1)
        w1 = global_min_freq;
    else
    {
        contfreq2digfreq(freq, freq, global_prm->Tm);
        double w;
        if (XX(dir) > 0.1)
            w = XX(freq) / XX(dir);
        else
            w = YY(freq) / YY(dir);
        w1 = XMIPP_MAX(global_min_freq, XMIPP_MIN(w1, w));
    }

    // Detect fifth zero
    global_ctfmodel.lookFor(5, dir, freq, 0);
    if (XX(freq) == -1 && YY(freq) == -1)
        w2 = global_max_freq;
    else
    {
        double w;
        contfreq2digfreq(freq, freq, global_prm->Tm);
        if (XX(dir) > 0.1)
            w = XX(freq) / XX(dir);
        else
            w = YY(freq) / YY(dir);
        w2 = XMIPP_MIN(global_max_freq, XMIPP_MAX(w2, w));
    }
}

/* Center focus ----------------------------------------------------------- */
void center_optimization_focus(bool adjust_freq, bool adjust_th, double margin =
                                   1)
{
    if (global_prm->show_optimization)
        std::cout << "Freq frame before focusing=" << global_min_freq << ","
        << global_max_freq << std::endl << "Value_th before focusing="
        << global_value_th << std::endl;

    double w1 = global_min_freq, w2 = global_max_freq;
    if (adjust_freq)
    {
        double w1U, w2U, w1V, w2V;
        compute_central_region(w1U, w2U, global_ctfmodel.azimuthal_angle);
        compute_central_region(w1V, w2V, global_ctfmodel.azimuthal_angle + 90);
        w1 = XMIPP_MIN(w1U, w1V);
        w2 = XMIPP_MAX(w2U, w2V);
        global_min_freq = XMIPP_MAX(global_min_freq, w1 - 0.05);
        global_max_freq = XMIPP_MIN(global_max_freq, w2 + 0.01);
    }

    // Compute maximum value within central region
    if (adjust_th)
    {
        Image<double> save;
        generateModelSoFar(save);
        double max_val = 0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
        {
            double w = global_w_digfreq(i, j);
            if (w >= w1 && w <= w2)
                max_val = XMIPP_MAX(max_val, save()(i, j));
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
void estimate_background_sqrt_parameters()
{
    if (global_prm->show_optimization)
        std::cout << "Computing first sqrt background ...\n";

    // Estimate the base line taking the value of the CTF
    // for the maximum X and Y frequencies
    double base_line = 0;
    int N = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    if (global_w_digfreq(i, j) > 0.4)
    {
        N++;
        base_line += (*f)(i, j);
    }
    global_ctfmodel.base_line = base_line / N;

    // Find the linear least squares solution for the sqrt part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_mask(i, j) <= 0)
            continue;

        // Compute weight for this point
        double weight = 1 + global_max_freq - global_w_digfreq(i, j);

        // Compute error
        global_ctfmodel.precomputeValues(global_x_contfreq(i, j),
                                         global_y_contfreq(i, j));
        double explained = global_ctfmodel.getValueNoiseAt();
        double unexplained = (*f)(i, j) - explained;
        if (unexplained <= 0)
            continue;
        unexplained = log(unexplained);

        double X = -sqrt(global_w_contfreq(i, j));
        A(0, 0) += weight * X * X;
        A(0, 1) += weight * X;
        A(1, 1) += weight * 1;
        b(0) += X * weight * unexplained;
        b(1) += weight * unexplained;
    }
    A(1, 0) = A(0, 1);

    b = A.inv() * b;

    global_ctfmodel.sqU = global_ctfmodel.sqV = b(0);
    global_ctfmodel.sqrt_K = exp(b(1));
    global_ctfmodel.sqrt_angle = 0;

    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (global_prm->show_optimization)
    {
        std::cout << "First SQRT Fit:\n" << global_ctfmodel << std::endl;
        saveIntermediateResults("step01a_first_sqrt_fit");
    }

    // Now optimize .........................................................
    double fitness;
    Matrix1D<double> steps;
    steps.resize(SQRT_CTF_PARAMETERS);
    steps.initConstant(0);
    steps(0) = steps(1) = steps(2) = 1;
    steps(3) = steps(4) = 1;

    // Optimize without penalization
    if (global_prm->show_optimization)
        std::cout << "Looking for best fitting sqrt ...\n";
    global_penalize = false;
    int iter;
    powellOptimizer(*global_adjust, FIRST_SQRT_PARAMETER + 1,
                    SQRT_CTF_PARAMETERS, &CTF_fitness, NULL, 0.05, fitness, iter, steps,
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
                        SQRT_CTF_PARAMETERS, &CTF_fitness, NULL, 0.05, fitness, iter,
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
        saveIntermediateResults("step01b_best_penalized_sqrt_fit");
    }

    center_optimization_focus(false, true, 1.5);
}

// Estimate gaussian parameters --------------------------------------------
//#define DEBUG
void estimate_background_gauss_parameters()
{

    if (global_prm->show_optimization)
        std::cout << "Computing first background Gaussian parameters ...\n";

    // Compute radial averages
    MultidimArray<double> radial_CTFmodel_avg(YSIZE(*f) / 2);
    MultidimArray<double> radial_CTFampl_avg(YSIZE(*f) / 2);
    MultidimArray<int> radial_N(YSIZE(*f) / 2);
    double w_max_gauss = 0.25;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_mask(i, j) <= 0)
            continue;
        double w = global_w_digfreq(i, j);
        if (w > w_max_gauss)
            continue;

        int r = FLOOR(w * (double)YSIZE(*f));
        global_ctfmodel.precomputeValues(global_x_contfreq(i, j),
                                         global_y_contfreq(i, j));
        radial_CTFmodel_avg(r) += global_ctfmodel.getValueNoiseAt();
        radial_CTFampl_avg(r) += (*f)(i, j);
        radial_N(r)++;
    }

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
    double error2_min = 0, wmin;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        double w = global_w_digfreq(i, 0);

        if (error(i) < 0 && first)
            continue;
        else if (error(i) < 0)
            break;
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
    double error2_max = 0, wmax, fmax;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        double w = global_w_digfreq(i, 0);
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
    fmax = global_ctfmodel.cV = global_ctfmodel.cU = wmax / global_prm->Tm;
#ifdef DEBUG

    std::cout << "Freq of the maximum error: " << wmax << " " << fmax << std::endl;
#endif

    // Find the linear least squares solution for the gauss part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_mask(i, j) <= 0)
            continue;
        if (global_w_digfreq(i, j) > wmin)
            continue;
        double fmod = global_w_contfreq(i, j);

        // Compute weight for this point
        double weight = 1 + global_max_freq - global_w_digfreq(i, j);

        // Compute error
        global_ctfmodel.precomputeValues(global_x_contfreq(i, j),
                                         global_y_contfreq(i, j));
        double explained = global_ctfmodel.getValueNoiseAt();

        double unexplained = (*f)(i, j) - explained;
        if (unexplained <= 0)
            continue;
        unexplained = log(unexplained);
        double F = -(fmod - fmax) * (fmod - fmax);

        A(0, 0) += weight * 1;
        A(0, 1) += weight * F;
        A(1, 1) += weight * F * F;
        b(0) += weight * unexplained;
        b(1) += F * weight * unexplained;
    }

    A(1, 0) = A(0, 1);
    if ( (A(0, 0)== 0) && (A(1, 0)== 0) && (A(1, 1)== 0))
    {
        std::cout << "A matrix es zero" << std::endl;
    }
    else
    {
        b = A.inv() * b;
        global_ctfmodel.sigmaU = XMIPP_MIN(fabs(b(1)), 95e3); // This value should be
        global_ctfmodel.sigmaV = XMIPP_MIN(fabs(b(1)), 95e3); // conformant with the physical
        // meaning routine in CTF.cc
        global_ctfmodel.gaussian_K = exp(b(0));
        // Store the CTF values in global_prm->adjust
        global_ctfmodel.forcePhysicalMeaning();
        COPY_ctfmodel_TO_CURRENT_GUESS;

        if (global_prm->show_optimization)
        {
            std::cout << "First Background Fit:\n" << global_ctfmodel << std::endl;
            saveIntermediateResults("step01c_first_background_fit");
        }
        center_optimization_focus(false, true, 1.5);
    }
}
#undef DEBUG

// Estimate second gaussian parameters -------------------------------------
//#define DEBUG
void estimate_background_gauss_parameters2()
{
    if (global_prm->show_optimization)
        std::cout << "Computing first background Gaussian2 parameters ...\n";

    // Compute radial averages
    MultidimArray<double> radial_CTFmodel_avg(YSIZE(*f) / 2);
    MultidimArray<double> radial_CTFampl_avg(YSIZE(*f) / 2);
    MultidimArray<int> radial_N(YSIZE(*f) / 2);
    double w_max_gauss = 0.25;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_mask(i, j) <= 0)
            continue;
        double w = global_w_digfreq(i, j);
        if (w > w_max_gauss)
            continue;

        int r = FLOOR(w * (double)YSIZE(*f));
        double f_x = DIRECT_A2D_ELEM(global_x_contfreq, i, j);
        double f_y = DIRECT_A2D_ELEM(global_y_contfreq, i, j);
        global_ctfmodel.precomputeValues(f_x, f_y);
        double bg = global_ctfmodel.getValueNoiseAt();
        double envelope = global_ctfmodel.getValueDampingAt();
        double ctf_without_damping =
            global_ctfmodel.getValuePureWithoutDampingAt();
        double ctf_with_damping = envelope * ctf_without_damping;
        double ctf2_th = bg + ctf_with_damping * ctf_with_damping;
        radial_CTFmodel_avg(r) += ctf2_th;
        radial_CTFampl_avg(r) += (*f)(i, j);
        radial_N(r)++;
    }

    // Compute the average radial error
    MultidimArray <double> error;
    error.initZeros(radial_CTFmodel_avg);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        error(i) = (radial_CTFampl_avg(i) - radial_CTFmodel_avg(i))
                   / radial_N(i);
    }
#ifdef DEBUG
    std::cout << "Error:\n" << error << std::endl;
#endif

    // Compute the frequency of the minimum error
    double wmin = 0.15;

    // Compute the maximum (negative) radial error
    double error_max = 0, wmax, fmax;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        double w = global_w_digfreq(i, 0);
        if (w > wmin)
            break;
        if (error(i) < error_max)
        {
            wmax = w;
            error_max = error(i);
        }
    }
    fmax = global_ctfmodel.cV2 = global_ctfmodel.cU2 = wmax / global_prm->Tm;
#ifdef DEBUG

    std::cout << "Freq of the maximum error: " << wmax << " " << fmax << std::endl;
#endif

    // Find the linear least squares solution for the gauss part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();
    int N = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_mask(i, j) <= 0)
            continue;
        if (global_w_digfreq(i, j) > wmin)
            continue;
        double fmod = global_w_contfreq(i, j);

        // Compute the zero on the direction of this point
        Matrix1D<double> u(2), fzero(2);
        XX(u) = global_x_contfreq(i, j) / fmod;
        YY(u) = global_y_contfreq(i, j) / fmod;
        global_ctfmodel.lookFor(1, u, fzero, 0);
        if (fmod > fzero.module())
            continue;

        // Compute weight for this point
        double weight = 1 + global_max_freq - global_w_digfreq(i, j);

        // Compute error
        double f_x = DIRECT_A2D_ELEM(global_x_contfreq, i, j);
        double f_y = DIRECT_A2D_ELEM(global_y_contfreq, i, j);
        global_ctfmodel.precomputeValues(f_x, f_y);
        double bg = global_ctfmodel.getValueNoiseAt();
        double envelope = global_ctfmodel.getValueDampingAt();
        double ctf_without_damping =
            global_ctfmodel.getValuePureWithoutDampingAt();
        double ctf_with_damping = envelope * ctf_without_damping;
        double ctf2_th = bg + ctf_with_damping * ctf_with_damping;
        double explained = ctf2_th;
        double unexplained = explained - (*f)(i, j);

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
    if (N != 0)
    {
        A(1, 0) = A(0, 1);

        double det=A.det();
        if (fabs(det)>1e-9)
        {
            b = A.inv() * b;
            global_ctfmodel.sigmaU2 = XMIPP_MIN(fabs(b(1)), 95e3); // This value should be
            global_ctfmodel.sigmaV2 = XMIPP_MIN(fabs(b(1)), 95e3); // conformant with the physical
            // meaning routine in CTF.cc
            global_ctfmodel.gaussian_K2 = exp(b(0));
        }
        else
        {
            global_ctfmodel.sigmaU2 = global_ctfmodel.sigmaV2 = 0;
            global_ctfmodel.gaussian_K2 = 0;
        }
    }
    else
    {
        global_ctfmodel.sigmaU2 = global_ctfmodel.sigmaV2 = 0;
        global_ctfmodel.gaussian_K2 = 0;
    }

    // Store the CTF values in global_prm->adjust
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

#ifdef DEBUG
    // Check
    FOR_ALL_ELEMENTS_IN_ARRAY2D(global_w_digfreq)
    {
        if (global_mask(i, j)<=0)
            continue;
        if (global_w_digfreq(i, j) > wmin)
            continue;
        double fmod = global_w_contfreq(i, j);

        // Compute the zero on the direction of this point
        Matrix1D<double> u(2), fzero(2);
        XX(u) = global_x_contfreq(i, j) / fmod;
        YY(u) = global_y_contfreq(i, j) / fmod;
        global_ctfmodel.zero(1, u, fzero);
        if (fmod > fzero.module())
            continue;

        // Compute error
        double f_x = DIRECT_A2D_ELEM(global_x_contfreq, i, j);
        double f_y = DIRECT_A2D_ELEM(global_y_contfreq, i, j);
        double bg = global_ctfmodel.getValueNoiseAt(f_x, f_y);
        double envelope = global_ctfmodel.getValueDampingAt(f_x, f_y);
        double ctf_without_damping = global_ctfmodel.getValuePureWithoutDampingAt(f_x, f_y);
        double ctf_with_damping = envelope * ctf_without_damping;
        double ctf2_th = bg + ctf_with_damping * ctf_with_damping;
        double explained = ctf2_th;
        double unexplained = explained - (*f)(i, j);

        if (unexplained <= 0)
            continue;
        std::cout << fmod << " " << unexplained << " "
        << global_ctfmodel.gaussian_K2*exp(-global_ctfmodel.sigmaU2*
                                           (fmod - fmax)*(fmod - fmax)) << std::endl;
    }
#endif

    if (global_prm->show_optimization)
    {
        std::cout << "First Background Gaussian 2 Fit:\n" << global_ctfmodel
        << std::endl;
        saveIntermediateResults("step04a_first_background2_fit");
    }
}
#undef DEBUG

// Estimate envelope parameters --------------------------------------------
//#define DEBUG
void estimate_envelope_parameters()
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
    COPY_ctfmodel_TO_CURRENT_GUESS;

    // Now optimize the envelope
    global_penalize = false;
    int iter;
    double fitness;
    Matrix1D<double> steps;
    steps.resize(ENVELOPE_PARAMETERS);
    steps.initConstant(1);
    steps(1) = 0; // Do not optimize Cs
    steps(5) = 0; // Do not optimize for alpha, since Ealpha depends on the
    // defocus
    if (global_prm->modelSimplification >= 1)
        steps(6) = steps(7) = 0; // Do not optimize DeltaF and DeltaR
    powellOptimizer(*global_adjust, FIRST_ENVELOPE_PARAMETER + 1,
                    ENVELOPE_PARAMETERS, &CTF_fitness, NULL, 0.05, fitness, iter, steps,
                    global_prm->show_optimization);

    // Keep the result in global_prm->adjust
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (global_prm->show_optimization)
    {
        std::cout << "Best envelope Fit:\n" << global_ctfmodel << std::endl;
        saveIntermediateResults("step02a_best_envelope_fit");
    }

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
                        ENVELOPE_PARAMETERS, &CTF_fitness, NULL, 0.05, fitness, iter,
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
        std::cout << "Best envelope Fit:\n" << global_ctfmodel << std::endl;
        saveIntermediateResults("step02b_best_penalized_envelope_fit");
    }
}
#undef DEBUG

// Estimate defoci ---------------------------------------------------------
void showFirstDefoci()
{
    if (global_prm->show_optimization)
    {
        std::cout << "First defocus Fit:\n" << global_ctfmodel << std::endl;
        saveIntermediateResults("step03a_first_defocus_fit");
        global_prm->enhanced_ctftomodel.write("step03a_enhanced_PSD.xmp");
        Image<double> save, save2, save3;
        save().resize(YSIZE(global_w_digfreq), XSIZE(global_w_digfreq));
        save2().resize(save());
        save3().resize(save());
        FOR_ALL_ELEMENTS_IN_ARRAY2D(save())
        {
            save()(i, j) = global_prm->enhanced_ctftomodel()(i, j);
            double f_x = DIRECT_A2D_ELEM(global_x_contfreq, i, j);
            double f_y = DIRECT_A2D_ELEM(global_y_contfreq, i, j);
            global_ctfmodel.precomputeValues(f_x, f_y);
            double ctf_without_damping =
                global_ctfmodel.getValuePureWithoutDampingAt();
            save2()(i, j) = ctf_without_damping * ctf_without_damping;
            save3()(i, j) = -global_prm->enhanced_ctftomodel()(i, j)
                            * ctf_without_damping * ctf_without_damping;
        }
        save.write("step03a_enhanced_PSD.xmp");
        save2.write("step03a_fitted_CTF.xmp");
        save3.write("step03a_superposition.xmp");
    }
}

//#define DEBUG
void estimate_defoci()
{
    if (global_prm->show_optimization)
        std::cout << "Looking for first defoci ...\n";
    double best_defocusU, best_defocusV, best_angle, best_K;
    double best_error = global_heavy_penalization * 1.1;
    bool first = true;
    int i, j;
    double defocusV, defocusU;

    double defocusV0 = 1e3, defocusU0 = 1e3;
    double defocusVF = 100e3, defocusUF = 100e3;
    double initial_defocusStep = 8e3;
    MultidimArray<double> error;

    // Check if there is no initial guess
    double min_allowed_defocusU = 1e3, max_allowed_defocusU = 100e3;
    double min_allowed_defocusV = 1e3, max_allowed_defocusV = 100e3;
    if (global_prm->initial_ctfmodel.DeltafU != 0)
    {
        initial_defocusStep = std::min(global_prm->defocus_range,20000.0);
        defocusU0 = std::max(
                        1e3,
                        global_prm->initial_ctfmodel.DeltafU
                        - global_prm->defocus_range);
        double maxDeviation = std::max(global_prm->defocus_range,
                                       0.25 * global_prm->initial_ctfmodel.DeltafU);
        max_allowed_defocusU = std::min(100e3,
                                        global_prm->initial_ctfmodel.DeltafU + maxDeviation);
        defocusUF = std::min(
                        150e3,
                        global_prm->initial_ctfmodel.DeltafU
                        + global_prm->defocus_range);
        min_allowed_defocusU = std::max(1e3,
                                        global_prm->initial_ctfmodel.DeltafU - maxDeviation);
        if (global_prm->initial_ctfmodel.DeltafV == 0)
        {
            defocusV0 = defocusU0;
            min_allowed_defocusV = min_allowed_defocusU;
            defocusVF = defocusUF;
            max_allowed_defocusV = max_allowed_defocusU;
        }
        else
        {
            defocusV0 = std::max(
                            1e3,
                            global_prm->initial_ctfmodel.DeltafV
                            - global_prm->defocus_range);
            max_allowed_defocusV = std::max(100e3,
                                            global_prm->initial_ctfmodel.DeltafV + maxDeviation);
            defocusVF = std::min(
                            150e3,
                            global_prm->initial_ctfmodel.DeltafV
                            + global_prm->defocus_range);
            min_allowed_defocusV = std::max(1e3,
                                            global_prm->initial_ctfmodel.DeltafV - maxDeviation);
        }
    }

    double K_so_far = global_ctfmodel.K;
    Matrix1D<double> steps(DEFOCUS_PARAMETERS);
    steps.initConstant(1);
    steps(3) = 0; // Do not optimize kV
    steps(4) = 0; // Do not optimize K
    for (double defocusStep = initial_defocusStep;
         defocusStep >= std::min(5000., global_prm->defocus_range / 2);
         defocusStep /= 2)
    {
        error.resize(CEIL((defocusVF - defocusV0) / defocusStep + 1),
                     CEIL((defocusUF - defocusU0) / defocusStep + 1));
        error.initConstant(global_heavy_penalization);
        if (global_prm->show_optimization)
            std::cout << "V=[" << defocusV0 << "," << defocusVF << "]\n"
            << "U=[" << defocusU0 << "," << defocusUF << "]\n"
            << "Defocus step=" << defocusStep << std::endl;
        for (defocusV = defocusV0, i = 0; defocusV <= defocusVF; defocusV +=
                 defocusStep, i++)
        {
            for (defocusU = defocusU0, j = 0; defocusU <= defocusUF; defocusU +=
                     defocusStep, j++)
            {
                bool first_angle = true;
                if (fabs(defocusU - defocusV) > 30e3)
                {
                    error(i, j) = global_heavy_penalization;
                    continue;
                }
                for (double angle = 0; angle < 180; angle += 45)
                {
                    int iter;
                    double fitness;

                    (*global_adjust)(0) = defocusU;
                    (*global_adjust)(1) = defocusV;
                    (*global_adjust)(2) = angle;
                    (*global_adjust)(4) = K_so_far;

                    powellOptimizer(*global_adjust, FIRST_DEFOCUS_PARAMETER + 1,
                                    DEFOCUS_PARAMETERS, &CTF_fitness, NULL, 0.05,
                                    fitness, iter, steps, false);

                    if ((first_angle || fitness < error(i, j))
                        && (global_ctfmodel.DeltafU >= min_allowed_defocusU
                            && global_ctfmodel.DeltafU
                            <= max_allowed_defocusU
                            && global_ctfmodel.DeltafV
                            >= min_allowed_defocusV
                            && global_ctfmodel.DeltafV
                            <= max_allowed_defocusV))
                    {
                        error(i, j) = fitness;
                        first_angle = false;
                        if (error(i, j) < best_error || first)
                        {
                            best_error = error(i, j);
                            best_defocusU = global_ctfmodel.DeltafU;
                            best_defocusV = global_ctfmodel.DeltafV;
                            best_angle = global_ctfmodel.azimuthal_angle;
                            best_K = global_ctfmodel.K;
                            first = false;
                            if (global_prm->show_optimization)
                            {
                                std::cout << "    (DefocusU,DefocusV)=("
                                << defocusU << "," << defocusV
                                << "), ang=" << angle << " --> ("
                                << global_ctfmodel.DeltafU << ","
                                << global_ctfmodel.DeltafV << "),"
                                << global_ctfmodel.azimuthal_angle
                                << " K=" << global_ctfmodel.K
                                << " error=" << error(i, j)
                                << std::endl;
#ifdef DEBUG

                                global_show=3;
                                CTF_fitness(global_adjust->vdata-1,NULL);
                                global_show=0;

                                Image<double> save;
                                save() = global_prm->enhanced_ctftomodel();
                                save.write("PPPenhanced.xmp");
                                for (int i = 0; i < YSIZE(global_w_digfreq); i += 1)
                                    for (int j = 0; j < XSIZE(global_w_digfreq); j += 1)
                                    {
                                        if (DIRECT_A2D_ELEM(global_mask, i, j)<=0)
                                            continue;
                                        double f_x = DIRECT_A2D_ELEM(global_x_contfreq, i, j);
                                        double f_y = DIRECT_A2D_ELEM(global_y_contfreq, i, j);
                                        double envelope = global_ctfmodel.getValueDampingAt(f_x, f_y);
                                        double ctf_without_damping = global_ctfmodel.getValuePureWithoutDampingAt(f_x, f_y);
                                        double ctf_with_damping = envelope * ctf_without_damping;
                                        double ctf_with_damping2 = ctf_with_damping * ctf_with_damping;
                                        save(i, j) = ctf_with_damping2;
                                    }
                                save.write("PPPctf2_with_damping2.xmp");
                                std::cout << "Press any key\n";
                                char c;
                                std::cin >> c;
#endif

                            }
                        }
                    }
                }
            }
        }

        // Compute the range of the errors
        double errmin = error(0, 0), errmax = error(0, 0);
        bool aValidErrorHasBeenFound=false;
        for (int ii = STARTINGY(error); ii <= FINISHINGY(error); ii++)
            for (int jj = STARTINGX(error); jj <= FINISHINGX(error); jj++)
            {
                if (error(ii, jj) != global_heavy_penalization)
                {
                	aValidErrorHasBeenFound=true;
                    if (error(ii, jj) < errmin)
                        errmin = error(ii, jj);
                    else if (errmax == global_heavy_penalization)
                        errmax = error(ii, jj);
                    else if (error(ii, jj) > errmax)
                        errmax = error(ii, jj);
                }
            }
        if (global_prm->show_optimization)
            std::cout << "Error matrix\n" << error << std::endl;
        if (!aValidErrorHasBeenFound)
        	REPORT_ERROR(ERR_NUMERICAL,"Cannot find any good defocus within the given range");

        // Find those defoci which are within a 10% of the maximum
        if (global_show >= 2)
            std::cout << "Range=" << errmax - errmin << std::endl;
        double best_defocusVmin = best_defocusV, best_defocusVmax =
                                      best_defocusV;
        double best_defocusUmin = best_defocusU, best_defocusUmax =
                                      best_defocusU;
        for (defocusV = defocusV0, i = 0; defocusV <= defocusVF; defocusV +=
                 defocusStep, i++)
        {
            for (defocusU = defocusU0, j = 0; defocusU <= defocusUF; defocusU +=
                     defocusStep, j++)
            {
                if (global_show >= 2)
                    std::cout << i << "," << j << " " << error(i, j) << " "
                    << defocusU << " " << defocusV << std::endl
                    << best_defocusUmin << " " << best_defocusUmax
                    << std::endl << best_defocusVmin << " "
                    << best_defocusVmax << std::endl;
                if (fabs((error(i, j) - errmin) / (errmax - errmin)) <= 0.1)
                {
                    if (defocusV < best_defocusVmin)
                        best_defocusVmin = defocusV;
                    if (defocusU < best_defocusUmin)
                        best_defocusUmin = defocusU;
                    if (defocusV > best_defocusVmax)
                        best_defocusVmax = defocusV;
                    if (defocusU > best_defocusUmax)
                        best_defocusUmax = defocusU;
                }
            }
        }

        defocusVF = std::min(max_allowed_defocusV,
                             best_defocusVmax + defocusStep);
        defocusV0 = std::max(min_allowed_defocusV,
                             best_defocusVmin - defocusStep);
        defocusUF = std::min(max_allowed_defocusU,
                             best_defocusUmax + defocusStep);
        defocusU0 = std::max(min_allowed_defocusU,
                             best_defocusUmin - defocusStep);
        i = j = 0;
        if (global_show >= 2)
        {
            Image<double> save;
            save() = error;
            save.write("error.xmp");
            std::cout << "Press any key: Error saved\n";
            char c;
            std::cin >> c;
        }
    }

    global_ctfmodel.DeltafU = best_defocusU;
    global_ctfmodel.DeltafV = best_defocusV;
    global_ctfmodel.azimuthal_angle = best_angle;
    global_ctfmodel.K = best_K;

    // Keep the result in global_prm->adjust
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;
    global_ctfmodel_defoci = global_ctfmodel;

    showFirstDefoci();
}
#undef DEBUG

// Estimate defoci with Zernike and SPTH transform---------------------------------------------
void estimate_defoci_Zernike(MultidimArray<double> &psdToModelFullSize, double min_freq, double max_freq, double Tm,
                             double kV, double lambdaPhase, int sizeWindowPhase,
                             double &defocusU, double &defocusV, double &ellipseAngle, int verbose)
{
	// Center enhanced PSD
    MultidimArray<double> centeredEnhancedPSD=psdToModelFullSize;
    CenterFFT(centeredEnhancedPSD,true);

    // Estimate phase, modulation and Zernikes
    MultidimArray<double> mod, phase;
    Matrix1D<double> coefs(13);
    coefs.initConstant(0);
    VEC_ELEM(coefs,0) = 1;
    VEC_ELEM(coefs,3) = 1;
    VEC_ELEM(coefs,4) = 1;
    VEC_ELEM(coefs,5) = 1;
    VEC_ELEM(coefs,12) =1;

    int x=(int)((0.3*max_freq+0.7*min_freq)*std::cos(PI/4)*XSIZE(centeredEnhancedPSD)+XSIZE(centeredEnhancedPSD)/2);
    DEBUG_TEXTFILE(formatString("Zernike1 %d",x));
    DEBUG_TEXTFILE(formatString("centeredEnhancedPSD80x80 %f",centeredEnhancedPSD(80,80)));
    DEBUG_TEXTFILE(formatString("centeredEnhancedPSD120x120 %f",centeredEnhancedPSD(120,120)));
    DEBUG_TEXTFILE(formatString("centeredEnhancedPSD160x160 %f",centeredEnhancedPSD(160,160)));

    int numElem = 10;
    kV = kV*1000;
    double K_so_far = global_ctfmodel.K;
    double lambda=12.2643247/std::sqrt(kV*(1.+0.978466e-6*kV));
    double Z8;
    double Z3;
    double Z4;
    double Z5;
    double eAngle;
    double deFocusAvg;
    double deFocusDiff;

    max_freq = max_freq*0.8;
    double fmax = max_freq;

    Matrix1D<double> arrayDefocusDiff(numElem);
    arrayDefocusDiff.initZeros();

    Matrix1D<double> arrayDefocusAvg(numElem);
    arrayDefocusAvg.initZeros();

    Matrix1D<double> arrayError(numElem);
    arrayError.initConstant(-1);

    int iter;
    double fitness;
    Matrix1D<double> steps(DEFOCUS_PARAMETERS);
    steps.initConstant(1);
    steps(3) = 0; // Do not optimize kV
    steps(4) = 0; // Do not optimize K
    double fmaxStep = (max_freq-min_freq)/numElem;

    lambdaPhase = 0.8;
    sizeWindowPhase = 10;

    Matrix1D<double> initialGlobalAdjust = (*global_adjust);

    for (int i = 1; i < numElem; i++)
    {
        if ( ( ((fmax - min_freq)/min_freq) > 0.5))
        {
            demodulate(centeredEnhancedPSD,lambdaPhase,sizeWindowPhase,
                          x,x,
                          (int)(min_freq*XSIZE(centeredEnhancedPSD)),
                          (int)(fmax*XSIZE(centeredEnhancedPSD)),
                          phase, mod, coefs, 0);

            Z8=VEC_ELEM(coefs,4);
            Z3=VEC_ELEM(coefs,12);
            Z4=VEC_ELEM(coefs,3);
            Z5=VEC_ELEM(coefs,5);

            eAngle = 0.5*RAD2DEG(std::atan2(Z5,Z4))+90.0;
            deFocusAvg  =  fabs(2*Tm*Tm*(2*Z3-6*Z8)/(PI*lambda));
            deFocusDiff =  fabs(2*Tm*Tm*(std::sqrt(Z4*Z4+Z5*Z5))/(PI*lambda));

            coefs.initConstant(0);
            VEC_ELEM(coefs,0) = 1;
            VEC_ELEM(coefs,3) = 1;
            VEC_ELEM(coefs,4) = 1;
            VEC_ELEM(coefs,5) = 1;
            VEC_ELEM(coefs,12) =1;

            fmax -= fmaxStep;

            (*global_adjust)(0) = deFocusAvg+deFocusDiff;
            (*global_adjust)(1) = deFocusAvg-deFocusDiff;
            (*global_adjust)(2) = eAngle;
            (*global_adjust)(4) = K_so_far;
            (*global_adjust)(6) = 2;


            fitness =0;
            powellOptimizer(*global_adjust, FIRST_DEFOCUS_PARAMETER + 1,
                            DEFOCUS_PARAMETERS, &CTF_fitness, NULL, 0.05,
                            fitness, iter, steps, false);

            VEC_ELEM(arrayDefocusAvg,i)  = ((*global_adjust)(0) +(*global_adjust)(1))/2;
            VEC_ELEM(arrayDefocusDiff,i) = ((*global_adjust)(0) -(*global_adjust)(1))/2;
            VEC_ELEM(arrayError,i) = (-1)*fitness;

        }
    }

    int maxInd=arrayError.maxIndex();

    while ( (VEC_ELEM(arrayDefocusAvg,maxInd) < 3000) || ((VEC_ELEM(arrayDefocusAvg,maxInd) > 50000) && VEC_ELEM(arrayError,maxInd)>-1e3 ))
    {
        VEC_ELEM(arrayError,maxInd) = -1e3;
        VEC_ELEM(arrayDefocusAvg,maxInd) = global_prm->initial_ctfmodel.DeltafU;
        VEC_ELEM(arrayDefocusDiff,maxInd) = global_prm->initial_ctfmodel.DeltafV;
        maxInd=arrayError.maxIndex();
    }
    if (VEC_ELEM(arrayError,maxInd)<=-1e3)
    {
    	estimate_defoci();
    	return;
    }

    Matrix1D<double> arrayDefocusU(3);
    arrayDefocusU.initZeros();

    Matrix1D<double> arrayDefocusV(3);
    arrayDefocusV.initZeros();

    Matrix1D<double> arrayError2(3);
    arrayError2.initConstant(-1);

    //We want to take care about more parameters
    // We optimize for deltaU, deltaV
    (*global_adjust)(0) = VEC_ELEM(arrayDefocusAvg,maxInd)+VEC_ELEM(arrayDefocusDiff,maxInd);
    (*global_adjust)(1) = VEC_ELEM(arrayDefocusAvg,maxInd)-VEC_ELEM(arrayDefocusDiff,maxInd);
    (*global_adjust)(2) = eAngle;
    (*global_adjust)(4) = K_so_far;
    (*global_adjust)(6) = 2;

    fitness =0;
    powellOptimizer(*global_adjust, FIRST_DEFOCUS_PARAMETER + 1,
                    DEFOCUS_PARAMETERS, &CTF_fitness, NULL, 0.05,
                    fitness, iter, steps, false);

    VEC_ELEM(arrayDefocusU,0) = (*global_adjust)(0);
    VEC_ELEM(arrayDefocusV,0) = (*global_adjust)(1);
    VEC_ELEM(arrayError2,0) = (-1)*fitness;

    // We optimize for deltaU, deltaU
    (*global_adjust)(0) = VEC_ELEM(arrayDefocusAvg,maxInd)+VEC_ELEM(arrayDefocusDiff,maxInd);
    (*global_adjust)(1) = VEC_ELEM(arrayDefocusAvg,maxInd)+VEC_ELEM(arrayDefocusDiff,maxInd);
    (*global_adjust)(2) = eAngle;
    (*global_adjust)(4) = K_so_far;
    (*global_adjust)(6) = 2;

    fitness =0;
    powellOptimizer(*global_adjust, FIRST_DEFOCUS_PARAMETER + 1,
                    DEFOCUS_PARAMETERS, &CTF_fitness, NULL, 0.05,
                    fitness, iter, steps, false);

    VEC_ELEM(arrayDefocusU,1) = (*global_adjust)(0);
    VEC_ELEM(arrayDefocusV,1) = (*global_adjust)(1);
    VEC_ELEM(arrayError2,1) = (-1)*fitness;

    // We optimize for deltaV, deltaV
    (*global_adjust)(0) = VEC_ELEM(arrayDefocusAvg,maxInd)-VEC_ELEM(arrayDefocusDiff,maxInd);
    (*global_adjust)(1) = VEC_ELEM(arrayDefocusAvg,maxInd)-VEC_ELEM(arrayDefocusDiff,maxInd);
    (*global_adjust)(2) = eAngle;
    (*global_adjust)(4) = K_so_far;
    (*global_adjust)(6) = 2;

    fitness =0;
    powellOptimizer(*global_adjust, FIRST_DEFOCUS_PARAMETER + 1,
                    DEFOCUS_PARAMETERS, &CTF_fitness, NULL, 0.05,
                    fitness, iter, steps, false);

    VEC_ELEM(arrayDefocusU,2) = (*global_adjust)(0);
    VEC_ELEM(arrayDefocusV,2) = (*global_adjust)(1);
    VEC_ELEM(arrayError2,2) = (-1)*fitness;

    //Here we select the best one
    maxInd=arrayError2.maxIndex();
    defocusU = VEC_ELEM(arrayDefocusU,maxInd);
    defocusV = VEC_ELEM(arrayDefocusV,maxInd);

    (*global_adjust)(0) = defocusU;
    (*global_adjust)(1) = defocusV;
    (*global_adjust)(2) = eAngle;
    (*global_adjust)(4) = K_so_far;
    (*global_adjust)(6) = 2;

    while ( (0.5*(defocusU+defocusV) < 2500) || (0.5*(defocusU+defocusV) > 60000) )
    {
        VEC_ELEM(arrayError2,maxInd) = -1e3;
        VEC_ELEM(arrayDefocusU,maxInd) = global_prm->initial_ctfmodel.DeltafU;
        VEC_ELEM(arrayDefocusV,maxInd) = global_prm->initial_ctfmodel.DeltafV;

        maxInd=arrayError2.maxIndex();
        defocusU = VEC_ELEM(arrayDefocusU,maxInd);
        defocusV = VEC_ELEM(arrayDefocusV,maxInd);
        (*global_adjust)(0) = defocusU;
        (*global_adjust)(1) = defocusV;
        (*global_adjust)(2) = eAngle;
        (*global_adjust)(4) = K_so_far;
        (*global_adjust)(6) = 2;

    }


    if (VEC_ELEM(arrayError2,maxInd) <= 0)
    {
        COPY_ctfmodel_TO_CURRENT_GUESS;
        global_ctfmodel_defoci = global_ctfmodel;

        global_action = 5;

        steps.resize(ALL_CTF_PARAMETERS);
        steps.initConstant(1);
        steps(3) = 0; // kV
        steps(5) = 0; // The spherical aberration (Cs) is not optimized
        if (global_ctfmodel.Q0!=0)
        	steps(12)=0;

        COPY_ctfmodel_TO_CURRENT_GUESS;

        global_evaluation_reduction = 2;
        powellOptimizer(*global_adjust, 0 + 1, ALL_CTF_PARAMETERS, &CTF_fitness,
                        NULL, 0.01, fitness, iter, steps, global_prm->show_optimization);
        COPY_ctfmodel_TO_CURRENT_GUESS;

        global_show=0;
        global_action = 3;
        global_evaluation_reduction = 1;

        double error = -CTF_fitness(global_adjust->vdata-1,NULL);
        if ( error <= -0.1)
        {
            *global_adjust = initialGlobalAdjust;
            COPY_ctfmodel_TO_CURRENT_GUESS;
            //There is nothing to do and we have to perform an exhaustive search
#ifndef RELEASE_MODE
            std::cout << " Entering in estimate_defoci, Performing exhaustive defocus search (SLOW)" << std::endl;
#endif
            estimate_defoci();
        }
    }
}

void estimate_defoci_Zernike()
{
    if (global_prm->show_optimization)
        std::cout << "Looking for first defoci ...\n";

    DEBUG_TEXTFILE("Step 6.1");
    DEBUG_MODEL_TEXTFILE;
    estimate_defoci_Zernike(global_prm->enhanced_ctftomodel_fullsize(),
                            global_prm->min_freq,global_prm->max_freq,global_prm->Tm,
                            global_prm->initial_ctfmodel.kV,
                            global_prm->lambdaPhase,global_prm->sizeWindowPhase,
                            global_ctfmodel.DeltafU, global_ctfmodel.DeltafV, global_ctfmodel.azimuthal_angle, 0);
    DEBUG_TEXTFILE("Step 6.2");
    DEBUG_MODEL_TEXTFILE;

    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;
    global_ctfmodel_defoci = global_ctfmodel;

    showFirstDefoci();
}

/* Main routine ------------------------------------------------------------ */
//#define DEBUG
double ROUT_Adjust_CTF(ProgCTFEstimateFromPSD &prm,
                       CTFDescription &output_ctfmodel, bool standalone)
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
    global_heavy_penalization = f->computeMax() * XSIZE(*f) * YSIZE(*f);
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
        estimate_background_sqrt_parameters();
        estimate_background_gauss_parameters();
    }

    // Optimize the current background
    global_action = 1;
    global_penalize = true;
    global_current_penalty = global_penalty;
    steps.resize(BACKGROUND_CTF_PARAMETERS);
    steps.initConstant(1);
    if (!global_prm->modelSimplification >= 3)
        steps(7) = steps(8) = steps(10) = 0;
    powellOptimizer(*global_adjust, FIRST_SQRT_PARAMETER + 1,
                    BACKGROUND_CTF_PARAMETERS, &CTF_fitness, NULL, 0.01, fitness, iter,
                    steps, global_prm->show_optimization);

    // Make sure that the model has physical meaning
    // (In some machines due to numerical imprecission this check is necessary
    // at the end)
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (global_prm->show_optimization)
    {
        std::cout << "Best background Fit:\n" << global_ctfmodel << std::endl;
        saveIntermediateResults("step01d_best_background_fit");
    }
    DEBUG_TEXTFILE(formatString("Step 4: CTF_fitness=%f",CTF_fitness));
    DEBUG_MODEL_TEXTFILE;

    /************************************************************************/
    /* STEPs 5 and 6:  Find envelope which best fits the CTF                */
    /************************************************************************/
    global_action = 2;
    global_ctfmodel.enable_CTF = true;
    if (prm.initial_ctfmodel.K == 0)
    {
        global_ctfmodel.kV = prm.initial_ctfmodel.kV;
        global_ctfmodel.Cs = prm.initial_ctfmodel.Cs;
        if (prm.initial_ctfmodel.Q0 != 0)
            global_ctfmodel.Q0 = prm.initial_ctfmodel.Q0;
        estimate_envelope_parameters();
    }
    else
    {
        global_ctfmodel.K = prm.initial_ctfmodel.K;
        global_ctfmodel.kV = prm.initial_ctfmodel.kV;
        global_ctfmodel.DeltafU = prm.initial_ctfmodel.DeltafU;
        global_ctfmodel.DeltafV = prm.initial_ctfmodel.DeltafV;
        global_ctfmodel.azimuthal_angle = prm.initial_ctfmodel.azimuthal_angle;
        global_ctfmodel.Cs = prm.initial_ctfmodel.Cs;
        global_ctfmodel.Ca = prm.initial_ctfmodel.Ca;
        global_ctfmodel.espr = prm.initial_ctfmodel.espr;
        global_ctfmodel.ispr = prm.initial_ctfmodel.ispr;
        global_ctfmodel.alpha = prm.initial_ctfmodel.alpha;
        global_ctfmodel.DeltaF = prm.initial_ctfmodel.DeltaF;
        global_ctfmodel.DeltaR = prm.initial_ctfmodel.DeltaR;
        global_ctfmodel.Q0 = prm.initial_ctfmodel.Q0;
        COPY_ctfmodel_TO_CURRENT_GUESS;
    }
    DEBUG_TEXTFILE(formatString("Step 6: espr=%f",global_ctfmodel.espr));
    DEBUG_MODEL_TEXTFILE;
    /************************************************************************/
    /* STEP 7:  the defocus and angular parameters                          */
    /************************************************************************/

    global_action = 3;
    global_evaluation_reduction = 1;
    if (prm.fastDefocusEstimate)
        estimate_defoci_Zernike();
    else
        estimate_defoci();

    DEBUG_TEXTFILE(formatString("Step 7: DeltafU=%f",global_ctfmodel.DeltafU));
    DEBUG_TEXTFILE(formatString("Step 7: DeltafV=%f",global_ctfmodel.DeltafV));
    DEBUG_TEXTFILE(formatString("Step 7: azimutalAngle=%f",global_ctfmodel.azimuthal_angle));
    DEBUG_MODEL_TEXTFILE;

    //This line is to test the results obtained
    //exit(1);

    /************************************************************************/
    /* STEPs 9, 10 and 11: all parameters included second Gaussian          */
    /************************************************************************/
    global_action = 5;
    if (prm.modelSimplification < 2)
        estimate_background_gauss_parameters2();

    steps.resize(ALL_CTF_PARAMETERS);
    steps.initConstant(1);
    steps(3) = 0; // kV
    steps(5) = 0; // The spherical aberration (Cs) is not optimized
    if (prm.initial_ctfmodel.Q0 != 0)
        steps(12) = 0; // Q0
    if (prm.modelSimplification >= 3)
        steps(20) = steps(21) = steps(23) = 0;
    if (prm.modelSimplification >= 2)
        steps(24) = steps(25) = steps(26) = steps(27) = steps(28) = steps(29) =
                                                0;
    if (prm.modelSimplification >= 1)
        steps(10) = steps(11) = 0;

    powellOptimizer(*global_adjust, 0 + 1, ALL_CTF_PARAMETERS, &CTF_fitness,
                    NULL, 0.01, fitness, iter, steps, global_prm->show_optimization);

    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (global_prm->show_optimization)
    {
        std::cout << "Best fit with Gaussian2:\n" << global_ctfmodel
        << std::endl;
        saveIntermediateResults("step04b_best_fit_with_gaussian2");
    }

    global_evaluation_reduction = 2;
    powellOptimizer(*global_adjust, 0 + 1, ALL_CTF_PARAMETERS, &CTF_fitness,
                    NULL, 0.01, fitness, iter, steps, global_prm->show_optimization);
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    global_evaluation_reduction = 1;
    powellOptimizer(*global_adjust, 0 + 1, ALL_CTF_PARAMETERS, &CTF_fitness,
                    NULL, 0.005, fitness, iter, steps, global_prm->show_optimization);
    global_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (global_prm->show_optimization)
    {
        std::cout << "Best fit:\n" << global_ctfmodel << std::endl;
        saveIntermediateResults("step04c_best_fit");
    }
    DEBUG_TEXTFILE(formatString("Step 11: DeltafU=%f fitness=%f",global_ctfmodel.DeltafU,fitness));
    DEBUG_MODEL_TEXTFILE;

    //We adopt that always  DeltafU > DeltafV so if this is not the case we change the values and the angle
    if ( global_ctfmodel.DeltafV > global_ctfmodel.DeltafU)
    {
        double temp;
        temp = global_ctfmodel.DeltafU;
        global_ctfmodel.DeltafU = global_ctfmodel.DeltafV;
        global_ctfmodel.DeltafV = temp;
        global_ctfmodel.azimuthal_angle -= 90;
        COPY_ctfmodel_TO_CURRENT_GUESS;
    }

    /************************************************************************/
    /* STEP 12: Produce output                                              */
    /************************************************************************/
    global_action = 6;

    if (prm.fn_psd != "")
    {
        // Define mask between first and third zero
        global_mask_between_zeroes.initZeros(global_mask);
        Matrix1D<double> u(2), z1(2), z3(2);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(global_mask_between_zeroes)
        {
            VECTOR_R2(u, global_x_digfreq(i, j), global_y_digfreq(i, j));
            u /= u.module();
            global_ctfmodel.lookFor(1, u, z1, 0);
            global_ctfmodel.lookFor(3, u, z3, 0);
            if (z1.module() < global_w_contfreq(i, j)
                && global_w_contfreq(i, j) < z3.module())
                global_mask_between_zeroes(i, j) = 1;
        }

        // Evaluate the correlation in this region
        CTF_fitness(global_adjust->adaptForNumericalRecipes(), NULL);

        // Save results
        FileName fn_rootCTFPARAM = prm.fn_psd.withoutExtension();

        FileName fn_rootMODEL = fn_rootCTFPARAM;
        size_t atPosition=fn_rootCTFPARAM.find('@');

        if (atPosition!=std::string::npos)
        {
            fn_rootMODEL=formatString("%03d@%s",textToInteger(fn_rootCTFPARAM.substr(0, atPosition)),
                                      fn_rootCTFPARAM.substr(atPosition+1).c_str());
            fn_rootCTFPARAM=formatString("region%03d@%s",textToInteger(fn_rootCTFPARAM.substr(0, atPosition)),
                                         fn_rootCTFPARAM.substr(atPosition+1).c_str());
        }
        else
            fn_rootCTFPARAM=(String)"fullMicrograph@"+fn_rootCTFPARAM;

        saveIntermediateResults(fn_rootMODEL, false);
        global_ctfmodel.Tm /= prm.downsampleFactor;
        global_ctfmodel.azimuthal_angle = std::fmod(global_ctfmodel.azimuthal_angle,360.);
        global_ctfmodel.write(fn_rootCTFPARAM + ".ctfparam_tmp");
        MetaData MD;
        MD.read(fn_rootCTFPARAM + ".ctfparam_tmp");
        size_t id = MD.firstObject();
        MD.setValue(MDL_CTF_X0, (double)output_ctfmodel.x0*prm.Tm, id);
        MD.setValue(MDL_CTF_XF, (double)output_ctfmodel.xF*prm.Tm, id);
        MD.setValue(MDL_CTF_Y0, (double)output_ctfmodel.y0*prm.Tm, id);
        MD.setValue(MDL_CTF_YF, (double)output_ctfmodel.yF*prm.Tm, id);
        MD.setValue(MDL_CTF_CRIT_FITTINGSCORE, fitness, id);
        MD.setValue(MDL_CTF_CRIT_FITTINGCORR13, global_corr13, id);
        MD.setValue(MDL_CTF_DOWNSAMPLE_PERFORMED, prm.downsampleFactor, id);
        MD.write(fn_rootCTFPARAM + ".ctfparam",MD_APPEND);
        fn_rootCTFPARAM = fn_rootCTFPARAM + ".ctfparam_tmp";
        fn_rootCTFPARAM.deleteFile();
    }
    output_ctfmodel = global_ctfmodel;

    DEBUG_CLOSE_TEXTFILE;
    return fitness;
}

void ProgCTFEstimateFromPSD::run()
{
    CTFDescription ctfmodel;
    ROUT_Adjust_CTF(*this, ctfmodel);
}
