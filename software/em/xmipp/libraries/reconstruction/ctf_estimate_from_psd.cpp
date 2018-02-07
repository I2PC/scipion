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
#define ALL_CTF_PARAMETERS           38
#define CTF_PARAMETERS               30
#define PARAMETRIC_CTF_PARAMETERS    16
#define BACKGROUND_CTF_PARAMETERS    14
#define SQRT_CTF_PARAMETERS           8
#define ENVELOPE_PARAMETERS          11
#define DEFOCUS_PARAMETERS            5
#define FIRST_SQRT_PARAMETER         16
#define FIRST_ENVELOPE_PARAMETER      4
#define FIRST_DEFOCUS_PARAMETER       0

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
// CTF model and noise model
//CTFDescription ctfmodel_defoci;

}

using namespace AdjustCTF;

#define ASSIGN_CTF_PARAM(index, paramName) if (ia <= index && l > 0) { ctfmodel.paramName = p[index]; --l; }

/* Assign ctfmodel from a vector and viceversa ----------------------------- */
void ProgCTFEstimateFromPSD::assignCTFfromParameters(double *p, CTFDescription &ctfmodel, int ia,
                             int l, int modelSimplification)
{

	ctfmodel.Tm = Tm;

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
    ASSIGN_CTF_PARAM(12, envR0);
    ASSIGN_CTF_PARAM(13, envR1);
    ASSIGN_CTF_PARAM(14, envR2);
    ASSIGN_CTF_PARAM(15, Q0);
    ASSIGN_CTF_PARAM(16, base_line);
    ASSIGN_CTF_PARAM(17, sqrt_K);
    ASSIGN_CTF_PARAM(18, sqU);
    ASSIGN_CTF_PARAM(19, sqV);
    ASSIGN_CTF_PARAM(20, sqrt_angle);
    ASSIGN_CTF_PARAM(21, bgR1);
    ASSIGN_CTF_PARAM(22, bgR2);
    ASSIGN_CTF_PARAM(23, bgR3);
    ASSIGN_CTF_PARAM(24, gaussian_K);
    ASSIGN_CTF_PARAM(25, sigmaU);

    if (ia <= 26 && l > 0)
    {
        if (modelSimplification < 3)
        {
            ctfmodel.sigmaV = p[26];
            l--;
        } //     7 *
        else
        {
            ctfmodel.sigmaV = p[25];
            l--;
        }
    }
    if (ia <= 27 && l > 0)
    {
        if (modelSimplification < 3)
        {
            ctfmodel.gaussian_angle = p[27];
            l--;
        } //     8 *
        else
        {
            ctfmodel.gaussian_angle = 0;
            l--;
        }
    }

    ASSIGN_CTF_PARAM(28, cU);

    if (ia <= 29 && l > 0)
    {
        if (modelSimplification < 3)
        {
            ctfmodel.cV = p[29];
            l--;
        } //    10 *
        else
        {
            ctfmodel.cV = p[28];
            l--;
        }
    }

    ASSIGN_CTF_PARAM(30, gaussian_K2);
    ASSIGN_CTF_PARAM(31, sigmaU2);
    ASSIGN_CTF_PARAM(32, sigmaV2);
    ASSIGN_CTF_PARAM(33, gaussian_angle2);
    ASSIGN_CTF_PARAM(34, cU2);
    ASSIGN_CTF_PARAM(35, cV2);
    ASSIGN_CTF_PARAM(36, phase_shift);
    ASSIGN_CTF_PARAM(37, VPP_radius);

}//function assignCTFfromParameters


#define ASSIGN_PARAM_CTF(index, paramName) if (ia <= index && l > 0) { p[index] = ctfmodel.paramName; --l; }

void ProgCTFEstimateFromPSD::assignParametersFromCTF(CTFDescription &ctfmodel, double *p, int ia,
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
    ASSIGN_PARAM_CTF(12, envR0);
    ASSIGN_PARAM_CTF(13, envR1);
    ASSIGN_PARAM_CTF(14, envR2);
    ASSIGN_PARAM_CTF(15, Q0);
    ASSIGN_PARAM_CTF(16, base_line);
    ASSIGN_PARAM_CTF(17, sqrt_K);
    ASSIGN_PARAM_CTF(18, sqU);
    ASSIGN_PARAM_CTF(19, sqV);
    ASSIGN_PARAM_CTF(20, sqrt_angle);
    ASSIGN_PARAM_CTF(21, bgR1);
    ASSIGN_PARAM_CTF(22, bgR2);
    ASSIGN_PARAM_CTF(23, bgR3);
    ASSIGN_PARAM_CTF(24, gaussian_K);
    ASSIGN_PARAM_CTF(25, sigmaU);

    if (ia <= 26 && l > 0)
    {
        if (modelSimplification < 3)
        {
            p[26] = ctfmodel.sigmaV;
            l--;
        }
        else
        {
            p[26] = ctfmodel.sigmaU;
            l--;
        }
    }
    if (ia <= 27 && l > 0)
    {
        if (modelSimplification < 3)
        {
            p[27] = ctfmodel.gaussian_angle;
            l--;
        }
        else
        {
            p[27] = 0;
            l--;
        }
    }

    ASSIGN_PARAM_CTF(28, cU);

    if (ia <= 29 && l > 0)
    {
        if (modelSimplification < 3)
        {
            p[29] = ctfmodel.cV;
            l--;
        }
        else
        {
            p[29] = ctfmodel.cU;
            l--;
        }
    }

    ASSIGN_PARAM_CTF(30, gaussian_K2);
    ASSIGN_PARAM_CTF(31, sigmaU2);
    ASSIGN_PARAM_CTF(32, sigmaV2);
    ASSIGN_PARAM_CTF(33, gaussian_angle2);
    ASSIGN_PARAM_CTF(34, cU2);
    ASSIGN_PARAM_CTF(35, cV2);
    ASSIGN_PARAM_CTF(36, phase_shift);
    ASSIGN_PARAM_CTF(37, VPP_radius);
}

#define COPY_ctfmodel_TO_CURRENT_GUESS \
		global_prm->assignParametersFromCTF(global_prm->current_ctfmodel, \
                               MATRIX1D_ARRAY(*global_prm->adjust_params),0,ALL_CTF_PARAMETERS, \
                               global_prm->modelSimplification);

/* Read parameters --------------------------------------------------------- */
void ProgCTFEstimateFromPSD::readBasicParams(XmippProgram *program)
{
	//ProgCTFBasicParams::readBasicParams(program);

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

/* Usage ------------------------------------------------------------------- */
void ProgCTFEstimateFromPSD::defineBasicParams(XmippProgram * program)
{
	ProgCTFBasicParams::defineBasicParams(program);
	CTFDescription::defineParams(program);

}

void ProgCTFEstimateFromPSD::defineParams()
{
    defineBasicParams(this);
}

/* Produce side information ------------------------------------------------ */
void ProgCTFEstimateFromPSD::produceSideInfo()
{
    adjust.resize(ALL_CTF_PARAMETERS);
    adjust.initZeros();
    current_ctfmodel.clear();
    ctfmodel_defoci.clear();
    assignParametersFromCTF(initial_ctfmodel, MATRIX1D_ARRAY(adjust), 0, ALL_CTF_PARAMETERS, true);

    // Read the CTF file, supposed to be the uncentered squared amplitudes
    if (fn_psd != "")
        ctftomodel.read(fn_psd);
    f = &(ctftomodel());

    ProgCTFBasicParams::produceSideInfo();
    current_ctfmodel.precomputeValues(x_contfreq, y_contfreq);
}

/* Generate model so far ---------------------------------------------------- */
/* The model is taken from global_adjust and global_ctfmodel is modified */
void ProgCTFEstimateFromPSD::generateModelSoFar(Image<double> &I, bool apply_log = false)
{
    Matrix1D<int> idx(2); // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    assignCTFfromParameters(MATRIX1D_ARRAY(*adjust_params), current_ctfmodel,
                            0, ALL_CTF_PARAMETERS, modelSimplification);
    current_ctfmodel.produceSideInfo();
    I().initZeros(*f);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I())
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(*f, idx, freq);
        double w=freq.module();
        if (w>max_freq_psd)
        	continue;
        digfreq2contfreq(freq, freq, Tm);

        // Decide what to save
        current_ctfmodel.precomputeValues(XX(freq), YY(freq));
        if (action <= 1)
            I()(i, j) = current_ctfmodel.getValueNoiseAt();
        else if (action == 2)
        {
            double E = current_ctfmodel.getValueDampingAt();
            I()(i, j) = current_ctfmodel.getValueNoiseAt() + E * E;
        }
        else if (action >= 3 && action <= 6)
        {
            double ctf = current_ctfmodel.getValuePureAt();
            I()(i, j) = current_ctfmodel.getValueNoiseAt() + ctf * ctf;
        }
        else
        {
            double ctf = current_ctfmodel.getValuePureAt();
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
void ProgCTFEstimateFromPSD::saveIntermediateResults(const FileName &fn_root, bool generate_profiles = true)
{
    std::ofstream plotX, plotY, plot_radial;
    Image<double> save;
    generateModelSoFar(save, false);

    Image<double> save_ctf;
    generate_model_halfplane(ctfmodelSize,
                                         ctfmodelSize, save_ctf());
    if (fn_root.find("@")==std::string::npos)
        save_ctf.write(fn_root + "_ctfmodel_halfplane.xmp");
    else
        save_ctf.write(fn_root + "_ctfmodel_halfplane.stk");
    generate_model_quadrant(ctfmodelSize,
                                        ctfmodelSize, save_ctf());
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
    plotX << "# freq_dig freq_angstrom model psd enhanced logModel logPsd\n";
    plotY << "# freq_dig freq_angstrom model psd enhanced logModel logPsd\n";
    //plot_radial << "# freq_dig freq_angstrom model psd enhanced logModel logPsd\n";

    // Generate cut along X
    for (int i = STARTINGY(save()); i <= FINISHINGY(save()) / 2; i++)
    {
        if (mask(i, 0) <= 0)
            continue;
        plotY << w_digfreq(i, 0) << " " << w_contfreq(i, 0) << " "
        << save()(i, 0) << " " << (*f)(i, 0) << " "
        << enhanced_ctftomodel()(i, 0) << " "
		<< log10(save()(i, 0)) << " " << log10((*f)(i, 0))
		<< std::endl;
    }

    // Generate cut along Y
    for (int j = STARTINGX(save()); j <= FINISHINGX(save()) / 2; j++)
    {
        if (mask(0, j) <= 0)
            continue;
        plotX << w_digfreq(0, j) << " " << w_contfreq(0, j) << " "
        << save()(0, j) << " " << (*f)(0, j) << " "
        << enhanced_ctftomodel()(0, j)
		<< log10(save()(0, j)) << " " << log10((*f)(0, j))
		<< std::endl;
    }

    // Generate radial average
    MultidimArray<double> radial_CTFmodel_avg, radial_CTFampl_avg,
    radial_enhanced_avg;
    MultidimArray<int> radial_N;
    radial_CTFmodel_avg.initZeros(YSIZE(save()) / 2);
    radial_CTFampl_avg.initZeros(YSIZE(save()) / 2);
    radial_enhanced_avg.initZeros(YSIZE(save()) / 2);
    radial_N.initZeros(YSIZE(save()) / 2);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    {
        if (mask(i, j) <= 0)
            continue;
        double model2 = save()(i, j);

        int r = w_digfreq_r(i, j);
        radial_CTFmodel_avg(r) += model2;
        radial_CTFampl_avg(r) += (*f)(i, j);
        radial_enhanced_avg(r) += enhanced_ctftomodel()(i, j);
        radial_N(r)++;
    }

    FOR_ALL_ELEMENTS_IN_ARRAY1D (radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        plot_radial << w_digfreq(i, 0) << " " << w_contfreq(i, 0)
        << " " << radial_CTFmodel_avg(i) / radial_N(i) << " " << radial_CTFampl_avg(i) / radial_N(i) << " "
		<< radial_enhanced_avg(i) / radial_N(i)
        << " " << log10(radial_CTFmodel_avg(i) / radial_N(i)) << " " << log10(radial_CTFampl_avg(i) / radial_N(i)) << " "
		<< std::endl;
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
    enhancedPSD = enhanced_ctftomodel_fullsize();
    CenterFFT(enhancedPSD, false);
    selfScaleToSize(BSPLINE3, enhancedPSD, Ydim, Xdim);
    CenterFFT(enhancedPSD, true);

    // Generate the CTF model
    assignCTFfromParameters(MATRIX1D_ARRAY(*adjust_params), current_ctfmodel,
                            0, ALL_CTF_PARAMETERS, modelSimplification);
    current_ctfmodel.produceSideInfo();

    // Write the two model quadrants
    MultidimArray<int> mask_norm;
    STARTINGX(enhancedPSD)=STARTINGY(enhancedPSD)=0;
    mask_norm.initZeros(enhancedPSD);
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
                mask_norm(i,j)=(int)mask(i,j);
            digfreq2contfreq(freq, freq, Tm);

            current_ctfmodel.precomputeValues(XX(freq), YY(freq));
            model(i, j) = current_ctfmodel.getValuePureAt();
            model(i, j) *= model(i, j);
        }
    }
    // Normalize the left part so that it has similar values to
    // the enhanced PSD
    model.rangeAdjust(enhancedPSD, &mask_norm);
    // Copy the part of the enhancedPSD
    FOR_ALL_ELEMENTS_IN_ARRAY2D(model)
    {
		if (!((j >= Xdim / 2 && i >= Ydim / 2) || (j < Xdim / 2 && i < Ydim / 2)))
			model(i, j) = enhancedPSD(i, j);
		else
		{
            XX(idx) = j;
            YY(idx) = i;
            FFT_idx2digfreq(model, idx, freq);
            double w=freq.module();
            if (w>max_freq_psd)
            	model(i,j)=0;

		}
    }

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
    enhancedPSD = enhanced_ctftomodel_fullsize();
    CenterFFT(enhancedPSD, false);
    selfScaleToSize(BSPLINE3, enhancedPSD, Ydim, Xdim);
    CenterFFT(enhancedPSD, true);

    // The left part is the CTF model
    assignCTFfromParameters(MATRIX1D_ARRAY(*adjust_params), current_ctfmodel,
                            0, CTF_PARAMETERS, modelSimplification);
    current_ctfmodel.produceSideInfo();

    MultidimArray<int> mask_norm;
    STARTINGX(enhancedPSD)=STARTINGY(enhancedPSD)=0;
    mask_norm.initZeros(enhancedPSD);
    model.resizeNoCopy(enhancedPSD);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(model)
    {
        if (j <= Xdim / 2)
        {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(model, idx, freq);
        double w=freq.module();
        if (w>max_freq_psd)
        	continue;
        if (fabs(XX(freq))>0.03 && fabs(YY(freq))>0.03)
            mask_norm(i,j)=(int)mask(i,j);
        digfreq2contfreq(freq, freq, Tm);

        current_ctfmodel.precomputeValues(XX(freq), YY(freq));
        model(i, j) = current_ctfmodel.getValuePureAt();
        model(i, j) *= model(i, j);
        }
    }
    // Normalize the left part so that it has similar values to
    // the enhanced PSD
    model.rangeAdjust(enhancedPSD, &mask_norm);
    // Copy the part of the enhancedPSD
    FOR_ALL_ELEMENTS_IN_ARRAY2D(model)
    {
    if (!(j <= Xdim / 2))
        model(i, j) = enhancedPSD(i, j);
    else
    {
    	XX(idx) = j;
		YY(idx) = i;
		FFT_idx2digfreq(model, idx, freq);
		double w=freq.module();
		if (w>max_freq_psd)
			model(i,j)=0;
    }
    }
    CenterFFT(model, true);
}

/* CTF fitness ------------------------------------------------------------- */
/* This function measures the distance between the estimated CTF and the
 measured CTF */
double ProgCTFEstimateFromPSD::CTF_fitness_object(double *p)
{
    double retval;
    // Generate CTF model
    switch (action)
    {
        // Remind that p is a vector whose first element is at index 1
    case 0:
        assignCTFfromParameters(p - FIRST_SQRT_PARAMETER + 1,
                                current_ctfmodel, FIRST_SQRT_PARAMETER, SQRT_CTF_PARAMETERS,
                                modelSimplification);
        if (show_inf >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= SQRT_CTF_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 1:
            assignCTFfromParameters(p - FIRST_SQRT_PARAMETER + 1,
            					current_ctfmodel, FIRST_SQRT_PARAMETER,
                                    BACKGROUND_CTF_PARAMETERS, modelSimplification);
        if (show_inf >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= BACKGROUND_CTF_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 2:
            assignCTFfromParameters(p - FIRST_ENVELOPE_PARAMETER + 1,
            					current_ctfmodel, FIRST_ENVELOPE_PARAMETER, ENVELOPE_PARAMETERS,
                                    modelSimplification);
        if (show_inf >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= ENVELOPE_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 3:
            assignCTFfromParameters(p - FIRST_DEFOCUS_PARAMETER + 1,
            					current_ctfmodel, FIRST_DEFOCUS_PARAMETER, DEFOCUS_PARAMETERS,
                                    modelSimplification);
            psd_theo_radial_derivative.initZeros();
        if (show_inf >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= DEFOCUS_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }
        break;
    case 4:
            assignCTFfromParameters(p - 0 + 1, current_ctfmodel, 0,
                                    CTF_PARAMETERS, modelSimplification);
            psd_theo_radial.initZeros();
        if (show_inf >= 2)
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
		assignCTFfromParameters(p - 0 + 1, current_ctfmodel, 0,
								ALL_CTF_PARAMETERS, modelSimplification);
		psd_theo_radial.initZeros();
        if (show_inf >= 2)
        {
            std::cout << "Input vector:";
            for (int i = 1; i <= ALL_CTF_PARAMETERS; i++)
                std::cout << p[i] << " ";
            std::cout << std::endl;
        }

        break;
    }

    current_ctfmodel.produceSideInfo();

    if (show_inf >= 2)
        std::cout << "Model:\n" << current_ctfmodel << std::endl;
    if (!current_ctfmodel.hasPhysicalMeaning())
    {
        if (show_inf >= 2)
            std::cout << "Does not have physical meaning\n";
        return heavy_penalization;
    }

    if (action > 3
        && (fabs(
                (current_ctfmodel.DeltafU - ctfmodel_defoci.DeltafU)
                / ctfmodel_defoci.DeltafU) > 0.2
            || fabs(
                (current_ctfmodel.DeltafV
                 - ctfmodel_defoci.DeltafV)
                / ctfmodel_defoci.DeltafU) > 0.2))
    {
        if (show_inf >= 2)
            std::cout << "Too large defocus\n";
        return heavy_penalization;
    }
    if (initial_ctfmodel.DeltafU != 0 && action >= 3)
    {
        // If there is an initial model, the true solution
        // cannot be too far
        if (fabs(initial_ctfmodel.DeltafU - current_ctfmodel.DeltafU) > defocus_range ||
            fabs(initial_ctfmodel.DeltafV - current_ctfmodel.DeltafV) > defocus_range)
        {
            if (show_inf >= 2)
            {
                std::cout << "Too far from hint: Initial (" << initial_ctfmodel.DeltafU << "," << initial_ctfmodel.DeltafV << ")"
                << " current guess (" << current_ctfmodel.DeltafU << "," << current_ctfmodel.DeltafV << ") max allowed difference: "
                << defocus_range << std::endl;
            }
            return heavy_penalization;
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
    double lowerLimit = 1.1 * min_freq_psd;
    double upperLimit = 0.9 * max_freq_psd;
    const MultidimArray<double>& local_enhanced_ctf = enhanced_ctftomodel();
    int XdimW=XSIZE(w_digfreq);
    int YdimW=YSIZE(w_digfreq);
    corr13=0;

    for (int i = 0; i < YdimW; i += evaluation_reduction)
        for (int j = 0; j < XdimW; j += evaluation_reduction)
        {
            if (DIRECT_A2D_ELEM(mask, i, j) <= 0)
                continue;

            // Compute each component
            current_ctfmodel.precomputeValues(i, j);
            double bg = current_ctfmodel.getValueNoiseAt();
            double envelope=0, ctf_without_damping, ctf_with_damping=0;
            double ctf2_th=0;

            switch (action)
            {
            case 0:
            case 1:
                ctf2_th = bg;
                break;
            case 2:
                envelope = current_ctfmodel.getValueDampingAt();
                ctf2_th = bg + envelope * envelope;
                break;
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
                envelope = current_ctfmodel.getValueDampingAt();
                ctf_without_damping =
                		current_ctfmodel.getValuePureWithoutDampingAt();
                ctf_with_damping = envelope * ctf_without_damping;
                ctf2_th = bg + ctf_with_damping * ctf_with_damping;
                break;
            }
            // Compute distance
            double ctf2 = DIRECT_A2D_ELEM(*f, i, j);
            double dist = 0;
            double ctf_with_damping2;
            switch (action)
            {
            case 0:
            case 1:
                dist = fabs(ctf2 - bg);
                if (penalize && bg > ctf2
                    && DIRECT_A2D_ELEM(w_digfreq, i, j)
                    > max_gauss_freq)
                    dist *= current_penalty;
                break;
            case 2:
                dist = fabs(ctf2 - ctf2_th);
                if (penalize && ctf2_th < ctf2
                    && DIRECT_A2D_ELEM(w_digfreq, i, j)
                    > max_gauss_freq)
                    dist *= current_penalty;
                break;
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
                if (DIRECT_A2D_ELEM(w_digfreq,i, j) < upperLimit
                    && DIRECT_A2D_ELEM(w_digfreq,i, j) > lowerLimit)
                {
                    if  (action == 3 ||
                         (action == 4 && DIRECT_A2D_ELEM(mask_between_zeroes,i,j) == 1) ||
                         (action == 7 && DIRECT_A2D_ELEM(mask_between_zeroes,i,j) == 1))
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

                        if (action==3)
                        {
                            int r = A2D_ELEM(w_digfreq_r,i, j);
                            A1D_ELEM(psd_theo_radial,r) += ctf2_th;
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
                // If we subtract two of this scaled values
                //    x-bg      y-bg       x-y
                //   ------- - ------- = -------
                //    env^2     env^2     env^2
                break;
            }
            distsum += dist * DIRECT_A2D_ELEM(mask,i,j);
            N++;
        }

    if (N > 0)
        retval = distsum / N;
    else
        retval = heavy_penalization;

    if (show_inf >=2)
        std::cout << "Fitness1=" << retval << std::endl;

    if ( (((action >= 3) && (action <= 4)) || (action == 7))
         && (Ncorr > 0) && (enhanced_weight != 0) )
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
            || (fabs(sigma1 - sigma2) / maxSigma > 0.9 && action>=5))
        {
            retval = heavy_penalization;
            if (show_inf>=2)
                std::cout << "Fitness2=" << heavy_penalization << " sigma1=" << sigma1 << " sigma2=" << sigma2 << std::endl;
        }
        else
        {
            correlation_coeff /= sigma1 * sigma2;
            if (action == 7)
            	corr13 = correlation_coeff;
            else
                retval -= enhanced_weight * correlation_coeff;
            if (show_inf >= 2)
            {
                std::cout << "model_avg=" << model_avg << std::endl;
                std::cout << "enhanced_avg=" << enhanced_avg << std::endl;
                std::cout << "enhanced_model=" << enhanced_model / Ncorr
                << std::endl;
                std::cout << "sigma1=" << sigma1 << std::endl;
                std::cout << "sigma2=" << sigma2 << std::endl;
                std::cout << "Fitness2="
                << -(enhanced_weight * correlation_coeff)
                << " (" << correlation_coeff << ")" << std::endl;
            }
        }
        // Correlation of the derivative of the radial profile
        if (action==3 || evaluation_reduction==1)
        {
            int state=0;
            double maxDiff=0;
            psd_theo_radial_derivative.initZeros();
            double lowerlimt=1.1*min_freq;
            double upperlimit=0.9*max_freq;
            FOR_ALL_ELEMENTS_IN_ARRAY1D(psd_theo_radial)
            if (A1D_ELEM(w_digfreq_r_iN,i)>0)
            {
                A1D_ELEM(psd_theo_radial,i)*=A1D_ELEM(w_digfreq_r_iN,i);
                double freq=A2D_ELEM(w_digfreq,i,0);
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
                        double diff=A1D_ELEM(psd_theo_radial,i)-A1D_ELEM(psd_theo_radial,i-1);
                        A1D_ELEM(psd_theo_radial_derivative,i)=diff;
                        maxDiff=std::max(maxDiff,fabs(diff));
                    }
                    break;
                }
            }

            double corrRadialDerivative=0,mux=0, muy=0, Ncorr=0, sigmax=0, sigmay=0;
            double iMaxDiff=1.0/maxDiff;
            FOR_ALL_ELEMENTS_IN_ARRAY1D(psd_theo_radial)
            {
                A1D_ELEM(psd_theo_radial_derivative,i)*=iMaxDiff;
                double x=A1D_ELEM(psd_exp_enhanced_radial_derivative,i);
                double y=A1D_ELEM(psd_theo_radial_derivative,i);
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

            if (show_inf>=2)
            {
                std::cout << "Fitness3=" << -corrRadialDerivative << std::endl;
                if (show_inf==3)
                {
                	psd_exp_radial.write("PPPexpRadial.txt");
                	psd_theo_radial.write("PPPtheoRadial.txt");
                	psd_exp_radial_derivative.write("PPPexpRadialDerivative.txt");
                	psd_theo_radial_derivative.write("PPPtheoRadialDerivative.txt");
                }
            }

        }

    }

    // Show some debugging information
    if (show_inf >= 2)
    {
        std::cout << "Fitness=" << retval << std::endl;
        if (show_inf == 3)
        {
            saveIntermediateResults("PPP");
            std::cout << "Press any key\n";
            char c;
            std::cin >> c;
        }
    }

    return retval;
}

double CTF_fitness(double *p, void *vprm)
{
	ProgCTFEstimateFromPSD *prm=(ProgCTFEstimateFromPSD *) vprm;
	return prm->CTF_fitness_object(p);
}

/* Compute central region -------------------------------------------------- */
void ProgCTFEstimateFromPSD::compute_central_region(double &w1, double &w2, double ang)
{
    w1 = max_freq_psd;
    w2 = min_freq_psd;
    Matrix1D<double> freq(2), dir(2);

    // Compute first and third zero in the given direction
    VECTOR_R2(dir, COSD(ang), SIND(ang));

    // Detect first zero
    current_ctfmodel.lookFor(1, dir, freq, 0);
    if (XX(freq) == -1 && YY(freq) == -1)
        w1 = min_freq_psd;
    else
    {
        contfreq2digfreq(freq, freq, Tm);
        double w;
        if (XX(dir) > 0.1)
            w = XX(freq) / XX(dir);
        else
            w = YY(freq) / YY(dir);
        w1 = XMIPP_MAX(min_freq_psd, XMIPP_MIN(w1, w));
    }

    // Detect fifth zero
    current_ctfmodel.lookFor(5, dir, freq, 0);
    if (XX(freq) == -1 && YY(freq) == -1)
        w2 = max_freq_psd;
    else
    {
        double w;
        contfreq2digfreq(freq, freq, Tm);
        if (XX(dir) > 0.1)
            w = XX(freq) / XX(dir);
        else
            w = YY(freq) / YY(dir);
        w2 = XMIPP_MIN(max_freq_psd, XMIPP_MAX(w2, w));
    }
}

/* Center focus ----------------------------------------------------------- */
void ProgCTFEstimateFromPSD::center_optimization_focus(bool adjust_freq, bool adjust_th, double margin = 1)
{
    if (show_optimization)
        std::cout << "Freq frame before focusing=" << min_freq_psd << ","
        << max_freq_psd << std::endl << "Value_th before focusing="
        << value_th << std::endl;

    double w1 = min_freq_psd, w2 = max_freq_psd;
    if (adjust_freq)
    {
        double w1U, w2U, w1V, w2V;
        compute_central_region(w1U, w2U, current_ctfmodel.azimuthal_angle);
        compute_central_region(w1V, w2V, current_ctfmodel.azimuthal_angle + 90);
        w1 = XMIPP_MIN(w1U, w1V);
        w2 = XMIPP_MAX(w2U, w2V);
        min_freq_psd = XMIPP_MAX(min_freq_psd, w1 - 0.05);
        max_freq_psd = XMIPP_MIN(max_freq_psd, w2 + 0.01);
    }

    // Compute maximum value within central region
    if (adjust_th)
    {
        Image<double> save;
        generateModelSoFar(save);
        double max_val = 0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
        {
            double w = w_digfreq(i, j);
            if (w >= w1 && w <= w2)
                max_val = XMIPP_MAX(max_val, save()(i, j));
        }
        if (value_th != -1)
        	value_th = XMIPP_MIN(value_th, max_val * margin);
        else
        	value_th = max_val * margin;
    }

    if (show_optimization)
        std::cout << "Freq frame after focusing=" << min_freq_psd << ","
        << max_freq_psd << std::endl << "Value_th after focusing="
        << value_th << std::endl;
}

// Estimate sqrt parameters ------------------------------------------------
// Results are written in global_ctfmodel
void ProgCTFEstimateFromPSD::estimate_background_sqrt_parameters()
{
    if (show_optimization)
        std::cout << "Computing first sqrt background ...\n";

    // Estimate the base line taking the value of the CTF
    // for the maximum X and Y frequencies
    double base_line = 0;
    int N = 0;

    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    if (w_digfreq(i, j) > 0.4)
    {
        N++;
        base_line += (*f)(i, j);

    }
    current_ctfmodel.base_line = base_line / N;
    // Find the linear least squares solution for the sqrt part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    {
        if (mask(i, j) <= 0)
            continue;

        // Compute weight for this point
        double weight = 1 + max_freq_psd - w_digfreq(i, j);

        // Compute error
        current_ctfmodel.precomputeValues(x_contfreq(i, j),y_contfreq(i, j));
        double explained = current_ctfmodel.getValueNoiseAt();
        double unexplained = (*f)(i, j) - explained;
        if (unexplained <= 0)
            continue;
        unexplained = log(unexplained);

        double X = -sqrt(w_contfreq(i, j));
        A(0, 0) += weight * X * X;
        A(0, 1) += weight * X;
        A(1, 1) += weight * 1;
        b(0) += X * weight * unexplained;
        b(1) += weight * unexplained;
    }
    A(1, 0) = A(0, 1);

    b = A.inv() * b;

    current_ctfmodel.sqU = current_ctfmodel.sqV = b(0);
    current_ctfmodel.sqrt_K = exp(b(1));
    current_ctfmodel.sqrt_angle = 0;

    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (show_optimization)
    {
        std::cout << "First SQRT Fit:\n" << current_ctfmodel << std::endl;
        saveIntermediateResults("step01a_first_sqrt_fit");
    }

    // Now optimize .........................................................
    double fitness;
    Matrix1D<double> steps;
    steps.resize(SQRT_CTF_PARAMETERS);
    steps.initConstant(1);
    // Optimize without penalization
    if (show_optimization)
        std::cout << "Looking for best fitting sqrt ...\n";
    penalize = false;
    int iter;
    powellOptimizer(*adjust_params, FIRST_SQRT_PARAMETER + 1,
                    SQRT_CTF_PARAMETERS, CTF_fitness, global_prm, 0.05, fitness, iter, steps,
                    show_optimization);
    // Optimize with penalization
    if (show_optimization)
        std::cout << "Penalizing best fitting sqrt ...\n";
    penalize = true;
    current_penalty = 2;
    int imax = CEIL(log(penalty) / log(2.0));
    for (int i = 1; i <= imax; i++)
    {
        if (show_optimization)
            std::cout << "     Iteration " << i << " penalty="
            << current_penalty << std::endl;
        powellOptimizer(*adjust_params, FIRST_SQRT_PARAMETER + 1,
                        SQRT_CTF_PARAMETERS, CTF_fitness, global_prm, 0.05, fitness, iter,
                        steps, show_optimization);
        current_penalty *= 2;
        current_penalty =
            XMIPP_MIN(current_penalty, penalty);
    }
    // Keep the result in global_prm->adjust
    current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (show_optimization)
    {
        std::cout << "Best penalized SQRT Fit:\n" << current_ctfmodel
        << std::endl;
        saveIntermediateResults("step01b_best_penalized_sqrt_fit");
    }

    center_optimization_focus(false, true, 1.5);
}

// Estimate gaussian parameters --------------------------------------------
//#define DEBUG
void ProgCTFEstimateFromPSD::estimate_background_gauss_parameters()
{

    if (show_optimization)
        std::cout << "Computing first background Gaussian parameters ...\n";

    // Compute radial averages
    MultidimArray<double> radial_CTFmodel_avg(YSIZE(*f) / 2);
    MultidimArray<double> radial_CTFampl_avg(YSIZE(*f) / 2);
    MultidimArray<int> radial_N(YSIZE(*f) / 2);
    double w_max_gauss = 0.25;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    {
        if (mask(i, j) <= 0)
            continue;
        double w = w_digfreq(i, j);
        if (w > w_max_gauss)
            continue;

        int r = FLOOR(w * (double)YSIZE(*f));
        current_ctfmodel.precomputeValues(x_contfreq(i, j),y_contfreq(i, j));
        radial_CTFmodel_avg(r) += current_ctfmodel.getValueNoiseAt();
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
    double error2_min = 0, wmin=0;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        double w = w_digfreq(i, 0);

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
    max_gauss_freq = wmin;
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
        double w = w_digfreq(i, 0);
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
    fmax = current_ctfmodel.cV = current_ctfmodel.cU = wmax / Tm;
#ifdef DEBUG

    std::cout << "Freq of the maximum error: " << wmax << " " << fmax << std::endl;
#endif

    // Find the linear least squares solution for the gauss part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    {
        if (mask(i, j) <= 0)
            continue;
        if (w_digfreq(i, j) > wmin)
            continue;
        double fmod = w_contfreq(i, j);

        // Compute weight for this point
        double weight = 1 + max_freq_psd - w_digfreq(i, j);
        // Compute error
        current_ctfmodel.precomputeValues(x_contfreq(i, j),y_contfreq(i, j));
        double explained = current_ctfmodel.getValueNoiseAt();

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
        current_ctfmodel.sigmaU = XMIPP_MIN(fabs(b(1)), 95e3); // This value should be
        current_ctfmodel.sigmaV = XMIPP_MIN(fabs(b(1)), 95e3); // conformant with the physical
        // meaning routine in CTF.cc
        current_ctfmodel.gaussian_K = exp(b(0));
        // Store the CTF values in global_prm->adjust
        current_ctfmodel.forcePhysicalMeaning();
        COPY_ctfmodel_TO_CURRENT_GUESS;

        if (show_optimization)
        {
            std::cout << "First Background Fit:\n" << current_ctfmodel << std::endl;
            saveIntermediateResults("step01c_first_background_fit");
        }
        center_optimization_focus(false, true, 1.5);
    }
}
#undef DEBUG

// Estimate second gaussian parameters -------------------------------------
//#define DEBUG
void ProgCTFEstimateFromPSD::estimate_background_gauss_parameters2()
{
    if (show_optimization)
        std::cout << "Computing first background Gaussian2 parameters ...\n";

    // Compute radial averages
    MultidimArray<double> radial_CTFmodel_avg(YSIZE(*f) / 2);
    MultidimArray<double> radial_CTFampl_avg(YSIZE(*f) / 2);
    MultidimArray<int> radial_N(YSIZE(*f) / 2);
    double w_max_gauss = 0.25;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    {
        if (mask(i, j) <= 0)
            continue;
        double w = w_digfreq(i, j);
        if (w > w_max_gauss)
            continue;

        int r = FLOOR(w * (double)YSIZE(*f));
        double f_x = DIRECT_A2D_ELEM(x_contfreq, i, j);
        double f_y = DIRECT_A2D_ELEM(y_contfreq, i, j);
        current_ctfmodel.precomputeValues(f_x, f_y);
        double bg = current_ctfmodel.getValueNoiseAt();
        double envelope = current_ctfmodel.getValueDampingAt();
        double ctf_without_damping =
        		current_ctfmodel.getValuePureWithoutDampingAt();
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
    double error_max = 0, wmax=0, fmax;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_CTFmodel_avg)
    {
        if (radial_N(i) == 0)
            continue;
        double w = w_digfreq(i, 0);
        if (w > wmin)
            break;
        if (error(i) < error_max)
        {
            wmax = w;
            error_max = error(i);
        }
    }
    fmax = current_ctfmodel.cV2 = current_ctfmodel.cU2 = wmax / Tm;
#ifdef DEBUG

    std::cout << "Freq of the maximum error: " << wmax << " " << fmax << std::endl;
#endif

    // Find the linear least squares solution for the gauss part
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();
    int N = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    {
        if (mask(i, j) <= 0)
            continue;
        if (w_digfreq(i, j) > wmin)
            continue;
        double fmod = w_contfreq(i, j);

        // Compute the zero on the direction of this point
        Matrix1D<double> u(2), fzero(2);
        XX(u) = x_contfreq(i, j) / fmod;
        YY(u) = y_contfreq(i, j) / fmod;
        current_ctfmodel.lookFor(1, u, fzero, 0);
        if (fmod > fzero.module())
            continue;

        // Compute weight for this point
        double weight = 1 + max_freq_psd - w_digfreq(i, j);

        // Compute error
        double f_x = DIRECT_A2D_ELEM(x_contfreq, i, j);
        double f_y = DIRECT_A2D_ELEM(y_contfreq, i, j);
        current_ctfmodel.precomputeValues(f_x, f_y);
        double bg = current_ctfmodel.getValueNoiseAt();
        double envelope = current_ctfmodel.getValueDampingAt();
        double ctf_without_damping =
        		current_ctfmodel.getValuePureWithoutDampingAt();
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
            current_ctfmodel.sigmaU2 = XMIPP_MIN(fabs(b(1)), 95e3); // This value should be
            current_ctfmodel.sigmaV2 = XMIPP_MIN(fabs(b(1)), 95e3); // conformant with the physical
            // meaning routine in CTF.cc
            current_ctfmodel.gaussian_K2 = exp(b(0));
        }
        else
        {
        	current_ctfmodel.sigmaU2 = current_ctfmodel.sigmaV2 = 0;
        	current_ctfmodel.gaussian_K2 = 0;
        }
    }
    else
    {
    	current_ctfmodel.sigmaU2 = current_ctfmodel.sigmaV2 = 0;
    	current_ctfmodel.gaussian_K2 = 0;
    }

    // Store the CTF values in global_prm->adjust
    current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

#ifdef DEBUG
    // Check
    FOR_ALL_ELEMENTS_IN_ARRAY2D(w_digfreq)
    {
        if (mask(i, j)<=0)
            continue;
        if (w_digfreq(i, j) > wmin)
            continue;
        double fmod = w_contfreq(i, j);

        // Compute the zero on the direction of this point
        Matrix1D<double> u(2), fzero(2);
        XX(u) = x_contfreq(i, j) / fmod;
        YY(u) = y_contfreq(i, j) / fmod;
        global_ctfmodel.zero(1, u, fzero);
        if (fmod > fzero.module())
            continue;

        // Compute error
        double f_x = DIRECT_A2D_ELEM(x_contfreq, i, j);
        double f_y = DIRECT_A2D_ELEM(y_contfreq, i, j);
        double bg = current_ctfmodel.getValueNoiseAt(f_x, f_y);
        double envelope = current_ctfmodel.getValueDampingAt(f_x, f_y);
        double ctf_without_damping = current_ctfmodel.getValuePureWithoutDampingAt(f_x, f_y);
        double ctf_with_damping = envelope * ctf_without_damping;
        double ctf2_th = bg + ctf_with_damping * ctf_with_damping;
        double explained = ctf2_th;
        double unexplained = explained - (*f)(i, j);

        if (unexplained <= 0)
            continue;
        std::cout << fmod << " " << unexplained << " "
        << current_ctfmodel.gaussian_K2*exp(-current_ctfmodel.sigmaU2*
                                           (fmod - fmax)*(fmod - fmax)) << std::endl;
    }
#endif

    if (show_optimization)
    {
        std::cout << "First Background Gaussian 2 Fit:\n" << current_ctfmodel
        << std::endl;
        saveIntermediateResults("step04a_first_background2_fit");
    }
}
#undef DEBUG

// Estimate envelope parameters --------------------------------------------
//#define DEBUG
void ProgCTFEstimateFromPSD::estimate_envelope_parameters()
{
    if (show_optimization)
        std::cout << "Looking for best fitting envelope ...\n";

    // Set the envelope
    current_ctfmodel.Ca = initial_ctfmodel.Ca;
    current_ctfmodel.K = 1.0;
    current_ctfmodel.espr = 0.0;
    current_ctfmodel.ispr = 0.0;
    current_ctfmodel.alpha = 0.0;
    current_ctfmodel.DeltaF = 0.0;
    current_ctfmodel.DeltaR = 0.0;
    current_ctfmodel.Q0 = initial_ctfmodel.Q0;
    current_ctfmodel.envR0 = 0.0;
    current_ctfmodel.envR1 = 0.0;
    COPY_ctfmodel_TO_CURRENT_GUESS;

    // Now optimize the envelope
    penalize = false;
    int iter;
    double fitness;
    Matrix1D<double> steps;
    steps.resize(ENVELOPE_PARAMETERS);
    steps.initConstant(1);
    steps(1) = 0; // Do not optimize Cs
    steps(5) = 0; // Do not optimize for alpha, since Ealpha depends on the
    // defocus
    if (modelSimplification >= 1)
        steps(6) = steps(7) = 0; // Do not optimize DeltaF and DeltaR
    powellOptimizer(*adjust_params, FIRST_ENVELOPE_PARAMETER + 1,
                    ENVELOPE_PARAMETERS, CTF_fitness, global_prm, 0.05, fitness, iter, steps,
                    show_optimization);

    // Keep the result in global_prm->adjust
    current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (show_optimization)
    {
        std::cout << "Best envelope Fit:\n" << current_ctfmodel << std::endl;
        saveIntermediateResults("step02a_best_envelope_fit");
    }

    // Optimize with penalization
    if (show_optimization)
        std::cout << "Penalizing best fitting envelope ...\n";
    penalize = true;
    current_penalty = 2;
    int imax = CEIL(log(penalty) / log(2.0));
    for (int i = 1; i <= imax; i++)
    {
        if (show_optimization)
            std::cout << "     Iteration " << i << " penalty="
            << current_penalty << std::endl;
        powellOptimizer(*adjust_params, FIRST_ENVELOPE_PARAMETER + 1,
                        ENVELOPE_PARAMETERS, CTF_fitness, global_prm, 0.05, fitness, iter,
                        steps, show_optimization);
        current_penalty *= 2;
        current_penalty =
            XMIPP_MIN(current_penalty,penalty);
    }
    // Keep the result in global_prm->adjust
    current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (show_optimization)
    {
        std::cout << "Best envelope Fit:\n" << current_ctfmodel << std::endl;
        saveIntermediateResults("step02b_best_penalized_envelope_fit");
    }
}
#undef DEBUG

// Estimate defoci ---------------------------------------------------------
void ProgCTFEstimateFromPSD::showFirstDefoci()
{
    if (show_optimization)
    {
        std::cout << "First defocus Fit:\n" << current_ctfmodel << std::endl;
        saveIntermediateResults("step03a_first_defocus_fit");
        enhanced_ctftomodel.write("step03a_enhanced_PSD.xmp");
        Image<double> save, save2, save3;
        save().resize(YSIZE(w_digfreq), XSIZE(w_digfreq));
        save2().resize(save());
        save3().resize(save());
        FOR_ALL_ELEMENTS_IN_ARRAY2D(save())
        {
            save()(i, j) = enhanced_ctftomodel()(i, j);
            double f_x = DIRECT_A2D_ELEM(x_contfreq, i, j);
            double f_y = DIRECT_A2D_ELEM(y_contfreq, i, j);
            current_ctfmodel.precomputeValues(f_x, f_y);
            double ctf_without_damping =
            		current_ctfmodel.getValuePureWithoutDampingAt();
            save2()(i, j) = ctf_without_damping * ctf_without_damping;
            save3()(i, j) = -enhanced_ctftomodel()(i, j)
                            * ctf_without_damping * ctf_without_damping;
        }
        save.write("step03a_enhanced_PSD.xmp");
        save2.write("step03a_fitted_CTF.xmp");
        save3.write("step03a_superposition.xmp");
    }
}

//#define DEBUG
void ProgCTFEstimateFromPSD::estimate_defoci()
{
    if (show_optimization)
        std::cout << "Looking for first defoci ...\n";
    double best_defocusU=0, best_defocusV=0, best_angle=0, best_K=1;
    double best_error = heavy_penalization * 1.1;
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
    if (initial_ctfmodel.DeltafU != 0)
    {
        initial_defocusStep = std::min(defocus_range,20000.0);
        defocusU0 = std::max(
                        1e3,
                        initial_ctfmodel.DeltafU
                        - defocus_range);
        double maxDeviation = std::max(defocus_range,
                                       0.25 * initial_ctfmodel.DeltafU);
        max_allowed_defocusU = std::min(100e3,
                                        initial_ctfmodel.DeltafU + maxDeviation);
        defocusUF = std::min(
                        150e3,
                        initial_ctfmodel.DeltafU
                        + defocus_range);
        min_allowed_defocusU = std::max(1e3,
                                        initial_ctfmodel.DeltafU - maxDeviation);
        if (initial_ctfmodel.DeltafV == 0)
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
                            initial_ctfmodel.DeltafV
                            - defocus_range);
            max_allowed_defocusV = std::max(100e3,
                                            initial_ctfmodel.DeltafV + maxDeviation);
            defocusVF = std::min(
                            150e3,
                            initial_ctfmodel.DeltafV
                            + defocus_range);
            min_allowed_defocusV = std::max(1e3,
                                            initial_ctfmodel.DeltafV - maxDeviation);
        }
    }

    double K_so_far = current_ctfmodel.K;
    Matrix1D<double> steps(DEFOCUS_PARAMETERS);
    steps.initConstant(1);
    steps(3) = 0; // Do not optimize kV
    steps(4) = 0; // Do not optimize K
    for (double defocusStep = initial_defocusStep;
         defocusStep >= std::min(5000., defocus_range / 2);
         defocusStep /= 2)
    {
        error.resize(CEIL((defocusVF - defocusV0) / defocusStep + 1),
                     CEIL((defocusUF - defocusU0) / defocusStep + 1));
        error.initConstant(heavy_penalization);
        if (show_optimization)
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
                    error(i, j) = heavy_penalization;
                    continue;
                }
                for (double angle = 0; angle < 180; angle += 45)
                {
                    int iter;
                    double fitness;

                    (*adjust_params)(0) = defocusU;
                    (*adjust_params)(1) = defocusV;
                    (*adjust_params)(2) = angle;
                    (*adjust_params)(4) = K_so_far;

                    powellOptimizer(*adjust_params, FIRST_DEFOCUS_PARAMETER + 1,
                                    DEFOCUS_PARAMETERS, CTF_fitness, global_prm, 0.05,
                                    fitness, iter, steps, false);

                    if ((first_angle || fitness < error(i, j))
                        && (current_ctfmodel.DeltafU >= min_allowed_defocusU
                            && current_ctfmodel.DeltafU
                            <= max_allowed_defocusU
                            && current_ctfmodel.DeltafV
                            >= min_allowed_defocusV
                            && current_ctfmodel.DeltafV
                            <= max_allowed_defocusV))
                    {
                        error(i, j) = fitness;
                        first_angle = false;
                        if (error(i, j) < best_error || first)
                        {
                            best_error = error(i, j);
                            best_defocusU = current_ctfmodel.DeltafU;
                            best_defocusV = current_ctfmodel.DeltafV;
                            best_angle = current_ctfmodel.azimuthal_angle;
                            best_K = current_ctfmodel.K;
                            first = false;
                            if (show_optimization)
                            {
                                std::cout << "    (DefocusU,DefocusV)=("
                                << defocusU << "," << defocusV
                                << "), ang=" << angle << " --> ("
                                << current_ctfmodel.DeltafU << ","
                                << current_ctfmodel.DeltafV << "),"
                                << current_ctfmodel.azimuthal_angle
                                << " K=" << current_ctfmodel.K
                                << " error=" << error(i, j)
                                << std::endl;
#ifdef DEBUG

                                show_inf=3;
                                CTF_fitness(adjust_params->vdata-1,NULL);
                                show_inf=0;

                                Image<double> save;
                                save() = enhanced_ctftomodel();
                                save.write("PPPenhanced.xmp");
                                for (int i = 0; i < YSIZE(w_digfreq); i += 1)
                                    for (int j = 0; j < XSIZE(w_digfreq); j += 1)
                                    {
                                        if (DIRECT_A2D_ELEM(mask, i, j)<=0)
                                            continue;
                                        double f_x = DIRECT_A2D_ELEM(x_contfreq, i, j);
                                        double f_y = DIRECT_A2D_ELEM(y_contfreq, i, j);
                                        double envelope = current_ctfmodel.getValueDampingAt(f_x, f_y);
                                        double ctf_without_damping = current_ctfmodel.getValuePureWithoutDampingAt(f_x, f_y);
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
                if (error(ii, jj) != heavy_penalization)
                {
                	aValidErrorHasBeenFound=true;
                    if (error(ii, jj) < errmin)
                        errmin = error(ii, jj);
                    else if (errmax == heavy_penalization)
                        errmax = error(ii, jj);
                    else if (error(ii, jj) > errmax)
                        errmax = error(ii, jj);
                }
            }
        if (show_optimization)
            std::cout << "Error matrix\n" << error << std::endl;
        if (!aValidErrorHasBeenFound)
        	REPORT_ERROR(ERR_NUMERICAL,"Cannot find any good defocus within the given range");

        // Find those defoci which are within a 10% of the maximum
        if (show_inf >= 2)
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
                if (show_inf >= 2)
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
        if (show_inf >= 2)
        {
            Image<double> save;
            save() = error;
            save.write("error.xmp");
            std::cout << "Press any key: Error saved\n";
            char c;
            std::cin >> c;
        }
    }

    current_ctfmodel.DeltafU = best_defocusU;
    current_ctfmodel.DeltafV = best_defocusV;
    current_ctfmodel.azimuthal_angle = best_angle;
    current_ctfmodel.K = best_K;

    // Keep the result in global_prm->adjust
    current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;
    ctfmodel_defoci = current_ctfmodel;

    showFirstDefoci();
}
#undef DEBUG

// Estimate defoci with Zernike and SPTH transform---------------------------------------------
void ProgCTFEstimateFromPSD::estimate_defoci_Zernike(MultidimArray<double> &psdToModelFullSize, double min_freq, double max_freq, double Tm,
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
    double K_so_far = current_ctfmodel.K;
    double lambda=12.2643247/std::sqrt(kV*(1.+0.978466e-6*kV));
    double Z8;
    double Z3;
    double Z4;
    double Z5;
    double eAngle=0.;
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

    Matrix1D<double> initialGlobalAdjust = (*adjust_params);

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

            (*adjust_params)(0) = deFocusAvg+deFocusDiff;
            (*adjust_params)(1) = deFocusAvg-deFocusDiff;
            (*adjust_params)(2) = eAngle;
            (*adjust_params)(4) = K_so_far;
            (*adjust_params)(6) = 2;


            fitness =0;
            powellOptimizer(*adjust_params, FIRST_DEFOCUS_PARAMETER + 1,
                            DEFOCUS_PARAMETERS, CTF_fitness, global_prm, 0.05,
                            fitness, iter, steps, false);

            VEC_ELEM(arrayDefocusAvg,i)  = ((*adjust_params)(0) +(*adjust_params)(1))/2;
            VEC_ELEM(arrayDefocusDiff,i) = ((*adjust_params)(0) -(*adjust_params)(1))/2;
            VEC_ELEM(arrayError,i) = (-1)*fitness;

        }
    }

    int maxInd=arrayError.maxIndex();

    while ( (VEC_ELEM(arrayDefocusAvg,maxInd) < 3000) || ((VEC_ELEM(arrayDefocusAvg,maxInd) > 50000) && VEC_ELEM(arrayError,maxInd)>-1e3 ))
    {
        VEC_ELEM(arrayError,maxInd) = -1e3;
        VEC_ELEM(arrayDefocusAvg,maxInd) = initial_ctfmodel.DeltafU;
        VEC_ELEM(arrayDefocusDiff,maxInd) = initial_ctfmodel.DeltafV;
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
    (*adjust_params)(0) = VEC_ELEM(arrayDefocusAvg,maxInd)+VEC_ELEM(arrayDefocusDiff,maxInd);
    (*adjust_params)(1) = VEC_ELEM(arrayDefocusAvg,maxInd)-VEC_ELEM(arrayDefocusDiff,maxInd);
    (*adjust_params)(2) = eAngle;
    (*adjust_params)(4) = K_so_far;
    (*adjust_params)(6) = 2;

    fitness =0;
    powellOptimizer(*adjust_params, FIRST_DEFOCUS_PARAMETER + 1,
                    DEFOCUS_PARAMETERS, CTF_fitness, global_prm, 0.05,
                    fitness, iter, steps, false);

    VEC_ELEM(arrayDefocusU,0) = (*adjust_params)(0);
    VEC_ELEM(arrayDefocusV,0) = (*adjust_params)(1);
    VEC_ELEM(arrayError2,0) = (-1)*fitness;

    // We optimize for deltaU, deltaU
    (*adjust_params)(0) = VEC_ELEM(arrayDefocusAvg,maxInd)+VEC_ELEM(arrayDefocusDiff,maxInd);
    (*adjust_params)(1) = VEC_ELEM(arrayDefocusAvg,maxInd)+VEC_ELEM(arrayDefocusDiff,maxInd);
    (*adjust_params)(2) = eAngle;
    (*adjust_params)(4) = K_so_far;
    (*adjust_params)(6) = 2;

    fitness =0;
    powellOptimizer(*adjust_params, FIRST_DEFOCUS_PARAMETER + 1,
                    DEFOCUS_PARAMETERS, CTF_fitness, global_prm, 0.05,
                    fitness, iter, steps, false);

    VEC_ELEM(arrayDefocusU,1) = (*adjust_params)(0);
    VEC_ELEM(arrayDefocusV,1) = (*adjust_params)(1);
    VEC_ELEM(arrayError2,1) = (-1)*fitness;

    // We optimize for deltaV, deltaV
    (*adjust_params)(0) = VEC_ELEM(arrayDefocusAvg,maxInd)-VEC_ELEM(arrayDefocusDiff,maxInd);
    (*adjust_params)(1) = VEC_ELEM(arrayDefocusAvg,maxInd)-VEC_ELEM(arrayDefocusDiff,maxInd);
    (*adjust_params)(2) = eAngle;
    (*adjust_params)(4) = K_so_far;
    (*adjust_params)(6) = 2;

    fitness =0;
    powellOptimizer(*adjust_params, FIRST_DEFOCUS_PARAMETER + 1,
                    DEFOCUS_PARAMETERS, CTF_fitness, global_prm, 0.05,
                    fitness, iter, steps, false);

    VEC_ELEM(arrayDefocusU,2) = (*adjust_params)(0);
    VEC_ELEM(arrayDefocusV,2) = (*adjust_params)(1);
    VEC_ELEM(arrayError2,2) = (-1)*fitness;

    //Here we select the best one
    maxInd=arrayError2.maxIndex();
    defocusU = VEC_ELEM(arrayDefocusU,maxInd);
    defocusV = VEC_ELEM(arrayDefocusV,maxInd);

    (*adjust_params)(0) = defocusU;
    (*adjust_params)(1) = defocusV;
    (*adjust_params)(2) = eAngle;
    (*adjust_params)(4) = K_so_far;
    (*adjust_params)(6) = 2;

    while ( (0.5*(defocusU+defocusV) < 2500) || (0.5*(defocusU+defocusV) > 60000) )
    {
        VEC_ELEM(arrayError2,maxInd) = -1e3;
        VEC_ELEM(arrayDefocusU,maxInd) = initial_ctfmodel.DeltafU;
        VEC_ELEM(arrayDefocusV,maxInd) = initial_ctfmodel.DeltafV;

        maxInd=arrayError2.maxIndex();
        defocusU = VEC_ELEM(arrayDefocusU,maxInd);
        defocusV = VEC_ELEM(arrayDefocusV,maxInd);
        (*adjust_params)(0) = defocusU;
        (*adjust_params)(1) = defocusV;
        (*adjust_params)(2) = eAngle;
        (*adjust_params)(4) = K_so_far;
        (*adjust_params)(6) = 2;
    }

    if (VEC_ELEM(arrayError2,maxInd) <= 0)
    {
        COPY_ctfmodel_TO_CURRENT_GUESS;
        ctfmodel_defoci = current_ctfmodel;

        action = 5;

        steps.resize(ALL_CTF_PARAMETERS);
        steps.initConstant(1);
        steps(3) = 0; // kV
        steps(5) = 0; // The spherical aberration (Cs) is not optimized
        if (current_ctfmodel.Q0!=0)
        	steps(12)=0;

        COPY_ctfmodel_TO_CURRENT_GUESS;

        evaluation_reduction = 2;
        powellOptimizer(*adjust_params, 0 + 1, ALL_CTF_PARAMETERS, CTF_fitness,
                        global_prm, 0.01, fitness, iter, steps,show_optimization);
        COPY_ctfmodel_TO_CURRENT_GUESS;

        show_inf=0;
        action = 3;
        evaluation_reduction = 1;

        double error = -CTF_fitness(adjust_params->vdata-1,NULL);
        if ( error <= -0.1)
        {
            *adjust_params = initialGlobalAdjust;
            COPY_ctfmodel_TO_CURRENT_GUESS;
            //There is nothing to do and we have to perform an exhaustive search
#ifndef RELEASE_MODE
            std::cout << " Entering in estimate_defoci, Performing exhaustive defocus search (SLOW)" << std::endl;
#endif
            estimate_defoci();
        }
    }
}

void ProgCTFEstimateFromPSD::estimate_defoci_Zernike()
{
    if (show_optimization)
        std::cout << "Looking for first defoci ...\n";

    DEBUG_TEXTFILE("Step 6.1");
    DEBUG_MODEL_TEXTFILE;
    estimate_defoci_Zernike(enhanced_ctftomodel_fullsize(),
                            min_freq,max_freq,Tm,
                            initial_ctfmodel.kV,
                            lambdaPhase,sizeWindowPhase,
							current_ctfmodel.DeltafU, current_ctfmodel.DeltafV, current_ctfmodel.azimuthal_angle, 0);
    DEBUG_TEXTFILE("Step 6.2");
    DEBUG_MODEL_TEXTFILE;

    current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;
    ctfmodel_defoci = current_ctfmodel;

    showFirstDefoci();
}

/* Main routine ------------------------------------------------------------ */
//#define DEBUG
double ROUT_Adjust_CTF(ProgCTFEstimateFromPSD &prm, CTFDescription &output_ctfmodel, bool standalone)
{
    DEBUG_OPEN_TEXTFILE(prm.fn_psd.removeLastExtension());
    global_prm = &prm;
    if (standalone || prm.show_optimization)
        prm.show();
    prm.produceSideInfo();
    DEBUG_TEXTFILE(formatString("After producing side info: Avg=%f",prm.ctftomodel().computeAvg()));
    DEBUG_MODEL_TEXTFILE;

    // Build initial frequency mask
    prm.value_th = -1;
    prm.min_freq_psd = prm.min_freq;
    prm.max_freq_psd = prm.max_freq;

    // Set some global variables
    prm.adjust_params = &prm.adjust;
    prm.penalize = false;
    prm.max_gauss_freq = 0;
    prm.heavy_penalization = prm.f->computeMax() * XSIZE(*prm.f) * YSIZE(*prm.f);
    prm.show_inf = 0;

    // Some variables needed by all steps
    int iter;
    double fitness;
    Matrix1D<double> steps;

    /************************************************************************/
    /* STEPs 1, 2, 3 and 4:  Find background which best fits the CTF        */
    /************************************************************************/

    prm.current_ctfmodel.enable_CTFnoise = true;
    prm.current_ctfmodel.enable_CTF = false;
    prm.evaluation_reduction = 4;

    // If initial parameters were not supplied for the gaussian curve,
    // estimate them from the CTF file
    prm.action = 0;
    if (prm.adjust(FIRST_SQRT_PARAMETER) == 0)
    {
        prm.estimate_background_sqrt_parameters();
        prm.estimate_background_gauss_parameters();
    }

    // Optimize the current background
    prm.action = 1;
    prm.penalize = true;
    prm.current_penalty = prm.penalty;
    steps.resize(BACKGROUND_CTF_PARAMETERS);
    steps.initConstant(1);
    if (!prm.modelSimplification >= 3)
        steps(7) = steps(8) = steps(10) = 0;
    powellOptimizer(*prm.adjust_params, FIRST_SQRT_PARAMETER + 1,
                    BACKGROUND_CTF_PARAMETERS, CTF_fitness, global_prm, 0.01, fitness, iter,
                    steps, prm.show_optimization);

    // Make sure that the model has physical meaning
    // (In some machines due to numerical imprecission this check is necessary
    // at the end)
    prm.current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (prm.show_optimization)
    {
        std::cout << "Best background Fit:\n" << prm.current_ctfmodel << std::endl;
        prm.saveIntermediateResults("step01d_best_background_fit");
    }
    DEBUG_TEXTFILE(formatString("Step 4: CTF_fitness=%f",CTF_fitness));
    DEBUG_MODEL_TEXTFILE;

    /************************************************************************/
    /* STEPs 5 and 6:  Find envelope which best fits the CTF                */
    /************************************************************************/
    prm.action = 2;
    prm.current_ctfmodel.enable_CTF = true;
    prm.current_ctfmodel.phase_shift = prm.initial_ctfmodel.phase_shift;
    prm.current_ctfmodel.VPP_radius = prm.initial_ctfmodel.VPP_radius;
    if (prm.initial_ctfmodel.K == 0)
    {
    	prm.current_ctfmodel.kV = prm.initial_ctfmodel.kV;
    	prm.current_ctfmodel.Cs = prm.initial_ctfmodel.Cs;
        if (prm.initial_ctfmodel.Q0 != 0)
        	prm.current_ctfmodel.Q0 = prm.initial_ctfmodel.Q0;
        prm.estimate_envelope_parameters();
    }
    else
    {
    	prm.current_ctfmodel.K = prm.initial_ctfmodel.K;
    	prm.current_ctfmodel.kV = prm.initial_ctfmodel.kV;
    	prm.current_ctfmodel.DeltafU = prm.initial_ctfmodel.DeltafU;
    	prm.current_ctfmodel.DeltafV = prm.initial_ctfmodel.DeltafV;
    	prm.current_ctfmodel.azimuthal_angle = prm.initial_ctfmodel.azimuthal_angle;
    	prm.current_ctfmodel.Cs = prm.initial_ctfmodel.Cs;
    	prm.current_ctfmodel.Ca = prm.initial_ctfmodel.Ca;
    	prm.current_ctfmodel.espr = prm.initial_ctfmodel.espr;
    	prm.current_ctfmodel.ispr = prm.initial_ctfmodel.ispr;
    	prm.current_ctfmodel.alpha = prm.initial_ctfmodel.alpha;
    	prm.current_ctfmodel.DeltaF = prm.initial_ctfmodel.DeltaF;
    	prm.current_ctfmodel.DeltaR = prm.initial_ctfmodel.DeltaR;
    	prm.current_ctfmodel.Q0 = prm.initial_ctfmodel.Q0;
        COPY_ctfmodel_TO_CURRENT_GUESS;
    }
    DEBUG_TEXTFILE(formatString("Step 6: espr=%f",prm.current_ctfmodel.espr));
    DEBUG_MODEL_TEXTFILE;
    /************************************************************************/
    /* STEP 7:  the defocus and angular parameters                          */
    /************************************************************************/

    prm.action = 3;
    prm.evaluation_reduction = 1;
    if (prm.fastDefocusEstimate)
        prm.estimate_defoci_Zernike();
    else
        prm.estimate_defoci();

    DEBUG_TEXTFILE(formatString("Step 7: DeltafU=%f",prm.current_ctfmodel.DeltafU));
    DEBUG_TEXTFILE(formatString("Step 7: DeltafV=%f",prm.current_ctfmodel.DeltafV));
    DEBUG_TEXTFILE(formatString("Step 7: azimutalAngle=%f",prm.current_ctfmodel.azimuthal_angle));
    DEBUG_MODEL_TEXTFILE;

    //This line is to test the results obtained
    //exit(1);

    /************************************************************************/
    /* STEPs 9, 10 and 11: all parameters included second Gaussian          */
    /************************************************************************/
    prm.action = 5;
    if (prm.modelSimplification < 2)
        prm.estimate_background_gauss_parameters2();

    steps.resize(ALL_CTF_PARAMETERS);
    steps.initConstant(1);
    steps(3) = 0; // kV
    steps(5) = 0; // The spherical aberration (Cs) is not optimized
    steps(37) = 0; //VPP radius not optimized
    if (prm.initial_ctfmodel.Q0 != 0)
        steps(15) = 0; // Q0
    if (prm.modelSimplification >= 3)
        steps(20) = steps(21) = steps(23) = 0;
    if (prm.modelSimplification >= 2)
        steps(24) = steps(25) = steps(26) = steps(27) = steps(28) = steps(29) = 0;
    if (prm.modelSimplification >= 1)
        steps(10) = steps(11) = 0;
    if (std::floor(prm.initial_ctfmodel.VPP_radius) == 0)
    	steps(36) = 0; //VPP phase shift

    powellOptimizer(*prm.adjust_params, 0 + 1, ALL_CTF_PARAMETERS, CTF_fitness,
                    global_prm, 0.01, fitness, iter, steps, prm.show_optimization);

    prm.current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (prm.show_optimization)
    {
        std::cout << "Best fit with Gaussian2:\n" << prm.current_ctfmodel << std::endl;
        prm.saveIntermediateResults("step04b_best_fit_with_gaussian2");
    }

    prm.evaluation_reduction = 2;
    powellOptimizer(*prm.adjust_params, 0 + 1, ALL_CTF_PARAMETERS, CTF_fitness,
                    global_prm, 0.01, fitness, iter, steps, prm.show_optimization);
    prm.current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    prm.evaluation_reduction = 1;
    powellOptimizer(*prm.adjust_params, 0 + 1, ALL_CTF_PARAMETERS, CTF_fitness,
                    global_prm, 0.005, fitness, iter, steps, prm.show_optimization);
    prm.current_ctfmodel.forcePhysicalMeaning();
    COPY_ctfmodel_TO_CURRENT_GUESS;

    if (prm.show_optimization)
    {
        std::cout << "Best fit:\n" << prm.current_ctfmodel << std::endl;
        prm.saveIntermediateResults("step04c_best_fit");
    }
    DEBUG_TEXTFILE(formatString("Step 11: DeltafU=%f fitness=%f",prm.current_ctfmodel.DeltafU,fitness));
    DEBUG_MODEL_TEXTFILE;

    //We adopt that always  DeltafU > DeltafV so if this is not the case we change the values and the angle
    if ( prm.current_ctfmodel.DeltafV > prm.current_ctfmodel.DeltafU)
    {
        double temp;
        temp = prm.current_ctfmodel.DeltafU;
        prm.current_ctfmodel.DeltafU = prm.current_ctfmodel.DeltafV;
        prm.current_ctfmodel.DeltafV = temp;
        prm.current_ctfmodel.azimuthal_angle -= 90;
        COPY_ctfmodel_TO_CURRENT_GUESS;
    }

    /************************************************************************/
    /* STEP 12: Produce output                                              */
    /************************************************************************/

    prm.action = 7;

    if (prm.fn_psd != "")
    {
        // Define mask between first and third zero
    	prm.mask_between_zeroes.initZeros(prm.mask);
        Matrix1D<double> u(2), z1(2), z3(2);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(prm.mask_between_zeroes)
        {
            VECTOR_R2(u, prm.x_digfreq(i, j), prm.y_digfreq(i, j));
            u /= u.module();
            prm.current_ctfmodel.lookFor(1, u, z1, 0);
            prm.current_ctfmodel.lookFor(3, u, z3, 0);
            if (z1.module() < prm.w_contfreq(i, j)
                && prm.w_contfreq(i, j) < z3.module())
            	prm.mask_between_zeroes(i, j) = 1;
        }
        // Evaluate the correlation in this region
        CTF_fitness(prm.adjust_params->adaptForNumericalRecipes(), &prm);
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

        prm.saveIntermediateResults(fn_rootMODEL, false);
        prm.current_ctfmodel.Tm /= prm.downsampleFactor;
        prm.current_ctfmodel.azimuthal_angle = std::fmod(prm.current_ctfmodel.azimuthal_angle,360.);
        prm.current_ctfmodel.write(fn_rootCTFPARAM + ".ctfparam_tmp");
        MetaData MD;
        MD.read(fn_rootCTFPARAM + ".ctfparam_tmp");
        size_t id = MD.firstObject();
        MD.setValue(MDL_CTF_X0, (double)output_ctfmodel.x0*prm.Tm, id);
        MD.setValue(MDL_CTF_XF, (double)output_ctfmodel.xF*prm.Tm, id);
        MD.setValue(MDL_CTF_Y0, (double)output_ctfmodel.y0*prm.Tm, id);
        MD.setValue(MDL_CTF_YF, (double)output_ctfmodel.yF*prm.Tm, id);
        MD.setValue(MDL_CTF_CRIT_FITTINGSCORE, fitness, id);
        MD.setValue(MDL_CTF_CRIT_FITTINGCORR13, prm.corr13, id);
        MD.setValue(MDL_CTF_CRIT_ICENESS, evaluateIceness(prm.psd_exp_radial, prm.Tm), id);
        MD.setValue(MDL_CTF_DOWNSAMPLE_PERFORMED, prm.downsampleFactor, id);
        MD.write(fn_rootCTFPARAM + ".ctfparam",MD_APPEND);
        fn_rootCTFPARAM = fn_rootCTFPARAM + ".ctfparam_tmp";
        fn_rootCTFPARAM.deleteFile();
    }
    output_ctfmodel = prm.current_ctfmodel;

    DEBUG_CLOSE_TEXTFILE;
    return fitness;
}

void ProgCTFEstimateFromPSD::run()
{
    CTFDescription ctfmodel;
    ROUT_Adjust_CTF(*this, ctfmodel);
}
