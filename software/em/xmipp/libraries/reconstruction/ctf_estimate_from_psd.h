/***************************************************************************
 *
 * Authors:      Javier Angel Velazquez Muriel    javi@cnb.csic.es
 *               Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef _ADJUST_CTF_HH
#define _ADJUST_CTF_HH

#include "fourier_filter.h"
#include "ctf_estimate_from_psd_base.h"
#include "ctf_estimate_from_psd_fast.h"

/**@defgroup AdjustParametricCTF adjust_ctf (Adjust CTF parameters to PSD)
   @ingroup ReconsLibrary */
//@{
/** Adjust CTF parameters. */
class ProgCTFEstimateFromPSD: public ProgCTFBasicParams
{
public:

    /// CTF model
    CTFDescription       initial_ctfmodel, current_ctfmodel, ctfmodel_defoci;

    ProgCTFEstimateFromPSD()
    {
    }
    //Copy useful initial values in 2D taken from 1D.
    ProgCTFEstimateFromPSD(ProgCTFEstimateFromPSDFast *copy)
    {
    	action = copy->action;
    	x_contfreq = copy->x_contfreq;
        y_contfreq = copy->y_contfreq;
        w_contfreq = copy->w_contfreq;
    	x_digfreq = copy->x_digfreq;
    	y_digfreq = copy->y_digfreq;
    	w_digfreq = copy->w_digfreq;
    	///PSD data
    	psd_exp_radial = copy->psd_exp_radial;
    	psd_exp_enhanced_radial = copy->psd_exp_enhanced_radial;
    	psd_exp_enhanced_radial_derivative = copy->psd_exp_enhanced_radial_derivative;
    	psd_theo_radial_derivative = copy->psd_theo_radial_derivative;
    	psd_exp_radial_derivative = copy->psd_exp_radial_derivative;
    	psd_theo_radial = copy->psd_theo_radial;
    	w_digfreq_r_iN = copy->w_digfreq_r_iN;
    	w_digfreq_r = copy->w_digfreq_r;
    	///Masks
    	mask = copy->mask;
    	mask_between_zeroes = copy->mask_between_zeroes;
    	max_freq = copy->max_freq;
    	min_freq = copy->min_freq;
    	min_freq_psd = copy->min_freq_psd;
    	max_freq_psd = copy->max_freq_psd;
    	corr13 = copy->corr13;

    	show_inf = copy->show_inf;
    	heavy_penalization = copy->heavy_penalization;
    	penalize = copy->penalize;
    	evaluation_reduction = copy->evaluation_reduction;
    	modelSimplification = copy->modelSimplification;
    	defocus_range = copy->defocus_range;
    	downsampleFactor = copy->downsampleFactor;

    	enhanced_ctftomodel() = copy->enhanced_ctftomodel();
    	enhanced_ctftomodel_fullsize() = copy->enhanced_ctftomodel_fullsize();
    	enhanced_weight = copy->enhanced_weight;

    	current_ctfmodel = copy->current_ctfmodel;
    	initial_ctfmodel = copy->initial_ctfmodel;
    	ctfmodel_defoci = copy->ctfmodel_defoci;
    	current_ctfmodel.precomputeValues(x_contfreq,y_contfreq);

    	Tm = copy->Tm;
    	f = copy->f;
    	ctfmodelSize = copy->ctfmodelSize;
    	show_optimization = copy->show_optimization;
    	fn_psd = copy->fn_psd;

    }

public:
    
    /// Read parameters
    void readBasicParams(XmippProgram *program);

    void readParams();

    /// Define basic parameters
    static void defineBasicParams(XmippProgram * program);

    /// Define Parameters
    void defineParams();

    /// Produce side information
    void produceSideInfo();

    /** Generate half-plane model at a given size.
        It is assumed that ROUT_Adjust_CTF has been already run */
    void generate_model_halfplane(int Ydim, int Xdim, MultidimArray<double> &model);

    /** Generate quadrant model at a given size.
        It is assumed that ROUT_Adjust_CTF has been already run */
    void generate_model_quadrant(int Ydim, int Xdim, MultidimArray<double> &model);

    /** Run */
    void run();

    /* Assign ctfmodel from a vector and viceversa ----------------------------- */
    void assignCTFfromParameters(double *p, CTFDescription &ctfmodel, int ia,
                                 int l, int modelSimplification);

    void assignParametersFromCTF(CTFDescription &ctfmodel, double *p, int ia,
                                 int l, int modelSimplification);

    /* Center focus ----------------------------------------------------------- */
    void center_optimization_focus(bool adjust_freq, bool adjust_th, double margin);

    /* The model is taken from global_adjust and global_ctfmodel is modified */
    void generateModelSoFar(Image<double> &I, bool apply_log);

    /* Compute central region -------------------------------------------------- */
    void compute_central_region(double &w1, double &w2, double ang);

    /* Save intermediate results ----------------------------------------------- */
    /* First call to generate model so far and then save the image, and a couple
     of cuts along X and Y.
     This function returns the fitting error.*/
    void saveIntermediateResults(const FileName &fn_root, bool generate_profiles);

    /** CTF fitness */
    double CTF_fitness_object(double *p);

    // Estimate sqrt parameters
    void estimate_background_sqrt_parameters();

    // Estimate first gaussian parameters ------------------------------------- */
    void estimate_background_gauss_parameters();

    // Estimate second gaussian parameters ------------------------------------- */
    void estimate_background_gauss_parameters2();

    // Estimate envelope parameters -------------------------------------------- */
    void estimate_envelope_parameters();

    // Estimate defoci --------------------------------------------------------- */
    void showFirstDefoci();
    void estimate_defoci();

    // Estimate defoci with Zernike and SPTH transform--------------------------------------------- */
    void estimate_defoci_Zernike(MultidimArray<double> &psdToModelFullSize, double min_freq, double max_freq, double Tm,
                                 double kV, double lambdaPhase, int sizeWindowPhase,
                                 double &defocusU, double &defocusV, double &ellipseAngle, int verbose);
    void estimate_defoci_Zernike();
};

double evaluateIceness(MultidimArray<double> &enhanced_ctftomodel, double Tm);

/** Core of the Adjust CTF routine.
    This is the routine which does everything. It returns the fitting error
    committed in the best fit.*/
double ROUT_Adjust_CTF(ProgCTFEstimateFromPSD &prm, CTFDescription &output_ctfmodel, 
    bool standalone = true);
//@}
#endif
