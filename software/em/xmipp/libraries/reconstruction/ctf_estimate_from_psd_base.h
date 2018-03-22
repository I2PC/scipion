/*
 * ctf_estimate_from_psd_base.h
 *
 *  Created on: Aug 29, 2017
 *      Author: javier mota <jmota@cnb.csic.es>
 */
#ifndef _ADJUST_CTF_BASIC_PARAMS
#define _ADJUST_CTF_BASIC_PARAMS

#include "fourier_filter.h"

class ProgCTFBasicParams: public XmippProgram
{
public:
	/// CTF filename
	FileName             fn_psd;
	/// Downsample performed
	double				 downsampleFactor;
	/// CTF amplitude to model
	Image<double>        ctftomodel;
	/// CTF amplitude to model
	Image<double>        enhanced_ctftomodel;
	/// CTF amplitude to model
	Image<double>        enhanced_ctftomodel_fullsize;
	/// Show convergence values
	bool                 show_optimization;
	/// X dimension of particle projections (-1=the same as the psd)
	int                  ctfmodelSize;
	/// Bootstrap estimation
	bool                 bootstrap;
	/// Fast defocus estimate
	bool                 fastDefocusEstimate;
	/// Regularization factor for the phase direction and unwrapping estimates (used in Zernike estimate)
	double               lambdaPhase;
	/// Size of the average window used during phase direction and unwrapping estimates (used in Zernike estimate)
	int                  sizeWindowPhase;
	/// Minimum frequency to adjust
	double               min_freq;
	/// Maximum frequency to adjust
	double               max_freq;
	/// Sampling rate
	double               Tm;
	/// Defocus range
	double               defocus_range;
	/// Enhancement filter low freq
	double               f1;
	/// Enhancement filter high freq
	double               f2;
	/// Weight of the enhanced image
	double               enhanced_weight;
	/// Set of parameters for the complete adjustment of the CTF
	Matrix1D<double>     adjust;
	/// Model simplification
	int                  modelSimplification;
	///Frequencies in axes
	MultidimArray<double> x_contfreq;
	MultidimArray<double> y_contfreq;
	MultidimArray<double> w_contfreq;
	MultidimArray<double> x_digfreq;
	MultidimArray<double> y_digfreq;
	MultidimArray<double> w_digfreq;
	///PSD data
	MultidimArray<double> psd_exp_radial;
	MultidimArray<double> psd_exp_enhanced_radial;
	MultidimArray<double> psd_exp_enhanced_radial_derivative;
	MultidimArray<double> psd_theo_radial_derivative;
	MultidimArray<double> psd_exp_radial_derivative;
	MultidimArray<double> psd_theo_radial;
	MultidimArray<double> w_digfreq_r_iN;
	MultidimArray<double> w_digfreq_r;
	///Masks
	MultidimArray<double> mask;
	MultidimArray<double> mask_between_zeroes;
	MultidimArray<double> w_count;
	//Global values
	// Autofocus
	double value_th;
	double min_freq_psd;
	double max_freq_psd;

	// Maximum of the gaussian
	double max_gauss_freq;

	int show_inf; // 0: Do not show
	// 1: Partially detailed
	// 2: Very detailed

	// Program status
	int action;
	// 0: Computing the background (sqrt)
	// 1: Computing the full background
	// 2: Computing the envelope
	// 3: Computing defoci
	// 4: Computing all CTF parameters
	// 5: Computing all CTF parameters + Gaussian2
	// 6: Computing all CTF parameters + Gaussian2
	// 7: Produce output

	double corr13;

	// Penalization factor for the background
	bool penalize;

	// Speed up factor
	int evaluation_reduction;

	// Penalization for forbidden values of the parameters
	double heavy_penalization;
	double current_penalty;
	static const double penalty; // Maximum penalization
	//Some aliases
	MultidimArray<double> *f; // The CTF to model
	Matrix1D<double> *adjust_params; // Current theoretical adjustment

	ProgCTFBasicParams()
	{
	}

public:
    /// Read parameters
    void readParams();

    /// Read parameters
    void readBasicParams(XmippProgram *program);

    /// Show parameters
    void show();

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

};

double evaluateIceness(MultidimArray<double> &psd, double Tm);

#endif
