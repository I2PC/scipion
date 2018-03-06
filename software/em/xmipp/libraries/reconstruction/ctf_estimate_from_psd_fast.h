/***************************************************************************
 *
 * Authors:      Javier Mota Garcia    jmota@cnb.csic.es
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

#ifndef _ADJUST_CTF_FAST_HH
#define _ADJUST_CTF_FAST_HH

#include "fourier_filter.h"
#include "ctf_estimate_from_psd_base.h"

/**@defgroup AdjustParametricCTFFast adjust_ctf (Adjust CTF parameters to PSD)
   @ingroup ReconsLibrary */
//@{
/** Adjust CTF parameters. */
class ProgCTFEstimateFromPSDFast: public ProgCTFBasicParams
{
public:

	CTFDescription1D    initial_ctfmodel, current_ctfmodel, ctfmodel_defoci;

	ProgCTFEstimateFromPSDFast()
	{
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

	/** Run */
	void run();

	/* Assign ctfmodel from a vector and viceversa ----------------------------- */
	void assignCTFfromParameters_fast(double *p, CTFDescription1D &ctf1Dmodel, int ia, int l, int modelSimplification);

	void assignParametersFromCTF_fast(CTFDescription1D &ctfmodel, double *p, int ia, int l, int modelSimplification);

	/* Center focus ----------------------------------------------------------- */
	void center_optimization_focus_fast(bool adjust_freq, bool adjust_th, double margin);

	/* The model is taken from global_adjust and global_ctfmodel is modified */
	void generateModelSoFar_fast(MultidimArray<double> &I, bool apply_log);

	/* Compute central region -------------------------------------------------- */
	void compute_central_region_fast(double &w1, double &w2, double ang);

	/* Save intermediate results ----------------------------------------------- */
	/* First call to generate model so far and then save the image, and a couple
	 of cuts along X and Y.
	 This function returns the fitting error.*/
	void saveIntermediateResults_fast(const FileName &fn_root, bool generate_profiles);

	/** CTF fitness */
	double CTF_fitness_object_fast(double *p);

	// Estimate sqrt parameters
	void estimate_background_sqrt_parameters_fast();

	// Estimate first gaussian parameters ------------------------------------- */
	void estimate_background_gauss_parameters_fast();

	// Estimate second gaussian parameters ------------------------------------- */
	void estimate_background_gauss_parameters2_fast();

	// Estimate envelope parameters -------------------------------------------- */
	void estimate_envelope_parameters_fast();

	// Estimate defoci --------------------------------------------------------- */
	void showFirstDefoci_fast();
	void estimate_defoci_fast();

};

/** Core of the Adjust CTF routine.
    This is the routine which does everything. It returns the fitting error
    committed in the best fit.*/
double ROUT_Adjust_CTFFast(ProgCTFEstimateFromPSDFast &prm, CTFDescription1D &output_ctfmodel,
    bool standalone = true);
//@}
#endif
