/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
 *
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

#include "score_micrograph.h"
#include "fourier_filter.h"
#include <data/args.h>
#include <data/filters.h>
#include <data/xmipp_fft.h>

/* Read parameters --------------------------------------------------------- */
void ProgScoreMicrograph::readParams()
{

	pieceDim = 256;
	overlap = 0.5;
	skipBorders = 2;
	Nsubpiece = 1;
	bootstrapN = -1;

	fn_micrograph = getParam("--micrograph");
	particleSize = getIntParam("--particleSize");

}

/* Usage ------------------------------------------------------------------- */
void ProgScoreMicrograph::defineParams()
{
	addUsageLine("Estimate the score of a given Micrograph attending on different parameters.");
    addParamsLine("   --micrograph <file>         : File with the micrograph");
    addParamsLine("  [--particleSize <p=100>]       : Size of the particle");
}

/* Apply ------------------------------------------------------------------- */
void ProgScoreMicrograph::run()
{

	prmEstimateCTFFromMicrograph.pieceDim = pieceDim;
	prmEstimateCTFFromMicrograph.overlap = overlap;
	prmEstimateCTFFromMicrograph.skipBorders = skipBorders;
	prmEstimateCTFFromMicrograph.Nsubpiece = 1;
	prmEstimateCTFFromMicrograph.bootstrapN = bootstrapN;
	prmEstimateCTFFromMicrograph.fn_micrograph = fn_micrograph;
	prmEstimateCTFFromMicrograph.fn_root = "/home/jvargas/Linux/Proyectos/LUMC/set_001_challenge/kk";

	//prmEstimateCTFFromMicrograph.PSDEstimator_mode = ProgCTFEstimateFromMicrograph.TPSDEstimator_mode.Periodogram;
	//	PSDEstimator_mode = prmEstimateCTFFromMicrograph.TPSDEstimator_mode->	Periodogram;

/*
	psd_mode = ProgCTFEstimateFromMicrograph::TPSD_mode::OnePerMicrograph;
	prmEstimateCTFFromMicrograph.psd_mode = psd_mode;

	PSDEstimator_mode = ProgCTFEstimateFromMicrograph::TPSDEstimator_mode::Periodogram;
	prmEstimateCTFFromMicrograph.PSDEstimator_mode = PSDEstimator_mode;
*/
	std::cout << prmEstimateCTFFromMicrograph.PSDEstimator_mode << std::endl;
	prmEstimateCTFFromMicrograph.run();

}

#undef DEBUG
