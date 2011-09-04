/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

#ifndef _PROG_EVALUATE_CLASSES
#define _PROG_EVALUATE_CLASSES

#include <data/xmipp_fftw.h>
#include <data/metadata.h>
#include <data/xmipp_program.h>

/**@defgroup EvaluateClasses Evaluate Classes
   @ingroup ReconsLibrary */
//@{

/// Class containing evaluators
class ClassEvaluation {
public:
	/** Digital frequency  <0.5) at which the Fourier Ring Correlation drops below 0.5 */
	double FRC_05;

	/** Differential phase residual at the frequency at which FRC=0.5 */
	double DPR_05;

	/** Possible overfitting */
	bool overfitting;

	/** Empty constructor */
	ClassEvaluation();
};

/** Evaluate the class contained in the Metadata.
 * It is assumed that the Metadata contains alignment parameters.
 */
void evaluateClass(MetaData &MD, ClassEvaluation &eval);

class ProgEvaluateClass: public XmippProgram
{
public:
	/// Name of the classification metadata
    FileName fnClass;
    /// Name of the output metadata
    FileName fnOut;
public:
    /// Read input parameters
    void readParams();

    /// Define input parameters
    void defineParams();

    /// Show
    void show();

    /// Execute
    void run();
};

//@}
#endif

