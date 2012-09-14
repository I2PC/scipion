/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
#ifndef _PROG_THRESHOLD_HH
#define _PROG_THRESHOLD_HH

#include <data/xmipp_funcs.h>
#include <data/xmipp_image.h>
#include <data/mask.h>
#include <data/symmetries.h>
#include <data/xmipp_program.h>

/**@defgroup ThresholdProgram Threshold
   @ingroup ReconsLibrary */
//@{
/// Threshold Parameters
class ProgThreshold : public XmippMetadataProgram
{
public:
	// Which pixels should be selected
	String selectionMethod;

	// Threshold
	double threshold;

	// Substitution method
	String substitutionMethod;

	// New value
	double newValue;

	// Noise average
	double noiseAvg;

	// Noise stddev
	double noiseStddev;
public:
    /** Read parameters from command line. */
    void readParams();

    /** Define Parameters */
    void defineParams();

    /** Show parameters */
    void show();

    /// Process image or volume
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
public:
    // Selection method as integer
    int iSelectionMethod;
};
//@}
#endif
