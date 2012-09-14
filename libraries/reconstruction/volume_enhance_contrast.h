/***************************************************************************
*
* Authors:    Carlos Oscar            coss@cnb.csic.es (2010)
*             Fernando Fuentes
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

#ifndef VOLUME_ENHANCE_CONTRAST_H
#define VOLUME_ENHANCE_CONTRAST_H

#include <string>
#include <data/morphology.h>
#include <data/filters.h>
#include <data/transformations.h>
#include <data/xmipp_program.h>


#include <queue>
#include <vector>
#include <iostream>


/// @defgroup EnhanceContrast Enhance constrast
/// @ingroup ReconsLibrary

/// Parameters for enhance contrast program
/// @ingroup Denoise
class ProgVolumeEnhanceContrast: public XmippProgram
{
public:
	// Input/Output volumes
	FileName fnIn, fnOut;

	// Confidence level for background identification
	double alpha;

	// Lower intensity (%)
	double lowerIntensity;

	// Remove the background
	bool removeBg;

	// Save mask
	FileName fnMask;

protected:
	/** Params definition */
	void defineParams();

    /** Read parameters from command line
     */
    void readParams();

    void run();

public:
    /** Show parameters. This function calls show_specific.
     */
    void show();

    /** Enhance contrast of a volume.
     */
    void enhance(MultidimArray< double >& vol);
};

#endif
