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

#ifndef ENHANCE_CONTRAST_H
#define ENHANCE_CONTRAST_H

#include <string>
#include <data/progs.h>
#include <data/morphology.h>
#include <data/filters.h>
#include <data/transformations.h>


#include <queue>
#include <vector>
#include <iostream>


/// @defgroup EnhanceContrast Enhance constrast
/// @ingroup ReconsLibraryPrograms

/// Parameters for enhance contrast program
/// @ingroup Denoise
class EnhanceContrast_parameters: public Prog_parameters
{
public:
	// Confidence level for background identification
	double alpha;

	// Lower intensity (%)
	double lowerIntensity;

	// Remove the background
	bool removeBg;

	// Save mask
	FileName fnMask;
public:
    /** Empty constructor
     */
    EnhanceContrast_parameters();

    /** Read parameters from command line
     */
    void read(int argc, char** argv);

    /** Produce side info.
     *
     * The DWT type is translated and set
     */
    void produce_side_info();

    /** Show parameters. This function calls show_specific.
     */
    void show();

    /** Usage. This function calls usage_specific.
     */
    void usage();

    /** Enhance contrast of a volume.
     */
    void enhance(MultidimArray< double >& vol);

};

#endif
