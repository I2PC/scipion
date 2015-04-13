/***************************************************************************
*
* Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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

#ifndef MAKESPECTRA_H
#define MAKESPECTRA_H

#include <data/rotational_spectrum.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>

#include <vector>

/// @defgroup MakeSpectra Make spectra
/// @ingroup ReconsLibrary

/** Make spectra parameters.
 * @ingroup MakeSpectra
 */
class ProgMakeSpectra: public XmippProgram
{
public:
	/** Input file */
	FileName fn_in;

	/** Output file */
	FileName fn_out;

    /** Rotational spectrum.
     */
    Rotational_Spectrum rot_spt;

public:
    /** Define parameters */
    void defineParams();

    /** Read parameters from command line. */
    void readParams();

    /** Show parameters. */
    void show();

    /** Run.
     */
    void run();
};

#endif
