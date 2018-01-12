/***************************************************************************
 * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
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

#ifndef _PROG_ELIMINATE_EMPTY_PARTICLES
#define _PROG_ELIMINATE_EMPTY_PARTICLES

#include <data/xmipp_program.h>

class ProgEliminateEmptyParticles: public XmippProgram
{
public:
	/// Name of the input metadata
    FileName fnIn;

    /// Name of the output metadata
    FileName fnOut;

    /// Name of the eliminated particle metadata
    FileName fnElim;

    /// Threshold for variance of variance
    float threshold;

    /// Add features
    bool addFeatures;

    /// Turning off denoising
    bool noDenoising;
public:
    /// Read input parameters
    void readParams();

    /// Show
    void show();

    /// Define input parameters
    void defineParams();

    /// Execute
    void run();
};

#endif
