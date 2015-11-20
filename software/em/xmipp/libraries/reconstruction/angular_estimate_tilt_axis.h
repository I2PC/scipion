/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2002)
 *             Carlos Oscar Sorzano          (coss@cnb.csic.es)
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

#ifndef ANGULAR_ESTIMATE_TILT_AXIS_H
#define ANGULAR_ESTIMATE_TILT_AXIS_H

#include <data/xmipp_program.h>
#include <data/micrograph.h>

/**@defgroup AngularEstimateTiltAxis Estimate the tilt axis position in a pair of micrographs
   @ingroup ReconsLibrary */
//@{
/** Centilt parameters. */
class ProgAngularEstimateTiltAxis: public XmippProgram
{
public:
    /**  Filename input untilted coords */
    FileName fnUntilted;
    /**  Filename input tilted coords */
    FileName fnTilted;
    /**  Filename output document file */
    FileName fnOut;
public:
    /// Define parameters in the command line
    void defineParams();

    /// Read parameters from the command line
    void readParams();

    /// Run over the whole input metadata
    void run();
};
//@}
#endif
