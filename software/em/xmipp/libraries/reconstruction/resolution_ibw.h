/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Alvaro Capell
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
#ifndef _PROG_RESOLUTION_IBW_HH
#  define _PROG_RESOLUTION_IBW_HH

#include <data/xmipp_funcs.h>
#include <data/xmipp_program.h>

/**@defgroup ResolutionIBW Measure resolution through the inverse of border widths
   @ingroup ReconsLibrary */
//@{

/* Resolution IBW Program Parameters ------------------------------------------ */
/** Parameter class for the resolution IBW program.
 *  The program is based on P. Marziliano, F. Dufaux, S. Winkler, T. Ebrahimi.
 *  Perceptual blur and ringing metrics: application to JPEG2000. Signal Processing:
 *  Image Communication, 19: 163-172 (2004)
 */
class ProgResolutionIBW: public XmippProgram
{
public:
    /// Input volume
    FileName fnVol;

    /// Output volume with local widths
    FileName fnOut;

public:
    /// Volume to evaluate
    Image<double> V;

public:
    /** Read from a command line. */
    void readParams();

    /** Define parameters. */
    void defineParams();

    /** Show parameters. */
    void show() const;

    /** Compute the correlation for all micrographs */
    void run();

    /** Compute the edge width at the edge pixels in a given direction.
     * step controls how fine is the search of the border width.
     */
    void edgeWidth(const MultidimArray<double> &volCoeffs, const MultidimArray<double> &edges,
    				MultidimArray <double>& widths, const Matrix1D<double> &dir,
    				double step=0.25) const;

    /** Calculate the inverse border width.
     * The calculation is only performed at the border pixels.
     */
    double calculateIBW(MultidimArray <double>& widths) const;
};
//@}

#endif
