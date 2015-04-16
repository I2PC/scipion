/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

#ifndef ANGULAR_ALIGN_DUAL_TILT_SERIES_H_
#define ANGULAR_ALIGN_DUAL_TILT_SERIES_H_

#include <data/metadata.h>
#include <data/multidim_array.h>
#include <data/matrix2d.h>
#include <data/xmipp_program.h>

/**@defgroup AlignDualProgram Align two dual tilt series
   @ingroup ReconsLibrary */
//@{
/** Align dual tilt series parameters. */
class ProgAlignDualTiltSeries: public XmippProgram
{
public:
    /// Reference tilt series
    FileName fnRef;

    /// Dual tilt series
    FileName fnDual;

    /// Aligned tilt series
    FileName fnOut;

    /// Scale factor
    double scaleFactor;
public:
    /// Read parameters from command line
    void readParams();

    /// Show parameters
    void show();

    /// Usage
    void defineParams();

    /// Read dual series
    void readDual();

    /// Produce side info
    void produceSideInfo();

    /// Find parameters (shift+rotation) at 0 degrees
    void findParametersAt0degrees(bool rotateDual);

    /// Distance between a pair of common lines
    double distanceBetweenCommonLines(
        int refi, int dualj, const Matrix2D<double> &E,
        double X, double Y, double Z);

    /// Evaluate alignment
    double evaluateAlignment(const Matrix1D<double> &_alignment);

    /// Optimize current alignment
    double optimizeAlignment();

    /// Align dual
    void alignDual();

    /// Run
    void run();

    /// Shift the projection in Z
    void shiftProjectionInZ(MultidimArray<double> &I, int dualj, double Z) const;

public:
    // Show debugging information
    bool debugging;

    // Reference tilt series
    MetaData SFRef;

    // Dual tilt series
    MetaData SFDual;

    // List of angles in ref
    MultidimArray<double> tiltRef;

    // List of angles in dual
    MultidimArray<double> tiltDual;

    // Images in the reference tilt series
    std::vector< MultidimArray<double> > imgRef;

    // Images in the dual tilt series
    std::vector< MultidimArray<double> > imgDual;

    // Dual images are rotated 180 degrees
    bool rotatedDual;

    // Name of the images closer to 0 degrees
    FileName fnRef0, fnDual0;

    // Shift at 0 degrees
    Matrix1D<double> alignment;

    // Auxiliary vectors;
    Matrix1D<double> normali, normalj,
        commonline, commonlinejE, commonlinei, commonlinej;
    MultidimArray<double> profilei, profilej, I;
    Matrix2D<double> E, Et, A;

    // Image dimensions
    int Xdim, Ydim;
};
//@}
#endif /* ANGULAR_ALIGN_DUAL_TILT_SERIES_H_ */
