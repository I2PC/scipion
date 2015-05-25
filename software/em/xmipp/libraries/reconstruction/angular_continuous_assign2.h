/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#ifndef _PROG_ANGULAR_PREDICT_CONTINOUOS2
#define _PROG_ANGULAR_PREDICT_CONTINUOUS2

#include <data/xmipp_program.h>
#include <data/ctf.h>
#include "fourier_projection.h"
#include "fourier_filter.h"

/**@defgroup AngularPredictContinuous2 angular_continuous_assign2 (Continuous angular assignment)
   @ingroup ReconsLibrary */
//@{
/** Predict Continuous Parameters. */
class ProgAngularContinuousAssign2: public XmippMetadataProgram
{
public:
    /** Filename of the reference volume */
    FileName fnVol;
    /** Filename of residuals */
    FileName fnResiduals;
    /** Maximum shift allowed */
    double maxShift;
    /** Maximum scale allowed */
    double maxScale;
    /** Maximum angular change allowed */
    double maxAngularChange;
    /** Maximum frequency (A) */
    double maxResol;
    /** Maximum defocus change (A) */
    double maxDefocusChange;
    /** Sampling rate */
    double Ts;
    /** Maximum radius */
    int Rmax;
    /** Padding factor */
    int pad;
    // Optimize gray
    bool optimizeGrayValues;
    // Optimize shift
    bool optimizeShift;
    // Optimize scale
    bool optimizeScale;
    // Optimize angles
    bool optimizeAngles;
    // Optimize defocus
    bool optimizeDefocus;
    // Apply transformation to this image
    String originalImageLabel;
    // Phase Flipped
    bool phaseFlipped;
    // Penalization for the average
    double penalization;
public:
    // 2D mask in real space
    MultidimArray<int> mask2D;
    // Inverse of the sum of Mask2D
    double iMask2Dsum;
    // Fourier projector
    FourierProjector *projector;
    // Reference volume
    Image<double> V;
    // Volume size
    size_t Xdim;
    // Input image
	Image<double> I, Ip, E, Ifiltered, Ifilteredp;
	// Theoretical projection
	Projection P;
	// Filter
    FourierFilter filter;
    // Transformation matrix
    Matrix2D<double> A;
    // Original angles
    double old_rot, old_tilt, old_psi;
    // Original shift
	double old_shiftX, old_shiftY;
	// Has CTF
	bool hasCTF;
	// Original defocus
	double old_defocusU, old_defocusV, old_defocusAngle;
	// CTF
	CTFDescription ctf;
	// Covariance matrices
	Matrix2D<double> C0, C;
public:
    /// Empty constructor
    ProgAngularContinuousAssign2();

    /// Destructor
    ~ProgAngularContinuousAssign2();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /** Start processing */
    void startProcessing();

    /** Produce side info.
        An exception is thrown if any of the files is not found*/
    void preProcess();

    /** Predict angles and shift.
        At the input the pose parameters must have an initial guess of the
        parameters. At the output they have the estimated pose.*/
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /** Post process */
    void postProcess();
};
//@}
#endif
