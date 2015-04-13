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

#ifndef _PROG_SSNR
#define _PROG_SSNR

#include <iostream>
#include <data/xmipp_program.h>
#include <data/xmipp_image.h>
#include <data/metadata.h>
#include <data/xmipp_fft.h>

/**@defgroup SSNR resolution_ssnr (Spectral Signal to Noise Ratio)
   @ingroup ReconsLibrary */
//@{
/** SSNR parameters. */

class ProgSSNR : public XmippProgram
{
public:
    /// Signal reconstructed volume
    FileName fn_S;
    /// Noise reconstructed volume
    FileName fn_N;
    /// Selfile with all the experimental and noise images
    FileName fn_SNsel;
    FileName fn_S_sel, fn_N_sel;
    /// Filename of the Volumetric SSNR, used only for radial averaging
    FileName fn_VSSNR;
    /// Fourier projections
    bool fourierProjections;
    /// Ringwidth
    double ring_width;
    /// Sampling rate
    double Tm;
    /** Output filename.
        If empty, SSNR is inserted before the extension in fn_S */
    FileName fn_out;
    /** Output rootname for the individual estimations.
        If empty, SSNR is inserted before the extension in fn_S */
    FileName fn_out_images;
    /** Symmetry */
    String sym;
    /** Generate VSSNR.*/
    bool generate_VSSNR;
    /** Generate radial average.*/
    bool radial_avg;
    /** Min_power: Threshold for not dividing */
    double min_power;
    /** Number of threads for ART */
    int Nthreads;
public:
    /* Side info -------------------------------------------------------- */
    // Signal volume
    Image<double> S;
    // Noise volume
    Image<double> N;
    // Selfile with all experimental images
    MetaData SF_SN, SF_S, SF_N;
    // SSNR3D for the radial_avg
    Image<double> VSSNR;

public:

    void defineParams();
    void readParams();
    void show();
    void produceSideInfo();
    void run();

    /** Estimate SSNR 2D.
        Generate images with the particular SSNR. The output filename
        is used as a rootname */
    void estimateSSNR(int dim, Matrix2D<double> &output);

    /** Radial average of a Volumetric SSNR.
        The Volumetric SSNR is stored as 10*log10(VSSNR+1). To perform
        a radial average that is consistent with the one produced
        by the 1D estimation the +1 must be properly eliminated.

        The columns of output are the following:
        Column 0: sample number in Fourier Space,
        Column 1: corresponding frequency in continuous freq (1/A),
        Column 2: corrected radial_avg
    */
    void radialAverage(Matrix2D<double> &output);
};
//@}
#endif
