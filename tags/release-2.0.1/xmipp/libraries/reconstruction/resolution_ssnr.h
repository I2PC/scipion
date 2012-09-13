/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef _PROG_SSNR
#define _PROG_SSNR

#include <iostream>

#include <data/volume.h>
#include <data/selfile.h>

/**@defgroup SSNR resolution_ssnr (Spectral Signal to Noise Ratio)
   @ingroup ReconsLibraryPrograms */
//@{
/** SSNR parameters. */
class Prog_SSNR_prm
{
public:
    /// Signal reconstructed volume or a selfile with its projections
    FileName fn_S;

    /// Noise reconstructed volume or a selfile with its projections
    FileName fn_N;

    /// Selfile with all the experimental images
    FileName fn_Ssel;

    /// Selfile with all the noise images
    FileName fn_Nsel;

    /// Filename of the Volumetric SSNR, used only for radial averaging
    FileName fn_VSSNR;

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

    /** Generate VSSNR.*/
    bool generate_VSSNR;

    /** Generate radial average.*/
    bool radial_avg;

    /** Min_power: Threshold for not dividing */
    double min_power;
public:
    /* Side info -------------------------------------------------------- */
    // Signal volume
    VolumeXmipp S;

    // Noise volume
    VolumeXmipp N;

    // Selfile with all experimental images
    SelFile SF_S;

    // Selfile with all noisy images
    SelFile SF_N;

    // Selfile with all the images projected from the signal volume
    SelFile SF_Sth;

    // Selfile with all the images projected from the noise volume
    SelFile SF_Nth;

    // SSNR3D for the radial_avg
    VolumeXmipp VSSNR;

public:
    /// Read parameters from command line
    void read(int argc, char **argv);

    /// Show parameters
    friend ostream & operator << (ostream &out, const Prog_SSNR_prm &prm);

    /// Usage
    void usage() const;

    /// Produce side Info
    void produce_side_info();

    /** Estimate SSNR 2D.
        Generate images with the particular SSNR. The output filename
        is used as a rootname */
    void Estimate_SSNR(int dim, Matrix2D<double> &output);

    /** Radial average of a Volumetric SSNR.
        The Volumetric SSNR is stored as 10*log10(VSSNR+1). To perform
        a radial average that is consistent with the one produced
        by the 1D estimation the +1 must be properly eliminated.

        The columns of output are the following:
        Column 0: sample number in Fourier Space,
        Column 1: corresponding frequency in continuous freq (1/A),
        Column 2: corrected radial_avg
    */
    void Radial_average(Matrix2D<double> &output);
};

/** Perform all the work.
    For the meaning of the output matrix look at the documentation
    of the function Estimate_SSNR_1D of the class Prog_SSNR_prm. */
void ROUT_SSNR(Prog_SSNR_prm &prm, Matrix2D<double> &output);
//@}
#endif
