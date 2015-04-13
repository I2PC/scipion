/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
#ifndef _PROG_ADJUST_VOLUME_HH
#define _PROG_ADJUST_VOLUME_HH

#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/multidim_array.h>
#include <data/xmipp_program.h>

/**@defgroup AdjustVolumeProgram adjust_volume_grey_values (Adjust volume grey values to a set of projections)
   @ingroup ReconsLibrary */
//@{
/* Adjust volume program*/
class ProgAdjustVolume: public XmippProgram
{
protected:
    /// Filename with the input volume
    FileName fn_vol;
    /// Filename with the input projections
    FileName fn_sel;
    /** Filename of the output volume.
        If empty the input one is used. */
    FileName fn_out;
    /// Optimize
    bool optimize;
    /// Probability of being evaluated
    double probb_eval;
    // Create temp file
    bool tempFile;

public:
    // Input volume
    Image<double> ImIn;
    MultidimArray<double> V;
    // SelFile
    MetaData SF;

protected:
    void defineParams();
    void readParams();

    /** Show parameters. */
    void show();

    /** Apply.
        This is the function that really does the job */
    void apply(MultidimArray<float> &output_volume);

public:
    /** Run. Calls apply and save the result. */
    void run();

    /** Mismatching. This function returns the overall mismatiching between the
     * experimental projections and the theoretical projections of the current
     * volume. */
    double mismatching(double a, double b);
};
//@}
#endif
