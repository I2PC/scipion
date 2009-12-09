/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Slavica Jonic (slavica.jonic@impmc.jussieu.fr)
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
#ifndef _PROG_SORT_PSD_HH
#  define _PROG_SORT_PSD_HH

#include <data/selfile.h>
#include <data/funcs.h>

/**@defgroup SortPSD psd_Sort (Sort visualization of the PSD)
   @ingroup ReconsLibraryPrograms */
//@{
/* Sort PSD Program Parameters ------------------------------------------ */
/** Parameter class for the project program */
class Prog_Sort_PSD_Parameters
{
public:
    /// List of micrographs to classify
    FileName fnSel;

    /// Name of the output file
    FileName fnOut;
    /// Size of the windows for the PSD estimation
    int windowSize;

    /// Bandpass filter low frequency (in Fourier space, max 0.5)
    double filter_w1;

    /// Bandpass filter high frequency (in Fourier space, max 0.5)
    double filter_w2;

    /// Decay width (raised cosine)
    double decay_width;

    /// Lower frequency for the mask (in Fourier space, max 0.5)
    double mask_w1;

    /// Higher frequency for the mask (in Fourier space, max 0.5)
    double mask_w2;

public:
    /// Selfile with images (side info)
    SelFile SF;

    /// List of correlation values (side info)
    Matrix1D<double> correlation;

public:
    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introdustd::cing this parameters. */
    void usage() const;

    /** Show parameters. */
    void show() const;

    /** Produce side information */
    void produceSideInfo();

    /** Compute the correlation value for a single micrograph */
    double computeCorrelation(const FileName &fnMicrograph) const;

    /** Compute the correlation for all micrographs */
    void run();
};
//@}
#endif
