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
#ifndef _PROG_ENHANCE_PSD_HH
#  define _PROG_ENHANCE_PSD_HH

#include <data/xmipp_program.h>

/**@defgroup EnhancePSD psd_enhance (Enhance visualization of the PSD)
   @ingroup ReconsLibrary */
//@{
/* Enhance PSD Program Parameters ------------------------------------------ */
/** Parameter class for the project program */
class ProgCTFEnhancePSD: public XmippMetadataProgram
{
public:
	/// Method
	String method;

    /// Bandpass filter low frequency (in Fourier space, max 0.5)
    double filter_w1;

    /// Bandpass filter high frequency (in Fourier space, max 0.5)
    double filter_w2;

    /// Decay width (raised cosine)
    double decay_width;

    /// Minimum number of fringes
    int N0;

    /// Maximum number of fringes
    int NF;

    /// Lower frequency for the mask (in Fourier space, max 0.5)
    double mask_w1;

    /// Higher frequency for the mask (in Fourier space, max 0.5)
    double mask_w2;
public:
    /** Read from a command line.*/
    void readParams();

    /** Define parameters. */
    void defineParams();

    /** Show parameters. */
    void show();

    /** Process one image */
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /** Apply filter method to a single PSD.
        The steps are basically: outlier removal, band pass filtration, masking
        and normalization. */
    void applyFilter(MultidimArray<double> &PSD);

    /** Apply SPHT to a single PSD.*/
    void applySPHT(MultidimArray<double> &PSD);
};
//@}
#endif
