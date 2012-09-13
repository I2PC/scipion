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
#ifndef _PROG_IDR_ART_HH
#  define _PROG_IDR_ART_HH

#include <data/xmipp_image.h>
#include <data/xmipp_program.h>
#include <data/projection.h>
#include "fourier_filter.h"

/**@defgroup IDR ctf_correct_idr (Iterative Data Refinement for CTF amplitude correction)
   @ingroup ReconsLibrary */
//@{
/* IDR Parameters ---------------------------------------------------------- */
/** IDR Parameters. */
class ProgCtfCorrectIdr: public XmippMetadataProgram
{
public:
    /// Reference volume
    FileName fn_vol;

    /// Relaxation factor
    double mu;

    /// Output filename root
    FileName fnRoot;
public:
    // Input Volume
    Image<double> V;

    // Temporary projections
    Projection Ireal, Itheo;
    MultidimArray<double> Inorm,Itheo_CTF;

    // CTF filter
    FourierFilter ctf;

    // Last CTF file read
    FileName last_fn_ctf;

public:
    /// Read params
    void readParams();

    /// Preprocess
    void preProcess();

    /// Show parameters
    void show();

    /// Define parameters
    void defineParams();

    /// Process one image
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
};
//@}
#endif
