/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.csic.es)
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

#include "fourier_filter.h"

#include <data/xmipp_image.h>
#include <data/metadata.h>
#include <data/xmipp_program.h>

/**@defgroup CorrectAmplitude3D ctf_correct_amplitude3D (3D Wiener filtering)
   @ingroup ReconsLibrary */
//@{
/// Correct Amplitude3D parameters
class ProgCtfCorrectAmplitude3D: public XmippProgram
{
public:
    /// Metadata with volume, ctf and number of images in that volume
    FileName fnIn;
    /// Rootname for output files
    FileName fnRoot;

    /// Wiener filter constant
    double wienerConstant;

    /// Low resolution cutoff to apply Wiener filter
    double minFreq;

    /// Flag for phase flipped images
    bool isFlipped;
public:
    /// Dimensions of the volumes
    size_t Zdim, Ydim, Xdim;
	
    /// Side Info: CTF
    FourierFilter ctf;
    
    /// Side Info: ctfdat
    MetaData ctfdat;

    /// The 3D CTFs and Wiener filters
    std::vector< MultidimArray<double> > Vctfs1D, Vwien1D;

public:
    /** Read parameters */
    void readParams();

    /** Show. */
    void show();

    /** Define Parameters*/
    void defineParams();

    /** Produce side information.
        The CTFdat, nr_imgs docfile and selection file with envelopes are read. */
    void produceSideInfo();

    /** Generate 1D CTFs. */
    void generateCTF1D(const FileName &fnCTF, 
		       size_t nr_steps,
		       MultidimArray<double> &CTF1D);

    /** Generate Wiener filters */
    void generateWienerFilters();

    /** Generate deconvolved volume */
    void generateVolumes();

    /** Do the whole thing ... */
    void run();
};
//@}

