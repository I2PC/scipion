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

#include <data/selfile.h>
#include <data/docfile.h>
#include <data/volume.h>

#define OVERSAMPLE 8

/**@defgroup CorrectAmplitude3D ctf_correct_amplitude3D (3D Wiener filtering)
   @ingroup ReconsLibraryPrograms */
//@{
/// Correct Amplitude3D parameters
class CorrectAmplitude3DParams
{
public:

    /// Filename for CTF datfile: a 2-column ASCII file with name of
    /// the volume and CTF parameter file for each defocus group
    FileName fnCtfdat;
    /// Filename for selfile with envelope ASCII files
    FileName fnEnvs;
    /// Rootname for output files
    FileName fnOut;
    /// Filename for docfile with number of images per defocus group
    FileName fnNrImgs;

    /// Wiener filter constant
    double wienConst;

    /// Low resolution cutoff to apply Wiener filter
    double minResol;

    /// Flag for phase flipped images
    bool isFlipped;

    /// Dimensions of the volumes
    int Zdim, Ydim, Xdim;

    /// Side Info: CTF
    FourierMask ctf;

    /// Side Info: ctfdat
    CTFDat ctfdat;

    /// Side Info: selfile with envelopes
    SelFile SFenv;

    /// Side Info: docfile with number of images in each defocus group
    DocFile DFimgs;

    /// The 3D CTFs and Wiener filters
    std::vector< Matrix1D<double> > Vctfs1D, Vwien1D;

public:
    /** Empty constructor */
    CorrectAmplitude3DParams(): wienConst(0)
    {}

    /** Read parameters from command line. */
    void read(int argc, char **argv);

    /** Show. */
    void show();

    /** Usage. */
    void usage();

    /** Produce side information.
        The CTFdat, nr_imgs docfile and selection file with envelopes are read. */
    void produceSideInfo();

    /** Generate 1D CTFs. */
    void generateCTF1D(const FileName &fnCTF,
                       const double nr_steps,
                       Matrix1D<double> &CTF1D);

    /** Generate Wiener filters */
    void generateWienerFilters();

    /** Generate deconvolved volume */
    void generateVolumes();

    /** Do the whole thing ... */
    void run();

};
//@}

