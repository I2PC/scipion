/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2002)
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

#include <data/fftw.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/docfile.h>
#include <data/selfile.h>
#include <data/image.h>
#include <data/geometry.h>
#include <data/filters.h>

/**@defgroup Centilt align_tilt_pairs (Align tilted and untilted images in a random conical tilt experiment)
   @ingroup ReconsLibraryPrograms */
//@{
/** Centilt parameters. */
class Prog_centilt_prm
{
public:
    /** SelFiles for untilted and tilted images */
    SelFile SFu, SFt;
    /**  Filename output document file */
    FileName fn_doc;
    /**  Filename output extension */
    FileName oext;
    /** Discard images that shift more than max_shift*/
    double max_shift;
    /** Force x-shift to be zero */
    bool force_x_zero;
    /** Perform cosine stretching */
    bool do_stretch;
    /** Perform centering */
    bool do_center;

public:
    /// Read arguments from command line
    void read(int argc, char **argv);

    /// Show
    void show();

    /// Usage
    void usage();

    /// Center one tilted image
    bool center_tilted_image(const ImageXmipp &Iu, ImageXmipp &It, double &ccf);

    /// Main routine
    void centilt();
};
//@}
