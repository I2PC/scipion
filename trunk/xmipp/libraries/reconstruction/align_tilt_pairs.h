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
#include <data/xmipp_fftw.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>

#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/xmipp_program.h>

/**@defgroup Centilt align_tilt_pairs (Align tilted and untilted images in a random conical tilt experiment)
   @ingroup ReconsLibrary */
//@{
/** Centilt parameters. */
class ProgAlignTiltPairs: public XmippProgram
{
protected:
    /** MetaData for untilted and tilted images */
	MetaData mdU, mdT;
	/** Object id for current image pair */
	size_t idU, idT;
    MDRow              rowU, rowT;
    /**  Filename output document file */
    FileName mdOut;
    /** Discard images that shift more than max_shift*/
    double max_shift;
    /** Force x-shift to be zero */
    bool force_x_zero;
    /** Perform cosine stretching */
    bool do_stretch;
    /** Perform centering */
    bool do_center;

    void defineParams();
    void readParams();
    void processImage();
    void run();

    /// Show
    void show();

    /// Center one tilted image
    bool centerTiltedImage(const Image<double> &Iu, Image<double> &It);
};
//@}
