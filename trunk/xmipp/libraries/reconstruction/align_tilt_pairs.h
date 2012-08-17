/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2002)
 *             Carlos Oscar Sorzano          (coss@cnb.csic.es)
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
public:
    /**  Filename input document file */
    FileName fnIn;
    /**  Filename untilted average */
    FileName fnRef;
    /**  Filename output document file */
    FileName fnOut;
    /** Discard images that shift more than max_shift*/
    double max_shift;
    /** Do stretch */
    bool do_stretch;
    /** Do not align tilted */
    bool do_not_align_tilted;
public:
    /// Define parameters in the command line
    void defineParams();

    /// Read parameters from the command line
    void readParams();

    /// Show
    void show();

    /// Process a single image
    void processImage();

    /// Run over the whole input metadata
    void run();

    /// Center one tilted image
    bool centerTiltedImage(const MultidimArray<double> &imgU, bool flip,
    		double inPlaneU, double shiftXu, double shiftYu,
    		double alphaT, double alphaU, double tilt,
    		MultidimArray<double> &imgT,
    		double &shiftX, double &shiftY, CorrelationAux &auxCorr);
};
//@}
