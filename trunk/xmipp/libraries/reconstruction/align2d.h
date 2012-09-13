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
#ifndef _PROG_ALIGN2D
#define _PROG_ALIGN2D

#include <data/xmipp_fftw.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>
#include <data/metadata_extension.h>
#include <data/metadata.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>
#include <vector>

/**@defgroup Align2DProgram align2d (Align a set of 2D images)
   @ingroup ReconsLibrary */
//@{
/** Align2D parameters. */
class ProgAlign2d: public XmippProgram
{
public:
    /** Filename selection file containing the images */
    FileName fnSel;
    /** Filename reference image */
    FileName fnRef;
    /**  Filename output root */
    FileName fnRoot;
    /** Integer number of iterations to perform */
    int Niter;
    /** Do not check mirrors */
    bool dont_mirror;
    /** Do pspc */
    bool pspc;
public:
    // SelFile with the input images
    MetaData SF;
    // Image holding current reference
    Image<double> Iref;

public:
    /// Read argument
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /// Align pairs
    void alignPairs(MetaData &MDin, MetaData &MDout, int level);

    /// Pyramidal combination of images to construct a reference
    void do_pspc();

    /// Compute mean
    void computeMean();

    /// Alignment of all images by iterative refinement
    void refinement();

    /// Main routine
    void run();
};
//@}
#endif
