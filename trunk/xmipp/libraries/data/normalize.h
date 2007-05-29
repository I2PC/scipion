/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef NORMALIZE_H
#define NORMALIZE_H

#include "matrix2d.h"
#include "progs.h"
#include "mask.h"

/// @defgroup Normalize Normalize program.

/** Normalize parameters.
 * @ingroup Normalize
 */
class Normalize_parameters: public Prog_parameters
{
public:

// FIXME Make this an enum or similar
#define NONE 0
#define OLDXMIPP 1
#define NEAR_OLDXMIPP 2
#define NEWXMIPP 3
#define MICHAEL 4
#define NEWXMIPP2 5
#define RANDOM 6
#define RAMP 7

    /** Normalizing method.
     * Valid methods are OLDXMIPP, NEAR_OLDXMIPP, NEWXMIPP, NEWXMIPP2, MICHAEL,
     * NONE, RANDOM.
     */
    int method;

    /** Nomalization of volumes.
     */
    bool volume;

// TODO Same thing, an enum
#define NONE 0
#define FRAME 1
#define CIRCLE 2

    /** Background mode.
     * Valid modes are NONE, FRAME, CIRCLE.
     */
    int background_mode;

    /** Frame width or circle radius.
     */
    int r;

    /** Upper limit of a in y=ax+b.
     */
    double aF;

    /** Lower limit of a in y=ax+b.
     */
    double a0;

    /** Upper limit of b in y=ax+b.
     */
    double bF;

    /** Lower limit of b in y=ax+b.
     */
    double b0;

    /** Flags for remving balck/white spots due to dust.
     */
    bool remove_black_dust;
    bool remove_white_dust;

    /** Threshold for removing black/white (dust) spots.
     */
    double thresh_black_dust;
    double thresh_white_dust;

    Matrix2D< int  > bg_mask;
    bool apply_geo;
    bool enable_mask;
    Mask_Params mask_prm;

    /** Read parameters from command line.
     */
    void read(int argc, char** argv);

    /** Produce side information.
     */
    void produce_side_info();

    /** Show parameters.
     */
    void show();

    /** Usage.
     */
    void usage();

    /** Apply inverse geometric transformation.
     * As stored in image header to the mask.
     */
    void apply_geo_mask(ImageXmipp& img);

    /** Apply to an image.
     * The input image is modified.
     */
    void apply(Image* img);
};

#endif
