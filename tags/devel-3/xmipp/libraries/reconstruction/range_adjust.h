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

#ifndef RANGEADJUST_H
#define RANGEADJUST_H

#include <data/progs.h>
#include <data/image.h>
#include <data/mask.h>
#include <data/args.h>

/// @defgroup RangeAdjust Adjust grey level range of images and volumes
/// @ingroup DataLibraryPrograms

/** Parameter class for the project program.
 * @ingroup RangeAdjust
 */
class Prog_Range_adjust_Parameters: public Prog_parameters
{
public:
    /// min_val.
    double min_val;

    /// max_val.
    double max_val;

    /// noise in %.
    double sigma;

    /// Mask
	Mask_Params mask_prm;
	
	/** Empty constructor */
	Prog_Range_adjust_Parameters(): Prog_parameters(), mask_prm(INT_MASK)
	{}

    /** Read from a command line.
     *
     * An exception might be thrown by any of the internal conversions, this
     * would mean that there is an error in the command line and you might
     * show a usage message.
     */
    void read(int argc, char** argv);

    /** Usage message.
     *
     * This function shows the way of introdustd::cing this parameters.
     */
    void usage();

    /** Show parameters.
     */
    void show();

    /** Range adjust of an image.
     *
     * The input image is modified.
     */
    void apply(Matrix2D< double >& I);

    /** Range adjust of a volume.
     *
     * The input volume is modified.
     */
    void apply(Matrix3D< double >& V);
};

#endif
