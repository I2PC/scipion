/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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
#ifndef _PROG_PROJECTION_NEIGHBOURHOOD
#define _PROG_PROJECTION_NEIGHBOURHOOD

#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/xmipp_program.h>

#include <data/symmetries.h>

/**@defgroup ProjectionNeighbourhood angular_neighbourhood (Projection neighbourhood)
   @ingroup ReconsLibrary */
//@{
/** Projection Neighbourhood parameters. */
class ProgAngularNeighbourhood: public XmippProgram
{
public:
    /** Filename reference projection directions document file */
    FileName fn_ref;
    /** Filename selection file containing the images */
    FileName fn_sel;
    /** Double maximum neighbour distance */
    double maxdist;
    /** Filename symmetry file */
    FileName fn_sym;
    /** Filename of output file*/
    FileName fn_out;
    /** Check mirrors */
    bool check_mirrors;
public:
    // DocFile images
    MetaData DF1;
    // DocFile nbhoods
    MetaData DF2;
    // Symmetry List
    SymList SL;
public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /** Check symmetries.
        Given two sets of angles, this function modifies set 2 so that
        the difference with set 1 are minimized searching in the symmetry
        list and the second set. Return the angle distance.*/
    double check_symmetries(
        double rot1, double tilt1, double &rot2, double &tilt2);

    /** Compute distance.
        Compute the distance between the two document files loaded. The
        average distance is returned.*/
    void run();
};
//@}
#endif
