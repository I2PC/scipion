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

#ifndef PROJECTION_H
#define PROJECTION_H

#include "image.h"

/// @defgroup Projections Projections.

/** Projection class.
 * @ingroup Projections
 *
 * A projection is an ImageXmipp plus some information (about the direction
 * of prejection) which makes it suitable for 3D reconstruction. A projection
 * is supposed to have the point (0,0) at the center of the image and not in
 * the corners as usual matrices have.
 *
 * The normal use of a projection is:
 *
 * @code
 * Projection P; // Create variable
 * P.reset(65, 65); // Init with zeros and set right origin
 * P.set_angles(30, 45, -50); // Set Euler angles
 * @endcode
 *
 * From now on, the projection can be treated as any other Image. You can even
 * write it on disk as it inherits from ImageXmipp.
 */
class Projection: public ImageXmipp
{
public:
    /** Vector perpendicular to the projection plane.
     * It is calculated as a function of rot and tilt.
     */
    matrix1D< double > direction;

    /** Matrix used to pass from the Universal coordinate system to the
     * projection coordinate system.
     *
     * @code
     * Rp = euler * Ru
     * @endcode
     */
    matrix2D< double > euler;

    /** Just the opposite.
     *
     * @code
     * Ru = eulert * Rp.
     * @endcode
     */
    matrix2D< double > eulert;

    /** Init_zeros and move origin to center.
     *
     * This function initialises the projection plane with 0's, and then moves
     * the logical origin of the image to the physical center of the image
     * (using the Xmipp conception of image center).
     */
    void reset(int Ydim, int Xdim);

    /** Set Euler angles for this projection.
     *
     * The Euler angles are stored in the Xmipp header, then the pass matrices
     * (Universe <---> Projection coordinate system) are computed, and the
     * vector perpendicular to this projection plane is also calculated.
     */
    void set_angles(double _rot, double _tilt, double _psi);

    /** Read a Projection from file.
     *
     * When a projection is read, the Euler matrices and perpendicular
     * direction is computed and stored in the Projection structure.
     */
    void read(const FileName& fn, const bool& apply_shifts = false);

    /** Assignment.
     */
    Projection& operator=(const Projection& P);

    /** Another function for assigment.
     */
    void assign(const Projection& P);
};

#endif
