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
#ifndef _PROG_PROJECT_CRYSTAL_HH
#define _PROG_PROJECT_CRYSTAL_HH

class ParametersProjection;
class PROJECT_Side_Info;
class Crystal_Projection_Parameters;

#include "project.h"

/**@defgroup CrystalProject project_crystal (Project crystals)
   @ingroup ReconsLibrary */
//@{
/* Crystal Projection parameters ------------------------------------------- */
/** Crystal Projecting parameters.
    This class reads a set of projection parameters in a file (see
    the user guide of the projection program
    for more information about the file structure) and extract the
    useful information from it. This class allows also to write in
    the same file format, for doing so you must first introduce the
    right values in the class fields, and the call to the procedure
    write.

    A crystal is supposed to have two lattice vectors: a and b as follows
    @code
     ^ b
     |
     |
     |
     -----------> a
    @endcode
*/
class Crystal_Projection_Parameters
{
public:
    /// Crystal X dimension
    int crystal_Xdim;
    /// Crystal Y dimension
    int crystal_Ydim;
    /// Orthogonalize projections
    bool orthogonal;
    /// Crystal vector a
    Matrix1D<double> a;
    /// Crustal vector b
    Matrix1D<double> b;
    /// Bias to apply to the magnitude shift
    double    Nshift_avg;
    /// Standard deviation of the magnitude shift
    double    Nshift_dev;
    /// Disappearing threshold
    double    disappearing_th;
    /// file with shifts
    FileName fn_shift;
    /// Document File for shifts. Order: H K x_SHIFT y_SHIFT z_SHIFT
    MetaData DF_shift;
    /// is doc file with shifts available
    bool     DF_shift_bool;

public:
    /** Empty constructor*/
    Crystal_Projection_Parameters();

    /** Read projection parameters from a file.
        An exception is thrown if the file is not found or any of the
        parameters is not found in the right place.*/
    void read(const FileName &fn_crystal, double scale = 1.0);

    /** Write projection parameters to a file.
        The projection parameters are written into a file wth the same
        structure as always. If the file cannot be openned for output
        an exception is thrown. */
    void write(const FileName &fn_crystal);
};

/* Project as crystal ------------------------------------------------------ */
/** Project a Mathematical volume as a crystal.
*/
void project_crystal(Phantom &phantom, Projection &P,
                     const ParametersProjection &prm,
                     PROJECT_Side_Info &side, const Crystal_Projection_Parameters &prm_crystal,
                     float rot, float tilt, float psi);

/** Find crystal limits.
    This function returns imin and imax such that imin*v and imax*v are
    still inside the projection. v is supposed to be in R2.
    proj_Corner1 and proj_corner2 are the top-left and right-bottom corners
    of the space to fit (for instance the projection).
    a and b are the lattice vectors where the unit cell defined
    by the two cell corners will be copied. iamin, iamax, ibmin, ibmax are
    the minimum and maximum indexes that you can use with these lattice
    vectors, unit cell and projection size such that the unit cell intersects
    the projection. */
void find_crystal_limits(
    const Matrix1D<double> &proj_corner1, const Matrix1D<double> &proj_corner2,
    const Matrix1D<double> &corner1, const Matrix1D<double> &corner2,
    const Matrix1D<double> &a, const Matrix1D<double> &b,
    int &iamin, int &iamax, int &ibmin, int &ibmax);

/** Move following a spiral.
    Starting in a given position this function visits a matrix using
    an spiral like this

                        DCBA
                        E329
                        F418
                        G567
                        HIJKL...

   r is the actual position in the matrix.
*/
void move_following_spiral(Matrix1D<double> &r, const MultidimArray<int> &visited);

/** Fill cell positions.
    This function returns the random shifts corresponding to all cells that
    can be seen in a projection using the lattice vectors given. It also
    returns a matrix telling if a lattice combination falls inside the
    projection or not. You must provide the vectors in the undeformed
    space (for orthogonal projections), the lattice vectors in the deformed
    space (if it is not deformed, supply the same as in the undeformed),
    and the corners of the projection of a single cell. Corners must be the
    unit cell corners in the deformed space. Output shifts are in the
    deformed space. */
void fill_cell_positions(Projection &P,
                         Matrix1D<double> &aproj,   Matrix1D<double> &bproj,
                         Matrix1D<double> &aprojd,  Matrix1D<double> &bprojd,
                         Matrix1D<double> &corner1, Matrix1D<double> &corner2,
                         const Crystal_Projection_Parameters &prm_crystal,
                         MultidimArray<double> &cell_shiftX, MultidimArray<double> &cell_shiftY,
                         MultidimArray<double> &cell_shiftZ,
                         MultidimArray<int> &cell_inside,
                         MultidimArray<double> &exp_shifts_matrix_X,
                         MultidimArray<double> &exp_shifts_matrix_Y,
                         MultidimArray<double> &exp_shifts_matrix_Z);
/** Fill auxiliary matrix with experimental shifts.
    the values 3D shifts stored in the doc file DF_shift are
    transfer to two 2D matrices

*/
void init_shift_matrix(const Crystal_Projection_Parameters &prm_crystal,
                       MultidimArray<int>    &cell_inside,
                       MultidimArray<double> &exp_shifts_matrix_X,
                       MultidimArray<double> &exp_shifts_matrix_Y,
                       MultidimArray<double> &exp_shifts_matrix_Z,
                       MultidimArray<double> &exp_normal_shifts_matrix_X,
                       MultidimArray<double> &exp_normal_shifts_matrix_Y,
                       MultidimArray<double> &exp_normal_shifts_matrix_Z,
                       double param_file_scale);

//@}
#endif
