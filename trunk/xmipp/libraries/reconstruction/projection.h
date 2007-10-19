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
#ifndef _PROJECTION_HH
#  define _PROJECTION_HH

#include <data/image.h>
#include <data/projection.h>

#include "grids.h"
#include "basis.h"

template <class T>
void project_SimpleGrid(VolumeT<T> &vol, const SimpleGrid &grid,
                        const Basis &basis,
                        Projection &proj, Projection &norm_proj, int FORW, int eq_mode,
                        const VolumeT<int> *VNeq, Matrix2D<double> *M,
                        double ray_length = -1.0);


/*---------------------------------------------------------------------------*/
/* PROJECTION                                                                */
/*---------------------------------------------------------------------------*/
/**@name Projections */
//@{
// Projecting functions ====================================================
#define FORWARD  1
#define BACKWARD 0
/**@name Single particle projections */
//@{
#define ARTK     1
#define CAVK     2
#define COUNT_EQ 3
#define CAV      4
#define CAVARTK  5

/** From basis volumes (Forward & Backward(ART)).
    Project a grid volume with a basis.
    The Grid volume is projected onto a projection plane defined by
    (rot, tilt, psi) (1st, 2nd and 3rd Euler angles). The projection
    is previously is resized to Ydim x Xdim and initialized to 0.
    The projection itself, from now on, will keep the Euler angles.

    FORWARD process:
       Each volume of the grid is projected on to the projection plane.
       The output is the projection itself and a normalising image, the
       normalising image is the projection of the same grid supposing
       that all basis are of value 1. This normalising image is used by
       the ART process

    BACKWARD process:
       During the backward process the normalising projection contains
       the correction image to apply to the volume (in the ART sense).
       The output is the volume itself, the projection image is useless
       in this case, and the normalising projection is not modified at
       all.

    As for the mode, valid modes are ARTK, CAVK, COUNT_EQ, CAVARTK.

    M is the matrix corresponding to the projection process.
    */
template <class T>
void project_Volume(GridVolumeT<T> &vol, const Basis &basis,
                    Projection &proj, Projection &norm_proj, int Ydim, int Xdim,
                    double rot, double tilt, double psi, int FORW, int eq_mode = ARTK,
                    GridVolumeT<int> *GVNeq = NULL, Matrix2D<double> *M = NULL,
                    double ray_length = -1);

/** From voxel volumes.

    The voxel volume is projected onto a projection plane defined by
    (rot, tilt, psi) (1st, 2nd and 3rd Euler angles) . The projection
    is previously is resized to Ydim x Xdim and initialized to 0.
    The projection itself, from now on, will keep the Euler angles.
    
    The offset is a 3D vector specifying the offset that must be applied
    when going from the projection space to the universal space
    
    rproj=E*r+roffset => r=E^t (rproj-roffset)
    
    Set it to NULL if you don't want to use it
 */
void project_Volume(Matrix3D<double> &V, Projection &P, int Ydim, int Xdim,
                    double rot, double tilt, double psi,
		    const Matrix1D<double> *roffset=NULL);

/** From voxel volumes, off-centered tilt axis.
    This routine projects a volume that is rotating <angle> degrees
    around the axis defined by the two angles <axisRot,axisTilt> and
    that passes through the point raxis. The projection can be futher
    inplane rotated and shifted through the parameters
    <inplaneRot> and <rinplane>.
    
    All vectors involved must be 3D.
    
    The projection model is rproj=H Rinplane Raxis r + 
                                  Rinplane (I-Raxis) raxis + rinplane
				  
    Where Raxis is the 3D rotation matrix given by the axis and
    the angle.
*/
void project_Volume_offCentered(Matrix3D<double> &V, Projection &P,
   int Ydim, int Xdim, double axisRot, double axisTilt,
   const Matrix1D<double> &raxis, double angle, double inplaneRot,
   const Matrix1D<double> &rinplane);

/** Projects a single particle into a voxels volume by updating its components this way:

 Voxel(i,j,k) = Voxel(i,j,k) + Pixel( x,y) * Distance.

 Where:

 - Voxel(i,j,k) is the voxel the ray is crossing.
 - Pixel( y,z ) is the value of the pixel where the ray departs from.
 - Distance is the distance the ray goes through the Voxel.
*/
void singleWBP(Matrix3D<double> &V, Projection &P);

/** Count equations in volume.
   For Component AVeraing (CAV), the number of equations in which
   each basis is involved is needed. */
void count_eqs_in_projection(GridVolumeT<int> &GVNeq,
                             const Basis &basis, Projection &read_proj);
//@}

/**@name Crystal projections */
//@{
/** Project a crystal basis volume.
    This function projects a crystal deformed basis volume, ie, in the
    documentation volume g. However the angles given must be those for
    volume f, the undeformed one. You must supply the deformed lattice
    vectors, and the matrix to pass from the deformed to the undeformed
    vectors (D and Dinv). a=D*ai;

    Valid eq_modes are ARTK, CAVARTK and CAV.
*/
void project_Crystal_Volume(GridVolume &vol, const Basis &basis,
                            Projection &proj, Projection &norm_proj,
                            int Ydim, int Xdim,
                            double rot, double tilt, double psi, const Matrix1D<double> &shift,
                            const Matrix1D<double> &aint, const Matrix1D<double> &bint,
                            const Matrix2D<double> &D, const Matrix2D<double> &Dinv,
                            const Matrix2D<int> &mask, int FORW, int eq_mode = ARTK);
//@}

//@}

// Implementations =========================================================
// Some aliases
#define x0   STARTINGX(IMGMATRIX(proj))
#define xF   FINISHINGX(IMGMATRIX(proj))
#define y0   STARTINGY(IMGMATRIX(proj))
#define yF   FINISHINGY(IMGMATRIX(proj))
#define xDim XSIZE(IMGMATRIX(proj))
#define yDim YSIZE(IMGMATRIX(proj))

// Projections from single particles #######################################
// Projections from basis volumes ==========================================

/* Projection of a simple volume ------------------------------------------- */
// The projection is not precleaned (set to 0) before projecting and its
// angles are supposed to be already written (and all Euler matrices
// precalculated
// The projection plane is supposed to pass through the Universal coordinate
// origin

/* Time measures ...........................................................
   This function has been optimized in time, several approaches have been
   tried and here are the improvements in time
   Case: Base (unit time measure)
         The points within the footprint are calculated inside the
         inner loop instead of computing the foot coordinate for the
         first corner and knowing the sampling rate in the oversampled
         image go on the next points. The image access was done directly
         (physical access).

   Notation simplification:
   Case: VOLVOXEL(V,k,i,j) ----> V(k,i,j):   + 38% (Inacceptable)
   Case: DIRECT_IMGPIXEL   ----> IMGPIXEL:   +  5% (Acceptable)

   Algorithmic changes:
   Case: project each basis position ---->   +325% (Inacceptable)
         get rid of all beginZ, beginY, prjX, ... and project each
         basis position onto the projection plane
   Case: footprint coordinates outside ----> - 33% (Great!!)
         the footprint coordinates are computed outside the inner
         loop, but you have to use the sampling rate instead to
         move over the oversampled image.
*/

//#define DEBUG_LITTLE
//#define DEBUG
const int ART_PIXEL_SUBSAMPLING = 2;

template <class T>
void project_SimpleGrid(VolumeT<T> &vol, const SimpleGrid &grid,
                        const Basis &basis,
                        Projection &proj, Projection &norm_proj, int FORW, int eq_mode,
                        const VolumeT<int> *VNeq, Matrix2D<double> *M,
                        double ray_length)
{
    Matrix1D<double> zero(3);                // Origin (0,0,0)
    static Matrix1D<double> prjPix(3);       // Position of the pixel within the
    // projection
    static Matrix1D<double> prjX(3);         // Coordinate: Projection of the
    static Matrix1D<double> prjY(3);         // 3 grid vectors
    static Matrix1D<double> prjZ(3);
    static Matrix1D<double> prjOrigin(3);    // Coordinate: Where in the 2D
    // projection plane the origin of
    // the grid projects
    static Matrix1D<double> prjDir(3);       // Projection direction

    static Matrix1D<double> actprj(3);       // Coord: Actual position inside
    // the projection plane
    static Matrix1D<double> beginZ(3);       // Coord: Plane coordinates of the
    // projection of the 3D point
    // (z0,YY(lowest),XX(lowest))
    static Matrix1D<double> univ_beginY(3);  // Coord: coordinates of the
    // grid point
    // (z0,y0,XX(lowest))
    static Matrix1D<double> univ_beginZ(3);  // Coord: coordinates of the
    // grid point
    // (z0,YY(lowest),XX(lowest))
    static Matrix1D<double> beginY(3);       // Coord: Plane coordinates of the
    // projection of the 3D point
    // (z0,y0,XX(lowest))
    double XX_footprint_size;                // The footprint is supposed
    double YY_footprint_size;                // to be defined between
    // (-vmax,+vmax) in the Y axis,
    // and (-umax,+umax) in the X axis
    // This footprint size is the
    // R2 vector (umax,vmax).

    int XX_corner2, XX_corner1;              // Coord: Corners of the
    int YY_corner2, YY_corner1;              // footprint when it is projected
    // onto the projection plane
    int           foot_V1, foot_U1;          // Img Coord: coordinate (in
    // an image fashion, not in an
    // oversampled image fashion)
    // inside the blobprint of the
    // corner1
    int           foot_V, foot_U;            // Img Coord: coordinate
    // corresponding to the blobprint
    // point which matches with this
    // pixel position
    int           Vsampling, Usampling;      // Sampling rate in Y and X
    // directions respectively
    // inside the blobprint
    double        vol_corr;                  // Correction for a volum element
    int           N_eq;                      // Number of equations in which
    // a blob is involved
    int           i, j, k;                   // volume element indexes

    // Prepare system matrix for printing ...................................
    if (M != NULL)
        M->initZeros(YSIZE(proj())*XSIZE(proj()), grid.get_number_of_samples());

    // Project grid axis ....................................................
    // These vectors ((1,0,0),(0,1,0),...) are referred to the grid
    // coord. system and the result is a 2D vector in the image plane
    // The unit size within the image plane is 1, ie, the same as for
    // the universal grid.
    // actprj is reused with a different purpose
    VECTOR_R3(actprj, 1, 0, 0);
    grid.Gdir_project_to_plane(actprj, proj.euler, prjX);
    VECTOR_R3(actprj, 0, 1, 0);
    grid.Gdir_project_to_plane(actprj, proj.euler, prjY);
    VECTOR_R3(actprj, 0, 0, 1);
    grid.Gdir_project_to_plane(actprj, proj.euler, prjZ);

    // This time the origin of the grid is in the universal coord system
    Uproject_to_plane(grid.origin, proj.euler, prjOrigin);

    // Get the projection direction .........................................
    proj.euler.getRow(2, prjDir);

    // Footprint size .......................................................
    // The following vectors are integer valued vectors, but they are
    // stored as real ones to make easier operations with other vectors
    if (basis.type == Basis::blobs)
    {
        XX_footprint_size = basis.blobprint.umax();
        YY_footprint_size = basis.blobprint.vmax();
        Usampling         = basis.blobprint.Ustep();
        Vsampling         = basis.blobprint.Vstep();
    }
    else if (basis.type == Basis::voxels || basis.type == Basis::splines)
    {
        YY_footprint_size = XX_footprint_size = CEIL(basis.max_length());
        Usampling = Vsampling = 0;
    }
    XX_footprint_size += XMIPP_EQUAL_ACCURACY;
    YY_footprint_size += XMIPP_EQUAL_ACCURACY;

    // Project the whole grid ...............................................
    // Corner of the plane defined by Z. These coordinates try to be within
    // the valid indexes of the projection (defined between STARTING and
    // FINISHING values, but it could be that they may lie outside.
    // These coordinates need not to be integer, in general, they are
    // real vectors.
    // The vectors returned by the projection routines are R3 but we
    // are only interested in their first 2 components, ie, in the
    // in-plane components
    beginZ = XX(grid.lowest) * prjX + YY(grid.lowest) * prjY + ZZ(grid.lowest) * prjZ
             + prjOrigin;

    // This type conversion gives more speed
    int ZZ_lowest = (int) ZZ(grid.lowest);
    int YY_lowest = (int) YY(grid.lowest);
    int XX_lowest = (int) XX(grid.lowest);
    int ZZ_highest = (int) ZZ(grid.highest);
    int YY_highest = (int) YY(grid.highest);
    int XX_highest = (int) XX(grid.highest);

    // Check if in VSSNR
    bool VSSNR_mode = (ray_length == basis.max_length());

#ifdef DEBUG_LITTLE
    int condition;
    condition = true;
    if (condition)
    {
        cout << "Equation mode " << eq_mode << endl;
        cout << "Footprint size " << YY_footprint_size << "x"
        << XX_footprint_size << endl;
        cout << "rot=" << proj.Phi() << " tilt=" << proj.Theta()
        << " psi=" << proj.Psi() << endl;
        cout << "Euler matrix " << proj.euler;
        cout << "Projection direction " << prjDir << endl;
        cout << grid;
        cout << "image limits (" << x0 << "," << y0 << ") (" << xF << ","
        << yF << ")\n";
        cout << "prjX           " << prjX.transpose()      << endl;
        cout << "prjY           " << prjY.transpose()      << endl;
        cout << "prjZ           " << prjZ.transpose()      << endl;
        cout << "prjOrigin      " << prjOrigin.transpose() << endl;
        cout << "beginZ(coord)  " << beginZ.transpose()    << endl;
        cout << "lowest         " << XX_lowest  << " " << YY_lowest
        << " " << XX_lowest  << endl;
        cout << "highest        " << XX_highest << " " << YY_highest
        << " " << XX_highest << endl;
        cout << "Stats of input basis volume ";
        vol().print_stats();
        cout << endl;
        cout.flush();
    }
#endif

    // Compute the grid lattice vectors in space ............................
    static Matrix2D<double> grid_basis(3, 3);
    grid_basis = grid.basis * grid.relative_size;
    static Matrix1D<double> gridX(3);  // Direction of the grid lattice vectors
    static Matrix1D<double> gridY(3);  // in universal coordinates
    static Matrix1D<double> gridZ(3);
    grid_basis.getCol(0, gridX);
    grid_basis.getCol(1, gridY);
    grid_basis.getCol(2, gridZ);

    univ_beginZ = XX(grid.lowest) * gridX + YY(grid.lowest) * gridY +
                  ZZ(grid.lowest) * gridZ + grid.origin;

    static Matrix1D<double> univ_position(3);
    int number_of_basis = 0;
    for (k = ZZ_lowest; k <= ZZ_highest; k++)
    {
        // Corner of the row defined by Y
        beginY = beginZ;
        univ_beginY = univ_beginZ;
        for (i = YY_lowest; i <= YY_highest; i++)
        {
            // First point in the row
            actprj = beginY;
            univ_position = univ_beginY;
            for (j = XX_lowest; j <= XX_highest; j++)
            {
                // Ray length interesting
                bool ray_length_interesting = true;
                if (ray_length != -1)
                    ray_length_interesting =
                        ABS(point_plane_distance_3D(univ_position, zero,
                                                    proj.direction)) <= ray_length;

                if (grid.is_interesting(univ_position) &&
                    ray_length_interesting)
                {
                    // Be careful that you cannot skip any basis, although its
                    // value be 0, because it is useful for norm_proj
#ifdef DEBUG
                    condition = (k == 0);
                    if (condition)
                    {
                        printf("\nProjecting grid coord (%d,%d,%d) ", j, i, k);
                        cout << "Vol there = " << VOLVOXEL(vol, k, i, j) << endl;
                        printf(" 3D universal position (%f,%f,%f) \n",
                               XX(univ_position), YY(univ_position), ZZ(univ_position));
                        cout << " Center of the basis proj (2D) " << XX(actprj) << "," << YY(actprj) << endl;
                        Matrix1D<double> aux;
                        Uproject_to_plane(univ_position, proj.euler, aux);
                        cout << " Center of the basis proj (more accurate) " << aux.transpose() << endl;
                    }
#endif

                    // Search for integer corners for this basis
                    XX_corner1 = CEIL(XMIPP_MAX(x0, XX(actprj) - XX_footprint_size));
                    YY_corner1 = CEIL(XMIPP_MAX(y0, YY(actprj) - YY_footprint_size));
                    XX_corner2 = FLOOR(XMIPP_MIN(xF, XX(actprj) + XX_footprint_size));
                    YY_corner2 = FLOOR(XMIPP_MIN(yF, YY(actprj) + YY_footprint_size));

#ifdef DEBUG
                    if (condition)
                    {
                        cout << "Clipped and rounded Corner 1 " << XX_corner1
                        << " " << YY_corner1 << " " << endl;
                        cout << "Clipped and rounded Corner 2 " << XX_corner2
                        << " " << YY_corner2 << " " << endl;
                    }
#endif

                    // Check if the basis falls outside the projection plane
                    // COSS: I have changed here
                    if (XX_corner1 <= XX_corner2 && YY_corner1 <= YY_corner2)
                    {
                        // Compute the index of the basis for corner1
                        if (basis.type == Basis::blobs)
                            OVER2IMG(basis.blobprint, (double)YY_corner1 - YY(actprj),
                                     (double)XX_corner1 - XX(actprj), foot_V1, foot_U1);

                        if (!FORW) vol_corr = 0;

                        // Effectively project this basis
                        // N_eq=(YY_corner2-YY_corner1+1)*(XX_corner2-XX_corner1+1);
                        N_eq = 0;
                        foot_V = foot_V1;
                        for (int y = YY_corner1; y <= YY_corner2; y++)
                        {
                            foot_U = foot_U1;
                            for (int x = XX_corner1; x <= XX_corner2; x++)
                            {
#ifdef DEBUG
                                if (condition)
                                {
                                    cout << "Position in projection (" << x << ","
                                    << y << ") ";
                                    double y, x;
                                    if (basis.type == Basis::blobs)
                                    {
                                        cout << "in footprint ("
                                        << foot_U << "," << foot_V << ")";
                                        IMG2OVER(basis.blobprint, foot_V, foot_U, y, x);
                                        cout << " (d= " << sqrt(y*y + x*x) << ") ";
                                        fflush(stdout);
                                    }
                                }
#endif
                                double a, a2;
                                // Check if volumetric interpolation (i.e., SSNR)
                                if (VSSNR_mode)
                                {
                                    // This is the VSSNR case
                                    // Get the pixel position in the universal coordinate
                                    // system
                                    SPEED_UP_temps;
                                    VECTOR_R3(prjPix, x, y, 0);
                                    M3x3_BY_V3x1(prjPix, proj.eulert, prjPix);
                                    V3_MINUS_V3(prjPix, prjPix, univ_position);
                                    a = basis.value_at(prjPix);
                                    a2 = a * a;
                                }
                                else
                                {
                                    // This is normal reconstruction from projections
                                    if (basis.type == Basis::blobs)
                                    {
                                        // Projection of a blob
                                        a = IMGPIXEL(basis.blobprint, foot_V, foot_U);
                                        a2 = IMGPIXEL(basis.blobprint2, foot_V, foot_U);
                                    }
                                    else
                                    {
                                        // Projection of other bases
                                        // If the basis is big enough, then
                                        // it is not necessary to integrate at several
                                        // places. Big enough is being greater than
                                        // 1.41 which is the maximum separation
                                        // between two pixels
                                        if (XX_footprint_size > 1.41)
                                        {
                                            // Get the pixel in universal coordinates
                                            SPEED_UP_temps;
                                            VECTOR_R3(prjPix, x, y, 0);
                                            // Express the point in a local coordinate system
                                            M3x3_BY_V3x1(prjPix, proj.eulert, prjPix);
#ifdef DEBUG
                                            if (condition)
                                                cout << " in volume coord ("
                                                << prjPix.transpose() << ")";
#endif
                                            V3_MINUS_V3(prjPix, prjPix, univ_position);
#ifdef DEBUG
                                            if (condition)
                                                cout << " in voxel coord ("
                                                << prjPix.transpose() << ")";
#endif
                                            a = basis.projection_at(prjDir, prjPix);
                                            a2 = a * a;
                                        }
                                        else
                                        {
                                            // If the basis is too small (of the
                                            // range of the voxel), then it is
                                            // necessary to sample in a few places
                                            const double p0 = 1.0 / (2 * ART_PIXEL_SUBSAMPLING) - 0.5;
                                            const double pStep = 1.0 / ART_PIXEL_SUBSAMPLING;
                                            const double pAvg = 1.0 / (ART_PIXEL_SUBSAMPLING * ART_PIXEL_SUBSAMPLING);
                                            int ii, jj;
                                            double px, py;
                                            a = 0;
#ifdef DEBUG
                                            if (condition) cout << endl;
#endif
                                            for (ii = 0, px = p0; ii < ART_PIXEL_SUBSAMPLING; ii++, px += pStep)
                                                for (jj = 0, py = p0; jj < ART_PIXEL_SUBSAMPLING; jj++, py += pStep)
                                                {
#ifdef DEBUG
                                                    if (condition)
                                                        cout << "    subsampling (" << ii << ","
                                                        << jj << ") ";
#endif
                                                    SPEED_UP_temps;
                                                    // Get the pixel in universal coordinates
                                                    VECTOR_R3(prjPix, x + px, y + py, 0);
                                                    // Express the point in a local coordinate system
                                                    M3x3_BY_V3x1(prjPix, proj.eulert, prjPix);
#ifdef DEBUG
                                                    if (condition)
                                                        cout << " in volume coord ("
                                                        << prjPix.transpose() << ")";
#endif
                                                    V3_MINUS_V3(prjPix, prjPix, univ_position);
#ifdef DEBUG
                                                    if (condition)
                                                        cout << " in voxel coord ("
                                                        << prjPix.transpose() << ")";
#endif
                                                    a += basis.projection_at(prjDir, prjPix);
#ifdef DEBUG
                                                    if (condition)
                                                        cout << " partial a="
                                                        << basis.projection_at(prjDir, prjPix)
                                                        << endl;
#endif
                                                }
                                            a *= pAvg;
                                            a2 = a * a;
#ifdef DEBUG
                                            if (condition)
                                                cout << "   Finally ";
#endif
                                        }
                                    }
                                }
#ifdef DEBUG
                                if (condition) cout << "=" << a << " , " << a2;
#endif
                                if (FORW)
                                {
                                    switch (eq_mode)
                                    {
                                    case CAVARTK:
                                    case ARTK:
                                        IMGPIXEL(proj, y, x) += VOLVOXEL(vol, k, i, j) * a;
                                        IMGPIXEL(norm_proj, y, x) += a2;
                                        if (M != NULL)
                                        {
                                            int py, px;
                                            proj().toPhysical(y, x, py, px);
                                            int number_of_pixel = py * XSIZE(proj()) + px;
                                            (*M)(number_of_pixel, number_of_basis) = a;
                                        }
                                        break;
                                    case CAVK:
                                        IMGPIXEL(proj, y, x) += VOLVOXEL(vol, k, i, j) * a;
                                        IMGPIXEL(norm_proj, y, x) += a2 * N_eq;
                                        break;
                                    case COUNT_EQ:
                                        VOLVOXEL(vol, k, i, j)++;
                                        break;
                                    case CAV:
                                        IMGPIXEL(proj, y, x) += VOLVOXEL(vol, k, i, j) * a;
                                        IMGPIXEL(norm_proj, y, x) += a2 *
                                                                     VOLVOXEL(*VNeq, k, i, j);
                                        break;
                                    }

#ifdef DEBUG
                                    if (condition)
                                    {
                                        cout << " proj= " << IMGPIXEL(proj, y, x)
                                        << " norm_proj=" << IMGPIXEL(norm_proj, y, x) << endl;
                                        cout.flush();
                                    }
#endif
                                }
                                else
                                {
                                    vol_corr += IMGPIXEL(norm_proj, y, x) * a;
                                    if (a != 0) N_eq++;
#ifdef DEBUG
                                    if (condition)
                                    {
                                        cout << " corr_img= " << IMGPIXEL(norm_proj, y, x)
                                        << " correction=" << vol_corr << endl;
                                        cout.flush();
                                    }
#endif
                                }

                                // Prepare for next operation
                                foot_U += Usampling;
                            }
                            foot_V += Vsampling;
                        } // Project this basis

                        if (!FORW)
                        {
                            T correction = 0;
                            switch (eq_mode)
                            {
                            case ARTK:
                                correction = (T) vol_corr;
                                break;
                            case CAVARTK:
                                if (N_eq != 0) correction = (T)(vol_corr / N_eq);
                                break;
                            case CAVK:
                            case CAV:
                                correction = (T) vol_corr;
                                break;
                            }
                            VOLVOXEL(vol, k, i, j) += correction;

#ifdef DEBUG
                            if (condition)
                            {
                                printf("\nFinal value at (%d,%d,%d) ", j, i, k);
                                cout << " = " << VOLVOXEL(vol, k, i, j) << endl;
                                cout.flush();
                            }
#endif
                        }
                    } // If not collapsed
                    number_of_basis++;
                } // If interesting

                // Prepare for next iteration
                V2_PLUS_V2(actprj, actprj, prjX);
                V3_PLUS_V3(univ_position, univ_position, gridX);
            }
            V2_PLUS_V2(beginY, beginY, prjY);
            V3_PLUS_V3(univ_beginY, univ_beginY, gridY);
        }
        V2_PLUS_V2(beginZ, beginZ, prjZ);
        V3_PLUS_V3(univ_beginZ, univ_beginZ, gridZ);
    }
}

#undef DEBUG
#undef DEBUG_LITTLE

/* Project a Grid Volume --------------------------------------------------- */
//#define DEBUG
//#define DEBUG_LITTLE
template <class T>
void project_Volume(
    GridVolumeT<T> &vol,                  // Volume
    const Basis &basis,                   // Basis
    Projection       &proj,               // Projection
    Projection       &norm_proj,          // Projection of a unitary volume
    int              Ydim,                // Dimensions of the projection
    int              Xdim,
    double rot, double tilt, double psi,  // Euler angles
    int              FORW,                // 1 if we are projecting a volume
    //   norm_proj is calculated
    // 0 if we are backprojecting
    //   norm_proj must be valid
    int              eq_mode,             // ARTK, CAVARTK, CAVK or CAV
    GridVolumeT<int> *GVNeq,              // Number of equations per blob
    Matrix2D<double> *M,                  // System matrix
    double            ray_length)         // Ray length of the projection
{
    // If projecting forward initialise projections
    if (FORW)
    {
        proj.reset(Ydim, Xdim);
        proj.set_angles(rot, tilt, psi);
        norm_proj().copyShape(proj());
        norm_proj().initZeros();
    }

#ifdef DEBUG_LITTLE
    if (FORW)
    {
        cout << "Number of volumes: " << vol.VolumesNo() << endl
        << "YdimxXdim: " << Ydim << "x" << Xdim << endl;
        for (int i = 0; i < vol.VolumesNo(); i++)
            cout << "Volume " << i << endl << vol.grid(i) << endl;
    }

#endif

    // Project each subvolume
    for (int i = 0; i < vol.VolumesNo(); i++)
    {
        VolumeT<int> *VNeq;
        if (GVNeq != NULL)   VNeq = &((*GVNeq)(i));
        else VNeq = NULL;
        project_SimpleGrid(vol(i), vol.grid(i), basis,
                           proj, norm_proj, FORW, eq_mode, VNeq, M, ray_length);
#ifdef DEBUG
        ImageXmipp save;
        save = norm_proj;
        if (FORW) save.write((string)"PPPnorm_FORW" + (char)(48 + i));
        else      save.write((string)"PPPnorm_BACK" + (char)(48 + i));
#endif
    }
}
#undef DEBUG_LITTLE
#undef DEBUG
#undef x0
#undef xF
#undef y0
#undef yF
#undef xDim
#undef yDim

#endif
