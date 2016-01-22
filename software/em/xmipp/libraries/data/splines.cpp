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

#include "splines.h"
#include "integration.h"

/* Bspline03 by a LUT ------------------------------------------------------ */
double Bspline03LUT(double x)
{
    static bool firstCall = true;
    static MultidimArray<double> table(BSPLINE03_SUBSAMPLING);
    static const double deltax = 2.0 / BSPLINE03_SUBSAMPLING;
    static const double ideltax = 1.0 / deltax;
    if (firstCall)
    {
        FOR_ALL_ELEMENTS_IN_ARRAY1D(table)
        table(i) = Bspline03(i * deltax);
        firstCall = false;
    }
    size_t i = (size_t)round(fabs(x) * ideltax);
    if (i >= XSIZE(table)) return 0;
    else return table(i);
}

/* Sum spline on a grid ---------------------------------------------------- */
double sum_spatial_Bspline03_SimpleGrid(const SimpleGrid &grid)
{
    Matrix1D<double> gr(3), ur(3), corner1(3), corner2(3);
    int          i, j, k;
    double        sum = 0.0;

// Compute the limits of the spline in the grid coordinate system
    grid.universe2grid(vectorR3(-2.0, -2.0, -2.0), corner1);
    grid.universe2grid(vectorR3(2.0, 2.0, 2.0), corner2);

// Compute the sum in the points inside the grid
// The integer part of the vectors is taken for not picking points
// just in the border of the spline, which we know they are 0.
    for (i = (int)XX(corner1); i <= (int)XX(corner2); i++)
        for (j = (int)YY(corner1); j <= (int)YY(corner2); j++)
            for (k = (int)ZZ(corner1); k <= (int)ZZ(corner2); k++)
            {
                VECTOR_R3(gr, i, j, k);
                grid.grid2universe(gr, ur);
                sum += spatial_Bspline03(ur);
            }
    return sum;
}

double sum_spatial_Bspline03_Grid(const Grid &grid)
{
    double sum = 0;
    for (size_t i = 0; i < grid.GridsNo(); i++)
        sum += sum_spatial_Bspline03_SimpleGrid(grid(i));
    return sum;
}

/* Line integral through a spline ------------------------------------------ */
static Matrix1D<double> global_r;
static Matrix1D<double> global_u;
static Matrix1D<double> global_aux(3);
double spatial_Bspline03_integralf(double alpha)
{
    XX(global_aux) = XX(global_r) + alpha * XX(global_u);
    YY(global_aux) = YY(global_r) + alpha * YY(global_u);
    ZZ(global_aux) = ZZ(global_r) + alpha * ZZ(global_u);
    return spatial_Bspline03LUT(global_aux);
}

double spatial_Bspline03_integral(const Matrix1D<double> &r,
                                  const Matrix1D<double> &u, double alpha0, double alphaF)
{
    global_r = r;
    global_u = u;
    return integrateNewtonCotes(&spatial_Bspline03_integralf, alpha0, alphaF, 9);
}

/* This is very similar to the projection of a voxel volume from -2 to 2.
   The code is taken from the function project_Volume of Src/projection.cc */
//#define DEBUG
double spatial_Bspline03_proj(
    const Matrix1D<double> &r, const Matrix1D<double> &u)
{
    // Avoids divisions by zero and allows orthogonal rays computation
    static Matrix1D<double> ur(3);
    if (XX(u) == 0) XX(ur) = XMIPP_EQUAL_ACCURACY;
    else XX(ur) = XX(u);
    if (YY(u) == 0) YY(ur) = XMIPP_EQUAL_ACCURACY;
    else YY(ur) = YY(u);
    if (ZZ(u) == 0) ZZ(ur) = XMIPP_EQUAL_ACCURACY;
    else ZZ(ur) = ZZ(u);

    // Some precalculated variables
    double x_sign = SGN(XX(ur));
    double y_sign = SGN(YY(ur));
    double z_sign = SGN(ZZ(ur));

    // Compute the minimum and maximum alpha for the ray
    double alpha_xmin = (-2 - XX(r)) / XX(ur);
    double alpha_xmax = (2 - XX(r)) / XX(ur);
    double alpha_ymin = (-2 - YY(r)) / YY(ur);
    double alpha_ymax = (2 - YY(r)) / YY(ur);
    double alpha_zmin = (-2 - ZZ(r)) / ZZ(ur);
    double alpha_zmax = (2 - ZZ(r)) / ZZ(ur);

    double alpha_min = XMIPP_MAX(XMIPP_MIN(alpha_xmin, alpha_xmax),
                           XMIPP_MIN(alpha_ymin, alpha_ymax));
    alpha_min = XMIPP_MAX(alpha_min, XMIPP_MIN(alpha_zmin, alpha_zmax));
    double alpha_max = XMIPP_MIN(XMIPP_MAX(alpha_xmin, alpha_xmax),
                           XMIPP_MAX(alpha_ymin, alpha_ymax));
    alpha_max = XMIPP_MIN(alpha_max, XMIPP_MAX(alpha_zmin, alpha_zmax));
    if (alpha_max - alpha_min < XMIPP_EQUAL_ACCURACY) return 0.0;

#ifdef DEBUG
    std::cout << "Pixel:  " << r.transpose() << std::endl
    << "Dir:    " << ur.transpose() << std::endl
    << "Alpha x:" << alpha_xmin << " " << alpha_xmax << std::endl
    << "        " << (r + alpha_xmin*ur).transpose() << std::endl
    << "        " << (r + alpha_xmax*ur).transpose() << std::endl
    << "Alpha y:" << alpha_ymin << " " << alpha_ymax << std::endl
    << "        " << (r + alpha_ymin*ur).transpose() << std::endl
    << "        " << (r + alpha_ymax*ur).transpose() << std::endl
    << "Alpha z:" << alpha_zmin << " " << alpha_zmax << std::endl
    << "        " << (r + alpha_zmin*ur).transpose() << std::endl
    << "        " << (r + alpha_zmax*ur).transpose() << std::endl
    << "alpha  :" << alpha_min  << " " << alpha_max  << std::endl
    << std::endl;
#endif

    // Compute the first point in the volume intersecting the ray
    static Matrix1D<double> v(3);
    V3_BY_CT(v, ur, alpha_min);
    V3_PLUS_V3(v, r, v);

#ifdef DEBUG
    std::cout << "First entry point: " << v.transpose() << std::endl;
    std::cout << "   Alpha_min: " << alpha_min << std::endl;
#endif

    // Follow the ray
    double alpha = alpha_min;
    double ray_sum = 0;
    do
    {
        double alpha_x = (XX(v) + x_sign - XX(r)) / XX(ur);
        double alpha_y = (YY(v) + y_sign - YY(r)) / YY(ur);
        double alpha_z = (ZZ(v) + z_sign - ZZ(r)) / ZZ(ur);

        // Which dimension will ray move next step into?, it isn't necessary to be only
        // one.
        double diffx = ABS(alpha - alpha_x);
        double diffy = ABS(alpha - alpha_y);
        double diffz = ABS(alpha - alpha_z);
        double diff_alpha = XMIPP_MIN(XMIPP_MIN(diffx, diffy), diffz);
        ray_sum += spatial_Bspline03_integral(r, ur, alpha, alpha + diff_alpha);

        // Update alpha and the next entry point
        if (ABS(diff_alpha - diffx) <= XMIPP_EQUAL_ACCURACY) alpha = alpha_x;
        if (ABS(diff_alpha - diffy) <= XMIPP_EQUAL_ACCURACY) alpha = alpha_y;
        if (ABS(diff_alpha - diffz) <= XMIPP_EQUAL_ACCURACY) alpha = alpha_z;
        XX(v) += diff_alpha * XX(ur);
        YY(v) += diff_alpha * YY(ur);
        ZZ(v) += diff_alpha * ZZ(ur);

#ifdef DEBUG
        std::cout << "Alpha x,y,z: " << alpha_x << " " << alpha_y
        << " " << alpha_z << " ---> " << alpha << std::endl;

        std::cout << "    Next entry point: " << v.transpose() << std::endl
        << "    diff_alpha: " << diff_alpha << std::endl
        << "    ray_sum: " << ray_sum << std::endl
        << "    Alfa tot: " << alpha << "alpha_max: " << alpha_max <<
        std::endl;

#endif
    }
    while ((alpha_max - alpha) > XMIPP_EQUAL_ACCURACY);
    return ray_sum;
}
#undef DEBUG

/* Splines -> Voxels for a SimpleGrid -------------------------------------- */
void spatial_Bspline032voxels_SimpleGrid(const MultidimArray<double> &vol_splines,
        const SimpleGrid &grid,
        MultidimArray<double> *vol_voxels,
        const MultidimArray<double> *vol_mask = NULL)
{
    Matrix1D<double> act_coord(3);           // Coord: Actual position inside
    // the voxel volume without deforming
    Matrix1D<double> beginZ(3);              // Coord: Voxel coordinates of the
    // blob at the 3D point
    // (z0,YY(lowest),XX(lowest))
    Matrix1D<double> beginY(3);              // Coord: Voxel coordinates of the
    // blob at the 3D point
    // (z0,y0,XX(lowest))
    Matrix1D<double> corner2(3), corner1(3); // Coord: Corners of the
    // blob in the voxel volume
    Matrix1D<double> gcurrent(3);            // Position in g of current point
    double        intx, inty, intz;          // Nearest integer voxel
    int           i, j, k;                   // Index within the blob volume
    bool          process;                   // True if this blob has to be
    // processed
    double spline_radius = 2 * sqrt(3.0);

    // Some aliases
#define x0 STARTINGX(*vol_voxels)
#define xF FINISHINGX(*vol_voxels)
#define y0 STARTINGY(*vol_voxels)
#define yF FINISHINGY(*vol_voxels)
#define z0 STARTINGZ(*vol_voxels)
#define zF FINISHINGZ(*vol_voxels)

#ifdef DEBUG
    bool condition = true;
    (*vol_voxels)().printShape();
    std::cout << std::endl;
    std::cout << "x0= " << x0 << " xF= " << xF << std::endl;
    std::cout << "y0= " << y0 << " yF= " << yF << std::endl;
    std::cout << "z0= " << z0 << " zF= " << zF << std::endl;
    std::cout << grid;
#endif

    // Convert the whole grid ...............................................
    // Corner of the plane defined by Z. These coordinates are in the
    // universal coord. system
    grid.grid2universe(grid.lowest, beginZ);

    Matrix1D<double> grid_index(3);
    for (k = (int) ZZ(grid.lowest); k <= (int) ZZ(grid.highest); k++)
    {
        // Corner of the row defined by Y
        beginY = beginZ;
        for (i = (int) YY(grid.lowest); i <= (int) YY(grid.highest); i++)
        {
            // First point in the row
            act_coord = beginY;
            for (j = (int) XX(grid.lowest); j <= (int) XX(grid.highest); j++)
            {
                VECTOR_R3(grid_index, j, i, k);
#ifdef DEBUG
                if (condition)
                {
                    printf("Dealing spline at (%d,%d,%d) = %f\n", j, i, k, A3D_ELEM(vol_splines, k, i, j));
                    std::cout << "Center of the blob      "
                    << act_coord.transpose() << std::endl;
                }
#endif

                // These two corners are also real valued
                process = true;
                if (XX(act_coord) >= xF) process = false;
                if (YY(act_coord) >= yF) process = false;
                if (ZZ(act_coord) >= zF) process = false;
                if (XX(act_coord) <= x0) process = false;
                if (YY(act_coord) <= y0) process = false;
                if (ZZ(act_coord) <= z0) process = false;
#ifdef DEBUG
                if (!process && condition) std::cout << "   It is outside output volume\n";
#endif
                if (!grid.is_interesting(act_coord))
                {
#ifdef DEBUG
                    if (process && condition) std::cout << "   It is not interesting\n";
#endif
                    process = false;
                }

                if (process)
                {
                    V3_PLUS_CT(corner1, act_coord, -spline_radius);
                    V3_PLUS_CT(corner2, act_coord, spline_radius);
#ifdef DEBUG
                    if (condition)
                    {
                        std::cout << "Corner 1 for this point " << corner1.transpose() << std::endl;
                        std::cout << "Corner 2 for this point " << corner2.transpose() << std::endl;
                    }
#endif

                    // Clip the corners to the volume borders
                    XX(corner1) = ROUND(CLIP(XX(corner1), x0, xF));
                    YY(corner1) = ROUND(CLIP(YY(corner1), y0, yF));
                    ZZ(corner1) = ROUND(CLIP(ZZ(corner1), z0, zF));
                    XX(corner2) = ROUND(CLIP(XX(corner2), x0, xF));
                    YY(corner2) = ROUND(CLIP(YY(corner2), y0, yF));
                    ZZ(corner2) = ROUND(CLIP(ZZ(corner2), z0, zF));
#ifdef DEBUG
                    if (condition)
                    {
                        std::cout << "Clipped and rounded Corner 1 " << corner1.transpose() << std::endl;
                        std::cout << "Clipped and rounded Corner 2 " << corner2.transpose() << std::endl;
                    }
#endif

                    // Effectively convert
                    for (intz = ZZ(corner1); intz <= ZZ(corner2); intz++)
                        for (inty = YY(corner1); inty <= YY(corner2); inty++)
                            for (intx = XX(corner1); intx <= XX(corner2); intx++)
                            {
                                int iz = (int)intz, iy = (int)inty, ix = (int)intx;
                                if (vol_mask != NULL)
                                    if (!A3D_ELEM(*vol_mask, iz, iy, ix)) continue;

                                // Compute the spline value at this point
                                VECTOR_R3(gcurrent, intx, inty, intz);
                                V3_MINUS_V3(gcurrent, act_coord, gcurrent);
                                double spline_value = spatial_Bspline03(gcurrent);
#ifdef DEBUG_MORE
                                if (condition)
                                {
                                    std::cout << "At (" << intx << ","
                                    << inty << "," << intz << ") value="
                                    << spline_value;
                                    std::cout.flush();
                                }
#endif

                                // Add at that position the corresponding spline value
                                A3D_ELEM(*vol_voxels, iz, iy, ix) +=
                                    A3D_ELEM(vol_splines, k, i, j) * spline_value;
#ifdef DEBUG_MORE
                                if (condition)
                                {
                                    std::cout << " adding " << A3D_ELEM(vol_splines, k, i, j)
                                    << " * " << value_spline << " = "
                                    << A3D_ELEM(vol_splines, k, i, j)*
                                    value_spline << std::endl;
                                    std::cout.flush();
                                }
#endif
                            }
                }

                // Prepare for next iteration
                XX(act_coord) = XX(act_coord) + grid.relative_size * grid.basis(0, 0);
                YY(act_coord) = YY(act_coord) + grid.relative_size * grid.basis(1, 0);
                ZZ(act_coord) = ZZ(act_coord) + grid.relative_size * grid.basis(2, 0);
            }
            XX(beginY) = XX(beginY) + grid.relative_size * grid.basis(0, 1);
            YY(beginY) = YY(beginY) + grid.relative_size * grid.basis(1, 1);
            ZZ(beginY) = ZZ(beginY) + grid.relative_size * grid.basis(2, 1);
        }
        XX(beginZ) = XX(beginZ) + grid.relative_size * grid.basis(0, 2);
        YY(beginZ) = YY(beginZ) + grid.relative_size * grid.basis(1, 2);
        ZZ(beginZ) = ZZ(beginZ) + grid.relative_size * grid.basis(2, 2);
    }
}
#undef x0
#undef y0
#undef z0
#undef xF
#undef yF
#undef zF
#undef DEBUG
#undef DEBUG_MORE

/* Splines -> Voxels for a Grid -------------------------------------------- */
//#define DEBUG
void spatial_Bspline032voxels(const GridVolume &vol_splines,
                              MultidimArray<double> *vol_voxels, int Zdim, int Ydim, int Xdim)
{
    // Resize and set starting corner .......................................
    if (Zdim == 0 || Ydim == 0 || Xdim == 0)
    {
        Matrix1D<double> size = vol_splines.grid(0).highest -
                                vol_splines.grid(0).lowest;
        Matrix1D<double> corner = vol_splines.grid(0).lowest;
        (*vol_voxels).initZeros(CEIL(ZZ(size)), CEIL(YY(size)), CEIL(XX(size)));
        STARTINGX(*vol_voxels) = FLOOR(XX(corner));
        STARTINGY(*vol_voxels) = FLOOR(YY(corner));
        STARTINGZ(*vol_voxels) = FLOOR(ZZ(corner));
    }
    else
    {
        (*vol_voxels).initZeros(Zdim, Ydim, Xdim);
        (*vol_voxels).setXmippOrigin();
    }

    // Convert each subvolume ...............................................
    for (size_t i = 0; i < vol_splines.VolumesNo(); i++)
    {
        spatial_Bspline032voxels_SimpleGrid(vol_splines(i)(), vol_splines.grid(i),
                                            vol_voxels);
#ifdef DEBUG
        std::cout << "Spline grid no " << i << " stats: ";
        vol_splines(i)().printStats();
        std::cout << std::endl;
        std::cout << "So far vol stats: ";
        (*vol_voxels).printStats();
        std::cout << std::endl;
        VolumeXmipp save;
        save() = *vol_voxels;
        save.write((std::string)"PPPvoxels" + integerToString(i));
#endif
    }

    // Now normalise the resulting volume ..................................
    double inorm = 1.0 / sum_spatial_Bspline03_Grid(vol_splines.grid());
    FOR_ALL_ELEMENTS_IN_ARRAY3D(*vol_voxels)
    A3D_ELEM(*vol_voxels, k, i, j) *= inorm;

    // Set voxels outside interest region to minimum value .................
    double R = vol_splines.grid(0).get_interest_radius();
    if (R != -1)
    {
        double R2 = (R - 6) * (R - 6);

        // Compute minimum value within sphere
        double min_val = A3D_ELEM(*vol_voxels, 0, 0, 0);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(*vol_voxels)
        if (j*j + i*i + k*k <= R2 - 4)
            min_val = XMIPP_MIN(min_val, A3D_ELEM(*vol_voxels, k, i, j));

        // Substitute minimum value
        R2 = (R - 2) * (R - 2);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(*vol_voxels)
        if (j*j + i*i + k*k >= R2)
            A3D_ELEM(*vol_voxels, k, i, j) = min_val;
    }
}
#undef DEBUG

