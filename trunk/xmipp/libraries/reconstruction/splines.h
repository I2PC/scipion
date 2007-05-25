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

#ifndef _SPLINES_HH
#define _SPLINES_HH

#include "grids.h"

#include <external/bilib/headers/kernel.h>

/* ========================================================================= */
/* SPLINES                                                                   */
/* ========================================================================= */
/**@name Splines */
//@{
// Spline value ------------------------------------------------------------
/** Spline value.
    This function returns the value of a Bspline of order 3 at a given
    point */
inline double spatial_Bspline03(const Matrix1D<double> &r)
{
    if (-2 <= XX(r) && XX(r) < 2 &&
        -2 <= YY(r) && YY(r) < 2 &&
        -2 <= ZZ(r) && ZZ(r) < 2)
        return Bspline03(XX(r))*Bspline03(YY(r))*Bspline03(ZZ(r));
    else return 0.0;
}

const int BSPLINE03_SUBSAMPLING = 2000;
/** Value of a Bspline of order 3 in a Look-Up Table. */
double Bspline03LUT(double x);

/** Spline value with Look-Up Table.
    This function returns the value of a Bspline of order 3 at a given
    point */
inline double spatial_Bspline03LUT(const Matrix1D<double> &r)
{
    if (-2 <= XX(r) && XX(r) < 2 &&
        -2 <= YY(r) && YY(r) < 2 &&
        -2 <= ZZ(r) && ZZ(r) < 2)
        return Bspline03LUT(XX(r))*Bspline03LUT(YY(r))*Bspline03LUT(ZZ(r));
    else return 0.0;
}

/** Sum of a single spline over a grid.
    As a normalisation factor, the sum of the splines values over a given grid
    is needed. This function puts a spline at coordinate (0,0,0)
    and sums all the splines values at points of the grid which are inside the
    spline. It doesn't matter if the grid is compound or not.
*/
double sum_spatial_Bspline03_Grid(const Grid &grid);

// Blob projection ---------------------------------------------------------
/** Spline projection.
    This function returns the value of the spline line integral through a
    straight line which passes through the point r with direction u*/
double spatial_Bspline03_proj(const Matrix1D<double> &r, const Matrix1D<double> &u);

/** Splines ---> Voxels.
    The voxel size is defined in the same grid as the spline volume.

   However, you might give a size, usually set to 0, i.e., no external size.
   If no size is provided a size is produced such that all spline centers
   fit into the output volume.
*/
void spatial_Bspline032voxels(const GridVolume &vol_splines,
                              matrix3D<double> *vol_voxels,
                              int Zdim = 0, int Ydim = 0, int Xdim = 0);
//@}
#endif
