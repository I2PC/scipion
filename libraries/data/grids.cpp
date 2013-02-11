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
#include "grids.h"
#include "numerical_tools.h"

#include <stdio.h>
#include <string.h> // for memcpy

/*****************************************************************************/
/* Simple Grids                                                              */
/*****************************************************************************/
// Constructor -------------------------------------------------------------
SimpleGrid::SimpleGrid()
{
    basis.initIdentity(3);
    inv_basis.initIdentity(3);
    origin        = vectorR3(0., 0., 0.);
    lowest        = vectorR3(-5., -5., -5.);
    highest       = -lowest;
    relative_size = 1;
    R2            = -1;
}

SimpleGrid::SimpleGrid(const SimpleGrid &SG)
{
    basis         = SG.basis;
    inv_basis     = SG.inv_basis;
    lowest        = SG.lowest;
    highest       = SG.highest;
    relative_size = SG.relative_size;
    origin        = SG.origin;
    R2            = SG.R2;
}

// std::cout --------------------------------------------------------------------
std::ostream& operator <<(std::ostream& o, const SimpleGrid &grid)
{
    o << "   Simple Grid -----" << std::endl;
    Matrix1D<double> aux;
    (grid.basis).getCol(0,aux); o << "   Vector 1: " << aux.transpose() << std::endl;
    (grid.basis).getCol(1,aux); o << "   Vector 2: " << aux.transpose() << std::endl;
    (grid.basis).getCol(2,aux); o << "   Vector 3: " << aux.transpose() << std::endl;
    o << "   Relative size:         " << grid.relative_size       << std::endl;
    o << "   Interest radius:       " << grid.R2                  << std::endl;
    o << "   Origin (univ.coords)   " << grid.origin.transpose()  << std::endl;
    o << "   Highest (grid. coord)  " << grid.highest.transpose() << std::endl;
    o << "   Lowest (grid. coord)   " << grid.lowest.transpose()  << std::endl;
    return o;
}

// Assignment --------------------------------------------------------------
SimpleGrid& SimpleGrid::operator = (const SimpleGrid &SG)
{
    if (&SG != this)
    {
        basis         = SG.basis;
        inv_basis     = SG.inv_basis;
        lowest        = SG.lowest;
        highest       = SG.highest;
        relative_size = SG.relative_size;
        R2            = SG.R2;
        origin        = SG.origin;
    }
    return *this;
}

// Another function for assignment -----------------------------------------
void SimpleGrid::assign(const SimpleGrid &SG)
{
    *this = SG;
}

// Number of samples -------------------------------------------------------
int SimpleGrid::get_number_of_samples() const
{
    if (R2 == -1)
    {
        int Zdim, Ydim, Xdim;
        getSize(Zdim, Ydim, Xdim);
        return Zdim*Ydim*Xdim;
    }
    else
    {
        int ZZ_lowest = (int) ZZ(lowest);
        int YY_lowest = (int) YY(lowest);
        int XX_lowest = (int) XX(lowest);
        int ZZ_highest = (int) ZZ(highest);
        int YY_highest = (int) YY(highest);
        int XX_highest = (int) XX(highest);
        Matrix1D<double> grid_index(3), univ_position(3);
        int N = 0;
        for (int k = ZZ_lowest; k <= ZZ_highest; k++)
            for (int i = YY_lowest; i <= YY_highest; i++)
                for (int j = XX_lowest; j <= XX_highest; j++)
                {
                    VECTOR_R3(grid_index, j, i, k);
                    grid2universe(grid_index, univ_position);
                    if (is_interesting(univ_position)) N++;
                }
        return N;
    }
}

// Prepare Grid ------------------------------------------------------------
void SimpleGrid::prepare_grid()
{
    // Compute matrix for inverse basis conversion
    try
    {
        inv_basis = basis.inv();
    }
    catch (XmippError &error)
    {
        REPORT_ERROR(ERR_UNCLASSIFIED, "The grid vectors are not a true 3D coordinate system");
    }
    lowest.setCol();
    highest.setCol();
    origin.setCol();
}

// Minimum size ------------------------------------------------------------
void Grid::voxel_corners(Matrix1D<double> &Gcorner1, Matrix1D<double> &Gcorner2,
                         const Matrix2D<double> *V) const
{
    Matrix1D<double> SGcorner1(3), SGcorner2(3);     // Subgrid corners
    SPEED_UP_temps012;

    // Look for the lowest and highest volume coordinate
    Gcorner1.resize(3);  // lowest and highest coord.
    Gcorner2.resize(3);
    for (size_t n = 0; n < GridsNo(); n++)
    {
        // Find box for this grid
        bool first;
        first = true;
        for (int k = (int)ZZ(LG[n].lowest); k <= ZZ(LG[n].highest); k++)
            for (int i = (int)YY(LG[n].lowest); i <= YY(LG[n].highest); i++)
                for (int j = (int)XX(LG[n].lowest); j <= XX(LG[n].highest); j++)
                {
                    Matrix1D<double> grid_index(3), univ_position(3);
                    VECTOR_R3(grid_index, j, i, k);
                    LG[n].grid2universe(grid_index, univ_position);
                    if (V != NULL)
                    {
                        M3x3_BY_V3x1(univ_position, *V, univ_position);
                    }
                    if (!LG[n].is_interesting(univ_position)) continue;
                    if (!first)
                    {
                        XX(SGcorner1) = XMIPP_MIN(XX(SGcorner1), XX(univ_position));
                        YY(SGcorner1) = XMIPP_MIN(YY(SGcorner1), YY(univ_position));
                        ZZ(SGcorner1) = XMIPP_MIN(ZZ(SGcorner1), ZZ(univ_position));

                        XX(SGcorner2) = XMIPP_MAX(XX(SGcorner2), XX(univ_position));
                        YY(SGcorner2) = XMIPP_MAX(YY(SGcorner2), YY(univ_position));
                        ZZ(SGcorner2) = XMIPP_MAX(ZZ(SGcorner2), ZZ(univ_position));
                    }
                    else
                    {
                        SGcorner2 = SGcorner1 = univ_position;
                        first = false;
                    }
                }

        // Compare with the rest of the grids
        if (n != 0)
        {
            XX(Gcorner1) = XMIPP_MIN(XX(Gcorner1), XX(SGcorner1));
            YY(Gcorner1) = XMIPP_MIN(YY(Gcorner1), YY(SGcorner1));
            ZZ(Gcorner1) = XMIPP_MIN(ZZ(Gcorner1), ZZ(SGcorner1));

            XX(Gcorner2) = XMIPP_MAX(XX(Gcorner2), XX(SGcorner2));
            YY(Gcorner2) = XMIPP_MAX(YY(Gcorner2), YY(SGcorner2));
            ZZ(Gcorner2) = XMIPP_MAX(ZZ(Gcorner2), ZZ(SGcorner2));
        }
        else
        {
            Gcorner1 = SGcorner1;
            Gcorner2 = SGcorner2;
        }

#ifdef DEBUG
        std::cout << LG[n];
        std::cout << "SGcorner1 " << SGcorner1.transpose() << std::endl;
        std::cout << "SGcorner2 " << SGcorner2.transpose() << std::endl;
        std::cout << "Gcorner1  " << Gcorner1.transpose() << std::endl;
        std::cout << "Gcorner2  " << Gcorner2.transpose() << std::endl;
#endif
    }
}

/*****************************************************************************/
/* Some useful Grids                                                         */
/*****************************************************************************/
/* Create CC Simple grid with a given origin ------------------------------- */
SimpleGrid Create_CC_grid(double relative_size, const Matrix1D<double> &corner1,
                          const Matrix1D<double> &corner2, const Matrix1D<double> &origin)
{
    SimpleGrid    grid;

    // The vectors of the grid are the default ones of (1,0,0), (0,1,0),
    // and (0,0,1), and its inverse matrix is already computed
    grid.relative_size = relative_size;
    grid.origin        = origin;

    // Compute the lowest and highest indexes inside the grid
    grid.universe2grid(corner1, grid.lowest);
    grid.lowest.selfFLOOR();
    grid.universe2grid(corner2, grid.highest);
    grid.highest.selfCEIL();

    grid.R2 = -1;
    return grid;
}

/* Create CC grid ---------------------------------------------------------- */
Grid Create_CC_grid(double relative_size, const Matrix1D<double> &corner1,
                    const Matrix1D<double> &corner2)
{
    Grid            result;
    SimpleGrid      aux_grid;

    Matrix1D<double> origin = (corner1 + corner2) / 2;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(origin)
        if (ABS(ROUND(origin(i))-origin(i))<0.45)
            origin(i)=ROUND(origin(i));
        else
            origin(i)=CEIL(origin(i));
    aux_grid = Create_CC_grid(relative_size, corner1, corner2, origin);
    result.add_grid(aux_grid);
    return result;
}

Grid Create_CC_grid(double relative_size, int Zdim, int Ydim, int Xdim)
{
    Grid            result;
    SimpleGrid      aux_grid;

    Matrix1D<double> origin =
        vectorR3((double)FLOOR(Xdim / 2.0), (double)FLOOR(Ydim / 2.0),
                  (double)FLOOR(Zdim / 2.0));
    aux_grid = Create_CC_grid(relative_size, -origin,
                              vectorR3((double)Xdim, (double)Ydim, (double)Zdim) - origin - 1,
                              vectorR3(0.0,0.0,0.0));
    result.add_grid(aux_grid);

    return result;
}

/* Create BCC grid --------------------------------------------------------- */
Grid Create_BCC_grid(double relative_size, const Matrix1D<double> &corner1,
                     const Matrix1D<double> &corner2)
{
    Grid             result;
    SimpleGrid       aux_grid;
    Matrix1D<double> origin = (corner1 + corner2) / 2;
    origin.selfROUND();

    //Even Slice
    //    0 1 2 3 4 5 6 7 8 9 10 11 12 (Col)
    //  0 A   A   A   A   A   A     A
    //  1
    //  2 A   A   A   A   A   A     A
    //  3
    //  4 A   A   A   A   A   A     A
    //  5
    //  6 A   A   A   A   A   A     A
    //  7
    //  8 A   A   A   A   A   A     A
    //  9
    // 10 A   A   A   A   A   A     A
    // 11
    // 12 A   A   A   A   A   A     A
    //(Row)
    //
    //Odd Slice
    //    0 1 2 3 4 5 6 7 8 9 10 11 12 (Col)
    //  0
    //  1   B   B   B   B   B    B
    //  2
    //  3   B   B   B   B   B    B
    //  4
    //  5   B   B   B   B   B    B
    //  6
    //  7   B   B   B   B   B    B
    //  8
    //  9   B   B   B   B   B    B
    // 10
    // 11   B   B   B   B   B    B
    // 12
    //(Row)

    // Grid A
    aux_grid = Create_CC_grid(relative_size, corner1, corner2, origin);
    result.add_grid(aux_grid);

    // Grid B
    origin = origin + relative_size / 2 * vectorR3(1., 1., 1.);
    aux_grid = Create_CC_grid(relative_size, corner1, corner2, origin);
    result.add_grid(aux_grid);

    return result;
}

/* Create FCC grid --------------------------------------------------------- */
Grid Create_FCC_grid(double relative_size, const Matrix1D<double> &corner1,
                     const Matrix1D<double> &corner2)
{

    Grid             result;
    SimpleGrid       aux_grid;
    Matrix1D<double> aux_origin;
    Matrix1D<double> cornerb;
    Matrix1D<double> origin = (corner1 + corner2) / 2;
    origin.selfROUND();

    //Even Slice
    //    0 1 2 3 4 5 6 7 8 9 10 11 12 (Col)
    //  0 A   A   A   A   A   A     A
    //  1   B   B   B   B   B    B
    //  2 A   A   A   A   A   A     A
    //  3   B   B   B   B   B    B
    //  4 A   A   A   A   A   A     A
    //  5   B   B   B   B   B    B
    //  6 A   A   A   A   A   A     A
    //  7   B   B   B   B   B    B
    //  8 A   A   A   A   A   A     A
    //  9   B   B   B   B   B    B
    // 10 A   A   A   A   A   A     A
    // 11   B   B   B   B   B    B
    // 12 A   A   A   A   A   A     A
    //(Row)
    //
    //Odd Slice
    //    0 1 2 3 4 5 6 7 8 9 10 11 12 (Col)
    //  0   C   C   C   C   C    C
    //  1 D   D   D   D   D   D     D
    //  2   C   C   C   C   C    C
    //  3 D   D   D   D   D   D     D
    //  4   C   C   C   C   C    C
    //  5 D   D   D   D   D   D     D
    //  6   C   C   C   C   C    C
    //  7 D   D   D   D   D   D     D
    //  8   C   C   C   C   C    C
    //  9 D   D   D   D   D   D     D
    // 10   C   C   C   C   C    C
    // 11 D   D   D   D   D   D     D
    // 12   C   C   C   C   C    C
    //(Row)

    // Grid A
    aux_grid = Create_CC_grid(relative_size, corner1, corner2, origin);
    result.add_grid(aux_grid);
    // Grid D
    aux_origin = origin + relative_size / 2 * vectorR3(0., 1., 1.);
    aux_grid = Create_CC_grid(relative_size, corner1, corner2, aux_origin);
    result.add_grid(aux_grid);
    // Grid C
    aux_origin = origin + relative_size / 2 * vectorR3(1., 0., 1.);
    aux_grid = Create_CC_grid(relative_size, corner1, corner2, aux_origin);
    result.add_grid(aux_grid);
    // Grid B
    cornerb = corner2;
    cornerb(0) = cornerb(0) - 1;
    cornerb(1) = cornerb(1) - 1;
    aux_origin = origin + relative_size / 2 * vectorR3(1., 1., 0.);
    aux_grid = Create_CC_grid(relative_size, corner1, cornerb, aux_origin);
    result.add_grid(aux_grid);

    return result;
}
#undef MULTIPLY_CC_GRID_BY_TWO

/* CC grid with region of interest ----------------------------------------- */
//#define DEBUG
SimpleGrid Create_grid_within_sphere(double relative_size,
                                     const Matrix1D<double> &origin,
                                     const Matrix1D<double> &X, const Matrix1D<double> &Y,
                                     const Matrix1D<double> &Z, double R2)
{
    SimpleGrid    grid;
    double R = sqrt(R2);

    grid.set_X(X);
    grid.set_Y(Y);
    grid.set_Z(Z);
    grid.relative_size = relative_size;
    grid.origin        = origin;
    grid.R2 = R2;
    grid.lowest.initZeros(3);
    grid.highest.initZeros(3);
    grid.prepare_grid();

    // Find grid limits
    int iR = CEIL(R);
    Matrix1D<double> univ_position(3), grid_position(3);
    for (int k = -iR; k <= iR; k++)
        for (int i = -iR; i <= iR; i++)
            for (int j = -iR; j <= iR; j++)
            {
                VECTOR_R3(univ_position, j, i, k);
                if (univ_position.module() > R) continue;
                grid.universe2grid(univ_position, grid_position);
                XX(grid.lowest) = XMIPP_MIN(XX(grid.lowest), FLOOR(XX(grid_position)));
                YY(grid.lowest) = XMIPP_MIN(YY(grid.lowest), FLOOR(YY(grid_position)));
                ZZ(grid.lowest) = XMIPP_MIN(ZZ(grid.lowest), FLOOR(ZZ(grid_position)));
                XX(grid.highest) = XMIPP_MAX(XX(grid.highest), CEIL(XX(grid_position)));
                YY(grid.highest) = XMIPP_MAX(YY(grid.highest), CEIL(YY(grid_position)));
                ZZ(grid.highest) = XMIPP_MAX(ZZ(grid.highest), CEIL(ZZ(grid_position)));
            }

#ifdef DEBUG
    std::cout << "Sphere radius = " << R << std::endl
    << "relative size = " << relative_size << std::endl
    << "X module      = " << X.module() << std::endl
    << "Y module      = " << Y.module() << std::endl
    << "Z module      = " << Z.module() << std::endl
    << grid
    ;
#endif
    return grid;
}
#undef DEBUG

/* Create CC grid ---------------------------------------------------------- */
Grid Create_CC_grid(double relative_size, double R)
{
    Grid            result;
    SimpleGrid      aux_grid;

    Matrix1D<double> origin(3);
    origin.initZeros();
    Matrix1D<double> x(3), y(3), z(3);
    VECTOR_R3(x, 1, 0, 0);
    VECTOR_R3(y, 0, 1, 0);
    VECTOR_R3(z, 0, 0, 1);
    aux_grid = Create_grid_within_sphere(relative_size, origin, x, y, z, R * R);
    result.add_grid(aux_grid);
    return result;
}

/* Create BCC grid --------------------------------------------------------- */
Grid Create_BCC_grid(double relative_size, double R)
{
    Grid            result;
    SimpleGrid      aux_grid;

    Matrix1D<double> origin(3);
    origin.initZeros();
    Matrix1D<double> x(3), y(3), z(3);
    VECTOR_R3(x, 0.5, 0.5, -0.5);
    VECTOR_R3(y, 0.5, -0.5, 0.5);
    VECTOR_R3(z, -0.5, 0.5, 0.5);
    aux_grid = Create_grid_within_sphere(relative_size, origin, x, y, z, R * R);
    result.add_grid(aux_grid);
    return result;
}

/* Create FCC grid --------------------------------------------------------- */
Grid Create_FCC_grid(double relative_size, double R)
{
    Grid            result;
    SimpleGrid      aux_grid;

    Matrix1D<double> origin(3);
    origin.initZeros();
    Matrix1D<double> x(3), y(3), z(3);
    VECTOR_R3(x, 0.5, 0.5, 0);
    VECTOR_R3(y, 0.5, 0, 0.5);
    VECTOR_R3(z, 0, 0.5, 0.5);
    aux_grid = Create_grid_within_sphere(relative_size, origin, x, y, z, R * R);
    result.add_grid(aux_grid);
    return result;
}
