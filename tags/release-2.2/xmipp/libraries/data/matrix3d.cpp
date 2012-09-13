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

/* ------------------------------------------------------------------------- */
/* MATRICES 3D                                                               */
/* ------------------------------------------------------------------------- */
#include "matrix3d.h"

/* Interpolated element ---------------------------------------------------- */
template <>
std::complex<double> Matrix3D< std::complex<double> >::interpolatedElement(
    double x, double y, double z, std::complex<double> outside_value)
{
    int x0 = FLOOR(x);
    double fx = x - x0;
    int x1 = x0 + 1;
    int y0 = FLOOR(y);
    double fy = y - y0;
    int y1 = y0 + 1;
    int z0 = FLOOR(z);
    double fz = z - z0;
    int z1 = z0 + 1;

    std::complex<double> d000 = outside(z0, y0, x0) ? outside_value : dVkij(*this, z0, y0, x0);
    std::complex<double> d001 = outside(z0, y0, x1) ? outside_value : dVkij(*this, z0, y0, x1);
    std::complex<double> d010 = outside(z0, y1, x0) ? outside_value : dVkij(*this, z0, y1, x0);
    std::complex<double> d011 = outside(z0, y1, x1) ? outside_value : dVkij(*this, z0, y1, x1);
    std::complex<double> d100 = outside(z1, y0, x0) ? outside_value : dVkij(*this, z1, y0, x0);
    std::complex<double> d101 = outside(z1, y0, x1) ? outside_value : dVkij(*this, z1, y0, x1);
    std::complex<double> d110 = outside(z1, y1, x0) ? outside_value : dVkij(*this, z1, y1, x0);
    std::complex<double> d111 = outside(z1, y1, x1) ? outside_value : dVkij(*this, z1, y1, x1);

    std::complex<double> dx00 = LIN_INTERP(fx, d000, d001);
    std::complex<double> dx01 = LIN_INTERP(fx, d100, d101);
    std::complex<double> dx10 = LIN_INTERP(fx, d010, d011);
    std::complex<double> dx11 = LIN_INTERP(fx, d110, d111);
    std::complex<double> dxy0 = LIN_INTERP(fy, dx00, dx10);
    std::complex<double> dxy1 = LIN_INTERP(fy, dx01, dx11);

    return LIN_INTERP(fz, dxy0, dxy1);
}

// Apply geom --------------------------------------------------------------
template <>
void applyGeometryBSpline(Matrix3D< std::complex<double> > &M2,
                        Matrix2D<double> A, const Matrix3D< std::complex<double> > &M1,
                        int Splinedegree, bool inv, bool wrap, std::complex<double> outside)
{
    REPORT_ERROR(1, "applyGeometryBSpline: Not yet implemented for complex matrices\n");
}
