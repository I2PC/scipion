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

#include "multidimensional_array.h"



// Special case for complex numbers
template <>
void applyGeometryBSpline2D(Matrix2D< std::complex<double> > &M2,
                            const Matrix2D<double> &A, const Matrix2D< std::complex<double> > &M1,
                            int Splinedegree, bool inv, bool wrap, std::complex<double> outside,
                            unsigned long n)
{
    Matrix2D<double> re, im, rotre, rotim;
    MultidimArray<std::complex<double> > oneImg;
    double outre, outim;
    re.resize(YSIZE(M1), XSIZE(M1));
    im.resize(YSIZE(M1), XSIZE(M1));
    outre = outside.real();
    outim = outside.imag();
    M1.getImage(n, oneImg);
    Complex2RealImag(MULTIDIM_ARRAY(oneImg),
                     MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                     MULTIDIM_SIZE(oneImg));
    applyGeometryBSpline(rotre, A, re, Splinedegree, inv, wrap, outre);
    applyGeometryBSpline(rotim, A, im, Splinedegree, inv, wrap, outim);
    M2.resize(oneImg);
    RealImag2Complex(MULTIDIM_ARRAY(rotre), MULTIDIM_ARRAY(rotim),
                     MULTIDIM_ARRAY(M2), MULTIDIM_SIZE(re));
}

