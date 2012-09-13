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

#ifndef ROTATIONALSPECTRUM_H
#define ROTATIONALSPECTRUM_H

#include "matrix2d.h"

#include <iostream>

/// @defgroup RotationalSpectrum Rotational spectrum
/// @ingroup DataLibrary

/** Cylindrical_Wave_Decomposition class.
 * @ingroup RotationalSpectrum
 */
class Cylindrical_Wave_Decomposition
{
public:
    /// Minimum harmonics
    int numin;

    /// Maximum harmonics
    int numax;

    /// Center of symmetry (x)
    double x0;

    /// Center of symmetry (x)
    double y0;

    /// Minimum integration radius
    double r1;

    /// Maximum integration radius
    double r2;

    /// Integration increment
    double r3;

    ///  Ir
    int ir;

    ///  Ampcos
    Matrix1D< double > out_ampcos;

    ///  Ampsin
    Matrix1D< double > out_ampsin;

    /// Show this object
    friend std::ostream& operator<<(std::ostream& _out,
                                    const Cylindrical_Wave_Decomposition& _cwd);

    /// Interpolate image value (bilinear)
    double interpolate(Matrix2D< double >& img, double y, double x);

    /// Compute the Cylindrical Wave decomposition of an image
    void compute_cwd(Matrix2D< double >& img);
};

/** Rotational spectrum.
 * @ingroup RotationalSpectrum
 *
 * Example of use:
 *
 * @code
 * int main(int argc, char** argv)
 * {
 *     ImageXmipp I(argv[1]);
 *
 *     int rl = 0;
 *     int rh = 22;
 *     int dr = 1;
 *
 *     Rotational_Spectrum spt;
 *     spt.rl = rl;
 *     spt.rh = rh;
 *     spt.dr = dr;
 *     spt.numin = 1;
 *     spt.numax = 15;
 *     spt.x0 = (double) XSIZE(I()) / 2;
 *     spt.y0 = (double) YSIZE(I()) / 2;
 *
 *     spt.compute_rotational_spectrum (I(), rl, rh, dr, rh - rl);
 *
 *     std::cout << spt.rot_spectrum << std::endl;
 *
 *     return 0;
 * }
 * @endcode
 */
class Rotational_Spectrum
{
public:
    /** Ir.
     *
     * This value is given by the Cylindrical Wave Decomposition. It is not
     * necessary to fill this field if the rotational spectrum is computed on
     * an image.
     */
    int ir;

    /// Minimum harmonics.
    int numin;

    /// Maximum harmonics.
    int numax;

    /// Center of symmetry (x).
    double x0;

    /// Center of symmetry (y).
    double y0;

    /// Minimum integration radius.
    double rl;

    /// Maximum integration radius.
    double rh;

    /// Integration increment.
    double dr;

    ///  Rotational spectrum.
    Matrix1D< double > rot_spectrum;

    /// Show.
    friend std::ostream& operator<<(std::ostream& _out,
                                    const Rotational_Spectrum& _spt);

    /** Compute rotational spectrum of an image.
     *
     * xr1 is the minimum integration radius. xr2 the maximum integration
     * radius. xdr the increment, and xr the length of integration. Usually,
     * xr=xr2-xr1.
     */
    void compute_rotational_spectrum(Matrix2D< double >& img,
                                     double xr1,
                                     double xr2,
                                     double xdr,
                                     double xr);

    /// Compute rotational spectrum using CWD.
    void compute_rotational_spectrum(Cylindrical_Wave_Decomposition& cwd,
                                     double xr1,
                                     double xr2,
                                     double xdr,
                                     double xr);

    /// Read parameters from command line.
    void read(int argc, char* *argv);

    /// Usage.
    void usage();
};

#endif
