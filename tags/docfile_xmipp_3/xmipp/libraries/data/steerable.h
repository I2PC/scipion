/***************************************************************************
 *
 * Authors:     Manuel Sanchez Pau 
 *              Carlos Oscar Sanchez Sorzano
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

#ifndef STEERABLE_H
#define STEERABLE_H

#include "matrix3d.h"
#include <vector>

/// @defgroup MissingWedge Missing wedge
/// @ingroup DataLibrary
///@{
class MissingWedge
{
public:
    /// Rot of the positive plane
    double rotPos;
    /// Tilt of the positive plane
    double tiltPos;
    /// Rot of the negative plane
    double rotNeg;
    /// Tilt of the negative plane
    double tiltNeg;
public:

    /// Empty constructor
    MissingWedge()
    {
        rotPos = 0.;
        tiltPos = 0.;
        rotNeg = 0.;
        tiltNeg = 0.;
    }

    /// Remove wedge
    void removeWedge(Matrix3D<double> &V) const;
};
///@}

/// @defgroup Steerable Steerable filters
/// @ingroup DataLibrary
///@{

/** Class for performing steerable filters */
class Steerable
{
public:
    // Basis functions for the steerability
    std::vector< Matrix3D<double> > basis;

    // Missing wedge
    const MissingWedge *MW;
public:
    /** Constructor.
       Sigma controls the width of the filter,
       deltaAng controls the accuracy of the final filtering.
       Vtomograph is the volume to filter.
       filterType is wall or filament. */
    Steerable(double sigma, Matrix3D<double> &Vtomograph, 
        double deltaAng, const std::string &filterType,
        const MissingWedge *_MW);
    
    /** This function is the one really filtering */
    void buildBasis(const Matrix3D<double> &Vtomograph, double sigma);

    /** Internal function for the generation of 1D filters. */
    void generate1DFilters(double sigma,
        const Matrix3D<double> &Vtomograph,
        std::vector< Matrix1D<double> > &hx,
        std::vector< Matrix1D<double> > &hy,
        std::vector< Matrix1D<double> > &hz);

    /** Internal function for the generation of 3D filters. */
    void generate3DFilter(Matrix3D<double>& h3D,
	std::vector< Matrix1D<double> > &hx,
	std::vector< Matrix1D<double> > &hy,
	std::vector< Matrix1D<double> > &hz);

    /** Internal function for filtering */
    void singleFilter(const Matrix3D<double>& Vin,
        Matrix1D<double> &hx, Matrix1D<double> &hy, 
        Matrix1D<double> &hz, Matrix3D<double> &Vout);
};
///@}
#endif
