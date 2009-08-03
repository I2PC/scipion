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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef STEERABLE_H
#define STEERABLE_H

#include "matrix3d.h"
#include <vector>

/// @defgroup Steerable Steerable filters
/// @ingroup DataLibrary
//@{

/** Class for performing steerable filters */
class Steerable
{
public:
    // Basis functions for the steerability
    std::vector< Matrix3D<double> > basis;

public:
    /** Constructor.
       Sigma controls the width of the filter,
       deltaAng controls the accuracy of the final filtering.
       Vtomograph is the volume to filter.
       filterType is wall or filament. */
    Steerable(double sigma, Matrix3D<double> &Vtomograph, 
        double deltaAng, const std::string &filterType);
    
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
        const Matrix1D<double> &hx, const Matrix1D<double> &hy, 
        const Matrix1D<double> &hz, Matrix3D<double> &Vout);
};

/** Class with the parameters for detect structures */
class Prog_Detect_Structures_Param {
public:
    /// Volume to filter
    FileName fnIn;
    
    /// Volume to filter
    FileName fnOut;
    
    /// Filter type
    std::string filterType;
    
    /// Initial width
    double sigma0;
    
    /// Final width
    double sigmaF;
    
    /// Width step
    double sigmaStep;
    
    /// Angular step
    double angStep;

    /// Remove background
    bool removeBackground;
public:
    /// Read parameters from command line
    void read(int argc, char **argv);
    
    /// Show parameters
    void show() const;
    
    /// Usage
    void usage() const;
    
    /// Run
    void run();
};
///@}
#endif
