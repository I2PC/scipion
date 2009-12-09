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

#ifndef DETECT_STRUCTURE_H
#define DETECT_STRUCTURE_H

#include <data/matrix3d.h>
#include <data/steerable.h>
#include <vector>

/// @defgroup DetectStructures Detect structures
/// @ingroup ReconsLibrary
//@{

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
    
    /// Remove missing wedge
    bool removeMissingWedge;

    /// Plane 1 rot
    double rot1;

    /// Plane 1 tilt
    double tilt1;

    /// Plane 2 rot
    double rot2;

    /// Plane 2 tilt
    double tilt2;
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
