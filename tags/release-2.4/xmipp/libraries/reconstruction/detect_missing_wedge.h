/***************************************************************************
*
* Authors:    Carlos Oscar            coss@cnb.csic.es (2009)
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

#ifndef DETECT_MISSING_WEDGE_H
#define DETECT_MISSING_WEDGE_H

#include <data/matrix3d.h>
#include <data/volume.h>

/// @defgroup DetectMissingWedge Detect missing wedge
/// @ingroup ReconsLibraryPrograms

/// Parameters for the program detecting the missing wedge
/// @ingroup DetectMissingWedge
class DetectMissingWedge_parameters
{
public:
    /// Input volume
    FileName fn_vol;
    
    /// Plane width
    double planeWidth;

    /// Maximum frequency of the plane
    double maxFreq;

    /// Save marks
    bool saveMarks;
    
    /// Save mask
    bool saveMask;
public:
    // Input volume
    VolumeXmipp *V;

    // Magnitude of the input volume
    Matrix3D<double> *Vmag;
    
    // Angles of the first plane
    double rotPos, tiltPos;

    // Angles of the first plane
    double rotNeg, tiltNeg;

public:
    /// Read parameters from command line
    void read(int argc, char** argv);

    /// Produce side info.
    void produceSideInfo();

    /// Show parameters
    void show() const;

    /// Usage
    void usage() const;

    /// Run
    void run();
};
///@}
#endif
