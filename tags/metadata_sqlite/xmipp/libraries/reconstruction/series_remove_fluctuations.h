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

#ifndef SERIES_REMOVE_FLUCTUATIONS_H
#define SERIES_REMOVE_FLUCTUATIONS_H

#include <data/matrix3d.h>
#include <data/selfile.h>

/// @defgroup RemoveFluctuations Remove fluctuations in tilt series
/// @ingroup ReconsLibraryPrograms

/// Parameters for the program removing the fluctuations
/// @ingroup RemoveFluctuations
class Series_remove_fluctuations_parameters
{
public:
    /// Input images
    FileName fn_sel;
    
    /// Rootname for the output
    FileName fn_root;
    
    /// Cutoff frequency of the lowpass filter (<0.5)
    double maxFreq;
public:
    // Selfile with the input mages
    SelFile SF;
    
    // Volume for holding the whole tilt series
    Matrix3D<double> V;

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
#endif
