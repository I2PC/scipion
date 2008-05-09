/***************************************************************************
 *
 * Authors:     Alberto Pascual
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
//-----------------------------------------------------------------------------
// xmippPlanes.hh
//-----------------------------------------------------------------------------

#ifndef XMIPPPLANES_H
#define XMIPPPLANES_H

#include <stdexcept>

#include "map.h"

/**@defgroup Planes Planes
   @ingroup ClassificationLibrary */
//@{
/// Alias for Plane Type (it uses a SOM map)
typedef xmippMap planeType;

/**
 * This class calculates the influence of each independent variable
 * in the SOM. It is used for visualizing the effect of each vector component on
 * the Self-Organizing Map.
 */
class xmippPlanes
{
public:
    typedef xmippMap In;
    typedef planeType Out;

    /** Planes default constructor
    */
    xmippPlanes()
    {};

    /** getPlane: gets the "plane" (a map representing the effect of the
    *   given variable.
    *   Parameter: _in: input SOM
    *   Parameter: _out: output plane
    *   Parameter: _plane: variable
    */
    void getPlane(const In& _in, Out& _out, const unsigned _plane) const;


};
//@}

#endif//XMIPPPLANES_H
