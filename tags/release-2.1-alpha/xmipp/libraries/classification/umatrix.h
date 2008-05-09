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
// xmippUmatrix.hh
//-----------------------------------------------------------------------------

#ifndef XMIPPUMATRIX_H
#define XMIPPUMATRIX_H

#include <vector>
#include <stdexcept>

#include "distance.h"
#include "map.h"

//-----------------------------------------------------------------------------
/// xmippUmatrix: Unified Distance Matrix (Umatrix)
//-----------------------------------------------------------------------------

/**@defgroup Umatrix Umatrix
   @ingroup ClassificationLibrary */
//@{

/** Alias for Umatrix Type (it uses a SOM map)
*/
typedef xmippMap umatType;

/**
 * This class implements the Unified Distance Matrix algorithm.
 * It is used for visualizing the clustering structure of Self-Organizing
 * Maps
 */
class xmippUmatrix
{
public:
    typedef xmippMap In;
    typedef xmippMap Out;
    typedef enum { NONE = 0, MEDIAN = 1, AVERAGE = 2 } smoothMode;

    /** Umatrix constructor
    *   Parameter: _sm: Defines the smoothing mode: NONE, MEDIAN or AVERAGE
    *
    */
    xmippUmatrix(smoothMode _sm = NONE): smooth(_sm)
    {};

    /**
    *   () operator.
    *    It takes a SOM as input and returns the calculated Umatrix
    */
    void operator()(const In& in, Out& out) const;

    /**
    *   Same as () operator.
    *   It takes a SOM as input and returns the calculated Umatrix
    */
    void getUmatrix(const In& in, Out& out) const;

    /**
    *   Gets the Umatrix but removing a variable from the analysis.
    *   It takes a SOM as input and returns the calculated Umatrix
    */
    void getUmatrix(const In& in, Out& out, const unsigned& _varOut) const;

    /**
    *   Gets the Umatrix but removing a list of variables from the analysis.
    *   It takes a SOM as input and returns the calculated Umatrix
    */

    void getUmatrix(const In& in, Out& out, const std::vector<unsigned>& _varsOut) const;


private:

    smoothMode smooth;          // Type of smoothing filter

    // Smoothing functions
    void medianSmoothing(Out& out) const;
    void averageSmoothing(Out& out) const;

};

//@}
#endif//XMIPPUMATRIX_H
