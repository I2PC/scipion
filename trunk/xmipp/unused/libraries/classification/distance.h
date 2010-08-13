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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

//-----------------------------------------------------------------------------
// xmippDistances.hh
//-----------------------------------------------------------------------------

#ifndef XMIPP_DISTANCES_H
#define XMIPP_DISTANCES_H

#pragma warning(disable:4786)

#include <cmath>              // sqrt, fabs
#include <stdexcept>           // runtime_error
#include <numeric>             // accumulate

#include "data_types.h"
#include "vector_ops.h"

/**@defgroup Distances Distances
   @ingroup ClassificationLibrary */
//@{
   
/**
 * Euclidean distance.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
Feature eDist(const FeatureVector& _v1, const FeatureVector& _v2);

//-----------------------------------------------------------------------------

/**
 * Manhattan distance.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
Feature mDist(const FeatureVector& _v1, const FeatureVector& _v2);

//-----------------------------------------------------------------------------

/**
 * This class defines an abstract interface for distances
 * Type T must have all arithmetic operations, specially - and conversion
 * to Feature
 */
class Distance
{
public:
    Distance()
    {}
    virtual ~Distance()
    {}

    /// Returns the distance between the two vectors
    virtual Feature operator()(const FeatureVector& _v1,
                                    const FeatureVector& _v2) = 0;
};

//-----------------------------------------------------------------------------

///Euclidean distance class
class EuclideanDistance: public Distance
{
public:
    EuclideanDistance(): Distance()
    {}
    ~EuclideanDistance()
    {}

    /**
     * Euclidean distance.
     * Parameter: _v1  First argument
     * Parameter: _v2  Second argument
     * @exception DifferentSize if _v1 and _v2  hasn't the same size
     */
    Feature operator()(const FeatureVector& _v1, const FeatureVector& _v2)
    {
        return eDist(_v1, _v2);
    }
};

//-----------------------------------------------------------------------------

///Manhattan distance
class ManhattanDistance: public Distance
{
public:
    ManhattanDistance(): Distance()
    {}
    ~ManhattanDistance()
    {}

    /**
     * Manhattan distance.
     * Parameter: _v1  First argument
     * Parameter: _v2  Second argument
     * @exception DifferentSize if _v1 and _v2  hasn't the same size
     */
    Feature operator()(const FeatureVector& _v1, const FeatureVector& _v2)
    {
        return mDist(_v1, _v2);
    }
};

//@}
#endif
