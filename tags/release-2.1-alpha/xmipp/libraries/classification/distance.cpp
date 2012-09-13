/***************************************************************************
 *
 * Authors:     Roberto Marabini
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
// xmippDistances.cpp
//-----------------------------------------------------------------------------
#pragma warning(disable:4786)

#include "distance.h"

//-----------------------------------------------------------------------------

/**
 * Euclidean distance.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
xmippFeature eDist(const xmippVector& _v1, const xmippVector& _v2)
{
    if (_v1.size() != _v2.size())
        throw std::runtime_error("vector of different size in eDist");

    double dist = 0;
    xmippVector::const_iterator i, j;
    for (i = _v1.begin(), j = _v2.begin() ; i < _v1.end(); i++, j++)
    {
        double tmp = (double)(*i) - (double)(*j);
        dist += tmp * tmp;
    }

    return (xmippFeature) sqrt(dist);
};

//-----------------------------------------------------------------------------

/**
 * Manhattan distance.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
xmippFeature mDist(const xmippVector& _v1, const xmippVector& _v2)
{
    if (_v1.size() != _v2.size())
        throw std::runtime_error("vector of different size in mDist");

    double dist = 0;
    xmippVector::const_iterator i, j;
    for (i = _v1.begin(), j = _v2.begin() ; i < _v1.end(); i++, j++)
        dist += fabs(((double)(*i) - (double)(*j)));

    return (xmippFeature) dist;
};

//-----------------------------------------------------------------------------
// Norm: norm of a vector (euclidean distance to origin)
//-----------------------------------------------------------------------------

xmippFeature xmippNorm::operator()(const xmippVector& v)
{
    double sum = 0.0;
    for (xmippVector::const_iterator i = v.begin(); i != v.end(); i++)
        sum += (double)(*i) * (double)(*i);
    return (xmippFeature) sqrt(sum);
}

//-----------------------------------------------------------------------------

