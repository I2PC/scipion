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
// xmippDistances.hh
//-----------------------------------------------------------------------------

#ifndef XMIPP_DISTANCES_H
#define XMIPP_DISTANCES_H

#pragma warning(disable:4786)

//-----------------------------------------------------------------------------

#include <math.h>              // sqrt, fabs
#include <stdexcept>           // runtime_error
#include <numeric>             // accumulate
#include "xmippCDataTypes.hh"
#include "xmippVectorOps.hh"

//-----------------------------------------------------------------------------

/**
 * Euclidean distance.
 * @param _v1  First argument
 * @param _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
xmippFeature eDist(const xmippVector& _v1, const xmippVector& _v2);

//-----------------------------------------------------------------------------

/**
 * Manhattan distance.
 * @param _v1  First argument
 * @param _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
xmippFeature mDist(const xmippVector& _v1, const xmippVector& _v2);

//-----------------------------------------------------------------------------

/**@name Distance class*/
//@{
/** 
 * This class defines an abstract interface for distances
 * Type T must have all arithmetic operations, specially - and conversion
 * to xmippFeature
 */
class xmippDistance {
 public:
  xmippDistance() {}
  virtual ~xmippDistance(){}
  
  /// Returns the distance between the two vectors
  virtual xmippFeature operator()(const xmippVector& _v1, 
				const xmippVector& _v2) = 0;
};
//@}

//-----------------------------------------------------------------------------

/**@name Euclidean Distance*/
//@{
///Euclidean distance class
class xmippEDistance: public xmippDistance {
 public:
  xmippEDistance(): xmippDistance() {}
  ~xmippEDistance() {}

  /**
   * Euclidean distance.
   * @param _v1  First argument
   * @param _v2  Second argument
   * @exception DifferentSize if _v1 and _v2  hasn't the same size
   */
  xmippFeature operator()(const xmippVector& _v1, const xmippVector& _v2) {
    return eDist( _v1, _v2 );
  }
};
//@}

//-----------------------------------------------------------------------------

/**@name Manhattan*/
//@{
///Manhattan distance
class xmippMDistance: public xmippDistance {
 public:
  xmippMDistance(): xmippDistance() {}
  ~xmippMDistance() {}

  /**
   * Manhattan distance.
   * @param _v1  First argument
   * @param _v2  Second argument
   * @exception DifferentSize if _v1 and _v2  hasn't the same size
   */
  xmippFeature operator()(const xmippVector& _v1, const xmippVector& _v2) {
    return mDist( _v1, _v2 );
  }
};
//@}


/**@name Norm class*/
//@{
/**
* Norm: norm of a vector (euclidean distance to origin)
*
*/
class xmippNorm
{
 public:
  xmippFeature operator()(const xmippVector& v);
};
//@}

//-----------------------------------------------------------------------------

#endif
