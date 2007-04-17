/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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
// xmippBatchSOM.hh
// Implements Kohonen Self-Organizing Feature Maps by using Batch SOM training
//-----------------------------------------------------------------------------

#ifndef XMIPPBATCHSOM_H
#define XMIPPBATCHSOM_H

#include "xmippBaseAlgo.hh"
#include "xmippMap.hh"
#include "xmippSOM.hh"



//---------------------------------------------------------------------------


/**@name Batch training of Kohonen SOM algorithm*/
//@{
/**
 * This class trains a Kohonen's Self Organizing Map using Batch SOM.
 */
class xmippBatchSOM : public xmippSOM
{
 public:


  /**
   * Constructs the algorithm
   * @param _radius     How is gonna decrease the radius of neighborhood
   * @param _nSteps     Number of training steps
   */
  xmippBatchSOM(Descent& _radius,  unsigned long _nSteps)
  : xmippSOM(_radius, _radius, BUBBLE, _nSteps) {};

  /**
   * Construct a BatchSOM from the code vectors in a stream
   * @param _is  The stream
   */
  xmippBatchSOM(istream& _is);


  /**
   * Virtual destructor
   */
  virtual ~xmippBatchSOM() {};

  /**
   * Trains the SOM
   * @param _som  The som to train
   * @param _ts   The training set
   */
  virtual void train (xmippMap& _som, const xmippCTVectors& _ts) const;

  /**
   * Trains the SOM
   * @param _som  The som to train
   * @param _ts   The training set
   */
//  virtual void train (xmippMap& _som, const TS& _ts) const {};


};

//@}



#endif
