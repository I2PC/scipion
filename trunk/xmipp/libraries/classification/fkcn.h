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
// xmippFKCN.hh
// Fuzzy Kohonen Clustering Network algorithm
//-----------------------------------------------------------------------------

#ifndef _XMIPPFKCN_H
#define _XMIPPFKCN_H

#pragma warning(disable:4786)

#include <iostream>
#include <vector>

#include "fcmeans.h"

/**@name Fuzzy Kohonen Clustering Network algorithm*/
//@{
/** This class implements Fuzzy Kohonen Clustering Network Algorithm (Bezdeck)
 *  an unsupervised clustering algorithm that is a combination of Fuzzy c-means
 *  and Kohonen SOM.
*/
class xmippFKCN:  public xmippFCMeans
{

public:


    /**  Big mega ctor. Creates a Fuzzy codebook, and initializes it
     * @param _m   Fuzzy constant
     * @param _epsilon  Stopping criterion
     * @param _epochs Number of epochs or iterations
    */
    xmippFKCN(double _m, double _epsilon, unsigned _epochs)
            : xmippFCMeans(_m, _epsilon, _epochs)
    {
        setID() = "Fuzzy Kohonen Clustering Network";
    };


    /// Virtual destructor
    virtual ~xmippFKCN()
    {};


    /**
     * Trains the Algorithm
     * @param _xmippDS Data structure to train, a codeBook in this case
     * @param _examples  A training set with the training examples
     */
    virtual void train(xmippFCB& _xmippDS,
                       const TS& _examples) const;

};

//@}

#endif//_XMIPPFKCN_H
