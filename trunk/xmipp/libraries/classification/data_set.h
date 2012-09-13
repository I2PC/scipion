/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
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

#ifndef XMIPPCDSET_H
#define XMIPPCDSET_H

#include "training_set.h"

/**@defgroup DataSets Data sets
   @ingroup ClassificationLibrary */
//@{
/**
 * Abstract Data Set class that should be used by all classification algorithms.
 * This is the parent class for all algorithms that can be trained and applied
 * in order to classify a set of data.
 * A xmipp classification object must have:
 * \\1) A train method that accepts a set of feature vectors
 * \\2) An apply method that accepts an example and classifies it
 */

template<class InClass, class OutClass>
class ClassificationDataSet
{
public:
    /// Class of input vectors. Usually a FeatureVector (vector of Feature)
    typedef InClass In;

    /// Class of the target. Can be a double, string, unsigned, even a vector ...
    typedef OutClass Out;

    /// Training set. Set of vectors (training vectors), probably classified.
    typedef ClassificationTrainingSet<In, Out> TS;

    /**
     * Constructor.
     * This constructor is empty.
     */
    ClassificationDataSet()
    {};

    /**
     * Destructor.
     * The default destructor
     */
    virtual ~ClassificationDataSet()
    {};

    /**
     * Method to classify a feature vector
     * It returns the 'class' to which the vector belongs
     * Parameter: _in  vcetor to test.
     * @return     The result of classification.
     */
    virtual Out apply(const In& _in) const = 0;

    /**
     * Method to classify an input vector.
     * This method returns an unsigned integer that would correspond to
     * the output neuron, the output codevector, or anything similar.
     * If it makes no sense, it should be declared as private. Although
     * it means many different things, it�s included here to have an uniform
     * representation
     * Parameter: _in  input vector to test.
     * @return     The result of classification.
     */
    virtual unsigned output(const In& _in) const = 0;
};

//@}
#endif//XMIPPCDSET_H
