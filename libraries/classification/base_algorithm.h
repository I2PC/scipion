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

//-----------------------------------------------------------------------------
// ClassificationAlgorithm.h
// Base class for all algorithms within Xmipp Library
//-----------------------------------------------------------------------------


#ifndef BASEALGORITHM_H
#define BASEALGORITHM_H

#include <iostream>

#include "training_set.h"
#include "data_set.h"

#include <data/xmipp_funcs.h>

/**@defgroup BaseAlgorithm Base Algorithm
   @ingroup ClassificationLibrary */
//@{
//-----------------------------------------------------------------------------
/**
 * This is the parent class for all algorithms that can be trained and applied
 * in order to classify a set of data.
 * A xmipp Algorithm object must have:
 * 1) a train method that accepts a set of examples
 * 2) an apply method that accepts an example and classifies it
 */
template<class DSClass>
class ClassificationAlgorithm
{
public:

    typedef DSClass DS;
    typedef typename DS::TS TS;

    /**
     * Constructor.
     * Parameter: _ID an ID string unique for each algorithm class
     */
    ClassificationAlgorithm(const std::string& _ID = ""): ID(_ID)
    {};

    /**
     * Destructor.
     * The default destructor
     */
    virtual ~ClassificationAlgorithm()
    {};

    /**
     * Method to train the algorithm.
     * Note that this method is pure virtual, so it ust be defined in every xmipp descendant class.
     * Parameter: _ds Data structure to train
     * Parameter: _examples  A training set with the training examples.
     */
    virtual void train(DS& _ds, const TS& _examples) const
        {};

    /**
     * Method to train the algorithm.
     * Note that this method is pure virtual, so it ust be defined in every xmipp descendant class.
     * Parameter: _ds Data structure to train
     * Parameter: _examples  A training set with the training examples.
     */
    virtual void train(DS& _ds, TS& _examples) const
        {};

    /**
     * Method to test the algorithm.
     * Note that this method is pure virtual, so it must be defined in every xmipp descendant class.
     * Parameter: _ds Data structure to train. Const because its not affected by test
     * Parameter: _examples  A training set with the training examples.
     * @return The test error, that is, error obtained for the test file. Might be
     * MSE error, or success in classification error, or whatever; just be consistent
     * for each algorithm
     */
    virtual double test(const DS& _ds, const TS& _examples) const = 0;

    /** Print itself on standard output
     */
    virtual void printSelf(std::ostream& _os) const
    {
        _os << "xmippBaseAlgo" << std::endl; // to identify it as an algorith
        _os << ID << std::endl;
    };

    /** Set ID (returns a const reference to the ID)
    */
    virtual const std::string& setID() const
    {
        return ID;
    };

    /** Set ID (returns a non-const reference to the ID)
    */
    virtual std::string& setID()
    {
        return ID;
    };


    /** Defines Listener class
    */
    void setListener(BaseListener* _listener)
    {
        listener = _listener;
    };

protected:
    std::string ID;// algorithm ID, an unique name to recognize the algorithm
    BaseListener* listener;   // Listener class

};

//-----------------------------------------------------------------------------
template<class DS>
std::ostream& operator << (std::ostream& _os, const ClassificationAlgorithm<DS>& _algo)
{
    _algo.printSelf(_os);
    return _os;
}
//@}
#endif//XMIPPBASEALGO_H
