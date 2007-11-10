/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia, CSIC
 *
 * Copyright (c) 2001 , CNB-CSIC.
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'pascual@cnb.uam.es'
 *
 *****************************************************************************/

#ifndef XMIPPRM_H
#define XMIPPRM_H

using namespace std;

#include "base_algorithm.h"
#include "data_types.h"
#include "training_vector.h"

#include <data/funcs.h>
#include <data/matrix2d.h>

/**@defgroup RMClasess RM classes
   @ingroup ClassificationLibrary */
//@{
/** Basic RM class */
class xmippRM
{
public:

    /**
    * Make an empty xmippRP
    */
    xmippRM(void)
    {
        matdist = "gaussian";
    }

    /**
    * Sets the matrix distribution:
    * "gaussian" or "sparse"
    * Parameter: _type The type of the matrix distribution
    */
    void setRMDistribution(string _type)
    {
        matdist = _type;
    }


    /**
    * Calculate the Random matrix
    * Parameter: ts The vectors.
    * Parameter: k Number of components
    */
    void calculateRM(xmippCTVectors const &ts, int _k);

    /** Clear.
    Clean the Random Matrix*/
    void clear();


    /** Project onto RM space.
    Given a vector of the same size as the Random matrix this function
    returns the projection of size k onto the first k random vectors.

    An exception is thrown if the input vectors are not of the same size
    as the RM ones.*/

    void Project(xmippCTVectors &input, xmippCTVectors &output, int _k);

    /** Defines Listener class
      */
    void setListener(xmippBaseListener* _listener)
    {
        listener = _listener;
    };


private:
    string matdist; // Matrix distribution (sparse or gaussian)
    Matrix2D<double> RM;
    xmippBaseListener* listener;   // Listener class


};
//@}
#endif

