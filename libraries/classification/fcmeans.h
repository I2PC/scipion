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
// FuzzyCMeans.hh
// Fuzzy c-means clustering algorithm
//-----------------------------------------------------------------------------

#ifndef _XMIPPFCMEANS_H
#define _XMIPPFCMEANS_H

#include <iostream>
#include <vector>

#include "base_algorithm.h"
#include "fuzzy_code_book.h"

/**@defgroup FuzzyCMeans Fuzzy c-means clustering algorithm
   @ingroup ClassificationLibrary */
//@{
/**
 *  This class implements Fuzzy c-means clustering method (Bezdeck)
 *  an unsupervised clustering algorithm.
*/
class FuzzyCMeans:  public ClassificationAlgorithm<FuzzyCodeBook >
{

public:

    /**
     * Big mega ctor. Creates a Fuzzy c-means codebook, and initializes it
     * Parameter: _m   Fuzzy constant
     * Parameter: _epsilon  Stopping criterion
     * Parameter: _epochs Number of epochs or iterations
    */
    FuzzyCMeans(double _m, double _epsilon, unsigned _epochs)
            : ClassificationAlgorithm< FuzzyCodeBook >("xmippFCMeans"),
            m(_m), epsilon(_epsilon),
            epochs(_epochs)
    {}
    ;


    /*
     * Ctor from stream
     * Parameter: _is Must have the parameters in the same order than the previous ctor.
     */
    //  FuzzyCMeans( std::istream& _is );


    /*
     * Virtual destructor
     */
    virtual ~FuzzyCMeans()
    {}
    ;


    /**
     * Trains the Algorithm
     * Parameter: _xmippDS Data structure to train, a codeBook in this case
     * Parameter: _examples  A training set with the training examples
     */
    virtual void train(FuzzyCodeBook& _xmippDS,
                       TS& _examples) const;


    /**
     * Tests with the training set using for training.
     * Fuzzy membership is used for testing
     * Parameter: _xmippDS Data structure to train, a codeBook in this case
     * Parameter: _examples  The training set
     * returns the quantization error
     */
    virtual double fuzzyTest(const FuzzyCodeBook& _xmippDS,
                             const TS& _examples) const;


    /**
     * Tests the Algorithm in a conventional way.
     * Parameter: _xmippDS Data structure to train, a codeBook in this case
     * Parameter: _examples  A training set with the training examples
     */
    virtual double test(const FuzzyCodeBook& _xmippDS,
                        const TS& _examples) const;


    /**
     * Calculates Partition Coefficient (F) validity functional
     * Parameter: _xmippDS Data structure to train, a codeBook in this case
     * (It should be maximum)
     * For more information see:
     *     J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
     *     Function Algorithms", Plenum Press, New York, 1981.
     *
     */
    double F(const FuzzyCodeBook& _xmippDS) const;

    /**
     * Calculates Partition Entropy (H) validity functional
     * Parameter: _xmippDS Data structure to train, a codeBook in this case
     * (It should be minimum)
     * For more information see:
     *     J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
     *     Function Algorithms", Plenum Press, New York, 1981.
     *
     */

    double H(const FuzzyCodeBook& _xmippDS) const;


    /**
     * Calculates Non-fuzzy index (NFI) validity functional
     * Parameter: _xmippDS Data structure to train, a codeBook in this case
     * (It should be maximum)
     * For more information see:
     *     M. Roubens, "Pattern Classification Problems and Fuzzy Sets",
     *     Fuzzy Sets and Systems, 1:239-253, 1978.
     *
     */

    double NFI(const FuzzyCodeBook& _xmippDS) const;

    /**
     * Calculates Compactness and separation index (S) validity functional
     * Parameter: _xmippDS Data structure to train, a codeBook in this case
     * Parameter: _examples  A training set with the training examples
     * (It should be minimum)
     * For more information see:
     *     X.L. Xie and G. Beni, "A Validity Measure for Fuzzy Clustering",
     *     IEEE Trans. PAMI, 13(8):841-847, 1991.
     *
     */
    double S(const FuzzyCodeBook& _xmippDS, const TS& _examples) const;

protected:

    /// print itself on standard output
    void printSelf(std::ostream& _os) const ;

    double m;       // Fuzzy constant
    double epsilon;      // Stopping criterion Error < epsilon
    unsigned epochs;                  // Number of presentations of the whole sample
};
//@}
#endif//_XMIPPFCMEANS_H
