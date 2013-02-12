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
// FuzzyCMeans.cc
// Fuzzy c-means clustering algorithm
//-----------------------------------------------------------------------------

#include "fcmeans.h"

/**  Ctor from stream
 * Parameter: _is Must have the parameters in the same order than the previous ctor.
   ****** check out this ************
 */
/*FuzzyCMeans::FuzzyCMeans( std::istream& _is )
  :ClassificationAlgorithm< FuzzyCodeBook >( "FuzzyCMeans") {
  _is >> m;
  _is >> epsilon;
  _is >> epochs;
};
*/

/**
 * Trains the Algorithm
 * Parameter: _xmippDS Data structure to train, a codeBook in this case
 * Parameter: _examples  A training set with the training examples
 */

void FuzzyCMeans::train(FuzzyCodeBook& _xmippDS, TS& _examples) const
{

    // Defines verbosity

    int verbosity = listener->getVerbosity();
    if (verbosity)
        listener->OnReportOperation((std::string) "Training....\n");
    if (verbosity == 1 || verbosity == 3)
        listener->OnInitOperation(epochs);

    // Create auxiliar Codebook

    FuzzyCodeBook auxCB;


    // Create auxiliar stuff

    unsigned numClusters = _xmippDS.size();
    unsigned numVectors = _examples.size();
    unsigned i, j, k;
    double stopError = 0, auxError = 0;
    double auxDist, auxProd, tmp, auxExp, auxSum;
    unsigned t = 0;  // Iteration index
    FeatureVector zero(_xmippDS.theItems[0].size()) ;
    fill(zero.begin(), zero.end(), 0.0);


    // Initialize auxiliary Codebook

    auxCB = _xmippDS;

    // Set auxExp

    auxExp = 2 / (m - 1);

    // This is the main code of the algorithm. Iterates "epochs" times


    stopError = epsilon + 1; // Initially must be higher

    while ((stopError > epsilon) && (t < epochs))
    {

        // Update Membership matrix

        for (k = 0; k < numVectors; k++)
        {

            auxProd = 1;
            for (j = 0; j < numClusters; j++)
                auxProd *= euclideanDistance(_xmippDS.theItems[j], _examples.theItems[k]);

            if (auxProd == 0.)
            { // Apply k-means criterion (Data-CB) must be > 0
                for (j = 0; j < numClusters; j ++)
                    if (euclideanDistance(_xmippDS.theItems[j], _examples.theItems[k]) == 0.)
                        _xmippDS.memb[k][j] = 1.0;
                    else
                        _xmippDS.memb[k][j] =  0.0;
            }
            else
            {
                for (i = 0; i < numClusters; i ++)
                {
                    auxDist = 0;
                    for (j = 0; j < numClusters; j ++)
                    {
                        tmp = euclideanDistance(_xmippDS.theItems[i], _examples.theItems[k]) /
                              euclideanDistance(_xmippDS.theItems[j], _examples.theItems[k]);
                        auxDist += pow(tmp, auxExp);
                    } // for j
                    _xmippDS.memb[k][i] = (floatFeature) 1.0 / auxDist;
                } // for i
            } // if auxProd
        } // for k


        // Update code vectors (Cluster Centers)

        for (i = 0; i < numClusters; i++)
        {
            _xmippDS.theItems[i] = zero;
            auxSum = 0;
            for (k = 0; k < numVectors; k++)
            {
                _xmippDS.theItems[i] += (floatFeature) pow((double)(_xmippDS.memb[k][i]), m) * _examples.theItems[k];
                auxSum += pow((double)(_xmippDS.memb[k][i]), m);
            } // for i
            _xmippDS.theItems[i] /= (floatFeature) auxSum;
        } // for k

        // Compute stopping criterion
        stopError = 0;
        for (i = 0; i < numClusters; i ++)
        {
            auxError = euclideanDistance(_xmippDS.theItems[i], auxCB.theItems[i]);
            stopError += auxError * auxError;
        } // for i


        // Update iteration index

        t++;
        auxCB = _xmippDS;

        if (verbosity == 1 || verbosity == 3)
            listener->OnProgress(t);
        if (verbosity >= 2)
        {
            char s[100];
            sprintf(s, "Iteration %d of %d. Code vectors variation: %g\n", t + 1, epochs, stopError);
            listener->OnReportOperation((std::string) s);
        }

    } // while

    if (verbosity == 1 || verbosity == 3)
        listener->OnProgress(epochs);

}// FuzzyCMeans::train

/**
 * Test the Algorithm in a conventional way
 * Parameter: _examples  A training set with the training examples
 */

double FuzzyCMeans::test(const FuzzyCodeBook& _xmippDS,
                         const TS& _examples) const
{

    // Defines verbosity

    int verbosity = listener->getVerbosity();
    if (verbosity)
    {
        listener->OnReportOperation((std::string) "Testing....\n");
        listener->OnInitOperation(_examples.size());
    }

    double distortion = 0;
    for (unsigned i = 0; i < _examples.size(); i ++)
    {
        const FeatureVector& auxS = _examples.theItems[i];
        unsigned best = _xmippDS.output(auxS);
        distortion += euclideanDistance(_xmippDS.theItems[best], _examples.theItems[i]);
        if (verbosity)
            listener->OnProgress(i);
    };

    if (verbosity)
        listener->OnProgress(_examples.size());

    return distortion / (double) _examples.size();
}

/**
 * Tests with the training set using for training.
 * Parameter: _examples  The training set
 */

double FuzzyCMeans::fuzzyTest(const FuzzyCodeBook& _xmippDS,
                              const TS& _examples) const
{

    // Defines verbosity

    int verbosity = listener->getVerbosity();
    if (verbosity)
    {
        listener->OnReportOperation((std::string) "Testing....\n");
        listener->OnInitOperation(_examples.size());
    }

    double distortion = 0;
    for (unsigned i = 0; i < _examples.size(); i ++)
    {
        unsigned best = _xmippDS.fuzzyOutput(i);
        distortion += euclideanDistance(_xmippDS.theItems[best], _examples.theItems[i]);
        if (verbosity)
            listener->OnProgress(i);
    };
    if (verbosity)
        listener->OnProgress(_examples.size());

    return distortion / (double)_examples.size();
}

/**
 * Calculates Partition Coefficient (F) validity functional
 * Parameter: _xmippDS Data structure to train, a codeBook in this case
 *
 * Notes on F:  For U in Mfc (fuzzy partition space)
 *              1/C <= F <= 1
 *              for F = 1, U is hard (zeros and ones only)
 *              for F = 1/C, U = 1/C*ones(C,n);
 *
 * (max)
 *
 * For more information see:
 *     J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
 *     Function Algorithms", Plenum Press, New York, 1981.
 *
 */

double FuzzyCMeans::F(const FuzzyCodeBook& _xmippDS) const
{
    double F = 0.;
    for (unsigned k = 0; k < _xmippDS.membVectors(); k++)
        for (unsigned i = 0; i < _xmippDS.membClusters(); i++)
            F += pow((double)(_xmippDS.memb[k][i]), 2);
    return (F / (double)(_xmippDS.membVectors()));
}

/**
 * Calculates Partition Entropy (H) validity functional
 * Parameter: _xmippDS Data structure to train, a codeBook in this case
 *
 * Notes on H:  For U in Mfc
 *              0 <= H <= log(C)
 *              for H = 0, U is hard
 *              for H = log(C), U = 1/C*ones(C,n);
 *              0 <= 1 - F <= H (strict inequality if U not hard)
 *
 * (min)
 *
 * For more information see:
 *     J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
 *     Function Algorithms", Plenum Press, New York, 1981.
 *
 */

double FuzzyCMeans::H(const FuzzyCodeBook& _xmippDS) const
{
    double H = 0.;
    for (unsigned k = 0; k < _xmippDS.membVectors(); k++)
        for (unsigned i = 0; i < _xmippDS.membClusters(); i++)
            if (_xmippDS.memb[k][i] != 0.)
                H += (double)(_xmippDS.memb[k][i]) * log((double)(_xmippDS.memb[k][i]));
    return (-H / (double)(_xmippDS.membVectors()));
}


/**
 * Calculates Non-fuzzy index (NFI) validity functional
 * Parameter: _xmippDS Data structure to train, a codeBook in this case
 *
 * (max)
 *
 * For more information see:
 *     M. Roubens, "Pattern Classification Problems and Fuzzy Sets",
 *     Fuzzy Sets and Systems, 1:239-253, 1978.
 *
 */

double FuzzyCMeans::NFI(const FuzzyCodeBook& _xmippDS) const
{
    double F = 0.;
    for (unsigned k = 0; k < _xmippDS.membVectors(); k++)
        for (unsigned i = 0; i < _xmippDS.membClusters(); i++)
            F += pow((double)(_xmippDS.memb[k][i]), 2);

    double NFI = (((double)(_xmippDS.membClusters()) * F - (double)(_xmippDS.membVectors())) / (double)(_xmippDS.membVectors())) / (double)(_xmippDS.membClusters() - 1);
    return NFI;
}


/**
 * Calculates Compactness and separation index (S) validity functional
 * Parameter: _xmippDS Data structure to train, a codeBook in this case
 * Parameter: _examples  A training set with the training examples
 *
 * (min)
 *
 * For more information see:
 *     X.L. Xie and G. Beni, "A Validity Measure for Fuzzy Clustering",
 *     IEEE Trans. PAMI, 13(8):841-847, 1991.
 *
 */

double FuzzyCMeans::S(const FuzzyCodeBook& _xmippDS,
                      const TS& _examples) const
{

    std::vector< std::vector< floatFeature > > ICD;       // Intercluster distance
    std::vector< std::vector< floatFeature > > D;         // Distance from each data to cluster centers

    unsigned i;
    D.resize(_xmippDS.membClusters());
    for (i = 0; i < _xmippDS.membClusters(); i++)
    {
        std::vector <floatFeature> d;
        d.resize(_xmippDS.membVectors());
        for (unsigned k = 0; k < _xmippDS.membVectors(); k++)
            d[k] = (floatFeature)euclideanDistance(_xmippDS.theItems[i], _examples.theItems[k]);
        D[i] = d;
    } // for i

    ICD.resize(_xmippDS.membClusters());
    for (i = 0; i < _xmippDS.membClusters(); i++)
    {
        std::vector <floatFeature> v;
        v.resize(_xmippDS.membVectors());
        for (unsigned j = 0; j < _xmippDS.membClusters(); j++)
            v[j] = (floatFeature)euclideanDistance(_xmippDS.theItems[i], _xmippDS.theItems[j]);
        ICD[i] = v;
    } // for i

    floatFeature auxSum = 0;
    for (i = 0; i < _xmippDS.membClusters(); i++)
        for (unsigned k = 0; k < _xmippDS.membVectors(); k++)
            auxSum += (floatFeature) pow((double)(D[i][k] * _xmippDS.memb[k][i]), (double)m);

    floatFeature auxMin = MAXFLOAT;
    for (i = 0; i < _xmippDS.membClusters(); i++)
        for (unsigned j = i + 1; j < _xmippDS.membClusters(); j++)
            if (auxMin > ICD[i][j])
                auxMin = ICD[i][j];

    double S = auxSum / (double)(_xmippDS.membVectors()) / (double)(auxMin);
    return S;

}

/// print itself on standard output
void FuzzyCMeans::printSelf(std::ostream& _os) const
{
    // Call base class, which will print ID
    _os << "Class (Algorithm): " << std::endl;
    ClassificationAlgorithm<FuzzyCodeBook>::printSelf(_os);
    _os << std::endl;
    // Print parameters in the same order they are declared
    _os << "Fuzzy constant m = " << m << std::endl;
    _os << "Epsilon eps = " << epsilon << std::endl;
    _os << "Iterations iter = " << epochs << std::endl << std::endl;
}
