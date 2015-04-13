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
// FuzzyKohonenCMeans.cc
// Fuzzy Kohonen Clustering Network Algorithm
//-----------------------------------------------------------------------------

#include "fkcn.h"

#ifdef __sun
#include <ieeefp.h>
#endif

/**
 * Trains the algorithm
 * Parameter: _xmippDS Data structure to train, a codeBook in this case
 * Parameter: _examples  A training set with the training examples
 */

void FuzzyKohonenCMeans::train(FuzzyCodeBook& _xmippDS, const TS& _examples) const
{
    using namespace std;
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
    unsigned i, j, k, cc, vv;
    double stopError = 0, auxError = 0;
    double auxDist, auxProd, auxExp, auxSum;
    unsigned t = 0;  // Iteration index

    // Initialize auxiliary Codebook
    auxCB = _xmippDS;
    std::vector<FeatureVector> alpha;
    alpha.resize(numVectors);
    for (vv = 0; vv < numVectors; vv++)
        alpha[vv].resize(numClusters, 0.);

    // Set auxiliar variables
    double Deltam = (double)(m - 1.) / (double) epochs;

    // This is the main code of the algorithm. Iterates "epochs" times
    stopError = epsilon + 1; // Initially must be higher

    while ((stopError > epsilon) && (t < epochs))
    {

        double mt = (m - t * Deltam);
        if ((mt - 1.) <= 1e-2)
            auxExp = 2. / 1e-2;
        else
            auxExp = 2. / (mt - 1.);

        // Update Membership matrix

        FeatureVector tmpD;
        tmpD.resize(numClusters);
        for (k = 0; k < numVectors; k++)
        {
            auxProd = 0;
            for (i = 0; i < numClusters; i ++)
            {
                auxDist = 0;
                for (size_t d = 0; d < _examples.theItems[0].size(); d++)
                    auxDist += ((double)(_examples.theItems[k][d]) - (double)(_xmippDS.theItems[i][d])) * ((double)(_examples.theItems[k][d]) - (double)(_xmippDS.theItems[i][d]));
                auxDist = (double) sqrt((double)auxDist);
                auxDist = (double) pow((double) auxDist, (double) auxExp);
                if (auxDist < MAXZERO) auxDist = MAXZERO;
                if (isnan(auxDist)) auxDist = MAXZERO;
                if (!finite(auxDist)) auxDist = 1e200;
                auxProd += 1. / auxDist;
                tmpD[i] = auxDist;
            }
            for (j = 0; j < numClusters; j ++)
            {
                _xmippDS.memb[k][j] = (floatFeature)(1. / (auxProd * tmpD[j]));
            }
        } // for k



        // Calculate Alpha
        for (cc = 0; cc < numClusters; cc++)
            for (vv = 0; vv < numVectors; vv++)
                alpha[vv][cc] = (floatFeature) pow((double)_xmippDS.memb[vv][cc], (double)mt);


        /* Step III: Update Code Vectors */


        for (cc = 0; cc < numClusters; cc++)
        {
            std::vector<double> tmpV;
            tmpV.resize(_examples.theItems[0].size(), 0.);
            auxSum = 0;
            for (vv = 0; vv < numVectors; vv++)
            {
                for (j = 0; j < _examples.theItems[0].size(); j++)
                    tmpV[j] += (double)((double)(alpha[vv][cc]) * ((double)(_examples.theItems[vv][j]) - (double)(_xmippDS.theItems[cc][j])));
                auxSum += (double)(alpha[vv][cc]);
            } // for vv
            if (auxSum != 0.)
                for (j = 0; j < _examples.theItems[0].size(); j++)
                    _xmippDS.theItems[cc][j] += (floatFeature)(tmpV[j] / auxSum);
        } // for cc


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
            sprintf(s, "Iteration %d of %d. Code vectors variation: %g\n", t, epochs, stopError);
            listener->OnReportOperation((std::string) s);
        }

    } // while

    if (verbosity == 1 || verbosity == 3)
        listener->OnProgress(epochs);

} // FuzzyKohonenCMeans::train

