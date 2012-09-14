/***************************************************************************
 *
 * Authors:     Vahid Abrishami (vabrishami@cnb.csic.es)
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

#ifndef KNN_CLASSIFIER_HH
#define KNN_CLASSIFIER_HH

/* Includes-----------------------------------------------------------------*/
#include <data/multidim_array.h>

/**@defgroup KNN Classifier
   @ingroup ClassificationLibrary */
//@{
/**
 * This class implements KNN (K Nearest Neighbors). It does the classification
 * by a voting of nearest neighbors of the point. It uses the exhaustive search
 * in order to find the nearest neighbors.
 */

class KNN
{
public:
    /// Type of distance
    typedef enum { EUCLIDEAN = 0, CITYBLOCK = 1} distType;

    /**
     * Constructs the algorithm
     * Parameter: k           The number of nearest neighbors
     */
    KNN(int k);

    /**
     * Now actually there is no training step.In future more
     * complex train step can be added.
     */
    void train(MultidimArray<double> &dataset,MultidimArray<double> &dataLabel,
               MultidimArray<double> &labelset);

    /**
     * Get an item an predict the related class of that by
     * doing the voting among its neighbors.
     */
    int predict(MultidimArray<double> &sample,double &score);

    /// Compute the K nearest neighbors to the sample
    void KNearestNeighbors(MultidimArray<double> &sample);

    /// Save the model for the classifier
    void saveModel(const FileName &fn);

    /// Load the model for the classifier
    void loadModel(const FileName &fn);

    /// Method for setting the K
    void setK(int k);

private:
    /// Compute the euclidean distance
    double euclideanDistance(MultidimArray<double> &sample,int index,double maximumDist);

    /// Compute the City Block (Manhattan) distance
    double cityBlockDistance(MultidimArray<double> &sample,int index,double maximumDist);

    /// This function find the index of maximum value in an 1D-array
    int findMaxIndex(MultidimArray<double> &inputArray);

    /// This function find the index of minimum value in an 1D-array
    int findMinIndex(MultidimArray<double> &inputArray);

private:
    /// The number of nearest neighbors
    int K;

    /// Pointer to the current dataset
    MultidimArray<double>  __dataset;

    /// Pointer to the current labelset
    MultidimArray<double>  __dataLabel;

    /// Pointer to the current labelset
    MultidimArray<double> __labelSet;

    /// Neighbors index in dataset
    MultidimArray<int> neighborsIndex;

    /// Neighbors distance in dataset
    MultidimArray<double> maxDist;
};
//@}
#endif
