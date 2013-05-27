/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
/* This file contains functions related to the SAMPLING Transform */

#ifndef _BASICPCA_HH
#define _BASICPCA_HH

#include <vector>
#include "multidim_array.h"


/**@defgroup BasicPCA Basic PCA class
   @ingroup DataLibrary */
//@{
/** Basic PCA class.
 *  The difference with PCAAnalyzer is that this one uses a different
 *  base class for the input vectors and the algorithm used to compute
 *  the PCA decomposition (see Roweis, "EM algorithms for PCA and SPCA",
 *  Neural Information Processing Systems 10 (NIPS'97) pp.626-632)
 *
 *  Example of use:
 *  @code
 *  PCAMahalanobisAnalyzer analyzer;
 *  FOR_ALL_VECTORS ...
 *     analyzer.addVector(v);
 *  analyzer.evaluateZScore(3,10); // 3 PCA components, 10 iterations
 *  @endcode
 *
 *  Another example:
 *  @code
 *  PCAMahalanobisAnalyzer analyzer;
 *  FOR_ALL_VECTORS ...
 *     analyzer.addVector(v);
 *  analyzer.subtractAvg();
 *  analyzer.learnPCABasis(NPCA, Niter);
    Matrix2D<double> proj;
    projectOnPCABasis(proj);
 *  @endcode
 *  */
class PCAMahalanobisAnalyzer
{
public:
    // Set of images assigned to the class
    std::vector< MultidimArray<float> > v;

    // Set of basis functions
    std::vector< MultidimArray<double> > PCAbasis;

    // The average of the vectors
    MultidimArray<double> avg;

    // Set of basis functions
    MultidimArray<double> Zscore;

    // Indexes to access in a sorted way
    MultidimArray<int> idx;

    // Matrix for the Eigen values
    Matrix1D<double> w;
public:
    /// Clear
    inline void clear()
    {
        v.clear();
        PCAbasis.clear();
        Zscore.clear();
        idx.clear();
        avg.clear();
    }

    /// Resize
    inline void reserve(int newSize)
    {
        v.reserve(newSize);
    }

    /// Add vector
    inline void addVector(const MultidimArray<float> &_v)
    {
        v.push_back(_v);
    }

    /// Compute Statistics
    void computeStatistics(MultidimArray<double> &avg,
                           MultidimArray<double> &stddev);

    /// Subtract average
    void subtractAvg();

    /**
     *This method computes the orthonormal basis for the vectors of an Eucldian
     *which is formed by the vectors on each column of the matrix. The output
     *of this method will be a matrix in which each column form an orthonormal
     *basis.
     */
    void gramSchmidt();

    /// Standardarize variables
    void standardarizeVariables();

    /// Learn basis
    void learnPCABasis(size_t NPCA, size_t Niter);

    /// Project on basis
    void projectOnPCABasis(Matrix2D<double> &CtY);

    /** Evaluate Zscore of the vectors stored with Mahalanobis distance.
     * NPCA is the dimension of the dimensionally reduced vectors before Mahalanobis
     * Niter is used to learn the PCA basis (typically, Niter=10).
     */
    void evaluateZScore(int NPCA, int Niter);

    /** Get the Zscore of vector n */
    inline double getZscore(int n)
    {
        return Zscore(n);
    }

    /** Get the Zscore of the n-th vector once sorted */
    inline double getSortedZscore(int n)
    {
        return Zscore(idx(n)-1);
    }

    /** Get the n-th index once sorted */
    inline int getSorted(int n)
    {
        return idx(n)-1;
    }
};
//@}
#endif
