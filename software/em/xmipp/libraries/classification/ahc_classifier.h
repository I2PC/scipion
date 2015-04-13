/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano   (coss@cnb.csic.es)
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

#ifndef XMIPP__AHC_CLASSIFIER_HH__
#define XMIPP__AHC_CLASSIFIER_HH__

/* Includes ---------------------------------------------------------------- */
#include <data/matrix2d.h>

/**@defgroup AHCClassifier Agglomerative Hierarchical Clustering
   @ingroup ClassificationLibrary */
//@{

/** AHC classifier class.
 */
class AHCClassifier
{
public:
	/** The i-th location is the number of the cluster of the i-th element in the input data.
	 * The length of this vector is the number of individuals.
	 */
	Matrix1D<int> clusterAssigned;

	/** There are as many elements in the vector as clusters.
	 * Each vector in the cluster contains the list of individuals assigned to this cluster.
	 */
	std::vector< std::vector<int> > cluster;
public:
	/** Cluster data.
	 * X is the data to classify, each row is an observation. Columns are features of that observation.
	 *
	 * For the distances see http://www.alglib.net/translator/man/manual.cpp.html#sub_clusterizersetpoints
	 *  0    Chebyshev distance  (L-inf norm)
	 *  1    city block distance (L1 norm)
	 *  2    Euclidean distance  (L2 norm)
	 * 10    Pearson correlation:
	 *       dist(a,b) = 1-corr(a,b)
	 * 11    Absolute Pearson correlation:
 	 *       dist(a,b) = 1-|corr(a,b)|
 	 * 12    Uncentered Pearson correlation (cosine of the angle):
	 *       dist(a,b) = a'*b/(|a|*|b|)
	 * 13    Absolute uncentered Pearson correlation
	 *       dist(a,b) = |a'*b|/(|a|*|b|)
	 * 20    Spearman rank correlation:
	 *       dist(a,b) = 1-rankcorr(a,b)
	 * 21    Absolute Spearman rank correlation
	 *	     dist(a,b) = 1-|rankcorr(a,b)|
	 *
	 * For linkage:
	 * 0     complete linkage (default algorithm)
	 * 1     single linkage
 	 * 2     unweighted average linkage
 	 * 3     weighted average linkage
	 */
	void clusterData(const Matrix2D<double> &X, int numberOfClusters=2, int distance=2, int linkageType=1);

	/** Cluster given a distance matrix.
	 * X is the data to classify, each row is an observation. Columns are features of that observation.
	 *
	 * For linkage:
	 * 0     complete linkage (default algorithm)
	 * 1     single linkage
 	 * 2     unweighted average linkage
 	 * 3     weighted average linkage
	 */
	void clusterWithDistance(const Matrix2D<double> &D, int numberOfClusters=2, int linkageType=1);
};
//@}
#endif
