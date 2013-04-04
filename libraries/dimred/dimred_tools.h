/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2013)
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
#ifndef _DIMRED_TOOLS
#define _DIMRED_TOOLS

#include <data/matrix2d.h>
#include <data/matrix1d.h>

/**@defgroup DimRedTools Tools for dimensionality reduction
   @ingroup DimRedLibrary */
//@{
/** Generate test data.
 * Translated from drtools/generate_data.m
 * http://homepage.tudelft.nl/19j49/Matlab_Toolbox_for_Dimensionality_Reduction.html
 *
 * Original code by Laurens van der Maaten, Delft University of Technology
 * */
class GenerateData
{
public:
	/** Generated data.
	  * Each row of the matrix is an individual observation */
	Matrix2D<double> X;

	/** Underlying manifold coordinates.
	 * Eeach row of the matrix corresponds to the manifold coordinates of the observation in X.
	 */
	Matrix2D<double> t;

	/** Vector of labels for the observations */
	Matrix1D<unsigned char> label;

public:
	/** Generate data with a given number of points and a given method.
	 * Generates an artificial dataset. Possible datasets are: 'swiss' for the Swiss roll
	 * dataset, 'helix' for the helix dataset, 'twinpeaks' for the twinpeaks dataset,
     * '3d_clusters' for the 3D clusters dataset, and 'intersect' for the intersecting
     * dataset. The variable n indicates the number of datapoints to generate
     * (default = 1000). The variable noise indicates the amount of noise that
     * is added to the data (default = 0.05). The function generates the
     * high-dimensional dataset in X, and corresponding labels in labels. In
     * addition, the function keeps the coordinates of the datapoints on the
     * underlying manifold in t.
	 */
	void generateNewDataset(const String& method, int N=1000, double noise=0.05);
};

/** Function type to compute the squared distance between individuals i1 and i2 of X */
typedef double (*DimRedDistance2)  (const Matrix2D<double> &X, int i1, int i2);

/** Compute the distance of all vs all elements in a matrix of observations.
 * Each observation is a row of the matrix X.
 */
void computeDistance(const Matrix2D<double> &X, Matrix2D<double> &distance, DimRedDistance2* f=NULL, bool computeSqrt=true);

/** Compute the distance of each observation to its K nearest neighbours.
 * Each observation is a row of the matrix X.
 * If there are N observations, the size of distance is NxN.
 */
void computeDistanceToNeighbours(const Matrix2D<double> &X, int K, Matrix2D<double> &distance, DimRedDistance2* f=NULL, bool computeSqrt=true);

/** Compute a similarity matrix from a squared distance matrix.
 * dij=exp(-dij/(2*sigma^2))
 * The distance matrix can be previously normalized so that the maximum distance is 1
 */
void computeSimilarityMatrix(Matrix2D<double> &D2, double sigma, bool skipZeros=false, bool normalize=false);

/** Compute graph laplacian.
 * L=D-G where D is a diagonal matrix with the row sums of G.
 */
void computeGraphLaplacian(const Matrix2D<double> &G, Matrix2D<double> &L);

/** Estimate the intrinsic dimensionality.
 * Performs an estimation of the intrinsic dimensionality of dataset X based
 * on the method specified by method. Possible values for method are 'CorrDim'
 * (based on correlation dimension), 'NearNbDim' (based on nearest neighbor
 * dimension), 'GMST' (based on the analysis of the geodesic minimum spanning
 * tree), 'PackingNumbers' (based on the analysis of data packing numbers),
 * 'EigValue' (based on analysis of PCA eigenvalues), and 'MLE' (maximum
 * likelihood estimator). The default method is 'MLE'. All methods are
 * parameterless.
 *
 * Columns of the input matrix may be normalized to have zero mean and standard deviation 1.
 *
 * Translated from drtools/intrinsic_dimensionality.m
 * http://homepage.tudelft.nl/19j49/Matlab_Toolbox_for_Dimensionality_Reduction.html
 *
 * Original code by Laurens van der Maaten, Delft University of Technology
 */
double intrinsicDimensionality(Matrix2D<double> &X, const String &method="MLE", bool normalize=true, DimRedDistance2* f=NULL);

/** k-Nearest neighbours.
 * Given a data matrix (each row is a sample, each column a variable), this function
 * returns a matrix of the indexes of the K nearest neighbours to each one of the input samples sorted by distance.
 * It also returns the corresponding distance.
 *
 * The element i,j of the output matrices is the index(distance) of the j-th nearest neighbor to the i-th sample.
 *
 * You can provide a distance function of your own. If not, Euclidean distance is used.
 */
void kNearestNeighbours(const Matrix2D<double> &X, int K, Matrix2D<int> &idx, Matrix2D<double> &distance, DimRedDistance2* f=NULL, bool computeSqrt=true);

/** Extract k-nearest neighbours.
 * This function extracts from the matrix X, the neighbours given by idx for the i-th observation.
 */
void extractNearestNeighbours(const Matrix2D<double> &X, Matrix2D<int> &idx, int i, Matrix2D<double> &Xi);

/** Generic class for dimensionality reduction */
class DimRedAlgorithm
{
public:
	/// Pointer to input data
	Matrix2D<double> *X;

	/// Output dim
	size_t outputDim;

	/// Output data
	Matrix2D<double> Y;

	/// Distance function
	DimRedDistance2 *distance;
public:
	/// Empty constructor
	DimRedAlgorithm();

	/// Set input data
	void setInputData(Matrix2D<double> &X);

	/// Set output dimensionality
	void setOutputDimensionality(size_t outputDim);

	/// Reduce dimensionality
	virtual void reduceDimensionality()=0;

	/// Get reduced data
	const Matrix2D<double> &getReducedData();
};
//@}
#endif
