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
double intrinsicDimensionality(Matrix2D<double> &X, const String &method="MLE", bool normalize=true);
//@}
#endif
