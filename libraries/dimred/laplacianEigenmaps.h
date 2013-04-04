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

#ifndef AUTOENCODER_H_
#define AUTOENCODER_H_

#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include "dimred_tools.h"

/**@defgroup LaplacianEigenmap Laplacian Eigenmap
   @ingroup DimRedLibrary */
//@{
/** Class for making a Kernel PCA dimensionality reduction */
class LaplacianEigenmap: public DimRedAlgorithm
{
public:
	double sigma;
	size_t numberOfNeighbours;
public:
	/// Set specific parameters
	void setSpecificParameters(double sigma=1.0, size_t numberOfNeighbours=12);

	/// Reduce dimensionality
	void reduceDimensionality();
};
//@}
#endif
