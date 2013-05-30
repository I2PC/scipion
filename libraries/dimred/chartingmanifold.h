/***************************************************************************
 *
 * Authors:    Clara M. Mejias Perez          clrmejias7@gmail.com (2013)
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
#ifndef _CHARTING_MANIFOLD
#define _CHARTING_MANIFOLD

#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include <data/multidim_array.h>
#include "dimred_tools.h"

/**@defgroup ChartingManifold
   @ingroup DimRedLibrary */
//@{
/** Class for making a Charting Manifold dimensionality reduction */
class ChartingManifold: public DimRedAlgorithm
{
public:
	int max_iterations;
	int no_analyzers;
public:
	/// Set specific parameters
	void setSpecificParameters(int max_iterations=200, int no_analyzers=40);

	/// Reduce dimensionality
	void reduceDimensionality();
};
//@}
#endif
