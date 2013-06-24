/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2013)
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

#include "spe.h"
#include <data/numerical_tools.h>

void SPE::setSpecificParameters(bool global, int k)
{
	this->global = global;
	this->k = k;
}

void SPE::reduceDimensionality()
{
	int lambda = 1;  // lambda
	int s = 100; // number of updates per iteration
	int max_iter = 20000 + round(0.04 * MAT_YSIZE(*X) * MAT_YSIZE(*X)); // number of iterations
	double tol = 1e-5;
	int n = MAT_YSIZE(*X);
	if (global)
		max_iter *= 3;

	MultidimArray<int> J;
	int ind1_0=0, ind1_F=s;
	int ind2_0=s+1, ind2_F=std::min(2*s,n-1);
	for(int i=1;i<max_iter;i++){
		randomPermutation(n,J);

		// Update lambda
		lambda = lambda - (lambda / max_iter);
	}
}
