/***************************************************************************
 *
 * Authors:    Emma Sesmero            emmasesmero@gmail.com (2013)
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
#ifndef _HESSIAN_LLE
#define _HESSIAN_LLE

#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include "dimred_tools.h"

/**@defgroup HessianLLE Hessian Local Linear Embedding
   @ingroup DimRedLibrary */
//@{
/** Class for making a Hessian LLE dimensionality reduction */
class HessianLLE: public DimRedAlgorithm
{
public:
    int kNeighbours;
public:
    /// Set specific parameters
    void setSpecificParameters(int kNeighbours=12);

    /// Reduce dimensionality
    void reduceDimensionality();

private:
    void modifiedGramSchmidtOrtogonalization(Matrix2D<double> *matrix);
    void completeYt(const Matrix2D<double> &V, const Matrix2D<double> &Yi, Matrix2D<double> &Yt_complete);
    void buildYiHessianEstimator(const Matrix2D<double> &V, Matrix2D<double> &Yi, size_t no_dim, size_t dp);
};
//@}
#endif
