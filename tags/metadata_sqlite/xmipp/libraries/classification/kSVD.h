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

#ifndef XMIPP__KSVD_HH__
#define XMIPP__KSVD_HH__

/* Includes ---------------------------------------------------------------- */
#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include <vector>

/**@defgroup kSVD
   @ingroup ClassificationLibrary */
//@{

/** Orthogonal matching pursuit (OMP)
 * @ingroup kSVD
 *  This function looks for the best approximation with only S components
 *  of a given dictionary D. The representation is alpha, i.e., x=D*alpha.
 *
 * The approximation error is returned.
 *
 * This implementation is based on Karl Skretting
 * http://www.ux.uis.no/~karlsk/proj02/index.html
 */
double orthogonalMatchingPursuit(const Matrix1D<double> &x,
    const Matrix2D<double> &D, int S, Matrix1D<double> &alpha);

/** Lasso projection
 * @ingroup kSVD
 *  This function looks for the best approximation with as few components
 *  as possible (regularized through a L1 norm)
 *  of a given dictionary D. The representation is alpha, i.e., x=D*alpha.
 * 
 * DtD is D^t*D and must be precomputed before calling the dictionary
 *
 * The approximation error is returned.
 *
 * This implementation is based on Mark Schmidt
 * http://www.cs.ubc.ca/~schmidtm/Software/lasso.html
 */
double lasso(const Matrix1D<double> &x,
    const Matrix2D<double> &D, 
    const Matrix2D<double> &DtD,
    const Matrix2D<double> &DtDlambdaInv,
    double lambda, Matrix1D<double> &alpha,
    const int maxIter=20, const double tol=0.005);

#define OMP_PROJECTION 1
#define LASSO_PROJECTION 2
/** kSVD
 * @ingroup kSVD
 *  This function optimizes the dictionary input and the representation
 *  of a set of vectors in this dictionary using at most S atoms.
 *
 * Valid projection methods are OMP_PROJECTION and LASSO_PROJECTION
 *
 * This implementation is based on Michael Elad
 * http://www.cs.technion.ac.il/~elad/software
 */
double kSVD(const std::vector< Matrix1D<double> > &X, int S,
    Matrix2D<double> &D, std::vector< Matrix1D<double> > &Alpha,
    bool keepFirstColumn=true, int maxIter=20, double minChange=0.005,
    int projectionMethod=OMP_PROJECTION, double lambda=0.2);
//@}
#endif
