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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
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
 */
double orthogonalMatchingPursuit(const Matrix1D<double> &x,
    const Matrix2D<double> &D, int S, Matrix1D<double> &alpha);

/** kSVD
 * @ingroup kSVD
 *  This function optimizes the dictionary input and the representation
 *  of a set of vectors in this dictionary using at most S atoms.
 */
double kSVD(const std::vector< Matrix1D<double> > &X, int S,
    Matrix2D<double> &D, std::vector< Matrix1D<double> > &Alpha,
    bool keepFirstColumn=true);
//@}
#endif
