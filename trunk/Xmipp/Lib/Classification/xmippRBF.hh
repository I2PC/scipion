/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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

#ifndef _XMIPP_RBF
   #define _XMIPP_RBF

#include "xmippCTVectors.hh"
#include <XmippData/xmippMatrices2D.hh>

class xmippRBF;

/**@name Radial Basis Functions. */
//@{
/** Compute design matrix.
    Given the set of centers (C), the set of radius and the set of points (X)
    this function computes the design matrix of the experiment. Element (i,j)
    is the value of the gaussian centered at center j on point i. It throws
    an exception if the centers and points are not of the same size.
    
    Look the Mark Orr's function rbf_dm */
void RBF_design_matrix(xmippCTVectors &C, matrix1D<double> &r,
   xmippCTVectors &X, matrix2D<double> &H) _THROW;

    
/** Compute the best scale.
    Range over a set of scales looking for the best one.
    The results are the Centers are radius that better fit the function.
    The input vectors are modified.
    
    Exceptions are thrown if any of the input dimensions do not
    meet what expected.*/
void RBF_train_best_scale(xmippCTVectors &candidate_C,  xmippCTVectors &X,
   vector<double> &y, double minscale,
   double maxscale, double scalestep, xmippRBF &RBF,
   double &error, vector<double> &y_predicted) _THROW;

/** Compute model given a certain scale.
    The selected centers are chosen via the index_out variable.
    */
void RBF_train(xmippCTVectors &C,  xmippCTVectors &X,
   vector<double> &y, matrix1D<double> &r, double scale,
   vector<int> &idx_out, matrix1D<double> &r_out,
   matrix1D<double> &w_out, double &error);

/** Predict the values at given points.
    Given the centers, weight and radii this function returns the predicted
    values at given X. */
void RBF_predict(xmippRBF &RBF,  xmippCTVectors &X, vector<double> &y_predicted);

/** RBF class.
    This class holds a RBF. */
class xmippRBF {
public:
   /// Radii
   matrix1D<double> r;
   /// Vector of weights
   matrix1D<double> w;
   /// Vector of Centers
   xmippCTVectors C;
public:
   /** Clear model. */
   void clear() {r.clear(); w.clear(); C.clear();}
   
   /** Model size. Number of centers*/
   int model_size() {return XSIZE(w);}
   
   /** Train model.
       Given the candidate centers, the set of pairs (X,y) this function
       returns the error of prediction (defined as the standard deviation
       of the predicted value with respect to y) and vector of predicted 
       values for the training set. */
   void train(xmippCTVectors &candidate_C, xmippCTVectors &X, vector<double> &y, 
      double &error, vector<double> &y_predicted, double minscale=0.1,
      double maxscale=1.0, double scalestep=0.15) {
      RBF_train_best_scale(candidate_C, X, y, minscale, maxscale, scalestep,
         *this, error, y_predicted);}

   /** Predict the values at given points.
       Given the points X this function returns the predicted points.*/
   void predict(xmippCTVectors &X, vector<double> &y_predicted)
      {RBF_predict(*this,X,y_predicted);}

   /** Show an RBF. */
   friend ostream & operator << (ostream &out, xmippRBF &rbf);

   /** Read an RBF. */
   friend istream & operator >> (istream &in, xmippRBF &rbf);
};

//@}
#endif
