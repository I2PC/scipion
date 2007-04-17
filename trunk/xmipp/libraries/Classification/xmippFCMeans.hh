/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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
//-----------------------------------------------------------------------------
// xmippFCMeans.hh
// Fuzzy c-means clustering algorithm
//-----------------------------------------------------------------------------

#ifndef _XMIPPFCMEANS_H
#define _XMIPPFCMEANS_H

#pragma warning(disable:4786)

#include <iostream>		// for cout
#include <vector>               // vectors
#include <Classification/xmippDistances.hh>
#include <Classification/xmippBaseAlgo.hh>
#include <Classification/xmippFuzzyCodeBook.hh>


/**@name Fuzzy c-means clustering algorithm*/
//@{
/** 
 *  This class implements Fuzzy c-means clustering method (Bezdeck) 
 *  an unsupervised clustering algorithm. 
*/
class xmippFCMeans:  public xmippBaseAlgo<xmippFCB > {
  
public:

  /**  
   * Big mega ctor. Creates a Fuzzy c-means codebook, and initializes it
   * @param _m   Fuzzy constant
   * @param _epsilon  Stopping criterion
   * @param _epochs Number of epochs or iterations
  */  
  xmippFCMeans( double _m, double _epsilon, unsigned _epochs)
    :xmippBaseAlgo< xmippFCB >( "xmippFCMeans"), 
     m(_m), epsilon (_epsilon),
     epochs (_epochs ) {};


  /*  
   * Ctor from stream
   * @param _is Must have the parameters in the same order than the previous ctor.
   */
//  xmippFCMeans( istream& _is );


  /*  
   * Virtual destructor
   */
  virtual ~xmippFCMeans() {};
    

  /**
   * Trains the Algorithm
   * @param _xmippDS Data structure to train, a codeBook in this case
   * @param _examples  A training set with the training examples
   */   
  virtual void train(xmippFCB& _xmippDS, 
		     TS& _examples) const;


  /**
   * Tests with the training set using for training.
   * Fuzzy membership is used for testing
   * @param _xmippDS Data structure to train, a codeBook in this case
   * @param _examples  The training set
   * returns the quantization error
   */
  virtual double fuzzyTest(const xmippFCB& _xmippDS, 
		      const TS& _examples) const;


  /**
   * Tests the Algorithm in a conventional way.
   * @param _xmippDS Data structure to train, a codeBook in this case
   * @param _examples  A training set with the training examples
   */
  virtual double test(const xmippFCB& _xmippDS, 
		      const TS& _examples) const;


  /**
   * Calculates Partition Coefficient (F) validity functional
   * @param _xmippDS Data structure to train, a codeBook in this case
   * (It should be maximum)
   * For more information see:
   *     J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
   *     Function Algorithms", Plenum Press, New York, 1981.
   * 
   */
  double F(const xmippFCB& _xmippDS) const;

  /**
   * Calculates Partition Entropy (H) validity functional
   * @param _xmippDS Data structure to train, a codeBook in this case
   * (It should be minimum)
   * For more information see:
   *     J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
   *     Function Algorithms", Plenum Press, New York, 1981.
   * 
   */

  double H(const xmippFCB& _xmippDS) const; 


  /**
   * Calculates Non-fuzzy index (NFI) validity functional
   * @param _xmippDS Data structure to train, a codeBook in this case
   * (It should be maximum)
   * For more information see:
   *     M. Roubens, "Pattern Classification Problems and Fuzzy Sets",
   *     Fuzzy Sets and Systems, 1:239-253, 1978.
   * 
   */

  double NFI(const xmippFCB& _xmippDS) const;


  /**
   * Calculates Compactness and separation index (S) validity functional
   * @param _xmippDS Data structure to train, a codeBook in this case
   * @param _examples  A training set with the training examples
   * (It should be minimum)
   * For more information see:
   *     X.L. Xie and G. Beni, "A Validity Measure for Fuzzy Clustering",
   *     IEEE Trans. PAMI, 13(8):841-847, 1991.
   * 
   */
  double S(const xmippFCB& _xmippDS, const TS& _examples) const;

protected:

  /// print itself on standard output
  void printSelf( ostream& _os ) const ;

  double m;			    // Fuzzy constant
  double epsilon;		    // Stopping criterion Error < epsilon
  unsigned epochs;                  // Number of presentations of the whole sample
};

//@}

#endif//_XMIPPFCMEANS_H
