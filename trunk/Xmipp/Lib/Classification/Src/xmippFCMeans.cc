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
// xmippFCMeans.cc
// Fuzzy c-means clustering algorithm
//-----------------------------------------------------------------------------

#include "../xmippFCMeans.hh"

#pragma warning(disable:4786)

/**  Ctor from stream
 * @param _is Must have the parameters in the same order than the previous ctor.
   ****** check out this ************
 */
/*xmippFCMeans::xmippFCMeans( istream& _is )
  :xmippBaseAlgo< xmippFCB >( "xmippFCMeans") {
  _is >> m;
  _is >> epsilon;
  _is >> epochs;
};
*/

/**
 * Trains the Algorithm
 * @param _xmippDS Data structure to train, a codeBook in this case
 * @param _examples  A training set with the training examples
 */

void xmippFCMeans::train(xmippFCB& _xmippDS, TS& _examples) const {

    // Defines verbosity
    
    int verbosity = listener->getVerbosity();
    if (verbosity)
  	listener->OnReportOperation((string) "Training....\n");  
    if (verbosity == 1 || verbosity == 3)    
      listener->OnInitOperation(epochs);

    // Create auxiliar Codebook
    
     xmippFCB auxCB;  

  
    // Create auxiliar stuff

     unsigned numClusters = _xmippDS.size();
     unsigned numVectors = _examples.size();
     unsigned i, j, k;
     double stopError = 0, auxError = 0;
     double auxDist, auxProd, tmp, auxExp, auxSum;    
     unsigned t = 0;		// Iteration index
     xmippVector zero( _xmippDS.theItems[0].size() ) ;
     fill( zero.begin(), zero.end(), 0.0 );

     
    // Initialize auxiliary Codebook
     
     auxCB = _xmippDS;
     
    // Set auxExp
    
    auxExp = 2/(m-1);  
    
    // This is the main code of the algorithm. Iterates "epochs" times

  
    stopError = epsilon + 1; // Initially must be higher
    
    while (( stopError > epsilon ) && (t < epochs)) {

      // Update Membership matrix

      for (k = 0; k < numVectors; k++ ) {     

        auxProd = 1;
	for(j=0; j< numClusters; j++) 
	  auxProd *= (double) eDist( _xmippDS.theItems[j], _examples.theItems[k] );
	  
	if (auxProd == 0.) { // Apply k-means criterion (Data-CB) must be > 0 
	  for ( j = 0; j < numClusters; j ++ ) 
	    if (eDist( _xmippDS.theItems[j], _examples.theItems[k] ) == 0.) _xmippDS.memb[k][j] = 1.0;
	    else _xmippDS.memb[k][j] =  0.0;
	} else {
	  for ( i = 0; i < numClusters; i ++ ) {
            auxDist = 0;
	    for ( j = 0; j < numClusters; j ++ ) {
	      tmp = (double) eDist( _xmippDS.theItems[i], _examples.theItems[k] ) / 
	      	    (double) eDist( _xmippDS.theItems[j], _examples.theItems[k] );
	      auxDist += pow(tmp, auxExp);
	    } // for j
	    _xmippDS.memb[k][i] =  (xmippFeature) 1.0/auxDist; 
	  } // for i	
	} // if auxProd
	 

      } // for k
        
    
      // Update code vectors (Cluster Centers)

       for(i=0; i<numClusters; i++) {
	 _xmippDS.theItems[i] = zero;
	 auxSum = 0;
         for(k=0; k<numVectors; k++) {
	   _xmippDS.theItems[i] += (xmippFeature) pow((double)(_xmippDS.memb[k][i]), m) * _examples.theItems[k];
           auxSum += pow((double)(_xmippDS.memb[k][i]), m);
         } // for i
	 _xmippDS.theItems[i] /= (xmippFeature) auxSum;
       } // for k
		          
       
    
      // Compute stopping criterion

        stopError = 0; 
        for ( i = 0; i < numClusters; i ++ ) {
	  auxError = (double) eDist( _xmippDS.theItems[i], auxCB.theItems[i] );
          stopError += auxError * auxError;
        } // for i    

    
      // Update iteration index 

        t++;
	auxCB = _xmippDS; 

      if (verbosity == 1 || verbosity == 3)    
      	listener->OnProgress(t);
      if (verbosity >= 2) {
	char s[100];
	sprintf(s, "Iteration %d of %d. Code vectors variation: %g\n", t+1, epochs, stopError);
  	listener->OnReportOperation((string) s);
      }	
              
    } // while

    if (verbosity == 1 || verbosity == 3)    
    	listener->OnProgress(epochs);
    
}; // xmippFCMeans::train


/**
 * Test the Algorithm in a conventional way
 * @param _examples  A training set with the training examples
 */
  
double xmippFCMeans::test(const xmippFCB& _xmippDS, 
			const TS& _examples) const { 

  // Defines verbosity
    
  int verbosity = listener->getVerbosity();
  if (verbosity) {
      listener->OnReportOperation((string) "Testing....\n");  
      listener->OnInitOperation(_examples.size());
  }
  
  double distortion = 0;
  for ( unsigned i = 0; i < _examples.size(); i ++ ) { 
    const xmippVector& auxS = _examples.theItems[i];
    unsigned best = _xmippDS.output( auxS );
    distortion += (double) eDist( _xmippDS.theItems[best], _examples.theItems[i] );
    if (verbosity)    
      	listener->OnProgress(i);                  
  };
   
  if (verbosity)    
    listener->OnProgress(_examples.size());
 
  return distortion / (double) _examples.size();
}

  /**
   * Tests with the training set using for training.
   * @param _examples  The training set
   */

double xmippFCMeans::fuzzyTest(const xmippFCB& _xmippDS, 
			const TS& _examples) const { 

  // Defines verbosity
    
  int verbosity = listener->getVerbosity();
  if (verbosity) {
      listener->OnReportOperation((string) "Testing....\n");  
      listener->OnInitOperation(_examples.size());
  }

  double distortion = 0;
  for ( unsigned i = 0; i < _examples.size(); i ++ ) {
    unsigned best = _xmippDS.fuzzyOutput( i );
    distortion += (double) eDist( _xmippDS.theItems[best], _examples.theItems[i] );
    if (verbosity)    
      	listener->OnProgress(i);                  
  };
  if (verbosity)    
    listener->OnProgress(_examples.size());

  return distortion / (double)_examples.size();
}

  /**
   * Calculates Partition Coefficient (F) validity functional
   * @param _xmippDS Data structure to train, a codeBook in this case
   * 
   * Notes on F:  For U in Mfc (fuzzy partition space)
   *              1/C <= F <= 1
   *              for F = 1, U is hard (zeros and ones only)
   *              for F = 1/C, U = 1/C*ones(C,n);
   *
   * (max)
   *
   * For more information see:
   *     J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
   *     Function Algorithms", Plenum Press, New York, 1981.
   * 
   */

double xmippFCMeans::F(const xmippFCB& _xmippDS) const { 
     double F = 0.;
     for(unsigned k=0; k<_xmippDS.membVectors(); k++)
       for(unsigned i=0; i<_xmippDS.membClusters(); i++)
	   F += pow((double)(_xmippDS.memb[k][i]), 2);
     return (F/(double)(_xmippDS.membVectors()));   
}

  /**
   * Calculates Partition Entropy (H) validity functional
   * @param _xmippDS Data structure to train, a codeBook in this case
   *
   * Notes on H:  For U in Mfc
   *              0 <= H <= log(C)
   *              for H = 0, U is hard
   *              for H = log(C), U = 1/C*ones(C,n);
   *              0 <= 1 - F <= H (strict inequality if U not hard)
   *
   * (min)
   *
   * For more information see:
   *     J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
   *     Function Algorithms", Plenum Press, New York, 1981.
   * 
   */

double xmippFCMeans::H(const xmippFCB& _xmippDS) const { 
     double H = 0.;
     for(unsigned k=0; k<_xmippDS.membVectors(); k++)
       for(unsigned i=0; i<_xmippDS.membClusters(); i++)
         if (_xmippDS.memb[k][i] != 0.)
	   H += (double)(_xmippDS.memb[k][i])*log((double)(_xmippDS.memb[k][i]));
     return (-H/(double)(_xmippDS.membVectors()));   
}


  /**
   * Calculates Non-fuzzy index (NFI) validity functional
   * @param _xmippDS Data structure to train, a codeBook in this case
   * 
   * (max)
   * 
   * For more information see:
   *     M. Roubens, "Pattern Classification Problems and Fuzzy Sets",
   *     Fuzzy Sets and Systems, 1:239-253, 1978.
   * 
   */

double xmippFCMeans::NFI(const xmippFCB& _xmippDS) const { 
     double F = 0.;
     for(unsigned k=0; k<_xmippDS.membVectors(); k++)
       for(unsigned i=0; i<_xmippDS.membClusters(); i++)
	   F += pow((double)(_xmippDS.memb[k][i]), 2);
     
     double NFI = (((double)(_xmippDS.membClusters())* F - (double)(_xmippDS.membVectors())) / (double)(_xmippDS.membVectors())) / (double)(_xmippDS.membClusters()-1);
     return NFI;   
}


  /**
   * Calculates Compactness and separation index (S) validity functional
   * @param _xmippDS Data structure to train, a codeBook in this case
   * @param _examples  A training set with the training examples
   *
   * (min)
   *
   * For more information see:
   *     X.L. Xie and G. Beni, "A Validity Measure for Fuzzy Clustering",
   *     IEEE Trans. PAMI, 13(8):841-847, 1991.
   * 
   */

double xmippFCMeans::S(const xmippFCB& _xmippDS, 
			const TS& _examples) const {
     
     vector< vector< xmippFeature > > ICD;       // Intercluster distance
     vector< vector< xmippFeature > > D;         // Distance from each data to cluster centers
     
     unsigned i;
     D.resize(_xmippDS.membClusters());
     for(i=0; i<_xmippDS.membClusters(); i++) {
       vector <xmippFeature> d;
       d.resize(_xmippDS.membVectors());
       for(unsigned k=0; k<_xmippDS.membVectors(); k++)
           d[k] = eDist( _xmippDS.theItems[i], _examples.theItems[k] );
       D[i] = d; 
     } // for i     

     ICD.resize(_xmippDS.membClusters());
     for(i=0; i<_xmippDS.membClusters(); i++) {
       vector <xmippFeature> v;
       v.resize(_xmippDS.membVectors());
       for(unsigned j=0; j<_xmippDS.membClusters(); j++)
           v[j] = eDist( _xmippDS.theItems[i], _xmippDS.theItems[j] );
       ICD[i] = v;
     } // for i  	   

     
     xmippFeature auxSum = 0;
	 
     for( i=0; i<_xmippDS.membClusters(); i++) 
       for(unsigned k=0; k<_xmippDS.membVectors(); k++)
         auxSum += (xmippFeature) pow((double)(D[i][k] * _xmippDS.memb[k][i]), (double)m);
     
     xmippFeature auxMin = MAXFLOAT;
     for(i=0; i<_xmippDS.membClusters(); i++)
       for(unsigned j=i+1; j<_xmippDS.membClusters(); j++) 
         if (auxMin > ICD[i][j]) auxMin = ICD[i][j];
     
     double S = auxSum / (double) (_xmippDS.membVectors()) / (double) (auxMin);
     return S;   
     
}


/// print itself on standard output
void xmippFCMeans::printSelf( ostream& _os ) const {
  // Call base class, which will print ID
  _os << "Class (Algorithm): " << endl;
  xmippBaseAlgo<xmippFCB>::printSelf( _os );
  _os << endl;
  // Print parameters in the same order they are declared
  _os << "Fuzzy constant m = " << m << endl;
  _os << "Epsilon eps = " << epsilon << endl;
  _os << "Iterations iter = " << epochs << endl << endl;
} ;
