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
// xmippSOM.cc
// Implements Smoothly Distributed Fuzzy c-means Self-Organizing Map
//-----------------------------------------------------------------------------

#include "../xmippFuzzySOM.hh"



//-----------------------------------------------------------------------------
  /**
   * Sets the number of training steps
   * @param _nSteps  Number of training steps
   */
  void xmippFuzzySOM::nSteps(const unsigned long& _nSteps) { somNSteps = _nSteps; };

//-----------------------------------------------------------------------------
  /**
   * Sets the Intial Fuzzy membership
   * @param _m0  
   */
  void xmippFuzzySOM::initialFuzzzyMembership(const double& _m0) { m0 = _m0; };


//-----------------------------------------------------------------------------
  /**
   * Sets the Final Fuzzy membership
   * @param _m1  
   */
  void xmippFuzzySOM::finalFuzzzyMembership(const double& _m1) { m1 = _m1; };


//-----------------------------------------------------------------------------
  /**
   * Sets the number of deterministic annealing training steps
   * @param _annSteps  Number of steps
   */
  void xmippFuzzySOM::setAnnSteps(const unsigned long& _annSteps) {annSteps = _annSteps;};


//-----------------------------------------------------------------------------
  /**
   * Sets the Regularization Constant
   * @param _reg  
   */
  void xmippFuzzySOM::regularization(const double& _reg) { reg = _reg;};


//-----------------------------------------------------------------------------
  /**
   * Trains the Fuzzy SOM
   * @param _som  The fuzzy som to train           
   * @param _ts   The training set
   */
  void xmippFuzzySOM::train (xmippFuzzyMap& _som, const TS& _examples)
  {

    numNeurons = _som.size();
    numVectors = _examples.size();
    dim = _examples.theItems[0].size();
    tmpV.resize(dim, 0.);
    tmpDens.resize(numNeurons, 0.);
    tmpMap.resize(numNeurons);
    for (int i = 0; i < numNeurons; i++ )
        tmpMap[i].resize(dim, 0.);

    tmpD.resize(numNeurons);
    double stopError;

    int verbosity = listener->getVerbosity();
    if (verbosity)
  	listener->OnReportOperation((string) "\nTraining Fuzzy SOM....\n");  

    // Deterministic annealing step
    if (annSteps) {

      for (int i = 0; i <  _examples.size(); i++) 
        for (int j = 0; j < _som.size(); j++) 
          _som.memb[i][j] = 0.;

      if (verbosity)
  	listener->OnReportOperation((string) "\nDeterministic annealing steps....\n");  
      if (verbosity == 1 || verbosity == 3)    
        listener->OnInitOperation(annSteps);
    
      double fuzz;
      for (int iter =0; iter < annSteps; iter++) {
      
         fuzz = m0-iter*(m0-m1)/(annSteps-1);
         stopError = updateU(_som, _examples, fuzz);
         updateV(_som, _examples, fuzz);	
	 if (verbosity == 1 || verbosity == 3)    
      	  listener->OnProgress(iter);
         if (verbosity >= 2) {
	   char s[100];
	   sprintf(s, "Iteration %d of %d. Code vectors variation: %g\n", iter + 1, annSteps, stopError);
  	   listener->OnReportOperation((string) s);
         }	
      }
      if (verbosity == 1 || verbosity == 3)    
      	 listener->OnProgress(annSteps);
    }
 
    
    // Fixed steps
    

    for (int i = 0; i <  _examples.size(); i++) 
      for (int j = 0; j < _som.size(); j++) 
         _som.memb[i][j] = 0.;
    

    if (verbosity)
      listener->OnReportOperation((string) "\nFinal steps....\n");  
    if (verbosity == 1 || verbosity == 3)    
      listener->OnInitOperation(somNSteps);

    unsigned t = 0;  // Iteration index
    stopError = epsilon + 1; // Initially must be higher

    while (( stopError > epsilon ) && (t < somNSteps)) {
      stopError = updateU(_som, _examples, m1);
      updateV(_som, _examples, m1);
      t++;
      if (verbosity == 1 || verbosity == 3)    
      	listener->OnProgress(t);
      if (verbosity >= 2) {
	char s[100];
	sprintf(s, "Iteration %d of %d. Code vectors variation: %g\n", t, somNSteps, stopError);
  	listener->OnReportOperation((string) s);
      }	
    } // while
    if (verbosity == 1 || verbosity == 3)    
      listener->OnProgress(somNSteps);


    tmpV.clear();
    tmpDens.clear();
    tmpMap.clear();
    tmpD.clear();

  }; 
  

//-----------------------------------------------------------------------------
  /**
   * Tests the Fuzzy SOM
   * @param _som        The fuzzy som to test
   * @param _examples   The training set of examples
   */
  double xmippFuzzySOM::test (const xmippFuzzyMap& _som, const TS& _examples) const
  {

	// Defines verbosity level
        int verbosity = listener->getVerbosity();
        if (verbosity) {
  	  listener->OnReportOperation((string) "\nEstimating quantization error....\n");  
          listener->OnInitOperation(_examples.size());
	}


  	/* Scan all data entries */
  	double qerror = 0.0;
        for (int i = 0; i < _examples.size(); i++) {
           SomIn& theBest = _som.fuzzyTest(i); // get the best
   	   qerror += (double) eDist(theBest, _examples.theItems[i]);
	   if (verbosity) {
	   	if ((i % (int)((_examples.size()*5)/100)) == 0)
	   		listener->OnProgress(i);
	   }
	}
	if (verbosity)listener->OnProgress(_examples.size());
	return (qerror/(double) _examples.size());
	
  };


//-----------------------------------------------------------------------------
  /**
   * Update Fuzzy Memberships
   */
  double xmippFuzzySOM::updateU(xmippFuzzyMap& _som, const TS& _examples, const double& _m)
  {

    // Create auxiliar stuff

    unsigned i, j, k;
    double var = 0;
    double auxDist, auxProd, auxExp;    

    auxExp = 1./(_m-1.);  
    // Update Membership matrix

    for (k = 0; k < numVectors; k++ ) {     
      auxProd = 0;
      for ( i = 0; i < numNeurons; i ++ ) {
          auxDist = 0;
          for (j = 0; j < dim; j++)
	     auxDist += (double) ( (double)(_examples.theItems[k][j]) - (double)(_som.theItems[i][j]))*((double)(_examples.theItems[k][j]) - (double)(_som.theItems[i][j]));
          auxDist = pow(auxDist, auxExp);
          if (!finite(auxDist)) auxDist = MAXFLOAT;
          if (auxDist < MAXZERO) auxDist = MAXZERO;
          auxProd += (double) 1./auxDist;
          tmpD[i] = auxDist;
      }    
      for (j = 0; j < numNeurons; j ++ ) {
          double tmp =  1./(auxProd*tmpD[j]);
          var += fabs((double)(_som.memb[k][j]) - tmp);
         _som.memb[k][j] = (xmippFeature) tmp;
      }


    } // for k

    var /= (double) numNeurons*numVectors;

    return var;

  }


//-----------------------------------------------------------------------------
  /**
   * Update Fuzzy Code vectors
   */
  void xmippFuzzySOM::updateV(xmippFuzzyMap& _som, const TS& _examples, const double& _m)
  {
     unsigned j, cc, vv;
     unsigned t2 = 0;  // Iteration index

     // Calculate Temporal scratch values
     for (cc = 0; cc < numNeurons; cc++ ) {
	for (j = 0; j < dim; j++) tmpMap[cc][j] = 0.;
	tmpDens[cc] = reg;
     	for(vv = 0; vv < numVectors; vv++) {
     	   double tmpU = (double) pow((double)(_som.memb[vv][cc]), (double)_m);
     	   tmpDens[cc] += tmpU; 
     	   for (j = 0; j < dim; j++) {
	     tmpMap[cc][j] += (double)((double) tmpU*(double)(_examples.theItems[vv][j]));
	   }
     	}    
     }


     // Update Code vectors using a sort of Gauss-Seidel iterative algorithm.
     // Usually 100 iterations are enough.


     double stopError2 = 1;
     while (( stopError2 > 1e-5 ) && (t2 < 100)) {
        t2++; stopError2 = 0;
        for(cc = 0; cc < numNeurons; cc++) {
	    if (reg != 0)
               _som.localAve(_som.indexToPos(cc), tmpV);
             for (j = 0; j < dim; j++) {
	        tmpV[j] *= reg; 
		double tmpU = (tmpMap[cc][j] + tmpV[j])/tmpDens[cc];
		stopError2 += fabs((double)(_som.theItems[cc][j]) - tmpU);
		_som.theItems[cc][j] = (xmippFeature) tmpU;  
	     }
	} // for
	stopError2 /= (double) (numNeurons*dim);
     } // while	

  }


//-----------------------------------------------------------------------------

  /**
   * Determines the functional value
   * Returns the fidelity to the data and penalty parts of the functional
   */
  double xmippFuzzySOM::functional(const TS& _examples, const xmippFuzzyMap& _som, double _m, double _reg, double& _fidelity, double& _penalty)
  {
    unsigned j, vv, cc;
    double t, t1, t2 = 0;
    _fidelity = 0;
    for (vv = 0; vv < numVectors; vv++ ) {
    	for (cc = 0; cc < numNeurons; cc++ ) {
	   t = 0;
	   for (int j = 0; j < dim; j++)
	     t += ((double)(_examples.theItems[vv][j]) - (double)(_som.theItems[cc][j]))*((double)(_examples.theItems[vv][j]) - (double)(_som.theItems[cc][j]));
           t1 = (double) pow((double)_som.memb[vv][cc], (double) _m);
	   t2 += t1;
	   _fidelity += t1*t;
	}
    }
    _fidelity /= t2; 
    _penalty = 0;
    
    if (reg != 0) {
      for (cc = 0; cc < numNeurons; cc++ ) {
         _som.localAve(_som.indexToPos(cc), tmpV);
         for (j = 0; j < dim; j++) {
            _penalty += (_som.theItems[cc][j] - tmpV[j])*(_som.theItems[cc][j] - tmpV[j]); 
         }
      }
      _penalty /= (double) (numNeurons);
    }  

    return (_fidelity - _reg*_penalty);
  }


//-----------------------------------------------------------------------------
void xmippFuzzySOM::showX(const TS& _ts) {
  cout << "Data (1..nd, 1..nv) "  << endl;
  for (int i = 0; i < _ts.size(); i++) {
     for (int j = 0; j < _ts.theItems[0].size(); j++) {
       cout << i + 1 << "  " << j + 1 << "  " << _ts.theItems[i][j] << endl;
     }  
  }
}


//-----------------------------------------------------------------------------
void xmippFuzzySOM::showV(xmippFuzzyMap& _som) {
  cout << "Code vectors (1..ni, 1..nj, 1..nv) "  << endl;
  for (int i = 0; i < _som.size(); i++) {
     int tmpj = _som.indexToPos(i).first;
     int tmpi = _som.indexToPos(i).second;
     for (int j = 0; j < _som.theItems[0].size(); j++)
       cout << tmpi + 1 << "  " << tmpj + 1 << "  " << j + 1 << "  " << _som.theItems[i][j] << endl;
  }
}

//-----------------------------------------------------------------------------
void xmippFuzzySOM::showU(xmippFuzzyMap& _som, const TS& _ts) {
  cout << " Memberships (1..nd,1..ni,1..nj)" << endl;
  for (int i = 0; i <  _ts.size(); i++) {
     for (int j = 0; j < _som.size(); j++) {
       int tmpj = _som.indexToPos(j).first;
       int tmpi = _som.indexToPos(j).second;
       cout << i +1 << "  " << tmpi + 1 << "  " << tmpj + 1 << "  " << _som.memb[i][j] << endl;
     }
  }


}
