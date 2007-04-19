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
// xmippFuzzySOM.hh
// Implements Smoothly Distributed Fuzzy c-means Self-Organizing Map
//-----------------------------------------------------------------------------

#ifndef XMIPPFUZZYSOM_H
#define XMIPPFUZZYSOM_H

#include "base_algorithm.h"
#include "map.h"

/**@name Smoothly Distributed Fuzzy c-means Self-Organizing Map algorithm*/
//@{
/**
 * This class trains a Smoothly Distributed Fuzzy c-Means Self Organizing Map
 */
class xmippFuzzySOM : public xmippBaseAlgo<xmippFuzzyMap>
{
 public:

  /**
   * Constructs the algorithm
   * @param _m0         Initial Fuzzy Membership constant
   * @param _m1         Final Fuzzy Membership constant
   * @param _annSteps   Number of steps in deterministic annealing
   * @param _reg        Regularization constant
   * @param _epsilon    Stopping criterion
   * @param _nSteps     Number of training steps
   */
  xmippFuzzySOM(double _m0, double _m1, unsigned long _annSteps,
  		double _reg, double _epsilon, unsigned long _nSteps)
  : xmippBaseAlgo<xmippFuzzyMap>(), m0(_m0), m1(_m1), annSteps(_annSteps),
    reg(_reg), epsilon(_epsilon), somNSteps(_nSteps) {};

  /**
   * Virtual destructor
   */
  virtual ~xmippFuzzySOM() {};

  /**
   * Sets the number of training steps
   * @param _nSteps  Number of training steps
   */
  void nSteps(const unsigned long& _nSteps);

  /**
   * Sets the Initial Fuzzy membership
   * @param _m0
   */
  void initialFuzzzyMembership(const double& _m0);


  /**
   * Sets the Final Fuzzy membership
   * @param _m1
   */
  void finalFuzzzyMembership(const double& _m1);


  /**
   * Sets the number of deterministic annealing training steps
   * @param _annSteps  Number of steps
   */
  void setAnnSteps(const unsigned long& _annSteps);



  /**
   * Sets the Regularization Constant
   * @param _reg
   */
  void regularization(const double& _reg);


  /**
   * Trains the Fuzzy SOM
   * @param _som  The som to train
   * @param _examples   The training set
   */
  virtual void train (xmippFuzzyMap& _som, const TS& _examples);

  /**
   * Tests the Fuzzy SOM
   * @param _som        The fuzzy som to test
   * @param _examples   The training set of examples
   */
  virtual double test (const xmippFuzzyMap& _som, const TS& _examples) const;

  /**
   * Determines the functional value
   * Returns the fidelity to the data and penalty parts of the functional
   */
  virtual double functional(const TS& _examples, const xmippFuzzyMap& _som, double _m, double _reg, double& _fidelity, double& _penalty);



 private:
   double m0;          	      // Initial Fuzzy Membership
   double m1;          	      // Final Fuzzy Membership
   unsigned long annSteps;    	      // number of deterministic annealing steps
   double reg;          	      // Regularization constant
   double epsilon;	      // Stopping criterion Error < epsilon
   unsigned long somNSteps;           // number of steps

  /** Declaration of virtual method */
   virtual void train (xmippFuzzyMap& _som, const TS& _examples) const {};

   double updateU(xmippFuzzyMap& _som, const TS& _examples, const double& _m);   	

   // Update Fuzzy Code vectors
   void updateV(xmippFuzzyMap& _som, const TS& _examples, const double& _m);

   // Local Average of the code vectors
   void LocalAve(xmippFuzzyMap& _som, xmippVector& tmpV, int tmpi, int tmpj);

   // Internal Scratch

   unsigned numNeurons;
   unsigned numVectors;
   unsigned dim;
   vector<double> tmpV, tmpD, tmpDens;
   vector < vector<double> > tmpMap;

void showX(const TS& _ts);

void showV(xmippFuzzyMap& _som);

void showU(xmippFuzzyMap& _som, const TS& _ts);


};
//@}



#endif
