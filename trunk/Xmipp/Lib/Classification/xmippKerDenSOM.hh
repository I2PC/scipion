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
// xmippKerDenSOM.hh
// Implements Smoothly Distributed Kernel Probability Density Estimator Self-Organizing Map
// This is an abstract base class for different variants of the KerDenSOM algorithm
//-----------------------------------------------------------------------------

#ifndef XMIPPKERDENSOM_H
#define XMIPPKERDENSOM_H

#include "xmippBaseAlgo.hh"
#include "xmippMap.hh"


//---------------------------------------------------------------------------


/**@name Smoothly Distributed Kernel Probability Density Estimator Self Organizing Map*/
//@{
/**
 * This class trains a Smoothly Distributed Kernel Probability Density Estimator Self Organizing Map
 */
class xmippKerDenSOM : public xmippBaseAlgo<xmippFuzzyMap>
{
 public:

  /**
   * Constructs the algorithm
   * @param _reg0       Initial regularization factor
   * @param _reg1       Final regularization factor
   * @param _annSteps   Number of steps in deterministic annealing
   * @param _epsilon    Stopping criterion
   * @param _nSteps     Number of training steps
   */
  xmippKerDenSOM(double _reg0, double _reg1, unsigned long _annSteps, 
  		double _epsilon, unsigned long _nSteps)
  : xmippBaseAlgo<xmippFuzzyMap>(), reg0(_reg0), reg1(_reg1), annSteps(_annSteps), 
    epsilon(_epsilon), somNSteps(_nSteps) {};

  /**
   * Virtual destructor
   */
  virtual ~xmippKerDenSOM() {};

  /**
   * Sets the number of training steps
   * @param _nSteps  Number of training steps
   */
  void nSteps(const unsigned long& _nSteps);

  /**
   * Gets the Kernel Width
   */
  virtual double getSigma();


  /**
   * Sets the number of deterministic annealing training steps
   * @param _annSteps  Number of steps
   */
  void setAnnSteps(const unsigned long& _annSteps);


  /**
   * Trains the KerDenSOM
   * @param _som  The KerDenSom to train           
   * @param _ts   The training set
   * @param _update True if uses _som as starting point for training.
   * @param _sigma If update = true, uses this sigma for the training.
   */
  virtual void train (xmippFuzzyMap& _som, TS& _examples, FileName& _fn_vectors, bool _update = false, double _sigma = 0) = 0;

  /**
   * Tests the KerDenSOM
   * @param _som        The KerDenSom to test
   * @param _examples   The testing set
   */
  virtual double test (const xmippFuzzyMap& _som, const TS& _examples) const;


  /**
   * Determines the functional value.
   * Returns the likelihood and penalty parts of the functional
   */
  virtual double functional(const TS* _examples, const xmippFuzzyMap* _som, double _sigma, double _reg, double& _likelihood, double& _penalty) = 0;

  /**
   * Determines the Random Approximation of GVC
   * (For determining the optimal Regularization Factor)
   */

  virtual double randApproxGVC(const TS* _examples, const xmippFuzzyMap* _som, double _dataSD, double _reg);


 protected:

   double sigma;        // Optimum Kernel Width
   unsigned long annSteps;    // number of deterministic annealing steps
   double reg1, reg0;   // Regularization factors
   double epsilon;      // Stopping criterion Error < epsilon
   unsigned long somNSteps;   // number of steps

   
   // Internal Scratch

   unsigned numNeurons;
   unsigned numVectors;
   unsigned dim;
   vector < vector<double> > tmpMap;
   vector<double> tmpD, tmpD1, tmpDens, tmpV;

    
  /** Declaration of virtual method */
   virtual void train (xmippFuzzyMap& _som, const TS& _examples) const {};

  /* Declaration of abstract methods */

   // Update Us
   virtual double updateU(xmippFuzzyMap* _som, const TS* _examples, const double& _sigma, double& _alpha) = 0;   	
   
   // Estimate Sigma II
   virtual double updateSigmaII(xmippFuzzyMap* _som, const TS* _examples, const double& _reg, const double& _alpha) = 0;

   // Estimate the PD (Method 1: Using the code vectors)
   virtual double codeDens(const xmippFuzzyMap* _som, const xmippVector* _example, double _sigma) const = 0;

   // Estimate the PD (Method 2: Using the data)
   virtual double dataDens(const TS* _examples, const xmippVector* _example, double _sigma) const = 0;


   /* Some other common methods */
   
   // Update Code vectors   
   virtual void updateV(xmippFuzzyMap* _som, const TS* _examples, const double& _sigma);

   // Main iterations
   virtual double mainIterations(xmippFuzzyMap* _som, const TS* _examples, double& _sigma, const double& _reg);   	

   // Special Initialization of the Us
   virtual void xmippKerDenSOM::initU(xmippFuzzyMap* _som);   

   // Special initialization of Code vectors   
   virtual void updateV1(xmippFuzzyMap* _som, const TS* _examples);

   // Special initialization of Membership Matrix   
   virtual void updateU1(xmippFuzzyMap* _som, const TS* _examples);

   // Estimate Sigma I
   virtual double updateSigmaI(xmippFuzzyMap* _som, const TS* _examples);
 
   // Some printing methods.
   void showX(const TS* _ts);
   void showV(xmippFuzzyMap* _som);
   void showU(xmippFuzzyMap* _som, const TS* _ts);
   void printV(xmippFuzzyMap* _som, const TS* _ts, FileName& _fname);

};

//@}


#endif
