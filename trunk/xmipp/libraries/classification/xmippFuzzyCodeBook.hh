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
// xmippFuzzyCodeBook.hh
// Implements a set of fuzzy code vectors, that is a Fuzzy code book.
//-----------------------------------------------------------------------------

#ifndef XMIPPFUZZYCB_H
#define XMIPPFUZZYCB_H

// To avoid problems with long template names
#pragma warning(disable:4786)

//-----------------------------------------------------------------------------

#include "xmippCodeBook.hh"

//-----------------------------------------------------------------------------

/**@name Fuzzy Code Book*/
//@{
/**
 * This class implements an specific Data Structure for Fuzzy 
 * algorithms. It inherits fron xmippCB, so this Data Structure is basically 
 * a codebook with a Fuzzy-membership Matrix.
 * Some of the methods from xmippCB are redefined. 
 */
class xmippFCB : public xmippCB
{
 public:
  
  typedef vector< vector< xmippFeature > > MM;       
  /// Alias for Membership Matrix
  //typedef xmippCTSet<vector<xmippFeature>,xmippLabel> TS;  // Alias for a Training set 
  typedef xmippCTVectors TS;  
  /// Alias for a Training set 
  typedef xmippCTSet<vector<xmippFeature>, xmippFeature> FV;    
  /// Alias for Fuzzy vectors 
   
   // Fuzzy membership matrix
      MM memb;  
  
  /**
   * Default constructor
   * @param _calib     Calibrated or not, that is, a CB with class labels or not
   */
  xmippFCB(const bool& _calib = false) : xmippCB(_calib) {};


   
  /**
   * Constructor.
   * Constructs a fuzzy codebook with initial code vectors filled with zero.
   * @param _n       Number of code vectors
   * @param _size    Size of code vectors
   * @param _data    Size of the training set   
   * @param _cal     Calibrated or not, that is, a CB with class labels or not
     It calls Base Class constructor (xmippCB) 
   */
   

  xmippFCB (unsigned _n, unsigned _size, unsigned _data, bool _cal = false );

   
  /**
   * Constructor.
   * Constructs a fuzzy codebook with random initial code vectors.
   * @param _n       Number of code vectors
   * @param _size    Size of code vectors
   * @param _data    Size of the training set
   * @param _lower   Lower value for random elements
   * @param _upper   Upper value for random elements
   * @param _cal     Calibrated or not, that is, a CB with class labels or not
     It calls Base Class constructor (xmippCB) 
   */
   
  xmippFCB (unsigned _n, unsigned _size, unsigned _data, double _lower = 0, double _upper = 1, bool _cal = false );

  /**
   * Constructor.
   * Constructs a fuzzy codebook with initial code vectors taken randomly from
   * the training file.
   * @param _n       Number of vectors
   * @param _ts	     Training set; will be used to get initial values
   * @param _use_rand_cvs  Use random code vectors (inherited from base class)
     It calls Base Class constructor (xmippCB) 
   */

  xmippFCB (unsigned _n, const xmippCTVectors& _ts, const bool _use_rand_cvs=false);

  /*
   * Constructs a fuzzy code book given a stream
   * @param _is  The input stream
   * @param _size Size of code vectors (number of data points)
   */

  xmippFCB(istream& _is, const unsigned _size = 0);
  
  /**
   * Virtual destructor 
   */

  virtual ~xmippFCB() {};

  
 /** Fuctions to access the Fuzzy Membership Matrix
 */ 
  
  /**
   * Returns a const reference to the specified item
   * @param _ci  cluster index
   * @param _di  data index
   * @exception out_of_range If _i is out of range
   */
  xmippFeature membAt ( unsigned _di, unsigned _ci) const;

  /**
   * Returns a  reference to the specified item
   * @param _ci  cluster index
   * @param _di  data index
   * @exception out_of_range If _i is out of range
   */
  xmippFeature& membAt ( unsigned _di, unsigned _ci); 


  /**
   * Returns dimensions of the Membership matrix
   */
  unsigned membClusters () const;
  unsigned membVectors () const;


  /**
   * Returns the code vector that represents the input in the codebook
   * @param _in  Sample to classify
     Note: The difference between Fuzzy codevector and non-Fuzzy 
     codevector is that the best (winner) is estimated using the 
     fuzzy membership matrix.
   */
  virtual xmippVector& fuzzyTest(unsigned _in) const;

  /**
   * Returns the index to the code vector that represents the input in the codebook
   * @param _in  Sample to classify
     Note: The difference between Fuzzy codevector and non-Fuzzy 
     codevector is that the best (winner) is estimated using the 
     fuzzy membership matrix.
   */
  virtual unsigned fuzzyTestIndex(unsigned _in) const;

  /**
   * Returns the label associated to an input
   * @param _in  Index to the sample to be classified
   */
  virtual xmippLabel fuzzyApply(unsigned _in) const;

  /**
   * Calibrates the code book
   * @param _ts   The calibrated training set
   * @param _def  Default target for non-calibrated vectors
   * @exception runtime_error  If the training set is not calibrated
   */
  virtual void fuzzyCalibrate (xmippCTVectors& _ts, xmippLabel _def = xmippLabel());

  /**
   * Returns the index of the codevector closest to an input.
   * This is the method used to classify inputs
   * @param _in  Index to the Sample to be classified
   */
  virtual unsigned fuzzyWinner(unsigned _in) const;

   /**
   * Returns the index of the codevector closest to an input.
   * This is the method used to classify inputs
   * @param _in  Index to the Sample to be classified
   */
  virtual unsigned fuzzyOutput(unsigned _in) const;

  /**
   * Fills the classifVectors with the list of the best input vectors associated to it.
   * In this case, it uses the Fuzzy Memberships to make the assignments
   * @param _ts  Sample list to classify
   */
  virtual void classify(const xmippCTVectors* _ts); 

   /**
   * Hard partition the Fuzzy membership matrix
   * using the Nearest Maximum Membership conversion
   */

  virtual void hardPartition();


  /**
   * Returns the alpha-core set (also called alpha-level set or simply "core")
   * @param _ts       The training set
   * @param _alpha    A threshold to identify the core.
   * @param _cluster  The cluster or partition
   * Returns a training vector "contained" in the cluster
   */
  virtual TS alphaCore(TS _ts, double _alpha, unsigned _cluster) const;

  /**
   * Writes the membership values
   * @param _os  The output stream
   */
  virtual void writeMembership (ostream& _os) const;

  /**
   * Reads the membership values
   * @param _is  The input stream
   */
  virtual void readMembership (istream& _is);


  /**
   * Saves the xmippFuzzyCodeBook class into a stream. 
   * this method can be used to save the status of the class.
   * @param _os The output stream
   */
  virtual void saveObject(ostream& _os) const;


  /**
   * Loads the xmippFuzzyCodeBook class from a stream. 
   * this method can be used to load the status of the class.
   * @param _is The output stream
   */
  virtual void loadObject(istream& _is);


  /**
   * Constructs a fuzzy code book given a stream
   * @param _is  The input stream
   * @param _size Size of code vectors (number of data points)
   * @exception  runtime_error  If there are problems with the stream
   */
  virtual void readSelf(istream& _is, const unsigned _size = 0);

  /**
   * Prints the density values of each Fuzzy codevector. 
   * @param _os  The the output stream
   */
  virtual void printDensity(ostream& _os) const;


 private:
 
   // Dimensions of the membership matrix   
      unsigned numClusters, numVectors;
    
};
//@}


//-----------------------------------------------------------------------------

#endif//XMIPPFUZZYCB_H
