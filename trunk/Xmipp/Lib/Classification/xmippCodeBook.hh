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
// xmippCB.h
//-----------------------------------------------------------------------------

#ifndef XMIPPCODEBOOK_H
#define XMIPPCODEBOOK_H

//-----------------------------------------------------------------------------

#include <strstream>
#include <map>

#include "xmippCDSet.hh"
#include "xmippCDataTypes.hh"
#include "xmippCTVectors.hh"
#include "xmippVectorOps.hh"
#include "xmippDistances.hh"
#include "xmippUniform.hh"

//-----------------------------------------------------------------------------


/**@name Code Book*/
//@{
/**
 * This class implements a codebook.
 * A codebook is a set of examples (each of them being a vector). These
 * examples are usually labeled, ie, classified (the codebook is calibrated),
 * but it is not necessary (the codebook is NOT calibrated). A codebook is the
 * result obtained after one has trained some kind of algorithms. The way to
 * classify a data once the algorithm has been trained is to look for the
 * example in the code book that best matches the data to classify. Then, the
 * same label of the example in the codebook is associated to the data wanted
 * to be classified (if it is a calibrated codebook), or the example itself is
 * returned (if it is a NO calibrated codebook) indicating with this that the
 * data belongs to the same 'class' that the returned example.
 */
class xmippCB : public xmippCDSet<xmippVector, xmippLabel>, public xmippCTSet<xmippVector, xmippLabel>
{
 public:
  
    vector< vector <unsigned> > classifVectors;	
    vector< double > aveDistances;	


  /**
   * Default constructor
   * @param _calib   Calibrated or not, that is, a CB with class labels or not
   */
  xmippCB(const bool& _calib = false) : xmippCDSet<xmippVector, xmippLabel>(),
                                      xmippCTSet<xmippVector, xmippLabel>(_calib) {};


  /**
   * Constructor.
   * Constructs a codebook with initial code vectors at zero.
   * from an unsigned integer to instantiate the template
   * @param _n       Number of vectors
   * @param _size    Size of code vectors
   * @param _cal     Calibrated or not, that is, a CB with class labels or not
   */
  xmippCB (unsigned _n, unsigned _size, bool _cal = false );


  /**
   * Constructor.
   * Constructs a codebook with random initial code vectors.
   * from an unsigned integer to instantiate the template
   * @param _n       Number of vectors
   * @param _size    Size of code vectors
   * @param _lower   Lower value for random elements
   * @param _upper   Upper value for random elements
   * @param _cal     Calibrated or not, that is, a CB with class labels or not
   */
  xmippCB (unsigned _n, unsigned _size, xmippFeature _lower = 0, xmippFeature _upper = 1,
         bool _cal = false );

  /**
   * Constructor.
   * Constructs a codebook with initial code vectors taken randomly from
   * the training file.
   * from an unsigned integer to instantiate the template
   * @param _n       Number of vectors
   * @param _ts	     Training set; will be used to get initial values
   * @param _use_rand_cvs  Use random code vector values
   */
  xmippCB (unsigned _n, const xmippCTVectors& _ts, const bool _use_rand_cvs);

  /**
   * Constructs a code book given a stream
   * @param _is  The input stream
   * @exception  runtime_error  If there are problems with the stream
   */
  xmippCB(istream& _is); 

  /**
   * Virtual destructor needed
   */
  virtual ~xmippCB() {};

  /**
   * Returns the code vector that represents the input in the codebook
   * @param _in  Sample to classify
   */
  virtual xmippVector& test(const xmippVector& _in) const;

  /**
   * Returns the index to the code vector that represents the input in the codebook
   * @param _in  Sample to classify
   */
  virtual unsigned testIndex(const xmippVector& _in) const;


  /**
   * Returns the index of the codevector closest to an input.
   * This is the method used to classify inputs
   * @param _ts  Training set
   * @param _in  Index to the Sample to be classified
   */
  virtual unsigned winner(const xmippCTVectors& _ts, unsigned _in) const;


  /**
   * Fills the classifVectors with the list of the best input vectors associated to it.
   * @param _ts  Sample list to classify
   */

  virtual void classify(const xmippCTVectors* _ts); 


  /**
   * Returns the list of input vectors associated to this code vector.
   * @_index  code vector index
   */
  virtual const vector< unsigned>& classifAt(const unsigned& _index) const;

   /**
   * Returns the number of input vectors associated to this code vector.
   * @_index  code vector index
   */
  virtual unsigned classifSizeAt(const unsigned& _index) const;

  /**
   * Returns the label associated to an input
   * @param _in  Sample to classify
   */
  virtual xmippLabel apply(const xmippVector& _in) const;

  /**
   * Calibrates the code book
   * @param _ts   The calibrated training set
   * @param _def  Default target for non-calibrated vectors
   * @exception runtime_error  If the training set is not calibrated
   */
  virtual void calibrate (xmippCTVectors& _ts, xmippLabel _def = "");

   /**
   * Returns the index of the codevector closest to an input.
   * This is the method used to classify inputs
   * @param _in  Sample to classify.
   */
  virtual unsigned output(const xmippVector& _in) const;


  /**
   * Standard output for a code book
   * @param _os The output stream
   */
  virtual void printSelf(ostream& _os) const; 

  /**
   * Standard input for a code book
   * @param _is The input stream
   */
  virtual void readSelf (istream& _is);

  /**
   * Saves the xmippCodeBook class into a stream. 
   * this method can be used to save the status of the class.
   * @param _os The output stream
   */
  virtual void saveObject(ostream& _os) const;


  /**
   * Loads the xmippCodeBook class from a stream. 
   * this method can be used to load the status of the class.
   * @param _is The output stream
   */
  virtual void loadObject(istream& _is);

  /** 
   *	Normalize all features in the codebook 
   * 	@param _varStats The normalization information
   */

  virtual void xmippCB::Normalize(const vector <xmippCTVectors::statsStruct>&  _varStats);


  /** 
   *	UnNormalize all features in the codebook 
   * 	@param _varStats The normalization information
   */

  virtual void xmippCB::unNormalize(const vector <xmippCTVectors::statsStruct>&  _varStats);


  /**
   * Prints the histogram values of each codevector. 
   * @param _os  The the output stream
   */
  virtual void printHistogram(ostream& _os) const;

  /**
   * Prints the Average Quantization Error of each codevector. 
   * @param _os  The the output stream
   */
  virtual void printQuantError(ostream& _os) const;


protected:

  /** 
   * Reads the classif vectors from a stream.
   * @param _is  The input stream
   */
    void readClassifVectors(istream& _is);

  /** 
   * Writes the classif vectors to a stream
   * @param _os  The output stream
   */
    void writeClassifVectors(ostream& _os) const;


};
//@}
//-----------------------------------------------------------------------------

#endif//XMIPPCODEBOOK_H
