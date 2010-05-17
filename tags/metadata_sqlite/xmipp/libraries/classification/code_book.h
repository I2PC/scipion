/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
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

//-----------------------------------------------------------------------------
// xmippCB.h
//-----------------------------------------------------------------------------

#ifndef XMIPPCODEBOOK_H
#define XMIPPCODEBOOK_H

#include <sstream>
#include <map>

#include "data_set.h"
#include "data_types.h"
#include "training_vector.h"
#include "vector_ops.h"
#include "distance.h"
#include "uniform.h"

/**@defgroup CodeBook Code book
   @ingroup ClassificationLibrary */
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

    std::vector< std::vector <unsigned> > classifVectors;
    std::vector< double > aveDistances;


    /**
     * Default constructor
     * Parameter: _calib   Calibrated or not, that is, a CB with class labels or not
     */
    xmippCB(const bool& _calib = false) : xmippCDSet<xmippVector, xmippLabel>(),
            xmippCTSet<xmippVector, xmippLabel>(_calib)
    {};


    /**
     * Constructor.
     * Constructs a codebook with initial code vectors at zero.
     * from an unsigned integer to instantiate the template
     * Parameter: _n       Number of vectors
     * Parameter: _size    Size of code vectors
     * Parameter: _cal     Calibrated or not, that is, a CB with class labels or not
     */
    xmippCB(unsigned _n, unsigned _size, bool _cal = false);


    /**
     * Constructor.
     * Constructs a codebook with random initial code vectors.
     * from an unsigned integer to instantiate the template
     * Parameter: _n       Number of vectors
     * Parameter: _size    Size of code vectors
     * Parameter: _lower   Lower value for random elements
     * Parameter: _upper   Upper value for random elements
     * Parameter: _cal     Calibrated or not, that is, a CB with class labels or not
     */
    xmippCB(unsigned _n, unsigned _size, xmippFeature _lower = 0, xmippFeature _upper = 1,
            bool _cal = false);

    /**
     * Constructor.
     * Constructs a codebook with initial code vectors taken randomly from
     * the training file.
     * from an unsigned integer to instantiate the template
     * Parameter: _n       Number of vectors
     * Parameter: _ts      Training set; will be used to get initial values
     * Parameter: _use_rand_cvs  Use random code vector values
     */
    xmippCB(unsigned _n, const xmippCTVectors& _ts, const bool _use_rand_cvs);

    /**
     * Constructs a code book given a stream
     * Parameter: _is  The input stream
     * @exception  runtime_error  If there are problems with the stream
     */
    xmippCB(std::istream& _is);

    /**
     * Virtual destructor needed
     */
    virtual ~xmippCB()
    {};

    /**
     * Returns the code vector that represents the input in the codebook
     * Parameter: _in  Sample to classify
     */
    virtual xmippVector& test(const xmippVector& _in) const;

    /**
     * Returns the index to the code vector that represents the input in the codebook
     * Parameter: _in  Sample to classify
     */
    virtual unsigned testIndex(const xmippVector& _in) const;


    /**
     * Returns the index of the codevector closest to an input.
     * This is the method used to classify inputs
     * Parameter: _ts  Training set
     * Parameter: _in  Index to the Sample to be classified
     */
    virtual unsigned winner(const xmippCTVectors& _ts, unsigned _in) const;


    /**
     * Fills the classifVectors with the list of the best input vectors associated to it.
     * Parameter: _ts  Sample list to classify
     */

    virtual void classify(const xmippCTVectors* _ts);


    /**
     * Returns the list of input vectors associated to this code vector.
     */
    virtual const std::vector< unsigned>& classifAt(const unsigned& _index) const;

    /**
    * Returns the number of input vectors associated to this code vector.
    */
    virtual unsigned classifSizeAt(const unsigned& _index) const;

    /**
     * Returns the label associated to an input
     * Parameter: _in  Sample to classify
     */
    virtual xmippLabel apply(const xmippVector& _in) const;

    /**
     * Calibrates the code book
     * Parameter: _ts   The calibrated training set
     * Parameter: _def  Default target for non-calibrated vectors
     * @exception runtime_error  If the training set is not calibrated
     */
    virtual void calibrate(xmippCTVectors& _ts, xmippLabel _def = "");

    /**
    * Returns the index of the codevector closest to an input.
    * This is the method used to classify inputs
    * Parameter: _in  Sample to classify.
    */
    virtual unsigned output(const xmippVector& _in) const;


    /**
     * Standard output for a code book
     * Parameter: _os The output stream
     */
    virtual void printSelf(std::ostream& _os) const;

    /**
     * Standard input for a code book
     * Parameter: _is The input stream
     */
    virtual void readSelf(std::istream& _is, long _dim = -1, long _size = -1);

    /**
     * Saves the xmippCodeBook class into a stream.
     * this method can be used to save the status of the class.
     * Parameter: _os The output stream
     */
    virtual void saveObject(std::ostream& _os) const;


    /**
     * Loads the xmippCodeBook class from a stream.
     * this method can be used to load the status of the class.
     * Parameter: _is The output stream
     */
    virtual void loadObject(std::istream& _is);

    /**
     * Normalize all features in the codebook
     *  Parameter: _varStats The normalization information
     */

    virtual void Normalize(const std::vector <xmippCTVectors::statsStruct>&  _varStats);


    /**
     * UnNormalize all features in the codebook
     *  Parameter: _varStats The normalization information
     */

    virtual void unNormalize(const std::vector <xmippCTVectors::statsStruct>&  _varStats);


    /**
     * Prints the histogram values of each codevector.
     * Parameter: _os  The the output stream
     */
    virtual void printHistogram(std::ostream& _os) const;

    /**
     * Prints the Average Quantization Error of each codevector.
     * Parameter: _os  The the output stream
     */
    virtual void printQuantError(std::ostream& _os) const;


protected:

    /**
     * Reads the classif vectors from a stream.
     * Parameter: _is  The input stream
     */
    void readClassifVectors(std::istream& _is);

    /**
     * Writes the classif vectors to a stream
     * Parameter: _os  The output stream
     */
    void writeClassifVectors(std::ostream& _os) const;
};
//@}
#endif//XMIPPCODEBOOK_H
