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
// xmippFuzzyCodeBook.hh
// Implements a set of fuzzy code vectors, that is a Fuzzy code book.
//-----------------------------------------------------------------------------

#ifndef XMIPPFUZZYCB_H
#define XMIPPFUZZYCB_H

#include "code_book.h"

/**@defgroup FuzzyCodeBook Fuzzy code book
   @ingroup ClassificationLibrary */
//@{
/**
 * This class implements an specific Data Structure for Fuzzy
 * algorithms. It inherits fron CodeBook, so this Data Structure is basically
 * a codebook with a Fuzzy-membership Matrix.
 * Some of the methods from CodeBook are redefined.
 */
class FuzzyCodeBook : public CodeBook
{
public:

    typedef std::vector< std::vector< floatFeature > > MM;
    /// Alias for Membership Matrix
    //typedef ClassificationTrainingSet<std::vector<Feature>,Label> TS;  // Alias for a Training set
    typedef ClassicTrainingVectors TS;
    /// Alias for a Training set
    typedef ClassificationTrainingSet<std::vector<floatFeature>, floatFeature > FV;
    /// Alias for Fuzzy vectors

    // Fuzzy membership matrix
    MM memb;

    /**
     * Default constructor
     * Parameter: _calib     Calibrated or not, that is, a CB with class labels or not
     */
    FuzzyCodeBook(const bool& _calib = false) : CodeBook(_calib)
    {};



    /**
     * Constructor.
     * Constructs a fuzzy codebook with initial code vectors filled with zero.
     * Parameter: _n       Number of code vectors
     * Parameter: _size    Size of code vectors
     * Parameter: _data    Size of the training set
     * Parameter: _cal     Calibrated or not, that is, a CB with class labels or not
       It calls Base Class constructor (CodeBook)
     */


    FuzzyCodeBook(unsigned _n, unsigned _size, unsigned _data, bool _cal = false);


    /**
     * Constructor.
     * Constructs a fuzzy codebook with random initial code vectors.
     * Parameter: _n       Number of code vectors
     * Parameter: _size    Size of code vectors
     * Parameter: _data    Size of the training set
     * Parameter: _lower   Lower value for random elements
     * Parameter: _upper   Upper value for random elements
     * Parameter: _cal     Calibrated or not, that is, a CB with class labels or not
       It calls Base Class constructor (CodeBook)
     */

    FuzzyCodeBook(unsigned _n, unsigned _size, unsigned _data, double _lower = 0, double _upper = 1, bool _cal = false);

    /**
     * Constructor.
     * Constructs a fuzzy codebook with initial code vectors taken randomly from
     * the training file.
     * Parameter: _n       Number of vectors
     * Parameter: _ts      Training set; will be used to get initial values
     * Parameter: _use_rand_cvs  Use random code vectors (inherited from base class)
       It calls Base Class constructor (CodeBook)
     */

    FuzzyCodeBook(unsigned _n, const ClassicTrainingVectors& _ts, const bool _use_rand_cvs = false);

    /*
     * Constructs a fuzzy code book given a stream
     * Parameter: _is  The input stream
     * Parameter: _size Size of code vectors (number of data points)
     */

    FuzzyCodeBook(std::istream& _is, const unsigned _size = 0);

    /**
     * Virtual destructor
     */

    virtual ~FuzzyCodeBook()
    {};


    /** Fuctions to access the Fuzzy Membership Matrix
    */

    /**
     * Returns a const reference to the specified item
     * Parameter: _ci  cluster index
     * Parameter: _di  data index
     * @exception out_of_range If _i is out of range
     */
    floatFeature membAt(unsigned _di, unsigned _ci) const;

    /**
     * Returns a  reference to the specified item
     * Parameter: _ci  cluster index
     * Parameter: _di  data index
     * @exception out_of_range If _i is out of range
     */
    floatFeature& membAt(unsigned _di, unsigned _ci);


    /**
     * Returns dimensions of the Membership matrix
     */
    unsigned membClusters() const;
    unsigned membVectors() const;


    /**
     * Returns the code vector that represents the input in the codebook
     * Parameter: _in  Sample to classify
       Note: The difference between Fuzzy codevector and non-Fuzzy
       codevector is that the best (winner) is estimated using the
       fuzzy membership matrix.
     */
    virtual FeatureVector& fuzzyTest(unsigned _in) const;

    /**
     * Returns the index to the code vector that represents the input in the codebook
     * Parameter: _in  Sample to classify
       Note: The difference between Fuzzy codevector and non-Fuzzy
       codevector is that the best (winner) is estimated using the
       fuzzy membership matrix.
     */
    virtual unsigned fuzzyTestIndex(unsigned _in) const;

    /**
     * Returns the label associated to an input
     * Parameter: _in  Index to the sample to be classified
     */
    virtual Label fuzzyApply(unsigned _in) const;

    /**
     * Calibrates the code book
     * Parameter: _ts   The calibrated training set
     * Parameter: _def  Default target for non-calibrated vectors
     * @exception runtime_error  If the training set is not calibrated
     */
    virtual void fuzzyCalibrate(ClassicTrainingVectors& _ts, Label _def = Label());

    /**
     * Returns the index of the codevector closest to an input.
     * This is the method used to classify inputs
     * Parameter: _in  Index to the Sample to be classified
     */
    virtual unsigned fuzzyWinner(unsigned _in) const;

    /**
    * Returns the index of the codevector closest to an input.
    * This is the method used to classify inputs
    * Parameter: _in  Index to the Sample to be classified
    */
    virtual unsigned fuzzyOutput(unsigned _in) const;

    /**
     * Fills the classifVectors with the list of the best input vectors associated to it.
     * In this case, it uses the Fuzzy Memberships to make the assignments
     * Parameter: _ts  Sample list to classify
     */
    virtual void classify(const ClassicTrainingVectors* _ts);

    /**
    * Hard partition the Fuzzy membership matrix
    * using the Nearest Maximum Membership conversion
    */

    virtual void hardPartition();


    /**
     * Returns the alpha-core set (also called alpha-level set or simply "core")
     * Parameter: _ts       The training set
     * Parameter: _alpha    A threshold to identify the core.
     * Parameter: _cluster  The cluster or partition
     * Returns a training vector "contained" in the cluster
     */
    virtual TS alphaCore(TS _ts, double _alpha, unsigned _cluster) const;

    /**
     * Writes the membership values
     * Parameter: _os  The output stream
     */
    virtual void writeMembership(std::ostream& _os) const;

    /**
     * Reads the membership values
     * Parameter: _is  The input stream
     */
    virtual void readMembership(std::istream& _is);


    /**
     * Saves the xmippFuzzyCodeBook class into a stream.
     * this method can be used to save the status of the class.
     * Parameter: _os The output stream
     */
    virtual void saveObject(std::ostream& _os) const;


    /**
     * Loads the xmippFuzzyCodeBook class from a stream.
     * this method can be used to load the status of the class.
     * Parameter: _is The output stream
     */
    virtual void loadObject(std::istream& _is);


    /**
     * Constructs a fuzzy code book given a stream
     * Parameter: _is  The input stream
     * Parameter: _size Size of code vectors (number of data points)
     * @exception  runtime_error  If there are problems with the stream
     */
    virtual void readSelf(std::istream& _is, const unsigned _size = 0);

    /**
     * Prints the density values of each Fuzzy codevector.
     * Parameter: _os  The the output stream
     */
    virtual void printDensity(std::ostream& _os) const;


private:

    // Dimensions of the membership matrix
    unsigned numClusters, numVectors;

};
//@}


//-----------------------------------------------------------------------------

#endif//XMIPPFUZZYCB_H
