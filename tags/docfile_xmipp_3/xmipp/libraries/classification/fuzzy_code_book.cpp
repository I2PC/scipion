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
// xmippFuzzyCodeBook.cc
// Implements a set of fuzzy code vectors, that is a Fuzzy code book.
//-----------------------------------------------------------------------------

#include "fuzzy_code_book.h"

/**
 * Constructor.
 * Constructs a codebook with initial code vectors filled with zero.
 * Parameter: _n       Number of vectors (clusters)
 * Parameter: _size    Size of code vectors
 * Parameter: _cal     Calibrated or not, that is, a CB with class labels or not
   It calls Base Class constructor (xmippCB)
 */

xmippFCB::xmippFCB(unsigned _n, unsigned _size, unsigned _data,
                   bool _cal): xmippCB(_n, _size, _cal)
{

    // Initialize Fuzzy membership Matrix

    numClusters = _n;
    numVectors = _data;
    memb.resize(numVectors);
    for (unsigned k = 0; k < numVectors; k++)
    {
        std::vector <xmippFeature> v;
        v.resize(numClusters, 0);
        memb[k] = v;
    } // for k

};


/**
 * Constructor.
 * Constructs a codebook with random initial code vectors.
 * Parameter: _n       Number of vectors
 * Parameter: _size    Size of code vectors
 * Parameter: _lower   Lower value for random elements
 * Parameter: _upper   Upper value for random elements
 * Parameter: _cal     Calibrated or not, that is, a CB with class labels or not
   It calls Base Class constructor (xmippCB)
 */


xmippFCB::xmippFCB(unsigned _n, unsigned _size, unsigned _data, double _lower, double _upper,
                   bool _cal): xmippCB(_n, _size, _lower, _upper, _cal)
{

    // Initialize Fuzzy membership Matrix

    numClusters = _n;
    numVectors = _data;
    memb.resize(numVectors);
    for (unsigned k = 0; k < numVectors; k++)
    {
        std::vector <xmippFeature> v;
        v.resize(numClusters, 0);
        memb[k] = v;
    } // for k

};

/**
 * Constructor.
 * Constructs a codebook with initial code vectors taken randomly from
 * the training file.
 * Parameter: _n       Number of vectors
 * Parameter: _ts      Training set; will be used to get initial values
 * Parameter: _use_rand_cvs  Use random code vectors (inherited from base class)
   It calls Base Class constructor (xmippCB)
 */

/* Part of this code were developed by Lorenzo Zampighi and Nelson Tang
   of the department of Physiology of the David Geffen School of Medistd::cine,
   University of California, Los Angeles
*/

xmippFCB::xmippFCB(unsigned _n, const xmippCTVectors& _ts, const bool _use_rand_cvs) : xmippCB(_n, _ts, _use_rand_cvs)
{

    // Initialize Fuzzy membership Matrix

    numClusters = _n;
    numVectors = _ts.size();
    memb.resize(numVectors);
    for (unsigned k = 0; k < numVectors; k++)
    {
        std::vector <xmippFeature> v;
        v.resize(numClusters, 0);
        memb[k] = v;
    } // for k

};


/**
 * Constructs a fuzzy code book given a stream
 * Parameter: _is  The input stream
 * Parameter: _size Size of code vectors (number of data points)
 * @exception  runtime_error  If there are problems with the stream
 */
xmippFCB::xmippFCB(std::istream& _is, const unsigned _size)
{
    readSelf(_is);
    // Initialize Fuzzy membership Matrix

    numClusters = theItems.size();
    numVectors = _size;
    memb.clear();
    memb.resize(numVectors);
    for (unsigned k = 0; k < numVectors; k++)
    {
        std::vector <xmippFeature> v;
        v.resize(numClusters, 0);
        memb[k] = v;
    } // for k
};


/** Fuctions to access the Fuzzy Membership Matrix
*/

/**
 * Returns a const reference to the specified item
 * Parameter: _ci  cluster index
 * Parameter: _di  data index
 * @exception out_of_range If _i is out of range
 */
xmippFeature xmippFCB::membAt(unsigned _di, unsigned _ci) const
{
    std::ostringstream msg;
    if ((_di >= membVectors()) || (_ci >= membClusters()))
    {
        msg << "Out of range. No item at position " << _di << "," << _ci << std::endl;
        throw std::out_of_range(msg.str());
    }

    return memb[_di][_ci];
};

/**
 * Returns a  reference to the specified item
 * Parameter: _ci  cluster index
 * Parameter: _di  data index
 * @exception out_of_range If _i is out of range
 */
xmippFeature& xmippFCB::membAt(unsigned _di, unsigned _ci)
{
    std::ostringstream msg;
    if ((_di >= membVectors()) || (_ci >= membClusters()))
    {
        msg << "Out of range. No item at position " << _di << "," << _ci << std::endl;
        throw std::out_of_range(msg.str());
    }

    return memb[_di][_ci];
};


/**
 * Returns dimensions of the Membership matrix
 */

unsigned xmippFCB::membClusters() const
{
    return numClusters;
}
unsigned xmippFCB::membVectors() const
{
    return numVectors;
}


/**
 * Returns the code vector that represents the input in the codebook
 * Parameter: _in  Sample to classify
   Note: The difference between Fuzzy codevector and non-Fuzzy
   codevector is that the best (winner) is estimated using the
   fuzzy membership matrix.
 */

xmippVector& xmippFCB::fuzzyTest(unsigned _in) const
{
    double maxMemb = 0;
    unsigned best = 0;

    for (unsigned c = 0; c < membClusters(); c++)
    {
        if (maxMemb < memb[_in][c])
        {
            maxMemb = (double) memb[_in][c];
            best = c;
        } //if
    } // for i

    return (xmippVector&) theItems[best];
};

/**
 * Returns the index to the code vector that represents the input in the codebook
 * Parameter: _in  Sample to classify
   Note: The difference between Fuzzy codevector and non-Fuzzy
   codevector is that the best (winner) is estimated using the
   fuzzy membership matrix.
 */
unsigned xmippFCB::fuzzyTestIndex(unsigned _in) const
{
    double maxMemb = 0;
    unsigned best = 0;

    for (unsigned c = 0; c < membClusters(); c++)
    {
        if (maxMemb < memb[_in][c])
        {
            maxMemb = (double) memb[_in][c];
            best = c;
        } //if
    } // for i

    return best;
};

/**
 * Returns the label associated to an input
 * Parameter: _in  Index to the sample to be classified
 */
xmippLabel xmippFCB::fuzzyApply(unsigned _in) const
{
    return theTargets[fuzzyTestIndex(_in)];
};

/**
 * Calibrates the code book
 * Parameter: _ts   The calibrated training set
 * Parameter: _def  Default target for non-calibrated vectors
 * @exception runtime_error  If the training set is not calibrated
 */
void xmippFCB::fuzzyCalibrate(xmippCTVectors& _ts, xmippLabel _def)
{
    // set the default label
    for (std::vector<xmippVector>::const_iterator i = itemsBegin();
         i < itemsEnd() ; i++)
        theTargets[i - itemsBegin()] = _def;

    if (_ts.calibrated())
    {
        for (unsigned j = 0 ; j < _ts.size() ; j++)
            theTargets[fuzzyTestIndex(j)] = _ts.theTargets[j];
        calibrated(true);
    }
    else
        calibrated(false);
};

/**
 * Returns the index of the codevector closest to an input.
 * This is the method used to classify inputs
 * Parameter: _in  Index to the Sample to be classified
 */
unsigned xmippFCB::fuzzyWinner(unsigned _in) const
{
    return fuzzyTestIndex(_in);
};

/**
* Returns the index of the codevector closest to an input.
* This is the method used to classify inputs
* Parameter: _in  Index to the Sample to be classified
*/
unsigned xmippFCB::fuzzyOutput(unsigned _in) const
{
    return fuzzyTestIndex(_in);
};


/**
 * Fills the classifVectors with the list of the best input vectors associated to it.
 * In this case, it uses the Fuzzy Memberships to make the assignments
 * Parameter: _ts  Sample list to classify
 */
void xmippFCB::classify(const xmippCTVectors* _ts)
{
    classifVectors.clear(); // clear previous classification.
    classifVectors.resize(size());
    aveDistances.clear(); // clear previous classification.
    aveDistances.resize(size());

    for (unsigned j = 0 ; j < _ts->size() ; j++)
        classifVectors[fuzzyTestIndex(j)].push_back(j);

    for (unsigned i = 0 ; i < size() ; i++)
    {
        double aveDist = 0;
        for (unsigned j = 0 ; j < classifVectors[i].size() ; j++)
            aveDist += (double) eDist(theItems[i], _ts->theItems[classifVectors[i][j]]);
        if (classifVectors[i].size() != 0)
            aveDist /= (double) classifVectors[i].size();
        aveDistances[i] = (double) aveDist;
    }

};


/**
* Hard partition the Fuzzy membership matrix
* using the Nearest Maximum Membership conversion
*/

void xmippFCB::hardPartition()
{

    for (unsigned k = 1; k < membVectors(); k++)
    {    // Number of input vectors
        double maxMemb = 0;
        unsigned maxIndex = membClusters();         // Impossible cluster index
        for (unsigned i = 0; i < membClusters(); i++)
        { // Number of clusters
            if (maxMemb < memb[k][i])
            {     // Always find the maximum
                if (maxIndex != membClusters())
                    memb[k][maxIndex] = 0.0;    // Previous maximum index set to 0
                maxMemb = memb[k][i];
                maxIndex = i;
                memb[k][i] = 1.0;  // Assume this is the maximum
            }
            else // if maxMemb
                memb[k][i] = 0.0;
        } // for i
    } // for k
}


/**
 * Returns the alpha-core set (also called alpha-level set or simply "core")
 * Parameter: _ts       The training set
 * Parameter: _alpha    A threshold to identify the core.
 * Parameter: _cluster  The cluster or partition
 */

xmippFCB::TS xmippFCB::alphaCore(TS _ts, double _alpha, unsigned _cluster) const
{
    xmippFCB::TS _alphaSet(0, _ts.calibrated());

    _alphaSet.theItems.resize(membVectors());
    if (_ts.calibrated())
        _alphaSet.theTargets.resize(membVectors());

    if ((_alpha < 0) || (_alpha > 1))
    {
        std::ostringstream msg;
        msg << "threshold must be in [0,1]";
        throw std::runtime_error(msg.str());
    }

    if ((_cluster < 0) || (_cluster >= membClusters()))
    {
        std::ostringstream msg;
        msg << "Invalid cluster number";
        throw std::runtime_error(msg.str());
    }
    for (unsigned k = 1; k < membVectors(); k++)
    {     // Number of input vectors

        double maxMemb = memb[k][_cluster];
        unsigned maxIndex = _cluster;                 // Impossible cluster index
        for (unsigned i = 0; i < membClusters(); i++)
        {   // Number of clusters
            if (i != _cluster)
                if (maxMemb < memb[k][i])
                {      // Always find the maximum
                    maxMemb = memb[k][i];
                    maxIndex = i;
                }
        } // for i

        if (maxIndex == _cluster && maxMemb >= _alpha)        // If above threshold
            _alphaSet.theItems[k] = _ts.theItems[k];
        if (_ts.calibrated())
            _alphaSet.theTargets[k] =  _ts.theTargets[k];

    } // for k

    return _alphaSet;
}


/**
 * Writes the membership values
 * Parameter: _os  The output stream
 */
void xmippFCB::writeMembership(std::ostream& _os) const
{
    _os << membVectors() << std::endl;
    for (unsigned k = 0; k < numVectors; k++)   // Number of input vectors
        _os << memb[k] << std::endl;
};


/**
 * Reads the membership values
 * Parameter: _is  The input stream
 */
void xmippFCB::readMembership(std::istream& _is)
{
    _is >> numVectors;
    memb.resize(numVectors);
    for (unsigned k = 0; k < numVectors; k++)    // Number of input vectors
        _is >> memb[k];
    numClusters = memb[0].size();
};


/**
 * Saves the xmippFCB class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void xmippFCB::saveObject(std::ostream& _os) const
{
    writeMembership(_os);
    xmippCB::saveObject(_os);
};


/**
 * Constructs a fuzzy code book given a stream
 * Parameter: _is  The input stream
 * Parameter: _size Size of code vectors (number of data points)
 * @exception  runtime_error  If there are problems with the stream
 */
void xmippFCB::readSelf(std::istream& _is, const unsigned _size)
{
    xmippCB::readSelf(_is);
    // Initialize Fuzzy membership Matrix
    numClusters = theItems.size();
    numVectors = _size;
    memb.clear();
    memb.resize(numVectors);
    for (int k = 0; k < numVectors; k++)
    {
        std::vector <xmippFeature> v;
        v.resize(numClusters, 0);
        memb[k] = v;
    } // for k
};



/**
 * Loads the xmippFCB class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void xmippFCB::loadObject(std::istream& _is)
{
    clear();
    readMembership(_is);
    xmippCB::loadObject(_is);
};


/**
 * Prints the density values of each Fuzzy codevector.
 * Parameter: _os  The the output stream
 */
void xmippFCB::printDensity(std::ostream& _os) const
{
    _os << "1" << std::endl;
    for (int j = 0; j < numClusters; j++)
    {
        double sumDens = 0;
        for (int i = 0; i < numVectors; i++)
            sumDens += (double) memb[i][j];
        _os << j << " " << sumDens << std::endl;
    }
};


//-----------------------------------------------------------------------------

