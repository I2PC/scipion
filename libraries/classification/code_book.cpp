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
// CodeBook.cc
//-----------------------------------------------------------------------------


#include "code_book.h"
#include <data/xmipp_funcs.h>

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


/**
 * Constructor.
 * Constructs a codebook with initial code vectors at zero.
 * from an unsigned integer to instantiate the template
 * Parameter: _n       Number of vectors
 * Parameter: _size    Size of code vectors
 * Parameter: _lower   Lower value for random elements
 * Parameter: _upper   Upper value for random elements
 * Parameter: _cal     Calibrated or not, that is, a CB with class labels or not
 */
CodeBook::CodeBook(unsigned _n, unsigned _size, bool _cal)
        : ClassificationDataSet<FeatureVector, Label>(),
        ClassificationTrainingSet<FeatureVector, Label>(_cal)
{
    // Fill vectors with zeros
    theItems.resize(_n);
    //if (calibrated())
    theTargets.resize(_n, "");
    for (unsigned i = 0 ; i < _n ; i++)
    {
        FeatureVector v;
        v.resize(_size, 0);
        theItems[i] = v;
    }
}



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
CodeBook::CodeBook(unsigned _n, unsigned _size, floatFeature _lower,
                   floatFeature _upper, bool _cal)
        : ClassificationDataSet<FeatureVector, Label>(),
        ClassificationTrainingSet<FeatureVector, Label>(_cal)
{
    theItems.resize(_n);
    //if (calibrated())
    theTargets.resize(_n, "");
    // Assign random vectors
    for (unsigned i = 0 ; i < _n ; i++)
    {
        std::vector<floatFeature> v;
        v = randomVector(_size, _lower, _upper);
        theItems[i] = v;
    }
}

/**
 * Constructor.
 * Constructs a codebook with initial code vectors taken randomly from
 * the training file.
 * from an unsigned integer to instantiate the template
 * Parameter: _n       Number of vectors
 * Parameter: _ts      Training set; will be used to get initial values
 * Parameter: _use_rand_cvs  Use random code vector values
 */

/* Part of this code were developed by Lorenzo Zampighi and Nelson Tang
   of the department of Physiology of the David Geffen School of Medistd::cine,
   University of California, Los Angeles
*/

CodeBook::CodeBook(unsigned _n, const ClassicTrainingVectors& _ts, const bool _use_rand_cvs)
        : ClassificationDataSet<FeatureVector, Label>(),
        ClassificationTrainingSet<FeatureVector, Label>(_ts.calibrated())
{
    // Take random samples
    //    RandomUniformGenerator<unsigned> chosen( 0, _ts.size() -1 );
    randomize_random_generator();

    theItems.resize(_n);
    //if (_ts.calibrated())
    theTargets.resize(_n, "");
    // Assign random vectors
    std::vector<floatFeature> v;
    for (unsigned i = 0 ; i < _n ; i++)
    {
        int index = (int) rnd_unif(0, _ts.size() - 1);
        v = _ts.theItems[index];
        if (_use_rand_cvs)
        {
            // NT: Scan this vector for the range of pixel values
            floatFeature minval, maxval;
            std::vector<floatFeature>::const_iterator viter = v.begin();
            minval = maxval = *viter;
            for (viter++; viter != v.end(); viter++)
            {
                if (*viter < minval)
                    minval = *viter;
                else if (*viter > maxval)
                    maxval = *viter;
            }
            v = randomVector(_ts.theItems[0].size(), minval, maxval);
        }
        theItems[i] = v;
    }

}

/**
 * Constructs a code book given a stream
 * Parameter: _is  The input stream
 * @exception  runtime_error  If there are problems with the stream
 */
CodeBook::CodeBook(std::istream& _is) : ClassificationDataSet<FeatureVector, Label>(),
        ClassificationTrainingSet<FeatureVector, Label>(_is)
{
    readSelf(_is);
}

/**
 * Returns the code vector that represents the input in the codebook
 * Parameter: _in    Sample to classify
 */
FeatureVector& CodeBook::test(const FeatureVector& _in) const
{
    // eval the first one to init best & bestDist
    std::vector<FeatureVector>::const_iterator i = itemsBegin();
    std::vector<FeatureVector>::const_iterator best = i;
    double bestDist = euclideanDistance(*i, _in);

    // eval the rest
    for (i++ ; i < itemsEnd() ; i++)
    {
        double dist = euclideanDistance(*i, _in);
        if (dist < bestDist)
        {
            bestDist = dist;
            best = i;
        }
    }

    return (FeatureVector&)*best;
}

/**
 * Returns the index to the code vector that represents the input in the codebook
 * Parameter: _in  Sample to classify
 */
unsigned CodeBook::testIndex(const FeatureVector& _in) const
{
    // eval the first one to init best & bestDist
    std::vector<FeatureVector>::const_iterator i = itemsBegin();
    std::vector<FeatureVector>::const_iterator best = i;
    double bestDist = euclideanDistance(*i, _in);

    // eval the rest
    unsigned bestIndex = 0, index = 1;
    for (i++ ; i < itemsEnd() ; i++)
    {
        double dist = euclideanDistance(*i, _in);
        if (dist < bestDist)
        {
            bestDist = dist;
            best = i;
            bestIndex = index;
        }
        index++;
    }

    return bestIndex;
}

/**
 * Returns the index of the codevector closest to an input.
 * This is the method used to classify inputs
 * Parameter: _ts  Training set
 * Parameter: _in  Index to the Sample to be classified
 */
unsigned CodeBook::winner(const ClassicTrainingVectors& _ts, unsigned _in) const
{
    return testIndex(_ts.theItems[_in]);
}


/**
 * Fills the classifVectors with the list of the best input vectors associated to it.
 * Parameter: _ts  Sample list to classify
 */
void CodeBook::classify(const ClassicTrainingVectors* _ts)
{
    classifVectors.clear(); // clear previous classification.
    classifVectors.resize(size());
    aveDistances.clear(); // clear previous classification.
    aveDistances.resize(size());
    for (unsigned j = 0 ; j < _ts->size() ; j++)
        classifVectors[testIndex(_ts->theItems[j])].push_back(j);

    for (unsigned i = 0 ; i < size() ; i++)
    {
        double aveDist = 0;
        for (unsigned j = 0 ; j < classifVectors[i].size() ; j++)
            aveDist += euclideanDistance(theItems[i], _ts->theItems[classifVectors[i][j]]);
        if (classifVectors[i].size() != 0)
            aveDist /= (double) classifVectors[i].size();
        aveDistances[i] = (double) aveDist;
    }

}


/**
 * Prints the histogram values of each Fuzzy codevector.
 * Parameter: _os  The the output stream
 */
void CodeBook::printHistogram(std::ostream& _os) const
{
    _os << "1 " << size() << std::endl;
    for (size_t j = 0; j < size(); j++)
        _os << j << " " << classifSizeAt(j) << std::endl;
}


/**
 * Prints the Average Quantization Error of each codevector.
 * Parameter: _os  The the output stream
 */
void CodeBook::printQuantError(std::ostream& _os) const
{
    _os << "1 " << size() << std::endl;
    for (size_t j = 0; j < size(); j++)
        _os << j << " " << aveDistances[j] << std::endl;
}

/**
 * Returns the list of input vectors associated to this code vector.
 */
const std::vector< unsigned>& CodeBook::classifAt(const unsigned& _index) const
{
    if (_index < 0 || _index > classifVectors.size())
    {
        std::ostringstream msg;
        msg << "index out of range";
        throw std::runtime_error(msg.str());
    }
    return classifVectors[_index];
}

/**
* Returns the number of input vectors associated to this code vector.
*/
unsigned CodeBook::classifSizeAt(const unsigned& _index) const
{
    if (_index < 0 || _index > classifVectors.size())
    {
        std::ostringstream msg;
        msg << "index out of range";
        throw std::runtime_error(msg.str());
    }
    return classifVectors[_index].size();
}


/**
 * Returns the label associated to an input
 * Parameter: _in  Sample to classify
 */
Label CodeBook::apply(const FeatureVector& _in) const
{
    return theTargets[testIndex(_in)];
}


/**
 * Calibrates the code book
 * Parameter: _ts   The calibrated training set
 * Parameter: _def  Default target for non-calibrated vectors
 * @exception runtime_error  If the training set is not calibrated
 */
void CodeBook::calibrate(ClassicTrainingVectors& _ts,
                         Label _def)
{
    // set the default label
    for (std::vector<FeatureVector>::const_iterator i = itemsBegin() ;
         i < itemsEnd() ; i++)
        theTargets[i - itemsBegin()] = _def;
    if (_ts.calibrated())
    {
        for (unsigned j = 0 ; j < _ts.size() ; j++)
            theTargets[testIndex(_ts.theItems[j])] = _ts.theTargets[j];
        calibrated(true);
    }
    else
        calibrated(false);
}

/**
* Returns the index of the codevector closest to an input.
* This is the method used to classify inputs
* Parameter: _in  Sample to classify.
*/
unsigned CodeBook::output(const FeatureVector& _in)
const
{
    return testIndex(_in);
}


/**
 * Standard output for a codebook
 * Parameter: _os The output stream
 */
void CodeBook::printSelf(std::ostream& _os) const
{
    ClassificationTrainingSet<FeatureVector, Label>::printSelf(_os);
}

/**
 * Standard input for a codebook
 * Parameter: _is The input stream
 */
void CodeBook::readSelf(std::istream& _is, long _dim, long _size)
{

#ifndef _NO_EXCEPTION
    try
    {
#endif
        clear();
        std::string line;

        // Determines the number of rows and columns in the training set

        long dim, size;
        if (_dim == -1)
        {
            _is >> dim;
        }
        else
            dim = _dim;
        if (_size == -1)
        {
            _is >> line;
            if (!sscanf(line.c_str(), "%ld", &size))
            {
                int x, y;
                _is >> x;
                _is >> y;
                size = x * y;
            }
        }
        else
            size = _size;
        getline(_is, line);
        theItems.resize(size);
        theTargets.resize(size);

        for (int i = 0; i < size; i++)
        {
            std::vector<floatFeature> v;
            v.resize(dim);
            for (int j = 0; j < dim; j++)
            {
                floatFeature var;
                _is >> var;
                v[j] = var;
            }
            getline(_is, line);
            theItems[i] = v;
            theTargets[i] = line;
        }
#ifndef _NO_EXCEPTION

    }
    catch (std::exception& e)
    {
        std::ostringstream msg;
        msg << e.what() << std::endl << "Error reading the code book";
        throw std::runtime_error(msg.str());
    }
#endif
}


/**
 * Reads the classif vectors from a stream.
 * Parameter: _is  The input stream
 */
void CodeBook::readClassifVectors(std::istream& _is)
{
    int dim;
    _is >> dim;
    classifVectors.resize(dim);
    for (size_t i = 0; i < classifVectors.size(); i++)
        _is >> classifVectors[i];
}


/**
 * Writes the classif vectors to a stream
 * Parameter: _os  The output stream
 */
void CodeBook::writeClassifVectors(std::ostream& _os) const
{
    _os << classifVectors.size() << std::endl;
    for (size_t i = 0; i < classifVectors.size(); i++)
        _os << classifVectors[i] << std::endl;
}


/**
 * Saves the CodeBook class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void CodeBook::saveObject(std::ostream& _os) const
{
    writeClassifVectors(_os);
    ClassificationTrainingSet<FeatureVector, Label>::saveObject(_os);
}


/**
 * Loads the CodeBook class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void CodeBook::loadObject(std::istream& _is)
{
    clear();
    readClassifVectors(_is);
    ClassificationTrainingSet<FeatureVector, Label>::loadObject(_is);
}


/**
 * UnNormalize all features in the codebook
 *  Parameter: _varStats The normalization information
 */

void CodeBook::unNormalize(const std::vector<ClassicTrainingVectors::statsStruct>&  _varStats)
{
    using namespace std;
    if (_varStats.size() != theItems[0].size())
    {
        std::ostringstream msg;
        msg << "Normalization information does not coincide with codebook structure";
        throw std::runtime_error(msg.str());
    }
    for (unsigned it = 0; it < size(); it++)
    {
        for (unsigned i = 0; i < theItems[0].size(); i++)
        {
            if (!isnan(theItems[it][i]))
                theItems[it][i] = theItems[it][i] * _varStats[i].sd + _varStats[i].mean;
        }
    }
}


/**
 * Normalize all features in the codebook
 *  Parameter: _varStats The normalization information
 */

void CodeBook::Normalize(const std::vector<ClassicTrainingVectors::statsStruct>&  _varStats)
{
    using namespace std;
    if (_varStats.size() != theItems[0].size())
    {
        std::ostringstream msg;
        msg << "Normalization information does not coincide with codebook structure";
        throw std::runtime_error(msg.str());
    }
    for (unsigned it = 0; it < size(); it++)
    {
        for (unsigned i = 0; i < theItems[0].size(); i++)
        {
            if (!isnan(theItems[it][i]))
            {
                if (_varStats[i].sd != 0)
                    theItems[it][i] = (theItems[it][i] - _varStats[i].mean) / _varStats[i].sd;
            }
        }
    }
}


//-----------------------------------------------------------------------------
