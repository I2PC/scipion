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
// ClassificationMap.cc
// Implements Self-Organizing Maps of the type used by Kohonen algorithms.
//-----------------------------------------------------------------------------

#include "map.h"

#include <data/args.h>


/**
 * Constructs a SOM with initial code vectors filled with zero.
 * Parameter: _layout  Type of layout
 * Parameter: _width   Width of the output plane
 * Parameter: _height  Height of the output plane
 * Parameter: _size    Size of code vectors
 */
ClassificationMap::ClassificationMap(const std::string& _layout,  unsigned _width,
                                     const unsigned& _height, const unsigned& _size)
        : CodeBook(_width*_height, _size, false), somWidth(_width), somHeight(_height)
{
    if (_layout == "HEXA")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
}




/**
 * Constructs a SOM with random initial code vectors
 * Parameter: _layout  Type of layout
 * Parameter: _width   Width of the output plane
 * Parameter: _height  Height of the output plane
 * Parameter: _size    Size of code vectors
 * Parameter: _lower   Lower value for random elements
 * Parameter: _upper   Upper value for random elements
 */
ClassificationMap::ClassificationMap(const std::string& _layout,  unsigned _width,
                                     const unsigned& _height, const unsigned& _size, const double& _lower,
                                     const double& _upper)
        : CodeBook(_width*_height, _size, _lower, _upper, false), somWidth(_width), somHeight(_height)
{
    if (_layout == "HEXA")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
}


/**
 * Constructs a SOM with initial code vectors taken randomly from
 * the training file.
 * Parameter: _layout  Type of layout
 * Parameter: _width   Width of the output plane
 * Parameter: _height  Height of the output plane
 * Parameter: _ts      Training set; will be used to get initial values
 * Parameter: _use_rand_cvs  Use random code vector pixel values
 */
/* Part of this code were developed by Lorenzo Zampighi and Nelson Tang
   of the department of Physiology of the David Geffen School of Medistd::cine,
   University of California, Los Angeles
*/
ClassificationMap::ClassificationMap(const std::string& _layout,  unsigned _width,
                                     const unsigned& _height, const ClassicTrainingVectors& _ts,
                                     const bool _use_rand_cvs)
        : CodeBook(_width*_height, _ts, _use_rand_cvs),
        somWidth(_width), somHeight(_height)
{
    if (_layout == "HEXA")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
}



/**
 * Construct a SOM from the code vectors in a stream
 * Parameter: _is  The stream
 * @exception  runtime_error  If there are problems with the stream
 */
ClassificationMap::ClassificationMap(std::istream& _is, bool _cv) : CodeBook(false)
{
    somLayout = NULL;
    if (_cv)
        readSelf(_is);
    else
        loadObject(_is);
}


/**
 * This method throws an exception if called. There is no sense in adding
 * vectors to a som
 * @exception range_error  If this method is called
 */
void ClassificationMap::add(const FeatureVector& _v, const Label& _l)
{
    throw std::runtime_error("You can't add vectors to a S.O.M.");
}

/**
 * Returns the id of layout that som has
 */
const std::string& ClassificationMap::layout() const
{
    return somLayout->id();
}

/**
 * Returns the neighborhood of a neuron
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
std::vector<unsigned> ClassificationMap::neighborhood(const SomPos& _center, double _radius) const
{
    return somLayout->neighborhood(this, _center, _radius);
}

/**
 * Returns the distance between two neurons according to the Layout
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _v       Position of the code vector
 */
double ClassificationMap::neighDist(const SomPos& _center, const SomPos& _v) const
{
    return somLayout->dist(_center, _v);
}


/**
 * Returns the width of the SOM
 */
unsigned ClassificationMap::width() const
{
    return somWidth;
}

/**
 * Returns the height of the SOM
 */
unsigned ClassificationMap::height() const
{
    return somHeight;
}

/**
 * Returns a code vector given its position
 * Parameter: _pos  The position of the code vector
 * @exception out_of _range   If _pos is out of range
 */
SomIn& ClassificationMap::itemAtPos(const SomPos& _pos)
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return itemAt((somWidth * _pos.second) + _pos.first);
}



/**
 * Returns a const reference to a code vector given its position
 * Parameter: _pos  The position of the code vector
 * @exception out_of _range   If _pos is out of range
 */
const SomIn& ClassificationMap::itemAtPos(const SomPos& _pos) const
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return itemAt((somWidth * _pos.second) + _pos.first);
}


/**
 * Returns the target of a code vector given its position
 * Parameter: _pos  The position of the code vector
 */
Label& ClassificationMap::targetAtPos(const SomPos& _pos)
{
    if (!calibrated())
    {
        std::ostringstream msg;
        msg << "The S.O.M. is not calibrated. No target at position ("
        << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No target at position (" << _pos.first << ", "
        << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return targetAt((somWidth * _pos.second) + _pos.first);
}


/**
 * Returns a const target of a code vector given its position
 * Parameter: _pos  The position of the code vector
 */
const Label& ClassificationMap::targetAtPos(const SomPos& _pos) const
{
    if (!calibrated())
    {
        std::ostringstream msg;
        msg << "The S.O.M. is not calibrated. No target at position ("
        << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No target at position (" << _pos.first << ", "
        << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return targetAt((somWidth * _pos.second) + _pos.first);
}



/**
 * Clears the Som
 */
void ClassificationMap::clear()
{
    CodeBook::clear();
    somLayout = NULL;
    if (somLayout)
        delete somLayout;
    somWidth = 0;
    somHeight = 0;
}

/**
 * Return the position associated to an index
 * Parameter: _i  Index of the code vector
 * @exception out_of _range   If _i is out of range
 */
SomPos ClassificationMap::indexToPos(const unsigned& _i) const
{
    if (_i >= somWidth * somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position " << _i;

        throw std::out_of_range(msg.str());
    }

    return SomPos(_i % somWidth, _i / somWidth);
}

/**
 * Return the index associated to a position
 * Parameter: _pos  Position of the code vector
 * @exception out_of _range   If _i is out of range
 */
unsigned ClassificationMap::PosToIndex(const SomPos& _pos) const
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return (unsigned)((somWidth * _pos.second) + _pos.first);
}


/**
 * Returns the position in the som of a code vector
 * Parameter: _v  Reference to the code vector
 */
SomPos ClassificationMap::codVecPos(SomIn& _v)
{
    return indexToPos(&_v - &(itemAt(0)));
}

/**
 * Returns the position of the code vector that represents the input in the
 * som
 * Parameter: _in  Sample to classify
 */
SomPos ClassificationMap::applyPos(const SomIn& _in)
{
    return codVecPos(test(_in));
}


/**
 * Standard output for a SOM
 * Parameter: _os The output stream
 */
void ClassificationMap::printSelf(std::ostream& _os) const
{
    _os << itemAt(0).size() << " " <<
    somLayout->id() << " " << somWidth << " " << somHeight << " gaussian" << std::endl;
    writeItems(_os);
}


/**
 * Standard input for a SOM
 * Parameter: _is The input stream
 */
void ClassificationMap::readSelf(std::istream& _is)
{
    clear();
    int dim;
    std::string layout, str;
    _is >> dim;
    _is >> layout;
    toLower(layout);

    if (layout == "hexa")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
    _is >> somWidth;
    _is >> somHeight;
    /* IT DOESN'T WORK PROPERLY
    str = integerToString(dim);
    str += " ";
    str += integerToString(somWidth*somHeight);
    str += " ";
    for (int i = str.size() - 1; i >= 0; i--)
    if (_is) _is.putback((char) str[i]);
    */
    CodeBook::readSelf(_is, dim, somWidth*somHeight);

    /*  IT DOESN'T WORK PROPERLY

         _is >> somWidth;
         _is >> somHeight;

       char aux[128];
         strstream ostr(aux,sizeof(aux));
         ostr << dim << " " << (somWidth*somHeight) << " " <<
            somWidth << " " << somHeight << std::endl;
         while (_is.good())
            ostr.put (_is.get());
         CodeBook::readSelf(ostr);*/
}


/**
 * Saves the ClassificationMap class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void ClassificationMap::saveObject(std::ostream& _os) const
{
    _os << somLayout->id() << " " << somWidth << " " << somHeight << std::endl;
    writeClassifVectors(_os);
    ClassificationTrainingSet<FeatureVector, Label>::saveObject(_os);
}


/**
 * Loads the ClassificationMap class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void ClassificationMap::loadObject(std::istream& _is)
{
    clear();
    std::string layout;
    _is >> layout;
    if (layout == "HEXA")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
    _is >> somWidth;
    _is >> somHeight;
    readClassifVectors(_is);
    ClassificationTrainingSet<FeatureVector, Label>::loadObject(_is);
}


/**
 * Constructs a neighborhood
 * Parameter: _som     The som
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
std::vector<unsigned> Layout::neighborhood(const ClassificationMap* _som, const SomPos& _center,
        double _radius) const
{
    std::vector<unsigned> neig;

    // try to find the neighbors
    for (unsigned i = 0 ; i < _som->size() ; i++)
    {
        if (isIn(_center, _som->indexToPos(i), _radius))
            neig.push_back(i);
    }

    return neig;
}


/**
 * Constructs a neighborhood
 * Parameter: _som     The Fuzzy Som
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
std::vector<unsigned> Layout::neighborhood(const FuzzyMap* _som, const SomPos& _center,
        double _radius) const
{
    std::vector<unsigned> neig;

    // try to find the neighbors
    for (unsigned i = 0 ; i < _som->size() ; i++)
    {
        if (isIn(_center, _som->indexToPos(i), _radius))
        {
            neig.push_back(i);
        }
    }

    return neig;
}


/**
 * Returns true if the vector in the given position is in the neighborhood
 * or false otherwise
 * Parameter: _center  Position of the center of neighborhood
 * Parameter: _v       Position of the code vector
 * Parameter: _radius  Radius of neighbohood
 */
bool Layout::isIn(const SomPos& _center, const SomPos& _v,
                  double _radius) const
{
    return (dist(_center, _v) <= _radius);
}


/**
 * Returns the id of the layout
 */
const std::string& Layout::id() const
{
    return theId;
}


//---------------------------------------------------------------------------


/**
 * Returns the distance between two vectors in their given position
 * Parameter: _center  Position of the center of neighborhood
 * Parameter: _v       Position of the code vector
 */

double RECTLayout::dist(const SomPos& _center, const SomPos& _v) const
{
    return ((double) sqrt((double)(_center.first - _v.first)*(_center.first - _v.first) +
                          (_center.second - _v.second)*(_center.second - _v.second)));
}

/**
 * Returns the local average of a neuron in a non-const reference.
 *     (average of the sourounding vectors)
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _aveVector: returns the average vector
 */
void RECTLayout::localAve(const FuzzyMap* _som, const SomPos& _center, std::vector<double>& _aveVector) const
{
    int j;
    int dim = _som->itemAt(0).size();
    double *ptrAveVector=&(_aveVector[0]);
    memset(ptrAveVector,0,dim*sizeof(double));
    int tmpi = _center.first;
    int tmpj = _center.second;
    int kk = 0;
    if ((tmpi - 1) >= 0)
    {
        kk++;
        const floatFeature *codevector=&(_som->itemAt(_som->PosToIndex(SomPos(tmpi - 1, tmpj)))[0]);
        for (j = 0; j < dim; j++)
        	ptrAveVector[j] += codevector[j];
    }
    if ((tmpi + 1) < (int)_som->width())
    {
        kk++;
        const floatFeature *codevector=&(_som->itemAt(_som->PosToIndex(SomPos(tmpi + 1, tmpj)))[0]);
        for (j = 0; j < dim; j++)
        	ptrAveVector[j] += codevector[j];
    }
    if ((tmpj - 1) >= 0)
    {
        kk++;
        const floatFeature *codevector=&(_som->itemAt(_som->PosToIndex(SomPos(tmpi, tmpj-1)))[0]);
        for (j = 0; j < dim; j++)
        	ptrAveVector[j] += codevector[j];
    }
    if ((tmpj + 1) < (int)_som->height())
    {
        kk++;
        const floatFeature *codevector=&(_som->itemAt(_som->PosToIndex(SomPos(tmpi, tmpj+1)))[0]);
        for (j = 0; j < dim; j++)
        	ptrAveVector[j] += codevector[j];
    }
    if (_som->height() == 1 || _som->width() == 1)
        for (j = 0; j < dim; j++)
        	ptrAveVector[j] *= 1.0/2.0;
    else
        for (j = 0; j < dim; j++)
        	ptrAveVector[j] *= 1.0/4.0;
}

/**
 * Returns the average number of intermediate neighbors.
 * Parameter: _center  Reference to the center of neighborhood
 */
double RECTLayout::numNeig(const FuzzyMap* _som, const SomPos& _center) const
{
    int tmpi = _center.first;
    int tmpj = _center.second;
    double kk = 0;
    if ((tmpi - 1) >= 0)
    {
        kk++;
    }
    if ((tmpi + 1) < (int)_som->width())
    {
        kk++;
    }
    if ((tmpj - 1) >= 0)
    {
        kk++;
    }
    if ((tmpj + 1) < (int)_som->height())
    {
        kk++;
    }
    if (_som->height() == 1 || _som->width() == 1)
        return (kk / 2.0);
    else
        return (kk / 4.0);
}


/**
 * Returns the distance between two vectors in their given position
 * Parameter: _center  Position of the center of neighborhood
 * Parameter: _v       Position of the code vector
 */

double HEXALayout::dist(const SomPos& _center, const SomPos& _v) const
{
    double ret, diff;
    diff = _center.first - _v.first;

    if (((_center.second - _v.second) % 2) != 0)
    {
        if ((_center.second % 2) == 0)
        {
            diff -= 0.5;
        }
        else
        {
            diff += 0.5;
        }
    }
    ret = diff * diff;
    diff = _center.second - _v.second;
    ret += 0.75 * diff * diff;
    ret = (double) sqrt(ret);
    return(ret);
}


/**
 * Returns the local average of a neuron in a non-const reference.
 *     (average of the sourounding vectors)
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _aveVector: returns the average vector
 */
void HEXALayout::localAve(const FuzzyMap* _som, const SomPos& _center, std::vector<double>& _aveVector) const
{

    int j;
    int dim = _som->itemAt(0).size();
    for (j = 0; j < dim; j++)
        _aveVector[j] = 0.0;
    int tmpi = _center.first;
    int tmpj = _center.second;
    int kk = 0;
    if ((tmpi - 1) >= 0)
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi - 1, tmpj)))[j]);
    }
    if ((tmpi + 1) < (int)_som->width())
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi + 1, tmpj)))[j]);
    }
    if ((tmpj - 1) >= 0)
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi, tmpj - 1)))[j]);
    }
    if ((tmpj + 1) < (int)_som->height())
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi, tmpj + 1)))[j]);
    }
    if (((tmpj - 1) >= 0) && ((tmpi - 1) >= 0))
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi - 1, tmpj - 1)))[j]);
    }
    if (((tmpj + 1) < (int)_som->height()) && ((tmpi - 1) >= 0))
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi - 1, tmpj + 1)))[j]);
    }
    /*         for (j = 0; j < dim; j++)
         _aveVector[j] /= kk;*/
    if (_som->height() == 1 || _som->width() == 1)
        for (j = 0; j < dim; j++)
            _aveVector[j] /= 2.0;
    else
        for (j = 0; j < dim; j++)
            _aveVector[j] /= 6.0;
}


/**
 * Returns the average number of intermediate neighbors.
 * Parameter: _center  Reference to the center of neighborhood
 */
double HEXALayout::numNeig(const FuzzyMap* _som, const SomPos& _center) const
{
    int tmpi = _center.first;
    int tmpj = _center.second;
    double kk = 0;
    if ((tmpi - 1) >= 0)
    {
        kk++;
    }
    if ((tmpi + 1) < (int)_som->width())
    {
        kk++;
    }
    if ((tmpj - 1) >= 0)
    {
        kk++;
    }
    if ((tmpj + 1) < (int)_som->height())
    {
        kk++;
    }
    if (((tmpj - 1) >= 0) && ((tmpi - 1) >= 0))
    {
        kk++;
    }
    if (((tmpj + 1) < (int)_som->height()) && ((tmpi - 1) >= 0))
    {
        kk++;
    }
    if (_som->width() == 1 || _som->height() == 1)
        return (kk / 2.0);
    else
        return (kk / 6.0);
}


/************** Fuzzy Map ******************************************/

/**
 * Constructs a Fuzzy SOM with random initial code vectors
 * Parameter: _layout  Type of layout
 * Parameter: _width   Width of the output plane
 * Parameter: _height  Height of the output plane
 * Parameter: _size    Size of code vectors
 * Parameter: _lower   Lower value for random elements
 * Parameter: _upper   Upper value for random elements
 */
FuzzyMap::FuzzyMap(const std::string& _layout,  unsigned _width,
                   const unsigned& _height, const unsigned& _size, const double& _lower,
                   const double& _upper)
        : FuzzyCodeBook(_width*_height, _size, 0, _lower, _upper, false), somWidth(_width), somHeight(_height)
{
    if (_layout == "HEXA")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
}

/**
 * Constructs a Fuzzy SOM with initial code vectors taken randomly from
 * the training file.
 * Parameter: _layout  Type of layout
 * Parameter: _width   Width of the output plane
 * Parameter: _height  Height of the output plane
 * Parameter: _ts      Training set; will be used to get initial values
 * Parameter: _use_rand_cvs  Use random code vector pixel values
 */
/* Part of this code were developed by Lorenzo Zampighi and Nelson Tang
   of the department of Physiology of the David Geffen School of Medistd::cine,
   University of California, Los Angeles
*/
FuzzyMap::FuzzyMap(const std::string& _layout,  unsigned _width,
                   const unsigned& _height, const ClassicTrainingVectors& _ts,
                   const bool _use_rand_cvs)
        : FuzzyCodeBook(_width*_height, _ts, _use_rand_cvs),
        somWidth(_width), somHeight(_height)
{
    if (_layout == "HEXA")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
}

/**
 * Construct a SOM from the code vectors in a stream
 * Parameter: _is  The stream
 * Parameter: _size Size of code vectors (number of data points)
 * @exception  runtime_error  If there are problems with the stream
 */
FuzzyMap::FuzzyMap(std::istream& _is, const unsigned _size, bool _cv) : FuzzyCodeBook(false)
{
    somLayout = NULL;
    if (_cv)
        readSelf(_is, _size);
    else
        loadObject(_is);
}


/**
 * This method throws an exception if called. There is no sense in adding
 * vectors to a som
 * @exception range_error  If this method is called
 */
void FuzzyMap::add(const FeatureVector& _v, const Label& _l)
{
    throw std::runtime_error("You can't add vectors to a S.O.M.");
}

/**
 * Returns the id of layout that som has
 */
const std::string& FuzzyMap::layout() const
{
    return somLayout->id();
}

/**
 * Returns the neighborhood of a neuron
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
std::vector<unsigned> FuzzyMap::neighborhood(const SomPos& _center, double _radius) const
{
    return somLayout->neighborhood(this, _center, _radius);
}

/**
 * Returns the neighborhood of a neuron in a non-const reference.
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
void FuzzyMap::neighborhood(const SomPos& _center, double _radius, std::vector<unsigned>& _neig) const
{
    _neig = somLayout->neighborhood(this, _center, _radius);
}

/**
 * Returns the local average of a neuron in a non-const reference.
 *     (average of the sourounding vectors)
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _aveVector: returns the average vector
 */
void FuzzyMap::localAve(const SomPos& _center, std::vector<double>& _aveVector) const
{
    somLayout->localAve(this, _center, _aveVector);
}


/**
 * Returns the distance between two neurons according to the Layout
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _v       Position of the code vector
 */
double FuzzyMap::neighDist(const SomPos& _center, const SomPos& _v) const
{
    return somLayout->dist(_center, _v);
}


/**
 * Returns the width of the SOM
 */
unsigned FuzzyMap::width() const
{
    return somWidth;
}

/**
 * Returns the height of the SOM
 */
unsigned FuzzyMap::height() const
{
    return somHeight;
}

/**
 * Returns a code vector given its position
 * Parameter: _pos  The position of the code vector
 * @exception out_of _range   If _pos is out of range
 */
SomIn& FuzzyMap::itemAtPos(const SomPos& _pos)
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return itemAt((somWidth * _pos.second) + _pos.first);
}


/**
 * Returns a const reference to a code vector given its position
 * Parameter: _pos  The position of the code vector
 * @exception out_of _range   If _pos is out of range
 */
const SomIn& FuzzyMap::itemAtPos(const SomPos& _pos) const
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return itemAt((somWidth * _pos.second) + _pos.first);
}


/**
 * Returns the target of a code vector given its position
 * Parameter: _pos  The position of the code vector
 */
Label& FuzzyMap::targetAtPos(const SomPos& _pos)
{
    if (!calibrated())
    {
        std::ostringstream msg;
        msg << "The S.O.M. is not calibrated. No target at position ("
        << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No target at position (" << _pos.first << ", "
        << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return targetAt((somWidth * _pos.second) + _pos.first);
}


/**
 * Returns a const target of a code vector given its position
 * Parameter: _pos  The position of the code vector
 */
const Label& FuzzyMap::targetAtPos(const SomPos& _pos) const
{
    if (!calibrated())
    {
        std::ostringstream msg;
        msg << "The S.O.M. is not calibrated. No target at position ("
        << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No target at position (" << _pos.first << ", "
        << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return targetAt((somWidth * _pos.second) + _pos.first);
}



/**
 * Clears the Som
 */
void FuzzyMap::clear()
{
    FuzzyCodeBook::clear();
    if (somLayout)
        delete somLayout;
    somWidth = 0;
    somHeight = 0;
}

/**
 * Return the position associated to an index
 * Parameter: _i  Index of the code vector
 * @exception out_of _range   If _i is out of range
 */
SomPos FuzzyMap::indexToPos(const unsigned& _i) const
{
    if (_i >= somWidth * somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position " << _i;

        throw std::out_of_range(msg.str());
    }

    return SomPos(_i % somWidth, _i / somWidth);
}


/**
 * Return the index associated to a position
 * Parameter: _pos  Position of the code vector
 * @exception out_of _range   If _i is out of range
 */
unsigned FuzzyMap::PosToIndex(const SomPos& _pos) const
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw std::out_of_range(msg.str());
    }

    return (unsigned)((somWidth * _pos.second) + _pos.first);
}


/**
 * Returns the position in the som of a code vector
 * Parameter: _v  Reference to the code vector
 */
SomPos FuzzyMap::codVecPos(SomIn& _v)
{
    return indexToPos(&_v - &(itemAt(0)));
}

/**
 * Returns the position of the code vector that represents the input in the
 * som
 * Parameter: _in  Sample to classify (index to the sample)
 */
SomPos FuzzyMap::applyPos(const unsigned& _in)
{
    return codVecPos(fuzzyTest(_in));
}


/**
 * Standard output for a Fuzzy SOM
 * Parameter: _os The output stream
 */
void FuzzyMap::printSelf(std::ostream& _os) const
{
    _os << itemAt(0).size() << " " <<
    somLayout->id() << " " << somWidth << " " << somHeight << " gaussian" << std::endl;
    writeItems(_os);
}


/**
 * Standard input for a Fuzzy SOM
 * Parameter: _is The input stream
 */
void FuzzyMap::readSelf(std::istream& _is, const unsigned _size)
{
    clear();
    int dim;
    std::string layout, str;
    _is >> dim;
    _is >> layout;
    if (layout == "HEXA")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
    _is >> somWidth;
    _is >> somHeight;
    str = integerToString(dim);
    str += " ";
    str += integerToString(somWidth * somHeight);
    str += " ";
    for (int i = str.size() - 1; i >= 0; i--)
        if (_is)
            _is.putback((char) str[i]);
    FuzzyCodeBook::readSelf(_is, _size);

}


/**
 * Saves the FuzzyMap class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void FuzzyMap::saveObject(std::ostream& _os) const
{
    _os << somLayout->id() << " " << somWidth << " " << somHeight << std::endl;
    writeClassifVectors(_os);
    writeMembership(_os);
    ClassificationTrainingSet<FeatureVector, Label>::saveObject(_os);
}


/**
 * Loads the FuzzyMap class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void FuzzyMap::loadObject(std::istream& _is)
{
    clear();
    std::string layout;
    _is >> layout;
    if (layout == "HEXA")
    {
        HEXALayout *tmpLayout = new HEXALayout();
        somLayout = tmpLayout;
    }
    else
    {
        RECTLayout *tmpLayout = new RECTLayout();
        somLayout = tmpLayout;
    }
    _is >> somWidth;
    _is >> somHeight;
    readClassifVectors(_is);
    readMembership(_is);
    ClassificationTrainingSet<FeatureVector, Label>::loadObject(_is);
}
