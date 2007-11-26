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
// xmippMap.cc
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
xmippMap::xmippMap(const string& _layout,  unsigned _width,
                   const unsigned& _height, const unsigned& _size)
        : xmippCB(_width*_height, _size, false), somWidth(_width), somHeight(_height)
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
};




/**
 * Constructs a SOM with random initial code vectors
 * Parameter: _layout  Type of layout
 * Parameter: _width   Width of the output plane
 * Parameter: _height  Height of the output plane
 * Parameter: _size    Size of code vectors
 * Parameter: _lower   Lower value for random elements
 * Parameter: _upper   Upper value for random elements
 */
xmippMap::xmippMap(const string& _layout,  unsigned _width,
                   const unsigned& _height, const unsigned& _size, const double& _lower,
                   const double& _upper)
        : xmippCB(_width*_height, _size, _lower, _upper, false), somWidth(_width), somHeight(_height)
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
};


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
   of the department of Physiology of the David Geffen School of Medicine,
   University of California, Los Angeles
*/
xmippMap::xmippMap(const string& _layout,  unsigned _width,
                   const unsigned& _height, const xmippCTVectors& _ts,
                   const bool _use_rand_cvs)
        : xmippCB(_width*_height, _ts, _use_rand_cvs),
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
};



/**
 * Construct a SOM from the code vectors in a stream
 * Parameter: _is  The stream
 * @exception  runtime_error  If there are problems with the stream
 */
xmippMap::xmippMap(istream& _is, bool _cv) : xmippCB(false)
{
    somLayout = NULL;
    if (_cv)
        readSelf(_is);
    else
        loadObject(_is);
};


/**
 * This method throws an exception if called. There is no sense in adding
 * vectors to a som
 * @exception range_error  If this method is called
 */
void xmippMap::add(const xmippVector& _v, const xmippLabel& _l)
{
    throw runtime_error("You can't add vectors to a S.O.M.");
};

/**
 * Returns the id of layout that som has
 */
const string& xmippMap::layout() const
{
    return somLayout->id();
};

/**
 * Returns the neighborhood of a neuron
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
vector<unsigned> xmippMap::neighborhood(const SomPos& _center, double _radius) const
{
    return somLayout->neighborhood(this, _center, _radius);
};

/**
 * Returns the distance between two neurons according to the Layout
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _v       Position of the code vector
 */
double xmippMap::neighDist(const SomPos& _center, const SomPos& _v) const
{
    return somLayout->dist(_center, _v);
};


/**
 * Returns the width of the SOM
 */
unsigned xmippMap::width() const
{
    return somWidth;
};

/**
 * Returns the height of the SOM
 */
unsigned xmippMap::height() const
{
    return somHeight;
};

/**
 * Returns a code vector given its position
 * Parameter: _pos  The position of the code vector
 * @exception out_of _range   If _pos is out of range
 */
SomIn& xmippMap::itemAtPos(const SomPos& _pos)
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return itemAt((somWidth * _pos.second) + _pos.first);
};



/**
 * Returns a const reference to a code vector given its position
 * Parameter: _pos  The position of the code vector
 * @exception out_of _range   If _pos is out of range
 */
const SomIn& xmippMap::itemAtPos(const SomPos& _pos) const
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return itemAt((somWidth * _pos.second) + _pos.first);
};


/**
 * Returns the target of a code vector given its position
 * Parameter: _pos  The position of the code vector
 */
xmippLabel& xmippMap::targetAtPos(const SomPos& _pos)
{
    if (!calibrated())
    {
        ostringstream msg;
        msg << "The S.O.M. is not calibrated. No target at position ("
        << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No target at position (" << _pos.first << ", "
        << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return targetAt((somWidth * _pos.second) + _pos.first);
};


/**
 * Returns a const target of a code vector given its position
 * Parameter: _pos  The position of the code vector
 */
const xmippLabel& xmippMap::targetAtPos(const SomPos& _pos) const
{
    if (!calibrated())
    {
        ostringstream msg;
        msg << "The S.O.M. is not calibrated. No target at position ("
        << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No target at position (" << _pos.first << ", "
        << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return targetAt((somWidth * _pos.second) + _pos.first);
};



/**
 * Clears the Som
 */
void xmippMap::clear()
{
    xmippCB::clear();
    somLayout = NULL;
    if (somLayout)
        delete somLayout;
    somWidth = 0;
    somHeight = 0;
};

/**
 * Return the position associated to an index
 * Parameter: _i  Index of the code vector
 * @exception out_of _range   If _i is out of range
 */
SomPos xmippMap::indexToPos(const unsigned& _i) const
{
    if (_i >= somWidth * somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No item at position " << _i;

        throw out_of_range(msg.str());
    }

    return SomPos(_i % somWidth, _i / somWidth);
};

/**
 * Return the index associated to a position
 * Parameter: _pos  Position of the code vector
 * @exception out_of _range   If _i is out of range
 */
unsigned xmippMap::PosToIndex(const SomPos& _pos) const
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return (unsigned)((somWidth * _pos.second) + _pos.first);
};


/**
 * Returns the position in the som of a code vector
 * Parameter: _v  Reference to the code vector
 */
SomPos xmippMap::codVecPos(SomIn& _v)
{
    return indexToPos(&_v - &(itemAt(0)));
};

/**
 * Returns the position of the code vector that represents the input in the
 * som
 * Parameter: _in  Sample to classify
 */
SomPos xmippMap::applyPos(const SomIn& _in)
{
    return codVecPos(test(_in));
};


/**
 * Standard output for a SOM
 * Parameter: _os The output stream
 */
void xmippMap::printSelf(ostream& _os) const
{
    _os << itemAt(0).size() << " " <<
    somLayout->id() << " " << somWidth << " " << somHeight << " gaussian" << endl;
    writeItems(_os);
};


/**
 * Standard input for a SOM
 * Parameter: _is The input stream
 */
void xmippMap::readSelf(istream& _is)
{
    clear();
    int dim;
    string layout, str;
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
    xmippCB::readSelf(_is, dim, somWidth*somHeight);

    /*  IT DOESN'T WORK PROPERLY

         _is >> somWidth;
         _is >> somHeight;

       char aux[128];
         strstream ostr(aux,sizeof(aux));
         ostr << dim << " " << (somWidth*somHeight) << " " <<
            somWidth << " " << somHeight << endl;
         while (_is.good())
            ostr.put (_is.get());
         xmippCB::readSelf(ostr);*/
};


/**
 * Saves the xmippMap class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void xmippMap::saveObject(ostream& _os) const
{
    _os << somLayout->id() << " " << somWidth << " " << somHeight << endl;
    writeClassifVectors(_os);
    xmippCTSet<xmippVector, xmippLabel>::saveObject(_os);
};


/**
 * Loads the xmippMap class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void xmippMap::loadObject(istream& _is)
{
    clear();
    string layout;
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
    xmippCTSet<xmippVector, xmippLabel>::loadObject(_is);
};


/**
 * Constructs a neighborhood
 * Parameter: _som     The som
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
vector<unsigned> Layout::neighborhood(const xmippMap* _som, const SomPos& _center,
                                      double _radius) const
{
    vector<unsigned> neig;

    // try to find the neighbors
    for (unsigned i = 0 ; i < _som->size() ; i++)
    {
        if (isIn(_center, _som->indexToPos(i), _radius))
            neig.push_back(i);
    }

    return neig;
};


/**
 * Constructs a neighborhood
 * Parameter: _som     The Fuzzy Som
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
vector<unsigned> Layout::neighborhood(const xmippFuzzyMap* _som, const SomPos& _center,
                                      double _radius) const
{
    vector<unsigned> neig;

    // try to find the neighbors
    for (unsigned i = 0 ; i < _som->size() ; i++)
    {
        if (isIn(_center, _som->indexToPos(i), _radius))
        {
            neig.push_back(i);
        }
    }

    return neig;
};


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
};


/**
 * Returns the id of the layout
 */
const string& Layout::id() const
{
    return theId;
};


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
};

/**
 * Returns the local average of a neuron in a non-const reference.
 *     (average of the sourounding vectors)
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _aveVector: returns the average vector
 */
void RECTLayout::localAve(const xmippFuzzyMap* _som, const SomPos& _center, vector<double>& _aveVector) const
{
    int j;
    int dim = _som->itemAt(0).size();
    for (j = 0; j < dim; j++) _aveVector[j] = 0.0;
    int tmpi = _center.first;
    int tmpj = _center.second;
    int kk = 0;
    if ((tmpi - 1) >= 0)
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi - 1, tmpj)))[j]);
    }
    if ((tmpi + 1) < _som->width())
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
    if ((tmpj + 1) < _som->height())
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi, tmpj + 1)))[j]);
    }
    if (_som->height() == 1 || _som->width() == 1)
        for (j = 0; j < dim; j++)
            _aveVector[j] /= 2.0;
    else
        for (j = 0; j < dim; j++)
            _aveVector[j] /= 4.0;
}


/**
 * Returns the average number of intermediate neighbors.
 * Parameter: _center  Reference to the center of neighborhood
 */
double RECTLayout::numNeig(const xmippFuzzyMap* _som, const SomPos& _center) const
{
    int dim = _som->itemAt(0).size();
    int tmpi = _center.first;
    int tmpj = _center.second;
    double kk = 0;
    if ((tmpi - 1) >= 0)
    {
        kk++;
    }
    if ((tmpi + 1) < _som->width())
    {
        kk++;
    }
    if ((tmpj - 1) >= 0)
    {
        kk++;
    }
    if ((tmpj + 1) < _som->height())
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
};


/**
 * Returns the local average of a neuron in a non-const reference.
 *     (average of the sourounding vectors)
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _aveVector: returns the average vector
 */
void HEXALayout::localAve(const xmippFuzzyMap* _som, const SomPos& _center, vector<double>& _aveVector) const
{

    int j;
    int dim = _som->itemAt(0).size();
    for (j = 0; j < dim; j++) _aveVector[j] = 0.0;
    int tmpi = _center.first;
    int tmpj = _center.second;
    int kk = 0;
    if ((tmpi - 1) >= 0)
    {
        kk++;
        for (j = 0; j < dim; j++)
            _aveVector[j] += (double)(_som->itemAt(_som->PosToIndex(SomPos(tmpi - 1, tmpj)))[j]);
    }
    if ((tmpi + 1) < _som->width())
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
    if ((tmpj + 1) < _som->height())
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
    if (((tmpj + 1) < _som->height()) && ((tmpi - 1) >= 0))
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
double HEXALayout::numNeig(const xmippFuzzyMap* _som, const SomPos& _center) const
{
    int dim = _som->itemAt(0).size();
    int tmpi = _center.first;
    int tmpj = _center.second;
    double kk = 0;
    if ((tmpi - 1) >= 0)
    {
        kk++;
    }
    if ((tmpi + 1) < _som->width())
    {
        kk++;
    }
    if ((tmpj - 1) >= 0)
    {
        kk++;
    }
    if ((tmpj + 1) < _som->height())
    {
        kk++;
    }
    if (((tmpj - 1) >= 0) && ((tmpi - 1) >= 0))
    {
        kk++;
    }
    if (((tmpj + 1) < _som->height()) && ((tmpi - 1) >= 0))
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
xmippFuzzyMap::xmippFuzzyMap(const string& _layout,  unsigned _width,
                             const unsigned& _height, const unsigned& _size, const double& _lower,
                             const double& _upper)
        : xmippFCB(_width*_height, _size, 0, _lower, _upper, false), somWidth(_width), somHeight(_height)
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
};

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
   of the department of Physiology of the David Geffen School of Medicine,
   University of California, Los Angeles
*/
xmippFuzzyMap::xmippFuzzyMap(const string& _layout,  unsigned _width,
                             const unsigned& _height, const xmippCTVectors& _ts,
                             const bool _use_rand_cvs)
        : xmippFCB(_width*_height, _ts, _use_rand_cvs), somWidth(_width), somHeight(_height)
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
};

/**
 * Construct a SOM from the code vectors in a stream
 * Parameter: _is  The stream
 * Parameter: _size Size of code vectors (number of data points)
 * @exception  runtime_error  If there are problems with the stream
 */
xmippFuzzyMap::xmippFuzzyMap(istream& _is, const unsigned _size, bool _cv) : xmippFCB(false)
{
    somLayout = NULL;
    if (_cv)
        readSelf(_is, _size);
    else
        loadObject(_is);
};


/**
 * This method throws an exception if called. There is no sense in adding
 * vectors to a som
 * @exception range_error  If this method is called
 */
void xmippFuzzyMap::add(const xmippVector& _v, const xmippLabel& _l)
{
    throw runtime_error("You can't add vectors to a S.O.M.");
};

/**
 * Returns the id of layout that som has
 */
const string& xmippFuzzyMap::layout() const
{
    return somLayout->id();
};

/**
 * Returns the neighborhood of a neuron
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
vector<unsigned> xmippFuzzyMap::neighborhood(const SomPos& _center, double _radius) const
{
    return somLayout->neighborhood(this, _center, _radius);
};

/**
 * Returns the neighborhood of a neuron in a non-const reference.
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _radius  Radius of neighbohood
 */
void xmippFuzzyMap::neighborhood(const SomPos& _center, double _radius, vector<unsigned>& _neig) const
{
    _neig = somLayout->neighborhood(this, _center, _radius);
};

/**
 * Returns the local average of a neuron in a non-const reference.
 *     (average of the sourounding vectors)
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _aveVector: returns the average vector
 */
void xmippFuzzyMap::localAve(const SomPos& _center, vector<double>& _aveVector) const
{
    somLayout->localAve(this, _center, _aveVector);
};


/**
 * Returns the distance between two neurons according to the Layout
 * Parameter: _center  Reference to the center of neighborhood
 * Parameter: _v       Position of the code vector
 */
double xmippFuzzyMap::neighDist(const SomPos& _center, const SomPos& _v) const
{
    return somLayout->dist(_center, _v);
};


/**
 * Returns the width of the SOM
 */
unsigned xmippFuzzyMap::width() const
{
    return somWidth;
};

/**
 * Returns the height of the SOM
 */
unsigned xmippFuzzyMap::height() const
{
    return somHeight;
};

/**
 * Returns a code vector given its position
 * Parameter: _pos  The position of the code vector
 * @exception out_of _range   If _pos is out of range
 */
SomIn& xmippFuzzyMap::itemAtPos(const SomPos& _pos)
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return itemAt((somWidth * _pos.second) + _pos.first);
};


/**
 * Returns a const reference to a code vector given its position
 * Parameter: _pos  The position of the code vector
 * @exception out_of _range   If _pos is out of range
 */
const SomIn& xmippFuzzyMap::itemAtPos(const SomPos& _pos) const
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return itemAt((somWidth * _pos.second) + _pos.first);
};


/**
 * Returns the target of a code vector given its position
 * Parameter: _pos  The position of the code vector
 */
xmippLabel& xmippFuzzyMap::targetAtPos(const SomPos& _pos)
{
    if (!calibrated())
    {
        ostringstream msg;
        msg << "The S.O.M. is not calibrated. No target at position ("
        << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No target at position (" << _pos.first << ", "
        << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return targetAt((somWidth * _pos.second) + _pos.first);
};


/**
 * Returns a const target of a code vector given its position
 * Parameter: _pos  The position of the code vector
 */
const xmippLabel& xmippFuzzyMap::targetAtPos(const SomPos& _pos) const
{
    if (!calibrated())
    {
        ostringstream msg;
        msg << "The S.O.M. is not calibrated. No target at position ("
        << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No target at position (" << _pos.first << ", "
        << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return targetAt((somWidth * _pos.second) + _pos.first);
};



/**
 * Clears the Som
 */
void xmippFuzzyMap::clear()
{
    xmippFCB::clear();
    if (somLayout)
        delete somLayout;
    somWidth = 0;
    somHeight = 0;
};

/**
 * Return the position associated to an index
 * Parameter: _i  Index of the code vector
 * @exception out_of _range   If _i is out of range
 */
SomPos xmippFuzzyMap::indexToPos(const unsigned& _i) const
{
    if (_i >= somWidth * somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No item at position " << _i;

        throw out_of_range(msg.str());
    }

    return SomPos(_i % somWidth, _i / somWidth);
};


/**
 * Return the index associated to a position
 * Parameter: _pos  Position of the code vector
 * @exception out_of _range   If _i is out of range
 */
unsigned xmippFuzzyMap::PosToIndex(const SomPos& _pos) const
{
    if (_pos.first >= (signed)somWidth || _pos.second >= (signed)somHeight)
    {
        ostringstream msg;
        msg << "Out of range. No item at position (" << _pos.first << ", " << _pos.second << ")";
        throw out_of_range(msg.str());
    }

    return (unsigned)((somWidth * _pos.second) + _pos.first);
};


/**
 * Returns the position in the som of a code vector
 * Parameter: _v  Reference to the code vector
 */
SomPos xmippFuzzyMap::codVecPos(SomIn& _v)
{
    return indexToPos(&_v - &(itemAt(0)));
};

/**
 * Returns the position of the code vector that represents the input in the
 * som
 * Parameter: _in  Sample to classify (index to the sample)
 */
SomPos xmippFuzzyMap::applyPos(const unsigned& _in)
{
    return codVecPos(fuzzyTest(_in));
};


/**
 * Standard output for a Fuzzy SOM
 * Parameter: _os The output stream
 */
void xmippFuzzyMap::printSelf(ostream& _os) const
{
    _os << itemAt(0).size() << " " <<
    somLayout->id() << " " << somWidth << " " << somHeight << " gaussian" << endl;
    writeItems(_os);
};


/**
 * Standard input for a Fuzzy SOM
 * Parameter: _is The input stream
 */
void xmippFuzzyMap::readSelf(istream& _is, const unsigned _size)
{
    clear();
    int dim;
    string layout, str;
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
        if (_is) _is.putback((char) str[i]);
    xmippFCB::readSelf(_is, _size);

};


/**
 * Saves the xmippFuzzyMap class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void xmippFuzzyMap::saveObject(ostream& _os) const
{
    _os << somLayout->id() << " " << somWidth << " " << somHeight << endl;
    writeClassifVectors(_os);
    writeMembership(_os);
    xmippCTSet<xmippVector, xmippLabel>::saveObject(_os);
};


/**
 * Loads the xmippFuzzyMap class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void xmippFuzzyMap::loadObject(istream& _is)
{
    clear();
    string layout;
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
    xmippCTSet<xmippVector, xmippLabel>::loadObject(_is);
};
