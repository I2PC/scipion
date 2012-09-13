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
// xmippMap.hh
// Implements a Map of the type used by Kohonen algorithms.
//-----------------------------------------------------------------------------

#ifndef XMIPPMAP_H
#define XMIPPMAP_H

// STL includes
#include <algorithm>
#include <vector>
#include <strstream>

// xmipp includes
#include "code_book.h"
#include "fuzzy_code_book.h"
#include "vector_ops.h"

typedef xmippVector SomIn;
typedef pair<long, long> SomPos;

// Forward declarations

class Layout;

/**@defgroup SomMap Self Organizing Map
   @ingroup ClassificationLibrary */
//@{
/**
 * This class implements a Map of the type used by Kohonen Algorithms
 */
class xmippMap : public xmippCB
{
public:


    //---------------------------------------------------------------------------



    /**
     * Constructs a SOM with initial code vectors filled with zero.
     * Parameter: _layout  Type of layout
     * Parameter: _width   Width of the output plane
     * Parameter: _height  Height of the output plane
     * Parameter: _size    Size of code vectors
     */
    xmippMap(const string& _layout,  unsigned _width,
             const unsigned& _height, const unsigned& _size);


    /**
     * Constructs a SOM with random initial code vectors
     * Parameter: _layout  Type of layout
     * Parameter: _width   Width of the output plane
     * Parameter: _height  Height of the output plane
     * Parameter: _size    Size of code vectors
     * Parameter: _lower   Lower value for random elements
     * Parameter: _upper   Upper value for random elements
     */
    xmippMap(const string& _layout,  unsigned _width,
             const unsigned& _height, const unsigned& _size, const double& _lower,
             const double& _upper);


    /**
     * Constructs a SOM with initial code vectors taken randomly from
     * the training file.
     * Parameter: _layout  Type of layout
     * Parameter: _width   Width of the output plane
     * Parameter: _height  Height of the output plane
     * Parameter: _ts      Training set; will be used to get initial values
     * Parameter: _use_rand_cvs  Use random code vector values
     */
    xmippMap(const string& _layout,  unsigned _width,
             const unsigned& _height, const xmippCTVectors& _ts,
             const bool _use_rand_cvs = false);

    /**
     * Construct a SOM from the code vectors in a stream
     * Parameter: _is  The stream
     * Parameter: _cv  If the stream holds a codevector file or a whole codebook file
     * @exception  runtime_error  If there are problems with the stream
     */
    xmippMap(istream& _is, bool _cv = true);

    /**
     * Virtual destructor is needed
     */
    virtual ~xmippMap()
    {};

    /**
     * This method throws an exception if called. There is no sense in adding
     * vectors to a som
     * @exception range_error  If this method is called
     */
    virtual void add(const xmippVector& _v, const xmippLabel& _l = xmippLabel());

    /**
     * Returns the id of layout that som has
     */
    const string& layout() const;

    /**
     * Returns the neighborhood of a neuron
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _radius  Radius of neighbohood
     */
    vector<unsigned> neighborhood(const SomPos& _center, double _radius) const;

    /**
     * Returns the distance between two neurons according to the Layout
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _v       Position of the code vector
     */
    double neighDist(const SomPos& _center, const SomPos& _v) const;


    /**
     * Returns the width of the SOM
     */
    unsigned width() const;

    /**
     * Returns the height of the SOM
     */
    unsigned height() const;

    /**
     * Returns a code vector given its position
     * Parameter: _pos  The position of the code vector
     * @exception out_of _range   If _pos is out of range
     */
    SomIn& itemAtPos(const SomPos& _pos);

    /**
     * Returns a const reference to a code vector given its position
     * Parameter: _pos  The position of the code vector
     * @exception out_of _range   If _pos is out of range
     */
    const SomIn& itemAtPos(const SomPos& _pos) const;

    /**
     * Returns the target of a code vector given its position
     * Parameter: _pos  The position of the code vector
     */
    xmippLabel& targetAtPos(const SomPos& _pos);

    /**
     * Returns a const target of a code vector given its position
     * Parameter: _pos  The position of the code vector
     */
    const xmippLabel& targetAtPos(const SomPos& _pos) const;


    /**
     * Clears the Som
     */
    void clear();

    /**
     * Return the position associated to an index
     * Parameter: _i  Index of the code vector
     * @exception out_of _range   If _i is out of range
     */
    SomPos indexToPos(const unsigned& _i) const;


    /**
     * Return the index associated to a position
     * Parameter: _pos  Position of the code vector
     * @exception out_of _range   If _i is out of range
     */

    unsigned PosToIndex(const SomPos& _pos) const;

    /**
     * Returns the position in the som of a code vector
     * Parameter: _v  Reference to the code vector
     */
    SomPos codVecPos(SomIn& _v);

    /**
     * Returns the position of the code vector that represents the input in the
     * som
     * Parameter: _in  Sample to classify
     */
    SomPos applyPos(const SomIn& _in);

    /*
     * Standard output for a som
     * Parameter: _os   The output stream
     * Parameter: _som  The som to be printed
     */
    friend ostream& operator << (ostream& _os, const xmippMap& _som)
    {
        _som.printSelf(_os);
    };

    /*
     * Standard input for a som
     * Parameter: _is The input stream
     * Parameter: _som The code book to be read
     * @exception  runtime_error  If there are problems with the stream
     */
    friend istream& operator >> (istream& _is, xmippMap& _som)
    {
        _som.readSelf(_is);
    };


    /**
    * Operator "="
    * Parameter: op1 xmippMap
    */
    xmippMap& operator= (const xmippMap &op1)
    {
        strstreambuf bf;
        ;
        iostream _str(&bf);
        op1.printSelf(_str);
        readSelf(_str);
        return *this;
    }


    /**
     * Standard output for a map
     * Parameter: _os The output stream
     */
    virtual void printSelf(ostream& _os) const;


    /**
     * Standard input for a map
     * Parameter: _is The input stream
     */
    virtual void readSelf(istream& _is);

    /**
     * Saves the xmippMap class into a stream.
     * this method can be used to save the status of the class.
     * Parameter: _os The output stream
     */
    virtual void saveObject(ostream& _os) const;


    /**
     * Loads the xmippMap class from a stream.
     * this method can be used to load the status of the class.
     * Parameter: _is The output stream
     */
    virtual void loadObject(istream& _is);


    /**
     * Const Reference to the Map Layout
     */
    virtual const Layout& getLayout()  const
    {
        return (Layout&) *somLayout;
    };


protected:

    const Layout* somLayout;   // type of layout
    unsigned somWidth;         // SOM's width
    unsigned somHeight;        // SOM's height

};
//@}


/**@defgroup FuzySOM Fuzzy Self Organizing Map
   @ingroup ClassificationLibrary */
//@{
/**
 * This class implements a Fuzzy Map of the type used by Fuzzy mapping Algorithms
 */
class xmippFuzzyMap : public xmippFCB
{
public:


    //---------------------------------------------------------------------------

    /**
     * Constructs a Fuzzy SOM with random initial code vectors
     * Parameter: _layout  Type of layout
     * Parameter: _width   Width of the output plane
     * Parameter: _height  Height of the output plane
     * Parameter: _size    Size of code vectors
     * Parameter: _lower   Lower value for random elements
     * Parameter: _upper   Upper value for random elements
     */
    xmippFuzzyMap(const string& _layout,  unsigned _width,
                  const unsigned& _height, const unsigned& _size, const double& _lower,
                  const double& _upper);


    /**
     * Constructs a Fuzzy SOM with initial code vectors taken randomly from
     * the training file.
     * Parameter: _layout  Type of layout
     * Parameter: _width   Width of the output plane
     * Parameter: _height  Height of the output plane
     * Parameter: _ts      Training set; will be used to get initial values
     * Parameter: _use_rand_cvs  Use random code vector pixel values
     */
    xmippFuzzyMap(const string& _layout,  unsigned _width,
                  const unsigned& _height, const xmippCTVectors& _ts,
                  const bool _use_rand_cvs = false);


    /**
     * Construct a Fuzzy SOM from the code vectors in a stream
     * Parameter: _is   The stream
     * Parameter: _size Size of code vectors (number of data points)
     * Parameter: _cv   If the stream holds a codevector file or a whole codebook file
     * @exception   runtime_error  If there are problems with the stream
     */
    xmippFuzzyMap(istream& _is, const unsigned _size = 0, bool _cv = true);

    /**
     * Virtual destructor is needed
     */
    virtual ~xmippFuzzyMap()
    {};

    /**
     * This method throws an exception if called. There is no sense in adding
     * vectors to a som
     * @exception range_error  If this method is called
     */
    virtual void add(const xmippVector& _v, const xmippLabel& _l = xmippLabel());

    /**
     * Returns the id of layout that som has
     */
    const string& layout() const;

    /**
     * Returns the neighborhood of a neuron
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _radius  Radius of neighbohood
     */
    vector<unsigned> neighborhood(const SomPos& _center, double _radius) const;

    /**
     * Returns the neighborhood of a neuron in a non-const reference.
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _radius  Radius of neighbohood
     */
    void neighborhood(const SomPos& _center, double _radius, vector<unsigned>& _neig) const;


    /**
     * Returns the local average of a neuron in a non-const reference.
     *     (average of the sourounding vectors)
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _aveVector: returns the average vector
     */
    void localAve(const SomPos& _center, vector<double>& _aveVector) const;


    /**
     * Returns the distance between two neurons according to the Layout
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _v       Position of the code vector
     */
    double neighDist(const SomPos& _center, const SomPos& _v) const;


    /**
     * Returns the width of the SOM
     */
    unsigned width() const;

    /**
     * Returns the height of the SOM
     */
    unsigned height() const;

    /**
     * Returns a code vector given its position
     * Parameter: _pos  The position of the code vector
     * @exception out_of _range   If _pos is out of range
     */
    SomIn& itemAtPos(const SomPos& _pos);

    /**
     * Returns a const reference to a code vector given its position
     * Parameter: _pos  The position of the code vector
     * @exception out_of _range   If _pos is out of range
     */
    const SomIn& itemAtPos(const SomPos& _pos) const;


    /**
     * Returns the target of a code vector given its position
     * Parameter: _pos  The position of the code vector
     */
    xmippLabel& targetAtPos(const SomPos& _pos);

    /**
     * Returns a const target of a code vector given its position
     * Parameter: _pos  The position of the code vector
     */
    const xmippLabel& targetAtPos(const SomPos& _pos) const;


    /**
     * Clears the Fuzzy Som
     */
    void clear();

    /**
     * Return the position associated to an index
     * Parameter: _i  Index of the code vector
     * @exception out_of _range   If _i is out of range
     */
    SomPos indexToPos(const unsigned& _i) const;

    /**
     * Return the index associated to a position
     * Parameter: _pos  Position of the code vector
     * @exception out_of _range   If _i is out of range
     */

    unsigned PosToIndex(const SomPos& _pos) const;

    /**
     * Returns the position in the som of a code vector
     * Parameter: _v  Reference to the code vector
     */
    SomPos codVecPos(SomIn& _v);

    /**
     * Returns the position of the code vector that represents the input in the
     * som
     * Parameter: _in  Sample to classify (index to the sample)
     */
    SomPos applyPos(const unsigned& _in);

    /*
     * Standard output for a fuzzy som
     * Parameter: _os   The output stream
     * Parameter: _som  The som to be printed
     */
    friend ostream& operator << (ostream& _os, const xmippFuzzyMap& _fsom)
    {
        _fsom.printSelf(_os);
    };

    /*
     * Standard input for a som
     * Parameter: _is The input stream
     * Parameter: _som The code book to be read
     * @exception  runtime_error  If there are problems with the stream
     */
    friend istream& operator >> (istream& _is, xmippFuzzyMap& _fsom)
    {
        _fsom.readSelf(_is);
    };


    /**
    * Operator "="
    * Parameter: op1 xmippFuzzyMap
    */
    xmippFuzzyMap& operator= (const xmippFuzzyMap &op1)
    {
        strstreambuf bf;
        ;
        iostream _str(&bf);
        op1.printSelf(_str);
        readSelf(_str);
        return *this;
    }


    /**
     * Standard output for a fuzzy map
     * Parameter: _os The output stream
     */
    virtual void printSelf(ostream& _os) const;


    /**
     * Standard input for a fuzzy map
     * Parameter: _size Size of code vectors (number of data points)
     * Parameter: _is The input stream
     */
    virtual void readSelf(istream& _is, const unsigned _size = 0);

    /**
     * Saves the xmippFuzzyMap class into a stream.
     * this method can be used to save the status of the class.
     * Parameter: _os The output stream
     */
    virtual void saveObject(ostream& _os) const;


    /**
     * Loads the xmippFuzzyMap class from a stream.
     * this method can be used to load the status of the class.
     * Parameter: _is The output stream
     */
    virtual void loadObject(istream& _is);


    /**
     * Const Reference to the Map Layout
     */
    virtual const Layout& getLayout()  const
    {
        return (Layout&) *somLayout;
    };


protected:

    const Layout* somLayout;   // type of layout
    unsigned somWidth;         // SOM's width
    unsigned somHeight;        // SOM's height

};
//@}



//-----------------------------------------------------------------------------
// Different layout managers
//-----------------------------------------------------------------------------



//---------------------------------------------------------------------------
// Class to manage different Layouts
//---------------------------------------------------------------------------

/**@defgroup SOMLayouts SOM Layouts
   @ingroup ClassificationLibrary */
//@{
/**
 * This class creates a SOM Layout. The neighborhood of a point depends
 * on the type of the layout chosen.
 */
class Layout
{
public:

    /**
     * Generic constructor
     * Parameter:  _id  The identification
     */
    Layout(const string& _id = "") : theId(_id)
    {};

    /**
     * Virtual destructor
     */
    virtual ~Layout()
    {};


    /**
     * Constructs a neighborhood for the SOM
     * Parameter: _som     The som
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _radius  Radius of neighbohood
     */
    vector<unsigned> neighborhood(const xmippMap* _som, const SomPos& _center,
                                  double _radius) const;


    /**
     * Constructs a neighborhood for the Fuzzy SOM
     * Parameter: _som     The Fuzzy Som
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _radius  Radius of neighbohood
     */
    vector<unsigned> neighborhood(const xmippFuzzyMap* _som, const SomPos& _center,
                                  double _radius) const;



    /**
     * Returns the local average of a neuron in a non-const reference.
     * (average of the sourounding vectors) in a Fuzzy SOM
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _aveVector: returns the average vector
     */
    virtual void localAve(const xmippFuzzyMap* _som, const SomPos& _center, vector<double>& _aveVector) const = 0;


    /**
     * Returns the average number of neighbors in a Fuzzy SOM
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _radius: Radius
     */
    virtual double numNeig(const xmippFuzzyMap* _som, const SomPos& _center) const = 0;


    /**
     * Returns true if the vector in the given position is in the neighborhood
     * or false otherwise
     * Parameter: _center  Position of the center of neighborhood
     * Parameter: _v       Position of the code vector
     * Parameter: _radius  Radius of neighbohood
     */
    virtual bool isIn(const SomPos& _center, const SomPos& _v,
                      double _radius) const;


    /**
     * Returns the distance between two vectors in their given position
     * Parameter: _center  Position of the center of neighborhood
     * Parameter: _v       Position of the code vector
     */

    virtual double dist(const SomPos& _center, const SomPos& _v) const = 0;


    /**
     * Returns the id of the layout
     */
    const string& id() const;

private:

    string theId;
};

/**
 * This class manages a layout in this way:
 * @code
 *
 *     RECTANGULAR:
 *  O O O O O O O O O
 *  O O O O & O O O O
 *  O O O & @ & O O O
 *  O O & @ + @ & O O
 *  O O O & @ & O O O
 *  O O O O & O O O O
 *  O O O O O O O O O
 * @endcode
 */
class RECTLayout : public Layout
{
public:

    /**
     * Generic constructor
     */
    RECTLayout() : Layout("rect")
    {};


    /**
     * Returns the distance between two vectors in their given position
     * Parameter: _center  Position of the center of neighborhood
     * Parameter: _v       Position of the code vector
     */
    virtual double dist(const SomPos& _center, const SomPos& _v) const;

    /**
     * Returns the local average of a neuron in a non-const reference.
     *     (average of the sourounding vectors)
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _aveVector: returns the average vector
     */
    virtual void localAve(const xmippFuzzyMap* _som, const SomPos& _center, vector<double>& _aveVector) const;

    /**
     * Returns the average number of neighbors.
     * Parameter: _center  Reference to the center of neighborhood
     */
    virtual double numNeig(const xmippFuzzyMap* _som, const SomPos& _center) const;


protected:

    RECTLayout(const string& _id) : Layout(_id)
    {};
};

/**
 * This class manages a layout in this way:
 * (Neurons are in the center of the hexagons, so distance to the 6
 * immediate neighbors are 1)
 * \\(Xdim is ------>)
 * @code
 *
 *               HEXAGONAL:
 *          O O O O O O O O O
 *         O O O & & & O O O
 *          O O & @ @ & O O O
 *  O O & @ + @ & O O
 *   O O & @ @ & O O O
 *  O O O & & & O O O
 *   O O O O O O O O O
 * @endcode
 */
class HEXALayout : public Layout
{
public:

    /**
     * Generic constructor
     */
    HEXALayout() : Layout("hexa")
    {};

    virtual double dist(const SomPos& _center, const SomPos& _v) const;

    /**
     * Returns the local average of a neuron in a non-const reference.
     *     (average of the sourounding vectors)
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _aveVector: returns the average vector
     */
    virtual void localAve(const xmippFuzzyMap* _som, const SomPos& _center, vector<double>& _aveVector) const;

    /**
     * Returns the average number of neighbors.
     * Parameter: _center  Reference to the center of neighborhood
     * Parameter: _radius: Radius
     */
    virtual double numNeig(const xmippFuzzyMap* _som, const SomPos& _center) const;

protected:

    HEXALayout(const string& _id) : Layout(_id)
    {};
};

//@}

#endif
