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
// xmippSOM.hh
// Implements Kohonen Self-Organizing Feature Maps
//-----------------------------------------------------------------------------

#ifndef XMIPPSOM_H
#define XMIPPSOM_H

#include "base_algorithm.h"
#include "map.h"

//---------------------------------------------------------------------------
// Class Descent
//---------------------------------------------------------------------------

/**@name Descent class*/
//@{
/**
 * This class implements a descent from an initial value to a final one in a
 * number of steps. The transition is linear by default. It's used to
 * decrease alpha and the radius of neighborhood during training of SOM.
 */

class Descent
{
public:
    /**
     * Constructor
     * @param _initial
     * @param _final
     */
    Descent(const double _initial = 1, const double _final = 0)
            : initial(_initial), final(_final)
    {};

    /**
     * Returns the function value associated a step if the transition from
     * the initial value to the final value es made in _nSteps steps
     * @param _step    The actual step
     * @param _nsteps  The number of steps to reach the final val from the
     *                 initial one
     */
    virtual double operator()(const unsigned _step, const unsigned _nSteps)
    const;


    /**
    * Standard output for a Descent class
    * @param _os The output stream
    */
    virtual void printSelf(ostream& _os) const;

    /**
    * Standard input for a Descent class
    * @param _is The input stream
    */
    virtual void readSelf(istream& _is);

    /**
     * Saves the Descent class into a stream.
     * this method can be used to save the status of the class.
     * @param _os The output stream
     */
    virtual void saveObject(ostream& _os) const;


    /**
     * Loads the Descent class from a stream.
     * this method can be used to load the status of the class.
     * @param _is The output stream
     */
    virtual void loadObject(istream& _is);


    /**
    * Standard output for a Descent class
    * @param _os   The output stream
    * @param _desc  The class
    */
    friend ostream& operator << (ostream& _os, const Descent& _desc)
    {
        _desc.printSelf(_os);
        return _os;
    };

    /**
    * Standard input for a Descent class
    * @param _is The input stream
    * @param _desc The class
    */
    friend istream& operator >> (istream& _is, Descent& _desc)
    {
        _desc.readSelf(_is);
        return _is;
    };


protected:
    double initial;
    double final;
};

//@}

//---------------------------------------------------------------------------

/**@name Kohonen Self-Organizing Feature Map (SOM)*/
//@{
/**
 * This class trains a Kohonen's Self Organizing Feature Map
 */
class xmippSOM : public xmippBaseAlgo<xmippMap>
{
public:
    /// Type of neighborhood function
    typedef enum { GAUSSIAN = 0, BUBBLE = 1} neighType;

    /**
     * Constructs the algorithm
     * @param _alpha      How is gonna decrease alpha
     * @param _radius     How is gonna decrease the radius of neighborhood
     * @param _radius     How is gonna decrease the radius of neighborhood
       @param _neighType  Type of neighborhood function
     * @param _nSteps     Number of training steps
     */
    xmippSOM(Descent& _alpha, Descent& _radius,  neighType _neighType, unsigned long _nSteps)
            : xmippBaseAlgo<xmippMap>(), somAlpha(_alpha), somRadius(_radius), somNeigh(_neighType), somNSteps(_nSteps)
    {};

    /**
     * Construct a SOM from the code vectors in a stream
     * @param _is  The stream
     */
    xmippSOM(istream& _is);


    /**
     * Virtual destructor
     */
    virtual ~xmippSOM()
    {};

    /**
     * Sets the alpha descent function
     * @param _alpha  alpha(t)
     */
    void alpha(Descent _alpha);

    /**
     * Sets the radius descent function
     * @param _radius  radius(t)
     */
    void radius(Descent _radius);

    /**
     * Sets the number of training steps
     * @param _nSteps  Number of training steps
     */
    void nSteps(const unsigned long& _nSteps);


    /**
     * Trains the SOM
     * @param _som  The som to train
     * @param _ts   The training set
     */
    virtual void train(xmippMap& _som, xmippCTVectors& _ts) const;

    /**
     * Tests the SOM
     * @param _som        The som to test
     * @param _examples   The training set of examples
     * returns the quantization error
     */
    virtual double test(const xmippMap& _som, const TS& _examples) const;


    /**
    * Clears the Algorithm
    */

    virtual void clear();

    /**
    * Standard output for a SOM algorithm
    * @param _os The output stream
    */
    virtual void printSelf(ostream& _os) const;

    /**
    * Standard input for a SOM algorithm
    * @param _is The input stream
    */
    virtual void readSelf(istream& _is);

    /**
     * Saves the xmippSOM class into a stream.
     * this method can be used to save the status of the class.
     * @param _os The output stream
     */
    virtual void saveObject(ostream& _os) const;


    /**
     * Loads the xmippSOM class from a stream.
     * this method can be used to load the status of the class.
     * @param _is The output stream
     */
    virtual void loadObject(istream& _is);

    /**
    * Standard output for a som algorithm
    * @param _os   The output stream
    * @param _som  The som algorithm to be printed
    */
    friend ostream& operator << (ostream& _os, const xmippSOM& _som)
    {
        _som.printSelf(_os);
        return _os;
    };

    /**
    * Standard input for a som algorithm
    * @param _is The input stream
    * @param _som The som algorithm
    */
    friend istream& operator >> (istream& _is, xmippSOM& _som)
    {
        _som.readSelf(_is);
        return _is;
    };


    /**
    * Operator "="
    * @param op1 xmippSOM
    */
    xmippSOM& operator= (const xmippSOM &op1)
    {
        strstreambuf bf;
        ;
        iostream _str(&bf);
        op1.printSelf(_str);
        readSelf(_str);
        return *this;
    }


protected:
    Descent somAlpha;          /// alpha(t)
    Descent somRadius;         /// radius(t)
    neighType somNeigh;       /// Neighborhood type for training (Bubble or Gaussian)
    unsigned long somNSteps;   /// number of steps

private:
    /*
      * Trains the SOM (never used)
      * @param _som  The som to train
      * @param _ts   The training set
      */
    virtual void train(xmippMap& _som, const TS& _ts) const
        {};


};

//@}


#endif
