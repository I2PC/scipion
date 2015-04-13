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
// ClassificationTrainingSet.hh
// Xmipp Classification Training Sets
//-----------------------------------------------------------------------------

#ifndef XMIPPTS_H
#define XMIPPTS_H

#include <cctype>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>

#include "data_types.h"
#include "vector_ops.h"

#include <data/args.h>

/**@defgroup TrainingSets Training Sets
   @ingroup ClassificationLibrary */
//@{
/**
 * Generic training sets.
 * This is the parent class for all possible training sets that can be used by
 * the classification algorithms. It provides the library with the basic
 * functionality of training vectors.
 * This a template class that takes as input the type of the training vectors
 * and the type of the target variable (for supervised analysis)
 * \\The following example shows the declaration of a training set
 * formed by a vector of doubles and a target variable of type string.
 * \\
 * \\ClassificationTrainingSet<double, string> myTS;
 */
template<class Item, class Target>
class ClassificationTrainingSet
{
public:

    /// Training sets mode
    typedef std::multimap<unsigned, unsigned, std::less<unsigned> > splitTS;

    /// iterator
    typedef splitTS::iterator splitIt;

    ///  Ways the training set can be used
    typedef enum { CROSSVAL, SPLITSAMPLE, JACKKNIFE, BOOTSTRAP } splitMode;

    /// use of samples
    typedef enum { TRAIN = 0, VALIDATION = 1, TEST = 2 } useMode;

protected:
    bool isCalibrated;
    splitTS splitTrainingSet;
    /*  std::vector<Item> theItems;
      std::vector<Target> theTargets;*/
    unsigned nTargets;  // Number of targets in the training set
public:
    std::vector<Item> theItems;
    std::vector<Target> theTargets;

    /**
     * Default constructor
     * Parameter: _calib: True if the trainign set is calibrated (labeled)
     * Parameter: _n: Number of items that the trained set should contain
     */
    ClassificationTrainingSet(const bool& _calib = true, unsigned _n = 0)
            :isCalibrated(_calib), /* splitTrainingSet(), */nTargets(0), theItems(_n) /*theTargets(_n), */
    {};

    /**
     * Constructs a training set given a stream
     * Parameter: _is  The input stream
     */

    ClassificationTrainingSet(std::istream & _is) :  isCalibrated(false), splitTrainingSet(), theItems(), theTargets()
    {
        loadObject(_is);
        computeNumTargets();
    };

    /**
     * Virtual destructor
     */
    virtual ~ClassificationTrainingSet()
    {};

    /** Sets the proportion of training, validation and testing and
        creates the map. Used for split-sample training and
        validation. Parameters must be between 0 and 1 and must add to
        less than 1; the rest will be used for testing
        Parameter: _tp proportion over 1 of samples devoted to training
        Parameter: _vp proportion over 1 of samples devoted to validation
     */
    void setSplit(float _tp, float _vp)
    {
        if ((_tp > 1) || (_tp <= 0) || (_vp > 1) || (_vp <= 0) || (_tp + _vp > 1))
        {
            throw std::invalid_argument("Split mode proportions must be < 1");
        }
        std::vector<float> acc(2); // 3 modes
        acc[TRAIN] = _tp;
        acc[VALIDATION] = _tp + _vp;
        RandomUniformGenerator<float> p(0.0, 1.0);
        for (unsigned i = 0; i < size(); i ++)
        {
            float aRnd = p();
            unsigned rw = TEST;
            for (unsigned j = TRAIN; j < TEST; j ++)
            {
                if (aRnd < acc[j])
                {
                    rw = j;
                    break;
                }
            }
            splitTrainingSet.insert(std::pair<unsigned, unsigned>(rw, i));
        }
    }



    /** Returns an iterator to the begining of the subset
        Parameter: _um is the mode: training, testing or validation
        @return an iterator to the begining of the subset;
    */
    splitIt beginSubset(unsigned _um)
    {
        return splitTrainingSet.lower_bound(_um);
    }

    /// Returns an iterator to the end of the subset
    splitIt  endSubset(unsigned _um)
    {
        return splitTrainingSet.upper_bound(_um);
    }


    /**
     * Adds an item to the training set
     * Parameter: _i   The item
     * Parameter: _tg  The target (Used only if the TS is calibrated).
     */
    virtual void add(const Item& _i, const Target& _tg/* = Target()*/)
    {
        theItems.push_back(_i);
        theTargets.push_back(_tg);
    };



    /**
     * Adds an item to a non calibrated training set
     * Parameter: _i   The item
     */
    virtual void add(const Item& _i)
    {
        theItems.push_back(_i);
    };


    /**
     * Deletes an item from the training set
     * Parameter: _idx   The index
     */
    virtual bool remove(unsigned int _idx)
    {

        if (_idx > theItems.size())
            return false;

        theItems.erase(theItems.begin() + _idx);
        if (isCalibrated)
            theTargets.erase(theTargets.begin() + _idx);

        return true;
    };


    /**
     * Returns the number of items in the training set
     */
    size_t size() const
    {
        return theItems.size();
    };


    /**
     * Returns a const reference to the specified target
     * Parameter: _i  The index
     * @exception out_of_range   If _i is out of range
     */
    const Target& targetAt(unsigned _i) const
    {
        if (!isCalibrated)
        {
            std::string msg;
            msg = "The training set is not calibrated.";
            throw std::out_of_range(msg);
        }
        if (_i >= size())
        {
            std::string msg;
            msg = "Out of range. No target at position " + integerToString(_i);
            throw std::out_of_range(msg);
        }

        return theTargets[_i];
    };


    /**
     * Returns a reference to a target.
     * Parameter: _i  The index
     * @exception out_of_range   If _i is out of range
     */
    Target& targetAt(unsigned _i)
    {
        if (_i >= size())
        {
            std::string msg;
            msg = "Out of range. No target at position " + integerToString(_i);
            throw std::out_of_range(msg);
        }

        return theTargets[_i];
    };



    /**
     * Returns a const reference to the specified item
     * Parameter: _i  The index
     * @exception out_of_range If _i is out of range
     */
    const Item& itemAt(unsigned _i) const
    {
        if (_i >= size())
        {
            std::string msg;
            msg = "Out of range. No item at position " + integerToString(_i);
            throw std::out_of_range(msg);
        }

        return theItems[_i];
    };



    /**
     * Returns a recerence to an item.
     * Parameter: _i     The index
     * Parameter: _item  The item
     * @exception out_of_range If _i is out of range
     */
    Item& itemAt(unsigned _i)
    {
        if (_i >= size())
        {
            std::string msg;
            msg = "Out of range. No item at position " + integerToString(_i);
            throw std::out_of_range(msg);
        }

        return theItems[_i];
    };



    /**
     * Returns true if the training set is calibreted or not
     */
    bool calibrated() const
    {
        return isCalibrated;
    }


    /**
     * Sets the calibrated flag.
     */
    void calibrated(const bool& _calib)
    {
        isCalibrated = _calib;
    };

    /**
     * Clears the training set: all items and targets are suppresed; all other variables remain the same
     */
    void clear()
    {
        theItems.clear();
        theTargets.clear();
        nTargets = 0;
    };


    /**
     * Standard output for a training set
     * Parameter: _os The output stream
     */
    virtual void printSelf(std::ostream& _os) const
    {
        writeItems(_os);
    }

    /**
     * Standard input for a training set
     * Parameter: _is The input stream
     * @exception  runtime_error  If there are problems with the stream
     * NOTE: This method is empty, it should be defined.
     */
    virtual void readSelf(std::istream& _is)
    {}

    /**
     * Saves the class into a stream.
     * this method can be used to save the status of the class.
     * Parameter: _os The output stream
     */
    virtual void saveObject(std::ostream& _os) const
    {
        writeCalibrated(_os);
        writeItems(_os, true);
    }

    /**
     * Loads the class from a stream.
     * this method can be used to load the status of the class.
     * Parameter: _is The output stream
     */
    virtual void loadObject(std::istream& _is)
    {
        clear();
        // first of all, check if the training set is calibrated or not
        checkCalibrated(_is);
        //afterwards, we have to read the item
        readItems(_is);
    }

    /**
     * Returns the number of different classes in the (calibrated) training set
     * @exception runtime_error  If the training set is not calibrated
     */
    unsigned numTargets() const
    {
        if (!calibrated())
            throw std::runtime_error("TS not calibrated");
        return nTargets;
    }

    /**
     * Swaps two items in the vector
     */
    virtual bool swapItems(unsigned _i, unsigned _j)
    {
        if (_i > size() || _j > size() || _i == _j)
            return false;
        swap(theItems[_i], theItems[_j]);
        if (isCalibrated)
            swap(theTargets[_i], theTargets[_j]);
        return true;
    }

    /** Standard output for a training set
      * Parameter: _os The output stream
      * Parameter: _ts  The training set to be printed
      */
    friend std::ostream& operator << (std::ostream& _os, const ClassificationTrainingSet& _ts)
    {
        _ts.printSelf(_os);
        return _os;
    }

    /**
     * Standard input for a training set
     * Parameter: _is The input stream
     * Parameter: _ts  The training set to be read
     * @exception  runtime_error  If there are problems with the stream
     */
    friend std::istream& operator >> (std::istream& _is, ClassificationTrainingSet& _ts)
    {
        _ts.readSelf(_is);
        return _is;
    }

protected:

    /**
     * Compute the number of different targets in the training set
     */
    void computeNumTargets()
    {
        if (calibrated())
        {
            typedef std::set< Target, std::less<Target> > STB;
            STB targetSet;
            for (unsigned i = 0; i < size(); i ++)
            {
                targetSet.insert(targetAt(i));
            }
            // Assign the number of targets
            nTargets = targetSet.size();
        }
        else
            nTargets = 0;
    }

    /**
     * Checks if the training set is calibrated or not
     * Parameter: _is  The input stream
     * @exception  runtime_error  If there are problems with the stream
     */
    void checkCalibrated(std::istream& _is)
    {
        skipComments(_is);

        if (_is)
        {
            std::string s;
            _is >> s;

            // Comments skipped, read calibrated
            if (_is)
            {
                std::string s2 = "";

                // uppercase s
                for (std::string::iterator i = s.begin() ; i != s.end() ; i++)
                    s2 += toupper(*i);

                if (s2 == "CALIBRATED")
                    isCalibrated = true;
                else
                    for (std::string::iterator i = s.end() ; i > s.begin() ; _is.putback(*(--i)));
            }
        }

        if (!_is)
        {
            std::string msg;
            msg = "Error reading the file";
            throw std::runtime_error(msg);
        }
    };


    /**
    * Read the items
    * Parameter: _is  The input stream
    * @exception  runtime_error  If there are problems with the stream
    */
    void readItems(std::istream& _is)
    {
        while (_is)
        {
            skipComments(_is);
            if (_is)
            {
                Item item;
                try
                {
                    _is >> item;
                }
                catch (std::exception&)
                {
                    std::string msg;
                    msg = "Error reading the item";
                    throw std::runtime_error(msg);
                }

                if (_is)
                {
                    theItems.push_back(item);
                    if (isCalibrated)
                    {
                        Target target;
                        try
                        {
                            char c;
                            if (_is) _is >> c;
                            // go back to the beginning of the line
                            if (_is) _is.putback(c);
                            // check for next line "<"
                            if (c == '<')
                                target = Target();
                            else
                                if (_is)
                                    _is >> target;
                                else
                                    target = Target();
                        }
                        catch (std::exception&)
                        {
                            std::string msg;
                            msg = "Error reading the item";
                            throw std::runtime_error(msg);
                        }

                        theTargets.push_back(target);
                    }
                    else
                        theTargets.push_back(Target());
                }
                else
                {
                    std::string msg;
                    msg = "Error reading the item";
                    throw std::runtime_error(msg);
                }
            }
        }
    };



    /**
     * Writes if the training set is calibrated or not
     * Parameter: _os  The output stream
     */
    void writeCalibrated(std::ostream& _os) const
    {
        if (isCalibrated)
            _os << "calibrated" << std::endl;
    };


    /**
     * Writes the items
     * Parameter: _os  The output stream
     * Parameter: _delim Flag to use "<" delimiters (false by default)
     */
    void writeItems(std::ostream& _os, bool _delim = false) const
    {
        typename std::vector<Item>::const_iterator i;
        typename std::vector<Target>::const_iterator j;
        for (i = theItems.begin(), j = theTargets.begin() ; i < theItems.end() ;
             i++, j++)
        {

            if (_delim)
                _os << *i;
            else
            {
                for (size_t d = 0; d < (*i).size(); d++)
                {
                    _os << (*i)[d];
                    if (d != (*i).size() - 1) _os << " ";
                }
            }

            /*      if (isCalibrated) {
                    if ((*j == " " || *j == "") && _delim)
                _os << "|";
                    else
                      _os << " " << *j;
                  }*/
            if (isCalibrated)
                _os << " " << *j;

            _os << std::endl;
        }
    };


    /**
     * Skip the comments if the stream
     * Parameter: _is  The input stream
     */
    void skipComments(std::istream& _is) const
    {
        char c;
        if (_is)
            _is >> c;

        // check for comments
        while (_is && c == '#')
        {
            // skip the comment
            while (_is && _is.get() != '\n');

            // read beginning of next line
            if (_is)
                _is >> c;
        }

        // go back to the beginning of the line
        if (_is)
            _is.putback(c);
    };


    /** Protected interface to the items
        @return a const_iterator that points to the first Item in the set
     */
    typename std::vector<Item>::const_iterator itemsBegin() const
    {
        return theItems.begin();
    }


    /** Protected interface to the items
     @return a const_iterator that points to the last Item in the set
     */
    typename std::vector<Item>::const_iterator itemsEnd() const
    {
        return theItems.end();
    }


    /** Protected interface to the targets
        @return a const_iterator that points to the first Target in the set
    */

    typename std::vector<Target>::const_iterator targetsBegin() const
    {
        return theTargets.begin();
    }

    /** Protected interface to the targets
        @return a const_iterator that points to the last Target in the set
    */
    typename std::vector<Target>::const_iterator targetsEnd() const
    {
        return theTargets.end();
    }
};
//@}
#endif
