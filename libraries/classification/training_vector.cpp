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
// ClassicTrainingVectors.cc
//-----------------------------------------------------------------------------

#include <cmath>

#include "training_vector.h"
#include <data/args.h>
#include <data/metadata.h>

/**
 * TrainingSet for ClassicTrainingVectors
 */

/**
 * Constructs a training set given a stream
 * Parameter: _is  The input stream
 * @exception  runtime_error  If there are problems with the stream
 */
ClassicTrainingVectors::ClassicTrainingVectors(std::istream & _is)
{
    try
    {
        clear();
        readSelf(_is);
    }
    catch (std::exception& e)
    {
        std::ostringstream msg;
        msg << e.what() << std::endl << "Error reading the training vector";
        throw std::runtime_error(msg.str());
    }
}



/**
 * Copy Constructor. Useful when returning a ClassicTrainingVectors Class.
 * Parameter: op1 ClassicTrainingVectors
 */
ClassicTrainingVectors::ClassicTrainingVectors(const ClassicTrainingVectors &op1)
{

    calibrated(op1.calibrated());

    for (size_t i = 0; i < op1.size(); i++)
        if (calibrated())
            add(op1.itemAt(i), op1.targetAt(i));
        else
            add(op1.itemAt(i));

    normalized = op1.normalized;
    varStats = op1.varStats;
}




/**
 * Returns amount of features
 */
unsigned ClassicTrainingVectors::featureSize() const
{
    return itemAt(0).size();
}

/**
 * Returns dimension (the same as above)
 */
unsigned ClassicTrainingVectors::dimension() const
{
    return itemAt(0).size();
}

/**
 * Clears the training set
 */
void ClassicTrainingVectors::clear()
{
    ClassificationTrainingSet<FeatureVector, Label>::clear();
    varStats.clear();
    normalized = false;
}

/**
 * Standard output for a training set
 * Parameter: _os The output stream
 * Parameter: _ts  The training set to be printed
 */
void ClassicTrainingVectors::printSelf(std::ostream& _os) const
{
    _os << dimension() << " " << theItems.size() << std::endl;
    ClassificationTrainingSet<FeatureVector, Label>::printSelf(_os);
}

/**
 * Standard input for a training set
 * Parameter: _is The input stream
 * Parameter: _ts  The training set to be read
 * @exception  runtime_error  If there are problems with the stream
 */
void ClassicTrainingVectors::readSelf(std::istream& _is)
{
#ifndef _NO_EXCEPTION
    try
    {
#endif
        clear();
        std::string line;

        // Determines the number of rows and columns in the training set

        long dim, size;
        _is >> dim;
        _is >> line;
        if (!sscanf(line.c_str(), "%ld", &size))
        {
            int x, y;
            _is >> x;
            _is >> y;
            size = x * y;
        }
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
            theTargets[i] = removeSpaces(line);
        }

#ifndef _NO_EXCEPTION
    }
    catch (std::exception& e)
    {
        std::ostringstream msg;
        msg << e.what() << std::endl << "Error reading the training set";
        throw std::runtime_error(msg.str());
    }
#endif
}

void ClassicTrainingVectors::read(const FileName& fnIn)
{
    clear();

    // Read header and content
    MetaData vectorHeader(formatString("vectorHeader@%s",fnIn.c_str()));
    MetaData vectorContent(formatString("vectorContent@%s",fnIn.c_str()));
    size_t Nvectors;
    size_t vectorSize;
    size_t id=vectorHeader.firstObject();
    vectorHeader.getValue(MDL_CLASSIFICATION_DATA_SIZE,vectorSize,id);
    vectorHeader.getValue(MDL_COUNT,Nvectors,id);
    theItems.reserve(Nvectors);
    theTargets.reserve(Nvectors);

    // Read the data
    FileName fnInRaw=formatString("%s.vec",fnIn.withoutExtension().c_str());
    std::ifstream fhInRaw(fnInRaw.c_str(),std::ios::binary);
    if (!fhInRaw)
        REPORT_ERROR(ERR_IO_NOTEXIST,fnInRaw);
    std::vector<floatFeature> v;
    v.resize(vectorSize);
    float *buffer=new float[vectorSize];
    String fnImg;
    size_t order;
    FOR_ALL_OBJECTS_IN_METADATA(vectorContent)
    {
        vectorContent.getValue(MDL_IMAGE,fnImg,__iter.objId);
        vectorContent.getValue(MDL_ORDER,order,__iter.objId);

        // Read raw values
        fhInRaw.seekg(order*vectorSize*sizeof(float));
        fhInRaw.read((char*)buffer,vectorSize*sizeof(float));
        if (!fhInRaw)
            REPORT_ERROR(ERR_IO_NOREAD,
                         formatString("Could not read image %lu from %s",
                                      order,fnInRaw.c_str()));
        for (size_t i=0; i<vectorSize; ++i)
        	v[i]=buffer[i];
        theTargets.push_back(fnImg);
        theItems.push_back(v);
    }
    delete []buffer;
}

/**
 * Saves the class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void ClassicTrainingVectors::saveObject(std::ostream& _os) const
{
    _os << dimension() << std::endl;
    _os << normalized << std::endl;
    if (normalized)
        for (size_t i = 0; i < varStats.size(); i++)
        {
            _os << varStats[i].mean << std::endl;
            _os << varStats[i].sd << std::endl;
        }
    ClassificationTrainingSet<FeatureVector, Label>::saveObject(_os);
}


/**
 * Loads the class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void ClassicTrainingVectors::loadObject(std::istream& _is)
{
    clear();
    int dim;
    _is >> dim;
    _is >> normalized;
    if (normalized)
        varStats.clear();
    varStats.resize(dim);
    for (size_t i = 0; i < varStats.size(); i++)
    {
        _is >> varStats[i].mean;
        _is >> varStats[i].sd;
    }
    ClassificationTrainingSet<FeatureVector, Label>::loadObject((std::istream&)_is);
}



/**
 * Deletes variable from Training set
 * Parameter: _var variable index
 */

void ClassicTrainingVectors::deleteVariable(int _var)
{
    for (unsigned int it = 0; it < size(); it++)
        itemAt(it).erase(itemAt(it).begin() + _var);
}


/**
 * Operator "="
 * Parameter: op1 ClassicTrainingVectors
 */
ClassicTrainingVectors& ClassicTrainingVectors::operator= (const ClassicTrainingVectors &op1)
{

    // This avoids memory leakage in assignments like v=v
    if (&op1 != this)
    {

        calibrated(op1.calibrated());

        for (size_t i = 0; i < op1.size(); i++)
            if (calibrated())
                add(op1.itemAt(i), op1.targetAt(i));
            else
                add(op1.itemAt(i));

        normalized = op1.normalized;
        varStats = op1.varStats;
    }
    return *this;
}


/** Copy the structure from another TS but leave it empty.
* Parameter: _ts ClassicTrainingVectors
* @note  Just the structure is copied, not the items or targets.
*/

bool ClassicTrainingVectors::copyStructure(ClassicTrainingVectors& _ts)
{

    // check if set is just initialized but empty

    if ((&_ts == this) || (size() + itemAt(0).size() != 0)) return false;
    calibrated(_ts.calibrated());
    normalized = _ts.normalized;
    varStats = _ts.varStats;
    return true;
}

/** Copy a row from an identical TS.
* Parameter: _ts ClassicTrainingVectors
* Parameter: _idx   row to be copied
* @note  No complete validation is done.
*/

bool ClassicTrainingVectors::insertRowFrom(ClassicTrainingVectors& _ts, unsigned int _idx)
{

    // just some validation, but not complete

    if (((&_ts == this) || (_idx > _ts.size())) ||
        (itemAt(0).size() != _ts.itemAt(0).size()))
        return false;

    if (calibrated())
        add(_ts.itemAt(_idx), _ts.targetAt(_idx));
    else
        add(_ts.itemAt(_idx));
    return true;

}

/** Delete a row from a TS.
* Parameter: _idx   row to be deleted
*/
bool ClassicTrainingVectors::deleteRow(unsigned int _idx)
{
    return remove(_idx);
}


/**
 * Normalize all features in the training set
 * Parameter: _i  The index to the feature
 */
void ClassicTrainingVectors::normalizeFeature(unsigned _i)
{
    using namespace std;

    // Do some validation

    if (_i > itemAt(0).size())
    {
        std::ostringstream msg;
        msg << "Out of range. No variable at position " << _i;
        throw std::out_of_range(msg.str());
    }

    // first calculates the mean
    floatFeature mean = 0;
    int nn = 0;
    for (size_t it = 0; it < size(); it++)
    {
        if (!isnan(itemAt(it)[_i]))
        {
            mean += itemAt(it)[_i];
            nn++;
        }

    }
    mean /= (floatFeature) nn;

    // Then calculates SD
    floatFeature sd = 0;
    for (size_t it = 0; it < size(); it++)
    {
        if (!isnan(itemAt(it)[_i]))
            sd += (itemAt(it)[_i] - mean) * (itemAt(it)[_i] - mean);
    }
    sd = sqrt(sd / (floatFeature)(nn - 1));

    // Now normalize the variable
    if (sd != 0)
    {
        for (size_t it = 0; it < size(); it++)
        {
            if (!isnan(itemAt(it)[_i]))
                itemAt(it)[_i] = (itemAt(it)[_i] - mean) / sd;
        }
    }

    varStats[_i].mean = mean;
    varStats[_i].sd = sd;
}


/**
 * Normalize all features in the training set
 */

void ClassicTrainingVectors::normalize()
{
    varStats.clear();
    varStats.resize(itemAt(0).size());
    for (unsigned i = 0; i < itemAt(0).size(); i++)
        normalizeFeature(i);
    normalized = true;
}


/**
 * UnNormalize all features in the training set
 */

void ClassicTrainingVectors::unNormalize()
{
    using namespace std;
    for (unsigned it = 0; it < size(); it++)
    {
        for (unsigned i = 0; i < itemAt(0).size(); i++)
        {
            if (!isnan(itemAt(it)[i]))
                itemAt(it)[i] = itemAt(it)[i] * varStats[i].sd + varStats[i].mean;
        }
    }
    varStats.clear();
    normalized = false;
}


/**
 * returns normalized variable in the original scale
 */

double ClassicTrainingVectors::getUnormalizedVar(unsigned _item, unsigned _var) const
{
    using namespace std;
    if (!normalized)
    {
        std::ostringstream msg;
        msg << "Variable is not normalized" << _var;
        throw std::runtime_error(msg.str());
    }

    if (_var > itemAt(0).size())
    {
        std::ostringstream msg;
        msg << "Out of range. No variable at position " << _var;
        throw std::out_of_range(msg.str());
    }

    if (_item > size())
    {
        std::ostringstream msg;
        msg << "Out of range. No item at position " << _var;
        throw std::out_of_range(msg.str());
    }

    double t;
    if (!isnan(itemAt(_item)[_var]))
        t = (double) itemAt(_item)[_var] * varStats[_var].sd + varStats[_var].mean;

    return t;
}

/**
 * Returns TRUE if recordset is normalized.
 */

bool ClassicTrainingVectors::isNormalized() const
{
    return normalized;
}


/**
  * Returns a const reference to the normalization vector
*/
/*  const std::vector<ClassicTrainingVectors::statsStruct>& ClassicTrainingVectors::getNormalizationInfo() const {
   return varStats;
  }*/


/**
 * Calcualtes the average and SD of a feature in the training set
 * Parameter: _i  The index to the feature
 */
void ClassicTrainingVectors::getFeatureStats(unsigned _i, floatFeature& _mean, floatFeature& _sd)
{
    using namespace std;

    // Do some validation
    if (_i > itemAt(0).size())
    {
        std::ostringstream msg;
        msg << "Out of range. No variable at position " << _i;
        throw std::out_of_range(msg.str());
    }

    // first calculates the mean
    _mean = 0;
    int nn = 0;
    for (size_t it = 0; it < size(); it++)
    {
        if (!isnan(itemAt(it)[_i]))
        {
            _mean += itemAt(it)[_i];
            nn++;
        }

    }
    _mean /= (floatFeature) nn;

    // Then calculates SD
    _sd = 0;
    for (size_t it = 0; it < size(); it++)
    {
        if (!isnan(itemAt(it)[_i]))
            _sd += (itemAt(it)[_i] - _mean) * (itemAt(it)[_i] - _mean);
    }
    _sd = sqrt(_sd / (floatFeature)(nn - 1));
}

/**
 * Returns a vector containing the average (item 0) and SD (item 1)
 */

ClassicTrainingVectors ClassicTrainingVectors::getStatVector()
{
    ClassicTrainingVectors myStatVector;
    myStatVector.theItems.resize(2);
    myStatVector.theItems[0].resize(itemAt(0).size(), 0);
    myStatVector.theItems[1].resize(itemAt(0).size(), 0);
    myStatVector.theTargets.resize(2);
    for (unsigned i = 0; i < itemAt(0).size(); i++)
        getFeatureStats(i, myStatVector.theItems[0][i], myStatVector.theItems[1][i]);
    myStatVector.theTargets[0] = "Average ";
    myStatVector.theTargets[1] = "SD ";
    return myStatVector;
}



//-----------------------------------------------------------------------------

