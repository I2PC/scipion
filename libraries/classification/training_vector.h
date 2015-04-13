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
// ClassicTrainingVectors.hh
//-----------------------------------------------------------------------------

#ifndef XMIPPCTVECTORS_H
#define XMIPPCTVECTORS_H

#include <sstream>
#include <data/xmipp_filename.h>
#include "data_types.h"
#include "training_set.h"
#include "vector_ops.h"

/**@defgroup TrainingVectors Training Vectors (Feature vectors) class
   @ingroup ClassificationLibrary */
//@{
/**
 * This class implements all the necessary functionality for classic
 * training vectors.
 */
class ClassicTrainingVectors : public ClassificationTrainingSet<FeatureVector, Label>
{
public:

    typedef struct stats
    {
        floatFeature mean;
        floatFeature sd;
    }
    statsStruct;


    /**
     * Default constructor
     Parameter: _vecSize Vector dimension; required to dim the feature and types vector
     Parameter: _calib calibration which should be true if the data set has labels
    */
    ClassicTrainingVectors(unsigned _vecSize = 0, bool _calib = true)
            : ClassificationTrainingSet<FeatureVector, Label>(_calib),  /*varStats(_vecSize),*/
            normalized(false)
    {}
    ;

    /**
     * Constructs a training set given a stream
     * Parameter: _is  The input stream
     * @exception  runtime_error  If there are problems with the stream
     */
    ClassicTrainingVectors(std::istream & _is);

    /**
     * Copy Constructor. Useful when returning a ClassicTrainingVectors Class.
     * Parameter: op1 ClassicTrainingVectors
     */
    ClassicTrainingVectors(const ClassicTrainingVectors &op1);

    /**
     * Virtual destructor
     */
    virtual ~ClassicTrainingVectors()
    {}
    ;

    /**
     * Returns the amount of features (cases)
     */
    unsigned featureSize() const;

    /**
     * Returns the dimension of the vectors (number of features)
     */
    unsigned dimension() const;

    /**
     * Clears the training set
     */
    void clear();

    /**
     * Standard output for a training set
     * Parameter: _os The output stream
     */
    virtual void printSelf(std::ostream& _os) const;

    /**
     * Standard input for a training set
     * Parameter: _is The input stream
     * @exception  runtime_error  If there are problems with the stream
     */
    virtual void readSelf(std::istream& _is);

    /** Read data from a metadata.
     */
    void read(const FileName & fnIn);

    /**
     * Saves the class into a stream.
     * this method can be used to save the status of the class.
     * Parameter: _os The output stream
     */
    virtual void saveObject(std::ostream& _os) const;


    /**
     * Loads the class from a stream.
     * this method can be used to load the status of the class.
     * Parameter: _is The output stream
     */
    virtual void loadObject(std::istream& _is);


    /**
     * Deletes a variable (feature) from Training set
     * Parameter: _var variable index
     */
    void deleteVariable(int _var);

    /**
     * Operator "="
     * Parameter: op1 ClassicTrainingVectors
     */
    ClassicTrainingVectors& operator= (const ClassicTrainingVectors &op1);


    /** Copy the structure from another TS but leave it empty.
    * Parameter: _ts ClassicTrainingVectors
    * @note  Just the structure is copied, not the items or targets.
    */
    bool copyStructure(ClassicTrainingVectors& _ts);


    /** Copy a row from an identical TS.
    * Parameter: _ts ClassicTrainingVectors
    * Parameter: _idx   row to be copied
    * @note  No complete validation is done.
    */
    bool insertRowFrom(ClassicTrainingVectors& _ts, unsigned int _idx);

    /** Delete a row from a TS.
    * Parameter: _idx   row to be deleted
    */
    bool deleteRow(unsigned int _idx);


    /**
    * Normalize a feature in the training set
    * Parameter: _i  The index to the feature
    */
    virtual void normalizeFeature(unsigned _i);


    /**
    * Normalize all features in the training set
    */
    void normalize();

    /**
    * UnNormalize all features in the training set
    */
    virtual void unNormalize();

    /**
    * returns normalized variable in the original scale
    */
    virtual double getUnormalizedVar(unsigned _item, unsigned _var) const;

    /**
    * Returns TRUE if recordset is normalized.
    */
    bool isNormalized() const;

    /**
    * Returns a const reference to the normalization vector
    */
    virtual const std::vector<statsStruct>& getNormalizationInfo() const
    {
        return varStats;
    };

    /**
     * Calcualtes the average and SD of a feature in the training set
     * Parameter: _i  The index to the feature
     */
    void getFeatureStats(unsigned _i, floatFeature& _mean, floatFeature& _sd);

    /**
    * Returns a vector containing the average (item 0) and SD (item 1)
    */

    ClassicTrainingVectors getStatVector();


protected:
    std::vector <statsStruct>  varStats;
    bool    normalized;

};
//@}
#endif//XMIPPCTVECTORS_H
