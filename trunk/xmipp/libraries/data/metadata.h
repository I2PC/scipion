/***************************************************************************
 * 
 * Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

#ifndef METADATA_H
#define METADATA_H

#include <map>
#include <vector>
#include <iostream>
#include <iterator>
#include <sstream>
#include "funcs.h"
#include "strings.h"
#include "metadata_sql.h"
#include <time.h>
#include <stdarg.h>
#include <stdio.h>
#include <sstream>

//empty vector used for inizialization

/// @defgroup MetaData Metadata management
/// @ingroup DataLibrary

/// @defgroup MetaDataClass Metadata class management
/// @ingroup MetaData

/** MetaData Manager.
 * @ingroup MetaDataClass
 *
 * The MetaData class manages all procedures related to
 * metadata. MetaData is intended to group toghether old
 * Xmipp specific files like Docfiles, Selfiles, etc..
 * 
 */

//useful to init values to zero
static double zeroD=0.;
static double oneD=1.;
static bool  falseb=false;

enum AggregateMode
{
    KEEP_OLD,
    KEEP_NEW,
    SUM
};

#define FOR_ALL_OBJECTS_IN_METADATA(kkkk_metadata) \
        for(long int kkkk = (kkkk_metadata).iteratorBegin(); \
             !(kkkk_metadata).iteratorEnd(); \
             kkkk=(kkkk_metadata).iteratorNext())

#define FOR_ALL_OBJECTS_IN_METADATA2(kkkk_metadata, jjjj_metadata) \
        for(long int kkkk = (kkkk_metadata).iteratorBegin(), jjjj = (jjjj_metadata).iteratorBegin(); \
             !(kkkk_metadata).iteratorEnd() && !(jjjj_metadata).iteratorEnd();\
             kkkk=(kkkk_metadata).iteratorNext(), jjjj=(jjjj_metadata).iteratorNext())



class MDQuery;
class MDSql;

class MetaData
{
private:
    //std::map<long int, MetaDataContainer *> objects; ///< Effectively stores all metadata

    // Used by firstObject, nextObject and lastObject to keep a pointer
    // to the "active" object. This way when you call setValue without
    // an objectId, this one is chosen
    //std::map<long int, MetaDataContainer *>::iterator objectsIterator;

    // Allows a fast search for pairs where the value is
    // a string, i.e. looking for filenames which is quite
    // usual
    std::map<std::string, long int> fastStringSearch;
    MDLabel fastStringSearchLabel;

    std::string path; ///< A parameter stored on MetaData Files
    std::string comment; ///< A general comment for the MetaData file

    bool isColumnFormat; ///< Format for the file, column or row formatted

    /**Input file name
     * Where does this MetaData come from/go to be stored?
     */
    FileName inFile;

    /** What labels have been read from a docfile/metadata file
     *   and/or will be stored on a new metadata file when "save" is
     *   called
     **/
    std::vector<MDLabel> activeLabels;

    /** When reading a column formated file, if a label is found that
     *   does not exists as a MDLabel, it is ignored. For further
     *   file processing, such columns must be ignored and this structure
     *   allows to do that
     **/
    std::vector<unsigned int> ignoreLabels;

    /** This variables should only be used by MDSql
     * for handling db status of metadata
     */
    /** The table id to do db operations */
    MDSql * myMDSql;

    /** The id of the object that is active
     * usefull for calling 'setValue' and 'getValue'
     * without specifying object id
     * when activeObjId = -1 means that aren't active object
     */
    int activeObjId;

    /** This are for iteration */
    int iterIndex;
    std::vector<long int> *iterObjectsId;

    /** Read, read data from an input stream and fill metadata
     * @ingroup MetaDataIO
     */
    //void read(std::ifstream *infile, std::vector<MDLabel> *labelsVector = NULL);

    /** Init, do some initializations tasks, used in constructors
     * @ingroup MetaDataConstructors
     */
    void init(const std::vector<MDLabel> *labelsVector = NULL);

    /** Copy info variables from another metadata
     * @ingroup MetaDataConstructors
     */
    void copyInfo(const MetaData &md);

    /** Copy all data from another metadata
     * @ingroup MetaDataConstructors
     */
    void copyMetadata(const MetaData &md);

    /** This private functions are for set real values
     * there is an explicit function signature
     * foreach type supported in Metadata.
     * This is done for some type checking of Metadata labels
     * and values
     */

    bool _setValue(long int objId, const MDValue &mdValueIn);
    bool _getValue(long int objId, MDValue &mdValueOut) const;

    /** This have the same logic of the public one,
     * but doesn't perform any range(wich implies do a size()) checks.
     */
    void _selectSplitPart(const MetaData &mdIn,
                          int n, int part, long int mdSize,
                          const MDLabel sortLabel);

    /** This function is for generalizate the sets operations
     * of unionDistinct, intersection, substraction
     * wich can be expressed in terms of
     * ADD, SUBSTRACT of intersection part
     */
    void _setOperates(const MetaData &mdIn, const MDLabel label, SetOperation operation);
    void _setOperates(const MetaData &mdInLeft, const MetaData &mdInRight, const MDLabel label, SetOperation operation);
    /** clear data and table structure */
    void _clear(bool onlyData=false);

    long int _iteratorBegin(const MDQuery *query = NULL);

    /** Some private reading functions */
    void _readColumns(std::istream& is, std::vector<MDValue>& columnValues,
                               std::vector<MDLabel>* desiredLabels = NULL);
    void _readRows(std::istream& is, std::vector<MDValue>& columnValues, bool useCommentAsImage);
    void _readRowFormat(std::istream& is);

public:

    /// @defgroup MetaDataConstructors Constructors for MetaData objects
    /// @ingroup MetaDataClass

    /** Empty Constructor.
     * @ingroup MetaDataConstructors
     *
     * The MetaData is created with no data stored on it. You can fill in it programmatically
     * or by a later reading from a MetaData file or old Xmipp formatted type.
     * if labels vectors is passed this labels are created on metadata
     */
    MetaData();
    MetaData(const std::vector<MDLabel> *labelsVector);

    /** From File Constructor.
     * @ingroup MetaDataConstructors
     *
     * The MetaData is created and data is read from provided fileName. Optionally, a vector
     * of labels can be provided to read just those required labels
     */
    MetaData(const FileName &fileName, const std::vector<MDLabel> *labelsVector = NULL);

    /** Copy constructor
     * @ingroup MetaDataConstructors
     *
     * Created a new metadata by copying all data from an existing MetaData object.
     */
    MetaData(const MetaData &md);

    /** Assignment operator
     * @ingroup MetaDataConstructors
     *
     * Copies MetaData from an existing MetaData object.
     */
    MetaData& operator =(const MetaData &md);

    /** Destructor
     * @ingroup MetaDataConstructors
     *
     * Frees all used memory and destroys object.
     */
    ~MetaData();

    /**Clear clear all data
     * @ingroup MetaDataConstructors
     */
    void clear();

    /// @defgroup GettersAndSetters some getters and setters functions
    /// @ingroup MetaDataClass
    /// all getter should be 'const'

    /**Get column format info.
     *@ingroup GettersAndSetters
     */
    bool getColumnFormat() const;

    /** Set to false for row format (parameter files)
     *  @ingroup GettersAndSetters
     *  set to true  for column format (this is the default) (docfiles)
     *
     */
    void setColumnFormat(bool column);
    /**Get path info.
     *@ingroup GettersAndSetters
     */
    std::string getPath()   const ;
    /**Set Path
     * @ingroup GettersAndSetters
     * will appear in first line
     */
    void setPath(std::string newPath = "");
    /**Get Header Comment
     *@ingroup GettersAndSetters
     *will appear in second line
     */
    std::string getComment() const;
    /**Set Header Comment
     *@ingroup GettersAndSetters
     *will appear in second line
     */
    void setComment(const std::string newComment = "No comment");

    /**Return metadata filename
     *@ingroup GettersAndSetters
     */
    FileName getFilename();

    /**Get safe access to active labels
     * @ingroup GettersAndSetters
     */
    std::vector<MDLabel> getActiveLabels() const;

    /**Get safe access to active labels
     * @ingroup GettersAndSetters
     */
    std::vector<MDLabel> *geActiveLabelsAddress() const;

    /**Get maximum string length of this value
    * @ingroup GettersAndSetters
    */

    /**Get the active object id
     * -1 will be returned if no object is active
     */
    long int getActiveObject();

    int MaxStringLength( const MDLabel thisLabel) const;

    /// @defgroup DataAccess Access to MetaData data for read or write values
    /// @ingroup MetaDataClass

    /** Set the value for some label in the object that has id 'objectId'
     * @ingroup DataAccess
     */
    template<class T>
    bool setValue(const MDLabel label, const T &valueIn, long int objectId=-1)
    {
        _setValue(objectId, MDValue(label, valueIn));
    }

    bool setValues(int n, ...);
    // Set a new pair/value for an specified object. If no objectId is given, that
    // pointed by the class iterator is used
    bool setValueFromStr(const MDLabel label, const std::string &value, long int objectId = -1);
    bool getStrFromValue(const MDLabel label, std::string &strOut, long int objectId = -1);


    /**IsEmpty check whether the metadata is empty or not
     * @ingroup DataAccess
     */
    bool isEmpty() const;

    /**size return the number of objects contained in the metadata
     * @ingroup DataAccess
     */
    long int size() const;

    /** Check whether a label is contained
     * in the metadata
     */
    bool containsLabel(const MDLabel label) const;
    /** Add a new label to the metadata
     *
     */
    bool addLabel(const MDLabel label);

    /** Adds a new, empty object to the objects map. If objectId == -1
     * @ingroup DataAccess
     *   the new ID will be that for the last object inserted + 1, else
     *   the given objectId is used. If there is already an object whose
     *   objectId == input objectId, just removes it and creates an empty
     *   one
     **/
    long int addObject(long int objectId = -1);

    /** Import objects from another metadata
     * @ingroup DataAccess
     */
    void importObject(const MetaData &md, const long int objId, bool doClear=true);
    void importObjects(const MetaData &md, const std::vector<long int> &objectsToAdd, bool doClear=true);
    void importObjects(const MetaData &md, const MDQuery &query, bool doClear=true);

    /**Remove the object with this id
    * @ingroup DataAccess
    * Returns true if the object was removed or false if
    * the object did not exist
    */
    bool removeObject(long int objectId);

    /** Removes the collection of objects of given vector id's
     * @ingroup DataAccess
     * NOTE: The iterator will point to the first object after any of these
     * operations
     */
    void removeObjects(const std::vector<long int> &toRemove);

    /** Removes objects with pair <label, value>
     * return the number of deleted rows
     * if not query, all objectes are removed
     * @ingroup DataAccess
     */
    int removeObjects(const MDQuery &query);
    int removeObjects();
    //FIXME: organize this
    // Possible error codes for the metadata operations
    enum MDErrors
    {
        NO_OBJECTS_STORED = -1, // NOTE: Do not change this value (-1)
        NO_MORE_OBJECTS = -2,
        NO_OBJECT_FOUND = -3
    };

    /** Indexes add and remove for fast search
     * in other labels, but insert are more expensive
     */
    void addIndex(MDLabel label);
    void removeIndex(MDLabel label);

    ///@defgroup MetaDataIteration Some functions related with iteration over metadata objects
    ///@ingroup MetaDataClass

    /**firstObject, move the 'activeObject' to the first object in the metadata
     *@ingroup MetaDataIteration
     * */
    long int firstObject();
    long int lastObject();
    long int nextObject();
    long int previousObject();
    long int goToObject(long int objectId);

    /** This will start iterating over all objects
     * the first object id will be returned
     */
    long int iteratorBegin();

    /**Same as previous but iterating over a subset of
     * objects
     */
    long int iteratorBegin(const MDQuery &query);

    /** Check whether the iteration has finished */
    bool iteratorEnd();

    /** Move to next object on iteration
     * return nextObject id
     */
    long int iteratorNext();


    ///@defgroup MetaDataSearch Some functions for perform searches on metadata objects
    ///@ingroup MetaDataClass

    /** Find all objects that match the query
     * if called without query, all objects are returned
     * if limit is provided only return a maximun of 'limit'
     * @ingroup MetaDataSearch
     */
    void findObjects(std::vector<long int> &objectsOut, const MDQuery &query, int limit = -1);
    void findObjects(std::vector<long int> &objectsOut, int limit = -1);
    /**Count all objects with pairs <label, value>
     * @ingroup MetaDataSearch
     */
    int countObjects(const MDQuery &query);

    /** Find if the object with this id is present in the metadata
     * @ingroup MetaDataSearch
     */
    bool containsObject(long int objectId);

    /**Check if exists at least one object with pair <label, value>
     * @ingroup MetaDataSearch
     */
    bool containsObject(const MDQuery &query);

    /**Move active object to the first
     * object with pair <label, value in range> if exists
     * @ingroup MetaDataSearch
     */
    long int gotoFirstObject(const MDQuery &query);


    ///@defgroup MetaDataIO functions related to the I/O of metadata
    ///@ingroup MetaDataClass

    /** Write, write metadata to disk, using filename
     * @ingroup MetaDataIO
     */
    void write(const FileName &outFile);
    /** Write, write metadata to out stream
     * @ingroup MetaDataIO
     */
    void write(std::ostream &os);

    /** Read, read data from file, using filename
     * @ingroup MetaDataIO
     */
    void read(const FileName &inFile, std::vector<MDLabel> *labelsVector = NULL);
    /** Aggregate metadata objects,
     * @ingroup SetOperations
     * result in calling metadata object
     * thisLabel label is used for aggregation, second. Valid operations are:
     *
     * MDL_AVG:  The avg function returns the average value of all  operationLabel within a group.
      The result of avg is always a floating point value as long as at there
      is at least one non-NULL input even if all inputs are integers.
       The result of avg is NULL if and only if there are no non-NULL inputs.

      MDL_COUNT: The count function returns a count of the number of times that operationLabel is in a group.

      MDL_MAX       The max aggregate function returns the maximum value of all values in the group.

      MDL_MIN       The min aggregate function returns the minimum  value of all values in the group.

     MDL_SUM The total aggregate functions return sum of all values in the group.
     If there are no non-NULL input rows then returns 0.0.

     The result of total() is always a floating point value.
     */
    /** Aggregate modes */

    void aggregate(const MetaData &mdIn, AggregateOperation op,
                   MDLabel aggregateLabel, MDLabel operateLabel, MDLabel resultLabel);
    void aggregate(const MetaData &mdIn, const std::vector<AggregateOperation> &ops,
                             MDLabel operateLabel, const std::vector<MDLabel> &resultLabels);

    /// @defgroup SetOperations Set operations on MetaData's
    /// @ingroup MetaDataClass
    /** union of  metadata objects,
     *@ingroup SetOperations
     * result in calling metadata object
     * union is a reserved word so I called this method unionDistinct
     */
    void unionDistinct(const MetaData &mdIn, const MDLabel label=MDL_OBJID);

    /** union of  metadata objects,
     * @ingroup SetOperations
     * result in calling metadata object
     * Repetion are allowed
     */
    void unionAll(const MetaData &mdIn);

    /** merge of a metadata
     * @ingroup SetOperations
     * This function reads another metadata and makes a union to this one
     */
    void merge(const FileName &fn);

    /** intersects two metadata objects,
     *@ingroup SetOperations
     * result in "calling" metadata
     */
    void intersection(const MetaData &mdIn, const MDLabel label);

    /** substract two metadata objects,
     * @ingroup SetOperations
     * result in "calling" metadata
     */
    void substraction(const MetaData &mdIn, const MDLabel label);

    /** join two metadata
     * @ingroup SetOperations
     * result in "calling" metadata
     */
    void join(const MetaData &mdInLeft, const MetaData &mdInRight, const MDLabel label, JoinType type=LEFT);

    /** basic operations on columns data
     * Mainly perform replacements on string values and
     * basic algebraic operations on numerical ones.
     */
    void operate(const std::string &expression);

    //---------TO ORGANIZE-------------------------
    /** Randomize this metadata, MDin is input
    */
    void randomize(MetaData &MDin);
    /*
    * Sort this metadata, by label
    * dirty implementation using sqlite
    */
    void sort(MetaData &MDin, const MDLabel sortLabel);


    /** Select only a part of this MetaData
     * the metadata will be divided in 'n'
     * almost equally parts and will take the '
    *
    */
    void split(int n, std::vector<MetaData> &results,
               const MDLabel sortLabel=MDL_OBJID);
    /** SelectSplitPart have the same logic than
     * split, but only one of the parts is taken
     * another difference is that result will be
     * in the calling object
     */
    void selectSplitPart(const MetaData &mdIn,
                         int n, int part,
                         const MDLabel sortLabel=MDL_OBJID);

    /** Select a part of the metadata in
     * from an starting position and
     * selecting an specified number of objects
     * if the numberOfObjects is -1, all objects
     * will be returned from startPosition
    */
    void selectPart (const MetaData &mdIn, long int startPosition, long int numberOfObjects,
                     const MDLabel sortLabel=MDL_OBJID);


    template<class T>
    bool getValue(const MDLabel label, T &valueOut, long int objectId = -1) const
    {
        try
        {
            MDValue mdValueOut(label);
            _getValue(objectId, mdValueOut);
            mdValueOut.getValue(valueOut);
            return true;
        }
        catch (Xmipp_error xe)
        {
            return false;
        }
    }

    /**Get all values of a column as a vector
     *
     */
    template<class T>
    void getColumnValues(const MDLabel label, std::vector<T> &valuesOut)
    {
        T value;
        MDValue mdValueOut(label);
        std::vector<long int> objectsId;
        findObjects(objectsId);
        int n = objectsId.size();
        valuesOut.resize(n);
        for (int i = 0; i < n; i++)
        {
            _getValue(objectsId[i], mdValueOut);
            mdValueOut.getValue(value);
            valuesOut[i] = value;
        }

    }
    friend class MDSql;
};

#endif
