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
#include <strings.h>

/** @defgroup MetaData Metadata Stuff
 * @ingroup DataLibrary
 * @{
 */

static double zeroD=0.;
static double oneD=1.;
static bool  falseb=false;

/** Write mode
 */
typedef enum
{
    OVERWRITE, //forget about the old file and overwrite it
    APPEND,    //append a data_ at the file end or replace an existing one
} WriteModeMetaData;
/** Iterate over all elements in MetaData
 *
 * This macro is used to generate loops over all elements in the MetaData.
 * At each iteration the 'active object' is changed so you can perform
 * the set and get default method on the MetaData.
 *
 * @code
 * MetaData md;
 *   //...
 * FOR_ALL_OBJECTS_IN_METADATA(md)
 * {
 *     std::string imageFile;
 *     md.getValue(MDL_IMAGE, imageFile);
 *     std::cout << "Image file: " << imageFile << " ";
 * }
 * @endcode
 */
#define FOR_ALL_OBJECTS_IN_METADATA(kkkk_metadata) \
        for(long int kkkk = (kkkk_metadata).iteratorBegin(); \
             !(kkkk_metadata).iteratorEnd(); \
             kkkk=(kkkk_metadata).iteratorNext())

/** Iterate over all elements of two MetaData at same time.
 *
 * This macro is useful to iterate over two MetaData with the same
 * number of elements and performs operations to elements in both of them.
 * At each iteration the 'active objects' in both MetaData are changed.
 *
 * @code
 * MetaData mdA, mdB, mdC;
 *  //Iterate over MetaData mdA and mdB
 *  //take image from the first and tilt angle from the second
 *  //and create a new MetaData.
 * FOR_ALL_OBJECTS_IN_METADATA2(mdA, mdB)
 * {
 *     std::string imageFile;
 *     double angle;
 *     mdA.getValue(MDL_IMAGE, imageFile);
 *     mdB.getValue(MDL_ANGLE_TILT, angle);
 *     mdC.addObject();
 *     mdC.setValue(MDL_IMAGE, imageFile);
 *     mdC.setValue(MDL_ANGLE_TILT, angle);
 * }
 * @endcode
 */
#define FOR_ALL_OBJECTS_IN_METADATA2(kkkk_metadata, jjjj_metadata) \
        for(long int kkkk = (kkkk_metadata).iteratorBegin(), jjjj = (jjjj_metadata).iteratorBegin(); \
             !(kkkk_metadata).iteratorEnd() && !(jjjj_metadata).iteratorEnd();\
             kkkk=(kkkk_metadata).iteratorNext(), jjjj=(jjjj_metadata).iteratorNext())



class MDQuery;
class MDSql;

/** Class to manage data files.
 *
 * The MetaData class manages all procedures related to
 * metadata. MetaData is intended to group toghether old
 * Xmipp specific files like Docfiles, Selfiles, etc..
 *
 */
class MetaData
{
protected:
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
     * and/or will be stored on a new metadata file when "save" is
     * called
     **/
    std::vector<MDLabel> activeLabels;

    /** When reading a column formated file, if a label is found that
     * does not exists as a MDLabel, it is ignored. For further
     * file processing, such columns must be ignored and this structure
     * allows to do that
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

    //Set the value of all objects in an specified column (both value and colum are specified in mdValueIn)
    bool _setValueCol(const MDObject &mdValueIn);

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
    void _readColumns(std::istream& is, MDRow & columnValues,
                      const std::vector<MDLabel>* desiredLabels = NULL);
    void _readRows(std::istream& is, MDRow & columnValues, bool useCommentAsImage);
    /** This function will be used to parse the rows data in START format
     * @param[out] columnValues MDRow with values to fill in
     * @param pchStar pointer to the position of '_loop' in memory
     * @param pEnd  pointer to the position of the next '_data' in memory
     */
    void _readRowsStar(MDRow & columnValues, char * pchStart, char * pEnd);
    void _readRowFormat(std::istream& is);

public:

    /** @name Constructors
     *  @{
     */

    /** Empty Constructor.
     *
     * The MetaData is created with no data stored on it. You can fill in it programmatically
     * or by a later reading from a MetaData file or old Xmipp formatted type.
     * if labels vectors is passed this labels are created on metadata
     */
    MetaData();
    MetaData(const std::vector<MDLabel> *labelsVector);

    /** From File Constructor.
     *
     * The MetaData is created and data is read from provided FileName. Optionally, a vector
     * of labels can be provided to read just those required labels
     */
    MetaData(const FileName &fileName, const std::vector<MDLabel> *desiredLabels = NULL);

    /** Copy constructor
     *
     * Created a new metadata by copying all data from an existing MetaData object.
     */
    MetaData(const MetaData &md);

    /** Assignment operator
     *
     * Copies MetaData from an existing MetaData object.
     */
    MetaData& operator =(const MetaData &md);

    /** Destructor
     *
     * Frees all used memory and destroys object.
     */
    ~MetaData();

    /**Clear all data
     */
    void clear();
    /** @} */

    /** @name Getters and setters
     * @{
     */

    /**Get column format info.
     */
    bool getColumnFormat() const;

    /** Set to false for row format (parameter files).
     *  set to true  for column format (this is the default) (docfiles)
     */
    void setColumnFormat(bool column);

    /** Check if the file (not the object) is in column format
     *  returns pointer do first two data_entries and firts loop
     */
    bool isColumnFormatFile(char * map,
                            char ** firstData,
                            char ** secondData,
                            char ** firstloop,
                            const char * blockName);

    /* Helper function to parse an MDObject and set its value.
     * The parsing will be from an input stream(istream)
     * and if parsing fails, an error will be raised
     */
    void _parseObject(std::istream &is, MDObject &object);

    /** Get Metadata labels for the block defined by start
     * and end loop pointers. Return pointer to newline after last label
     */
    char * _readColumnsStar(char * start,
                            char * end,
                            MDRow & columnValues,
                            const std::vector<MDLabel>* desiredLabels);
#ifdef NEVERDEFINED
    /** This function will read the possible columns and values from the file
     * in ROW format
     * and mark as MDL_UNDEFINED those who aren't valid labels
     * or those who appears in the IgnoreLabels vector
     * also set the activeLabels (for new STAR files)
     */

    char * _readRowFormatStar(char * pStart,
                              char * pEnd,
                              MDRow & columnValues,
                              const std::vector<MDLabel>* desiredLabels);
#endif
    /**Get path.
     */
    std::string getPath() const ;

    /**Set Path.
     * the path will appear in first line
     */
    void setPath(std::string newPath = "");

    /**Get Header Comment.
     * the comment will appear in second line.
     */
    std::string getComment() const;

    /**Set Header Comment.
     * the comment will appear in second line
     */
    void setComment(const std::string &newComment = "No comment");

    /**Get metadata filename.
     */
    FileName getFilename() const;

    /**Set metadata filename.
     */
    void setFilename(const FileName _filename);

    /**Get safe access to active labels.
     */
    std::vector<MDLabel> getActiveLabels() const;

    /**Get unsafe pointer to active labels.
     */
    std::vector<MDLabel> *geActiveLabelsAddress() const;

    /**Get the active object id.
     * -1 will be returned if no object is active
     */
    long int getActiveObject() const;

    /**Get maximum string length of column values.
    */
    int MaxStringLength( const MDLabel thisLabel) const;

    /** @} */

    /** @name MetaData Manipulation
     * @{
     */

    /**Set the value of all objects in an specified column.
     * @code
     * MetaData md;
     * md.setValueCol(MDL_IMAGE, "images/image00011.xmp");
     * @endcode
     */

    template<class T>
    bool setValueCol(const MDLabel label, const T &valueIn)
    {
        _setValueCol(MDObject(label, valueIn));
    }
    /** Set the value for some label.
     * to the object that has id 'objectId'
     * or to 'activeObject' if is objectId=-1.
     * This is one of the most used functions to programatically
     * fill a metadata.
     * @code
     * MetaData md;
     * md.addObject();
     * md.setValue(MDL_IMAGE, "images/image00011.xmp");
     * md.setValue(MDL_ANGLE_ROT, 0.);
     * @endcode
     */
    template<class T>
    bool setValue(const MDLabel label, const T &valueIn, long int objectId=-1)
    {
        setValue(MDObject(label, valueIn), objectId);
    }

    /** This functions are using MDObject for set real values
     * there is an explicit function signature
     * foreach type supported in Metadata.
     * This is done for some type checking of Metadata labels
     * and values
     */
    bool setValue(const MDObject &mdValueIn, long int objId = -1);
    bool getValue(MDObject &mdValueOut, long int objId = -1) const;

    /** Get the value of some label.
     * from the object that has id 'objectId'
     * or from 'activeObject' if objectId=-1.
     * @code
     * MetaData md;
     * md.read("images.xmd");
     * FileName imageFn;     *
     * FOR_ALL_OBJECTS_IN_METADATA(md)
     * {
     *      md.getValue(MDL_IMAGE, imageFn);
     *      std::out << "Image: " << imageFn);
     * }
     * @endcode
     */
    template<class T>
    bool getValue(const MDLabel label, T &valueOut, long int objectId = -1) const
    {
        try
        {
            MDObject mdValueOut(label);
            if (!getValue(mdValueOut,objectId))
                return false;
            mdValueOut.getValue(valueOut);
            return true;
        }
        catch (XmippError xe)
        {
            return false;
        }
    }

    /** Get all values of a column as a vector.
     */
    template<class T>
    void getColumnValues(const MDLabel label, std::vector<T> &valuesOut)
    {
        T value;
        MDObject mdValueOut(label);
        std::vector<long int> objectsId;
        findObjects(objectsId);
        int n = objectsId.size();
        valuesOut.resize(n);
        for (int i = 0; i < n; i++)
        {
            getValue(mdValueOut, objectsId[i]);
            mdValueOut.getValue(value);
            valuesOut[i] = value;
        }

    }

    /** Get all values of an MetaData row of an specified objId*/
    bool getRow(MDRow &row, long int objId = -1);

    /** Set label values from string representation.
     */
    bool setValueFromStr(const MDLabel label, const std::string &value, long int objectId = -1);

    /** Get string representation from label value.
     */
    bool getStrFromValue(const MDLabel label, std::string &strOut, long int objectId = -1);

    /**Check whether the metadata is empty.
     */
    bool isEmpty() const;

    /**Number of objects contained in the metadata.
     */
    long int size() const;

    /** Check whether a label is contained in metadata.
     */
    bool containsLabel(const MDLabel label) const;

    /** Add a new label to the metadata.
     * By default the label is added at the end,
     * if the position is specified and is between 0 and n-1
     * the new label is inserted at that position.
     */
    bool addLabel(const MDLabel label, int pos = -1);

    /** Adds a new, empty object to the objects map. If objectId == -1
     * the new ID will be that for the last object inserted + 1, else
     * the given objectId is used. If there is already an object whose
     * objectId == input objectId, just removes it and creates an empty
     * one
     */
    long int addObject(long int objectId = -1);

    /** Import objects from another metadata.
     * @code
     * //Import object 1000 from metadata B into metadata A
     * A.importObject(B, 1000);
     * //Import all objects with rotational angle greater that 60
     * A.importObjects(B, MDValuesGT(MDL_ANGLE_ROT, 60));
     * //Import all objects
     * A.importObjects(B);     *
     * @endcode
     */
    void importObject(const MetaData &md, const long int objId, bool doClear=true);
    void importObjects(const MetaData &md, const std::vector<long int> &objectsToAdd, bool doClear=true);
    void importObjects(const MetaData &md, const MDQuery &query, bool doClear=true);

    /** Remove the object with this id.
    * Returns true if the object was removed or false if
    * the object did not exist
    */
    bool removeObject(long int objectId);

    /** Removes the collection of objects of given vector id's
     * NOTE: The iterator will point to the first object after any of these
     * operations
     */
    void removeObjects(const std::vector<long int> &toRemove);

    /** Removes objects from metadata.
     * return the number of deleted rows
     * if not query, all objectes are removed
     * Queries can be used in the same way
     * as in the importObjects function
     */
    int removeObjects(const MDQuery &query);
    int removeObjects();

    /** Add and remove indexes for fast search
     * in other labels, but insert are more expensive
     */
    void addIndex(MDLabel label);
    void removeIndex(MDLabel label);

    /** @} */

    /** @name Iteration functions
     * @{
     */

    /** Goto first metadata object.
     * This function changes he 'activeObject'
     * to the first object in the metadata
     * so all calls to setValue() and getValue() will
     * be performed to this object.
     */
    long int firstObject();

    /** Goto last metadata object.*/
    long int lastObject();

    /** Goto next object to the active object.*/
    long int nextObject();

    /** Goto previous object to the active object.*/
    long int previousObject();

    /** Goto to specific object */
    long int goToObject(long int objectId);

    /** Goto to the first object that match the query.
     * This function will change the active object
     * to the first object that match the query if
     * the query results are non-empty
     */
    long int gotoFirstObject(const MDQuery &query);

    /** Starts iterating over all objects.
     * the first object id will be returned.
     * Also the active object will be changed
     * to the first object.
     */
    long int iteratorBegin();

    /**Starts iterating over a subset of objects.
     * the first object of the subset will
     * be returned and active object changed.
     */
    long int iteratorBegin(const MDQuery &query);

    /** Check whether the iteration has finished.
     */
    bool iteratorEnd() const;

    /** Move to next object on iteration.
     * return nextObject id
     */
    long int iteratorNext();
    /** @} */

    /** @name Search operations
     * @{
     */

    /** Find all objects that match a query.
     * if called without query, all objects are returned
     * if limit is provided only return a maximun of 'limit'
     */
    void findObjects(std::vector<long int> &objectsOut, const MDQuery &query);
    void findObjects(std::vector<long int> &objectsOut, int limit = -1);

    /**Count all objects that match a query.
     */
    int countObjects(const MDQuery &query);

    /** Find if the object with this id is present in the metadata
     */
    bool containsObject(long int objectId);

    /**Check if exists at least one object that match query.
     */
    bool containsObject(const MDQuery &query);
    /** @} */

    /** @name I/O functions
     * @{
     */

    /** Write rows data to disk. */
    void _writeRows(std::ostream &os);

    /** Write metadata to disk.
     * This will write the metadata content to disk.
     */
    void write(const FileName &outFile,const std::string & blockName="", WriteModeMetaData mode=OVERWRITE);

    /** Write metadata to out stream
     */
    void write(std::ostream &os, const std::string & blockName="",WriteModeMetaData mode=OVERWRITE);

    /** Append data lines to file.
     * This function can be used to add new data to
     * an existing metadata. Now should be used with
     * files with only one metadata, maybe can be extended later.
     * For now it will not check any compatibility beetween the
     * existent metadata and the new data to append.
     */
    void append(const FileName &outFile);

    /** Read data from file.
     */
    void read(const FileName &inFile, const std::vector<MDLabel> *desiredLabels = NULL, const std::string & blockName="", bool addStack=true);
    /** @} */

    /** @name Set Operations
     * @{
     */

    /** Aggregate metadata objects,
     * result in calling metadata object (except for aggregateSingle)
     * thisLabel label is used for aggregation, second. Valid operations are:
     *
     * MDL_AVG:  The avg function returns the average value of all  operationLabel within a group.
      The result of avg is always a floating point value as long as at there
      is at least one non-NULL input even if all inputs are integers.
       The result of avg is NULL if and only if there are no non-NULL inputs.

      AGGR_COUNT: The count function returns a count of the number of times that operationLabel is in a group.

      AGGR_MAX       The max aggregate function returns the maximum value of all values in the group.

      AGGR_MIN       The min aggregate function returns the minimum  value of all values in the group.

     AGGRL_SUM The total aggregate functions return sum of all values in the group.
     If there are no non-NULL input rows then returns 0.0.


     */

    void aggregate(const MetaData &mdIn, AggregateOperation op,
                   MDLabel aggregateLabel, MDLabel operateLabel, MDLabel resultLabel);
    void aggregate(const MetaData &mdIn, const std::vector<AggregateOperation> &ops,
                   MDLabel operateLabel, const std::vector<MDLabel> &resultLabels);
    double aggregateSingle(AggregateOperation op,
                           MDLabel aggregateLabel);
    /** Union of elements in two Metadatas, without duplicating.
     * Result in calling metadata object
     * union is a reserved word so I called this method unionDistinct
     */
    void unionDistinct(const MetaData &mdIn, const MDLabel label=MDL_OBJID);

    /** Union of all elements in two Metadata, duplicating common elements.
     * Result in calling metadata object
     * Repetion are allowed
     */
    void unionAll(const MetaData &mdIn);

    /** Merge of a Metadata.
     * This function reads another metadata and makes a union to this one
     */
    void merge(const FileName &fn);

    /** Intersects two Metadatas.
     * Result in "calling" Metadata
     */
    void intersection(const MetaData &mdIn, const MDLabel label);

    /** Substract two Metadatas.
     * Result in "calling" metadata
     */
    void subtraction(const MetaData &mdIn, const MDLabel label);

    /** Join two Metadatas
     * Result in "calling" metadata
     */
    void join(const MetaData &mdInLeft, const MetaData &mdInRight, const MDLabel label, JoinType type=LEFT);

    /** Basic operations on columns data.
     * Mainly perform replacements on string values and
     * basic algebraic operations on numerical ones.
     */
    void operate(const std::string &expression);

    /** Randomize a metadata.
     * MDin is input and the "randomized"
     * result will be in the "calling" Metadata.
    */
    void randomize(MetaData &MDin);

    /*
    * Sort a Metadata by a label.
    * Sort the content of MDin comparing
    * the label supplied, the result will
    * be in the "calling" MetaData.
    */
    void sort(MetaData &MDin, const MDLabel sortLabel);

    /** Split Metadata in several Metadatas.
     * The Metadata will be divided in 'n'
     * almost equally parts and the result will
     * be a vector of Metadatas. The "calling"
     * Metadata will no be modified.
     * @code
     *   // Divide the images metadata in 10 metadatas.
     *   std::vector<MetaData> imagesGroups;
     *   imageMD.split(10, imagesGroups);
     * @endcode
     */
    void split(int n, std::vector<MetaData> &results,
               const MDLabel sortLabel=MDL_OBJID);

    /** Take a part from MetaData.
     * This function is equivallent to divide
     * the input MetaData in n parts and take one.
     * The result will be in "calling" MetaData.
     */
    void selectSplitPart(const MetaData &mdIn,
                         int n, int part,
                         const MDLabel sortLabel=MDL_OBJID);

    /** Select some part from Metadata.
     * Select elements from input Metadata
     * at some starting position
     * if the numberOfObjects is -1, all objects
     * will be returned from startPosition to the end.
    */
    void selectPart (const MetaData &mdIn, long int startPosition, long int numberOfObjects,
                     const MDLabel sortLabel=MDL_OBJID);

    /** Makes filenames with absolute paths
    *
    */
    void makeAbsPath(const MDLabel label=MDL_IMAGE);


    /** @} */


    friend class MDSql;
};


/** Convert string to write mode metadata enum.
 *
 */
WriteModeMetaData metadataModeConvert (String mode);


/** @} */

#endif
