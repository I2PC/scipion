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
#include "metadata_container.h"
#include "image.h"
#include <time.h>
#include <stdio.h>
#include "cppsqlite3.h"

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
class MetaData
{
    std::map<long int, MetaDataContainer *> objects; ///< Effectively stores all metadata

    // Used by firstObject, nextObject and lastObject to keep a pointer
    // to the "active" object. This way when you call setValue without
    // an objectID, this one is chosen
    std::map<long int, MetaDataContainer *>::iterator objectsIterator;

    // Allows a fast search for pairs where the value is
    // a string, i.e. looking for filenames which is quite
    // usual
    std::map<std::string, long int> fastStringSearch;
    MetaDataLabel fastStringSearchLabel;

    std::vector<long int> randomOrderedObjects;

    std::string path; ///< A parameter stored on MetaData Files
    std::string comment; ///< A general comment for the MetaData file

    MetaDataContainer * getObject(long int objectID = -1);
    void read(std::ifstream *infile, std::vector<MetaDataLabel> * labelsVector);

    bool isColumnFormat; ///< Format for the file, column or row formatted

    /**Input file name
     * Where does this MetaData come from/go to be stored?
     */
    FileName inFile;

public:

    /// @defgroup MetaDataConstructors Constructors for MetaData objects
    /// @ingroup MetaDataClass

    /** Empty Constructor.
     * @ingroup MetaDataConstructors
     *
     * The MetaData is created with no data stored on it. You can fill in it programmatically
     * or by a later reading from a MetaData file or old Xmipp formatted type.
     */
    MetaData();

    /** From File Constructor.
     * @ingroup MetaDataConstructors
     *
     * The MetaData is created and data is read from provided fileName. Optionally, a vector
     * of labels can be provided to read just those required labels
     */
    MetaData(FileName fileName, std::vector<MetaDataLabel> * labelsVector =
                 NULL);

    /** Copy constructor
     * @ingroup MetaDataConstructors
     *
     * Created a new metadata by copying all data from an existing MetaData object.
     */
    MetaData(MetaData & c);

    /** Assignment operator
     * @ingroup MetaDataConstructors
     *
     * Copies MetaData from an existing MetaData object.
     */
    MetaData& operator =(MetaData &MD);

    /** union of  metadata objects, result in calling metadata object
     * union is a reserved word so I called this method union_
     */
    void union_(MetaData &MD);

    /** intersects two metadata objects, result in "calling" metadata
     */
    void intersection(MetaData & minuend, MetaData & subtrahend,
                      MetaDataLabel thisLabel);

    /** substract two metadata objects, result in "calling" metadata
     */
    void substraction(MetaData & minuend, MetaData & subtrahend,
                      MetaDataLabel thisLabel);

    /** Destructor
     * @ingroup MetaDataConstructors
     *
     * Frees all used memory and destroys object.
     */
    ~MetaData();

    /** Set to false for row format (parameter files)
     *  set to true  for column format (this is the default) (docfiles)
     *
     */
    void setColumnFormat(bool column);
    /**Get column format info.
     *
     */
    bool getColumnFormat()
    {
        return isColumnFormat;
    }
    // Set a new pair/value for an specified object. If no objectID is given, that
    // pointed by the class iterator is used
    bool setValue(MetaDataLabel name, double value, long int objectID = -1);
    bool setValue(MetaDataLabel name, int value, long int objectID = -1);
    bool setValue(MetaDataLabel name, bool value, long int objectID = -1);
    bool setValue(MetaDataLabel name, const std::string &value,
                  long int objectID = -1);
    bool setValue(MetaDataLabel name, const std::vector<double> &value,
                  long int objectID = -1);
    bool setValue(const std::string &name, const std::string &value,
                  long int objectID = -1);

    void readOldSelFile(std::ifstream *infile);
    void readOldDocFile(std::ifstream *infile,
                        std::vector<MetaDataLabel> * labelsVector);
    void read(FileName infile,
              std::vector<MetaDataLabel> * labelsVector = NULL);

    /**convert metafile to table in database
     *
     */
    void toDataBase(const FileName & DBname,
                    const std::string & tableName = "",
                    std::vector<MetaDataLabel> * labelsVector = NULL);
    /**Convert table from database in metadata
     *
     */
    void fromDataBase(const FileName & DBname, const std::string & tableName,
                      std::vector<MetaDataLabel> * labelsVector = NULL);
    std::vector<long int> & getRandomOrderedObjects();

    /** What labels have been read from a docfile/metadata file
     *   and/or will be stored on a new metadata file when "save" is
     *   called
     **/
    std::vector<MetaDataLabel> activeLabels;

    /** When reading a column formated file, if a label is found that
     *   does not exists as a MetaDataLabel, it is ignored. For further
     *   file processing, such columns must be ignored and this structure
     *   allows to do that
     **/
    std::vector<unsigned int> ignoreLabels;

    /** Adds a new, empty object to the objects map. If objectID == -1
     *   the new ID will be that for the last object inserted + 1, else
     *   the given objectID is used. If there is already an object whose
     *   objectID == input objectID, just removes it and creates an empty
     *   one
     **/
    long int addObject(long int objectID = -1);

    // Possible error codes for the map
    enum errors
    {
        NO_OBJECTS_STORED = -1, // NOTE: Do not change this value (-1)
        NO_MORE_OBJECTS = -2,
        NO_OBJECT_FOUND = -3
    };

    long int firstObject();
    long int nextObject();
    long int lastObject();
    long int goToObject(long int objectID);

    void write(const std::string &fileName);

    bool isEmpty();

    void clear();

    void writeValueToString(std::string & result,
                            const std::string & inputLabel);

    long int fastSearch(MetaDataLabel name, std::string value, bool recompute =
                            false);

    /**Create new metadata with selected objectId
     *
     */
    void fillMetaData(MetaData &base, std::vector<long int> objectsToAdd);

    void combine(MetaData & other, MetaDataLabel thisLabel = MDL_UNDEFINED);
    void combineWithFiles(MetaDataLabel thisLabel);

    // Removes the collection of objects whose pair label/value is given
    // NOTE: The iterator will point to the first object after any of these
    // operations.
    void removeObjects(std::vector<long int> &toRemove);

    // Returns true if the object was removed or false if
    // the object did not exist
    bool removeObject(long int objectID);

    /**Set Path
     * will appear in first line
     */
    void setPath(std::string newPath = "");

    /**Set Header Comment
     * will appear in second line
     */
    void setComment(std::string Comment = "");

    std::string getPath();
    std::string getComment();

    size_t size(void)
    {
        return objects.size();
    }
    /** Return metafile filename
     *
     */

    FileName getFilename()
    {
        return (inFile);
    }



    /*Detect if there is at least one entry with the given label-value pair
     * This can be much faster than 'countObjects' as it stops iterating once the first
     * object has been found.
     */
    bool detectObjects(MetaDataLabel name, int value);

    /**Add object with metadata label name in the range given by minvalue and maxvalue
     * This template function may be accessed from swig
     *
     * */
    template<class T>
    friend void addObjectsInRange(MetaData &MDin, MetaData &MDout,
                                  MetaDataLabel name, T minValue, T maxValue)
    {
        MDout.fillMetaData(MDin, MDin.findObjectsInRange(name, minValue, maxValue));
    }

    template<class T>
    std::vector<long int> findObjectsInRange(MetaDataLabel name, T minValue, T maxValue)
    {
        std::vector<long int> result;

        // Traverse all the structure looking for objects
        // that satisfy search criteria

        std::map<long int, MetaDataContainer *>::iterator It;

        MetaDataContainer * aux;

        for (It = objects.begin(); It != objects.end(); It++)
        {
            aux = It->second;

            if (aux->valueExists(name))
            {
                double value;
                aux->getValue(name, value);

                if (value >= minValue && value <= maxValue)
                {
                    result.push_back(It->first);
                }
            }
        }

        return result;
    }

    template<class T>
    std::vector<long int> findObjects(MetaDataLabel name,
                                      T value)
    {
        std::vector<long int> result;

        // Traverse all the structure looking for objects
        // that satisfy search criteria

        std::map<long int, MetaDataContainer *>::iterator It;

        MetaDataContainer * aux;

        for (It = objects.begin(); It != objects.end(); It++)
        {
            aux = It->second;

            if (aux->pairExists(name, value))
                result.push_back(It->first);
        }

        return result;
    }

    template<class T>
    void removeObjects(MetaDataLabel name, const T &value)
    {
        std::vector<long int> toRemove = findObjects(name, value);
        removeObjects(toRemove);
    }

    /**Count number of objects that satisfy a given label,entry pair
     */
    template<class T>
    long int countObjects(MetaDataLabel name, const T &value)
    {
        return ((findObjects(name, value)).size());
    }

    template<class T>
    bool getValue(MetaDataLabel name, T &value,
                  long int objectID = -1)
    {
        MetaDataContainer * aux = getObject(objectID);
        if (!aux->valueExists(name))
        {
            return false;
        }
        else
        {
            aux->getValue(name, value);
        }

        return true;
    }


};

/** Compute images metadata estatistics
 * This use to be part of Metadata but should not
 */

void get_statistics(MetaData MT, Image& _ave, Image& _sd, double& _min,
                    double& _max, bool apply_geo);

/** Get image size
 *
 */
void ImgSize(MetaData MD, int &Xdim, int &Ydim, int &Zdim, int &Ndim);

/** For all objects.
 @code
 FOR_ALL_OBJECTS_IN_METADATA(metadata) {
 double rot; DF.getValue( MDL_ANGLEROT, rot);
 }
 @endcode
 */
#define FOR_ALL_OBJECTS_IN_METADATA(kkkk_metadata) \
        for(long int kkkk = kkkk_metadata.firstObject(); \
             kkkk != MetaData::NO_MORE_OBJECTS; \
             kkkk=kkkk_metadata.nextObject())
#endif

