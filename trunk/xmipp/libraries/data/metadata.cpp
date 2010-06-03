/***************************************************************************
 * 
 * Authors:      J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

#include "metadata.h"

//-----Constructors and related functions ------------
void MetaData::_clear(bool onlyData)
{
    if (onlyData)
    {
        myMDSql->deleteObjects();
    }
    else
    {
        path.clear();
        comment.clear();
        fastStringSearch.clear();
        fastStringSearchLabel = MDL_UNDEFINED;

        activeLabels.clear();
        ignoreLabels.clear();
        isColumnFormat = true;
        inFile = FileName::FileName();
        activeObjId = -1;//no active object
        if (iterObjectsId != NULL)
        {
            delete iterObjectsId;
        }
        iterObjectsId = NULL;
        myMDSql->clearMd();
    }
}//close clear

void MetaData::clear()
{
    _clear(true);
}

void MetaData::init(const std::vector<MDLabel> *labelsVector)
{
    _clear();
    if (labelsVector != NULL)
        this->activeLabels = *labelsVector;
    //Create table in database
    myMDSql->createMd();
}//close init

void MetaData::copyInfo(const MetaData &md)
{
    if (this == &md) //not sense to copy same metadata
        return;
    this->setComment(md.getComment());
    this->setPath(md.getPath());
    this->isColumnFormat = md.isColumnFormat;
    this->inFile = md.inFile;
    this->fastStringSearchLabel = md.fastStringSearchLabel;
    this->activeLabels = md.activeLabels;
    this->ignoreLabels = md.ignoreLabels;

}//close copyInfo



void MetaData::copyMetadata(const MetaData &md)
{
    if (this == &md) //not sense to copy same metadata
        return;
    init(&(md.activeLabels));
    copyInfo(md);
    md.myMDSql->copyObjects(this);
}

bool MetaData::_setValue(long int objId, const MDValue &mdValueIn)
{
    if (objId == -1)
    {
        if (activeObjId != -1)
            objId = activeObjId;
        else
        {
            REPORT_ERROR(-1, "No active object, please provide objId for 'setValue'");
            exit(1);
        }
    }
    //add label if not exists, this is checked in addlabel
    addLabel(mdValueIn.label);
    //MDL::voidPtr2Value(label, valuePtr, mdValue);
    myMDSql->setObjectValue(objId, mdValueIn);
}

bool MetaData::_getValue(long int objId, MDValue &mdValueOut) const
{
    if (!containsLabel(mdValueOut.label))
        return false;

    if (objId == -1)
    {
        if (activeObjId != -1)
            objId = activeObjId;
        else
        {
            REPORT_ERROR(-1, "No active object, please provide objId for 'getValue'");
            exit(1);
        }
    }
    //MDValue mdValue;
    return myMDSql->getObjectValue(objId, mdValueOut);
    //MDL::value2VoidPtr(label, mdValue, valuePtrOut);
}

MetaData::MetaData()
{
    myMDSql = new MDSql(this);
    iterObjectsId = NULL;
    init(NULL);
}//close MetaData default Constructor

MetaData::MetaData(const std::vector<MDLabel> *labelsVector)
{
    myMDSql = new MDSql(this);
    iterObjectsId = NULL;
    init(labelsVector);
}//close MetaData default Constructor

MetaData::MetaData(const FileName &fileName, const std::vector<MDLabel> *labelsVector)
{
    myMDSql = new MDSql(this);
    iterObjectsId = NULL;
    init(labelsVector);
    //FIXME: what to do when labels vector is provided???
    read(fileName);
}//close MetaData from file Constructor

MetaData::MetaData(const MetaData &md)
{
    myMDSql = new MDSql(this);
    iterObjectsId = NULL;
    copyMetadata(md);
}//close MetaData copy Constructor

MetaData& MetaData::operator =(const MetaData &md)
{
    copyMetadata(md);
    return *this;
}//close metadata operator =

MetaData::~MetaData()
{
    _clear();
    delete myMDSql;
}//close MetaData Destructor

//-------- Getters and Setters ----------

bool MetaData::getColumnFormat() const
{
    return isColumnFormat;
}
/** Set to false for row format (parameter files)
 *  @ingroup GettersAndSetters
 *  set to true  for column format (this is the default) (docfiles)
 *
 */
void MetaData::setColumnFormat(bool column)
{
    isColumnFormat = column;
}
std::string MetaData::getPath()   const
{
    return path;
}

void MetaData::setPath(std::string newPath)
{
    const size_t length = 512;
    char _buffer[length];
    path = (newPath == "") ? std::string(getcwd(_buffer, length)) : newPath;
}

std::string MetaData::getComment() const
{
    return  comment;
}

void MetaData::setComment(const std::string newComment )
{
    comment = newComment;
}

FileName MetaData::getFilename()
{
    return inFile;
}

std::vector<MDLabel> MetaData::getActiveLabels() const
{
    return activeLabels;
}
long int  MetaData::getActiveObject()
{
    return activeObjId;
}

int MetaData::MaxStringLength(const MDLabel thisLabel) const
{
    if (!containsLabel(thisLabel))
        return -1;

    return myMDSql->columnMaxLength(thisLabel);
}

bool MetaData::setValueFromStr(const MDLabel label, const std::string &value, long int objectId)
{
    addLabel(label);

    if (objectId == -1)
    {
        if (activeObjId != -1)
            objectId = activeObjId;
        else
        {
            REPORT_ERROR(-1, "No active object, please provide objId for 'setValue'");
            exit(1);
        }
    }
    MDValue mdValue(label);
    mdValue.fromString(value);
    return myMDSql->setObjectValue(objectId, mdValue);
}

bool MetaData::getStrFromValue(const MDLabel label, std::string &strOut, long int objectId)
{
    MDValue mdValueOut(label);
    _getValue(objectId, mdValueOut);
    strOut = mdValueOut.toString();
}

bool MetaData::isEmpty() const
{
    return size() == 0;
}

long int MetaData::size() const
{
    std::vector<long int> objects;
    myMDSql->selectObjects(objects);

    return (long int)objects.size();
}

bool MetaData::containsLabel(const MDLabel label) const
{
    return vectorContainsLabel(activeLabels, label);
}

bool MetaData::addLabel(const MDLabel label)
{
    if (containsLabel(label))
        return false;
    activeLabels.push_back(label);
    myMDSql->addColumn(label);
    return true;
}

long int MetaData::addObject(long int objectId)
{
    activeObjId = myMDSql->addRow();
    return activeObjId;
}

void MetaData::importObject(const MetaData &md, const long int objId, bool doClear)
{
    md.myMDSql->copyObjects(this, new MDValueEqual(MDL_OBJID, objId));
}

void MetaData::importObjects(const MetaData &md, const std::vector<long int> &objectsToAdd, bool doClear)
{
    init(&(md.activeLabels));
    copyInfo(md);
    int size = objectsToAdd.size();
    for (int i = 0; i < size; i++)
        importObject(md, objectsToAdd[i]);
}

void MetaData::importObjects(const MetaData &md, const MDQuery &query, bool doClear)
{
    if (doClear)
    {
        //Copy all structure and info from the other metadata
        init(&(md.activeLabels));
        copyInfo(md);
    }
    else
    {
        //If not clear, ensure that the have the same labels
        for (int i = 0; i < md.activeLabels.size(); i++)
            addLabel(md.activeLabels[i]);
    }
    md.myMDSql->copyObjects(this, &query);
}

bool MetaData::removeObject(long int objectId)
{
    int removed = removeObjects(MDValueEqual(MDL_OBJID, objectId));
    return (removed > 0);
}

void MetaData::removeObjects(const std::vector<long int> &toRemove)
{
    int size = toRemove.size();
    for (int i = 0; i < size; i++)
        removeObject(toRemove[i]);
}

int MetaData::removeObjects(const MDQuery &query)
{
    int removed = myMDSql->deleteObjects(&query);
    firstObject(); //I prefer put active object to -1
    return removed;
}
int MetaData::removeObjects()
{
    int removed = myMDSql->deleteObjects();
    activeObjId = -1;
    return removed;
}

void MetaData::addIndex(MDLabel label)
{
    myMDSql->indexModify(label, true);
}

void MetaData::removeIndex(MDLabel label)
{
    myMDSql->indexModify(label, false);
}

long int MetaData::_iteratorBegin(const MDQuery *query)
{
    if (iterObjectsId == NULL)
        iterObjectsId = new std::vector<long int>;
    findObjects(*iterObjectsId);

    if (iterObjectsId->size() > 0)
    {
        iterIndex = 0;
        activeObjId = iterObjectsId->at(iterIndex);
    }
    else
    {
        activeObjId = iterIndex = -1;
    }
    return activeObjId;
}

long int MetaData::iteratorBegin()
{
    _iteratorBegin();
}

/**Same as previous but iterating over a subset of
 * objects
 */
long int MetaData::iteratorBegin(const MDQuery &query)
{
    _iteratorBegin(&query);
}

/** Check whether the iteration if finished */
bool MetaData::iteratorEnd()
{
    return iterIndex == -1;
}

/** Move to next object on iteration
 * return nextObject id
 */
long int MetaData::iteratorNext()
{
    if (iterObjectsId != NULL && iterIndex < iterObjectsId->size() - 1)
    {
        //The end not reached yet
        iterIndex++;
        activeObjId = iterObjectsId->at(iterIndex);
    }
    else
    {
        activeObjId = iterIndex = -1;
    }
    return activeObjId;
}

//----------Iteration functions -------------------
long int MetaData::firstObject()
{
    //std::vector<long int> objects = MDSql::selectObjects(this, 1);
    //return (objects.size() > 0) ? objects[0] : -1;
    activeObjId = myMDSql->firstRow();
    return activeObjId;
}

long int MetaData::lastObject()
{
    activeObjId = myMDSql->lastRow();
    return activeObjId;
}

long int MetaData::nextObject()
{
    if (activeObjId == -1)
        REPORT_ERROR(-55, "Couldn't call 'nextObject' when 'activeObject' is -1");
    activeObjId = myMDSql->nextRow(activeObjId);
    return activeObjId;
}

long int MetaData::previousObject()
{
    if (activeObjId == -1)
        REPORT_ERROR(-55, "Couldn't call 'previousObject' when 'activeObject' is -1");
    activeObjId = myMDSql->previousRow(activeObjId);
    return activeObjId;
}

long int MetaData::goToObject(long int objectId)
{
    if (containsObject(objectId))
        activeObjId = objectId;
    else
        REPORT_ERROR(-55, "Couldn't 'gotoObject', ID doesn't exist in metadata");
    return activeObjId;

}

//-------------Search functions-------------------
void MetaData::findObjects(std::vector<long int> &objectsOut, const MDQuery &query, int limit)
{
    objectsOut.clear();
    myMDSql->selectObjects(objectsOut, limit, &query);
}

void MetaData::findObjects(std::vector<long int> &objectsOut, int limit)
{
    objectsOut.clear();
    myMDSql->selectObjects(objectsOut, limit);
}

int MetaData::countObjects(MDQuery query)
{
    std::vector<long int> objects;
    findObjects(objects, query);
    return objects.size();
}

bool MetaData::containsObject(long int objectId)
{
    return containsObject(MDValueEqual(MDL_OBJID, objectId));
}

bool MetaData::containsObject(MDQuery query)
{
    std::vector<long int> objects;
    findObjects(objects, query, 1);
    return objects.size() > 0;
}

long int MetaData::gotoFirstObject(MDQuery query)
{
    std::vector<long int> objects;
    findObjects(objects, query, 1);

    activeObjId = objects.size() == 1 ? objects[0] : -1;
    return activeObjId;
}


//--------------IO functions -----------------------

void MetaData::write(const FileName &outFile)
{
    // Open file
    std::ofstream ofs(outFile.data(), std::ios_base::out);
    write(ofs);
}

void MetaData::write(std::ostream &os)
{
    os << "; XMIPP_3 * " << (isColumnFormat ? "column" : "row")
    << "_format * " << path << std::endl //write wich type of format (column or row) and the path
    << ";" << comment << std::endl; //write md comment in the 2nd comment line of header

    if (isColumnFormat)
    {
        //write md columns in 3rd comment line of the header
        os << "; ";
        for (int i = 0; i < activeLabels.size(); i++)
        {
            if (activeLabels.at(i) != MDL_COMMENT)
            {
                os.width(10);
                os << MDL::label2Str(activeLabels.at(i)) << " ";
            }
        }
        os << std::endl;
        //Write data
        FOR_ALL_OBJECTS_IN_METADATA(*this)
        {
            for (int i = 0; i < activeLabels.size(); i++)
            {
                if (activeLabels[i] != MDL_COMMENT)
                {
                    MDValue mdValue(activeLabels[i]);
                    os.width(10);
                    myMDSql->getObjectValue(activeObjId, mdValue);
                    mdValue.toStream(os);
                    os << " ";
                }
            }
            os << std::endl;
        }
    }
    else //rowFormat
    {
        // Get first object. In this case (row format) there is a single object
        int objId = firstObject();

        if (objId != -1)
        {
            int maxWidth=20;
            for (int i = 0; i < activeLabels.size(); i++)
            {
                if (activeLabels.at(i) != MDL_COMMENT)
                {
                    int w=MDL::label2Str(activeLabels.at(i)).length();
                    if (w>maxWidth)
                        maxWidth=w;
                }
            }

            for (int i = 0; i < activeLabels.size(); i++)
            {
                if (activeLabels[i] != MDL_COMMENT)
                {
                    MDValue mdValue(activeLabels[i]);
                    os.width(maxWidth + 1);
                    os << MDL::label2Str(activeLabels.at(i)) << " ";
                    myMDSql->getObjectValue(objId, mdValue);
                    mdValue.toStream(os);
                    os << std::endl;
                }
            }
        }

    }
}//write

/** This function will read the posible columns from the file
 * and mark as MDL_UNDEFINED those who aren't valid labels
 * or those who appears in the IgnoreLabels vector
 * also set the activeLabels
 */
void MetaData::_readColumns(std::istream& is, std::vector<MDValue>& columnValues,
                            std::vector<MDLabel>* desiredLabels)
{
    std::string token;
    MDLabel label;

    while (is >> token)
        if (token.find('(') == std::string::npos)
        {
            //label is not reconized, the MDValue will be created
            //with MDL_UNDEFINED, wich will be ignored while reading data
            label = MDL::str2Label(token);
            if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
                label = MDL_UNDEFINED; //ignore if not present in desiredLabels
            columnValues.push_back(MDValue(label));
            if (label != MDL_UNDEFINED)
                addLabel(label);
        }
}

/** This function will be used to parse the rows data
 * having read the columns labels before and setting wich are desired
 * the useCommentAsImage is for compatibility with old DocFile format
 * where the image were in comments
 */
void MetaData::_readRows(std::istream& is, std::vector<MDValue>& columnValues, bool useCommentAsImage)
{
    std::string line = "";
    int objId;
    while (!is.eof() && !is.fail())
    {
        //Move until the ';' or the first alphanumeric character
        while (is.peek() != ';' && !isalnum(is.peek()) && !is.eof())
            is.ignore(1);
        if (!is.eof())

            if (is.peek() == ';')//is a comment
            {
                is.ignore(1); //ignore the ';'
                getline(is, line);
                trim(line);
            }
            else if (isalnum(is.peek()))
            {
                objId = addObject();
                if (line != "")
                {
                    if (!useCommentAsImage)
                        setValue(MDL_COMMENT, line);
                    else
                        setValue(MDL_IMAGE, line);
                }
                for (int i = 0; i < columnValues.size(); i++)
                {
                    is >> columnValues[i];
                    if (is.fail())
                    {
                        REPORT_ERROR(-44, "MetaData::read2: Error parsing data column");
                    }
                    else
                        if (columnValues[i].label != MDL_UNDEFINED)
                            _setValue(activeObjId, columnValues[i]);
                }
            }

    }
}

/**This function will read the md data if is in row format */
void MetaData::_readRowFormat(std::istream& is)
{
    std::string line, token;
    MDLabel label;

    long int objectID = addObject();

    // Read data and fill structures accordingly
    while (getline(is, line, '\n'))
    {
        if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
            continue;

        // Parse labels
        std::stringstream os(line);

        os >> token;
        label = MDL::str2Label(token);
        MDValue value(label);
        os >> value;
        if (label != MDL_UNDEFINED)
            _setValue(objectID, value);
    }
}

void MetaData::read(const FileName &filename, std::vector<MDLabel> *desiredLabels)
{
    int pos;
    std::ifstream is(filename.data(), std::ios_base::in);
    std::stringstream ss;
    std::string line, token;
    std::vector<MDValue> columnValues;
    isColumnFormat = true;
    bool useCommentAsImage = false;

    getline(is, line); //get first line to identify the type of file

    if (is.fail())
    {
        REPORT_ERROR(-44, (std::string) "MetaData::read: File " + filename + " does not exists" );
    }

    _clear();
    myMDSql->createMd();

    if (pos = line.find("XMIPP_3 *") != std::string::npos)
    {
        //We have a new XMIPP MetaData format here, parse header
        is.seekg(0, std::ios::beg);
        is.ignore(256, '*') >> token; //Ignore all until first '*' and get md format in token
        isColumnFormat = token != "row_format";
        is.ignore(256, '*') >> token;
        if (token != ";") //there is path, need to ignore ';' of the next line
        {
            setPath(token);
            is.ignore(256, ';');
        }
        getline(is, line);
        setComment(line);
        if (isColumnFormat)
        {
            is.ignore(256, ';'); //ignore ';' to start parsing column labels
            getline(is, line);
            ss.str(line);
            //Read column labels
            _readColumns(ss, columnValues, desiredLabels);
        }
    }
    else if (pos = line.find("Headerinfo columns:") != std::string::npos)
    {
        is.seekg(0, std::ios::beg);
        //This looks like an old DocFile, parse header
        std::cerr << "WARNING: ** You are using an old file format (DOCFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        is.ignore(256, ':'); //ignore all until ':' to start parsing column labels
        getline(is, line);
        ss.str(line);
        columnValues.resize(2, MDValue(MDL_UNDEFINED)); //start with 2 undefined to avoid 2 first columns of old format
        addLabel(MDL_IMAGE);
        _readColumns(ss, columnValues, desiredLabels);
        useCommentAsImage = true;
    }
    else
    {
        std::cerr << "WARNING: ** You are using an old file format (SELFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        //I will assume that is an old SelFile, so only need to add two columns
        columnValues.push_back(MDValue(MDL_IMAGE));//addLabel(MDL_IMAGE);
        columnValues.push_back(MDValue(MDL_ENABLED));//addLabel(MDL_ENABLED);
    }

    if (isColumnFormat)
        _readRows(is, columnValues, useCommentAsImage);
    else
        _readRowFormat(is);

}

void MetaData::merge(const FileName &fn)
{
    MetaData md;
    md.read(fn);
    unionAll(md);
}

void MetaData::aggregate(const MetaData &mdIn, AggregateOperation op,
                         MDLabel aggregateLabel, MDLabel operateLabel, MDLabel resultLabel)

{
    std::vector<MDLabel> labels(2);
    labels[0] = aggregateLabel;
    labels[1] = resultLabel;
    init(&labels);
    std::vector<AggregateOperation> ops(1);
    ops[0] = op;
    mdIn.myMDSql->aggregateMd(this, ops, operateLabel);
}

void MetaData::aggregate(const MetaData &mdIn, const std::vector<AggregateOperation> &ops,
                         MDLabel operateLabel, const std::vector<MDLabel> &resultLabels)
{
    if (resultLabels.size() - ops.size() != 1)
        REPORT_ERROR(-55, "Labels vectors should contain one element more than operations");
    init(&resultLabels);
    mdIn.myMDSql->aggregateMd(this, ops, operateLabel);
}


//-------------Set Operations ----------------------
void MetaData::_setOperates(const MetaData &mdIn, const MDLabel label, int operation)
{
    if (this == &mdIn) //not sense to operate on same metadata
        return;
    //Add labels to be sure are present
    for (int i = 0; i < mdIn.activeLabels.size(); i++)
        addLabel(mdIn.activeLabels[i]);

    mdIn.myMDSql->setOperate(this, label, operation);
}

void MetaData::unionDistinct(const MetaData &mdIn, const MDLabel label)
{
    _setOperates(mdIn, label, 1);
}

void MetaData::unionAll(const MetaData &mdIn)
{
    if (this == &mdIn) //not sense to copy same metadata
        return;
    //Add labels to be sure are present
    for (int i = 0; i < mdIn.activeLabels.size(); i++)
        addLabel(mdIn.activeLabels[i]);
    mdIn.myMDSql->copyObjects(this);
}


void MetaData::intersection(const MetaData &mdIn, const MDLabel label)
{
    _setOperates(mdIn, label, 3);
}
void MetaData::substraction(const MetaData &mdIn, const MDLabel label)
{
    _setOperates(mdIn, label, 2);
}

void MetaData::randomize(MetaData &MDin)
{
    std::vector<long int> objects;
    MDin.myMDSql->selectObjects(objects);
    std::random_shuffle(objects.begin(), objects.end());
    importObjects(MDin, objects);
}

void MetaData::sort(MetaData &MDin, const MDLabel sortLabel)
{
    init(&(MDin.activeLabels));
    copyInfo(MDin);
    MDin.myMDSql->copyObjects(this, NULL, sortLabel);
}

void MetaData::split(int n, std::vector<MetaData> &results, const MDLabel sortLabel)
{
    long int mdSize = size();
    if (n > mdSize)
        REPORT_ERROR(-55, "split: Couldn't split a metadata in more parts than its size");

    results.clear();
    results.resize(n);
    for (int i = 0; i < n; i++)
    {
        MetaData &md = results.at(i);
        md._selectSplitPart(*this, n, i, mdSize, sortLabel);
    }
}

void MetaData::_selectSplitPart(const MetaData &mdIn,
                                int n, int part, long int mdSize,
                                const MDLabel sortLabel)
{
    int first, last, n_images;
    n_images = divide_equally(mdSize, n, part, first, last);
    init(&(mdIn.activeLabels));
    copyInfo(mdIn);
    mdIn,myMDSql->copyObjects(this, NULL, sortLabel, n_images, first);
}

void MetaData::selectSplitPart(const MetaData &mdIn, int n, int part, const MDLabel sortLabel)
{
    long int mdSize = mdIn.size();
    if (n > mdSize)
        REPORT_ERROR(-55, "selectSplitPart: Couldn't split a metadata in more parts than its size");
    if (part < 0 || part >= n)
        REPORT_ERROR(-55, "selectSplitPart: 'part' should be between 0 and n-1");
    _selectSplitPart(mdIn, n, part, mdSize, sortLabel);

}

void MetaData::selectPart (const MetaData &mdIn, long int startPosition, long int numberOfObjects,
                           const MDLabel sortLabel)
{
    long int mdSize = mdIn.size();
    if (startPosition < 0 || startPosition >= mdSize)
        REPORT_ERROR(-55, "selectPart: 'startPosition' should be between 0 and size()-1");
    init(&(mdIn.activeLabels));
    copyInfo(mdIn);
    mdIn.myMDSql->copyObjects(this, NULL, sortLabel, numberOfObjects, startPosition);
}

