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
#include "image.h"

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
        inFile = FileName();
        activeObjId = -1;//no active object
        delete iterObjectsId;
        iterObjectsId = NULL;
        myMDSql->clearMd();
    }
}//close clear

void MetaData::clear()
{
    //_clear(true);
    init();
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
    if (!md.activeLabels.empty())
    {
        md.myMDSql->copyObjects(this);
        firstObject(); //set first object as active
    }
    else
    {
        int n = md.size();
        for (int i = 0; i < n; i++)
            addObject();
    }
}

bool MetaData::setValue(const MDObject &mdValueIn, long int objId)
{
    if (objId == -1)
    {
        if (activeObjId != -1)
            objId = activeObjId;
        else
        {
            REPORT_ERROR(ERR_MD_NOACTIVE, "setValue: please provide objId other than -1");
            exit(1);
        }
    }
    //add label if not exists, this is checked in addlabel
    addLabel(mdValueIn.label);
    myMDSql->setObjectValue(objId, mdValueIn);
}

bool MetaData::setValueCol(const MDObject &mdValueIn)
{
    //add label if not exists, this is checked in addlabel
    addLabel(mdValueIn.label);
    myMDSql->setObjectValue(mdValueIn);
}

bool MetaData::getValue(MDObject &mdValueOut, long int objId) const
{
    if (!containsLabel(mdValueOut.label))
        return false;

    if (objId == -1)
    {
        if (activeObjId != -1)
            objId = activeObjId;
        else
        {
            REPORT_ERROR(ERR_MD_NOACTIVE, "getValue: please provide objId other than -1");
            exit(1);
        }
    }
    //MDValue mdValue;
    return myMDSql->getObjectValue(objId, mdValueOut);
}

bool MetaData::getRow(MDRow &row, long int objId)
{
    row.clear();
    MDObject * obj;
    for (std::vector<MDLabel>::const_iterator it = activeLabels.begin(); it != activeLabels.end(); ++it)
    {
        obj = new MDObject(*it);
        if (!getValue(*obj, objId))
            return false;
        row.push_back(obj);
    }
    return true;
}

void MetaData::setRow(const MDRow &row)
{
    int imax=row.size();
    for (int i=0; i<imax; ++i)
        setValue(*row[i]);
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

MetaData::MetaData(const FileName &fileName, const std::vector<MDLabel> *desiredLabels)
{
    myMDSql = new MDSql(this);
    iterObjectsId = NULL;
    init(desiredLabels);
    read(fileName, desiredLabels);
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

void MetaData::setComment(const std::string &newComment)
{
    comment = newComment;
}

FileName MetaData::getFilename() const
{
    return inFile;
}

void MetaData::setFilename(const FileName _fileName)
{
    inFile=_fileName;
}

std::vector<MDLabel> MetaData::getActiveLabels() const
{
    return activeLabels;
}

std::vector<MDLabel>* MetaData::geActiveLabelsAddress() const
{
    return (std::vector<MDLabel>*) (&activeLabels);
}

long int  MetaData::getActiveObject() const
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
            REPORT_ERROR(ERR_MD_NOACTIVE, "setValue: please provide objId other than -1");
            exit(1);
        }
    }
    MDObject mdValue(label);
    mdValue.fromString(value);
    return myMDSql->setObjectValue(objectId, mdValue);
}

bool MetaData::getStrFromValue(const MDLabel label, std::string &strOut, long int objectId)
{
    MDObject mdValueOut(label);
    if (!getValue(mdValueOut, objectId))
        return false;
    strOut = mdValueOut.toString();
    return true;
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

bool MetaData::addLabel(const MDLabel label, int pos)
{
    if (containsLabel(label))
        return false;
    if (pos < 0 || pos >= activeLabels.size())
        activeLabels.push_back(label);
    else
        activeLabels.insert(activeLabels.begin() + pos, label);
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
    md.myMDSql->copyObjects(this, new MDValueEQ(MDL_OBJID, objId));
}

void MetaData::importObjects(const MetaData &md, const std::vector<long int> &objectsToAdd, bool doClear)
{
    init(&(md.activeLabels));
    copyInfo(md);
    int size = objectsToAdd.size();
    for (int i = 0; i < size; i++)
        importObject(md, objectsToAdd[i]);
    firstObject();
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
    firstObject();
}

bool MetaData::removeObject(long int objectId)
{
    int removed = removeObjects(MDValueEQ(MDL_OBJID, objectId));
    firstObject();
    return (removed > 0);
}

void MetaData::removeObjects(const std::vector<long int> &toRemove)
{
    int size = toRemove.size();
    for (int i = 0; i < size; i++)
        removeObject(toRemove[i]);
    firstObject();
}

int MetaData::removeObjects(const MDQuery &query)
{
    int removed = myMDSql->deleteObjects(&query);
    firstObject();
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
    return _iteratorBegin();
}

/**Same as previous but iterating over a subset of
 * objects
 */
long int MetaData::iteratorBegin(const MDQuery &query)
{
    return _iteratorBegin(&query);
}

/** Check whether the iteration if finished */
bool MetaData::iteratorEnd() const
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
        REPORT_ERROR(ERR_MD_NOACTIVE, "nextObject: Couldn't perform this operation when 'activeObject' is -1");
    activeObjId = myMDSql->nextRow(activeObjId);
    return activeObjId;
}

long int MetaData::previousObject()
{
    if (activeObjId == -1)
        REPORT_ERROR(ERR_MD_NOACTIVE, "previousObject: Couldn't perform this operation when 'activeObject' is -1");
    activeObjId = myMDSql->previousRow(activeObjId);
    return activeObjId;
}

long int MetaData::goToObject(long int objectId)
{
    if (containsObject(objectId))
        activeObjId = objectId;
    else
        REPORT_ERROR(ERR_MD_NOOBJ, (std::string)"gotoObject: object with ID " + integerToString(objectId) + " doesn't exist in metadata");
    return activeObjId;

}

//-------------Search functions-------------------
void MetaData::findObjects(std::vector<long int> &objectsOut, const MDQuery &query)
{
    objectsOut.clear();
    myMDSql->selectObjects(objectsOut, &query);
}

void MetaData::findObjects(std::vector<long int> &objectsOut, int limit)
{
    objectsOut.clear();
    MDQuery query(limit);
    myMDSql->selectObjects(objectsOut, &query);
}

int MetaData::countObjects(const MDQuery &query)
{
    std::vector<long int> objects;
    findObjects(objects, query);
    return objects.size();
}

bool MetaData::containsObject(long int objectId)
{
    return containsObject(MDValueEQ(MDL_OBJID, objectId));
}

bool MetaData::containsObject(const MDQuery &query)
{
    std::vector<long int> objects;
    findObjects(objects, query);
    return objects.size() > 0;
}

long int MetaData::gotoFirstObject(const MDQuery &query)
{
    std::vector<long int> objects;
    findObjects(objects, query);

    activeObjId = objects.size() == 1 ? objects[0] : -1;
    return activeObjId;
}

//--------------IO functions -----------------------
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

void MetaData::write(const FileName &_outFile, WriteModeMetaData mode)
{
    std::string blockName;
    FileName outFile;

    blockName=_outFile.getBlockName();
    outFile = _outFile.removeBlockName();
    _write(outFile, blockName, mode);

}

void MetaData::_write(const FileName &outFile,const std::string &blockName, WriteModeMetaData mode)
{
    struct stat file_status;
    int fd;
    char *map;

    //check if file exists or not block name has been given
    //in our format no two identical data_xxx strings may exists
    if(mode==OVERWRITE)
        ;
    else if (blockName=="")
        mode=OVERWRITE;
    else if(!exists(outFile))
        mode=OVERWRITE;
    else
    {
        //does blockname exists?
        //remove it from file in this case
        // get length of file:
        if(stat(outFile.data(), &file_status) != 0)
            REPORT_ERROR(ERR_IO_NOPATH,"Metadata:write can not get filesize for file "+outFile);
        size_t size = file_status.st_size;
        if(size!=0)//size=0 for /dev/stderr
        {
            fd = open(outFile.data(),  O_RDWR, S_IREAD | S_IWRITE);
            if (fd == -1)
                REPORT_ERROR(ERR_IO_NOPATH,"Metadata:write can not read file named "+outFile);

            map = (char *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (map == MAP_FAILED)
                REPORT_ERROR(ERR_MEM_BADREQUEST,"Metadata:write can not map memory ");

            // Is this a START formatted FILE

            if(strncmp(map,"# XMIPP_STAR",12)!=0)
            {
                mode=OVERWRITE;
            }
            else
            {
                //block name
                std::string _szBlockName = (std::string)("data_") + blockName;

                //search for the string
                char * target, * target2;
                target = (char *) strstr(map,_szBlockName.data());
                if(target!=NULL)
                {
                    target2 = (char *) strstr(target+1,"data_");

                    if (target2==NULL)//truncate file at target
                        ftruncate(fd, target - map);
                    else//copy file from target2 to target and truncate
                    {
                        memmove(target,target2, (map + size) - target2);
                        ftruncate(fd, (target-map) + ((map + size) - target2) );
                    }
                }
            }
            if (munmap(map, size) == -1)
            {
                REPORT_ERROR(ERR_MEM_NOTDEALLOC,"metadata:write, Can not unmap memory");
            }
            close(fd);
        }
        else
            mode=OVERWRITE;

    }

    std::ios_base::openmode openMode;
    if(mode==OVERWRITE)
        openMode = std::ios_base::out;
    else if(mode=APPEND)
        openMode = std::ios_base::app;
    std::ofstream ofs(outFile.data(), openMode);

    write(ofs, blockName, mode);
    ofs.close();

}

void MetaData::append(const FileName &outFile)
{
    std::ofstream ofs(outFile.data(), std::ios_base::app);
    _writeRows(ofs);
    ofs.close();
}

void MetaData::_writeRows(std::ostream &os)
{

    FOR_ALL_OBJECTS_IN_METADATA(*this)
    {
        for (int i = 0; i < activeLabels.size(); i++)
        {
            if (activeLabels[i] != MDL_COMMENT)
            {
                MDObject mdValue(activeLabels[i]);
                os.width(1);
                myMDSql->getObjectValue(activeObjId, mdValue);
                mdValue.toStream(os);
                os << " ";
            }
        }
        os << std::endl;
    }
}

void MetaData::write(std::ostream &os,const std::string &blockName, WriteModeMetaData mode )
{
    if(mode==OVERWRITE)
        os << "# XMIPP_STAR_1 * "// << (isColumnFormat ? "column" : "row")
        << std::endl //write wich type of format (column or row) and the path
        << "# " << comment << std::endl;     //write md comment in the 2nd comment line of header
    //write data block
    std::string _szBlockName = (std::string)("data_") + blockName;

    if (isColumnFormat)
    {
        //write md columns in 3rd comment line of the header
        os << _szBlockName << std::endl;
        os << "loop_" << std::endl;

        for (int i = 0; i < activeLabels.size(); i++)
        {
            if (activeLabels.at(i) != MDL_COMMENT)
            {
                os << " _" << MDL::label2Str(activeLabels.at(i)) << std::endl;
            }
        }

        _writeRows(os);
        //Put the activeObject to the first, if exists
        firstObject();
    }
    else //rowFormat
    {
        os << _szBlockName << std::endl;

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
                    MDObject mdValue(activeLabels[i]);
                    os << " _" << MDL::label2Str(activeLabels.at(i)) << " ";
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
 * also set the activeLabels (for OLD doc files)
 */
void MetaData::_readColumns(std::istream& is, MDRow & columnValues,
                            const std::vector<MDLabel>* desiredLabels)
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
            columnValues.push_back(new MDObject(label));
            if (label != MDL_UNDEFINED)
                addLabel(label);

        }
}

/* Helper function to parse an MDObject and set its value.
 * The parsing will be from an input stream(istream)
 * and if parsing fails, an error will be raised
 */
void MetaData::_parseObject(std::istream &is, MDObject &object)
{
    is >> object;
    if (is.fail())
    {
        REPORT_ERROR(ERR_MD_BADLABEL, (std::string)"read: Error parsing data column, expecting " + MDL::label2Str(object.label));
    }
    else
        if (object.label != MDL_UNDEFINED)
            setValue(object, activeObjId);
}//end of function parseObject

#define END_OF_LINE() ((char*) memchr (pchStart, '\n', pEnd-pchStart+1))
/** This function will read the posible columns from the file
 * and mark as MDL_UNDEFINED those who aren't valid labels
 * or those who appears in the IgnoreLabels vector
 * also set the activeLabels (for new STAR files)
 */
char * MetaData::_readColumnsStar(char * pStart,
                                  char * pEnd,
                                  MDRow & columnValues,
                                  const std::vector<MDLabel>* desiredLabels)
{
    char * pchStart, *pchEnd;
    pchStart =pStart;
    pchStart = END_OF_LINE() + 1;//skip _loop line
    MDLabel label;
    while(1)
    {
        //Next char should be _ or we are done
        //if((pchStart+1)[0]!='_')
        //{
        if( isspace((pchStart)[0]))//trim spaces and newlines at the beginning
        {
            ++pchStart;
            continue;
        }
        if((pchStart)[0]=='#')//this is a comment
        {
            pchStart = END_OF_LINE() + 1;
            continue;
        }
        if( (pchStart)[0] !='_')
        {
            break;
        }

        ++pchStart;//skip '_'

        if(pchStart < pEnd /*&& pchStart !=NULL*/)
        {
            std::stringstream ss;
            pchEnd = END_OF_LINE();
            String s(pchStart, pchEnd-pchStart);//get string line
            ss.str(s);//set the string of the stream
            ss >> s; //get the first token, the label
            //will fail is string label has spaces, tabs or newlines
            label = MDL::str2Label(s);
            if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
            {
                label = MDL_UNDEFINED; //ignore if not present in desiredLabels
            }
            else
                addLabel(label);
            MDObject * _mdObject = new MDObject(label);
            columnValues.push_back(_mdObject);//add the value here with a char
            if(!isColumnFormat)
                _parseObject(ss, *_mdObject);
            pchStart = pchEnd + 1;//go to next line character
        }
        else
            break;
    }
    return pchStart;
}
#ifdef NEVERDEFINED
/** This function will read the possible columns and values from the file
 * in ROW format
 * and mark as MDL_UNDEFINED those who aren't valid labels
 * or those who appears in the IgnoreLabels vector
 * also set the activeLabels (for new STAR files)
 */

char * MetaData::_readRowFormatStar(char * pStart,
                                    char * pEnd,
                                    MDRow & columnValues,
                                    std::vector<MDLabel>* desiredLabels)
{
    //search for "_" and "\n"
    char * pchStart, *pchEnd;
    pchStart =pStart;
    pchEnd=pEnd;
    pchStart = (char*) memchr (pchStart, '\n', pchEnd-pchStart);//after nextline
    MDLabel label;
    char * szAux;
    szAux = (char *) calloc(1024,1);
    char toks[] = "\t \n"; // tabs, spaces and newlines
    addObject();
    while(1)
    {
        //Is next char '_'
        if((pchStart+1)[0]!='_')
        {
            if((pchStart+1)[0]=='#')//this is a comment
            {
                pchStart = (char*) memchr (pchStart+1, '\n', pchEnd-pchStart-1);//after nextline
                continue;
            }
            else if((pchStart+1)[0]==' '||(pchStart+1)[0]=='\t')//trim spaces at the beginning
            {
                ++pchStart;
                continue;
            }
            else
            {
                //pchEnd = pchStart;
                break;
            }
        }

        pchStart = (char*) memchr (pchStart, '_' , pEnd-pchStart) ;//begining label
        if(pchStart==NULL)
            break;
        ++pchStart;
        //pchEnd = (char*) memchr (pchStart, '\n', pEnd-pchStart);//after endlabel
        if(pchStart<pEnd && pchStart !=NULL)
        {
            //will fail is string label has spaces, tabs or newlines
            memtok(&pchStart,&pchEnd,toks);//get label
            std::string s(pchStart,pchEnd-pchStart);
            //trim(s);
            label = MDL::str2Label(s);
            if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
            {
                label = MDL_UNDEFINED; //ignore if not present in desiredLabels
            }
            //columnValues[column]->fromChar(szAux);
            //std::cerr << "rlabel " << MDL::label2Str(label) <<std::endl;
            pchStart=pchEnd;
            memtok(&pchStart,&pchEnd,toks);//get value
            if (label != MDL_UNDEFINED)
            {
                strncpy(szAux,pchStart,pchEnd-pchStart);
                szAux[pchEnd-pchStart]=NULL;
                MDObject * _mdObject = new MDObject(label);
                _mdObject->fromChar(szAux);
                addLabel(label);
                columnValues.push_back(_mdObject);//add the value here with a char

                setValue(activeObjId, *(_mdObject));

            }
            pchStart=pchEnd;
        }
        else
            break;
    }
    free(szAux);
    return pchEnd;
}
#endif

/** This function will be used to parse the rows data
 * having read the columns labels before and setting wich are desired
 * the useCommentAsImage is for compatibility with old DocFile format
 * where the image were in comments
 */
void MetaData::_readRows(std::istream& is, MDRow & columnValues, bool useCommentAsImage)
{
    std::string line = "";
    while (!is.eof() && !is.fail())
    {
        //Move until the ';' or the first alphanumeric character
        while (is.peek() != ';' && isspace(is.peek()) && !is.eof())
            is.ignore(1);
        if (!is.eof())

            if (is.peek() == ';')//is a comment
            {
                is.ignore(1); //ignore the ';'
                getline(is, line);
                trim(line);
            }
            else if (!isspace(is.peek()))
            {
                addObject();
                if (line != "")//this is for old format files
                {
                    if (!useCommentAsImage)
                        setValue(MDL_COMMENT, line);
                    else
                        setValue(MDL_IMAGE, line);
                }
                int nCol = columnValues.size();
                for (int i = 0; i < nCol; ++i)
                    _parseObject(is, *(columnValues[i]));
            }

    }
}

/** This function will be used to parse the rows data in START format
 * @param[out] columnValues MDRow with values to fill in
 * @param pchStar pointer to the position of '_loop' in memory
 * @param pEnd  pointer to the position of the next '_data' in memory
 */
void MetaData::_readRowsStar(MDRow & columnValues, char * pchStart, char * pEnd)
{
    char * pchEnd;
    String line;
    std::stringstream ss;
    int nCol = columnValues.size();

    while (pchStart < pEnd)//while there are data lines
    {
        pchEnd = END_OF_LINE();
        line.assign(pchStart, pchEnd-pchStart);
        trim(line);
        if (line != "")
        {
            std::stringstream ss(line);
            addObject();
            for (int i = 0; i < nCol; ++i)
                _parseObject(ss, *(columnValues[i]));
        }
        pchStart = pchEnd + 1;//go to next line
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
        MDObject value(label);
        os >> value;
        if (label != MDL_UNDEFINED)
            setValue(value, objectID);
    }
}
void MetaData::read(const FileName &_filename,
                    const std::vector<MDLabel> *desiredLabels,
                    bool decomposeStack)
{
    std::string BlockName;
    FileName filename;
    BlockName = _filename.getBlockName();
    filename  = _filename.removeBlockName();
    _read(filename,desiredLabels,BlockName,decomposeStack);

}

void MetaData::_read(const FileName &filename,
                     const std::vector<MDLabel> *desiredLabels,
                     const std::string & blockName, bool decomposeStack)
{
    //First try to open the file as a metadata
    _clear();
    myMDSql->createMd();
    isColumnFormat = true;

    std::string dataBlockName = (std::string) "data_" + blockName;

    if (!filename.isMetaData())//if not a metadata, try to read as image or stack
    {
        Image<char> image;
        image.read(filename, false);
        if (image().ndim == 1 || !decomposeStack) //single image
        {
            addObject();
            setValue(MDL_IMAGE, filename);
            setValue(MDL_ENABLED, 1);
        }
        else //stack
        {
            FileName fnTemp;
            for (size_t i = 0; i < image().ndim; ++i)
            {
                fnTemp.compose(i, filename);
                addObject();
                setValue(MDL_IMAGE, fnTemp);
                setValue(MDL_ENABLED, 1);
            }
        }
        return;
    }

    std::ifstream is(filename.data(), std::ios_base::in);
    std::stringstream ss;
    std::string line, token;
    MDRow columnValues;

    getline(is, line); //get first line to identify the type of file

    if (is.fail())
    {
        REPORT_ERROR(ERR_IO_NOTEXIST, (std::string) "MetaData::read: File " + filename + " does not exists" );
    }

    bool useCommentAsImage = false;
    this->inFile = filename;
    bool oldFormat=true;

    is.seekg(0, std::ios::beg);//reset the stream position to the beginning to start parsing

    if (line.find("XMIPP_STAR_1") != std::string::npos)
    {
        //map file
        struct stat file_status;
        int fd;
        char *map;
        oldFormat=false;

        if(stat(filename.data(), &file_status) != 0)
            REPORT_ERROR(ERR_IO_NOPATH,"Metadata:isColumnFormat can not get filesize for file "+filename);
        size_t size = file_status.st_size;
        if(size==0)
            REPORT_ERROR(ERR_IO_NOPATH,"Metadata:isColumnFormat: File size=0, can not read it ("+filename+")");

        fd = open(filename.data(),  O_RDWR, S_IREAD | S_IWRITE);
        if (fd == -1)
            REPORT_ERROR(ERR_IO_NOPATH,"Metadata:isColumnFormat can not read file named "+filename);

        map = (char *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (map == MAP_FAILED)
            REPORT_ERROR(ERR_MEM_BADREQUEST,"Metadata:write can not map memory ");

        char * firstData, * secondData, * firstloop;
        isColumnFormat=isColumnFormatFile(map,
                                          &firstData,
                                          &secondData,
                                          &firstloop,
                                          dataBlockName.data());
        if(secondData ==NULL)//this should not be necessary but you never know
            secondData = map + size;

        is.ignore(256,'#');
        is.ignore(256,'#');
        getline(is, line);
        trim(line);
        setComment(line);//this stream is no longer needed
        //        if (isColumnFormat)
        //        {
        //Read column labels from the datablock that starts at firstData
        //Label ends at firstloop
        char * aux;
        if (isColumnFormat)
        {
            aux = _readColumnsStar(firstloop, secondData, columnValues, desiredLabels);
            _readRowsStar(columnValues, aux, secondData);
        }
        else
        {
            addObject();
            _readColumnsStar(firstData,secondData, columnValues, desiredLabels);
        }
        //        }
        //        else//row format
        //        {
        //            _readRowFormatStar(firstData,secondData, columnValues, desiredLabels);
        //        }
        if (munmap(map, size) == -1)
        {
            REPORT_ERROR(ERR_MEM_NOTDEALLOC,"metadata:write, Can not unmap memory");
        }
        close(fd);
    }
    else if (line.find("Headerinfo columns:") != std::string::npos)
    {
        //This looks like an old DocFile, parse header
        std::cerr << "WARNING: ** You are using an old file format (DOCFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        is.ignore(256, ':'); //ignore all until ':' to start parsing column labels
        getline(is, line);
        ss.str(line);
        columnValues.push_back(new MDObject(MDL_UNDEFINED));
        columnValues.push_back(new MDObject(MDL_UNDEFINED));

        addLabel(MDL_IMAGE);
        _readColumns(ss, columnValues, desiredLabels);
        useCommentAsImage = true;
    }
    else
    {
        std::cerr << "WARNING: ** You are using an old file format (SELFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        //I will assume that is an old SelFile, so only need to add two columns
        columnValues.push_back(new MDObject(MDL_IMAGE));//addLabel(MDL_IMAGE);
        columnValues.push_back(new MDObject(MDL_ENABLED));//addLabel(MDL_ENABLED);
    }

    if (oldFormat)
        _readRows(is, columnValues, useCommentAsImage);
    //    else
    //        _readRowFormat(is);
    firstObject();
}

void MetaData::merge(const FileName &fn)
{
    MetaData md;
    md.read(fn);
    unionAll(md);
}

void MetaData::aggregateSingle(MDObject &mdValueOut, AggregateOperation op,
                               MDLabel aggregateLabel)

{
    mdValueOut.setValue(myMDSql->aggregateSingleDouble(op,aggregateLabel));
}

bool MetaData::isColumnFormatFile(char * map,
                                  char ** firstData,
                                  char ** secondData,
                                  char ** firstloop,
                                  const char * szBlockName)
{
    //search for the first data_XXX line
    //then for the next data_ line if any
    //rewind to the first hit and look for loop_
    //std::string _szBlockName = (std::string)("data_") + blockName;
    *firstData  = (char *)  strstr(map, szBlockName);
    if(*firstData==NULL)
    	REPORT_ERROR(ERR_MD_WRONGDATABLOCK,(std::string) "Block Named: " + szBlockName + " does not exist");
    *secondData = (char *)  strstr((*firstData+1),"data_");
    *firstloop  = (char *)  strstr((*firstData),"loop_");
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "in_map:  "         << (void *) map  << std::endl;
    std::cerr << "in_firstData:  "   << (void *) *firstData  << std::endl;
    std::cerr << "in_secondData: "   << (void *) *secondData << std::endl;
    std::cerr << "in_firstLoop:  "   << (void *) *firstloop  << std::endl;
#endif
#undef DEBUG

    if (*firstloop==NULL)
        return (false);
    else if(*secondData == NULL || (*secondData) > (*firstloop))
        return (true);
    return (false);
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
    firstObject();
}

void MetaData::aggregate(const MetaData &mdIn, const std::vector<AggregateOperation> &ops,
                         MDLabel operateLabel, const std::vector<MDLabel> &resultLabels)
{
    if (resultLabels.size() - ops.size() != 1)
        REPORT_ERROR(ERR_MD, "Labels vectors should contain one element more than operations");
    init(&resultLabels);
    mdIn.myMDSql->aggregateMd(this, ops, operateLabel);
    firstObject();
}


//-------------Set Operations ----------------------
void MetaData::_setOperates(const MetaData &mdIn, const MDLabel label, SetOperation operation)
{
    if (this == &mdIn) //not sense to operate on same metadata
        REPORT_ERROR(ERR_MD, "Couldn't perform this operation on input metadata");
    if (size() == 0 && mdIn.size() == 0)
        REPORT_ERROR(ERR_MD, "Couldn't perform this operation if both metadata are empty");
    //Add labels to be sure are present
    for (int i = 0; i < mdIn.activeLabels.size(); i++)
        addLabel(mdIn.activeLabels[i]);

    mdIn.myMDSql->setOperate(this, label, operation);
    firstObject();
}

void MetaData::_setOperates(const MetaData &mdInLeft, const MetaData &mdInRight, const MDLabel label, SetOperation operation)
{
    if (this == &mdInLeft || this == &mdInRight) //not sense to operate on same metadata
        REPORT_ERROR(ERR_MD, "Couldn't perform this operation on input metadata");
    //Add labels to be sure are present
    for (int i = 0; i < mdInLeft.activeLabels.size(); i++)
        addLabel(mdInLeft.activeLabels[i]);
    for (int i = 0; i < mdInRight.activeLabels.size(); i++)
        addLabel(mdInRight.activeLabels[i]);

    myMDSql->setOperate(&mdInLeft, &mdInRight, label, operation);
    firstObject();
}

void MetaData::unionDistinct(const MetaData &mdIn, const MDLabel label)
{
    if(mdIn.isEmpty())
        return;
    _setOperates(mdIn, label, UNION_DISTINCT);
}

void MetaData::unionAll(const MetaData &mdIn)
{
    if(mdIn.isEmpty())
        return;
    _setOperates(mdIn, MDL_UNDEFINED, UNION);//label not needed for unionAll operation
}


void MetaData::intersection(const MetaData &mdIn, const MDLabel label)
{
    if(mdIn.isEmpty())
        clear();
    else
        _setOperates(mdIn, label, INTERSECTION);
}
void MetaData::subtraction(const MetaData &mdIn, const MDLabel label)
{
    if(mdIn.isEmpty())
        return;
    _setOperates(mdIn, label, SUBSTRACTION);
}

void MetaData::join(const MetaData &mdInLeft, const MetaData &mdInRight, const MDLabel label, JoinType type)
{
    _setOperates(mdInLeft, mdInRight, label, (SetOperation)type);
}

void MetaData::operate(const std::string &expression)
{
    if (!myMDSql->operate(expression))
        REPORT_ERROR(ERR_MD, "MetaData::operate: error doing operation");
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
    if (MDin.containsLabel(sortLabel))
    {
        init(&(MDin.activeLabels));
        copyInfo(MDin);
        MDin.myMDSql->copyObjects(this, new MDQuery(-1, 0, sortLabel));
    }
    else
        *this=MDin;
    firstObject();
}

void MetaData::sort(MetaData &MDin, const std::string &sortLabel)
{
    // Check if the label has semicolon
    int ipos=sortLabel.find(':');
    if (ipos!=std::string::npos || MDL::labelType(sortLabel)==LABEL_VECTOR)
    {
        MDLabel label;
        int column;
        if (ipos!=std::string::npos)
        {
            // Check that the label is a vector field
            std::vector< std::string > results;
            splitString(sortLabel,":",results);
            column=textToInteger(results[1]);
            if (MDL::labelType(results[0])!=LABEL_VECTOR)
                REPORT_ERROR(ERR_ARG_INCORRECT,"Column specifications cannot be used with non-vector labels");
            label=MDL::str2Label(results[0]);
        }
        else
        {
            label=MDL::str2Label(sortLabel);
            column=0;
        }

        // Get the column values
        MultidimArray<double> v;
        v.resizeNoCopy(MDin.size());
        std::vector<double> vectorValues;
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(MDin)
        {
            MDin.getValue(label,vectorValues);
            if (column>=vectorValues.size())
                REPORT_ERROR(ERR_MULTIDIM_SIZE,"Trying to access to inexistent column in vector");
            DIRECT_A1D_ELEM(v,i)=vectorValues[column];
            i++;
        }

        // Sort
        MultidimArray<int> idx;
        v.indexSort(idx);

        // Construct output Metadata
        init(&(MDin.activeLabels));
        copyInfo(MDin);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(idx)
        {
            MDRow row;
            MDin.getRow(row,DIRECT_A1D_ELEM(idx,i));
            addObject();
            setRow(row);
        }

        firstObject();
    }
    else
        sort(MDin, MDL::str2Label(sortLabel));
}

void MetaData::split(int n, std::vector<MetaData> &results, const MDLabel sortLabel)
{
    long int mdSize = size();
    if (n > mdSize)
        REPORT_ERROR(ERR_MD, "MetaData::split: Couldn't split a metadata in more parts than its size");

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
    mdIn.myMDSql->copyObjects(this, new MDQuery(n_images, first, sortLabel));
    firstObject();
}

void MetaData::selectSplitPart(const MetaData &mdIn, int n, int part, const MDLabel sortLabel)
{
    long int mdSize = mdIn.size();
    if (n > mdSize)
        REPORT_ERROR(ERR_MD, "selectSplitPart: Couldn't split a metadata in more parts than its size");
    if (part < 0 || part >= n)
        REPORT_ERROR(ERR_MD, "selectSplitPart: 'part' should be between 0 and n-1");
    _selectSplitPart(mdIn, n, part, mdSize, sortLabel);

}

void MetaData::selectPart(const MetaData &mdIn, long int startPosition, long int numberOfObjects,
                          const MDLabel sortLabel)
{
    long int mdSize = mdIn.size();
    if (startPosition < 0 || startPosition >= mdSize)
        REPORT_ERROR(ERR_MD, "selectPart: 'startPosition' should be between 0 and size()-1");
    init(&(mdIn.activeLabels));
    copyInfo(mdIn);
    mdIn.myMDSql->copyObjects(this, new MDQuery(numberOfObjects, startPosition, sortLabel));
    firstObject();
}


WriteModeMetaData metadataModeConvert (String mode)
{
    toLower(mode);
    if (mode.npos != mode.find("overwrite"))
        return OVERWRITE;
    if (mode.npos != mode.find("append"))
        return APPEND;
    REPORT_ERROR(ERR_ARG_INCORRECT,"metadataModeConvert: Invalid mode");
}

void MetaData::makeAbsPath(const MDLabel label)
{

    std::string aux_string;
    std::string aux_string_path;
    char buffer[1024];

    getcwd(buffer, 1023);
    std::string path_str(buffer);
    path_str += "/";

    this->firstObject();
    getValue(label,aux_string);

    if (aux_string[0] == '/')
        return;

    FOR_ALL_OBJECTS_IN_METADATA(*this)
    {
        aux_string_path = path_str;
        getValue(label,aux_string);
        aux_string_path += aux_string;
        setValue(label,aux_string_path);
    }
}
