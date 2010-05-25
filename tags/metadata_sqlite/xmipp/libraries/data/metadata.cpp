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
void MetaData::clear()
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
    MDSql::clearMd(this);
}//close clear

void MetaData::init(const std::vector<MDLabel> *labelsVector)
{
    tableId = MDSql::getMdUniqueId();
    clear();
    if (labelsVector != NULL)
        this->activeLabels = *labelsVector;
    //Create table in database
    MDSql::createMd(this);
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
    MDSql::copyObjects(&md, this);
}

bool MetaData::setValueObject(long int objId, const MDValue &mdValueIn)
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
    MDSql::setObjectValue(this, objId, mdValueIn);
}

bool MetaData::getValueObject(long int objId, MDValue &mdValueOut) const
{
    REQUIRE_LABEL_EXISTS(mdValueOut.label, "getValue");
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
    MDSql::getObjectValue(this, objId, mdValueOut);
    //MDL::value2VoidPtr(label, mdValue, valuePtrOut);
}

MetaData::MetaData()
{
    init(NULL);
}//close MetaData default Constructor

MetaData::MetaData(const std::vector<MDLabel> *labelsVector)
{
    init(labelsVector);
}//close MetaData default Constructor

MetaData::MetaData(const FileName &fileName, const std::vector<MDLabel> *labelsVector)
{
    init(labelsVector);
    //TODO: Read values from File and store in MD
}//close MetaData from file Constructor

MetaData::MetaData(const MetaData &md)
{
    copyMetadata(md);
}//close MetaData copy Constructor

MetaData& MetaData::operator =(const MetaData &md)
{
    copyMetadata(md);
    return *this;
}//close metadata operator =

MetaData::~MetaData()
{
    clear();
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
    path = (newPath == "") ? std::string(getcwd(NULL, 0)) : newPath;
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

int MetaData::MaxStringLength(const MDLabel thisLabel) const
{
   REQUIRE_LABEL_EXISTS(thisLabel, "MaxStringLength");

    return MDSql::columnMaxLength(this, thisLabel);
}

bool MetaData::setValueFromStr(const MDLabel label, const std::string &value, long int objectId)
{
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
    MDSql::setObjectValue(this, objectId, mdValue);
}

bool MetaData::isEmpty() const
{
    return size() == 0;
}

size_t MetaData::size() const
{
    std::vector<long int> objects = MDSql::selectObjects(this);

    return objects.size();
}

bool MetaData::containsLabel(const MDLabel label) const
{
    int size = activeLabels.size();
    for (int i = 0; i < size; i++)
        if (label == activeLabels[i])
            return true;
    return false;
}

bool MetaData::addLabel(const MDLabel label)
{
    if (containsLabel(label))
        return false;
    activeLabels.push_back(label);
    MDSql::addColumn(this, label);
    return true;
}

long int MetaData::addObject(long int objectId)
{
    activeObjId = MDSql::addRow(this);
    return activeObjId;
}

void MetaData::importObject(const MetaData &md, const long int objId)
{
    importObjects(md, MDValueEqual(MDL_OBJID, objId));
}

void MetaData::importObjects(const MetaData &md, const std::vector<long int> &objectsToAdd)
{
    int size = objectsToAdd.size();
    for (int i = 0; i < size; i++)
        importObject(md, objectsToAdd[i]);
}

void MetaData::importObjects(const MetaData &md, MDQuery query)
{
    init(&(md.activeLabels));
    copyInfo(md);
    MDSql::copyObjects(&md, this, &query);
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

int MetaData::removeObjects(MDQuery query)
{
    int removed = MDSql::deleteObjects(this, &query);
    firstObject(); //I prefer put active object to -1
    return removed;
}


//----------Iteration functions -------------------
long int MetaData::firstObject()
{
    //std::vector<long int> objects = MDSql::selectObjects(this, 1);
    //return (objects.size() > 0) ? objects[0] : -1;
    activeObjId = MDSql::firstRow(this);
    return activeObjId;
}

long int MetaData::lastObject()
{
    activeObjId = MDSql::lastRow(this);
    return activeObjId;
}

long int MetaData::nextObject()
{
    if (activeObjId == -1)
        REPORT_ERROR(-55, "Couldn't call 'nextObject' when 'activeObject' is -1");
    activeObjId = MDSql::nextRow(this, activeObjId);
    return activeObjId;
}

long int MetaData::previousObject()
{
    if (activeObjId == -1)
        REPORT_ERROR(-55, "Couldn't call 'previousObject' when 'activeObject' is -1");
    activeObjId = MDSql::previousRow(this, activeObjId);
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
std::vector<long int> MetaData::findObjects(MDQuery query, int limit)
{
    //FIXME: VECTORS ARE RETURNED BY VALUE
    std::vector<long int> objects = MDSql::selectObjects(this, limit, &query);
    return objects;
}

int MetaData::countObjects(MDQuery query)
{
    return findObjects(query).size();
}

bool MetaData::containsObject(long int objectId)
{
    return containsObject(MDValueEqual(MDL_OBJID, objectId));
}

bool MetaData::containsObject(MDQuery query)
{
    return findObjects(query, 1).size() > 0;
}

long int MetaData::gotoFirstObject(MDQuery query)
{
    std::vector<long int> objects = findObjects(query, 1);

    activeObjId = objects.size() == 1 ? objects[0] : -1;
    return activeObjId;
}


//--------------IO functions -----------------------

void MetaData::write(const FileName &outFile)
{
    std::cerr << "writing to file: " << outFile <<std::endl;
    // Open file
    std::ofstream ofs(outFile.data(), std::ios_base::out);
    write(ofs);
}

void MetaData::write(std::ostream &os)
{
    os << "; XMIPP_3 * " << (isColumnFormat ? "column" : "row")
    << "_format *" << std::endl;
    //TODO: not saving the path
    os << "; " << comment << std::endl;

    if (isColumnFormat)
    {
        //Write metadata header
        //write columns
        int labelsSize = activeLabels.size();
        os << "; ";
        for (int i = 0; i < labelsSize; i++)
        {
            if (activeLabels.at(i) != MDL_COMMENT)
            {
                os.width(10);
                os << MDL::label2Str(activeLabels.at(i)) << " ";
            }
        }
        os << std::endl;
        //Write data
        std::vector<long int> objects = MDSql::selectObjects(this);
        int objsSize = objects.size();
        for (int o = 0; o < objsSize; o++)
        {
            for (int i = 0; i < labelsSize; i++)
            {
                if (activeLabels.at(i) != MDL_COMMENT)
                {
                    MDValue mdValue(activeLabels[i]);
                    os.width(10);
                    MDSql::getObjectValue(this, objects.at(o), mdValue);
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
                    MDSql::getObjectValue(this, objId, mdValue);
                    mdValue.toStream(os);
                    os << std::endl;
                }
            }
        }

    }
}//write

void MetaData::read(const FileName &inFile)
{
    std::ifstream ifs(inFile.data(), std::ios_base::in);
    read(ifs);
}

void MetaData::read(std::istream &is)
{
    std::vector<MDLabel> activeLabels;

    is.seekg(0, std::ios::beg);
    std::string line;
    getline(is, line, '\n');
    int pos = line.find("*");

    if (pos == std::string::npos)
    {
        REPORT_ERROR( 200, "End of string reached" );
    }
    else
    {
        line.erase(0, pos + 1);
        pos = line.find(" row_format ");

        if (pos != std::string::npos)
        {
            isColumnFormat = false;
        }
    }

    pos = line.find("*");
    line.erase(0, pos + 1);
    line = removeChar(line, ' ');
    setPath(line);
    getline(is, line, '\n');
    setComment(line.erase(0, 2));

    if (isColumnFormat)
    {
        // Get Labels line
        getline(is, line, '\n');
        // Remove ';'
        line.erase(0, line.find(";") + 1);
        // Parse labels
        std::stringstream os(line);
        std::string newLabel;
        int labelPosition = 0;

        while (os >> newLabel)
        {
            MDLabel label = MDL::str2Label(newLabel);

            if (label == MDL_UNDEFINED)
                ignoreLabels.push_back(labelPosition);
            else
                activeLabels.push_back(label);
            labelPosition++;
        }

        init(&activeLabels);
        // Read data and fill structures accordingly
        while (getline(is, line, '\n'))
        {
            if (line[0] == '\0' || line[0] == '#')
                continue;

            long int objectID = addObject();

            if (line[0] == ';')
            {
                line.erase(0, 1);
                line = simplify(line);
                setValue(MDL_COMMENT, line);
                getline(is, line, '\n');
            }

            // Parse labels
            std::stringstream os2(line);
            std::string value;

            int labelPosition = 0;
            int counterIgnored = 0;

            while (os2 >> value)
            {
                if (std::find(ignoreLabels.begin(), ignoreLabels.end(),
                              labelPosition) != ignoreLabels.end())
                {
                    // Ignore this column
                    counterIgnored++;
                    labelPosition++;
                    continue;
                }

                if (MDL::isVector(activeLabels[labelPosition - counterIgnored])
                    && value == "**")
                {
                    std::string aux;
                    while (os2 >> value)
                        if (value == "**")
                            break;
                        else
                            aux += value + " ";
                    value = aux;
                    std::cerr << "is vector value" << value << std::endl;
                }

                setValueFromStr(activeLabels[labelPosition - counterIgnored], value);
                labelPosition++;
            }
        }
    }
    else //RowFormat
    {
        std::string newLabel;
        std::string value;

        long int objectID = addObject();

        // Read data and fill structures accordingly
        while (getline(is, line, '\n'))
        {
            if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
                continue;

            // Parse labels
            std::stringstream os(line);

            os >> newLabel;
            MDLabel label = MDL::str2Label(newLabel);

            if(!MDL::isVector(label))
                os >> value;
            else
            {
                std::vector<std::string> v;
                Tokenize(line,v,(std::string)"**");
                value = v[1];
            }
            if (label != MDL_UNDEFINED)
            {
                activeLabels.push_back(label);
                setValueFromStr(label, value);
            }
        }
    }
}//close read

//-------------Set Operations ----------------------
void MetaData::union_(const MetaData &MD, const MDLabel thisLabel)
{}
void MetaData::unionAll(const MetaData &MD)
{}
void MetaData::aggregate(MetaData MDIn,
                         MDLabel aggregateLabel,
                         MDLabel entryLabel,
                         MDLabel operationLabel)
{}

void MetaData::merge(const FileName &fn)
{}
void MetaData::intersection(MetaData & minuend, MetaData & , MDLabel thisLabel)
{}
void MetaData::substraction(MetaData & minuend, MetaData & subtrahend,
                            MDLabel thisLabel)
{}


/*----------   Statistics --------------------------------------- */
#include "metadata_extension.h"
void get_statistics(MetaData MT_in, Image<double> & _ave, Image<double> & _sd, double& _min,
                    double& _max, bool apply_geo)
{
    REPORT_ERROR(-55, "get_statistics not yet implemented");
}

void ImgSize(MetaData &MD, int &Xdim, int &Ydim, int &Zdim, int &Ndim)
{
    if (!MD.isEmpty())
    {
        FileName fn_img;
        Image<double> img;
        MD.getValue(MDL_IMAGE, fn_img);
        img.read(fn_img, false);
        img.getDimensions(Xdim, Ydim, Zdim, Ndim);
    }
    else
        REPORT_ERROR(-1, "Can not read image size from empty metadata");
}
