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

MetaData::MetaData()
{
    clear();
}

MetaData::MetaData(MetaData &MD)
{
    clear();
    this->setComment(MD.getComment());
    this->setPath(MD.getPath());
    this->isColumnFormat = MD.isColumnFormat;
    this->inFile = MD.inFile;
    this->fastStringSearchLabel = MDL_UNDEFINED;
    this->activeLabels = MD.activeLabels;

    //objects, define iterator
    std::map<long int, MetaDataContainer *>::iterator objIt;
    for (objIt = MD.objects.begin(); objIt != MD.objects.end(); objIt++)
    {
        long int idx = this->addObject();

        this->objects[idx] = new MetaDataContainer(*(objIt->second));
    }
    this->objectsIterator = objects.begin();

}

MetaData& MetaData::operator =(MetaData &MD)
{
    clear();
    if (this != &MD)
    {
        this->setComment(MD.getComment());
        this->setPath(MD.getPath());
        this->isColumnFormat = MD.isColumnFormat;
        this->inFile = MD.inFile;
        this->fastStringSearchLabel = MDL_UNDEFINED;
        this->activeLabels = MD.activeLabels;

        //objects, define iterator
        std::map<long int, MetaDataContainer *>::iterator objIt;
        for (objIt = MD.objects.begin(); objIt != MD.objects.end(); objIt++)
        {
            long int idx = this->addObject();

            this->objects[idx] = new MetaDataContainer(*(objIt->second));
        }
        this->objectsIterator = objects.begin();
    }
    return *this;
}

MetaData& MetaData::add (MetaData &MD)
{
    std::map<long int, MetaDataContainer *>::iterator objIt;
    for (objIt = MD.objects.begin(); objIt != MD.objects.end(); objIt++)
    {
        long int idx = this->addObject();
        objects[idx] = new MetaDataContainer(*(objIt->second));
    }
    this->objectsIterator = objects.begin();
    return *this;
}

void MetaData::read(std::ifstream *infile,
                    std::vector<MetaDataLabel> * labelsVector)
{
    infile->seekg(0, std::ios::beg);
    std::string line;

    getline(*infile, line, '\n');

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
    getline(*infile, line, '\n');
    setComment(line.erase(0, 2));

    if (isColumnFormat)
    {
        // Get Labels line
        getline(*infile, line, '\n');

        // Remove ';'
        line.erase(0, line.find(";") + 1);

        // Parse labels
        std::stringstream os(line);
        std::string newLabel;

        int labelPosition = 0;

        while (os >> newLabel)
        {
            MetaDataLabel label = MetaDataContainer::codifyLabel(newLabel);

            if (label == MDL_UNDEFINED)
                ignoreLabels.push_back(labelPosition);
            else
                activeLabels.push_back(label);

            labelPosition++;
        }

        // Read data and fill structures accordingly
        while (getline(*infile, line, '\n'))
        {
            if (line[0] == '\0' || line[0] == '#')
                continue;

            long int objectID = addObject();

            if (line[0] == ';')
            {
                line.erase(0, 1);
                line = simplify(line);
                setValue(MDL_COMMENT, line);
                getline(*infile, line, '\n');
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

                if (isVector(activeLabels[labelPosition - counterIgnored])
                    && value == "**")
                {
                    std::string aux;
                    while (os2 >> value)
                        if (value == "**")
                            break;
                        else
                            aux += value + " ";
                    value = aux;
                }

                setValue(MetaDataContainer::decodeLabel(
                             activeLabels[labelPosition - counterIgnored]), value);
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
        while (getline(*infile, line, '\n'))
        {
            if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
                continue;

            // Parse labels
            std::stringstream os(line);

            os >> newLabel;
            os >> value;

            MetaDataLabel label = MetaDataContainer::codifyLabel(newLabel);
            if (label != MDL_UNDEFINED)
            {
                activeLabels.push_back(label);
                setValue(newLabel, value);
            }
        }
    }
}

void MetaData::readOldDocFile(std::ifstream *infile,
                              std::vector<MetaDataLabel> * labelsVector)
{
    bool saveName = true;
    infile->seekg(0, std::ios::beg);
    std::string line;

    // Search for Headerinfo, if present we are processing an old-styled docfile
    // else we are processing a new Xmipp MetaData file
    getline(*infile, line, '\n');

    int pos = line.find("Headerinfo");

    // Remove from the beginning to the end of "Headerinfo columns:"
    line = line.erase(0, line.find(":") + 1);

    // In the old docfile format the "image" label did not exist, it was
    // a ";" commented line containing the name of the projection. Therefore,
    // for the new format it must be added by hand, if necessary
    if (labelsVector != NULL)
    {
        std::vector<MetaDataLabel>::iterator location;

        location = std::find(labelsVector->begin(), labelsVector->end(),
                             MDL_IMAGE);

        if (location != labelsVector->end())
        {
            activeLabels.push_back(MDL_IMAGE);
            saveName = true;
        }
    }
    else
    {
        activeLabels.push_back(MDL_IMAGE);
    }

    // Extract labels until the string is empty
    while (line != "")
    {
        pos = line.find(")");
        std::string newLabel = line.substr(0, pos + 1);
        line.erase(0, pos + 1);

        // The token can now contain a ',', if so, remove it
        if ((pos = newLabel.find(",")) != std::string::npos)
            newLabel.erase(pos, 1);

        // Remove unneded parentheses and contents
        pos = newLabel.find("(");
        newLabel.erase(pos, newLabel.find(")") - pos + 1);

        // Remove white spaces
        newLabel = removeChar(newLabel, ' ');

        if (labelsVector != NULL)
        {
            std::vector<MetaDataLabel>::iterator location;

            location = std::find(labelsVector->begin(), labelsVector->end(),
                                 MetaDataContainer::codifyLabel(newLabel));

            if (location != labelsVector->end())
            {
                activeLabels.push_back(MetaDataContainer::codifyLabel(newLabel));
            }
        }
        else
        {
            activeLabels.push_back(MetaDataContainer::codifyLabel(newLabel));
        }
    }

    int isname = 0;
    while (getline(*infile, line, '\n'))
    {
        if (isname % 2 == 0)
        {
            long int objectID = addObject();
            line.erase(0, line.find(";") + 1);

            // Remove spaces from string
            line = removeChar(line, ' ');
            setValue(MDL_IMAGE, line);
        }
        else
        {
            // Parse labels
            std::stringstream os2(line);
            std::string value;

            int counter = 0;
            while (os2 >> value)
            {
                if (counter >= 2) // Discard two first numbers
                {
                    if (saveName)
                        setValue(MetaDataContainer::decodeLabel(
                                     activeLabels[counter - 1]), value);
                    else
                        setValue(MetaDataContainer::decodeLabel(
                                     activeLabels[counter - 2]), value);
                }
                counter++;
            }
        }

        isname++;
    }
}

void MetaData::readOldSelFile(std::ifstream *infile)
{
    infile->seekg(0, std::ios::beg);
    std::string line;

    activeLabels.push_back(MDL_IMAGE);
    activeLabels.push_back(MDL_ENABLED);

    while (getline(*infile, line, '\n'))
    {
        line = simplify(line);
        if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
            continue;
        else
        {
            int pos = line.find(" ");
            std::string name = line.substr(0, pos);
            line.erase(0, pos + 1);
            int i = atoi(line.c_str());
            addObject();
            setValue(MDL_IMAGE, name);
            setValue(MDL_ENABLED, i);
        }
    }
}

MetaData::MetaData(FileName fileName, std::vector<MetaDataLabel> * labelsVector)
{
    read(fileName, labelsVector);
}

void MetaData::setColumnFormat(bool column)
{
    isColumnFormat = column;
}

void MetaData::toDataBase(const FileName & DBname,
                          const std::string & tableName,
                          std::vector<MetaDataLabel> * labelsVector)
{
    try
    {
        CppSQLite3DB db;
        db.open(DBname,SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE);

        int labelsCounter = 0;
        std::string sqlCommandCreate = "create table ";

        if (tableName == "")
        {
            if (inFile == "")
                REPORT_ERROR( -1, "Please, supply a name for the table" );
            else
                sqlCommandCreate += inFile;
        }
        else
            sqlCommandCreate += tableName;
        sqlCommandCreate += "(objId int primary key,";
        std::vector<MetaDataLabel>::iterator strIt;
        std::string sqlCommandInsert = "insert into " + \
                                       tableName + \
                                       "(objId," ;
        for (strIt = activeLabels.begin();
             strIt != activeLabels.end();
             strIt++)
        {
            if (labelsVector==NULL || std::find(labelsVector->begin(),
                                                labelsVector->end(),
                                                *strIt) != labelsVector->end())
            {
                std::string label = MetaDataContainer::decodeLabel(*strIt);
                sqlCommandInsert += label + ",";
                if (isDouble(*strIt))
                {
                    sqlCommandCreate += label + " double,";
                    labelsCounter++;
                }
                else if (isString(*strIt))
                {
                    sqlCommandCreate += label + " text,";
                    labelsCounter++;
                }
                else if (isInt(*strIt))
                {
                    sqlCommandCreate += label + " int,";
                    labelsCounter++;
                }
                else if (isBool(*strIt))
                {
                    sqlCommandCreate += label + " int,";
                    labelsCounter++;
                }
                else if (isVector(*strIt))
                {
                    std::cerr
                    << "SQLLITE does not support vectors, skipping vector labeled "
                    << label << std::endl;
                }
            }
        }//end table creation
        // remove last ','
        sqlCommandCreate.erase(sqlCommandCreate.end() - 1);
        sqlCommandInsert.erase(sqlCommandInsert.end() - 1);
        sqlCommandCreate += ");";
        sqlCommandInsert += ") values(";
        for (int i = 0; i <= labelsCounter; i++)
        {
            sqlCommandInsert += "?,";
        }
        sqlCommandInsert.erase(sqlCommandInsert.end() - 1);
        sqlCommandInsert += ");";
        if (!db.tableExists(tableName))
            db.execDML(sqlCommandCreate);

        db.execDML("begin transaction;");
        //std::vector<MetaDataLabel>::iterator strIt;
        CppSQLite3Statement stmt = db.compileStatement(sqlCommandInsert);

        for (long int IDthis = firstObject(); IDthis != NO_MORE_OBJECTS; IDthis
             = nextObject())
        {
            stmt.bind(1,(int)IDthis);
            labelsCounter = 1;
            /*
             *
             for (strIt = activeLabels.begin();
                 strIt != activeLabels.end();
                 strIt++)
            */


            MetaDataContainer * aux = getObject();
            //save objID - primary key-
            //notice that you can not add two object with the same objID
            stmt.bind(labelsCounter,(int)IDthis);
            for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
            {
                if (labelsVector==NULL || std::find(labelsVector->begin(),
                                                    labelsVector->end(),
                                                    *strIt) != labelsVector->end())
                {
                    if (isDouble(*strIt))
                    {
                        labelsCounter++;
                        double value;
                        getValue(*strIt, value);
                        stmt.bind(labelsCounter, value);
                    }
                    else if (isString(*strIt))
                    {
                        labelsCounter++;
                        std::string value;
                        getValue(*strIt, value);
                        stmt.bind(labelsCounter, value.c_str());
                    }
                    else if (isInt(*strIt))
                    {
                        labelsCounter++;
                        int value;
                        getValue(*strIt, value);
                        stmt.bind(labelsCounter, value);
                    }
                    else if (isBool(*strIt))
                    {
                        labelsCounter++;
                        bool value;
                        getValue(*strIt, value);
                        stmt.bind(labelsCounter, value);
                    }
                }
            }
            stmt.execDML();
            stmt.reset();
        }
        db.execDML("commit transaction;");
        db.close();
    }
    catch (CppSQLite3Exception& e)
    {
        REPORT_ERROR( e.errorCode(), e.errorMessage());
    }
}
void MetaData::fromDataBase(const FileName & DBname,
                            const std::string & tableName,
                            std::vector<MetaDataLabel> * labelsVector)
{
    //check that database and table exits
    try
    {
        std::vector<MetaDataLabel> _labelsVector;
        std::string sqlCommand;
        std::vector<MetaDataLabel>::iterator strIt;
        if (!exists(DBname))
            REPORT_ERROR( 1, (std::string)"database "+ DBname + " does not exists");
        CppSQLite3DB db;
        db.open(DBname,SQLITE_OPEN_READONLY);
        if (!db.tableExists(tableName) )
            REPORT_ERROR( 1, (std::string)"table "+ tableName + " does not exists");
        clear();
        // we can use db.getTable(sqlCommand) but it is slow so let us go low level.
        //IF NULL get table column names
        /*
        PRAGMA table_info([table]) returns:

        [0] = column number
        [1] = column name
        [2] = column type
        [3] = flag: is the column NOT NULL?
        [4] = default value of the column
        [5] = flag: is this column part of the table's PRIMARY KEY?
        */
        if (labelsVector==NULL)
        {
            sqlCommand=(std::string)"PRAGMA table_info(" + tableName + ");";
            CppSQLite3Query q = db.execQuery(sqlCommand);
            while (!q.eof())
            {
                _labelsVector.push_back(MetaDataContainer::codifyLabel(q.fieldValue(1)));
                q.nextRow();
            }
        }
        else
        {
            _labelsVector.push_back(MDL_OBJID);
            for (strIt  = (*labelsVector).begin();
                 strIt != (*labelsVector).end();
                 strIt++)
            {
                _labelsVector.push_back(*strIt);
            }
        }
        sqlCommand = "select ";
        for (strIt  = (_labelsVector).begin();
             strIt != (_labelsVector).end();
             strIt++)
        {
            sqlCommand += MetaDataContainer::decodeLabel(*strIt) + ",";
        }
        sqlCommand.erase(sqlCommand.end() - 1);
        sqlCommand += (std::string)" from " + tableName;
        sqlCommand += " order by 1;"; //sort always by objId


        CppSQLite3Query q = db.execQuery(sqlCommand);
        int labelsCounter = 0;

        while (!q.eof())
        {
            labelsCounter = 0;
            int _objID=q.getIntField(labelsCounter);
            addObject(_objID);
            for (strIt  = (_labelsVector).begin()+1;
                 strIt != (_labelsVector).end();
                 strIt++)
            {
                if (isDouble(*strIt))
                {
                    labelsCounter++;
                    double value;
                    value=q.getFloatField(labelsCounter);
                    setValue(*strIt,value);
                }
                else if (isInt(*strIt))
                {
                    labelsCounter++;
                    int value;
                    value=q.getIntField(labelsCounter);
                    setValue(*strIt,value);
                }
                else if (isString(*strIt))
                {
                    labelsCounter++;
                    std::string value;
                    value.assign( q.getStringField(labelsCounter) );
                    setValue(*strIt,value);
                }
                else if (isBool(*strIt))
                {
                    labelsCounter++;
                    int value;
                    value=q.getIntField(labelsCounter);
                    setValue(*strIt,(bool)value);
                }
                else
                {
                    REPORT_ERROR( 1, (std::string)"unknown type in fromDataBase routine");
                }
            }
            q.nextRow();
        }

        //        CppSQLite3Table t = db.getTable(sqlCommand);
        db.close();
    }
    catch (CppSQLite3Exception& e)
    {
        REPORT_ERROR( e.errorCode(), e.errorMessage());
    }

}
void MetaData::combineWithFiles(MetaDataLabel thisLabel)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be combined" );
    }

    MetaData auxMetaData;
    MetaDataContainer * auxMetaDataContainer;

    for (long int IDthis = firstObject(); IDthis != NO_MORE_OBJECTS; IDthis
         = nextObject())
    {
        MetaDataContainer * aux = getObject();

        // Read fileName
        FileName fileName;
        aux->writeValueToString(fileName, thisLabel);

        // Read file
        auxMetaData.read(fileName);
        auxMetaDataContainer = auxMetaData.getObject();

        for (MetaDataLabel mdl = MDL_FIRST_LABEL; mdl <= MDL_LAST_LABEL; mdl
             = MetaDataLabel(mdl + 1))
        {
            if (auxMetaDataContainer->valueExists(mdl))
            {
                std::string value;

                auxMetaDataContainer->writeValueToString(value, mdl);

                setValue(MetaDataContainer::decodeLabel(mdl), value);
            }
        }
    }
}
void MetaData::combine(MetaData & other, MetaDataLabel thisLabel)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be combined" );
    }

    MetaDataContainer * aux, *aux2;
    std::string value1, value2;

    for (long int IDthis = firstObject(); IDthis != NO_MORE_OBJECTS; IDthis
         = nextObject())
    {
        aux = getObject();
        aux->writeValueToString(value1, thisLabel);
        for (long int IDother = other.firstObject(); IDother != NO_MORE_OBJECTS; IDother
             = other.nextObject())
        {
            aux2 = other.getObject();
            aux2->writeValueToString(value2, thisLabel);

            if (value2 == value1)
            {
                for (MetaDataLabel mdl =
                         static_cast<MetaDataLabel> (MDL_FIRST_LABEL); mdl
                     <= MDL_LAST_LABEL; mdl = MetaDataLabel(mdl + 1))
                {
                    if (aux2->
                        valueExists(mdl))
                    {
                        std::string value;

                        aux2->writeValueToString(value, mdl);

                        setValue(MetaDataContainer::decodeLabel(mdl), value);
                    }
                }
                break;
            }
        }
    }
}

void MetaData::read(FileName fileName,
                    std::vector<MetaDataLabel> * labelsVector)
{
    clear();
    inFile = fileName;

    // Open file
    std::ifstream infile(fileName.data(), std::ios_base::in);
    std::string line;

    if (infile.fail())
    {
        REPORT_ERROR( 200, (std::string) "File " + fileName + " does not exits" );
    }

    // Search for Headerinfo, if present we are processing an old-styled docfile
    // else we are processing a new Xmipp MetaData file
    getline(infile, line, '\n');

    int pos = line.find("Headerinfo");

    if (pos != std::string::npos) // Headerinfo token found
    {
        readOldDocFile(&infile, labelsVector);
        std::cerr
        << (std::string) "WARNING: ** You are using an old file format (DOCFILE) which is going "
        + "to be deprecated in next Xmipp release **"
        << std::endl;
    } else
    {
        pos = line.find("XMIPP_3 * ");

        if (pos != std::string::npos) // xmipp_3 token found
        {
            read(&infile, labelsVector);
        } else // We are reading an old selfile
        {
            readOldSelFile(&infile);
            std::cerr
            << (std::string) "WARNING: ** You are using an old file format (SELFILE) which is going "
            + "to be deprecated in next Xmipp release **"
            << std::endl;
        }
    }

    objectsIterator = objects.begin();

    infile.close();
}

void MetaData::writeValueToString(std::string & result,
                                  const std::string &inputLabel)
{
    MetaDataContainer * aux = getObject();
    aux->writeValueToString(result, MetaDataContainer::codifyLabel(inputLabel));
}

void MetaData::write(const std::string &fileName)
{
    // Open file
    std::ofstream outfile(fileName.data(), std::ios_base::out);

    outfile << "; ";
    if (isColumnFormat)
        outfile << "XMIPP_3 * column_format * ";
    else
        outfile << "XMIPP_3 * row_format * ";

    outfile << path << std::endl;
    outfile << "; ";
    outfile << comment;
    outfile << std::endl;

    std::map<long int, MetaDataContainer *>::iterator It;
    std::vector<MetaDataLabel>::iterator strIt;

    if (isColumnFormat)
    {
        outfile << "; ";
        for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            if (*strIt != MDL_COMMENT)
            {
                outfile << MetaDataContainer::decodeLabel(*strIt);
                outfile << " ";
            }
        }
        outfile << std::endl;

        for (It = objects.begin(); It != objects.end(); It++)
        {
            if ((It->second)->valueExists(MDL_COMMENT))
            {
                std::string entryComment;
                (It->second)->getValue(MDL_COMMENT, entryComment);
                if (entryComment != std::string(""))
                    outfile << "; " << entryComment << std::endl;
            }

            for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
            {
                if (*strIt != MDL_COMMENT)
                {
                    outfile.width(10);
                    (It->second)->writeValueToFile(outfile, *strIt);
                    outfile << " ";
                }
            }

            outfile << std::endl;
        }
    }
    else
    {
        // Get first object. In this case (row format) there is a single object
        MetaDataContainer * object = getObject();

        for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            if (*strIt != MDL_COMMENT)
            {
                outfile.width(20);
                outfile << MetaDataContainer::decodeLabel(*strIt);
                outfile << " ";
                object->writeValueToFile(outfile, *strIt);
                outfile << std::endl;
            }
        }
    }
}

MetaData::~MetaData()
{
    clear();
}

bool MetaData::isEmpty()
{
    return objects.empty();
}

void MetaData::clear()
{
    srand ( time(NULL) );
    path.clear();
    comment.clear();
    objects.clear();

    objectsIterator = objects.end();

    fastStringSearch.clear();
    ;
    fastStringSearchLabel = MDL_UNDEFINED;

    activeLabels.clear();
    ignoreLabels.clear();

    isColumnFormat = true;
    inFile = FileName::FileName();
}

void MetaData::setPath(std::string newPath)
{
    if (newPath == "")
    {
        path = std::string(getcwd(NULL, 0));
    }
    else
    {
        path = newPath;
    }
}

std::string MetaData::getPath()
{
    return path;
}

void MetaData::setComment(std::string newComment)
{
    if (newComment == "")
    {
        comment = (std::string) "No comment";
    }
    else
    {
        comment = newComment;
    }
}

std::string MetaData::getComment()
{
    return comment;
}

long int MetaData::addObject(long int objectID)
{
    typedef std::pair<long int, MetaDataContainer *> newPair;
    long int result;

    if (objectID == -1)
    {
        result = lastObject() + 1;

        objects.insert(newPair(result, new MetaDataContainer()));

        // Set iterator pointing to the newly added object
        objectsIterator = objects.end();
        objectsIterator--;
    }
    else
    {
        result = objectID;

        // Does the object already exist ?, remove
        objects.erase(objectID);

        // The object does not exist, create it
        objects.insert(newPair(objectID, new MetaDataContainer()));
        objectsIterator = objects.find(objectID);
    }

    // Set default values for the existing labels
    std::vector<MetaDataLabel>::iterator It;
    for (It = activeLabels.begin(); It != activeLabels.end(); It++)
    {
        (objectsIterator->second)->addValue(
            MetaDataContainer::decodeLabel(*It), std::string(""));
    }

#ifdef NEVERDEFINE
    //¿Que leches es esto?
    //randomizar siempre?
    // randomize no ha  sido inicializado
    std::cerr << "before randomize, print randomize size " << randomOrderedObjects.size() << std::endl;
    size_t randomIndex = rand() % randomOrderedObjects.size();
    randomOrderedObjects[randomIndex]=result;
    std::cerr << "after randomize" << std::endl;
#endif

    return result;
}

long int MetaData::firstObject()
{
    long int result = 0;

    if (!objects.empty())
    {
        objectsIterator = objects.begin();
        result = objectsIterator->first;
    }
    else
    {
        result = NO_OBJECTS_STORED; // Map is empty
    }

    return result;
}
;

long int MetaData::nextObject()
{
    long int result = 0;

    if (!objects.empty())
    {
        objectsIterator++;

        if (objectsIterator != objects.end())
        {
            result = objectsIterator->first;
        }
        else
        {
            result = NO_MORE_OBJECTS;
        }
    }
    else
    {
        result = NO_OBJECTS_STORED;
    }

    return result;
}
;

long int MetaData::lastObject()
{
    long int result = 0;

    if (!objects.empty())
    {
        objectsIterator = objects.end();
        objectsIterator--;
        result = objectsIterator->first;
    }
    else
    {
        result = NO_OBJECTS_STORED;
    }

    return result;
}
;

long int MetaData::goToObject(long int objectID)
{
    long int result = objectID;

    if (!objects.empty())
    {
        std::map<long int, MetaDataContainer *>::iterator It;

        It = objects.find(objectID);

        if (It == objects.end())
        {
            result = NO_OBJECT_FOUND;
        }
        else
        {
            objectsIterator = It;
        }
    }
    else
    {
        result = NO_OBJECTS_STORED;
    }

    return result;
}

bool MetaData::setValue(MetaDataLabel name, double value, long int objectID)
{
    long int auxID;

    if (!objects.empty() && MetaDataContainer::isValidLabel(name))
    {
        if (objectID == -1)
        {
            auxID = objectsIterator->first;
        }
        else
        {
            auxID = objectID;
        }

        MetaDataContainer * aux = objects[auxID];

        // Check whether label is correct (belongs to the enum in the metadata_container header
        // and whether it is present in the activeLabels vector. If not, add it to all the other
        // objects with default values
        std::vector<MetaDataLabel>::iterator location;
        std::map<long int, MetaDataContainer *>::iterator It;

        location = std::find(activeLabels.begin(), activeLabels.end(), name);

        if (location == activeLabels.end())
        {
            activeLabels.push_back(name);

            // Add this label to the rest of the objects in this class
            for (It = objects.begin(); It != objects.end(); It++)
            {
                if (It->second != aux)
                {
                    (It->second)->addValue(name, double());
                }
            }
        }

        aux->addValue(name, value);

        return true;
    }
    else
    {
        return false;
    }
}

bool MetaData::setValue(MetaDataLabel name, int value, long int objectID)
{
    long int auxID;

    if (!objects.empty() && MetaDataContainer::isValidLabel(name))
    {
        if (objectID == -1)
        {
            auxID = objectsIterator->first;
        }
        else
        {
            auxID = objectID;
        }

        MetaDataContainer * aux = objects[auxID];

        // Check whether label is correct (belongs to the enum in the metadata_container header
        // and whether it is present in the activeLabels vector. If not, add it to all the other
        // objects with default values
        std::vector<MetaDataLabel>::iterator location;
        std::map<long int, MetaDataContainer *>::iterator It;

        location = std::find(activeLabels.begin(), activeLabels.end(), name);

        if (location == activeLabels.end())
        {
            activeLabels.push_back(name);

            // Add this label to the rest of the objects in this class
            for (It = objects.begin(); It != objects.end(); It++)
            {
                if (It->second != aux)
                {
                    (It->second)->addValue(name, int());
                }
            }

        }

        aux->addValue(name, value);

        return true;
    }
    else
    {
        return false;
    }
}

bool MetaData::setValue(MetaDataLabel name, bool value, long int objectID)
{
    long int auxID;

    if (!objects.empty() && MetaDataContainer::isValidLabel(name))
    {
        if (objectID == -1)
        {
            auxID = objectsIterator->first;
        }
        else
        {
            auxID = objectID;
        }

        MetaDataContainer * aux = objects[auxID];

        // Check whether label is correct (belongs to the enum in the metadata_container header
        // and whether it is present in the activeLabels vector. If not, add it to all the other
        // objects with default values
        std::vector<MetaDataLabel>::iterator location;
        std::map<long int, MetaDataContainer *>::iterator It;

        location = std::find(activeLabels.begin(), activeLabels.end(), name);

        if (location == activeLabels.end())
        {
            activeLabels.push_back(name);

            // Add this label to the rest of the objects in this class
            for (It = objects.begin(); It != objects.end(); It++)
            {
                if (It->second != aux)
                {
                    (It->second)->addValue(name, bool());
                }
            }
        }

        aux->addValue(name, value);

        return true;
    }
    else
    {
        return false;
    }
}

bool MetaData::setValue(MetaDataLabel name, const std::string &value,
                        long int objectID)
{
    long int auxID;

    if (!objects.empty() && MetaDataContainer::isValidLabel(name))
    {
        if (objectID == -1)
        {
            auxID = objectsIterator->first;
        }
        else
        {
            auxID = objectID;
        }

        MetaDataContainer * aux = objects[auxID];

        // Check whether label is correct (belongs to the enum in the metadata_container header
        // and whether it is present in the activeLabels vector. If not, add it to all the other
        // objects with default values
        std::vector<MetaDataLabel>::iterator location;
        std::map<long int, MetaDataContainer *>::iterator It;

        location = std::find(activeLabels.begin(), activeLabels.end(), name);

        if (location == activeLabels.end())
        {
            activeLabels.push_back(name);

            // Add this label to the rest of the objects in this class
            for (It = objects.begin(); It != objects.end(); It++)
            {
                if (It->second != aux)
                {
                    (It->second)->addValue(name, std::string(""));
                }
            }
        }

        aux->addValue(name, value);

        return true;
    }
    else
    {
        return false;
    }
}

bool MetaData::setValue(MetaDataLabel name, const std::vector<double> &value,
                        long int objectID)
{
    long int auxID;

    if (!objects.empty() && MetaDataContainer::isValidLabel(name))
    {
        if (objectID == -1)
        {
            auxID = objectsIterator->first;
        }
        else
        {
            auxID = objectID;
        }

        MetaDataContainer * aux = objects[auxID];

        // Check whether label is correct (belongs to the enum in the metadata_container header
        // and whether it is present in the activeLabels vector. If not, add it to all the other
        // objects with default values
        std::vector<MetaDataLabel>::iterator location;
        std::map<long int, MetaDataContainer *>::iterator It;

        location = std::find(activeLabels.begin(), activeLabels.end(), name);

        if (location == activeLabels.end())
        {
            activeLabels.push_back(name);

            // Add this label to the rest of the objects in this class
            for (It = objects.begin(); It != objects.end(); It++)
            {
                if (It->second != aux)
                {
                    (It->second)->addValue(name, std::vector<double>());
                }
            }
        }

        aux->addValue(name, value);

        return true;
    }
    else
    {
        return false;
    }
}

bool MetaData::setValue(const std::string &name, const std::string &value,
                        long int objectID)
{
    long int auxID;

    MetaDataLabel label = MetaDataContainer::codifyLabel(name);

    if (!objects.empty() && MetaDataContainer::isValidLabel(label))
    {
        if (objectID == -1)
        {
            auxID = objectsIterator->first;
        }
        else
        {
            auxID = objectID;
        }

        MetaDataContainer * aux = objects[auxID];

        // Check whether label is correct (belongs to the enum in the metadata_container header
        // and whether it is present in the activeLabels vector. If not, add it to all the other
        // objects with default values
        std::vector<MetaDataLabel>::iterator location;
        std::map<long int, MetaDataContainer *>::iterator It;

        location = std::find(activeLabels.begin(), activeLabels.end(), label);

        if (location == activeLabels.end())
        {
            activeLabels.push_back(label);

            // Add this label to the rest of the objects in this class
            for (It = objects.begin(); It != objects.end(); It++)
            {
                if (It->second != aux)
                {
                    (It->second)->addValue(label, std::string(""));
                }
            }

        }

        aux->addValue(name, value);

        return true;
    }
    else
    {
        return false;
    }
}

std::vector<long int> & MetaData::getRandomOrderedObjects()
{
    return randomOrderedObjects;
}

bool MetaData::removeObject(long int objectID)
{
    int result = objects.erase(objectID);
    //randomOrderedObjects.erase(std::find(randomOrderedObjects.begin(),randomOrderedObjects.end(), objectID));
    objectsIterator = objects.begin();
    return result;
}

void MetaData::removeObjects(std::vector<long int> &toRemove)
{
    std::vector<long int>::iterator It;

    for (It = toRemove.begin(); It != toRemove.end(); It++)
    {
        objects.erase(*It);
        //randomOrderedObjects.erase(std::find(randomOrderedObjects.begin(),randomOrderedObjects.end(), *It));
    }

    // Reset iterator due to unknown iterator state after removing items
    objectsIterator = objects.begin();
}

void MetaData::removeObjects(MetaDataLabel name, double value)
{
    std::vector<long int> toRemove = findObjects(name, value);
    removeObjects(toRemove);
}

void MetaData::removeObjects(MetaDataLabel name, int value)
{
    std::vector<long int> toRemove = findObjects(name, value);
    removeObjects(toRemove);
}

void MetaData::removeObjects(MetaDataLabel name, bool value)
{
    std::vector<long int> toRemove = findObjects(name, value);
    removeObjects(toRemove);
}

void MetaData::removeObjects(MetaDataLabel name, std::string value)
{
    std::vector<long int> toRemove = findObjects(name, value);
    removeObjects(toRemove);
}

std::vector<long int> MetaData::findObjectsInRange(MetaDataLabel name,
        double minValue, double maxValue)
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

            std::cerr << "VAL_D: " << value << std::endl;
            if (value >= minValue && value <= maxValue)
            {
                result.push_back(It->first);
            }
        }
    }

    return result;
}

std::vector<long int> MetaData::findObjectsInRange(MetaDataLabel name,
        int minValue, int maxValue)
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
            int value;
            aux->getValue(name, value);
            std::cerr << "VAL_I: " << value << std::endl;

            if (value >= minValue && value <= maxValue)
            {
                result.push_back(It->first);
            }
        }
    }

    return result;

}

std::vector<long int> MetaData::findObjects(MetaDataLabel name, double value)
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

std::vector<long int> MetaData::findObjects(MetaDataLabel name, int value)
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

std::vector<long int> MetaData::findObjects(MetaDataLabel name, bool value)
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

std::vector<long int> MetaData::findObjects(MetaDataLabel name,
        std::string value)
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

bool MetaData::detectObjects(MetaDataLabel name, int value)
{
    bool result = false;
    // Traverse all the structure looking for objects
    // that satisfy search criteria
    MetaDataContainer * aux;
    std::map<long int, MetaDataContainer *>::iterator It;
    for (It = objects.begin(); It != objects.end(); It++)
    {
        aux = It->second;
        if (aux->pairExists(name, value))
        {
            result = true;
            break;
        }
    }
    return result;
}

long int MetaData::countObjects(MetaDataLabel name, int value)
{
    return ((findObjects(name, value)).size());
}

long int MetaData::countObjects(MetaDataLabel name, double value)
{
    return ((findObjects(name, value)).size());
}

long int MetaData::countObjects(MetaDataLabel name, bool value)
{
    return ((findObjects(name, value)).size());
}

long int MetaData::countObjects(MetaDataLabel name, const std::string &value)
{
    return ((findObjects(name, value)).size());
}

bool MetaData::getValue(MetaDataLabel name, double &value, long int objectID)
{
    MetaDataContainer * aux = getObject(objectID);
    if( !aux->valueExists(name) )
    {
        return false;
    }
    else
    {
        aux->getValue(name, value);
    }

    return true;
}

bool MetaData::getValue(MetaDataLabel name, int &value, long int objectID)
{
    MetaDataContainer * aux = getObject(objectID);
    if( !aux->valueExists(name) )
    {
        return false;
    }
    else
    {
        aux->getValue(name, value);
    }

    return true;
}

bool MetaData::getValue(MetaDataLabel name, bool &value, long int objectID)
{
    MetaDataContainer * aux = getObject(objectID);
    if( !aux->valueExists(name) )
    {
        return false;
    }
    else
    {
        aux->getValue(name, value);
    }

    return true;
}

bool MetaData::getValue(MetaDataLabel name, std::vector<double> &value,
                        long int objectID)
{
    MetaDataContainer * aux = getObject(objectID);
    if( !aux->valueExists(name) )
    {
        return false;
    }
    else
    {
        aux->getValue(name, value);
    }

    return true;
}

bool MetaData::getValue(MetaDataLabel name, std::string &value,
                        long int objectID)
{
    MetaDataContainer * aux = getObject(objectID);
    if( !aux->valueExists(name) )
    {
        return false;
    }
    else
    {
        aux->getValue(name, value);
    }

    return true;
}

MetaDataContainer * MetaData::getObject(long int objectID)
{
    if (isEmpty())
    {
        // The objects map is empty, error
        REPORT_ERROR( -1, "Requested objecID not found (no objects stored). Exiting... ");
    }

    MetaDataContainer * aux;

    if (objectID == -1)
        aux = objectsIterator->second;
    else
        aux = objects[objectID];

    if (aux == NULL)
    {
        // This objectID does not exist, finish execution
        REPORT_ERROR( -1, "Requested objecID not found. Exiting... ");
    }

    return aux;
}

long int MetaData::fastSearch(MetaDataLabel name, std::string value,
                              bool recompute)
{
    long int result;

    if (recompute || fastStringSearch.empty() || fastStringSearchLabel != name)
    {
        fastStringSearch.clear();
        fastStringSearchLabel = name;

        // Repopulate list
        // Traverse all the structure looking for objects
        // that satisfy search criteria

        std::map<long int, MetaDataContainer *>::iterator It;

        MetaDataContainer * aux;

        for (It = objects.begin(); It != objects.end(); It++)
        {
            aux = It->second;

            if (aux->valueExists(name))
            {
                std::string auxStr;
                (It->second)->getValue(name, auxStr);
                fastStringSearch[auxStr] = It->first;

                if (aux->pairExists(name, value))
                {
                    result = It->first;
                }
            }
        }
    }
    else
    {
        std::map<std::string, long int>::iterator It;

        if ((It = fastStringSearch.find(value)) != fastStringSearch.end())
        {
            result = It->second;
        }
        else
        {
            result = -1; // Not found
        }
    }

    return result;
}

/* Statistics -------------------------------------------------------------- */
void get_statistics(MetaData MT_in, Image& _ave, Image& _sd, double& _min,
                    double& _max, bool apply_geo)
{
    MetaData MT(MT_in); //copy constructor so original MT is not changed
    _min = MAXFLOAT;
    _max = 0;
    bool first = true;
    int n = 0;
    // Calculate Mean
    long int ret = MT.firstObject();
    if (ret == MetaData::NO_OBJECTS_STORED)
    {
        std::cerr << "Empty inputFile File\n";
        exit(1);
    }
    FileName image_name;
    int _enabled;
    do
    {
        MT.getValue(MDL_IMAGE, image_name);
        MT.getValue(MDL_ENABLED, _enabled);
        if (_enabled == (-1) || image_name == "")
            continue;
        Image *image = Image::LoadImage(image_name, apply_geo); // reads image
        double min, max, avg, stddev;
        (*image)().computeStats(avg, stddev, min, max);
        if (_min > min)
            _min = min;
        if (_max < max)
            _max = max;
        if (first)
        {
            _ave = *image;
            first = false;
        }
        else
        {
            _ave() += (*image)();
        }
        delete image;
        n++;
    }
    while (MT.nextObject() != MetaData::NO_MORE_OBJECTS);

    if (n > 0)
        _ave() /= n;
    _sd = _ave;
    _sd().initZeros();
    // Calculate SD
    MT.firstObject();
    do
    {
        MT.getValue(MDL_IMAGE, image_name);
        MT.getValue(MDL_ENABLED, _enabled);
        if (_enabled == (-1) || image_name == "")
            continue;
        Image *image = Image::LoadImage(image_name, apply_geo); // reads image
        Image tmpImg;
        tmpImg() = (((*image)() - _ave()));
        tmpImg() *= tmpImg();
        _sd() += tmpImg();
        delete image;
    }
    while (MT.nextObject() != MetaData::NO_MORE_OBJECTS);
    _sd() /= (n - 1);
    _sd().selfSQRTnD();
}

void ImgSize(int &Ydim, int &Xdim)
{
    /*
     std::vector<SelLine>::iterator aux = current_line;
     go_first_ACTIVE();
     FileName fn_img = (*current_line).text;
     if (fn_img.find("imagic:") != -1)
     {
     Image *img = Image::LoadImage(fn_img);
     Ydim = (*img)().ydim;
     Xdim = (*img)().xdim;
     delete img;
     }
     else if (Is_ImageXmipp(fn_img))
     {
     ImageXmipp img;
     img.read(fn_img);
     Ydim = img().ydim;
     Xdim = img().xdim;
     }
     else if (Is_FourierImageXmipp(fn_img))
     {
     FourierImageXmipp img;
     img.read(fn_img);
     Ydim = img().ydim;
     Xdim = img().xdim;
     }
     else
     REPORT_ERROR(1, "SelFile::ImgSize: First Active file is not an image");

     current_line = aux;
     */
}
