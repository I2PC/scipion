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

MetaData::MetaData(const MetaData &MD)
{
    clear();
    this->setComment(MD.getComment());
    this->setPath(MD.getPath());
    this->isColumnFormat = MD.isColumnFormat;
    this->inFile = MD.inFile;
    this->fastStringSearchLabel = MDL_UNDEFINED;
    this->activeLabels = MD.activeLabels;

    //objects, define iterator
    std::map<long int, MetaDataContainer *>::const_iterator objIt;
    for (objIt = MD.objects.begin(); objIt != MD.objects.end(); objIt++)
    {
        long int idx = this->addObject();

        this->objects[idx] = new MetaDataContainer(*(objIt->second));
    }
    this->objectsIterator = objects.begin();

}

void MetaData::fillWithNextNObjects (MetaData &MD, long int start, long int numberObjects)
{
    clear();
    this->setComment(MD.getComment());
    this->setPath(MD.getPath());
    this->isColumnFormat = MD.isColumnFormat;
    this->inFile = MD.inFile;
    this->fastStringSearchLabel = MDL_UNDEFINED;
    this->activeLabels = MD.activeLabels;

    //objects, define iterator
    int counter=0;
    std::map<long int, MetaDataContainer *>::iterator objIt;
    objIt = MD.objects.begin();
    advance( objIt, start );
    for (; objIt != MD.objects.end() && counter<numberObjects; objIt++,counter++)
    {
        long int idx = this->addObject(objIt->first);
        this->objects[idx] = new MetaDataContainer(*(objIt->second));
    }
    this->objectsIterator = objects.begin();
}

MetaData& MetaData::operator =(const MetaData &MD)
{
    if (this != &MD)
    {
        clear();
        this->setComment(MD.getComment());
        this->setPath(MD.getPath());
        this->isColumnFormat = MD.isColumnFormat;
        this->inFile = MD.inFile;
        this->fastStringSearchLabel = MDL_UNDEFINED;
        this->activeLabels = MD.activeLabels;

        //objects, define iterator
        std::map<long int, MetaDataContainer *>::const_iterator objIt;
        for (objIt = MD.objects.begin(); objIt != MD.objects.end(); objIt++)
        {
            long int idx = this->addObject();

            this->objects[idx] = new MetaDataContainer(*(objIt->second));
        }
        this->objectsIterator = objects.begin();
    }
    return *this;
}

void MetaData::union_(MetaData &MD, MDLabel thisLabel)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be added" );
    }

    if (this->size() == 0)
    {
        *this = MD;
        return;
    }
    else if (MD.size()==0)
        return;

    //this->activeLabels = minuend.activeLabels;
    MetaDataContainer * aux, *aux2;
    std::string value1, value2;
    std::map<long int, MetaDataContainer *>::iterator itMinuend;
    bool add;
    for (itMinuend = MD.objects.begin(); itMinuend
         != MD.objects.end(); itMinuend++)
    {
        aux = MD.getObject(itMinuend->first);
        aux->writeValueToString(value1, thisLabel);
        add=true;
        for (long int idSubtrahendr = this->firstObject(); idSubtrahendr
             != NO_MORE_OBJECTS; idSubtrahendr = this->nextObject())
        {
            aux2 = this->getObject(idSubtrahendr);
            aux2->writeValueToString(value2, thisLabel);

            if (value2 == value1)
            {
                add=false;
                break;
            }
        }
        if(add)
        {
            long int idx = this->addObject();
            this->objects[idx]
            = new MetaDataContainer(*(itMinuend->second));
        }
    }
    this->objectsIterator = this->objects.begin();
}

void MetaData::unionAll(MetaData &MD)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be added" );
    }

    if (this->size() == 0)
    {
        *this = MD;
        return;
    }
    else if (MD.size()==0)
        return;

    std::map<long int, MetaDataContainer *>::iterator itMinuend;
    for (itMinuend = MD.objects.begin(); itMinuend
         != MD.objects.end(); itMinuend++)
    {
        long int idx = this->addObject();
        this->objects[idx]
        = new MetaDataContainer(*(itMinuend->second));
    }
    this->objectsIterator = this->objects.begin();
}

void MetaData::merge(const FileName &fn)
{
    MetaData aux;
    aux.read(fn);
    union_(aux);
}

void MetaData::intersection(MetaData &minuend, MetaData &subtrahend,
                            MDLabel thisLabel)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be substracted" );
    }

    if (minuend.size() == 0)
    {
        *this = minuend;
        return;
    }
    if (subtrahend.size() == 0)
    {
        *this = subtrahend;
        return;
    }

    this->activeLabels = minuend.activeLabels;
    MetaDataContainer * aux, *aux2;
    std::string value1, value2;
    std::map<long int, MetaDataContainer *>::iterator itMinuend;
    for (itMinuend = minuend.objects.begin(); itMinuend
         != minuend.objects.end(); itMinuend++)
    {
        aux = minuend.getObject(itMinuend->first);
        aux->writeValueToString(value1, thisLabel);
        for (long int idSubtrahendr = subtrahend.firstObject(); idSubtrahendr
             != NO_MORE_OBJECTS; idSubtrahendr = subtrahend.nextObject())
        {
            aux2 = subtrahend.getObject(idSubtrahendr);
            aux2->writeValueToString(value2, thisLabel);

            if (value2 == value1)
            {
                long int idx = this->addObject();
                this->objects[idx]
                = new MetaDataContainer(*(itMinuend->second));
                break;
            }
        }
    }
    this->objectsIterator = this->objects.begin();
}

void MetaData::substraction(MetaData &minuend, MetaData &subtrahend,
                            MDLabel thisLabel)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be substracted" );
    }

    if (minuend.size() == 0 || subtrahend.size() == 0)
    {
        *this = minuend;
        return;
    }

    this->activeLabels = minuend.activeLabels;
    MetaDataContainer * aux, *aux2;
    std::string value1, value2;
    std::map<long int, MetaDataContainer *>::iterator itMinuend;
    for (itMinuend = minuend.objects.begin(); itMinuend
         != minuend.objects.end(); itMinuend++)
    {
        aux = minuend.getObject(itMinuend->first);
        aux->writeValueToString(value1, thisLabel);
        bool doSave = false;
        for (long int idSubtrahendr = subtrahend.firstObject(); idSubtrahendr
             != NO_MORE_OBJECTS; idSubtrahendr = subtrahend.nextObject())
        {
            aux2 = subtrahend.getObject(idSubtrahendr);
            aux2->writeValueToString(value2, thisLabel);

            if (value2 == value1)
            {
                doSave = false;
                break;
            }
            else
            {
                doSave = true;
            }
        }
        if (doSave)
        {
            long int idx = this->addObject();
            this->objects[idx] = new MetaDataContainer(*(itMinuend->second));
        }
    }
    this->objectsIterator = this->objects.begin();
}

void MetaData::read(std::ifstream *infile,
                    std::vector<MDLabel> * labelsVector)
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
            MDLabel label = MDL::str2Label(newLabel);

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

                setValue(MDL::label2Str(
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
                setValue(newLabel, value);
            }
        }
    }
}

void MetaData::readOldDocFile(std::ifstream *infile,
                              std::vector<MDLabel> * labelsVector)
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
        std::vector<MDLabel>::iterator location;

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
            newLabel.erase(
                pos, 1);

        // Remove unneded parentheses and contents
        pos = newLabel.find("(");
        newLabel.erase(pos, newLabel.find(")") - pos + 1);

        // Remove white spaces
        newLabel = removeChar(newLabel, ' ');

        if (labelsVector != NULL)
        {
            std::vector<MDLabel>::iterator location;

            location = std::find(labelsVector->begin(), labelsVector->end(),
                                 MDL::str2Label(newLabel));

            if (location != labelsVector->end())
            {
                activeLabels.push_back(MDL::str2Label(newLabel));
            }
        }
        else
        {
            activeLabels.push_back(MDL::str2Label(newLabel));
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
                        setValue(MDL::label2Str(
                                     activeLabels[counter - 1]), value);
                    else
                        setValue(MDL::label2Str(
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

MetaData::MetaData(FileName fileName, std::vector<MDLabel> * labelsVector)
{
    read(fileName, labelsVector);
}

void MetaData::setColumnFormat(bool column)
{
    isColumnFormat = column;
}

void MetaData::toDataBase(CppSQLite3DB &db,
                          const FileName & DBname,
                          const std::string & tableName,
                          const std::vector<MDLabel> * labelsVector,
                          bool OpenDb,
                          bool CloseDb
                         )
{
    try
    {
        if(OpenDb)
            db.open(DBname, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);

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
        std::vector<MDLabel>::iterator strIt;
        std::string sqlCommandInsert = "insert into " + tableName + "(objId,";
        for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            if (labelsVector == NULL || std::find(labelsVector->begin(),
                                                  labelsVector->end(), *strIt) != labelsVector->end())
            {
                std::string label = MDL::label2Str(*strIt);
                sqlCommandInsert += label + ",";
                if (MDL::isDouble(*strIt))
                {
                    sqlCommandCreate += label + " double,";
                    labelsCounter++;
                }
                else if (MDL::isString(*strIt))
                {
                    sqlCommandCreate += label + " text,";
                    labelsCounter++;
                }
                else if (MDL::isInt(*strIt))
                {
                    sqlCommandCreate += label + " int,";
                    labelsCounter++;
                }
                else if (MDL::isBool(*strIt))
                {
                    sqlCommandCreate += label + " int,";
                    labelsCounter++;
                }
                else if (MDL::isVector(*strIt))
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
        //std::vector<MDLabel>::iterator strIt;
        CppSQLite3Statement stmt = db.compileStatement(sqlCommandInsert);

        for (long int IDthis = firstObject(); IDthis != NO_MORE_OBJECTS; IDthis
             = nextObject())
        {
            stmt.bind(1, (int) IDthis);
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
            stmt.bind(labelsCounter, (int) IDthis);
            for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
            {
                if (labelsVector == NULL || std::find(labelsVector->begin(),
                                                      labelsVector->end(), *strIt) != labelsVector->end())
                {
                    if (MDL::isDouble(*strIt))
                    {
                        labelsCounter++;
                        double value;
                        getValue(*strIt, value);
                        stmt.bind(labelsCounter, value);
                    }
                    else if (MDL::isString(*strIt))
                    {
                        labelsCounter++;
                        std::string value;
                        getValue(*strIt, value);
                        stmt.bind(labelsCounter, value.c_str());
                    }
                    else if (MDL::isInt(*strIt))
                    {
                        labelsCounter++;
                        int value;
                        getValue(*strIt, value);
                        stmt.bind(labelsCounter, value);
                    }
                    else if (MDL::isBool(*strIt))
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
        if(CloseDb)
            db.close();
    }
    catch (CppSQLite3Exception& e)
    {
        REPORT_ERROR( e.errorCode(), e.errorMessage());
    }
}
void MetaData::fromDataBase(CppSQLite3DB &db,
                            const FileName & DBname,
                            const std::string & tableName,
                            const MDLabel sortLabel,
                            std::vector<MDLabel> * labelsVector,
                            bool OpenDb,
                            bool CloseDb
                           )
{
    //check that database and table exits
    try
    {
        std::ostringstream stm;
        std::vector<MDLabel> _labelsVector;
        std::string sqlCommand;
        std::vector<MDLabel>::iterator strIt;
        if(OpenDb)
        {
            if (!exists(DBname))
                REPORT_ERROR( 1, (std::string)"database "+ DBname + " does not exists");
            db.open(DBname, SQLITE_OPEN_READONLY);
        }
        if (!db.tableExists(tableName))
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
        if (labelsVector == NULL)
        {
            sqlCommand = (std::string) "PRAGMA table_info(" + tableName + ");";
            CppSQLite3Query q = db.execQuery(sqlCommand);
            while (!q.eof())
            {
                _labelsVector.push_back(MDL::str2Label(
                                            q.fieldValue(1)));
                q.nextRow();
            }
        }
        else
        {
            _labelsVector.push_back(MDL_OBJID);
            for (strIt = (*labelsVector).begin(); strIt
                 != (*labelsVector).end(); strIt++)
            {
                _labelsVector.push_back(*strIt);
            }
        }

        //Order by...
        strIt = std::find(_labelsVector.begin(), _labelsVector.end(), sortLabel);
        int _order=1;
        if (strIt == _labelsVector.end())
            REPORT_ERROR(1,"label" + MDL::label2Str(*strIt)+ "does not exists");
        else
            _order += std::distance(_labelsVector.begin(), strIt);

        sqlCommand = "select ";
        for (strIt = (_labelsVector).begin(); strIt != (_labelsVector).end(); strIt++)
        {
            sqlCommand += MDL::label2Str(*strIt) + ",";
        }
        sqlCommand.erase(sqlCommand.end() - 1);
        sqlCommand += (std::string) " from " + tableName;
        stm << _order;
        sqlCommand += " order by " + stm.str(); //sort always by objId

        CppSQLite3Query q = db.execQuery(sqlCommand);
        int labelsCounter = 0;

        int _objID=-1;
        while (!q.eof())
        {
            labelsCounter = 0;
            if(_order==1)
                _objID = q.getIntField(labelsCounter);
            else
                _objID++;
            addObject(_objID);
            for (strIt = (_labelsVector).begin() + 1; strIt
                 != (_labelsVector).end(); strIt++)
            {
                if (MDL::isDouble(*strIt))
                {
                    labelsCounter++;
                    double value;
                    value = q.getFloatField(labelsCounter);
                    setValue(*strIt, value);
                }
                else if (MDL::isInt(*strIt))
                {
                    labelsCounter++;
                    int value;
                    value = q.getIntField(labelsCounter);
                    setValue(*strIt, value);
                }
                else if (MDL::isString(*strIt))
                {
                    labelsCounter++;
                    std::string value;
                    value.assign(q.getStringField(labelsCounter));
                    setValue(*strIt, value);
                }
                else if (MDL::isBool(*strIt))
                {
                    labelsCounter++;
                    int value;
                    value = q.getIntField(labelsCounter);
                    bool value2 = (bool) value;
                    setValue(*strIt, value2);
                }
                else
                {
                    REPORT_ERROR( 1, (std::string)"unknown type in fromDataBase routine");
                }
            }
            q.nextRow();
        }

        //        CppSQLite3Table t = db.getTable(sqlCommand);
        if(CloseDb)
            db.close();
    }
    catch (CppSQLite3Exception& e)
    {
        REPORT_ERROR( e.errorCode(), e.errorMessage());
    }

}
void MetaData::combineWithFiles(MDLabel thisLabel)
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

        for (MDLabel mdl = MDL_FIRST_LABEL; mdl <= MDL_LAST_LABEL; mdl
             = MDLabel(mdl + 1))
        {
            if (auxMetaDataContainer->valueExists(mdl))
            {
                std::string value;

                auxMetaDataContainer->writeValueToString(value, mdl);

                setValue(MDL::label2Str(mdl), value);
            }
        }
    }
}
void MetaData::combine(MetaData & other, MDLabel thisLabel)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be combined" );
    }

    MetaDataContainer * aux, *aux2;
    std::string value1, value2;
    std::string value;

    std::vector<MDLabel>::iterator strIt;
    //init everything with zeros, this should not be needed by write is not properlly writting
    //FIXIT
    for (long int IDthis = firstObject(); IDthis != NO_MORE_OBJECTS; IDthis
         = nextObject())
    {
        aux = getObject();
        for (strIt = other.activeLabels.begin(); strIt != other.activeLabels.end(); strIt++)
        {
            MDLabel mdl = *strIt;
            value="0";
            if(thisLabel!=mdl)
                setValue(MDL::label2Str(mdl), value);
        }
    }
    bool notFound;
    for (long int IDthis = firstObject(); IDthis != NO_MORE_OBJECTS; IDthis
         = nextObject())
    {
        aux = getObject();
        aux->writeValueToString(value1, thisLabel);
        notFound=true;
        for (long int IDother = other.firstObject(); IDother != NO_MORE_OBJECTS; IDother
             = other.nextObject())
        {
            aux2 = other.getObject();
            aux2->writeValueToString(value2, thisLabel);

            if (value2 == value1)
            {
                for (strIt = other.activeLabels.begin(); strIt != other.activeLabels.end(); strIt++)
                {
                    MDLabel mdl = *strIt;
                    aux2->writeValueToString(value, mdl);
                    setValue(MDL::label2Str(mdl), value);
                }
                notFound=false;
                break;
            }
        }
        if(notFound)
            std::cerr << "WARNING: no match found for label '"
            << MDL::label2Str(thisLabel)
            << "' entry='" << value1
            << "' unknown values filled with 0" << std::endl;
    }
}

void MetaData::read(FileName fileName,
                    std::vector<MDLabel> * labelsVector)
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
    }
    else
    {
        pos = line.find("XMIPP_3 * ");

        if (pos != std::string::npos) // xmipp_3 token found
        {
            read(&infile, labelsVector);
        }
        else // We are reading an old selfile
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
    aux->writeValueToString(result, MDL::str2Label(inputLabel));
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
    std::vector<MDLabel>::iterator strIt;

    if (isColumnFormat)
    {
        outfile << "; ";
        for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            if (*strIt != MDL_COMMENT)
            {
                outfile << MDL::label2Str(*strIt);
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
                    outfile << "; "
                    << entryComment << std::endl;
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

        int maxWidth=20;
        for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            if (*strIt != MDL_COMMENT)
            {
                int w=MDL::label2Str(*strIt).length();
                if (w>maxWidth)
                    maxWidth=w;
            }
        }

        for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            if (*strIt != MDL_COMMENT)
            {
                outfile.width(maxWidth+1);
                outfile << MDL::label2Str(*strIt);
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

bool MetaData::isEmpty() const
{
    return objects.empty();
}

void MetaData::clear()
{
    //By default, random_shuffle uses rand, so do NOT remove srand next line
    srand(time(NULL));
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

std::string MetaData::getPath() const
{
    return path;
}

void MetaData::setComment(const std::string newComment)
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

std::string MetaData::getComment() const
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
    std::vector<MDLabel>::iterator It;
    for (It = activeLabels.begin(); It != activeLabels.end(); It++)
    {
        (objectsIterator->second)->addValue(
            MDL::label2Str(*It), std::string(""));
    }

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

bool MetaData::setValue(const std::string &name, const std::string &value,
                        long int objectID)
{
    long int auxID;

    MDLabel label = MDL::str2Label(name);

    if (!objects.empty() && MDL::isValidLabel(label))
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
        std::vector<MDLabel>::iterator location;
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

bool MetaData::removeObject(long int objectID)
{
    int result = objects.erase(objectID);
    objectsIterator = objects.begin();
    return result;
}

void MetaData::removeObjects(std::vector<long int> &toRemove)
{
    std::vector<long int>::iterator It;

    for (It = toRemove.begin(); It != toRemove.end(); It++)
    {
        objects.erase(*It);
    }

    // Reset iterator due to unknown iterator state after removing items
    objectsIterator = objects.begin();
}


void MetaData::fillMetaData(MetaData &MD, std::vector<long int> objectsToAdd)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be added" );
    }

    this->activeLabels = MD.activeLabels;
    for (int i = 0; i < objectsToAdd.size(); i++)
    {
        long int idx = this->addObject();
        this->objects[idx] = new MetaDataContainer(*(MD.objects[objectsToAdd[i]]));
    }
    this->objectsIterator = this->objects.begin();
}

void MetaData::importObjects(MetaData &base, int first, int last)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData can not be added" );
    }

    activeLabels = base.activeLabels;
    if (base.size()==0)
    {
        clear();
        return;
    }

    first=XMIPP_MIN(first,base.size());
    last=XMIPP_MIN(last,base.size());
    for (int i = first; i <=last; i++)
    {
        long int idx = addObject();
        this->objects[idx] = new MetaDataContainer(*(base.objects[i]));
    }
    objectsIterator = objects.begin();
}

MetaDataContainer * MetaData::getObject(const long int objectID) const
{
    if (isEmpty())
    {
        // The objects map is empty, error
        REPORT_ERROR( -1, "Requested objectID not found (no objects stored). Exiting... ");
    }

    MetaDataContainer * aux;

    if (objectID == -1)
        aux = objectsIterator->second;
    else
        aux = objects.at(objectID);

    if (aux == NULL)
    {
        // This objectID does not exist, finish execution
        REPORT_ERROR( -1, "Requested objectID not found. Exiting... ");
    }

    return aux;
}

// Allows a fast search for pairs where the value is
// a string, i.e. looking for filenames which is quite
// usual
long int MetaData::fastSearch(MDLabel name, std::string value,
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
//Randomize
void MetaData::randomize(MetaData &MDin)
{
    std::vector<long int> v;
    //v.reserve(objects.size());
    std::map<long int, MetaDataContainer *>::iterator it;
    for (it = MDin.objects.begin(); it
         != MDin.objects.end(); it++)
    {
        v.push_back(it->first);
    }
    random_shuffle(v.begin(), v.end());
    int i=0;

    this->activeLabels = MDin.activeLabels;

    for (it = MDin.objects.begin(); it != MDin.objects.end(); it++)
    {
        long int idx = this->addObject(v[i++]);
        objects[idx] = new MetaDataContainer(*(it->second));
    }
}

/* Statistics -------------------------------------------------------------- */
#include "metadata_extension.h"
void get_statistics(MetaData MT_in, Image<double> & _ave, Image<double> & _sd, double& _min,
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
        Image<double> image;
        MetaDataContainer * aux = MT.getObject();
        image.read(image_name,true,
                   -1,apply_geo,
                   false);
        double min, max, avg, stddev;
        image().computeStats(avg, stddev, min, max);
        if (_min > min)
            _min = min;
        if (_max < max)
            _max = max;
        if (first)
        {
            _ave = image;
            first = false;
        }
        else
        {
            _ave() += image();
        }
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

        Image<double> image, tmpImg;
        image.read(image_name);
        tmpImg() = ((image() - _ave()));
        tmpImg() *= tmpImg();
        _sd() += tmpImg();
    }
    while (MT.nextObject() != MetaData::NO_MORE_OBJECTS);
    _sd() /= (n - 1);
    _sd().selfSQRT();
}

void ImgSize(MetaData &MD, int &Xdim, int &Ydim, int &Zdim, int &Ndim)
{
    if (MD.firstObject() != MetaData::NO_OBJECTS_STORED)
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

void ImgSize(FileName &fn_img, int &Xdim, int &Ydim, int &Zdim, int &Ndim)
{
    Image<double> img;
    img.read(fn_img, false);
    img.getDimensions(Xdim, Ydim, Zdim, Ndim);
}

void MetaData::split_in_two(MetaData &SF1, MetaData &SF2,MDLabel label)
{
    if (!isColumnFormat)
    {
        REPORT_ERROR( -1, "Row formatted MetaData cannot be splitted in two" );
    }

    MetaData SFtmp;
    SFtmp.randomize(*this);
    SFtmp.firstObject();

    SF1.clear();
    SF2.clear();
    if(inFile.length()<1)
    {
        SF1.setComment("Splitted in half from file " + inFile);
        SF2.setComment("Splited in half from file " + inFile);
    }
    SF2.isColumnFormat = SF1.isColumnFormat = SFtmp.isColumnFormat;
    SF2.fastStringSearchLabel =SF1.fastStringSearchLabel = MDL_UNDEFINED;
    SF2.activeLabels =SF1.activeLabels = SFtmp.activeLabels;
    SF2.isColumnFormat = SFtmp.isColumnFormat;

    //objects, define iterator
    std::map<long int, MetaDataContainer *>::iterator objIt;
    for (objIt = SFtmp.objects.begin(); objIt != SFtmp.objects.end(); objIt++)
    {
        long int idx = SF1.addObject();
        SF1.objects[idx] = new MetaDataContainer(*(objIt->second));
        objIt++;
        if(objIt == SFtmp.objects.end())
            break;
        else
        {
            idx = SF2.addObject();
            SF2.objects[idx] = new MetaDataContainer(*(objIt->second));
        }
    }
    if(label!=MDL_UNDEFINED)
    {
        SF1.sort(SF1,label);
        SF2.sort(SF2,label);
    }
}

void MetaData::mpi_select_part(int rank, int size, int &num_img_tot)
{
    num_img_tot = (*this).size();
    int remaining = num_img_tot % size;
    int Npart = (int)(num_img_tot - remaining) / size;
    int myFirst, myLast;
    if (rank < remaining)
    {
        myFirst = rank * (Npart + 1);
        myLast = myFirst + Npart;
    }
    else
    {
        myFirst = rank * Npart + remaining;
        myLast = myFirst + Npart - 1;
    }

    // Now discard all images in Selfile that are outside myFirst-myLast
    for (int nr = num_img_tot-1; nr >=0; nr--)
    {
        if (nr<myFirst || nr>myLast)
            removeObject(nr);
    }
}

void MetaData::sort(MetaData & MDin, MDLabel sortlabel)
{
    CppSQLite3DB db;
    if(sortlabel==MDL_UNDEFINED)
        return;
    std::string tempTable;
    tempTable="tempTable";
    //"" means temporary in memory database
    MDin.toDataBase( db,"",tempTable,NULL,true,false);
    this->clear();
    this->fromDataBase( db,"", tempTable,sortlabel,NULL,false,true);
    db.close();
    //    if(remove ( tmpFileName ) == -1)
    //        std::cerr << "cannot remove file " << tmpFileName <<std::endl;
}

void MetaData::aggregate(MetaData MDIn,
                         MDLabel aggregateLabel,
                         MDLabel entryLabel,
                         MDLabel operationLabel)
{
    CppSQLite3DB db;
    if(aggregateLabel ==MDL_UNDEFINED ||
       operationLabel ==MDL_UNDEFINED ||
       entryLabel     ==MDL_UNDEFINED)
        return;
    //check is valid operation
    if(operationLabel  != MDL_AVG &&
       operationLabel != MDL_COUNT&&
       operationLabel != MDL_MAX &&
       operationLabel != MDL_MIN &&
       operationLabel != MDL_SUM)
        REPORT_ERROR(1,"invalid operation:" +MDL::label2Str(operationLabel) +" in aggregate");
    std::string tempTable;
    tempTable="tempTable";
    this->clear();
    try
    {
        std::vector<MDLabel> _labelsVector;
        std::string sqlCommand;
        std::vector<MDLabel>::iterator strIt;
        //convert to database in memory
        MDIn.toDataBase(db,"",tempTable,NULL,true,false);
        _labelsVector.push_back(aggregateLabel);
        _labelsVector.push_back(operationLabel);

        //Order by...

        sqlCommand  = (std::string)"select ";
        sqlCommand += MDL::label2Str(aggregateLabel);
        sqlCommand += (std::string)", ";
        if(operationLabel==MDL_AVG)
            sqlCommand += "avg(";
        else if(operationLabel==MDL_COUNT)
            sqlCommand += "count(";
        else if(operationLabel==MDL_MAX)
            sqlCommand += "max(";
        else if(operationLabel==MDL_MIN)
            sqlCommand += "min(";
        else if(operationLabel==MDL_SUM)
            sqlCommand += "total(";
        else
            REPORT_ERROR(1,"invalid operation selected in MetaData::aggregate");

        sqlCommand += MDL::label2Str(entryLabel);
        sqlCommand += (std::string)") ";
        sqlCommand += " as " + MDL::label2Str(operationLabel);
        sqlCommand += (std::string) " from " + tempTable;
        sqlCommand += " group by ";
        sqlCommand += MDL::label2Str(aggregateLabel) + " ";
        sqlCommand += " order by " ;
        sqlCommand += MDL::label2Str(aggregateLabel) +";";

        //std::cerr << "sqlCommand: " << sqlCommand <<std::endl;

        CppSQLite3Query q = db.execQuery(sqlCommand);
        int labelsCounter;

        while (!q.eof())
        {
            labelsCounter = -1;
            addObject();
            for (strIt = (_labelsVector).begin() ; strIt
                 != (_labelsVector).end(); strIt++)
            {
                if (MDL::isDouble(*strIt))
                {
                    labelsCounter++;
                    double value;
                    value = q.getFloatField(labelsCounter);
                    setValue(*strIt, value);
                }
                else if (MDL::isInt(*strIt))
                {
                    labelsCounter++;
                    int value;
                    value = q.getIntField(labelsCounter);
                    setValue(*strIt, value);
                }
                else if (MDL::isString(*strIt))
                {
                    labelsCounter++;
                    std::string value;
                    value.assign(q.getStringField(labelsCounter));
                    setValue(*strIt, value);
                }
                else if (MDL::isBool(*strIt))
                {
                    labelsCounter++;
                    int value;
                    value = q.getIntField(labelsCounter);
                    bool value2 = (bool) value;
                    setValue(*strIt, value2);
                }
                else
                {
                    REPORT_ERROR( 1, (std::string)"unknown type in fromDataBase routine");
                }
            }
            q.nextRow();
        }
        db.close();
    }
    catch (CppSQLite3Exception& e)
    {
        REPORT_ERROR( e.errorCode(), e.errorMessage());
    }
}
int MetaData::MaxStringLength( MDLabel thisLabel)
{
    if (!MDL::isString(thisLabel))
    {
        REPORT_ERROR(1,"MaxFileNameLength only works for strings");
    }
    std::string ss;
    int max_length =0,length;
    std::map<long int, MetaDataContainer *>::iterator it;
    for (it = objects.begin(); it != objects.end(); it++)
    {
        (it->second)->getValue(thisLabel, ss);
        length = ss.length();
        if ( length > max_length)
            max_length = length;
    }
    return max_length;
}

