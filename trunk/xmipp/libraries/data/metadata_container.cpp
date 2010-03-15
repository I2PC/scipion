/***************************************************************************
* 
* Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
*              R. Marabini        (roberto@cnb.csic.es)
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
#include "metadata_container.h"
//called by creator
void metaDataContainer::open(const FileName baseFileName,int flag)
{
	int nRet = sqlite3_open_v2(baseFileName.c_str(), &dB,flag, NULL);
//#define DEBUG_OPEN
#ifdef DEBUG_OPEN
    std::cerr<< "open database " <<  baseFileName
    		 << "with flags    " << flag
    		 << "result is     " << nRet
    		 << std::endl;

#endif
	if (nRet != SQLITE_OK)
	{
		REPORT_ERROR(1,sqlite3_errmsg(dB));
	}

	setBusyTimeout(mnBusyTimeoutMs);
	char *errmsg;
	//some optimization, no too effetctive
	sqlite3_exec (dB, "PRAGMA temp_store=MEMORY",NULL, NULL, &errmsg);
	sqlite3_exec (dB, "PRAGMA synchronous=OFF",NULL, NULL, &errmsg);
	sqlite3_exec (dB, "PRAGMA count_changes=OFF",NULL, NULL, &errmsg);
	sqlite3_exec (dB, "PRAGMA page_size=4092",NULL, NULL, &errmsg);
}

void metaDataContainer::close()
{
	if (dB)
	{
		sqlite3_close(dB);
		dB = 0;
	}
}

void metaDataContainer::init(const FileName baseFileName,int flag)
{
	dB = 0;
	mnBusyTimeoutMs = 60000; // 60 second
	open(baseFileName,flag);

	//Valid types and kind of variable
    //Angles
	validAttributes.insert(std::make_pair("angleRot" ,metaDataReal));
	validAttributes.insert(std::make_pair("angleTilt",metaDataReal));
	validAttributes.insert(std::make_pair("anglePsi" ,metaDataReal));
	//shifts
	validAttributes.insert(std::make_pair("shiftX"   ,metaDataReal));
	validAttributes.insert(std::make_pair("shiftY"   ,metaDataReal));
	validAttributes.insert(std::make_pair("shiftZ"   ,metaDataReal));
	//origin
	validAttributes.insert(std::make_pair("originX"  ,metaDataReal));
	validAttributes.insert(std::make_pair("originY"  ,metaDataReal));
	validAttributes.insert(std::make_pair("originZ"  ,metaDataReal));
	//filename
	validAttributes.insert(std::make_pair("fileName" ,metaDataString));
	//label (unique key) can work as primary key but it is not compulsory
	validAttributes.insert(std::make_pair("label",    metaDataStringUnique));
	//primary key, this is handled by the class not by the program
	//validAttributes.insert(std::make_pair("entryId",  metaDataPK));
	//belong to class
	validAttributes.insert(std::make_pair("belongToClass",metaDataInteger));
	//format
	validAttributes.insert(std::make_pair("formatId",metaDataInteger));
}

//import data from file, only attributeNames fields
metaDataContainer::metaDataContainer(FileName fileNameBase,
		                             std::string tableName,
		                             std::vector<std::string> &attributeNames,
		                             int flag)
{
	//part common for all creators
	init(fileNameBase,flag);
	//create table
	createTable(tableName);
    std::string invalidAttribute;
	//if read-only check if all asked attributes exist
	//otherwise create them if they di nit exists
	if(flag|SQLITE_OPEN_READONLY)
		if(checkAttributes(tableName,attributeNames,invalidAttribute)!=0)
			REPORT_ERROR(1,"Attribute \"" +invalidAttribute + "\" does not "
					+ "exist and data base open for reading only");
#ifdef NEVER
	else
		addAttributes(tableName,attributeNames);
#endif
}
#ifdef NEVER

//open database from program with attributeNames
metaDataContainer::metaDataContainer(FileName fileNameBase,
		                             std::vector<std::string> &attributeNames,
		                             int flag)
{
	//part common for all creators
	init(const FileName baseFileName,,flag);
	//TO BE DONE
	//create table (only one per file)
	//readCVSFile()
	//createTable()
	//addAttribute()
}
#endif

void metaDataContainer::createTable(const std::string tableName)
{

    /** 1. Create table if it does not exist **/
	std::string sqlCommand;
	sqlCommand = (std::string)"CREATE TABLE IF NOT EXISTS " +\
			     tableName +\
			     (std::string) "( entryId integer PRIMARY KEY AUTOINCREMENT) "
				 ;
	// Execute the query for   table creation
	if(sqlite3_exec(dB,sqlCommand.c_str(),0,0,0))
	{
    	REPORT_ERROR(1,sqlite3_errmsg(dB));
	}
}


int metaDataContainer::checkAttributes(const std::string tableName,
		                              std::vector<std::string> &attributeNames,
		                              std::string &invalidAttribute
		                              )
{
	//table should exist
	if(!tableExists(tableName))
    	REPORT_ERROR(1,sqlite3_errmsg(dB));

	std::vector<std::string>::iterator attributeNamesIterator;
    const char *datatype;
    const char *colseq;
    int    NotNull,PrimaryKey,Autoinc;
    int    numberColDoNOTexist=0;

    for(attributeNamesIterator =attributeNames.begin();
		attributeNamesIterator != attributeNames.end();
		attributeNamesIterator++)
	{
    	if(sqlite3_table_column_metadata(dB,
    			                        NULL,
    			                        tableName.c_str(),
    			                        (*attributeNamesIterator).c_str(),
    			                        &datatype,
    			                        &colseq,
    			                        &NotNull,
    			                        &PrimaryKey,
    			                        &Autoinc)!=SQLITE_OK)
    	{
    		numberColDoNOTexist++;
    		invalidAttribute = (*attributeNamesIterator);
    		break;
    	}
	}

//    int sqlite3_table_column_metadata(
//      sqlite3 *db,                /* Connection handle */
//      const char *zDbName,        /* Database name or NULL */
//      const char *zTableName,     /* Table name */
//      const char *zColumnName,    /* Column name */
//      char const **pzDataType,    /* OUTPUT: Declared data type */
//      char const **pzCollSeq,     /* OUTPUT: Collation sequence name */
//      int *pNotNull,              /* OUTPUT: True if NOT NULL constraint exists */
//      int *pPrimaryKey,           /* OUTPUT: True if column part of PK */
//      int *pAutoinc               /* OUTPUT: True if column is auto-increment */

    return (numberColDoNOTexist != 0);
}

#ifdef NEVER

void metaDataContainer::addAttributes(const std::string tableName,
		                              std::vector<std::string> &attributeNames
		                              )
{
	//table should exist
	if(!tableExists(tableName))
    	REPORT_ERROR(1,sqlite3_errmsg(dB));
    // get number of columns
	int nCols = sqlite3_column_count(statement);
	std::vector<std::string>::iterator attributeNamesIterator;
	std::string sqlCommand;
	//add all column, if there are already there they will fail
    for(attributeNamesIterator =attributeNames.begin();
		attributeNamesIterator != attributeNames.end();
		attributeNamesIterator++)
	{
    	sqlCommand = (std::string)"ALTER TABLE " +\
    			     tableName +\
    			     " ADD COLUMN " +\
    			     *attributeNamesIterator + \
    			     (std::string)" " +validAttributes[*attributeNamesIterator]\
    			     ;
    	if(sqlite3_exec(dB,sqlCommand,0,0,0))
    	{
    		std::string error = sqlite3_errmsg(dB);
            if(errot.find((std::string)"duplicate"))
            	std::cerr << "column "
            	          << *attributeNamesIterator
            	          << " already exist "
            	          << std::endl;
            else
            	REPORT_ERROR(1,sqlite3_errmsg(dB));
    	}
	}
}
#endif

metaDataContainer::~metaDataContainer()
{
	close();
}

void metaDataContainer::setBusyTimeout(int nMillisecs)
{
	mnBusyTimeoutMs = nMillisecs;
	sqlite3_busy_timeout(dB, mnBusyTimeoutMs);
}

bool metaDataContainer::tableExists(const std::string tableName)
{
	std::string sqlCommand;
	sqlCommand = (std::string)"select count(*) from sqlite_master where type='table' and name=?";
    if(sqlite3_prepare_v2(dB, sqlCommand.c_str(), -1, &statement, 0) != SQLITE_OK)
    	REPORT_ERROR(1,sqlite3_errmsg(dB));
	if(sqlite3_bind_text(statement, 1, tableName.c_str(),-1, SQLITE_TRANSIENT)!= SQLITE_OK)
    	REPORT_ERROR(1,sqlite3_errmsg(dB));
    int result = sqlite3_step(statement);
    if(result != SQLITE_ROW)
    	REPORT_ERROR(1,sqlite3_errmsg(dB));
    int numberAnswers = sqlite3_column_int(statement, 0);
	//sqlite3_reset(statement);
	sqlite3_finalize(statement);
    return (numberAnswers > 0);
}



//make prepared query and iterate trough results
//read cvs, skip #, process first colum in a different awya
//write cvs

//if statement open close it
//start_reading_data
//bind
//next
//end_reading_data

//if statement open close it
//start_writing_data
//bind
//next
//end_writing_data
//begin_transaction
//insert row


//executeSQLcommand
#ifdef NEVEREVER
bool metaDataContainer::tableTypeCompatible(const std::string tableName,
		                  std::vector<std::string> &attributeNames,
		                  std::vector<int> &attributeTypes)
{
    //if new, table is compatible
	if(!tableExists(tableName))
		return (1);
	std::string sqlCommand;
	sqlCommand = (std::string)"select * from " + tableName;
	if(sqlite3_prepare_v2(dB, sqlCommand.c_str(), -1, &statement, 0) != SQLITE_OK)
		REPORT_ERROR(1,sqlite3_errmsg(dB));
    int result = sqlite3_step(statement);
	int nCols = sqlite3_column_count(statement);
	std::vector<int>::iterator attributeTypesIterator;
	std::vector<std::string>::iterator attributeNamesIterator;
	bool returnValue=0;
	for(attributeTypesIterator = attributeTypes.begin(),
		attributeNamesIterator =attributeNames.begin();
		attributeTypesIterator != attributeTypes.end();
		attributeTypesIterator++,attributeNamesIterator++)
	{
		for (int nField = 0; nField < nCols; nField++)
		{
			const char* szTemp = sqlite3_column_name(statement, nField);
			int iType = sqlite3_column_type(statement, nField);
			std::string Temp=szTemp;
			if((*attributeNamesIterator).size() != Temp.size())
			 {
				continue;
			 }
			if((*attributeNamesIterator).compare(Temp)==0)
			{
				if(iType==(*attributeTypesIterator))
			    {
					returnValue=1;
					std::cout << " column name " << szTemp << "with type " << iType << "is a " << (*attributeTypesIterator) <<std::endl;
					break;
			    }
				else
					returnValue=0;
			}
		}
		if(returnValue==0)
			break;
	}
	sqlite3_finalize(statement);
    return (returnValue > 0);
}//FALTA CASO NO HAY DATOS y devuelve NULL, entonces borrar tabla y recrearla

int CppSQLite3Query::getIntField(int nField, int nNullValue/*=0*/)
{
	if (fieldDataType(nField) == SQLITE_NULL)
	{
		return nNullValue;
	}
	else
	{
		return sqlite3_column_int(mpVM, nField);
	}
}


int CppSQLite3Query::getIntField(const char* szField, int nNullValue/*=0*/)
{
	int nField = fieldIndex(szField);
	return getIntField(nField, nNullValue);
}


double CppSQLite3Query::getFloatField(int nField, double fNullValue/*=0.0*/)
{
	if (fieldDataType(nField) == SQLITE_NULL)
	{
		return fNullValue;
	}
	else
	{
		return sqlite3_column_double(mpVM, nField);
	}
}


double CppSQLite3Query::getFloatField(const char* szField, double fNullValue/*=0.0*/)
{
	int nField = fieldIndex(szField);
	return getFloatField(nField, fNullValue);
}


const char* CppSQLite3Query::getStringField(int nField, const char* szNullValue/*=""*/)
{
	if (fieldDataType(nField) == SQLITE_NULL)
	{
		return szNullValue;
	}
	else
	{
		return (const char*)sqlite3_column_text(mpVM, nField);
	}
}


const char* CppSQLite3Query::getStringField(const char* szField, const char* szNullValue/*=""*/)
{
	int nField = fieldIndex(szField);
	return getStringField(nField, szNullValue);
}

int CppSQLite3Query::fieldIndex(const char* szField)
{
	checkVM();

	if (szField)
	{
		for (int nField = 0; nField < mnCols; nField++)
		{
			const char* szTemp = sqlite3_column_name(mpVM, nField);

			if (strcmp(szField, szTemp) == 0)
			{
				return nField;
			}
		}
	}

	throw CppSQLite3Exception(CPPSQLITE_ERROR,
							"Invalid field name requested",
							DONT_DELETE_MSG);
}

const char* CppSQLite3Query::fieldName(int nCol)
{
	checkVM();

	if (nCol < 0 || nCol > mnCols-1)
	{
		throw CppSQLite3Exception(CPPSQLITE_ERROR,
								"Invalid field index requested",
								DONT_DELETE_MSG);
	}

	return sqlite3_column_name(mpVM, nCol);
}

bool CppSQLite3Query::eof()
void CppSQLite3Query::nextRow()
int CppSQLite3Table::numRows()
{
	checkResults();
	return mnRows;
}


const char* CppSQLite3Table::fieldValue(int nField)
{
	checkResults();

	if (nField < 0 || nField > mnCols-1)
	{
		throw CppSQLite3Exception(CPPSQLITE_ERROR,
								"Invalid field index requested",
								DONT_DELETE_MSG);
	}

	int nIndex = (mnCurrentRow*mnCols) + mnCols + nField;
	return mpaszResults[nIndex];
}



int metaDataContainer::addValue( std::string name, int value )
{
	void * newValue = (void *)(new int(value));
	return insertVoidPtr( name, newValue );
}

int metaDataContainer::addValue( std::string name, double value )
{
	void * newValue = (void *)(new double(value));
	return insertVoidPtr( name, newValue );
}

int metaDataContainer::insertVoidPtr( std::string name, void * value )
{
	// Return value for "insert" call
	std::pair<std::map<std::string, void *>::iterator,bool> ret;

	ret = values.insert( std::pair<std::string, void *>(name,value) );
	
	if (ret.second==false)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

void * metaDataContainer::getValue( std::string name )
{	
	std::map<std::string, void *>::iterator element; 

	element = values.find( name );

	if ( element == values.end( ) )
	{
		return NULL;
	}
	else
	{
		return element->second;
	}
}

bool metaDataContainer::valueExists( std::string name )
{
	if( values.find( name ) == values.end( ) )
	{
		return false;
	}
	else
	{
		return true;
	}
}

void metaDataContainer::deleteValue( std::string name )
{
	values.erase( name );
}
#endif
// Testing	
/*int main( )
{
	metaDataContainer * params = new metaDataContainer( );

	params->addValue( std::string("rot"), 5.);
	params->addValue( std::string("tilt"), 2.3);

	double rot = rot( params );
	double tilt = tilt( params );	

	std::cout << "ROT: " << rot << " TILT: " << tilt << std::endl;

	return 0;
}*/
