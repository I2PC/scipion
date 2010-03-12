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
#include "metadata_container.h"
//called by creator
void metaDataContainer::open(const FileName fileName,int flag)
{
	int nRet = sqlite3_open_v2(fileName.c_str(), &mpDB,flag, NULL);

	if (nRet != SQLITE_OK)
	{
		REPORT_ERROR(1,sqlite3_errmsg(mpDB));
	}
	/*
	> Suppose you are trying to perform a select, but there's an update
	> currently in progress. Without sqlite3_busy_timeout or explicit busy
	> handler, your operation will fail immediately with SQLITE_BUSY. But if
	> you set up sqlite3_busy_timeout, SQLite will automatically retry your
	> operation several times before erroring out. Hopefully the writing
	> transaction completes before the timeout has expired, and your select is
	> allowed to proceed.
	 */
	setBusyTimeout(mnBusyTimeoutMs);
	char *errmsg;
	//some optimization, do not look to efective
	sqlite3_exec (mpDB, "PRAGMA temp_store=MEMORY",NULL, NULL, &errmsg);
	sqlite3_exec (mpDB, "PRAGMA synchronous=OFF",NULL, NULL, &errmsg);
	sqlite3_exec (mpDB, "PRAGMA count_changes=OFF",NULL, NULL, &errmsg);
	sqlite3_exec (mpDB, "PRAGMA page_size=4092",NULL, NULL, &errmsg);
}

metaDataContainer::metaDataContainer(FileName fileName)
{
	mpDB = 0;
	mnBusyTimeoutMs = 1000; // 1 second
}
metaDataContainer::~metaDataContainer(){};


#ifdef NEVER
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
