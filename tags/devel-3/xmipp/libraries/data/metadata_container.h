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

#ifndef METADATACONTAINER_H
#define METADATACONTAINER_H

#include <external/sqlite-3.6.23/sqlite3.h>
#include "funcs.h"
#include <iostream>
#include <map>
#include <string>

/** valid data types
 *
 */
#define metaDataInteger "INTEGER"
#define metaDataReal    "FLOAT"
#define metaDataString  "TEXT"
//primary key integer
#define metaDataPK                "PRIMARY KEY AUTOINCREMENT"
#define metaDataStringUnique      "TEXT UNIQUE"

/** Some other variables that may be of interest
 *
 */
#define metaDataBlob  SQLITE_BLOB
#define metaDataNull  SQLITE_NULL


/**
 * Valid image formats // check Sjors,...
 */

#define CCP4      1
#define MRC       2
#define PIF       3
#define EM        4
#define IMAGIC    5
#define SPIDER    6
#define SUPRIM    7
#define TIFF      8
#define PNG       9
#define JPEG     10

/**
 * Whenever make sense follow this naming convention
 * Images
 */

class metaDataContainer
{
    private:
		/** database handler
		 *
		 */
	    sqlite3* dB;
	    /** prepared statement for database
	     *
	     */
		sqlite3_stmt* statement;
		/**lock time for database (see function setBusyTimeout)
		 *
		 */
	    int mnBusyTimeoutMs; //lock time
	    /**
	     * valid attribute names are in map
	     */
	    std::map<std::string, std::string> validAttributes;

	    /** iterator for validAttributes
	    *
	    */
	    std::map<std::string, std::string>::iterator validAttributesIterator;

        /** open metaContainer DataBase for readind/writing ...  */
	    void open(const FileName fileName,int flag);

        /** close metaContainer DataBase */
	    void close();

	    /** Check if table exists */
	    bool tableExists(const std::string tableName);
		/** Block common for all creators
		 *
		 */
	    void init(const FileName baseFileName,  int flag);

	    /** check if attributes exits
	     *
	     */
	    int checkAttributes(const std::string tableName,
	    		              std::vector<std::string> &attributeNames,
                              std::string &invalidAttribute
	    		              );

    public:

	    /** Check if types are compatible */
	    bool tableTypeCompatible(const std::string tableName,
	    		                  std::vector<std::string> &attributeNames,
	    		                  std::vector<int> &attributeTypes);

    public:
		/** Constructor open data base
            vector has column names that should be retrieved
            from table tablename
	    valid flags

	    #define SQLITE_OPEN_READONLY         0x00000001
	    #define SQLITE_OPEN_READWRITE        0x00000002
	    #define SQLITE_OPEN_CREATE           0x00000004
	    #define SQLITE_OPEN_NOMUTEX          0x00008000
	    #define SQLITE_OPEN_FULLMUTEX        0x00010000
	    #define SQLITE_OPEN_SHAREDCACHE      0x00020000
	    #define SQLITE_OPEN_PRIVATECACHE     0x00040000
	    */
        /** Read metaData from file
         *  check open for flags (read(write by default)
         */
	    metaDataContainer(FileName fileNameBase,\
	    		          std::string tableName,\
	    		          std::vector<std::string> &attributeNames,\
	    		          int flag=(SQLITE_OPEN_READONLY|SQLITE_OPEN_SHAREDCACHE));

	    /**
	     * check open for flags (read only by default)
	     */
	    metaDataContainer(FileName fileNameBase,\
	    		          std::vector<std::string> &attributeNames,\
	    		          int flag=SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);
	    /** Destructor */
		~metaDataContainer();

		/** 	set maximum waiting time
		Suppose you are trying to perform a select, but there's an update
		currently in progress. Without sqlite3_busy_timeout or explicit busy
		handler, your operation will fail immediately with SQLITE_BUSY. But if
		you set up sqlite3_busy_timeout, SQLite will automatically retry your
		operation several times before erroring out. Hopefully the writing
		transaction completes before the timeout has expired, and your select is
		allowed to proceed.
	    */
		void setBusyTimeout(int nMillisecs);
		/** Create table
		 * create (if needed) table with primary key
		 */
		void createTable(const std::string tableName);

#ifdef NEVER
	/** Container for pairs "name" and value. Note that void * allows to use
	    mixed types */
	std::map<std::string, void *> values;
	
	int insertVoidPtr( std::string name, void * value );

	public:
	
	/** Constructor open data base*/
	metaDataContainer(FileName fileName);
	

	
	/** Create a new pair name-value of integer type */
	int addValue( std::string name, int value );
	
	/** Create a new pair name-value of double type */
	int addValue( std::string name, double value );
	
	void * getValue( std::string name );
	bool valueExists( std::string name );
	void deleteValue( std::string name );
#endif
};

#endif
