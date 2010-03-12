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

#include <map>
#include <string>
#include <iostream>
#include <external/sqlite-3.6.23/sqlite3.h>

class metaDataContainer
{
    private:
		sqlite3* mpDB;      //database pointer
		sqlite3_stmt* mpVM; //prepared statment
	    int mnBusyTimeoutMs; //lock time
    public:
		/** Constructor open data base

	    valid flags

	    #define SQLITE_OPEN_READONLY         0x00000001   Ok
	    #define SQLITE_OPEN_READWRITE        0x00000002   Ok
	    #define SQLITE_OPEN_CREATE           0x00000004   Ok
	    #define SQLITE_OPEN_NOMUTEX          0x00008000   Ok
	    #define SQLITE_OPEN_FULLMUTEX        0x00010000   Ok
	    #define SQLITE_OPEN_SHAREDCACHE      0x00020000   Ok
	    #define SQLITE_OPEN_PRIVATECACHE     0x00040000   Ok
	    */
	    metaDataContainer::metaDataContainer(FileName fileName,
				          int flag=SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE)
		{
			open(fileName,flag);
		}
	    void metaDataContainer::open(const FileName fileName,int flag);

#ifdef NEVER
	/** Container for pairs "name" and value. Note that void * allows to use
	    mixed types */
	std::map<std::string, void *> values;
	
	int insertVoidPtr( std::string name, void * value );

	public:
	
	/** Constructor open data base*/
	metaDataContainer(FileName fileName);
	
	/** Destructor */
	~metaDataContainer();
	
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
