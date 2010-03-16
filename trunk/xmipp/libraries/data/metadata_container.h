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

class metaDataContainer
{
	/** Container for pairs "name" and value. Note that void * allows to use
	    mixed types */
	std::map<std::string, void *> values;
	
	void insertVoidPtr( std::string name, void * value );

	public:
	
	/** Constructor */
	metaDataContainer();
	
	/** Destructor */
	~metaDataContainer();
	
	/** Create a new pair name-value of integer type */
	void addValue( std::string name, double value );
	void addValue( std::string name, float value );
	void addValue( std::string name, int value );
	void addValue( std::string name, bool value );
	void addValue( std::string name, std::string value );
	
	void * getValue( std::string name );
	bool valueExists( std::string name );
	
	bool pairExists( std::string name, double value );
	bool pairExists( std::string name, float value );
	bool pairExists( std::string name, int value );
	bool pairExists( std::string name, bool value );
	bool pairExists( std::string name, std::string value );
	
	void deleteValue( std::string name );
};

#endif
