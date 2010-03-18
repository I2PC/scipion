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
#include <sstream>

enum label 
{ 
	ANGLEROT, 
	ANGLETILT, 
	ANGLEPSI,
	IMAGE,
	MICROGRAPH,
	CTFMODEL,
	SHIFTX, 
	SHIFTY, 
	SHIFTZ, 
	ENABLED, 
	ORIGINX, 
	ORIGINY, 
	ORIGINZ,
	WEIGHT,
	FLIP
};

class metaDataContainer
{
	/** Container for pairs "name" and value. Note that void * allows to use
	    mixed types */
	std::map<label, void *> values;
	
	void insertVoidPtr( label name, void * value );

	public:
	
	/** Constructor */
	metaDataContainer();
	
	/** Destructor */
	~metaDataContainer();
	
	/** Create a new pair name-value of integer type */
	void addValue( label name, double value );
	void addValue( label name, float value );
	void addValue( label name, int value );
	void addValue( label name, bool value );
	void addValue( label name, std::string value );
	void addValue( std::string name, std::string value );
	
	void * getValue( label name );
	bool valueExists( label name );
	
	bool pairExists( label name, double value );
	bool pairExists( label name, float value );
	bool pairExists( label name, int value );
	bool pairExists( label name, bool value );
	bool pairExists( label name, std::string value );
	
	void deleteValue( label name );
	
	static label codifyLabel( std::string strLabel );
	static std::string decodeLabel( label inputLabel );
};

#endif
