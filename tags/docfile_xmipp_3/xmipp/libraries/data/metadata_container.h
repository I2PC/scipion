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
#include "strings.h"
#include <iostream>
#include <sstream>
#include <fstream>

// This enum defines what MetaDataLabels this class can manage, if
// you need a new one add it here and modify affected methods:
//
// 	- static MetaDataLabel codifyLabel( std::string strLabel );
//	- static std::string decodeLabel( MetaDataLabel inputLabel );      
//	- void writeValuesToFile( std::ofstream &outfile, MetaDataLabel inputLabel );
//	- void addValue( std::string name, std::string value );
//
// The order these labels are defined in the "label" enum
// is that that will be used when writing a metadata file
enum MetaDataLabel 
{ 
	MDL_UNDEFINED = -1,
	MDL_ANGLEROT, 
	MDL_ANGLETILT, 
	MDL_ANGLEPSI,
	MDL_IMAGE,
	MDL_MICROGRAPH,
	MDL_CTFMODEL,
	MDL_SHIFTX, 
	MDL_SHIFTY, 
	MDL_SHIFTZ, 
	MDL_ENABLED, 
	MDL_ORIGINX, 
	MDL_ORIGINY, 
	MDL_ORIGINZ,
	MDL_WEIGHT,
	MDL_FLIP,
	MDL_REF,
	MDL_MAXCC,
	MDL_LAST_LABEL	// NOTE: Do keep this label always at the end
			 	// it is here for looping purposes  	
};


class MetaDataContainer
{
	/** Container for pairs "name" and value. Note that void * allows to use
	    mixed types */
	std::map<MetaDataLabel, void *> values;
	
	void insertVoidPtr( MetaDataLabel name, void * value );

	public:
	
	/** Constructor */
	MetaDataContainer();
	
	/** Destructor */
	~MetaDataContainer();
	
	/** Create a new pair name-value of integer type */
	void addValue( MetaDataLabel name, double value );
	void addValue( MetaDataLabel name, float value );
	void addValue( MetaDataLabel name, int value );
	void addValue( MetaDataLabel name, bool value );
	void addValue( MetaDataLabel name, std::string value );
	void addValue( std::string name, std::string value );
	
	void * getValue( MetaDataLabel name );
	bool valueExists( MetaDataLabel name );
	
	bool pairExists( MetaDataLabel name, double value );
	bool pairExists( MetaDataLabel name, float value );
	bool pairExists( MetaDataLabel name, int value );
	bool pairExists( MetaDataLabel name, bool value );
	bool pairExists( MetaDataLabel name, std::string value );
	
	void deleteValue( MetaDataLabel name );
	
	void writeValuesToFile( std::ofstream &outfile, MetaDataLabel inputLabel );

	static MetaDataLabel codifyLabel( std::string strLabel );
	static std::string decodeLabel( MetaDataLabel inputLabel );
};

#endif
