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

MetaDataContainer::MetaDataContainer(){};
MetaDataContainer::~MetaDataContainer(){};

void MetaDataContainer::addValue( MetaDataLabel name, int value )
{
	void * newValue = (void *)(new int(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( MetaDataLabel name, double value )
{
    void * newValue = (void *)(new double(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( MetaDataLabel name, float value )
{
	void * newValue = (void *)(new float(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( MetaDataLabel name, bool value )
{
	void * newValue = (void *)(new bool(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( MetaDataLabel name, std::string value )
{
	void * newValue = (void *)(new std::string(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( std::string name, std::string value )
{	
	MetaDataLabel lCode = codifyLabel( name );
	std::istringstream i( value );
	
	// Look for a double value
	if( lCode == MDL_ANGLEROT || lCode == MDL_ANGLETILT || lCode == MDL_ANGLEPSI ||
	   lCode == MDL_SHIFTX || lCode == MDL_SHIFTY || lCode == MDL_SHIFTZ ||
	   lCode == MDL_ORIGINX || lCode == MDL_ORIGINY || lCode == MDL_ORIGINZ || 
	   lCode == MDL_WEIGHT || lCode == MDL_MAXCC ) 
	{
		double doubleValue;
		
		i >> doubleValue;			
				
		addValue( lCode, doubleValue );
	}
	else if( lCode == MDL_IMAGE || lCode == MDL_MICROGRAPH || lCode == MDL_CTFMODEL )
	{
		addValue( lCode, value );
	}
	else if( lCode == MDL_REF || lCode == MDL_ENABLED )
	{
		int intValue;
		
		i >> intValue;
				
		addValue( lCode, intValue ); 
	}
	else if( lCode == MDL_FLIP )
	{
		bool boolValue;
		
		i >> boolValue;
				
		addValue( lCode, boolValue );
	}
	
}

void MetaDataContainer::insertVoidPtr( MetaDataLabel name, void * value )
{
	values[ name ] = value;
}

void * MetaDataContainer::getValue( MetaDataLabel name )
{	
	std::map<MetaDataLabel, void *>::iterator element; 

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

bool MetaDataContainer::valueExists( MetaDataLabel name )
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

bool MetaDataContainer::pairExists( MetaDataLabel name, double value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< MetaDataLabel, void *>::iterator It;
	
	It = values.find( name );
	
	if( It != values.end( ))
	{
		if( value == *((double *)(It->second)) )
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}

bool MetaDataContainer::pairExists( MetaDataLabel name, float value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< MetaDataLabel, void *>::iterator It;
	
	It = values.find( name );
	
	if( It != values.end( ))
	{
		if( value == *((float *)(It->second)) )
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}

bool MetaDataContainer::pairExists( MetaDataLabel name, int value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< MetaDataLabel, void *>::iterator It;
	
	It = values.find( name );
	
	if( It != values.end( ))
	{
		if( value == *((int *)(It->second)) )
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}

bool MetaDataContainer::pairExists( MetaDataLabel name, bool value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< MetaDataLabel, void *>::iterator It;
	
	It = values.find( name );
	
	if( It != values.end( ))
	{
		if( value == *((bool *)(It->second)) )
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}

bool MetaDataContainer::pairExists( MetaDataLabel name, std::string value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< MetaDataLabel, void *>::iterator It;
	
	It = values.find( name );
	
	if( It != values.end( ))
	{
		if( value == *((std::string *)It->second) )
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}

void MetaDataContainer::deleteValue( MetaDataLabel name )
{
	values.erase( name );
}
	
MetaDataLabel MetaDataContainer::codifyLabel( std::string strLabel )
{
	// NOTE: For the "||" cases, the left MetaDataLabel is the XMIPP_3 version's and
	// the right value is the old-styled one

	if( strLabel == "angleRot" || strLabel == "rot" )
	{
		return MDL_ANGLEROT;
	}
	else if( strLabel == "angleTilt" || strLabel == "tilt" )
	{
		return MDL_ANGLETILT;
	}
	else if( strLabel == "anglePsi" || strLabel == "psi" )
	{
		return MDL_ANGLEPSI;
	}
	else if( strLabel == "image" )
	{
		return MDL_IMAGE;
	}
	else if( strLabel == "micrograph" )
	{
		return MDL_MICROGRAPH;
	}
	else if( strLabel == "CTFModel" )
	{
		return MDL_CTFMODEL;
	}
	else if( strLabel == "shiftX" || strLabel == "Xoff" )
	{
		return MDL_SHIFTX;
	}
	else if( strLabel == "shiftY" || strLabel == "Yoff" )
	{
		return MDL_SHIFTY;
	}
	else if( strLabel == "shiftZ" || strLabel == "Zoff" )
	{
		return MDL_SHIFTZ;
	}
	else if( strLabel == "enabled" )
	{
		return MDL_ENABLED;
	}
	else if( strLabel == "originX" )
	{
		return MDL_ORIGINX;
	}
	else if( strLabel == "originY" )
	{
		return MDL_ORIGINY;
	}
	else if( strLabel == "originZ" )
	{
		return MDL_ORIGINZ;
	}
	else if( strLabel == "weight" || strLabel == "Weight"  )
	{
		return MDL_WEIGHT;
	}
	else if( strLabel == "flip" || strLabel == "Flip" )
	{
		return MDL_FLIP;
	}
	else if( strLabel == "ref" || strLabel == "Ref" )
	{
		return MDL_REF;
	}
	else if( strLabel == "maxCC" )
	{
		return MDL_MAXCC;
	}
	else
	{
		return MDL_UNDEFINED;
	}
}

std::string MetaDataContainer::decodeLabel( MetaDataLabel inputLabel )
{
	switch ( inputLabel ) 
	{
		case MDL_ANGLEROT:
			return std::string( "angleRot" );
			break;
		case MDL_ANGLETILT:
			return std::string( "angleTilt" );
			break;
		case MDL_ANGLEPSI:
			return std::string( "anglePsi" );
			break;
		case MDL_IMAGE:
			return std::string( "image" );
			break;
		case MDL_MICROGRAPH:
			return std::string( "micrograph" );
			break;
		case MDL_CTFMODEL:
			return std::string( "CTFModel" );
			break;
		case MDL_SHIFTX:
			return std::string( "shiftX" );
			break;
		case MDL_SHIFTY:
			return std::string( "shiftY" );
			break;
		case MDL_SHIFTZ:
			return std::string( "shiftZ" );
			break;
		case MDL_ENABLED:
			return std::string( "enabled" );
			break;
		case MDL_ORIGINX:
			return std::string( "originX" );
			break;
		case MDL_ORIGINY:
			return std::string( "originY" );
			break;
		case MDL_ORIGINZ:
			return std::string( "originZ" );
			break;
		case MDL_WEIGHT:
			return std::string( "weight" );
			break;
		case MDL_FLIP: 
			return std::string( "flip" );
			break;
		case MDL_REF: 
			return std::string( "ref" );
			break;
		case MDL_MAXCC: 
			return std::string( "maxCC" );
			break;
		default:
			return std::string( "" );
			break;
	}
}

void MetaDataContainer::writeValueToFile( std::ofstream &outfile, MetaDataLabel inputLabel )
{
	switch ( inputLabel ) 
	{
		case MDL_ANGLEROT:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_ANGLETILT:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_ANGLEPSI:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_IMAGE:
			outfile << *((std::string*)(getValue( inputLabel )));
			break;
		case MDL_MICROGRAPH:
			outfile << *((std::string*)(getValue( inputLabel )));
			break;
		case MDL_CTFMODEL:
			outfile << *((std::string*)(getValue( inputLabel )));
			break;
		case MDL_SHIFTX:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_SHIFTY:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_SHIFTZ:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_ENABLED:
			outfile << *((int*)(getValue( inputLabel )));
			break;
		case MDL_ORIGINX:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_ORIGINY:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_ORIGINZ:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_WEIGHT:
			outfile << *((double*)(getValue( inputLabel )));
			break;
		case MDL_FLIP: 
			outfile << *((bool*)(getValue( inputLabel )));
			break;
		case MDL_REF: 
			outfile << *((int*)(getValue( inputLabel )));
			break;
		case MDL_MAXCC: 
			outfile << *((double*)(getValue( inputLabel )));
			break;
		default:
			break;
	}
}

bool MetaDataContainer::isValidLabel( MetaDataLabel inputLabel )
{
	if( inputLabel > MDL_UNDEFINED && inputLabel < MDL_LAST_LABEL )
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool MetaDataContainer::isValidLabel( std::string inputLabel )
{
	MetaDataLabel label = codifyLabel( inputLabel );
	
	if( label > MDL_UNDEFINED && label < MDL_LAST_LABEL )
	{
		return true;
	}
	else
	{
		return false;
	}
}
