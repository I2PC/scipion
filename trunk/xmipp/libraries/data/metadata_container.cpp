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

void MetaDataContainer::addValue( label name, int value )
{
	void * newValue = (void *)(new int(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( label name, double value )
{
	void * newValue = (void *)(new double(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( label name, float value )
{
	void * newValue = (void *)(new float(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( label name, bool value )
{
	void * newValue = (void *)(new bool(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( label name, XmpString value )
{
	void * newValue = (void *)(new XmpString(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( XmpString name, XmpString value )
{	
	label lCode = codifyLabel( name );
	std::istringstream i( value );
	
	// Look for a double value
	if( lCode == ANGLEROT || lCode == ANGLETILT || lCode == ANGLEPSI ||
	   lCode == SHIFTX || lCode == SHIFTY || lCode == SHIFTZ ||
	   lCode == ORIGINX || lCode == ORIGINY || lCode == ORIGINZ ) 
	{
		double doubleValue;
		
		i >> doubleValue;			
		
		addValue( lCode, doubleValue );
	}
	else if( lCode == IMAGE || lCode == MICROGRAPH || lCode == CTFMODEL )
	{
		addValue( lCode, value );
	}
	else if( lCode == ENABLED )
	{
		bool boolValue;
		
		i >> boolValue;
		
		addValue( lCode, boolValue );
	}
	
}

void MetaDataContainer::insertVoidPtr( label name, void * value )
{
	values[ name ] = value;
}

void * MetaDataContainer::getValue( label name )
{	
	std::map<label, void *>::iterator element; 

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

bool MetaDataContainer::valueExists( label name )
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

bool MetaDataContainer::pairExists( label name, double value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< label, void *>::iterator It;
	
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

bool MetaDataContainer::pairExists( label name, float value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< label, void *>::iterator It;
	
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

bool MetaDataContainer::pairExists( label name, int value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< label, void *>::iterator It;
	
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

bool MetaDataContainer::pairExists( label name, bool value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< label, void *>::iterator It;
	
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

bool MetaDataContainer::pairExists( label name, XmpString value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< label, void *>::iterator It;
	
	It = values.find( name );
	
	if( It != values.end( ))
	{
		if( value == *((XmpString *)(It->second)) )
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}

void MetaDataContainer::deleteValue( label name )
{
	values.erase( name );
}
	
label MetaDataContainer::codifyLabel( XmpString strLabel )
{
	if( strLabel == "angleRot" || strLabel == "rot" )
	{
		return ANGLEROT;
	}
	else if( strLabel == "angleTilt" || strLabel == "tilt" )
	{
		return ANGLETILT;
	}
	else if( strLabel == "anglePsi" || strLabel == "psi" )
	{
		return ANGLEPSI;
	}
	else if( strLabel == "image" )
	{
		return IMAGE;
	}
	else if( strLabel == "micrograph" )
	{
		return MICROGRAPH;
	}
	else if( strLabel == "CTFModel" )
	{
		return CTFMODEL;
	}
	else if( strLabel == "shiftX" || strLabel == "Xoff" )
	{
		return SHIFTX;
	}
	else if( strLabel == "shiftY" || strLabel == "Yoff" )
	{
		return SHIFTY;
	}
	else if( strLabel == "shiftZ" || strLabel == "Zoff" )
	{
		return SHIFTZ;
	}
	else if( strLabel == "enabled" )
	{
		return ENABLED;
	}
	else if( strLabel == "originX" )
	{
		return ORIGINX;
	}
	else if( strLabel == "originY" )
	{
		return ORIGINY;
	}
	else if( strLabel == "originZ" )
	{
		return ORIGINZ;
	}
	else if( strLabel == "Weight" )
	{
		return WEIGHT;
	}
	else if( strLabel == "Flip" )
	{
		return FLIP;
	}
}

XmpString MetaDataContainer::decodeLabel( label inputLabel )
{
	switch ( inputLabel ) 
	{
		case ANGLEROT:
			return XmpString( "angleRot" );
			break;
		case ANGLETILT:
			return XmpString( "angleTilt" );
			break;
		case ANGLEPSI:
			return XmpString( "anglePsi" );
			break;
		case IMAGE:
			return XmpString( "image" );
			break;
		case MICROGRAPH:
			return XmpString( "micrograph" );
			break;
		case CTFMODEL:
			return XmpString( "CTFModel" );
			break;
		case SHIFTX:
			return XmpString( "shiftX" );
			break;
		case SHIFTY:
			return XmpString( "shiftY" );
			break;
		case SHIFTZ:
			return XmpString( "shiftZ" );
			break;
		case ENABLED:
			return XmpString( "enabled" );
			break;
		case ORIGINX:
			return XmpString( "originX" );
			break;
		case ORIGINY:
			return XmpString( "originY" );
			break;
		case ORIGINZ:
			return XmpString( "originZ" );
			break;
		case WEIGHT:
			return XmpString( "weight" );
			break;
		case FLIP: 
			return XmpString( "flip" );
			break;
		default:
			return XmpString( "" );
			break;
	}
}
