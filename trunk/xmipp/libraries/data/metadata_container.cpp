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

metaDataContainer::metaDataContainer(){};
metaDataContainer::~metaDataContainer(){};

void metaDataContainer::addValue( std::string name, int value )
{
	void * newValue = (void *)(new int(value));
	insertVoidPtr( name, newValue );
}

void metaDataContainer::addValue( std::string name, double value )
{
	void * newValue = (void *)(new double(value));
	insertVoidPtr( name, newValue );
}

void metaDataContainer::addValue( std::string name, float value )
{
	void * newValue = (void *)(new float(value));
	insertVoidPtr( name, newValue );
}

void metaDataContainer::addValue( std::string name, bool value )
{
	void * newValue = (void *)(new bool(value));
	insertVoidPtr( name, newValue );
}

void metaDataContainer::addValue( std::string name, std::string value )
{
	void * newValue = (void *)(new std::string(value));
	insertVoidPtr( name, newValue );
}

void metaDataContainer::insertVoidPtr( std::string name, void * value )
{
	values[ name ] = value;
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

bool metaDataContainer::pairExists( std::string name, double value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< std::string, void *>::iterator It;
	
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

bool metaDataContainer::pairExists( std::string name, float value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< std::string, void *>::iterator It;
	
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

bool metaDataContainer::pairExists( std::string name, int value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< std::string, void *>::iterator It;
	
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

bool metaDataContainer::pairExists( std::string name, bool value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< std::string, void *>::iterator It;
	
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

bool metaDataContainer::pairExists( std::string name, std::string value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< std::string, void *>::iterator It;
	
	It = values.find( name );
	
	if( It != values.end( ))
	{
		if( value == *((std::string *)(It->second)) )
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}

void metaDataContainer::deleteValue( std::string name )
{
	values.erase( name );
}
	
