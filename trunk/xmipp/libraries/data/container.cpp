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
#include "container.h"

xmpContainer::xmpContainer(){};
xmpContainer::~xmpContainer(){};

int xmpContainer::addValue( std::string name, int value )
{
	void * newValue = (void *)(new int(value));
	return insertVoidPtr( name, newValue );
}

int xmpContainer::addValue( std::string name, double value )
{
	void * newValue = (void *)(new double(value));
	return insertVoidPtr( name, newValue );
}

int xmpContainer::insertVoidPtr( std::string name, void * value )
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

void * xmpContainer::getValue( std::string name )
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

bool xmpContainer::valueExists( std::string name )
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

void xmpContainer::deleteValue( std::string name )
{
	values.erase( name );
}
	
// Testing	
/*int main( )
{
	xmpContainer * params = new xmpContainer( );

	params->addValue( std::string("rot"), 5.);
	params->addValue( std::string("tilt"), 2.3);

	double rot = rot( params );
	double tilt = tilt( params );	

	std::cout << "ROT: " << rot << " TILT: " << tilt << std::endl;

	return 0;
}*/
