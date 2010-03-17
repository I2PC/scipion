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

#include "metadata.h"

metaData::metaData()
{
	programName = std::string("");
	path = std::string("");
	objectsType = UNDEFINED_TYPE;
	objects.clear( );	
}

metaData::metaData( std::string fileName, std::string newObjectsType, std::vector<std::string> * labelsVector )
{
	programName = std::string("");
	path = std::string("");
	objectsType = newObjectsType;
}

metaData::~metaData()
{
	objects.clear( );
}

bool metaData::isEmpty( )
{
	return objects.empty( );
}

void metaData::clear( )
{
	programName = std::string("");
	path = std::string("");
	objectsType = UNDEFINED_TYPE;
	objects.clear( );		
}

void metaData::setProgram( std::string newProgName )
{
	programName = newProgName;
}

void metaData::setPath( std::string newPath )
{
	path = newPath;
}

void metaData::setType( unsigned int newObjectsType )
{
	objectsType = newObjectsType;
}

std::string metaData::getProgram( )
{
	return programName;
}

std::string metaData::getPath( )
{
	return path;
}

unsigned int metaData::getType( )
{
	return objectsType;
}

long int metaData::addObject( )
{
	long int result = lastObject();
	typedef std::pair<long int, metaDataContainer *> newPair;
	
	if( result == NO_OBJECTS_STORED )
	{
		result = 0;
	}
	else
	{
		result++;
	}

	objects.insert( newPair( result, NULL ) );
	
	return result;
}

size_t metaData::firstObject( ) 
{ 
	size_t result = 0;
	
	if( !objects.empty( ))
	{
		objectsIterator = objects.begin(); 
		result = objectsIterator->first;
	}
	else
	{
		result = NO_OBJECTS_STORED;	// Map is empty
	}
	
	return result;
};

size_t metaData::nextObject( ) 
{ 
	size_t result = 0;
	
	if( !objects.empty( ))
	{
		objectsIterator++;
		
		if( objectsIterator != objects.end( ) )
		{
			result = objectsIterator->first;
		}
		else
		{
			result = NO_MORE_OBJECTS;
		}
	}
	else
	{
		result = NO_OBJECTS_STORED;
	}
	
	return result;
};

long int metaData::lastObject( )
{
	size_t result = 0;
	
	if( !objects.empty( ))
	{
		objectsIterator = objects.end(); 
		objectsIterator--;
		
		result = objectsIterator->first;
	}
	else
	{
		result = NO_OBJECTS_STORED;
	}
	
	return result;
};

bool metaData::setValue( std::string label, double value, long int objectID )
{
	long int auxID;
	
	if( !objects.empty( ))
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		metaDataContainer * aux = objects[auxID];
		
		aux->setValue( label, value );
		
		return true;
	}	
	else
	{
		return false;
	}
}

bool metaData::setValue( std::string label, float value, long int objectID )
{
	long int auxID;
	
	if( !objects.empty( ))
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		metaDataContainer * aux = objects[auxID];
		
		aux->setValue( label, value );

		return true;
	}	
	else
	{
		return false;
	}	
}

bool metaData::setValue( std::string label, int value, long int objectID )
{
	long int auxID;
	
	if( !objects.empty( ))
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		metaDataContainer * aux = objects[auxID];
		
		aux->setValue( label, value );

		return true;
	}	
	else
	{
		return false;
	}	
}

bool metaData::setValue( std::string label, std::string value, long int objectID )
{
	long int auxID;
	
	if( !objects.empty( ))
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		metaDataContainer * aux = objects[auxID];
		
		aux->setValue( label, value );
		
		return true;
	}	
	else
	{
		return false;
	}	
}

std::vector<long int> metaData::findObjects( std::string label, double value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( label, value ) )
			result.push_back( It->first ):
	}
	
	return result;
}

std::vector<long int> metaData::findObjects( std::string label, float value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( label, value ) )
			result.push_back( It->first ):
	}
	
	return result;
}

std::vector<long int> metaData::findObjects( std::string label, int value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( label, value ) )
			result.push_back( It->first ):
	}
	
	return result;
}

std::vector<long int> metaData::findObjects( std::string label, bool value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( label, value ) )
			result.push_back( It->first ):
	}
	
	return result;
}

std::vector<long int> metaData::findObjects( std::string label, std::string value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( label, value ) )
			result.push_back( It->first ):
			}
	
	return result;
}

double metadata::rot( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("angleRot") );
		
		if( result == NULL )
		{
			std::cerr << "No 'angleRot' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metadata::tilt( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("angleTilt") );
		
		if( result == NULL )
		{
			std::cerr << "No 'angleTilt' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metadata::psi( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("anglePsi") );
		
		if( result == NULL )
		{
			std::cerr << "No 'anglePs' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

bool metadata::enabled( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		bool * result = aux->getValue( std::string("enabled") );
		
		if( result == NULL )
		{
			std::cerr << "No 'enabled' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

std::string metadata::fileName( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		std::string * result = aux->getValue( std::string("fileName") );
		
		if( result == NULL )
		{
			std::cerr << "No 'fileName' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metadata::shiftX( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("shiftX") );
		
		if( result == NULL )
		{
			std::cerr << "No 'shiftX' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metadata::shiftY( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("shiftY") );
		
		if( result == NULL )
		{
			std::cerr << "No 'shiftY' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metadata::shiftZ( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("shiftZ") );
		
		if( result == NULL )
		{
			std::cerr << "No 'shiftZ' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metadata::originX( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("originX") );
		
		if( result == NULL )
		{
			std::cerr << "No 'originX' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metadata::originY( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("originY") );
		
		if( result == NULL )
		{
			std::cerr << "No 'originY' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metadata::originZ( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = aux->getValue( std::string("originZ") );
		
		if( result == NULL )
		{
			std::cerr << "No 'originZ' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}