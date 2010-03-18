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
	path = std::string("");
	objects.clear( );	
}

std::string metaData::eraseChar( std::string origStr, char character )
{
	std::string temp;
			
	for( unsigned int i = 0 ; i < origStr.length( ) ; i++ )
		if ( origStr[ i ] != character ) 
			temp += origStr[ i ];

	return temp;
}

metaData::metaData( std::string fileName, std::vector<label> * labelsVector )
{
	path = std::string("");

	// Open file
	std::ifstream infile ( fileName.data(), std::ios_base::in );
	std::string line;
	
	// Search for Headerinfo, if present we are processing an old-styled docfile
	// else we are processing a new Xmipp metadata file
	getline( infile, line, '\n');
	
	int pos = line.find( "Headerinfo" );
	
	if( pos != std::string::npos ) // Headerinfo token founds
	{
		// Remove from the beginning to the end of "Headerinfo columns:"
		pos = line.find( ":" );
		line = line.erase( 0, pos + 1 );
		
		std::vector<std::string> labelsVector;
		
		// Extract labels until the string is empty
		while ( line != "" )
		{
			pos = line.find( ")" );
			std::string label = line.substr( 0, pos+1 );
			line.erase( 0,pos+1);
			
			// The token can now contain a ',', if so, remove it
			if( ( pos = label.find( "," ) ) != std::string::npos )
				label.erase( pos, 1 );
			
			// Remove unneded parentheses and contents
			pos = label.find( "(" );
			int pos2 = label.find( ")" );
			label.erase( pos, pos2-pos+1 );
			
			// Remove white spaces
			label = eraseChar( label, ' ' );
		
			labelsVector.push_back( label );
		}
		
		int isname=0;
		while ( getline( infile, line, '\n') )
		{
			if( isname % 2 == 0 )
			{
				long int objectID = addObject( );
				line.erase(0,line.find(";")+1);
				
				// Remove spaces from string
				line = eraseChar( line, ' ' );
				
				setValue( IMAGE, line );
			}
			else
			{
				// Parse labels
				std::stringstream os2( line );          
				std::string value;
				
				int counter = 0;
				while ( os2 >> value )
				{
					if( counter >= 2 ) // Discard two first numbers
						setValue( labelsVector[counter-2], value );
					counter++;
				}
			}
		}
	}
	else
	{
		pos = line.rfind( "*" );
		
		if( pos == std::string::npos )
		{
			REPORT_ERROR( 200, "End of string reached" );
		}
		else
		{
			line = eraseChar( line, ' ' );
		}
		
		setPath( line );
		
		// Get Labels line
		getline( infile, line, '\n');
		
		// Remove ';'
		line.erase(0,line.find(";")+1);
				
		// Parse labels
		std::stringstream os( line );          
		std::string label;                 
		
		std::vector<std::string> labelsVector;
		
		while ( os >> label )
		{
			labelsVector.push_back( label );
		}
		
		// Read data and fill structures accordingly
		while ( getline( infile, line, '\n') )
		{
			long int objectID = addObject( );
			
			// Parse labels
			std::stringstream os2( line );          
			std::string value;
			
			int counter = 0;
			while ( os2 >> value )
			{
				setValue( labelsVector[counter], value );
				counter++;
			}
		}
	}
	
	infile.close( );
}

metaData::~metaData( )
{
	objects.clear( );
}

bool metaData::isEmpty( )
{
	return objects.empty( );
}

void metaData::clear( )
{
	path = std::string("");
	objects.clear( );		
}

void metaData::setPath( std::string newPath )
{
	if( newPath == "" )
	{  
		path = std::string( getcwd( NULL, 0 ) );
	}
	else
	{
		path = newPath;
	}
}

std::string metaData::getPath( )
{
	return path;
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

long int metaData::firstObject( ) 
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

long int metaData::nextObject( ) 
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

bool metaData::setValue( label name, double value, long int objectID )
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
		
		aux->addValue( name, value );
		
		return true;
	}	
	else
	{
		return false;
	}
}

bool metaData::setValue( label name, float value, long int objectID )
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
		
		aux->addValue( name, value );

		return true;
	}	
	else
	{
		return false;
	}	
}

bool metaData::setValue( label name, int value, long int objectID )
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
		
		aux->addValue( name, value );

		return true;
	}	
	else
	{
		return false;
	}	
}

bool metaData::setValue( label name, std::string value, long int objectID )
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
		
		aux->addValue( name, value );
		
		return true;
	}	
	else
	{
		return false;
	}	
}

bool metaData::setValue( std::string name, std::string value, long int objectID )
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
		
		aux->addValue( name, value );
		
		return true;
	}	
	else
	{
		return false;
	}	
}

std::vector<long int> metaData::findObjects( label name, double value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
	}
	
	return result;
}

std::vector<long int> metaData::findObjects( label name, float value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
	}
	
	return result;
}

std::vector<long int> metaData::findObjects( label name, int value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
	}
	
	return result;
}

std::vector<long int> metaData::findObjects( label name, bool value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
	}
	
	return result;
}

std::vector<long int> metaData::findObjects( label name, std::string value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, metaDataContainer *>::iterator It;
	
	metaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
			}
	
	return result;
}

double metaData::angleRot( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( ANGLEROT );
		
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

double metaData::angleTilt( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( ANGLETILT );
		
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

double metaData::anglePsi( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( ANGLEPSI );
		
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

bool metaData::enabled( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		bool * result = (bool *)aux->getValue( ENABLED );
		
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

double metaData::shiftX( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( SHIFTX );
		
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

double metaData::shiftY( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( SHIFTY );
		
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

double metaData::shiftZ( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( SHIFTZ );
		
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

double metaData::originX( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( ORIGINX );
		
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

double metaData::originY( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( ORIGINY );
		
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

double metaData::originZ( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( ORIGINZ );
		
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

std::string metaData::image( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		std::string * result = (std::string *)aux->getValue( IMAGE );
		
		if( result == NULL )
		{
			std::cerr << "No 'image' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

std::string metaData::CTFModel( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		std::string * result = (std::string *)aux->getValue( CTFMODEL );
		
		if( result == NULL )
		{
			std::cerr << "No 'CTFModel' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

std::string metaData::micrograph( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		std::string * result = (std::string *)aux->getValue( MICROGRAPH );
		
		if( result == NULL )
		{
			std::cerr << "No 'micrograph' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metaData::weight( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( WEIGHT );
		
		if( result == NULL )
		{
			std::cerr << "No 'weight' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}

double metaData::flip( long int objectID )
{
	if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	}
	else
	{
		metaDataContainer * aux = objects[ objectID ];
		double * result = (double *)aux->getValue( FLIP );
		
		if( result == NULL )
		{
			std::cerr << "No 'flip' label found for objectID = " << objectID << " . Exiting... " << std::endl;
		}
		else
		{
			return (*result);
		}
	}
}
