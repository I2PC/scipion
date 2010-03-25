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

MetaData::MetaData()
{
	setPath();
	objects.clear( );
	fastStringSearchLabel = MDL_UNDEFINED;	
    objectsIterator = objects.begin();
}

void MetaData::read( std::ifstream *infile, bool skipDisabled, std::vector<MetaDataLabel> * labelsVector )
{
    infile->seekg(0, std::ios::beg); 
	std::string line;

	getline( *infile, line, '\n');
	
    int pos = line.find( "*" );
		
	if( pos == std::string::npos )
    {
		REPORT_ERROR( 200, "End of string reached" );
	}
	else
	{
		line.erase( 0, pos+1 );
		line = removeChar( line, ' ' );
	}
	
	setPath( line );
		
	// Get Labels line
	getline( *infile, line, '\n');
		
	// Remove ';'
	line.erase(0,line.find(";")+1);
				
	// Parse labels
	std::stringstream os( line );          
	std::string newLabel;                 
				
	while ( os >> newLabel )
	{
		activeLabels.push_back(  MetaDataContainer::codifyLabel(newLabel) );
	}
		
	// Read data and fill structures accordingly
	while ( getline( *infile, line, '\n') )
	{
		long int objectID = addObject( );
		
		// Parse labels
		std::stringstream os2( line );          
		std::string value;
		
		int counter = 0;
		while ( os2 >> value )
		{
			setValue( MetaDataContainer::decodeLabel(activeLabels[counter]), value );
			counter++;
		}
	}
}

void MetaData::readOldDocFile( std::ifstream *infile, bool skipDisabled, std::vector<MetaDataLabel> * labelsVector )
{
    bool saveName = true;
    infile->seekg(0, std::ios::beg);     
    std::string line;
    
	// Search for Headerinfo, if present we are processing an old-styled docfile
	// else we are processing a new Xmipp MetaData file
	getline( *infile, line, '\n');
	
	int pos = line.find( "Headerinfo" );
    
    // Remove from the beginning to the end of "Headerinfo columns:"
    line = line.erase( 0, line.find( ":" ) + 1 );
           
    // In the old docfile format the "image" label did not exist, it was
    // a ";" commented line containing the name of the projection. Therefore,
    // for the new format it must be added by hand, if necessary
    if( labelsVector != NULL )
    {   
        std::vector< MetaDataLabel >::iterator location;
            		
        location = std::find( labelsVector->begin(), labelsVector->end(), MDL_IMAGE );
            
        if ( location != labelsVector->end() )
        {
           	activeLabels.push_back( MDL_IMAGE );
           	saveName = true;
        }
    }   
    else
    {				
        activeLabels.push_back( MDL_IMAGE );
    }   
            				
    // Extract labels until the string is empty
    while ( line != "" )
    {   
        pos = line.find( ")" );
        std::string newLabel = line.substr( 0, pos+1 );
        line.erase( 0, pos+1 );
            
        // The token can now contain a ',', if so, remove it
        if( ( pos = newLabel.find( "," ) ) != std::string::npos )
         	newLabel.erase( pos, 1 );
            
        // Remove unneded parentheses and contents
        pos = newLabel.find( "(" );
        newLabel.erase( pos, newLabel.find( ")" )-pos+1 );
            
        // Remove white spaces
        newLabel = removeChar( newLabel, ' ' );
            
        if( labelsVector != NULL )
        {
            std::vector< MetaDataLabel >::iterator location;
            	
            location = std::find( labelsVector->begin(), labelsVector->end(), MetaDataContainer::codifyLabel( newLabel ) );

            if ( location != labelsVector->end() )
            {
            	activeLabels.push_back( MetaDataContainer::codifyLabel(newLabel) );
            }
        }
        else
        {
          	activeLabels.push_back( MetaDataContainer::codifyLabel(newLabel) );
        }
    }   
            	
    int isname=0;
    while ( getline( *infile, line, '\n') )
    {   
        if( isname % 2 == 0 )
        {
          	long int objectID = addObject( );
           	line.erase(0,line.find(";")+1);
            	
           	// Remove spaces from string
           	line = removeChar( line, ' ' );
           					
           	setValue( MDL_IMAGE, line );
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
           		{
           			if( saveName )
           				setValue( MetaDataContainer::decodeLabel(activeLabels[counter-1]), value );
           			else
           				setValue( MetaDataContainer::decodeLabel(activeLabels[counter-2]), value );
           		}
           		counter++;
           	}
        }
            
        isname++;
    }   
}

void MetaData::readOldSelFile( std::ifstream *infile, bool skipDisabled )
{	
    infile->seekg(0, std::ios::beg);     
    std::string line;

	getline( *infile, line, '\n');
	while ( getline( *infile, line, '\n') )
    {
	    line=simplify(line);
	    if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
	         continue;
	    else
	    {
	        int pos = line.find( " " );
	        std::string name=line.substr( 0, pos );
            line.erase(0,pos+1);
            int i = atoi (line.c_str());
            if (skipDisabled && i==(-1))
                continue;
            addObject();
            setValue(MDL_IMAGE,name);
            setValue(MDL_ENABLED,i);
        } 
    }	
}

MetaData::MetaData( std::string fileName, bool skipDisabled, std::vector<MetaDataLabel> * labelsVector )
{
	setPath();
	objects.clear( );
	fastStringSearchLabel = MDL_UNDEFINED;	
    objectsIterator = objects.begin();

	// Open file
	std::ifstream infile ( fileName.data(), std::ios_base::in );
	std::string line;
	
	// Search for Headerinfo, if present we are processing an old-styled docfile
	// else we are processing a new Xmipp MetaData file
	getline( infile, line, '\n');
	
	int pos = line.find( "Headerinfo" );
	
	if( pos != std::string::npos ) // Headerinfo token found
	{
        readOldDocFile( &infile , skipDisabled, labelsVector );
	}
	else
	{
        getline( infile, line, '\n');
        
        pos = line.find( "XMIPP_3 * " );
	    
        if( pos != std::string::npos ) // xmipp_3 token found
        {
            read( &infile, skipDisabled, labelsVector );
	    }
        else    // We are reading an old selfile
        {
            readOldSelFile( &infile, skipDisabled );   
        }
    }
	
	infile.close( );
}

void MetaData::write( std::string fileName )
{
	// Open file
	std::ofstream outfile ( fileName.data(), std::ios_base::out );

	outfile << "; ";
	outfile << "XMIPP_3 * ";
	outfile << path << std::endl;
	
	std::map< long int, MetaDataContainer *>::iterator It;
	std::vector< MetaDataLabel >::iterator strIt;
	
	outfile << "; ";	
	for( strIt = activeLabels.begin( ); strIt != activeLabels.end( ); strIt ++ )
	{
		outfile << MetaDataContainer::decodeLabel(*strIt);
		outfile << " ";
	}
	outfile << std::endl;
	
	for( It = objects.begin( ); It != objects.end(); It ++)
	{
		for( strIt = activeLabels.begin( ); strIt != activeLabels.end( ); strIt ++ )
		{
			(It->second)->writeValueToFile( outfile, *strIt );
			outfile << " ";
		}
		
		outfile << std::endl;
	} 
}

MetaData::~MetaData( )
{
	clear( );
}

bool MetaData::isEmpty( )
{
	return objects.empty( );
}

void MetaData::clear( )
{
	path = std::string("");
	objects.clear( );		
	
	objectsIterator = objects.end();
	
	fastStringSearch.clear( );;
	fastStringSearchLabel = MDL_UNDEFINED;

	activeLabels.clear( );

	path = std::string(""); 
}

void MetaData::setPath( std::string newPath )
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

std::string MetaData::getPath( )
{
	return path;
}

long int MetaData::addObject( )
{
	long int result = lastObject();
	typedef std::pair<long int, MetaDataContainer *> newPair;
		
	if( result == NO_OBJECTS_STORED )
	{
		result = 0;
	}
	else
	{
		result++;
	}

	objects.insert( newPair( result, new MetaDataContainer() ) );
	
	// Set iterator pointing to the newly added object
	objectsIterator = objects.end( );
	objectsIterator--;

	return result;
}

long int MetaData::firstObject( ) 
{ 
	long int result = 0;
	
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

long int MetaData::nextObject( ) 
{ 
	long int result = 0;
	
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

long int MetaData::lastObject( )
{
	long int result = 0;
	
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

bool MetaData::setValue( MetaDataLabel name, double value, long int objectID )
{
	long int auxID;
	    
	if( !objects.empty( ) && MetaDataContainer::isValidLabel( name ))
	{
    	if( objectID == -1 )
		{
			auxID = objectsIterator->first;
    	}
		else
		{
			auxID = objectID;
		}
		
		MetaDataContainer * aux = objects[auxID];
		
		// Check whether label is correct (belongs to the enum in the metadata_container header
		// and whether it is present in the activeLabels vector. If not, add it to all the other
		// objects with default values
		std::vector< MetaDataLabel >::iterator location;
		std::map< long int, MetaDataContainer *>::iterator It;
									
   		location = std::find( activeLabels.begin(), activeLabels.end(), name );
		
	   	if ( location == activeLabels.end() )
		{
			activeLabels.push_back( name );
			
			// Add this label to the rest of the objects in this class
			for( It = objects.begin( ); It != objects.end(); It ++)
			{
				if( It->second != aux )
				{
					(It->second)->addValue( name, double() );
				}		
			} 
			
		}
			
		aux->addValue( name, value );
		
		return true;
	}	
	else
	{
		return false;
	}
}

bool MetaData::setValue( MetaDataLabel name, float value, long int objectID )
{
	long int auxID;
	
	if( !objects.empty( ) && MetaDataContainer::isValidLabel( name ))
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		MetaDataContainer * aux = objects[auxID];
		
		// Check whether label is correct (belongs to the enum in the metadata_container header
		// and whether it is present in the activeLabels vector. If not, add it to all the other
		// objects with default values
		std::vector< MetaDataLabel >::iterator location;
		std::map< long int, MetaDataContainer *>::iterator It;
									
   		location = std::find( activeLabels.begin(), activeLabels.end(), name );
		
	   	if ( location == activeLabels.end() )
		{
			activeLabels.push_back( name );
			
			// Add this label to the rest of the objects in this class
			for( It = objects.begin( ); It != objects.end(); It ++)
			{
				if( It->second != aux )
				{
					(It->second)->addValue( name, float() );
				}		
			} 
			
		}
			
		aux->addValue( name, value );

		return true;
	}	
	else
	{
		return false;
	}	
}

bool MetaData::setValue( MetaDataLabel name, int value, long int objectID )
{
	long int auxID;
	
	if( !objects.empty( ) && MetaDataContainer::isValidLabel( name ))
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		MetaDataContainer * aux = objects[auxID];
		
		// Check whether label is correct (belongs to the enum in the metadata_container header
		// and whether it is present in the activeLabels vector. If not, add it to all the other
		// objects with default values
		std::vector< MetaDataLabel >::iterator location;
		std::map< long int, MetaDataContainer *>::iterator It;
									
   		location = std::find( activeLabels.begin(), activeLabels.end(), name );
		
	   	if ( location == activeLabels.end() )
		{
			activeLabels.push_back( name );
			
			// Add this label to the rest of the objects in this class
			for( It = objects.begin( ); It != objects.end(); It ++)
			{
				if( It->second != aux )
				{
					(It->second)->addValue( name, int() );
				}		
			} 
			
		}
		
		aux->addValue( name, value );

		return true;
	}	
	else
	{
		return false;
	}	
}

bool MetaData::setValue( MetaDataLabel name, bool value, long int objectID )
{
	long int auxID;

	if( !objects.empty( ) && MetaDataContainer::isValidLabel( name ))
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		MetaDataContainer * aux = objects[auxID];

		// Check whether label is correct (belongs to the enum in the metadata_container header
		// and whether it is present in the activeLabels vector. If not, add it to all the other
		// objects with default values
		std::vector< MetaDataLabel >::iterator location;
		std::map< long int, MetaDataContainer *>::iterator It;
									
   		location = std::find( activeLabels.begin(), activeLabels.end(), name );
		
	   	if ( location == activeLabels.end() )
		{
			activeLabels.push_back( name );
			
			// Add this label to the rest of the objects in this class
			for( It = objects.begin( ); It != objects.end(); It ++)
			{
				if( It->second != aux )
				{
					(It->second)->addValue( name, bool( ) );
				}		
			} 
			
		}
		
		aux->addValue( name, value );
		
		return true;
	}	
	else
	{
		return false;
	}	
}

bool MetaData::setValue( MetaDataLabel name, std::string value, long int objectID )
{
	long int auxID;

	if( !objects.empty( ) && MetaDataContainer::isValidLabel( name ))
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		MetaDataContainer * aux = objects[auxID];

		// Check whether label is correct (belongs to the enum in the metadata_container header
		// and whether it is present in the activeLabels vector. If not, add it to all the other
		// objects with default values
		std::vector< MetaDataLabel >::iterator location;
		std::map< long int, MetaDataContainer *>::iterator It;
									
   		location = std::find( activeLabels.begin(), activeLabels.end(), name );
		
	   	if ( location == activeLabels.end() )
		{
			activeLabels.push_back( name );
			
			// Add this label to the rest of the objects in this class
			for( It = objects.begin( ); It != objects.end(); It ++)
			{
				if( It->second != aux )
				{
					(It->second)->addValue( name, std::string( "" ) );
				}		
			} 
			
		}
		
		aux->addValue( name, value );
		
		return true;
	}	
	else
	{
		return false;
	}	
}

bool MetaData::setValue( std::string name, std::string value, long int objectID )
{
	long int auxID;
	
	MetaDataLabel label = MetaDataContainer::codifyLabel( name );
	
	if( !objects.empty( ) && MetaDataContainer::isValidLabel( label ) )
	{
		if( objectID == -1 )
		{
			auxID = objectsIterator->first;
		}
		else
		{
			auxID = objectID;
		}
		
		MetaDataContainer * aux = objects[auxID];

		// Check whether label is correct (belongs to the enum in the metadata_container header
		// and whether it is present in the activeLabels vector. If not, add it to all the other
		// objects with default values
		std::vector< MetaDataLabel >::iterator location;
		std::map< long int, MetaDataContainer *>::iterator It;
									
   		location = std::find( activeLabels.begin(), activeLabels.end(), label );
		
	   	if ( location == activeLabels.end() )
		{
			activeLabels.push_back( label );
			
			// Add this label to the rest of the objects in this class
			for( It = objects.begin( ); It != objects.end(); It ++)
			{
				if( It->second != aux )
				{
					(It->second)->addValue( label, std::string( "" ) );
				}		
			} 
			
		}
		
		aux->addValue( name, value );
		
		return true;
	}	
	else
	{
		return false;
	}	
}

std::vector<long int> MetaData::findObjects( MetaDataLabel name, double value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, MetaDataContainer *>::iterator It;
	
	MetaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
	}
	
	return result;
}

std::vector<long int> MetaData::findObjects( MetaDataLabel name, float value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, MetaDataContainer *>::iterator It;
	
	MetaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
	}
	
	return result;
}

std::vector<long int> MetaData::findObjects( MetaDataLabel name, int value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, MetaDataContainer *>::iterator It;
	
	MetaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
	}
	
	return result;
}

std::vector<long int> MetaData::findObjects( MetaDataLabel name, bool value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, MetaDataContainer *>::iterator It;
	
	MetaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
	}
	
	return result;
}

std::vector<long int> MetaData::findObjects( MetaDataLabel name, std::string value )
{
	std::vector<long int> result;
	
	// Traverse all the structure looking for objects
	// that satisfy search criteria
	
	std::map< long int, MetaDataContainer *>::iterator It;
	
	MetaDataContainer * aux;
	
	for( It = objects.begin( ) ; It != objects.end( ); It ++ )
	{
		aux = It->second;
		
		if( aux->pairExists( name, value ) )
			result.push_back( It->first );
			}
	
	return result;
}
double MetaData::angleRot( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_ANGLEROT );
		
	if( result == NULL )
	{
		std::cerr << "No 'angleRot' label found for objectID = " << objectID << " . Exiting... " << std::endl;
        exit( 1 );
	}

	return (*result);
}

double MetaData::angleTilt( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_ANGLETILT );
		
	if( result == NULL )
	{
		std::cerr << "No 'angleTilt' label found for objectID = " << objectID << " . Exiting... " << std::endl;
        exit( 1 );
	}

	return (*result);
}

double MetaData::anglePsi( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_ANGLEPSI );
		
	if( result == NULL )
	{
		std::cerr << "No 'anglePs' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}

	return (*result);
}

int MetaData::enabled( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    int * result = (int *)aux->getValue( MDL_ENABLED );
		
	if( result == NULL )
	{
		std::cerr << "No 'enabled' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}	
			
    return (*result);
}

double MetaData::shiftX( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_SHIFTX );
	
	if( result == NULL )		
    {
		std::cerr << "No 'siftX' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}
	
    return (*result);
}

double MetaData::shiftY( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_SHIFTY );
		
	if( result == NULL )
	{
		std::cerr << "No 'shiftY' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}

	return (*result);
}

double MetaData::shiftZ( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_SHIFTZ );
	
	if( result == NULL )
    {
		std::cerr << "No 'shiftZ' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}	

	return (*result);
}

double MetaData::originX( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_ORIGINX );
		
	if( result == NULL )
	{
		std::cerr << "No 'originX' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}
	
    return (*result);
}

double MetaData::originY( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_ORIGINY );
		
	if( result == NULL )
	{
		std::cerr << "No 'originY' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}
	
    return (*result);
}

double MetaData::originZ( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_ORIGINZ );
		
	if( result == NULL )
	{
		std::cerr << "No 'originZ' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}

	return (*result);
}

double MetaData::pMax( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_PMAX );
		
	if( result == NULL )
	{
		std::cerr << "No 'pMax' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}
	
    return (*result);
}


std::string MetaData::image( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
		
    std::string * result = (std::string *)aux->getValue( MDL_IMAGE );
		
	if( result == NULL )
	{
		std::cerr << "No 'image' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}

	return (*result);
}

std::string MetaData::CTFModel( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    std::string * result = (std::string *)aux->getValue( MDL_CTFMODEL );
		
	if( result == NULL )
	{
		std::cerr << "No 'CTFModel' label found for objectID = " << objectID << " . Exiting... " << std::endl;
       	exit( 1 );
	}
	
    return (*result);
}

std::string MetaData::micrograph( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    std::string * result = (std::string *)aux->getValue( MDL_MICROGRAPH );
		
	if( result == NULL )
	{
		std::cerr << "No 'micrograph' label found for objectID = " << objectID << " . Exiting... " << std::endl;
   	    exit( 1 );
	}

    return (*result);
}

std::string MetaData::CTFInputParams( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    std::string * result = (std::string *)aux->getValue( MDL_CTFINPUTPARAMS );
		
	if( result == NULL )
	{
		std::cerr << "No 'CTFInputParam' label found for objectID = " << objectID << " . Exiting... " << std::endl;
        exit( 1 );
	}

	return (*result);
}

std::string MetaData::periodogram( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    std::string * result = (std::string *)aux->getValue( MDL_PERIODOGRAM );
		
	if( result == NULL )
	{
		std::cerr << "No 'periodogram' label found for objectID = " << objectID << " . Exiting... " << std::endl;
	    exit( 1 );
	}

	return (*result);
}

std::string MetaData::serie( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    std::string * result = (std::string *)aux->getValue( MDL_SERIE );
		
	if( result == NULL )
	{
		std::cerr << "No 'serie' label found for objectID = " << objectID << " . Exiting... " << std::endl;
	    exit( 1 );
 	}
    
    return (*result);
}

double MetaData::weight( long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_WEIGHT );
	
	if( result == NULL )	
    {
		std::cerr << "No 'weight' label found for objectID = " << objectID << " . Exiting... " << std::endl;
	    exit( 1 );
	}

	return (*result);
}

double MetaData::flip( long int objectID )
{
    MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_FLIP );
        
	if( result == NULL )
	{
		std::cerr << "No 'flip' label found for objectID = " << objectID << " . Exiting... " << std::endl;
        exit( 1 );
	}

	return (*result);
}

double MetaData::maxCC( long int objectID )
{
    MetaDataContainer * aux = getObject( objectID );
    
    double * result = (double *)aux->getValue( MDL_MAXCC );
		
	if( result == NULL )
	{
		std::cerr << "No 'maxCC' label found for objectID = " << objectID << " . Exiting... " << std::endl;
	    exit( 1 );
    }

	return (*result);
}

int MetaData::ref( long int objectID )
{
    MetaDataContainer * aux = getObject( objectID );
	
    int * result = (int *)aux->getValue( MDL_REF );
		
	if( result == NULL )
	{
		std::cerr << "No 'ref' label found for objectID = " << objectID << " . Exiting... " << std::endl;
	    exit( 1 );
    }
	
    return (*result);
}

void MetaData::setAngleRot( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_ANGLEROT, value );
}

void MetaData::setAngleTilt( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_ANGLETILT, value );
}

void MetaData::setAnglePsi( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_ANGLEPSI, value );
}

void MetaData::setEnabled( int value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_ENABLED, value );
}

void MetaData::setShiftX( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_SHIFTX, value );
}

void MetaData::setShiftY( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_SHIFTY, value );
}

void MetaData::setShiftZ( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_SHIFTZ, value );
}

void MetaData::setOriginX( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_ORIGINX, value );
}

void MetaData::setOriginY( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_ORIGINY, value );
}

void MetaData::setOriginZ( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_ORIGINZ, value );
}

void MetaData::setPMax( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_PMAX, value );
}


void MetaData::setImage( std::string value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_IMAGE, value );
}

void MetaData::setCTFModel( std::string value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_CTFMODEL, value );
}

void MetaData::setMicrograph( std::string value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_MICROGRAPH, value );
}

void MetaData::setCTFInputParams( std::string value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_CTFINPUTPARAMS, value );
}

void MetaData::setPeriodogram( std::string value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_PERIODOGRAM, value );
}

void MetaData::setSerie( std::string value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_SERIE, value );
}

void MetaData::setWeight( double value, long int objectID )
{
	MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_WEIGHT, value );
}

void MetaData::setFlip( bool value, long int objectID )
{
    MetaDataContainer * aux = getObject( objectID );
    aux->addValue( MDL_FLIP, value );
}

void MetaData::setMaxCC( double value, long int objectID )
{
    MetaDataContainer * aux = getObject( objectID );   
    aux->addValue( MDL_MAXCC, value );
}

void MetaData::setRef( int value, long int objectID )
{
    MetaDataContainer * aux = getObject( objectID );	
    aux->addValue( MDL_REF, value );
}

MetaDataContainer * MetaData::getObject( long int objectID )
{
    MetaDataContainer * aux;
    
    if( objectID == -1 )
    {
        aux = objectsIterator->second;
    }
	else if( objects.find( objectID ) == objects.end( ) )
	{
		// This objectID does not exist, finish execution
		std::cerr << "No objectID = " << objectID << " found. Exiting... " << std::endl;
	    exit( 1 );
    }
	else
	{
		aux = objects[ objectID ];
	}
 
    return aux;  
}

long int MetaData::fastSearch( MetaDataLabel name, std::string value, bool recompute )
{
	long int result;
	
	if( recompute || fastStringSearch.empty( ) || fastStringSearchLabel != name )
	{
		fastStringSearch.clear( );
		fastStringSearchLabel = name;
		
		// Repopulate list
		// Traverse all the structure looking for objects
		// that satisfy search criteria
	
		std::map< long int, MetaDataContainer *>::iterator It;
	
		MetaDataContainer * aux;
	
		for( It = objects.begin( ) ; It != objects.end( ); It ++ )
		{
			aux = It->second;
		
			if( aux->valueExists( name ) )
			{
				fastStringSearch[ *((std::string *)((It->second)->getValue( name ))) ] = It->first ;
				
				if( aux->pairExists( name, value) )
				{
					result = It->first ;
				}
			}
		}
	}
	else
	{
		std::map< std::string, long int>::iterator It;
	
		if( ( It = fastStringSearch.find( value ) ) != fastStringSearch.end( ) )
		{
			result = It->second;
		}
		else
		{
			result = -1; // Not found
		}
	}
	
	return result;
}
		
