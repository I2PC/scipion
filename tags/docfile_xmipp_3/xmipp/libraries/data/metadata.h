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

#ifndef METADATA_H
#define METADATA_H

#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "metadata_container.h"

// Possible error codes for the map
#define NO_OBJECTS_STORED	-1
#define NO_MORE_OBJECTS		-2

// Types that can be stored as metaData
#define UNDEFINED_TYPE		0
#define IMAGE_TYPE			1
#define MICROGRAPH_TYPE		2
#define VOLUME_TYPE			3

class metaData
	{
		std::map< long int, metaDataContainer *> objects;
		std::map< long int, metaDataContainer *>::iterator objectsIterator;
		
		std::string programName;
		std::string path;
		unsigned int objectsType;
		
	public:
		
		metaData();
		metaData( std::string fileName, std::string type, std::vector<std::string> * labelsVector = NULL );
		~metaData();
		
		bool isEmpty( );
		
		void clear( );
		
		long int addObject( );
		long int firstObject( );
		long int nextObject( );
		long int lastObject( );
		
		// Set a new pair/value for an specified object. If no objectID is given, that
		// pointed by the class iterator is used 
		bool setValue( std::string label, double value, long int objectID = -1 );
		bool setValue( std::string label, float value, long int objectID = -1 );
		bool setValue( std::string label, int value, long int objectID = -1 );
		bool setValue( std::string label, bool, long int objectID = -1 );
		bool setValue( std::string label, std::string value, long int objectID = -1 );
		
		// Get the collection of objects whose pair label/value is given
		std::vector<long int> findObjects( std::string label, double value );
		std::vector<long int> findObjects( std::string label, float value );
		std::vector<long int> findObjects( std::string label, int value );
		std::vector<long int> findObjects( std::string label, bool value );
		std::vector<long int> findObjects( std::string label, std::string value );
		
		// Xmipp-specific, for new parameters add here.
		void setProgram( std::string newProgName );		
		void setPath( std::string newPath );	
		void setType( unsigned int newObjectsType );
		
		std::string getProgram( );
		std::string getPath( );
		unsigned int getType( );
		
		double rot( long int objectID );
		double tilt( long int objectID );
		double psi( long int objectID );
		double shiftX( long int objectID );
		double shiftY( long int objectID );
		double shiftZ( long int objectID );
		double originX( long int objectID );
		double originY( long int objectID );
		double originZ( long int objectID );
		bool enabled( long int objectID );
		std::string fileName( long int objectID );

	};

#endif
