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
#include <iostream>
#include <iterator>
#include <sstream>
#include "funcs.h"
#include "strings.h"

#include "metadata_container.h"

// Possible error codes for the map
#define NO_OBJECTS_STORED	-1
#define NO_MORE_OBJECTS		-2

class MetaData
{
	std::map< long int, MetaDataContainer *> objects;
	std::map< long int, MetaDataContainer *>::iterator objectsIterator;

	std::map< std::string, long int> fastStringSearch;
	label fastStringSearchLabel;

	std::vector< label > activeLabels;

	std::string path;

public:
	
	MetaData();
	MetaData( std::string fileName, std::vector<label> * labelsVector = NULL );
	~MetaData();
	
	void save( std::string fileName );
	
	bool isEmpty( );
	
	void clear( );
	
	long int addObject( );
	long int firstObject( );
	long int nextObject( );
	long int lastObject( );
	
	long int fastSearch( label name, std::string value, bool recompute = false );
	
	// Set a new pair/value for an specified object. If no objectID is given, that
	// pointed by the class iterator is used 
	bool setValue( label name, double value, long int objectID = -1 );
	bool setValue( label name, float value, long int objectID = -1 );
	bool setValue( label name, int value, long int objectID = -1 );
	bool setValue( label name, bool, long int objectID = -1 );
	bool setValue( label name, std::string value, long int objectID = -1 );
	
	bool setValue( std::string name, std::string value, long int objectID = -1 );
	
	// Get the collection of objects whose pair label/value is given
	std::vector<long int> findObjects( label name, double value );
	std::vector<long int> findObjects( label name, float value );
	std::vector<long int> findObjects( label name, int value );
	std::vector<long int> findObjects( label name, bool value );
	std::vector<long int> findObjects( label name, std::string value );
	
	// Xmipp-specific, for new parameters add here.
	void setPath( std::string newPath = "" );	
	
	std::string getPath( );
	
	double angleRot( long int objectID );
	double angleTilt( long int objectID );
	double anglePsi( long int objectID );
	double shiftX( long int objectID );
	double shiftY( long int objectID );
	double shiftZ( long int objectID );
	double originX( long int objectID );
	double originY( long int objectID );
	double originZ( long int objectID );
	double weight( long int objectID );
	double flip( long int objectID );
	double maxCC( long int objectID );
	int ref( long int objectID );
	bool enabled( long int objectID );
	std::string CTFModel( long int objectID );
	std::string image( long int objectID );
	std::string micrograph( long int objectID );
};

#endif
