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
#include "funcs.h"
#include "metadata_container.h"
#include "image.h"
#include <time.h>
#include <stdio.h>

class MetaData
{

	std::map< long int, MetaDataContainer *> objects;
	
	// Used by firstObject, nextObject and lastObject to keep a pointer
	// to the "active" object. This way when you call setValue without
	// an objectID, this one is chosen
	std::map< long int, MetaDataContainer *>::iterator objectsIterator;
	
	// Allows a fast search for pairs where the value is
	// a string, i.e. looking for filenames which is quite
	// usual
	std::map< std::string, long int> fastStringSearch;
	MetaDataLabel fastStringSearchLabel;


	std::string path;
	std::string comment;
    
    MetaDataContainer * getObject( long int objectID = -1 );

	// Set a new pair/value for an specified object. If no objectID is given, that
	// pointed by the class iterator is used 
	bool setValue( MetaDataLabel name, double value, long int objectID = -1 );
	bool setValue( MetaDataLabel name, float value, long int objectID = -1 );
	bool setValue( MetaDataLabel name, int value, long int objectID = -1 );
	bool setValue( MetaDataLabel name, bool, long int objectID = -1 );
	bool setValue( MetaDataLabel name, std::string value, long int objectID = -1 );
	bool setValue( std::string name, std::string value, long int objectID = -1 );
    

    void readOldSelFile( std::ifstream *infile );
    void readOldDocFile( std::ifstream *infile, std::vector<MetaDataLabel> * labelsVector );
    void read( std::ifstream *infile, std::vector<MetaDataLabel> * labelsVector );

    bool isColumnFormat;
    /**Input file name
     *
     */
    FileName infile;
public:
	/** What labels have been read from a docfile/metadata file
	* and/or will be stored on a new metadata file when "save" is
	* called
	* */
	std::vector< MetaDataLabel > activeLabels;

	long int addObject( );
    void read( FileName infile, std::vector<MetaDataLabel> * labelsVector=NULL );
    /** Overloading of the read function for python
     *
     */
    //void read( std::string infile){read(infile,NULL);}

	// Possible error codes for the map
    enum errors
    {
        NO_OBJECTS_STORED = -1, // NOTE: Do not change this value (-1)
        NO_MORE_OBJECTS = -2
    };
    
	long int firstObject( );
	long int nextObject( );
	long int lastObject( );
	
	MetaData();
	MetaData( FileName fileName, std::vector<MetaDataLabel> * labelsVector = NULL );
	~MetaData();
	
	void write( std::string fileName );
	
	bool isEmpty( );
	
	void clear( );
	
	long int fastSearch( MetaDataLabel name, std::string value, bool recompute = false );
	
	// Get the collection of objects whose pair label/value is given
	std::vector<long int> findObjects( MetaDataLabel name, double value );
	std::vector<long int> findObjects( MetaDataLabel name, float value );
	std::vector<long int> findObjects( MetaDataLabel name, int value );
	std::vector<long int> findObjects( MetaDataLabel name, bool value );
	std::vector<long int> findObjects( MetaDataLabel name, std::string value );
	
    // Removes the collection of objects whose pair label/value is given
    // NOTE: The iterator will point to the first object after any of these
    // operations.
	void removeObjects( MetaDataLabel name, double value );
	void removeObjects( MetaDataLabel name, float value );
	void removeObjects( MetaDataLabel name, int value );
	void removeObjects( MetaDataLabel name, bool value );
	void removeObjects( MetaDataLabel name, std::string value );
	
    // Xmipp-specific, for new parameters add here.
	void setPath( std::string newPath = "" );	
	void setComment( std::string Comment = "" );
	
	std::string getPath( );
	std::string getComment( );
	
//    void getAngles( double &rot, double &tilt, double &psi, long int objectID = -1 );
    
	double angleRot( long int objectID = -1 );
	double angleTilt( long int objectID = -1 );
	double anglePsi( long int objectID = -1 );
	double shiftX( long int objectID = -1 );
	double shiftY( long int objectID = -1 );
	double shiftZ( long int objectID = -1 );
	double originX( long int objectID = -1 );
	double originY( long int objectID = -1 );
	double originZ( long int objectID = -1 );
	double weight( long int objectID = -1 );
	double flip( long int objectID = -1 );
	double maxCC( long int objectID = -1);
    double pMax( long int objectID = -1 );
	int ref( long int objectID = -1 );
	int enabled( long int objectID = -1 );
	std::string CTFModel( long int objectID = -1 );
	std::string image( long int objectID = -1 );
	std::string micrograph( long int objectID = -1 );
	std::string CTFInputParams( long int objectID = -1 );
	std::string periodogram( long int objectID = -1 );
	std::string serie( long int objectID = -1 );
    
//    void setAngles( double rot, double tilt, double psi, long int objectID = -1 );
    
    void setAngleRot( double value, long int objectID = -1 );
    void setAngleTilt( double value, long int objectID = -1 );
    void setAnglePsi( double value, long int objectID = -1 );
    void setEnabled( int value, long int objectID = -1 );
    void setShiftX( double value, long int objectID = -1 );
    void setShiftY( double value, long int objectID = -1 );
    void setShiftZ( double value, long int objectID = -1 );
    void setOriginX( double value, long int objectID = -1 );
    void setOriginY( double value, long int objectID = -1 );
    void setOriginZ( double value, long int objectID = -1 );
    void setPMax( double value, long int objectID = -1 );
    void setImage( std::string value, long int objectID = -1 );
    void setCTFModel( std::string value, long int objectID = -1 );
    void setMicrograph( std::string value, long int objectID = -1 );
    void setCTFInputParams( std::string value, long int objectID = -1 );
    void setPeriodogram( std::string value, long int objectID = -1 );
    void setSerie( std::string value, long int objectID = -1 );
    void setWeight( double value, long int objectID = -1 );
    void setFlip( bool value, long int objectID = -1 );
    void setMaxCC( double value, long int objectID = -1 );
    void setRef( int value, long int objectID = -1 );
    size_t size(void){return objects.size();}
    /** Return metafile filename
     *
     */

    FileName getFilename(){return (infile);}
    /*Detect is there is at least one entry with the given label,entry pair
     * So far only implemented for int
     */
    bool detectObjects( MetaDataLabel name, int value );
    /**Count number of objects that satisfy a given label,entry pair
     * So far only implemented for int
     */
    long int countObjects( MetaDataLabel name, int value );

};
/** Compute images metadata estatistics
 * This use to be part of Metadata but should not
 */

void get_statistics(MetaData MT,Image& _ave, Image& _sd, double& _min,
                             double& _max, bool apply_geo);

#endif
