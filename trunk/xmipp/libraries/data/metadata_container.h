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

#ifndef METADATACONTAINER_H
#define METADATACONTAINER_H

#include <map>
#include "strings.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "funcs.h"

// This enum defines what MetaDataLabels this class can manage, if
// you need a new one add it here and modify affected methods:
//
// 	- static MetaDataLabel codifyLabel( std::string strLabel );
//	- static std::string decodeLabel( MetaDataLabel inputLabel );      
//	- void writeValuesToFile( std::ofstream &outfile, MetaDataLabel inputLabel );
//	- void addValue( std::string name, std::string value );
//
// Keep this special structure (using MDL_FIRSTLABEL and MDL_LAST_LABEL) so the
// programmer can iterate through it like this:
//
//  for( MetaDataLabel mdl = MDL_FIRST_LABEL ; mdl < MDL_LAST_LABEL ; MetaDataLabel( mdl+1 ) )
//

enum MetaDataLabel 
{ 
	MDL_UNDEFINED = -1,
    MDL_FIRST_LABEL,
	MDL_ANGLEROT = MDL_FIRST_LABEL,       // Rotation angle of an image (double)
    MDL_COMMENT,                          // A comment for this object /*** NOTE THIS IS A SPECIAL CASE AND SO IS TREATED ***/
	MDL_ANGLETILT,                        // Tilting angle of an image (double)
	MDL_ANGLEPSI,                         // Psi angle of an image (double)
	MDL_IMAGE,                            // Name of an image (std::string)
	MDL_MICROGRAPH,                       // Name of a micrograph (std::string)
	MDL_CTFMODEL,                         // Name for the CTF Model (std::string)
	MDL_SHIFTX,                           // Shift for the image in the X axis (double)
	MDL_SHIFTY,                           // Shift for the image in the Y axis (double)
	MDL_SHIFTZ,                           // Shift for the image in the Z axis (double)
	MDL_ENABLED,                          // Is this image enabled? (int [-1 or 1])
	MDL_ORIGINX,                          // Origin for the image in the X axis (double)
	MDL_ORIGINY,                          // Origin for the image in the Y axis (double)
	MDL_ORIGINZ,                          // Origin for the image in the Z axis (double)
	MDL_WEIGHT,                           // Weight assigned to the image (double)
	MDL_FLIP,                             // Flip the image? (bool)
	MDL_REF,                              // Class to which the image belongs (int)
	MDL_MAXCC,                            // Cross-correlation for the image (double)
        MDL_SERIE,                            // A collection of micrographs, e.g. a tilt serie (std::string)
        MDL_PMAX,                             // Maximum value of normalized probability function (now called "Pmax/sumP") (double)
        MDL_CTFINPUTPARAMS,                   // Parameters file for the CTF Model (std::string)
        MDL_PERIODOGRAM,                      // A periodogram's file name (std::string)
        MDL_NMA,                              // Normal mode displacements
	MDL_LAST_LABEL	                      // **** NOTE ****: Do keep this label always at the end
			 		      // it is here for looping purposes  	
};

inline bool isString(MetaDataLabel lCode)
{
    if(lCode == MDL_COMMENT     || lCode == MDL_IMAGE          || 
       lCode == MDL_MICROGRAPH  || 
       lCode == MDL_CTFMODEL    || lCode == MDL_CTFINPUTPARAMS || 
       lCode == MDL_PERIODOGRAM || lCode == MDL_SERIE)
        return true;
    else
        return false;
}

inline bool isDouble(MetaDataLabel lCode)
{
    if(lCode == MDL_ANGLEROT || lCode == MDL_ANGLETILT || lCode == MDL_ANGLEPSI ||\
       lCode == MDL_SHIFTX   || lCode == MDL_SHIFTY    || lCode == MDL_SHIFTZ   ||\
       lCode == MDL_ORIGINX  || lCode == MDL_ORIGINY   || lCode == MDL_ORIGINZ  ||\
       lCode == MDL_WEIGHT   || lCode == MDL_MAXCC     || lCode == MDL_PMAX)
        return true;
    else
        return false;
}

inline bool isVector(MetaDataLabel lCode)
{
    if(lCode==MDL_NMA)
        return true;
    else
        return false;
}

inline bool isBool(MetaDataLabel lCode)
{
    if(lCode==MDL_FLIP)
        return true;
    else
        return false;
}

inline bool isInt(MetaDataLabel lCode)
{
    if(lCode == MDL_REF || lCode == MDL_ENABLED)
        return true;
    else
        return false;
}

class MetaDataContainer
{
	/** Container for pairs "name" and value. Note that void * allows to use
	    mixed types */
	std::map<MetaDataLabel, void *> values;

	void insertVoidPtr( MetaDataLabel name, void * value );
	void * getVoidPtr( MetaDataLabel name );

	public:

	/**Assignment operator
	 *
	 */
	MetaDataContainer& operator = ( MetaDataContainer &MDc);
	
	/** Constructor */
	MetaDataContainer();
	/** Copy constructor
	 *
	 */
	MetaDataContainer( MetaDataContainer &MDc);

	/** Destructor */
	~MetaDataContainer();
	
	/** Create a new pair name-value of integer type */
	void addValue( MetaDataLabel name, double value );
	void addValue( MetaDataLabel name, int value );
	void addValue( MetaDataLabel name, bool value );
	void addValue( MetaDataLabel name, const std::string &value  );
	void addValue( MetaDataLabel name, const std::vector<double> &value  );
	void addValue( const std::string &name, const std::string &value );

	void getValue( MetaDataLabel name, int &value );
	void getValue( MetaDataLabel name, double &value );
	void getValue( MetaDataLabel name, std::string &value );
	void getValue( MetaDataLabel name, bool &value );
	void getValue( MetaDataLabel name, std::vector<double> &value );

	bool valueExists( MetaDataLabel name );

	bool pairExists( MetaDataLabel name, double value );
	bool pairExists( MetaDataLabel name, int value );
	bool pairExists( MetaDataLabel name, bool value );
	bool pairExists( MetaDataLabel name, const std::string &value );
        // pairExists for vectors makes nosense, not implemented

	void deleteValue( MetaDataLabel name );
	
	void writeValueToFile( std::ofstream &outfile, MetaDataLabel inputLabel );
    void writeValueToString( std::string &outString, MetaDataLabel inputLabel );

	static MetaDataLabel codifyLabel( std::string strLabel );
	static std::string decodeLabel( MetaDataLabel inputLabel );
	static bool isValidLabel( MetaDataLabel inputLabel );
	static bool isValidLabel( std::string inputLabel );
};

#endif
