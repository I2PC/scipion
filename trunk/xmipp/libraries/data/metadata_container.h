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

#ifndef CONTAINER_H
#define CONTAINER_H

#include <map>
#include <string>
#include <iostream>

#define angle_t double

#define rot(object) *((angle_t*)(object->getValue(std::string("rot"))))
#define tilt(object) *((angle_t*)(object->getValue(std::string("tilt"))))
#define psi(object) *((angle_t*)(object->getValue(std::string("psi"))))

class xmpContainer
{
	/** Container for pairs "name" and value. Note that void * allows to use
	    mixed types */
	std::map<std::string, void *> values;
	
	int insertVoidPtr( std::string name, void * value );

	public:
	
	/** Constructor */
	xmpContainer();
	
	/** Destructor */
	~xmpContainer();
	
	/** Create a new pair name-value of integer type */
	int addValue( std::string name, int value );
	
	/** Create a new pair name-value of double type */
	int addValue( std::string name, double value );
	
	void * getValue( std::string name );
	bool valueExists( std::string name );
	void deleteValue( std::string name );
};

//@}
#endif
