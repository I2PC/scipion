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

/// @defgroup Containers Pairs Container
/// @ingroup DataLibrary

/** A container for heterogeneous types ordered by label
 * @ingroup Containers
 *
 * This class provides a way to store heterogeneous data types on the
 * same data structure. This is a flexible way for extending Xmipp 
 * features and capabilities.
 *
 * Data is ordered by means of a pair "label"->value. Each "label" must
 * be unique and should have some useful meaning (i.e. "rot" for a rotation
 * angle).
 *
 * To hide type return conversions for certain well-known labels, the 
 * programmer should write macros in the container.h file as the next examples:
 *
 * @code
 * #define rot(object) *((double*)(object->getValue(std::string("rot"))))
 * #define dims(object) *((int*)(object->getValue(std::string("dims"))))
 * @endcode
 *
 * where "object" is the container itself.
 *
 * @code
 * xmpContainer * params = new xmpContainer( );
 *
 * // Adds an int value
 * params->addValue( std::string("dims"), 5);
 * // Adds a double value
 * params->addValue( std::string("tilt"), 2.3);
 *
 * // dims and tilt are macros defined in container.h
 * // they take care of types conversion
 * int rot = dims( params );
 * double tilt = tilt( params );	
 *
 * std::cout << "ROT: " << rot << " TILT: " << tilt << std::endl;
 *
 * @endcode
 */
class xmpContainer
{
    /** Container for pairs "name" and value. Note that void * allows to use
    heterogeneous types */
    std::map<std::string, void *> values;

    /** Adds a new pair to the container */
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

    /** Returns a void * with the requested value, if the label does
        not exist, returns NULL */
    void * getValue( std::string name );
    
    /** Says whether a label exists or not */
    bool valueExists( std::string name );
    
    /** Deletes a value whose label equals "name" */
    void deleteValue( std::string name );
};

//@}
#endif
