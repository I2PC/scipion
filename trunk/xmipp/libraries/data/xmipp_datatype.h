/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#ifndef DATATYPE_H_
#define DATATYPE_H_

#include <string>
/** @defgroup Datatypes Datatypes for MultidimArrays
 *  @ingroup DataLibrary
*/
//@{
/** Data type.
 * This class defines the datatype of the data inside this image.
 */
typedef enum
{
    Default = -1,           // For writing purposes
    Unknown_Type = 0,       // Undefined data type
    UChar = 1,              // Unsigned character or byte type
    SChar = 2,              // Signed character (for CCP4)
    UShort = 3,             // Unsigned integer (2-byte)
    Short = 4,              // Signed integer (2-byte)
    UInt = 5,               // Unsigned integer (4-byte)
    Int = 6,                // Signed integer (4-byte)
    Long = 7,               // Signed integer (4 or 8 byte, depending on system)
    Float = 8,              // Floating point (4-byte)
    Double = 9,             // Double precision floating point (8-byte)
    ComplexShort = 10,      // Complex two-byte integer (4-byte)
    ComplexInt = 11,        // Complex integer (8-byte)
    ComplexFloat = 12,      // Complex floating point (8-byte)
    ComplexDouble = 13,     // Complex floating point (16-byte)
    Bool = 14,              // Boolean (1-byte?)
    LastEntry = 15          // This must be the last entry
} DataType;


/// Returns memory size of datatype
size_t gettypesize(DataType type);

/** Convert datatype string to datatype enum */
DataType str2Datatype(const std::string & str);

/** Convert datatype to string */
std::string datatype2Str(DataType datatype);

//@}
#endif /* DATATYPE_H_ */
