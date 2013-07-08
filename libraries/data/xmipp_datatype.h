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
    DT_Default = -1,           // For writing purposes
    DT_Unknown = 0,       // Undefined data type
    DT_UChar = 1,              // Unsigned character or byte type
    DT_SChar = 2,              // Signed character (for CCP4)
    DT_UShort = 3,             // Unsigned integer (2-byte)
    DT_Short = 4,              // Signed integer (2-byte)
    DT_UInt = 5,               // Unsigned integer (4-byte)
    DT_Int = 6,                // Signed integer (4-byte)
    DT_ULong = 7,               // Unsigned integer (4 or 8 byte, depending on system)
    DT_Long = 8,               //  Signed integer (4 or 8 byte, depending on system)
    DT_Float = 9,              // Floating point (4-byte)
    DT_Double = 10,             // DT_Double precision floating point (8-byte)
    DT_CShort = 11,      // Complex two-byte integer (4-byte)
    DT_CInt = 12,        // Complex integer (8-byte)
    DT_CFloat = 13,      // Complex floating point (8-byte)
    DT_CDouble = 14,     // Complex floating point (16-byte)
    DT_Bool = 15,              // Boolean (1-byte?)
    DT_LastEntry = 16          // This must be the last entry
} DataType;


/// Returns memory size of datatype
size_t gettypesize(DataType type);

/** Convert datatype string to datatype enum */
DataType str2Datatype(const std::string & str);

/** Convert datatype to string */
std::string datatype2Str(DataType datatype);

/** Convert datatype to string in long format */
std::string datatype2StrLong(DataType datatype);

//@}
#endif /* DATATYPE_H_ */
