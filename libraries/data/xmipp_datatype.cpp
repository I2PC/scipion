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

#include <complex>
#include "xmipp_datatype.h"
#include "xmipp_error.h"


// Get size of datatype
size_t gettypesize(DataType type)
{
    size_t   size;

    switch ( type )
    {
    case DT_UChar:
    case DT_SChar:
        size = sizeof(char);
        break;
    case DT_UShort:
    case DT_Short:
        size = sizeof(short);
        break;
    case DT_UInt:
    case DT_Int:
        size = sizeof(int);
        break;
    case DT_Float:
        size = sizeof(float);
        break;
    case DT_Double:
        size = sizeof(double);
        break;
    case DT_CShort:
        size = sizeof(std::complex<short>);
        break;
    case DT_CInt:
        size = sizeof(std::complex<int>);
        break;
    case DT_CFloat:
        size = sizeof(std::complex<float>);
        break;
    case DT_CDouble:
        size = sizeof(std::complex<double>);
        break;
    case DT_Bool:
        size = sizeof(bool);
        break;
    default:
        size = 0;
    }

    return(size);
}

/** Convert datatype string to datatypr enun */
DataType str2Datatype(const std::string & str)
{
    DataType datatype;

    if(str=="uint8")
        datatype = DT_UChar;
    else if (str=="int8")
        datatype = DT_SChar;
    else if (str=="uint16")
        datatype = DT_UShort;
    else if (str=="int16")
        datatype = DT_Short;
    else if (str=="uint32")
        datatype = DT_UInt;
    else if (str=="int32")
        datatype = DT_Int;
    else if (str=="long")
        datatype = DT_Long;
    else if (str=="float")
        datatype = DT_Float;
    else if (str=="double")
        datatype = DT_Double;
    else if (str=="cint16")
        datatype = DT_CShort;
    else if (str=="cint32")
        datatype = DT_CInt;
    else if (str=="cfloat")
        datatype = DT_CFloat;
    else if (str=="cdouble")
        datatype = DT_CDouble;
    else if (str=="bool")
        datatype = DT_Bool;
    else
        REPORT_ERROR(ERR_TYPE_INCORRECT, "datatypeString2int; unknown datatype");

    return datatype;
}

/** Convert datatype to string */
std::string datatype2Str(DataType datatype)
{
    switch ( datatype )
    {
    case DT_UChar:
        return "uint8";
    case DT_SChar:
        return "int8";
    case DT_UShort:
        return "uint16";
    case DT_Short:
        return "int16";
    case DT_UInt:
        return "uint32";
    case DT_Int:
        return "int32";
    case DT_Long:
        return "int64";
    case DT_Float:
        return "float";
    case DT_Double:
        return "double";
    case DT_CShort:
        return "cint16";
    case DT_CInt:
        return "cint32";
    case DT_CFloat:
        return "cfloat";
    case DT_CDouble:
        return "cdouble";
    case DT_Bool:
        return "bool";
    default:
        return "unknown type";
    }
}

std::string datatype2StrLong(DataType datatype)
{
    switch (datatype)
    {
    case DT_UChar:
        return "Unsigned character or byte type (UInt8)";
        break;
    case DT_SChar:
        return "Signed character (Int8)";
        break;
    case DT_UShort:
        return "Unsigned short integer (UInt16)";
        break;
    case DT_Short:
        return "Signed short integer (Int16)";
        break;
    case DT_UInt:
        return "Unsigned integer (UInt32)";
        break;
    case DT_Int:
        return "Signed integer (Int32)";
        break;
    case DT_Long:
        return "Signed integer (4 or 8 byte, depending on system)";
        break;
    case DT_Float:
        return "Floating point (4-byte)";
        break;
    case DT_Double:
        return "Double precision floating point (8-byte)";
        break;
    case DT_CShort:
        return "Complex two-byte integer (4-byte)";
        break;
    case DT_CInt:
        return "Complex integer (8-byte)";
        break;
    case DT_CFloat:
        return "Complex floating point (8-byte)";
        break;
    case DT_CDouble:
        return "Complex floating point (16-byte)";
        break;
    case DT_Bool:
        return "Boolean (1-byte?)";
        break;
    case DT_Unknown:
        return "Unknown data type";
        break;
    default:
        return "Undefined data type";
        break;
    }
}

