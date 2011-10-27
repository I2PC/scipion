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
    case UChar:
    case SChar:
        size = sizeof(char);
        break;
    case UShort:
    case Short:
        size = sizeof(short);
        break;
    case UInt:
    case Int:
        size = sizeof(int);
        break;
    case Float:
        size = sizeof(float);
        break;
    case Double:
        size = sizeof(double);
        break;
    case ComplexShort:
        size = sizeof(std::complex<short>);
        break;
    case ComplexInt:
        size = sizeof(std::complex<int>);
        break;
    case ComplexFloat:
        size = sizeof(std::complex<float>);
        break;
    case ComplexDouble:
        size = sizeof(std::complex<double>);
        break;
    case Bool:
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
        datatype = UChar;
    else if (str=="int8")
        datatype = SChar;
    else if (str=="uint16")
        datatype = UShort;
    else if (str=="int16")
        datatype = Short;
    else if (str=="uint32")
        datatype = UInt;
    else if (str=="int32")
        datatype = Int;
    else if (str=="long")
        datatype = Long;
    else if (str=="float")
        datatype = Float;
    else if (str=="double")
        datatype = Double;
    else if (str=="cint16")
        datatype = ComplexShort;
    else if (str=="cint32")
        datatype = ComplexInt;
    else if (str=="cfloat")
        datatype = ComplexFloat;
    else if (str=="cdouble")
        datatype = ComplexDouble;
    else if (str=="bool")
        datatype = Bool;
    else
        REPORT_ERROR(ERR_TYPE_INCORRECT, "datatypeString2int; unknown datatype");

    return datatype;
}

/** Convert datatype to string */
std::string datatype2Str(DataType datatype)
{
    switch ( datatype )
    {
    case UChar:
        return "uint8";
    case SChar:
        return "int8";
    case UShort:
        return "uint16";
    case Short:
        return "int16";
    case UInt:
        return "uint32";
    case Int:
        return "int32";
    case Float:
        return "float";
    case Double:
        return "double";
    case ComplexShort:
        return "cint16";
    case ComplexInt:
        return "cint32";
    case ComplexFloat:
        return "cfloat";
    case ComplexDouble:
        return "cdouble";
    case Bool:
        return "bool";
    default:
        return "unknown type";
    }
}
