/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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

#include "image.h"

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
DataType datatypeString2Int(std::string str)
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


// Special cases for complex numbers
template<>
void Image< std::complex< double > >::castPage2T(char * page,
        std::complex<double> * ptrDest,
        DataType datatype,
        size_t pageSize)
{

    switch (datatype)
    {
    case ComplexShort:
        {
            std::complex<short> * ptr = (std::complex<short> *) page;
            for(int i=0; i<pageSize;i++)
                ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
            break;
        }
    case ComplexInt:
            {
                std::complex<int> * ptr = (std::complex<int> *) page;
                for(int i=0; i<pageSize;i++)
                ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
                break;
            }
        case ComplexFloat:
                {
                    std::complex<float> * ptr = (std::complex<float> *) page;
                    for(int i=0; i<pageSize;i++)
                        ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
                        break;
                    }
                case ComplexDouble:
                        {
                            memcpy(ptrDest, page, pageSize*sizeof(std::complex<double>));
                            break;
                        }
                    default:
                            {
                                std::cerr<<"Datatype= "<<datatype<<std::endl;
                                REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to std::complex<double>");
                                break;
                            }
                        }
}

template<>
void Image< std::complex< double > >::castPage2Datatype(std::complex<double> * srcPtr,
        char * page,
        DataType datatype,
        size_t pageSize)
{
    switch (datatype)
    {
    case ComplexShort:
        {
            std::complex<short> * ptr = (std::complex<short> *) page;
            for(int i=0; i<pageSize;i++)
                ptr[i] = (std::complex<short>)srcPtr[i];
            break;
        }
    case ComplexInt:
            {
                std::complex<int> * ptr = (std::complex<int> *) page;
                for(int i=0; i<pageSize;i++)
                ptr[i] = (std::complex<int>)srcPtr[i];
                break;
            }
        case ComplexFloat:
                {
                    std::complex<float> * ptr = (std::complex<float> *) page;
                    for(int i=0; i<pageSize;i++)
                        ptr[i] = (std::complex<float>)srcPtr[i];
                        break;
                    }
                case ComplexDouble:
                        {
                            memcpy(page, srcPtr, pageSize*sizeof(std::complex<double>));
                            break;
                        }
                    default:
                            {
                                std::cerr<<"Datatype= "<<datatype<<std::endl;
                                REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to complex<double>");
                                break;
                            }
                        }
}

template<>
void Image< std::complex< double > >::castConvertPage2Datatype(std::complex< double > * srcPtr,
        char * page, DataType datatype, size_t pageSize,double min0,double max0,CastWriteMode castMode)
{

    switch (datatype)
    {
    case ComplexFloat:
        {
            std::complex<float> * ptr = (std::complex<float> *) page;
            for(int i=0; i<pageSize;i++)
                ptr[i] = (std::complex<float>)srcPtr[i];
            break;
        }
    default:
            {
                std::cerr<<"Datatype= "<<datatype<<std::endl;
                REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast&convert datatype to complex<double>");
                break;
            }
        }
}


