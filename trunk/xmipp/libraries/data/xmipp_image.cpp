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

#include "xmipp_image.h"


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


