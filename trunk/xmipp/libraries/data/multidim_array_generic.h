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

#ifndef MULTIDIM_ARRAY_GENERIC_H_
#define MULTIDIM_ARRAY_GENERIC_H_


#include"multidim_array.h"
#include "image.h"


/// @addtogroup MultidimensionalArrays

//@{

/**
 * MultidimArrayGeneric class to handle arrays with independence of the data type
 */
class MultidimArrayGeneric
{
public:
    DataType       datatype;
    MultidimArrayBase *im;

public:

    /**
     * Constructor with pointer to array to be linked and datatype definition of
     * the linked array.
     */
    MultidimArrayGeneric(MultidimArrayBase* array, DataType _datatype);

    /**
     * Destructor.
     */
    ~MultidimArrayGeneric()
    {}

    /**
     * Link the internal array base to a specific multidimarray object.
     */
    void link(MultidimArrayBase* array);

    /**
     *  Call the resize function of the linked array.
     */
    void resize(unsigned long int Ndim, int Zdim, int Ydim, int Xdim, bool copy=true)
    {
        im->resize(Ndim,Zdim,Ydim,Xdim,copy);
    }

    /**
     *  Copy a specific slice of the linked array.
     */
    template <typename T>
    void getSliceT(int k, MultidimArray<T> &M, char axis = 'Z', unsigned long n = 0) const
    {
        switch(datatype)
        {
        case Float:
            ((MultidimArray<float>*) im)->getSlice(k, M, axis, n);
            break;
        case UInt:
            ((MultidimArray<unsigned int>*) im)->getSlice(k, M, axis, n);
            break;
        case Int:
            ((MultidimArray<int>*) im)->getSlice(k, M, axis, n);
            break;
        case Short:
            ((MultidimArray<short>*) im)->getSlice(k, M, axis, n);
            break;
        case UShort:
            ((MultidimArray<unsigned short>*) im)->getSlice(k, M, axis, n);
            break;
        case SChar:
            ((MultidimArray<char>*) im)->getSlice(k, M, axis, n);
            break;
        case UChar:
            ((MultidimArray<unsigned char>*) im)->getSlice(k, M, axis, n);
            break;
        }
    }

    /**
     *  Copy a specific slice of the linked array.
     */
    void getSlice(int k, MultidimArrayGeneric* M, char axis = 'Z', unsigned long n = 0) const
    {
        switch(M->datatype)
        {
        case Float:
            getSliceT(k, *(MultidimArray<float>*)M->im, axis, n);
            break;
        case UInt:
            getSliceT(k, *(MultidimArray<unsigned int>*)M->im, axis, n);
            break;
        case Int:
            getSliceT(k, *(MultidimArray<int>*)M->im, axis, n);
            break;
        case Short:
            getSliceT(k, *(MultidimArray<short>*)M->im, axis, n);
            break;
        case UShort:
            getSliceT(k, *(MultidimArray<unsigned short>*)M->im, axis, n);
            break;
        case SChar:
            getSliceT(k, *(MultidimArray<char>*)M->im, axis, n);
            break;
        case UChar:
            getSliceT(k, *(MultidimArray<unsigned char>*)M->im, axis, n);
            break;
        }
    }

    /**
     *  Set a specific slice of the linked array.
     */
    template <typename T1>
    void setSlice(int k, const MultidimArray <T1>& v, unsigned long n = 0)
    {
        switch(datatype)
        {
        case Float:
            ((MultidimArray<float>*) im)->setSlice(k, v, n);
            break;
        case UInt:
            ((MultidimArray<unsigned int>*) im)->setSlice(k, v, n);
            break;
        case Int:
            ((MultidimArray<int>*) im)->setSlice(k, v, n);
            break;
        case Short:
            ((MultidimArray<short>*) im)->setSlice(k, v, n);
            break;
        case UShort:
            ((MultidimArray<unsigned short>*) im)->setSlice(k, v, n);
            break;
        case SChar:
            ((MultidimArray<char>*) im)->setSlice(k, v, n);
            break;
        case UChar:
            ((MultidimArray<unsigned char>*) im)->setSlice(k, v, n);
            break;
        }
    }

    /**
     * Set a specific slice of the linked array.
     */
    void setSlice(int k, const MultidimArrayGeneric* v, unsigned long n = 0)
    {
        switch(v->datatype)
        {
        case Float:
            setSlice(k,*(MultidimArray<float>*) v->im, n);
            break;
        case UInt:
            setSlice(k,*(MultidimArray<unsigned int>*) v->im, n);
            break;
        case Int:
            setSlice(k,*(MultidimArray<int>*) v->im, n);
            break;
        case Short:
            setSlice(k,*(MultidimArray<short>*) v->im, n);
            break;
        case UShort:
            setSlice(k,*(MultidimArray<unsigned short>*) v->im, n);
            break;
        case SChar:
            setSlice(k,*(MultidimArray<char>*) v->im, n);
            break;
        case UChar:
            setSlice(k,*(MultidimArray<unsigned char>*) v->im, n);
            break;
        }
    }

    /**
     * Get the dimensions of the linked array.
     */
    void getDimensions(int& Xdim, int& Ydim, int& Zdim, unsigned long int &Ndim) const
    {
        im->getDimensions(Xdim,Ydim,Zdim,Ndim);
    }
    void getDimensions(int& Xdim, int& Ydim, int& Zdim) const
    {
        unsigned long
        int Ndim;
        im->getDimensions(Xdim,Ydim,Zdim,Ndim);
    }

    /**
      * Get the dimensions of the linked array.
      */
    void setXmippOrigin()
    {
        im->setXmippOrigin();
    }
}
;
//@}

#endif /* MULTIDIM_ARRAY_GENERIC_H_ */
