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

#include "datatype.h"
#include "multidim_array.h"

/* Switch among different datatypes.
 *
 * This macro replicates the code for the different data type options.
 *
 *@code
 *
 *#define MYFUNC(type)  getSlice(k, *(MultidimArray<type>*)image, axis, n)
 *
 *SWITCHDATATYPE(datatype,MYFUNC)
 *
 *@endcode
 */
#define SWITCHDATATYPE(datatype,OP) \
    switch (datatype)\
        {\
        case Float:\
            OP(float);\
            break;\
        case UInt:\
            OP(unsigned int);\
            break;\
        case Int:\
            OP(int);\
            break;\
        case Short:\
            OP(short);\
            break;\
        case UShort:\
            OP(unsigned short);\
            break;\
        case SChar:\
            OP(char);\
            break;\
        case UChar:\
            OP(unsigned char);\
            break;\
        }

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
    {
    }

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
    void getSlice(int k, MultidimArray<T> &M, char axis = 'Z', unsigned long n = 0) const
    {
#define GETSLICE(type) ((MultidimArray<type>*) im)->getSlice(k, M, axis, n)

      SWITCHDATATYPE(datatype,GETSLICE);

#undef GETSLICE
    }

    /**
     *  Copy a specific slice of the linked array.
     */
    void getSlice(int k, MultidimArrayGeneric* M, char axis = 'Z', unsigned long n = 0) const
    {
#define GETSLICE(type) getSlice(k, *(MultidimArray<type>*)M->im, axis, n)

      SWITCHDATATYPE(M->datatype,GETSLICE);

#undef GETSLICE
    }

    /**
     *  Set a specific slice of the linked array.
     */
    template <typename T1>
    void setSlice(int k, const MultidimArray <T1>& v, unsigned long n = 0)
    {
#define SETSLICE(type) ((MultidimArray<type>*) im)->setSlice(k, v, n)

      SWITCHDATATYPE(datatype,SETSLICE);

#undef SETSLICE
    }

    /**
     * Set a specific slice of the linked array.
     */
    void setSlice(int k, const MultidimArrayGeneric* v, unsigned long n = 0)
    {
#define SETSLICE(type) setSlice(k,*(MultidimArray<type>*) v->im, n);

      SWITCHDATATYPE(v->datatype,SETSLICE);

#undef SETSLICE
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
 inline   void setXmippOrigin()
    {
        im->setXmippOrigin();
    }
}
;
//@}

#endif /* MULTIDIM_ARRAY_GENERIC_H_ */
