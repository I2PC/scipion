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

#include "xmipp_datatype.h"
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
     case DT_Double:\
         {OP(double)};\
            break;\
        case DT_Float:\
            {OP(float)};\
            break;\
        case DT_UInt:\
            {OP(unsigned int)};\
            break;\
        case DT_Int:\
            {OP(int)};\
            break;\
        case DT_Short:\
            {OP(short)};\
            break;\
        case DT_UShort:\
            {OP(unsigned short)};\
            break;\
        case DT_SChar:\
            {OP(char)};\
            break;\
        case DT_UChar:\
            {OP(unsigned char)};\
            break;\
        default:\
			REPORT_ERROR(ERR_ARG_INCORRECT,"Do not know how to handle this type at this point");\
        }

/// @addtogroup MultidimensionalArrays

//@{

/** @name MultidimArrayGenericSpeedUp Speed up macros */
/** Array access.
 *
 * This macros gives you access to the array (T **)
 */
//@{
#ifndef MULTIDIM_ARRAY_BASE
#define MULTIDIM_ARRAY_BASE(v) (*((v).data->im))
#endif

#ifndef MULTIDIM_ARRAY_GENERIC
#define MULTIDIM_ARRAY_GENERIC(v) (*((v).data))
#endif

#ifndef MULTIDIM_ARRAY_TYPE
#define MULTIDIM_ARRAY_TYPE(v, type) (*((MultidimArray<type>*)(v.im)))
#endif

/**
 * MultidimArrayGeneric class to handle arrays with independence of the data type
 */
class MultidimArrayGeneric
{
public:
    DataType       datatype;
    MultidimArrayBase *im;

protected:
    // Flag to allow destroy MultidimArrayBase im data.
    bool destroyData;

public:

    /* Empty constructor */
    MultidimArrayGeneric()
    {
        init();
    }

    /**
     * Constructor with pointer to array to be linked and datatype definition of
     * the linked array.
     */
    MultidimArrayGeneric(MultidimArrayBase* array, DataType _datatype);

    /* Get an aliasSlice of the selected slice from the multidimarray
     */
    MultidimArrayGeneric(MultidimArrayGeneric &mdim, int select_slice);

    /**
     * Destructor.
     */
    ~MultidimArrayGeneric();

    /* Initialize
     */
    void init();

    /* Clear the MultidimArrayBase and others
     */
    void clear();

    /* Set the datatype of the multidimarray object
     */
    void setDatatype(DataType imgType);

    /* Return mdim pointing to a specific slice inside this mda generic
     */
    void aliasSlice(MultidimArrayGeneric &mdim, int select_slice);

    /**
     * Link the internal array base to a specific multidimarray object.
     */
    void link(MultidimArrayBase* array);

    /**
     *  Call the resize function of the linked array.
     */
    void resize(size_t Ndim, int Zdim, int Ydim, int Xdim, bool copy=true)
    {
        im->resize(Ndim,Zdim,Ydim,Xdim,copy);
    }

    void resize(ArrayDim &adim, bool copy=true)
    {
        im->resize(adim, copy);
    }

    void resize(MultidimArrayGeneric &mdim, bool copy=true)
    {
        ArrayDim adim;
        mdim.getDimensions(adim);
        im->resize(adim, copy);
    }

    /** Reverse matrix values over X axis, keep in this object. */
    void selfReverseX()
    {
        im->selfReverseX();
    }
    /** Reverse matrix values over Y axis, keep in this object. */
    void selfReverseY()
    {
        im->selfReverseY();
    }
    /** Reverse matrix values over Z axis, keep in this object. */
    void selfReverseZ()
    {
        im->selfReverseZ();
    }

    /**
     *  Return a pointer to internal multidimarray casted to template T.
     */
    template <typename T>
    void getArrayPointer(T* &M) const
    {
#define GETMULTIDIMARRAY(type) M = (T*) (((MultidimArray<type>*) im)->data);
        SWITCHDATATYPE(datatype,GETMULTIDIMARRAY)
#undef GETMULTIDIMARRAY

    }

    void * getArrayPointer() const
    {
#define GETMULTIDIMARRAY(type) return (void *) (((MultidimArray<type>*) im)->data);
        SWITCHDATATYPE(datatype,GETMULTIDIMARRAY)
#undef GETMULTIDIMARRAY

    }

    /**
     *  Return a pointer to internal multidimarray casted to template T.
     */
    template <typename T>
    void getMultidimArrayPointer(MultidimArray<T>* &M) const
    {
#define GETMULTIDIMARRAY(type) M = (MultidimArray<T>*) ((MultidimArray<type>*) im);
        SWITCHDATATYPE(datatype,GETMULTIDIMARRAY)
#undef GETMULTIDIMARRAY

    }

    /** Get a Window from the image*/
    void window(MultidimArrayGeneric &result, int z0, int y0, int x0,
                int zF, int yF, int xF,
                double init_value = 0.) const
    {
#define WINDOW(type) _window(*((MultidimArray<type>*)(result.im)), z0,y0,x0,zF,yF,xF,init_value);
        SWITCHDATATYPE(result.datatype,WINDOW)
#undef WINDOW

    }

    template <class T>
    void _window(MultidimArray<T> &result, int z0, int y0, int x0,
                 int zF, int yF, int xF,
                 double init_value = 0.) const
    {
#define WINDOW(type) ((MultidimArray<type>*)im)->window(result, z0,y0,x0,zF,yF,xF,(T)init_value);
        SWITCHDATATYPE(datatype,WINDOW)
#undef WINDOW

    }

    void selfWindow(int z0, int y0, int x0,
                    int zF, int yF, int xF,
                    double init_value = 0.)
    {
        if (im->mmapOn)
            REPORT_ERROR(ERR_MMAP, "Cannot resize the image when it is mapped to file.");

#define WINDOW(type) ((MultidimArray<type>*)im)->selfWindow(z0,y0,x0,zF,yF,xF,(type)init_value);

        SWITCHDATATYPE(datatype,WINDOW)
#undef WINDOW

    }

    void patch(MultidimArrayGeneric &patchArray, int x, int y)
    {
#define PATCH(type) ((MultidimArray<type>*)im)->patch(*((MultidimArray<type>*)(patchArray.im)), x, y);
        SWITCHDATATYPE(datatype, PATCH)
#undef PATCH

    }

    /**
     *  Copy a specific slice of the linked array.
     */
    template <typename T>
    void getSlice(int k, MultidimArray<T> &M, char axis = 'Z', bool reverse = false,  size_t n = 0) const
    {
#define GETSLICE(type) ((MultidimArray<type>*) im)->getSlice(k, M, axis, reverse, n);
        SWITCHDATATYPE(datatype,GETSLICE)
#undef GETSLICE

    }

    /**
     *  Copy a specific slice of the linked array.
     */
    void getSlice(int k, MultidimArrayGeneric* M, char axis = 'Z', bool reverse = false, size_t n = 0) const
    {
#define GETSLICE(type) getSlice(k, *(MultidimArray<type>*)M->im, axis, reverse, n);
        SWITCHDATATYPE(M->datatype,GETSLICE)
#undef GETSLICE

    }

    /**
     *  Set a specific slice of the linked array.
     */
    template <typename T1>
    void setSlice(int k, const MultidimArray <T1>& v, size_t n = 0)
    {
#define SETSLICE(type) ((MultidimArray<type>*) im)->setSlice(k, v, n);

        SWITCHDATATYPE(datatype,SETSLICE)

#undef SETSLICE

    }

    /**
     * Set a specific slice of the linked array.
     */
    void setSlice(int k, const MultidimArrayGeneric* v, size_t n = 0)
    {
#define SETSLICE(type) setSlice(k,*(MultidimArray<type>*) v->im, n);
        SWITCHDATATYPE(v->datatype,SETSLICE);
#undef SETSLICE

    }

    /**
     * Get the dimensions of the linked array.
     */
    void getDimensions(size_t& Xdim, size_t& Ydim, size_t& Zdim, size_t &Ndim) const
    {
        im->getDimensions(Xdim,Ydim,Zdim,Ndim);
    }

    void getDimensions(size_t& Xdim, size_t& Ydim, size_t& Zdim) const
    {
        size_t Ndim;
        im->getDimensions(Xdim,Ydim,Zdim,Ndim);
    }

    void getDimensions(size_t& Xdim, size_t& Ydim) const
    {
        size_t Zdim, Ndim;
        im->getDimensions(Xdim,Ydim,Zdim,Ndim);
    }

    void getDimensions(ArrayDim &adim)
    {
        im->getDimensions(adim);
    }

    /**
      * Set Xmipp origin.
      */
    inline void setXmippOrigin()
    {
        im->setXmippOrigin();
    }

    void maxIndex(ArrayCoord &pos)
    {
        im->maxIndex(pos);
    }

    /** Compute average */
    double computeAvg() const
    {
        return im->computeAvg();
    }

    /** Compute statistics.
         *
         * The average, standard deviation, minimum and maximum value are
         * returned.
         */
    void computeStats(double& avg, double& stddev, double& minval, double& maxval) const
    {
#define COMPUTESTATS(type) type Tminval; \
                           type Tmaxval; \
                           ((MultidimArray<type>*)(im))->computeStats(avg, stddev, Tminval, Tmaxval);\
                           minval = Tminval;\
                           maxval = Tmaxval;

        SWITCHDATATYPE(datatype,  COMPUTESTATS)
#undef COMPUTESTATS

    }

    /** Compute minimum and maximum as double values.
         */
    void computeDoubleMinMax(double& minval, double& maxval) const
    {
#define COMPUTESDOUBLEMINMAX(type) ((MultidimArray<type>*)(im))->computeDoubleMinMax(minval, maxval);
        SWITCHDATATYPE(datatype,  COMPUTESDOUBLEMINMAX)
#undef COMPUTESDOUBLEMINMAX

    }

    /**
     *  Range adjust using an example
     */
    void rangeAdjust(const MultidimArrayGeneric &example, const MultidimArray<int> *mask=NULL)
    {
#define RANGEADJUST(type) ((MultidimArray<type>*)(im))->rangeAdjust(*(MultidimArray<type>*)(example.im),mask);
        SWITCHDATATYPE(datatype,RANGEADJUST)
#undef RANGEADJUST
    }

    /** Assignment **/
    MultidimArrayGeneric& operator=(const MultidimArrayGeneric& input)
    {
        if (&input != this && input.datatype != DT_Unknown)
        {
            setDatatype(input.datatype);
            *im = *input.im;
        }
        return *this;
    }

    /** *= **/
    void operator*=(double op1)
    {
#define OPERATORTIMESEQUAL(type) *((MultidimArray<type>*)(im)) *= (type)op1;
        SWITCHDATATYPE(datatype,OPERATORTIMESEQUAL)
#undef OPERATORTIMESEQUAL
    }

    /** += **/
    void operator+=(const MultidimArrayGeneric& op1)
    {
#define OPERATORPLUSEQUAL(type) *((MultidimArray<type>*)(im)) += *((MultidimArray<type>*)(op1.im));
        SWITCHDATATYPE(datatype,OPERATORPLUSEQUAL)
#undef OPERATORPLUSEQUAL
    }

    /** Get constant access */
    double operator()(size_t n, int k, int i, int j) const
    {
#define GETVALUE(type) return NZYX_ELEM(*(MultidimArray<type>*)im,n,k,i,j);
        SWITCHDATATYPE(datatype,GETVALUE)
#undef GETVALUE

    }
    /**
     * equal operator
     */
    bool operator==(const MultidimArrayGeneric &mdA) const;

    /** Equality.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const MultidimArrayGeneric &op,
               double accuracy = XMIPP_EQUAL_ACCURACY) const;

    /** Get constant access */
    double operator()(int i, int j) const
    {
#define GETVALUE(type) return A2D_ELEM(*(MultidimArray<type>*)im,i,j);
        SWITCHDATATYPE(datatype,GETVALUE)
#undef GETVALUE

    }

    /** Get array */
    MultidimArrayBase & operator()()
    {
        return *im;
    }

    /** Get array */
    const MultidimArrayBase & operator()() const
    {
        return *im;
    }


    /** Copy the image in MultidimarrayGeneric to a specific T MultidimArray
     */
    template <typename T>
    void getImage(MultidimArray<T> &M) const
    {
#define TYPECAST(type) typeCast(*(MultidimArray<type>*)(im), M);
        SWITCHDATATYPE(datatype, TYPECAST)
#undef TYPECAST

    }

    /** Copy in MultidimarrayGeneric an image from a specific T MultidimArray
     */
    template <typename T>
    void setImage(MultidimArray<T> &M)
    {
#define TYPECAST(type) typeCast(M, *(MultidimArray<type>*)(im));
        SWITCHDATATYPE(datatype, TYPECAST)
#undef TYPECAST

    }

    /** Get window.
     */
    template <typename T1>
    void window(MultidimArray<T1> &result, int n0,int z0, int y0, int x0,
                int nF,int zF, int yF, int xF,
                T1 init_value = 0) const
    {
#define WINDOW(type) ((MultidimArray<type>*)(im))->window(result, n0, z0, y0, x0, nF, zF, yF, xF, init_value);
        SWITCHDATATYPE(datatype, WINDOW)
#undef WINDOW

    }
}
;
//@}
//@}
#endif /* MULTIDIM_ARRAY_GENERIC_H_ */
