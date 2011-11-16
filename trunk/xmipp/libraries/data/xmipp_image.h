/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
 *    Joaquin Oton       (joton@cnb.csic.es)
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

#ifndef IMAGE_H
#define IMAGE_H

#include "xmipp_image_base.h"
#include "xmipp_image_generic.h"
#include "xmipp_color.h"

/// @addtogroup Images
//@{

/** Size of the page used to read and write images from/to file */
const size_t rw_max_page_size = 4194304; // 4Mb

/** Template class for images.
 * The image class is the general image handling class.
 */
template<typename T>
class Image: public ImageBase
{

public:
    MultidimArray<T>    data;        // The image data array

public:
    /** Empty constructor
     *
     * An empty image is created.
     *
     * @code
     * Image<double> I;
     * @endcode
     */
    Image()
    {
        init();
    }

    /** Constructor with size and filename
      *
      * An image file, which name and format are given by filename,
      * is created with the given size. Then the image is mapped to this file.
      *
      * @code
      * Image I(64,64,1,1,"image.spi");
      * @endcode
      */
    Image(int Xdim, int Ydim, int Zdim, int Ndim, const FileName &_filename)
    {
        init();
        mmapOnWrite = true;
        data.setDimensions(Xdim, Ydim, Zdim, Ndim);
        MD.resize(Ndim);
        filename = _filename;
        ImageFHandler *hFile = openFile(_filename, WRITE_OVERWRITE);
        _write(_filename, hFile, ALL_IMAGES, false, WRITE_OVERWRITE);
        closeFile(hFile);
    }

    /** Constructor with size
     *
     * A blank image (0.0 filled) is created with the given size. Pay attention
     * to the dimension order: Y and then X. If _mmapOn is True then image is allocated
     * in a temporary file.
     *
     * @code
     * Image I(64,64);
     * @endcode
     */
    Image(int Xdim, int Ydim, int Zdim=1, int Ndim=1, bool _mmapOn=false)
    {
        init();
        data.setMmap(_mmapOn);
        data.coreAllocate(Ndim, Zdim, Ydim, Xdim);
        MD.resize(Ndim);
    }

    /** Constructor with MultidimArray alias
     *
     *  An image is created directly with its multidimarray aliased to im.
     *  This function is useful when debugging and you want to save multidimarrays.
     */
    Image(const MultidimArray<T> &im)
    {
        init();
        data.alias(im);
    }

    /** Destructor.
        */
    virtual ~Image()
    {
        if (mmapOnRead || mmapOnWrite)
            munmapFile();
        else
            data.clear();
    }

    /** Clear.
     * Initialize everything to 0
     */
    void clear()
    {
        if (mmapOnRead || mmapOnWrite)
            munmapFile();
        else
            data.clear();
        init();
    }

    /** Check whether image is complex based on T
     */
    bool isComplexT() const
    {
        return ( typeid(T) == typeid(std::complex<double>) ||
                 typeid(T) == typeid(std::complex<float>) );
    }


    /** Cast a page of data from type dataType to type Tdest
     *    input pointer  char *
     */
    void castPage2T(char * page, T * ptrDest, DataType datatype, size_t pageSize )
    {
        switch (datatype)
        {
        case Unknown_Type:
            REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR: datatype is Unknown_Type");
        case UChar:
            {
                if (typeid(T) == typeid(unsigned char))
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                else
                {
                    unsigned char * ptr = (unsigned char *) page;
                    for(int i=0; i<pageSize; ++i, ++ptr)
                        ptrDest[i]=(T) *ptr;
                }
                break;
            }
        case SChar:
                {
                    if (typeid(T) == typeid(signed char))
                {
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                    }
                    else
                    {
                        signed char * ptr = (signed char *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            ptrDest[i]=(T) *ptr;
                    }
                break;
            }
        case UShort:
                {
                    if (typeid(T) == typeid(unsigned short))
                {
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                    }
                    else
                    {
                        unsigned short * ptr = (unsigned short *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            ptrDest[i]=(T) *ptr;
                    }
                break;
            }
        case Short:
                {
                    if (typeid(T) == typeid(short))
                {
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                    }
                    else
                    {
                        short * ptr = (short *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            ptrDest[i]=(T) *ptr;
                    }
                break;
            }
        case UInt:
                {
                    if (typeid(T) == typeid(unsigned int))
                {
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                    }
                    else
                    {
                        unsigned int * ptr = (unsigned int *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            ptrDest[i]=(T) *ptr;
                    }
                break;
            }
        case Int:
                {
                    if (typeid(T) == typeid(int))
                {
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                    }
                    else
                    {
                        int * ptr = (int *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            ptrDest[i]=(T) *ptr;
                    }
                break;
            }
        case Long:
                {
                    if (typeid(T) == typeid(long))
                {
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                    }
                    else
                    {
                        long * ptr = (long *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            ptrDest[i]=(T) *ptr;
                    }
                break;
            }
        case Float:
                {
                    if (typeid(T) == typeid(float))
                {
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                    }
                    else
                    {
                        float * ptr = (float *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            ptrDest[i]=(T) *ptr;
                    }
                break;
            }
        case Double:
                {
                    if (typeid(T) == typeid(double))
                {
                    memcpy(ptrDest, page, pageSize*sizeof(T));
                    }
                    else
                    {
                        double * ptr = (double *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            ptrDest[i]=(T) *ptr;
                    }
                break;
            }
        default:
                {
                    std::cerr<<"Datatype= "<<datatype<<std::endl;
                    REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to T");
                    break;
                }
            }

    }

    /** Cast page from T to datatype
     *  input pointer char *
     */
    void castPage2Datatype(T * srcPtr, char * page, DataType datatype, size_t pageSize )
    {
        switch (datatype)
        {
        case Float:
            {
                if (typeid(T) == typeid(float))
                {
                    memcpy(page, srcPtr, pageSize*sizeof(T));
                }
                else
                {
                    float * ptr = (float *) page;
                    for(int i=0; i<pageSize; ++i, ++ptr)
                        *ptr = (float)srcPtr[i];
                }
                break;
            }
        case Double:
                {
                    if (typeid(T) == typeid(double))
                {
                    memcpy(page, srcPtr, pageSize*sizeof(T));
                    }
                    else
                    {
                        double * ptr = (double *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            *ptr = (double)srcPtr[i];
                    }
                break;
            }
        case UShort:
                {
                    if (typeid(T) == typeid(unsigned short))
                {
                    memcpy(page, srcPtr, pageSize*sizeof(T));
                    }
                    else
                    {
                        unsigned short * ptr = (unsigned short *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            *ptr = (unsigned short)srcPtr[i];
                    }
                break;
            }
        case Short:
                {
                    if (typeid(T) == typeid(short))
                {
                    memcpy(page, srcPtr, pageSize*sizeof(T));
                    }
                    else
                    {
                        short * ptr = (short *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            *ptr = (short)srcPtr[i];
                    }
                break;
            }
        case UChar:
                {
                    if (typeid(T) == typeid(unsigned char))
                {
                    memcpy(page, srcPtr, pageSize*sizeof(T));
                    }
                    else
                    {
                        unsigned char * ptr = (unsigned char *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            *ptr = (unsigned char)srcPtr[i];
                    }
                break;
            }
        case SChar:
                {
                    if (typeid(T) == typeid(char))
                {
                    memcpy(page, srcPtr, pageSize*sizeof(T));
                    }
                    else
                    {
                        char * ptr = (char *) page;
                        for(int i=0; i<pageSize; ++i, ++ptr)
                            *ptr = (char)srcPtr[i];
                    }
                break;
            }
        default:
                {
                    std::cerr<<"outputDatatype = " << datatype << std::endl;
                    REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast T to outputDatatype");
                    break;
                }
            }
    }

    /* Convert the pixels values from one datatype to another, taking into account for datatypes
     * of same bitdepth the shift of the minimum values. In other cases, the conversion is done
     * adjusting the input values in the range of output datatype.
     */
    void castConvertPage2Datatype(T * srcPtr, char * page, DataType datatype, size_t pageSize,
                                  double min0, double max0, CastWriteMode castMode=CW_CONVERT)
    {

        double minF, maxF;
        double slope;
        size_t n;
        DataType  myTypeId = myT();

        switch(datatype)
        {
        case UChar:
            {
                if (castMode==CW_CONVERT && myTypeId == SChar)
                {
                    slope = 1;
                    min0 -= CHAR_MIN;
                }
                else
                {
                    minF = 0;
                    maxF = UCHAR_MAX;
                    if (max0 != min0)
                        slope = static_cast< double >(maxF - minF) /
                                static_cast< double >(max0 - min0);
                    else
                        slope = 0;
                }
                unsigned char * ptr = (unsigned char *) page;

                for( n=0; n<pageSize; n++)
                    ptr[n] = minF + static_cast< unsigned char >(slope *
                             static_cast< double >(srcPtr[n] - min0));

                break;
            }
        case SChar:
                {
                    if (castMode==CW_CONVERT &&  myTypeId == UChar)
                {
                    slope = 1;
                    min0 += CHAR_MIN;
                }
                else
                {
                    minF = CHAR_MIN;
                    maxF = CHAR_MAX;
                    if (max0 != min0)
                            slope = static_cast< double >(maxF - minF) /
                                    static_cast< double >(max0 - min0);
                        else
                            slope = 0;
                    }
                char * ptr = (char *) page;

                for( n=0; n<pageSize; n++)
                ptr[n] = minF + static_cast< char >(slope *
                                                    static_cast< double >(srcPtr[n] - min0));

                break;
            }
        case UShort:
                {
                    if (castMode==CW_CONVERT && (myTypeId == SChar|| myTypeId == Short))
                    {
                        slope = 1;
                        min0 -= SHRT_MIN;
                    }
                    else if (castMode==CW_CONVERT && (myTypeId == UChar))
                    {
                        slope = 1;
                    }
                    else
                    {
                        minF = 0;
                        maxF = USHRT_MAX;
                        if (max0 != min0)
                                slope = static_cast<double >(maxF-minF)/static_cast<double >(max0-min0);
                            else
                                slope = 0;
                        }

                unsigned short * ptr = (unsigned short *) page;

                for( n=0; n<pageSize; n++)
                ptr[n] = minF + static_cast< unsigned short >(slope *
                         static_cast< double >(srcPtr[n] - min0));

                break;
            }
        case Short:
                {
                    if (castMode==CW_CONVERT && (myTypeId == UChar || myTypeId == UShort))
                    {
                        slope = 1;
                        min0 += SHRT_MIN;
                    }
                    else if (castMode==CW_CONVERT && (myTypeId == SChar))
                    {
                        slope = 1;
                    }
                    else
                    {
                        minF = SHRT_MIN;
                        maxF = SHRT_MAX;
                        if (max0 != min0)
                                slope = static_cast< double >(maxF - minF) /
                                        static_cast< double >(max0 - min0);
                            else
                                slope = 0;
                        }
                short * ptr = (short *) page;

                for( n=0; n<pageSize; n++)
                ptr[n] = minF + static_cast< short >(slope *
                                                     static_cast< double >(srcPtr[n] - min0));

                break;
            }
        case UInt:
                {
                    if (castMode==CW_CONVERT &&  (myTypeId == SChar
                                                      || myTypeId == Short
                                                      || myTypeId == Int))
                    {
                        slope = 1;
                        min0 -= INT_MIN;
                    }
                    else if (castMode==CW_CONVERT &&  (myTypeId == UShort || myTypeId == UChar))
                    {
                        slope = 1;
                    }
                    else
                    {
                        minF = 0;
                        maxF = UINT_MAX;
                        if (max0 != min0)
                                slope = static_cast< double >(maxF - minF) /
                                        static_cast< double >(max0 - min0);
                            else
                                slope = 0;
                        }
                unsigned int * ptr = (unsigned int *) page;

                for( n=0; n<pageSize; n++)
                ptr[n] = minF + static_cast< unsigned int >(slope *
                         static_cast< double >(srcPtr[n] - min0));
                break;
            }
        case Int:
                {
                    if (castMode==CW_CONVERT &&  (myTypeId == UChar
                                                      || myTypeId == UShort
                                                      || myTypeId == UInt))
                    {
                        slope = 1;
                        min0 += INT_MIN;
                    }
                    else if (castMode==CW_CONVERT &&  (myTypeId == Short || myTypeId == SChar))
                    {
                        slope = 1;
                    }
                    else
                    {
                        minF = INT_MIN;
                        maxF = INT_MAX;
                        if (max0 != min0)
                                slope = static_cast< double >(maxF - minF) /
                                        static_cast< double >(max0 - min0);
                            else
                                slope = 0;
                        }
                int * ptr = (int *) page;

                for( n=0; n<pageSize; n++)
                ptr[n] = minF + static_cast< int >(slope *
                                                   static_cast< double >(srcPtr[n] - min0));
                break;
            }
        default:
                castPage2Datatype(srcPtr,page,datatype,pageSize);
        }

    }
    /** Check if file Datatype is the same as the declared image object (T type) to use mmap.
     */
    bool checkMmapT(DataType datatype)
    {

        switch (datatype)
        {
        case Unknown_Type:
            REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR: datatype is Unknown_Type");
        case UChar:
            {
                if (typeid(T) == typeid(unsigned char))
                    return 1;
                else
                    return 0;
                break;
            }
        case SChar:
            {
                if (typeid(T) == typeid(char))
                    return 1;
                else
                    return 0;
                break;
            }
        case UShort:
            {
                if (typeid(T) == typeid(unsigned short))
                    return 1;
                else
                    return 0;
                break;
            }
        case Short:
            {
                if (typeid(T) == typeid(short))
                    return 1;
                else
                    return 0;
                break;
            }
        case UInt:
            {
                if (typeid(T) == typeid(unsigned int))
                    return 1;
                else
                    return 0;
                break;
            }
        case Int:
            {
                if (typeid(T) == typeid(int))
                    return 1;
                else
                    return 0;
                break;
            }
        case Long:
            {
                if (typeid(T) == typeid(long))
                    return 1;
                else
                    return 0;
                break;
            }
        case Float:
            {
                if (typeid(T) == typeid(float))
                    return 1;
                else
                    return 0;
                break;
            }
        case Double:
            {
                if (typeid(T) == typeid(double))
                    return 1;
                else
                    return 0;
                break;
            }
        default:
            {
                std::cerr<<"Datatype= "<<datatype<<std::endl;
                REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to T");
                break;
            }
        }
        //               int * iTemp = (int*) map;
        //                ptrDest = reinterpret_cast<T*> (iTemp);
    }

    /* Read an image with a lower resolution as a preview image.
     * If Zdim parameter is not passed, then all slices are rescaled.
     * If Ydim is not passed, then Ydim is rescaled same factor as Xdim.
     */
    int readPreview(const FileName &name, int Xdim, int Ydim = -1, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE)
    {
        // Zdim is used to choose the slices: -1 = CENTRAL_SLICE, 0 = ALL_SLICES, else This Slice

        ImageGeneric im;
        int imXdim, imYdim, imZdim, Zdim;
        int err;
        err = im.readMapped(name, select_img);
        im.getDimensions(imXdim, imYdim, imZdim);
        ImageInfo imgInfo;
        im.getInfo(imgInfo);

        //Set information from image file
        setName(name);
        setDatatype(imgInfo.datatype);
        aDimFile = imgInfo.adim;

        im().setXmippOrigin();

        double scale;

        // If only Xdim is passed, it is the higher allowable size, for any dimension
        if (Ydim == -1 && imXdim < imYdim)
        {
            Ydim = Xdim;
            scale = ((double) Ydim)/((double) imYdim);
            Xdim = imXdim * scale;
        }
        else
        {
            scale = ((double) Xdim)/((double) imXdim);
            if (Ydim == -1)
                Ydim = imYdim * scale;
        }

        int mode = (scale <= 1)? NEAREST : LINEAR; // If scale factor is higher than 1, LINEAR mode is used to avoid artifacts

        if (select_slice > ALL_SLICES) // In this case a specific slice number has been chosen (Not central slice)
        {
            MultidimArrayGeneric array(im(), select_slice - 1);
            array.setXmippOrigin();

            scaleToSize(mode, IMGMATRIX(*this), array ,Xdim, Ydim);
        }
        else // Otherwise, All slices or Central slice is selected
        {
            Zdim = (select_slice == ALL_SLICES)? imZdim: 1;
            scaleToSize(mode, IMGMATRIX(*this), im(), Xdim, Ydim, Zdim);
        }

        IMGMATRIX(*this).resetOrigin();
        return err;
    }

    /** It changes the behavior of the internal multidimarray so it points to a specific slice of
     *  the initial volume. No information is deallocated from memory, so it is also possible to
     *  repoint to the whole volume (passing select_slice = ALL_SLICES), or CENTRAL_SLICE.
     */
    void movePointerToSlice(int select_slice = ALL_SLICES)
    {
        if (select_slice > aDimFile.zdim)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, formatString("movePointerToSlice: Selected slice %4d cannot be higher than Z size %4d.",
                         select_slice,aDimFile.zdim));

        ArrayDim newDim = aDimFile;
        int phys_slice;

        switch (select_slice)
        {
        case CENTRAL_SLICE:
            phys_slice = aDimFile.zdim/2;
            newDim.zdim = 1;
            break;
        case ALL_SLICES:
            phys_slice = 0;
            break;
        default:
            phys_slice = select_slice - 1;
            newDim.zdim = 1;
            break;
        }

        VOLMATRIX(*this).setDimensions(newDim);
        MULTIDIM_ARRAY(VOLMATRIX(*this)) += YXSIZE(VOLMATRIX(*this)) * (phys_slice - SLICE_INDEX(mappedSlice));
        mappedSlice = (select_slice == CENTRAL_SLICE)? phys_slice + 1 : select_slice;
    }

    /** Write an entire page as datatype
     *
     * A page of datasize_n elements T is cast to datatype and written to fimg
     * The memory for the casted page is allocated and freed internally.
     */
    void writePageAsDatatype(FILE * fimg, DataType datatype, size_t datasize_n )
    {
        size_t datasize = datasize_n * gettypesize(datatype);
        char * fdata = (char *) askMemory(datasize);
        castPage2Datatype(MULTIDIM_ARRAY(data), fdata, datatype, datasize_n);
        fwrite( fdata, datasize, 1, fimg );
        freeMemory(fdata, datasize);
    }

    /** Data access
     *
     * This operator can be used to access the data multidimarray.
     * In this way we could resize an image just by
     * resizing its associated matrix or we could add two images by adding their
     * matrices.
     * @code
     * I().resize(128, 128);
     * I2() = I1() + I2();
     * @endcode
     */
    MultidimArray<T>& operator()()
    {
        return data;
    }
    const MultidimArray<T>& operator()() const
    {
        return data;
    }

    /** Pixel access
    *
    * This operator is used to access a pixel within a 2D image. This is a
    * logical access, so you could access to negative positions if the image
    * has been defined so (see the general explanation for the class).
    *
    * @code
    * std::cout << "Grey level of pixel (-3,-3) of the image = " << I(-3, -3)
    * << std::endl;
    *
    * I(-3, -3) = I(-3, -2);
    * @endcode
    */
    T& operator()(int i, int j) const
    {
        return A2D_ELEM(data, i, j);
    }
    /**
     * equal operator
     */
    bool operator==(const Image<T> &i1) const
    {
        return(this->data == i1.data);
    }

    /** Voxel access
     *
     * This operator is used to access a voxel within a 3D image. This is a
     * logical access, so you could access to negative positions if the image
     * has been defined so (see the general explanation for the class).
     *
     * @code
     * std::cout << "Grey level of pixel (-3,-3, 1) of the volume = " << I(-3, -3, 1)
     * << std::endl;
     *
     * I(-3, -3, 1) = I(-3, -2, 0);
     * @endcode
     */
    T& operator()(int k, int i, int j) const
    {
        return A3D_ELEM(data, k, i, j);
    }

    /** Get Image dimensions
     */
    void getDimensions(int &Xdim, int &Ydim, int &Zdim, size_t &Ndim) const
    {
        Xdim = XSIZE(data);
        Ydim = YSIZE(data);
        Zdim = ZSIZE(data);
        Ndim = NSIZE(data);
    }

    size_t getSize() const
    {
        return NZYXSIZE(data);
    }

    /** Get geometric transformation matrix from 2D-image header
      */
    void getTransformationMatrix(Matrix2D<double> &A,
                                 bool only_apply_shifts = false,
                                 const size_t n = 0)
    {
        // This has only been implemented for 2D images...
        MULTIDIM_ARRAY(*this).checkDimension(2);
        A.resizeNoCopy(3,3);
        geo2TransformationMatrix(MD[n], A, only_apply_shifts);
    }

    /** Sum this object with other file and keep in this object
      */
    void sumWithFile(const FileName &fn)
    {
        Image<T> aux;
        aux.read(fn, DATA, -1, true);
        (*this)()+=aux();
    }

    /**
     *  Specific read functions for different file formats
     */
#include "rwTIFF.h"

protected:

    /** Apply geometry in refering metadata to the image */
    void applyGeo(const MDRow &row, bool only_apply_shifts = false, bool wrap = WRAP)
    {
        //This implementation does not handle stacks,
        //read in a block
        if (data.ndim != 1)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "Hgeometric transformation cannot be applied to stacks!!!");
        if (MD.size()==0)
            MD.push_back(MDL::emptyHeader);
        MDRow &rowAux = MD[0];

        double aux;
        //origins
        if (row.getValue(MDL_ORIGINX, aux))
            rowAux.setValue(MDL_ORIGINX, aux);
        if (row.getValue(MDL_ORIGINY, aux))
            rowAux.setValue(MDL_ORIGINY, aux);
        if (row.getValue(MDL_ORIGINZ, aux))
            rowAux.setValue(MDL_ORIGINZ, aux);
        //shifts
        if (row.getValue(MDL_SHIFTX, aux))
            rowAux.setValue(MDL_SHIFTX, aux);
        if (row.getValue(MDL_SHIFTY, aux))
            rowAux.setValue(MDL_SHIFTY, aux);
        if (row.getValue(MDL_SHIFTZ, aux))
            rowAux.setValue(MDL_SHIFTZ, aux);
        //rotations
        if (row.getValue(MDL_ANGLEROT, aux))
            rowAux.setValue(MDL_ANGLEROT, aux);
        if (row.getValue(MDL_ANGLETILT, aux))
            rowAux.setValue(MDL_ANGLETILT, aux);
        if (row.getValue(MDL_ANGLEPSI, aux))
            rowAux.setValue(MDL_ANGLEPSI, aux);
        //scale
        if (row.getValue(MDL_SCALE, aux))
            rowAux.setValue(MDL_SCALE, aux);
        //weight
        if (row.getValue(MDL_WEIGHT, aux))
            rowAux.setValue(MDL_WEIGHT, aux);
        bool auxBool;
        if (row.getValue(MDL_FLIP, auxBool))
            rowAux.setValue(MDL_FLIP, auxBool);

        //apply geo has not been defined for volumes
        //and only make sense when reading data
        if (data.getDim() < 3 && dataMode >= DATA)
        {
            Matrix2D< double > A;
            getTransformationMatrix(A, only_apply_shifts);
            if (!A.isIdentity())
            {
                MultidimArray<T> tmp=MULTIDIM_ARRAY(*this);
                applyGeometry(BSPLINE3, MULTIDIM_ARRAY(*this), tmp,
                              A, IS_NOT_INV, wrap);
            }
        }
    }

    /** Set the image dimensions
     */
    void setDimensions(int Xdim, int Ydim, int Zdim, size_t Ndim)
    {
        data.setDimensions(Xdim,Ydim,Zdim,Ndim);
        data.getDimensions(aDimFile);
    }

private:

    /** Read the raw data
      */
    void readData(FILE* fimg, size_t select_img, DataType datatype, size_t pad)
    {
        //#define DEBUG
#ifdef DEBUG
        std::cerr<<"entering readdata"<<std::endl;
        std::cerr<<" readData flag= "<<dataMode<<std::endl;
#endif
#undef DEBUG

        if ( dataMode < DATA )
            return;

        // If only half of a transform is stored, it needs to be handled
        if (transform == Hermitian || transform == CentHerm )
            data.setXdim(XSIZE(data)/2 + 1);

        size_t selectImgOffset, readsize, readsize_n, pagemax = 4194304; //4Mb
        size_t datatypesize=gettypesize(datatype);
        size_t pagesize  =ZYXSIZE(data)*datatypesize;
        size_t haveread_n = 0;

        selectImgOffset = offset + IMG_INDEX(select_img) * (pagesize + pad);

        // Flag to know that data is not going to be mapped although mmapOn is true
        if (mmapOnRead && !checkMmapT(datatype))
        {
            reportWarning("Image::readData: File datatype and image declaration not "
                          "compatible with mmap. Loading into memory.");
            mmapOnRead = false;
            mFd = -1;
        }

        if (mmapOnRead)
        {
            // Image mmapOn is not compatible with Multidimarray mmapOn
            if(data.mmapOn)
                REPORT_ERROR(ERR_MULTIDIM_DIM,"Image Class::ReadData: mmap option can not be selected simultaneously\
                             for both Image class and its Multidimarray.");
            if ( NSIZE(data) > 1 )
            {
                REPORT_ERROR(ERR_MMAP,"Image Class::ReadData: mmap with multiple "
                             "images file not compatible. Try selecting a unique image.");
            }
            mappedOffset = selectImgOffset;
            mappedSize   = mappedOffset + pagesize;
            mmapFile();
        }
        else
        {
            char* page = NULL;

            // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
            //if memory already allocated use it (no resize allowed)
            data.coreAllocateReuse();
            //ROB
            //#define DEBUG
#ifdef DEBUG

            data.printShape();
            printf("DEBUG: Page size: %ld offset= %ld \n", pagesize, offset);
            printf("DEBUG: Swap = %d  Pad = %ld  Offset = %ld\n", swap, pad, offset);
            printf("DEBUG: myoffset = %ld select_img= %ld \n", selectImgOffset, selectImgSizeT);
#endif
#undef DEBUG

            if (pagesize > pagemax)
                page = (char *) askMemory(pagemax*sizeof(char));
            else
                page = (char *) askMemory(pagesize*sizeof(char));

            if(fseek( fimg, selectImgOffset, SEEK_SET ) == -1)
                REPORT_ERROR(ERR_IO_SIZE,"readData: can not seek the file pointer");
            for ( size_t myn = 0; myn < NSIZE(data); myn++ )
            {
                for (size_t myj = 0; myj < pagesize; myj += pagemax )//pagesize size of object
                {
                    // Read next page. Divide pages larger than pagemax
                    readsize = pagesize - myj;
                    if ( readsize > pagemax )
                        readsize = pagemax;
                    readsize_n = readsize/datatypesize;

                    //Read page from disc
                    if (fread( page, readsize, 1, fimg )!=1)
                        REPORT_ERROR(ERR_IO_NOREAD,"Cannot read the whole page");
                    //swap per page
                    if (swap)
                        swapPage(page, readsize, datatype, swap);
                    // cast to T per page
                    castPage2T(page, MULTIDIM_ARRAY(data) + haveread_n, datatype, readsize_n);
                    haveread_n += readsize_n;
                }
                if ( pad > 0 )
                    //fread( padpage, pad, 1, fimg);
                    if (fseek(fimg, pad, SEEK_CUR) == -1)
                        REPORT_ERROR(ERR_IO_SIZE,"readData: can not seek the file pointer");
            }
            //if ( pad > 0 )
            //    freeMemory(padpage, pad*sizeof(char));
            if ( page > 0 )
                freeMemory(page, pagesize*sizeof(char));

#ifdef DEBUG

            printf("DEBUG img_read_data: Finished reading and converting data\n");
#endif

        }
        return;
    }

    /* Write the raw date after a data type casting.
     */
    void writeData(FILE* fimg, size_t offset, DataType wDType, size_t datasize_n,
                   CastWriteMode castMode=CW_CAST)
    {
        size_t dTypeSize = gettypesize(wDType);
        size_t datasize  = datasize_n * dTypeSize;
        size_t ds2Write  = rw_max_page_size;
        size_t dsN2Write = rw_max_page_size/dTypeSize;
        size_t rw_max_n  = dsN2Write;

        char* fdata;
        double min0 = 0, max0 = 0;

        if (castMode != CW_CAST)
            data.computeDoubleMinMaxRange(min0, max0, offset, datasize_n);

        if (datasize > rw_max_page_size)
            fdata = (char *) askMemory(rw_max_page_size*sizeof(char));
        else
            fdata = (char *) askMemory(datasize*sizeof(char));


        for (size_t writtenDataN = 0; writtenDataN < datasize_n; writtenDataN += rw_max_n )
        {

            if (writtenDataN + rw_max_n > datasize_n )
            {
                dsN2Write = datasize_n - writtenDataN;
                ds2Write = dsN2Write * dTypeSize;
            }

            if (castMode == CW_CAST)
                castPage2Datatype(MULTIDIM_ARRAY(data)+offset+writtenDataN, fdata, wDType, dsN2Write);
            else
                castConvertPage2Datatype(MULTIDIM_ARRAY(data)+offset+writtenDataN,fdata,wDType,dsN2Write,min0,max0,castMode);

            //swap per page
            if (swapWrite)
                swapPage(fdata, ds2Write, wDType);

            fwrite( fdata, ds2Write, 1, fimg );
        }
        freeMemory(fdata, rw_max_page_size);
    }

    /* Mmap the Image class to an image file.
     */
    void mmapFile()
    {
        if (this->hFile->mode == WRITE_READONLY)
            mFd = open(dataFName.c_str(), O_RDONLY, S_IREAD);
        else
            mFd = open(dataFName.c_str(), O_RDWR, S_IREAD | S_IWRITE);

        if ( mFd == -1 )
        {
            if (errno == EACCES)
                REPORT_ERROR(ERR_IO_NOPERM,formatString("Image Class::mmapFile: permission denied when opening %s",dataFName.c_str()));
            else
                REPORT_ERROR(ERR_IO_NOTOPEN,"Image Class::mmapFile: Error opening the image file to be mapped.");
        }
        char * map;

        if (this->hFile->mode == WRITE_READONLY)
            map = (char*) mmap(0,mappedSize, PROT_READ, MAP_SHARED, mFd, 0);
        else
            map = (char*) mmap(0,mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED, mFd, 0);

        if ( map == (void*) -1 )
            REPORT_ERROR(ERR_MMAP_NOTADDR,"Image Class::ReadData: mmap of image file failed.");
        data.data = reinterpret_cast<T*> (map+mappedOffset);
        data.nzyxdimAlloc = XSIZE(data)*YSIZE(data)*ZSIZE(data)*NSIZE(data);
    }

    /* Munmap the image file.
     */
    void munmapFile()
    {
        munmap((char*)(data.data)-mappedOffset,mappedSize);
        close(mFd);
        data.data = NULL;
        mappedSize = mappedOffset = 0;
    }

    /* Return the datatype of the current image object
     */
    DataType myT()
    {
        if (typeid(T) == typeid(unsigned char))
            return UChar;
        else if (typeid(T) == typeid(char))
            return SChar;
        else if (typeid(T) == typeid(unsigned short))
            return UShort;
        else if (typeid(T) == typeid(short))
            return Short;
        else if (typeid(T) == typeid(unsigned int))
            return UInt;
        else if (typeid(T) == typeid(int))
            return Int;
        else if (typeid(T) == typeid(unsigned int))
            return UInt;
        else if (typeid(T) == typeid(int))
            return Int;
        else if (typeid(T) == typeid(long))
            return Long;
        else if (typeid(T) == typeid(float))
            return Float;
        else if (typeid(T) == typeid(double))
            return Double;
        else if (typeid(T) == typeid(std::complex<short>))
            return ComplexShort;
        else if (typeid(T) == typeid(std::complex<int>))
            return ComplexInt;
        else if (typeid(T) == typeid(std::complex<float>))
            return ComplexFloat;
        else if (typeid(T) == typeid(std::complex<double>))
            return ComplexDouble;
        else if (typeid(T) == typeid(bool))
            return Bool;
        else
            return Unknown_Type;
    }

    /* friend declaration for stacks handling purposes
    */
    friend class ImageCollection;
    template <typename TT>
    friend class Image;

}
;

// Special cases for complex numbers
template<>
void Image< std::complex< double > >::castPage2T(char * page,
        std::complex<double> * ptrDest,
        DataType datatype,
        size_t pageSize);
template<>
void Image< std::complex< double > >::castPage2Datatype(std::complex< double > * srcPtr,
        char * page,
        DataType datatype,
        size_t pageSize);
template<>
void Image< std::complex< double > >::castConvertPage2Datatype(std::complex< double > * srcPtr,
        char * page, DataType datatype, size_t pageSize,double min0,double max0,CastWriteMode castMode);

//@}
#endif
