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

#ifndef IMAGE_H
#define IMAGE_H

#include "image_base.h"
#include "datatype.h"
#include "image_generic.h"

/** Template class for images.
 * The image class is the general image handling class.
 */
template<typename T>
class Image: public ImageBase
{
public:
    MultidimArray<T>    data;        // The image data array
    std::vector<MDRow>  MD;                     // data for each subimage
    MDRow               MDMainHeader;           // data for the file

protected:
    FileName            filename;    // File name
    FileName            dataFName;   // Data File name without flags
    FILE*                fimg;        // Image File handler
    FILE*                fhed;        // Image File header handler
    TIFF*                tif;         // TIFF Image file hander
    bool                stayOpen;    // To maintain the image file open after read/write
    int                 dataflag;    // Flag to force reading of the data
    unsigned long       i;           // Current image number (may be > NSIZE)
    size_t       offset;      // Data offset
    int                 swap;        // Perform byte swapping upon reading
    TransformType       transform;   // Transform type
    int                 replaceNsize;// Stack size in the replace case
    bool                 _exists;     // does target file exists?
    // equal 0 is not exists or not a stack
    bool                mmapOnRead;  // Mapping when reading from file
    bool                mmapOnWrite; // Mapping when writing to file
    int                 mFd;         // Handle the file in reading method and mmap
    size_t              mappedSize;  // Size of the mapped file
    size_t              mappedOffset;// Offset for the mapped file

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

    /** Constructor with size and filename
     *
     * An image file, which name and format are given by filename,
     * is created with the given size. Then the image is mapped to this file.
     *
     * @code
     * Image I(64,64,1,1,"image.spi");
     * @endcode
     */
    Image(int Xdim, int Ydim, int Zdim, int Ndim, const FileName _filename)
    {
        init();
        mmapOnWrite = true;
        data.setDimensions(Xdim, Ydim, Zdim, Ndim);
        MD.resize(Ndim);
        filename = _filename;
        ImageFHandler *hFile = openFile(_filename, WRITE_OVERWRITE);
        _write(_filename, hFile, -1, false, WRITE_OVERWRITE);
        closeFile(hFile);
    }

    /** Init.
     * Initialize everything to 0
     */
    void init()
    {
        clearHeader();
        dataflag = -1;
        if (isComplexT())
            transform = Standard;
        else
            transform = NoTransform;
        i = 0;
        filename = "";
        offset = 0;
        swap = 0;
        replaceNsize=0;
        mmapOnRead = mmapOnWrite = false;
        mappedSize = 0;
        mFd    = NULL;
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

    /** Clear the header of the image
     */
    void clearHeader()
    {
        MDMainHeader.clear();
        MD.clear();
        //Just to ensure there is an empty MDRow
        MD.push_back(MDMainHeader);
    }

    /** Check whether image is complex based on T
     */
    bool isComplexT() const
    {
        return ( typeid(T) == typeid(std::complex<double>) ||
                 typeid(T) == typeid(std::complex<float>) );
    }

    /** Check whether image is complex based on transform
      */
    bool isComplex() const
    {
        return !(transform==NoTransform);
    }

    /** Destructor.
     */
    ~Image()
    {
        if (mmapOnRead || mmapOnWrite)
            munmapFile();
        else
            data.clear();
    }

    /** Is this file an image
     *
     *  Check whether a real-space image can be read
     *
     */
    bool isImage(const FileName &name)
    {
        return !read(name, false);
    }

    /** Is this file a real-valued image
     *
     *  Check whether a real-space image can be read
     *
     */
    bool isRealImage(const FileName &name)
    {
        return (isImage(name) && !isComplex());
    }

    /** Is this file a complex image
     *
     *  Check whether a fourier-space (complex) image can be read
     *
     */
    bool isComplexImage(const FileName &name)
    {
        return (isImage(name) && isComplex());
    }

    /** Rename the image
      */
    void rename (const FileName &name)
    {
        filename = name;
    }

    /** Create a mapped image file
     *
     * An image file, which name and format are given by filename,
     * is created with the given size. Then the image is mapped to this file.
     *
     * @code
     * Image I(64,64,1,1,"image.spi");
     * @endcode
     */
    void newMappedFile(int Xdim, int Ydim, int Zdim, int Ndim, const FileName _filename)
    {
        clear();
        mmapOnWrite = true;
        data.setDimensions(Xdim, Ydim, Zdim, Ndim);
        MD.resize(Ndim);
        filename = _filename;
        ImageFHandler *hFile = openFile(_filename, WRITE_OVERWRITE);
        _write(_filename, hFile, -1, false, WRITE_OVERWRITE);
        closeFile(hFile);
    }

    /** General read function
     * you can read a single image from a single image file
     * or a single image file from an stack, in the second case
     * the select slide may come in the image name or in the select_img parameter
     * file name takes precedence over select_img
     * If -1 is given the whole object is read
     *
     */
    int read(const FileName &name, bool readdata=true, int select_img = -1,
             bool apply_geo = false, bool only_apply_shifts = false,
             MDRow * row = NULL, bool mapData = false)
    {
        ImageFHandler* hFile = openFile(name);
        int err = _read(name, hFile, readdata, select_img, apply_geo, only_apply_shifts, row, mapData);
        closeFile(hFile);

        return err;
    }

    /** Read an image from metadata*/

    /* Read an image with a lower resolution as a preview image.
     * If Zdim parameter is not passed, then all slices are rescaled.
     */
    int readPreview(const FileName &name, int Xdim, int Ydim, int Zdim = -1, int select_img = 0)
    {
        ImageGeneric im;
        int imXdim, imYdim, imZdim;

        im.readMapped(name, select_img);
        im.getDimensions(imXdim, imYdim, imZdim);
        im().setXmippOrigin();

        scaleToSize(0,IMGMATRIX(*this),im(),Xdim,Ydim,(Zdim != -1)? Zdim:imZdim);
    }

    /** General write function
     * select_img= which slice should I replace
     * overwrite = 0, append slice
     * overwrite = 1 overwrite slice
     */
    void write(const FileName &name="", int select_img=-1, bool isStack=false,
               int mode=WRITE_OVERWRITE,bool adjust=false)
    {
        // If image is already mapped to file then close the file and clear.
        if (mmapOnWrite && mappedSize > 0)
        {
            munmapFile();
            return;
        }

        const FileName &fname = (name == "") ? filename : name;

        /* If the filename is in stack we will suppose you want to write this,
         * even if you have not set the flags to.
         */
        if (fname.isInStack() && isStack == false && mode == WRITE_OVERWRITE)
        {
            isStack = true;
            mode = WRITE_APPEND;
        }

        ImageFHandler* hFile = openFile(fname, mode);
        _write(fname, hFile, select_img, isStack, mode, adjust);
        closeFile(hFile);
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

    /* Convert the image to another datatype and save it
        *
        */
    void castConvertPage2Datatype(T * srcPtr, char * page, DataType datatype, size_t pageSize,
                                  double min0, double max0, CastWriteMode castMode=CONVERT)
    {

        double minF, maxF;
        double slope;
        unsigned long int n;
        DataType  myTypeId = myT();

        switch(datatype)
        {
        case UChar:
            {
                if (castMode==CONVERT && myTypeId == SChar)
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
                    if (castMode==CONVERT &&  myTypeId == UChar)
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
                    if (castMode==CONVERT && (myTypeId == SChar|| myTypeId == Short))
                    {
                        slope = 1;
                        min0 -= SHRT_MIN;
                    }
                    else if (castMode==CONVERT && (myTypeId == UChar))
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
                    if (castMode==CONVERT && (myTypeId == UChar || myTypeId == UShort))
                    {
                        slope = 1;
                        min0 += SHRT_MIN;
                    }
                    else if (castMode==CONVERT && (myTypeId == SChar))
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
                    if (castMode==CONVERT &&  (myTypeId == SChar
                                                   || myTypeId == Short
                                                   || myTypeId == Int))
                    {
                        slope = 1;
                        min0 -= INT_MIN;
                    }
                    else if (castMode==CONVERT &&  (myTypeId == UShort || myTypeId == UChar))
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
                    if (castMode==CONVERT &&  (myTypeId == UChar
                                                   || myTypeId == UShort
                                                   || myTypeId == UInt))
                    {
                        slope = 1;
                        min0 += INT_MIN;
                    }
                    else if (castMode==CONVERT &&  (myTypeId == Short || myTypeId == SChar))
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
    /** Check file Datatype is same as T type to use mmap.
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
                if (typeid(T) == typeid(signed char))
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

    /** Swap an entire page
      * input pointer char *
      */
    void swapPage(char * page, size_t pageNrElements, DataType datatype)
    {
        unsigned long datatypesize = gettypesize(datatype);
#ifdef DEBUG

        std::cerr<<"DEBUG swapPage: Swapping image data with swap= "
        << swap<<" datatypesize= "<<datatypesize
        << " pageNrElements " << pageNrElements
        << " datatype " << datatype
        <<std::endl;
        ;
#endif

        // Swap bytes if required
        if ( swap == 1 )
        {
            if ( datatype >= ComplexShort )
                datatypesize /= 2;
            for ( unsigned long i=0; i<pageNrElements; i+=datatypesize )
                swapbytes(page+i, datatypesize);
        }
        else if ( swap > 1 )
        {
            for ( unsigned long i=0; i<pageNrElements; i+=swap )
                swapbytes(page+i, swap);
        }
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
    /** Set pixel
     * (direct access) needed by swig
     */
    void setPixel(int i, int j, T v)
    {
        IMGPIXEL(*this,i,j)=v;
    }

    /** Get pixel
     * (direct acces) needed by swig
     */
    T getPixel(int i, int j) const
    {
        return IMGPIXEL(*this,i,j);
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

    /** Get file name
     *
     * @code
     * std::cout << "Image name = " << I.name() << std::endl;
     * @endcode
     */
    const FileName & name() const
    {
        return filename;
    }

    /** Get Image dimensions
     */
    void getDimensions(int &Xdim, int &Ydim, int &Zdim, unsigned long &Ndim) const
    {
        Xdim = XSIZE(data);
        Ydim = YSIZE(data);
        Zdim = ZSIZE(data);
        Ndim = NSIZE(data);
    }

    long unsigned int getSize() const
    {
        return NZYXSIZE(data);
    }

    /** Get Image offset and swap
     */
    void getOffsetAndSwap(size_t &_offset, int &_swap) const
    {
        _offset = offset;
        _swap = swap;
    }

    /* Is there label in the individual header */
    bool individualContainsLabel(MDLabel label) const
    {
        return (!MD.empty() && MD[0].containsLabel(label));
    }

    /* Is there label in the main header */
    bool mainContainsLabel(MDLabel label) const
    {
        return MDMainHeader.containsLabel(label);
    }

    /** Get Rot angle
    *
    * @code
    * std::cout << "First Euler angle " << I.rot() << std::endl;
    * @endcode
    */
    double rot(const long int n = 0) const
    {
        double dummy = 0;
        MD[n].getValue(MDL_ANGLEROT, dummy);
        return dummy;
    }

    /** Get Tilt angle
     *
     * @code
     * std::cout << "Second Euler angle " << I.tilt() << std::endl;
     * @endcode
     */
    double tilt(const long int n = 0) const
    {
        double dummy = 0;
        MD[n].getValue(MDL_ANGLETILT, dummy);
        return dummy;
    }

    /** Get Psi angle
     *
     * @code
     * std::cout << "Third Euler angle " << I.psi() << std::endl;
     * @endcode
     */
    double psi(const long int n = 0) const
    {
        double dummy = 0;
        MD[n].getValue(MDL_ANGLEPSI, dummy);
        return dummy;
    }

    /** Get Xoff
     *
     * @code
     * std::cout << "Origin offset in X " << I.Xoff() << std::endl;
     * @endcode
     */
    double Xoff(const long int n = 0) const
    {
        double dummy = 0;
        MD[n].getValue(MDL_ORIGINX, dummy);
        return dummy;
    }

    /** Get Yoff
     *
     * @code
     * std::cout << "Origin offset in Y " << I.Yoff() << std::endl;
     * @endcode
     */
    double Yoff(const long int n = 0) const
    {
        double dummy = 0;
        MD[n].getValue(MDL_ORIGINY, dummy);
        return dummy;
    }

    /** Get Zoff
     *
     * @code
     * std::cout << "Origin offset in Z " << I.Zoff() << std::endl;
     * @endcode
     */
    double Zoff(const long int n = 0) const
    {
        double dummy = 0;
        MD[n].getValue(MDL_ORIGINZ, dummy);
        return dummy;
    }

    /** Get Weight
    *
    * @code
    * std::cout << "weight= " << I.weight() << std::endl;
    * @endcode
    */
    double weight(const long int n = 0) const
    {
        double dummy = 1;
        MD[n].getValue(MDL_WEIGHT, dummy);
        return dummy;
    }

    /** Get Scale factor
    *
    * @code
    * std::cout << "scale= " << I.scale() << std::endl;
    * @endcode
    */
    double scale(const long int n = 0) const
    {
        double dummy = 1;
        MD[n].getValue(MDL_SCALE, dummy);
        return dummy;
    }


    /** Get Flip
    *
    * @code
    * std::cout << "flip= " << flip() << std::endl;
    * @endcode
    */
    bool flip(const long int n = 0) const
    {
        bool dummy = false;
        MD[n].getValue(MDL_FLIP, dummy);
        return dummy;
    }

    /** Data type
        *
        * @code
        * std::cout << "datatype= " << dataType() << std::endl;
        * @endcode
        */
    DataType dataType() const
    {
        int dummy;
        MDMainHeader.getValue(MDL_DATATYPE, dummy);
        return (DataType)dummy;
    }

    /** Sampling RateX
    *
    * @code
    * std::cout << "sampling= " << samplingRateX() << std::endl;
    * @endcode
    */
    double samplingRateX() const
    {
        double dummy = 1.;
        MDMainHeader.getValue(MDL_SAMPLINGRATEX, dummy);
        return dummy;
    }

    /** Set file name
     */
    void setName(const FileName &_filename)
    {
        filename = _filename;
    }

    /** Set Euler angles in image header
     */
    void setEulerAngles(double rot, double tilt, double psi,
                        long int n = 0)
    {
        MD[n].setValue(MDL_ANGLEROT, rot);
        MD[n].setValue(MDL_ANGLETILT, tilt);
        MD[n].setValue(MDL_ANGLEPSI, psi);
    }

    /** Get Euler angles from image header
     */
    void getEulerAngles(double &rot, double &tilt, double &psi,
                        long int n = 0)
    {
        MD[n].getValue(MDL_ANGLEROT, rot);
        MD[n].getValue(MDL_ANGLETILT, tilt);
        MD[n].getValue(MDL_ANGLEPSI, psi);
    }

    /** Set Rotation angle to image */
    void setRot(double rot, long int n = 0)
    {
        MD[n].setValue(MDL_ANGLEROT, rot);
    }

    /** Set Tilt angle to image */
    void setTilt(double tilt, long int n = 0)
    {
        MD[n].setValue(MDL_ANGLETILT, tilt);
    }

    /** Set Rotation angle to image */
    void setPsi(double psi, long int n = 0)
    {
        MD[n].setValue(MDL_ANGLEPSI, psi);
    }

    /** Set origin offsets in image header
     */
    void setShifts(double xoff, double yoff, double zoff = 0.,
                   long int n = 0)
    {
        MD[n].setValue(MDL_ORIGINX, xoff);
        MD[n].setValue(MDL_ORIGINY, yoff);
        MD[n].setValue(MDL_ORIGINZ, zoff);
    }
    /** Get origin offsets from image header
      */
    void getShifts(double &xoff, double &yoff, double &zoff,
                   long int n = 0)
    {
        MD[n].getValue(MDL_ORIGINX, xoff);
        MD[n].getValue(MDL_ORIGINY, yoff);
        MD[n].getValue(MDL_ORIGINZ, zoff);
    }

    /** Set X offset in image header
     */
    void setXoff(double xoff, long int n = 0)
    {
        MD[n].setValue(MDL_ORIGINX, xoff);
    }

    /** Set Y offset in image header
     */
    void setYoff(double yoff, long int n = 0)
    {
        MD[n].setValue(MDL_ORIGINY, yoff);
    }

    /** Set Z offset in image header
     */
    void setZoff(double zoff, long int n = 0)
    {
        MD[n].setValue(MDL_ORIGINZ, zoff);
    }

    /** Set scale in image header
     */
    void setScale(double scale, long int n = 0)
    {
        MD[n].setValue(MDL_SCALE, scale);
    }

    /** Get scale from image header
     */
    void getScale(double &scale, long int n = 0)
    {
        MD[n].getValue(MDL_SCALE, scale);
    }

    /** Set flip in image header
     */
    void setFlip(bool flip, long int n = 0)
    {
        MD[n].setValue(MDL_FLIP, flip);
    }

    /** Set Weight in image header
    */
    void setWeight(double weight, long int n = 0)
    {
        MD[n].setValue(MDL_WEIGHT, weight);
    }

    /** Get geometric transformation matrix from 2D-image header
      */
    void getTransformationMatrix(Matrix2D<double> &A,
        bool only_apply_shifts = false,
        long int n = 0)
    {
        // This has only been implemented for 2D images...
        MULTIDIM_ARRAY(*this).checkDimension(2);

        const MDRow &rowAux=MD[n];

        double phi,psi,theta,xoff,yoff,scale;
        bool flip;
        rowAux.getValue(MDL_ANGLEROT, phi);
        phi = realWRAP(phi, 0., 360.);
        rowAux.getValue(MDL_ANGLETILT, theta);
        theta = realWRAP(theta, 0., 360.);
        rowAux.getValue(MDL_ANGLEPSI, psi);
        psi = realWRAP(psi, 0., 360.);
        rowAux.getValue(MDL_ORIGINX, xoff);
        rowAux.getValue(MDL_ORIGINY, yoff);
        rowAux.getValue(MDL_SCALE, scale);

        A.initIdentity(3);
        if (only_apply_shifts)
        {
            Euler_angles2matrix(0., 0., 0., A);
            MAT_ELEM(A, 0, 2) = -xoff;
            MAT_ELEM(A, 1, 2) = -yoff;
        }
        else
        {
            if (theta == 0.)
            {
                // For untilted images: apply Euler matrix
                Euler_angles2matrix(phi, 0., psi, A);
            }
            else
            {
                // For tilted images: only apply Psi
                // Take another_set into account
                if (theta < 0.)
                {
                    theta = -theta;
                    psi = realWRAP(psi - 180., -180, 180);
                }
                Euler_angles2matrix(0., 0., psi, A);
            }
            MAT_ELEM(A, 0, 2) = -xoff;
            MAT_ELEM(A, 1, 2) = -yoff;
        }
        A *= scale;

        // Also for only_apply_shifts: mirror if necessary!
        rowAux.getValue(MDL_FLIP, flip);

        if (flip)
        {
            MAT_ELEM(A, 0, 0) *= -1.;
            MAT_ELEM(A, 0, 1) *= -1.;
        }
    }

    /** Show image properties
      */
    friend std::ostream& operator<<(std::ostream& o, const Image<T>& I)
    {
        o << "Image type   : ";
        if (I.isComplex())
            o << "Fourier-space image" << std::endl;
        else
            o << "Real-space image" << std::endl;

        o << "Reversed     : ";
        if (I.swap)
            o << "TRUE"  << std::endl;
        else
            o << "FALSE" << std::endl;

        o << "Data type    : ";
        switch (I.dataType())
        {
        case Unknown_Type:
            o << "Undefined data type";
            break;
        case UChar:
            o << "Unsigned character or byte type (UInt8)";
            break;
        case SChar:
            o << "Signed character (Int8)";
            break;
        case UShort:
            o << "Unsigned short integer (UInt16)";
            break;
        case Short:
            o << "Signed short integer (Int16)";
            break;
        case UInt:
            o << "Unsigned integer (UInt32)";
            break;
        case Int:
            o << "Signed integer (Int32)";
            break;
        case Long:
            o << "Signed integer (4 or 8 byte, depending on system)";
            break;
        case Float:
            o << "Floating point (4-byte)";
            break;
        case Double:
            o << "Double precision floating point (8-byte)";
            break;
        case ComplexShort:
            o << "Complex two-byte integer (4-byte)";
            break;
        case ComplexInt:
            o << "Complex integer (8-byte)";
            break;
        case ComplexFloat:
            o << "Complex floating point (8-byte)";
            break;
        case ComplexDouble:
            o << "Complex floating point (16-byte)";
            break;
        case Bool:
            o << "Boolean (1-byte?)";
            break;
        }
        o << std::endl;

        o << "dimensions   : " << NSIZE(I()) << " x " << ZSIZE(I()) << " x " << YSIZE(I()) << " x " << XSIZE(I());
        o << "  (noObjects x slices x rows x columns)" << std::endl;
        if (I.individualContainsLabel(MDL_ANGLEROT))
        {
            o << "Euler angles : " << std::endl;
            o << "  Phi   (rotation around Z axis) = " << I.rot() << std::endl;
            o << "  theta (tilt, second rotation around new Y axis) = " << I.tilt() << std::endl;
            o << "  Psi   (third rotation around new Z axis) = " << I.psi() << std::endl;
        }
        if (I.individualContainsLabel(MDL_ORIGINX))
        {
            o << "Origin Offsets : " << std::endl;
            o << "  Xoff  (origin offset in X-direction) = " << I.Xoff() << std::endl;
            o << "  Yoff  (origin offset in Y-direction) = " << I.Yoff() << std::endl;
            o << "  Zoff  (origin offset in Z-direction) = " << I.Zoff() << std::endl;
        }
        if(I.individualContainsLabel(MDL_SCALE))
            o << "Scale  : " <<I.scale() << std::endl;
        o << "Header size  : " << I.offset << std::endl;
        if (I.individualContainsLabel(MDL_WEIGHT))
            o << "Weight  : " << I.weight() << std::endl;
        if (I.individualContainsLabel(MDL_FLIP))
            o << "Flip    : " << I.flip() << std::endl;
        return o;
    }

    /** Sum this object with other file and keep in this object
      */
    void sumWithFile(const FileName &fn)
    {
        Image<T> aux;
        aux.read(fn);
        (*this)()+=aux();
    }


    /**
     *  Specific read functions for different file formats
     */
#include "rwDM3.h"
#include "rwIMAGIC.h"
#include "rwMRC.h"
#include "rwINF.h"
#include "rwRAW.h"
#include "rwSPIDER.h"
#include "rwSPE.h"
#include "rwTIA.h"
#include "rwTIFF.h"

private:

    /** Open file function
      * Open the image file and returns its file hander.
      */
    ImageFHandler* openFile(const FileName &name, int mode = WRITE_READONLY) const
    {
        ImageFHandler* hFile = new ImageFHandler;
        FileName fileName, headName = "";
        FileName ext_name = name.getFileFormat();

        int dump;
        name.decompose(dump, fileName);

        fileName = fileName.removeFileFormat();

        size_t found = fileName.find_first_of("%");
        if (found!=std::string::npos)
            fileName = fileName.substr(0, found) ;

        hFile->exist = exists(fileName);

        std::string wmChar;

        switch (mode)
        {
        case WRITE_READONLY:
            if (!hFile->exist)
                REPORT_ERROR(ERR_IO_NOTEXIST,(std::string) "Cannot read file "
                             + fileName + ". It does not exist" );
            wmChar = "r";
            break;
        case WRITE_OVERWRITE:
            wmChar = "w";
            break;
        case WRITE_APPEND:
        case WRITE_REPLACE:
            if (hFile->exist)
                wmChar = "r+";
            else
                wmChar = "w+";
            break;

        }

        if (ext_name.contains("tif"))
        {
            TIFFSetWarningHandler(NULL); // Switch off warning messages
            if ((hFile->tif = TIFFOpen(fileName.c_str(), wmChar.c_str())) == NULL)
                REPORT_ERROR(ERR_IO_NOTOPEN,"rwTIFF: There is a problem opening the TIFF file.");
            hFile->fimg = NULL;
            hFile->fhed = NULL;
        }
        else
        {
            hFile->tif = NULL;

            if (ext_name.contains("img") || ext_name.contains("hed"))
            {
                fileName = fileName.withoutExtension();
                headName = fileName.addExtension("hed");
                fileName = fileName.addExtension("img");
            }
            else if (ext_name.contains("raw"))
            {
                if (mode != WRITE_READONLY || exists(fileName.addExtension("inf")) )
                {
                    headName = fileName.addExtension("inf");
                    ext_name = "inf";
                }
                else
                    ext_name = "raw";
            }
            else if (ext_name.contains("inf"))
            {
                headName = fileName;
                fileName = fileName.withoutExtension();
                ext_name = "inf";
            }

            // Open image file
            if ( ( hFile->fimg = fopen(fileName.c_str(), wmChar.c_str()) ) == NULL )
                REPORT_ERROR(ERR_IO_NOTOPEN,(std::string)"Image::openFile cannot open: " + fileName);

            if (headName != "")
            {
                if ( ( hFile->fhed = fopen(headName.c_str(), wmChar.c_str()) ) == NULL )
                    REPORT_ERROR(ERR_IO_NOTOPEN,(std::string)"Image::openFile cannot open: " + headName);
            }
            else
                hFile->fhed = NULL;

        }
        hFile->fileName = fileName;
        hFile->headName = headName;
        hFile->ext_name = ext_name;

        return hFile;
    }

    /** Close file function.
      * Close the image file according to its name and file handler.
      */
    void closeFile(ImageFHandler* hFile = NULL)
    {
        FileName ext_name;
        FILE* fimg, *fhed;
        TIFF* tif;

        if (hFile != NULL)
        {
            ext_name = hFile->ext_name;
            fimg = hFile->fimg;
            fhed = hFile->fhed;
            tif  = hFile->tif;
        }
        else
        {
            ext_name = filename.getFileFormat();
            fimg = this->fimg;
            fhed = this->fhed;
            tif  = this->tif;
        }

        if (ext_name.contains("tif"))
            TIFFClose(tif);
        else
        {
            if (fclose(fimg) != 0 )
                REPORT_ERROR(ERR_IO_NOCLOSED,(std::string)"Can not close image file "+ filename);

            if (fhed != NULL &&  fclose(fhed) != 0 )
                REPORT_ERROR(ERR_IO_NOCLOSED,(std::string)"Can not close header file of "
                             + filename);
        }
        delete hFile;
    }

    /* Internal read image file method.
     */
    int _read(const FileName &name, ImageFHandler* hFile, bool readdata=true, int select_img = -1,
              bool apply_geo = false, bool only_apply_shifts = false,
              MDRow * row = NULL, bool mapData = false)
    {
        //const MetaData &docFile = *docFilePtr;
        //std::vector<MDLabel> &activeLabels = *activeLabelsPtr;

        int err = 0;

        // Check whether to read the data or only the header
        dataflag = ( readdata ) ? 1 : -1;

        // If Image has been previously used with mmap, then close the previous file
        if (mappedSize != 0)
            munmapFile();

        // Check whether to map the data or not
        mmapOnRead = mapData;

        FileName ext_name = hFile->ext_name;
        fimg = hFile->fimg;
        fhed = hFile->fhed;
        tif  = hFile->tif;

        int dump;
        name.decompose(dump, filename);
        filename = name;
        dataFName = hFile->fileName;

        if (dump != -1)
            select_img = dump;

#undef DEBUG
        //#define DEBUG
#ifdef DEBUG

        std::cerr << "READ\n" <<
        "name="<<name <<std::endl;
        std::cerr << "ext= "<<ext_name <<std::endl;
        std::cerr << " now reading: "<< filename <<" dataflag= "<<dataflag
        << " select_img "  << select_img << std::endl;
#endif
#undef DEBUG

        //Just clear the header before reading
        MDMainHeader.clear();
        //Set the file pointer at beginning
        if (fimg != NULL)
            fseek(fimg, 0, SEEK_SET);
        if (fhed != NULL)
            fseek(fhed, 0, SEEK_SET);

        if (ext_name.contains("spi") || ext_name.contains("xmp")  ||
            ext_name.contains("stk") || ext_name.contains("vol"))//mrc stack MUST go BEFORE plain MRC
            err = readSPIDER(select_img);
        else if (ext_name.contains("mrcs"))//mrc stack MUST go BEFORE plain MRC
            err = readMRC(select_img,true);
        else if (ext_name.contains("mrc"))//mrc
            err = readMRC(select_img,false);
        else if (ext_name.contains("img") || ext_name.contains("hed"))//
            err = readIMAGIC(select_img);//imagic is always an stack
        else if (ext_name.contains("ser"))//TIA
            err = readTIA(select_img,false);
        else if (ext_name.contains("dm3"))//DM3
            err = readDM3(select_img,false);
        else if (ext_name.contains("inf"))//RAW with INF file
            err = readINF(select_img,false);
        else if (ext_name.contains("raw"))//RAW without INF file
            err = readRAW(select_img,false);
        else if (ext_name.contains("tif"))//TIFF
            err = readTIFF(select_img,false);
        else if (ext_name.contains("spe"))//SPE
            err = readSPE(select_img,false);
        else
            err = readSPIDER(select_img);
        //This implementation does not handle stacks,
        //read in a block
        if (row != NULL)
        {
            if (data.ndim != 1)
                REPORT_ERROR(ERR_MULTIDIM_SIZE, "Header overwriting not available for stacks!!!");
            MDLabel label;
            MDRow &rowAux=MD[0];

            for (MDRow::const_iterator it = row->begin(); it != row->end(); ++it)
            {
                label = (*it)->label;
                if (rowAux.containsLabel(label))
                    *(rowAux.getObject(label)) = *(*it);
                else
                  rowAux.push_back(new MDObject(*(*it)));
            }
        }

        //apply geo has not been defined for volumes
        if(this->data.getDim()>2)
            apply_geo=false;

        if (readdata && (apply_geo || only_apply_shifts))
        {
            Matrix2D< double > A;
            getTransformationMatrix(A,only_apply_shifts);
            if (!A.isIdentity())
            {
                MultidimArray<T> tmp=MULTIDIM_ARRAY(*this);
                applyGeometry(BSPLINE3, MULTIDIM_ARRAY(*this), tmp,
                    A, IS_INV, WRAP);
            }
        }

        // Negative errors are bad.
        return err;
    }

    /* Internal write image file method.
     */
    void _write(const FileName &name, ImageFHandler* hFile, int select_img=-1,
                bool isStack=false, int mode=WRITE_OVERWRITE, bool adjust=false)
    {
        int err = 0;

        // if image is mapped to file then close the file and clear

        if (mmapOnWrite && mappedSize > 0)
        {
            munmapFile();
            return;
        }

        filename = name;
        dataFName = hFile->fileName;
        _exists = hFile->exist;
        fimg = hFile->fimg;
        fhed = hFile->fhed;
        tif  = hFile->tif;

        FileName ext_name = hFile->ext_name;

        int aux;
        FileName filNamePlusExt;
        name.decompose(aux, filNamePlusExt);

        if (select_img == -1)
            select_img = aux;

        size_t found = filNamePlusExt.find_first_of("%");

        std::string imParam = "";

        if (found!=std::string::npos)
        {
            imParam =  filNamePlusExt.substr(found+1).c_str();
            filNamePlusExt = filNamePlusExt.substr(0, found) ;
        }

        found = filNamePlusExt.find_first_of(":");
        if ( found!=std::string::npos)
            filNamePlusExt   = filNamePlusExt.substr(0, found);


        //#define DEBUG
#ifdef DEBUG

        std::cerr << "write" <<std::endl;
        std::cerr<<"extension for write= "<<ext_name<<std::endl;
        std::cerr<<"filename= "<<filename<<std::endl;
        std::cerr<<"mode= "<<mode<<std::endl;
        std::cerr<<"isStack= "<<isStack<<std::endl;
        std::cerr<<"select_img= "<<select_img<<std::endl;
#endif
    #undef DEBUG
        // Check that image is not empty
        if (getSize() < 1)
            REPORT_ERROR(ERR_MULTIDIM_EMPTY,"write Image ERROR: image is empty!");

        // CHECK FOR INCONSISTENCIES BETWEEN data.xdim and x, etc???
        int Xdim, Ydim, Zdim;
        unsigned long Ndim;
        Xdim=Ydim=Zdim=Ndim=0;
        if (_exists)
            this->getDimensions(Xdim,Ydim, Zdim, Ndim);

        Image<T> auxI;
        replaceNsize=0;//reset replaceNsize in case image is reused
        if(select_img == -1 && mode == WRITE_REPLACE)
            REPORT_ERROR(ERR_VALUE_INCORRECT,"Please specify object to be replaced");
        else if (_exists && (mode == WRITE_REPLACE || mode == WRITE_APPEND))
        {
            auxI.dataflag = -2;
            auxI._read(filNamePlusExt, hFile, false);
            int _Xdim, _Ydim, _Zdim;
            unsigned long _Ndim;
            auxI.getDimensions(_Xdim,_Ydim, _Zdim, _Ndim);
            replaceNsize=_Ndim;
            if(Xdim!=_Xdim ||
               Ydim!=_Ydim ||
               Zdim!=_Zdim)
            {
                std::cerr << "(x,y,z) " << Xdim << " " << Ydim << " " << Zdim << " "<< Ndim << std::endl;
                std::cerr << "(_x,_y,_z) " << _Xdim << " " << _Ydim << " " << _Zdim << " " <<_Ndim <<std::endl;
                REPORT_ERROR(ERR_MULTIDIM_SIZE,"write: target and source objects have different size");
            }
            if(mode==WRITE_REPLACE && select_img>_Ndim)
                replaceNsize = select_img;
            if(auxI.replaceNsize <1 &&
               (mode==WRITE_REPLACE || mode==WRITE_APPEND))
                REPORT_ERROR(ERR_IO,"write: output file is not an stack");
        }
        else if(!_exists && mode==WRITE_APPEND)
        {
            ;
        }
        else if (mode == WRITE_READONLY)//If new file we are in the WRITE_OVERWRITE mode
        {
            REPORT_ERROR(ERR_ARG_INCORRECT, (std::string) "File " + name
                         + " opened in read-only mode. Cannot write.");
        }
        /*
         * SELECT FORMAT
         */
        //Set the file pointer at beginning
        if (fimg != NULL)
            fseek(fimg, 0, SEEK_SET);
        if (fhed != NULL)
            fseek(fhed, 0, SEEK_SET);

        if(ext_name.contains("spi") || ext_name.contains("xmp") ||
           ext_name.contains("stk") || ext_name.contains("vol"))
            err = writeSPIDER(select_img,isStack,mode);
        else if (ext_name.contains("mrcs"))
            writeMRC(select_img,true,mode);
        else if (ext_name.contains("mrc"))
            writeMRC(select_img,false,mode,imParam,adjust);
        else if (ext_name.contains("img") || ext_name.contains("hed"))
            writeIMAGIC(select_img,mode,imParam,adjust);
        else if (ext_name.contains("dm3"))
            writeDM3(select_img,false,mode);
        else if (ext_name.contains("ser"))
            writeTIA(select_img,false,mode);
        else if (ext_name.contains("raw") || ext_name.contains("inf"))
            writeINF(select_img,false,mode,imParam,adjust);
        else if (ext_name.contains("tif"))
            writeTIFF(select_img,isStack,mode,imParam,adjust);
        else if (ext_name.contains("spe"))
            writeSPE(select_img,isStack,mode);
        else
            err = writeSPIDER(select_img,isStack,mode);

        if ( err < 0 )
        {
            std::cerr << " Filename = " << filename << " Extension= " << ext_name << std::endl;
            REPORT_ERROR(ERR_IO_NOWRITE, "Error writing file");
        }

        /* If initially the file did not existed, once the first image is written,
         * then the file exists
         */
        if (!_exists)
            hFile->exist = _exists = true;
    }

    /** Read the raw data
      */
    void readData(FILE* fimg, int select_img, DataType datatype, unsigned long pad)
    {
        //#define DEBUG
#ifdef DEBUG
        std::cerr<<"entering readdata"<<std::endl;
        std::cerr<<" readData flag= "<<dataflag<<std::endl;
#endif
#undef DEBUG

        if ( dataflag < 1 )
            return;

        // If only half of a transform is stored, it needs to be handled
        if (transform == Hermitian || transform == CentHerm )
            data.setXdim(XSIZE(data)/2 + 1);

        size_t selectImgOffset, readsize, readsize_n, pagemax = 1073741824; //1Gb
        size_t datatypesize=gettypesize(datatype);
        size_t pagesize  =ZYXSIZE(data)*datatypesize;
        size_t haveread_n=0;
        size_t selectImgSizeT=0;

        // Reset select to get the correct offset
        selectImgSizeT = ( select_img < 0 )? 0 : (size_t) select_img;

        selectImgOffset = offset + selectImgSizeT*(pagesize + pad);

        // Flag to know that data is not going to be mapped although mmapOn is true
        if (mmapOnRead && !checkMmapT(datatype))
        {
            std::cout << "WARNING: Image Class. File datatype and image declaration not compatible with mmap. Loading into memory." <<std::endl;
            mmapOnRead = false;
            mFd = -1;
        }

        if (mmapOnRead)
        {
            // Image mmapOn is not compatible with Multidimarray mmapOn
            if(data.mmapOn)
                REPORT_ERROR(ERR_MULTIDIM_DIM,"Image Class::ReadData: mmap option can not be selected simoultanesouslly\
                             for both Image class and its Multidimarray.");
            if ( NSIZE(data) > 1 )
            {
                REPORT_ERROR(ERR_MMAP,"Image Class::ReadData: mmap with multiple "
                             "images file not compatible. Try selecting a unique image.");
            }
            //            fclose(fimg);
            mappedOffset = selectImgOffset;
            mappedSize = mappedOffset + pagesize;
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

            if(fseek( fimg, selectImgOffset, SEEK_SET )==-1)
                REPORT_ERROR(ERR_IO_SIZE,"readData: can not seek the file pointer");
            for ( size_t myn=0; myn<NSIZE(data); myn++ )
            {
                for (size_t myj=0; myj<pagesize; myj+=pagemax )//pagesize size of object
                {
                    // Read next page. Divide pages larger than pagemax
                    readsize = pagesize - myj;
                    if ( readsize > pagemax )
                        readsize = pagemax;
                    readsize_n = readsize/datatypesize;

                    //Read page from disc
                    fread( page, readsize, 1, fimg );
                    //swap per page
                    if (swap)
                        swapPage(page, readsize, datatype);
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
                   CastWriteMode castMode=CAST)
    {
        size_t datasize = datasize_n * gettypesize(wDType);
        char* fdata = (char *) askMemory(datasize);

        switch(castMode)
        {
        case CAST:
            {
                castPage2Datatype(MULTIDIM_ARRAY(data)+offset, fdata, wDType, datasize_n);
                break;
            }
        default:
            {
                double min0, max0;
                data.computeDoubleMinMaxRange(min0, max0, offset, datasize_n);
                castConvertPage2Datatype(MULTIDIM_ARRAY(data)+offset, fdata, wDType, datasize_n,min0 ,max0, castMode);
            }
        }

        fwrite( fdata, datasize, 1, fimg );
        freeMemory(fdata, datasize);
    }

    /* Mmap the Image class to an image file.
     */
    void mmapFile()
    {
        if ( ( mFd = open(dataFName.c_str(), O_RDWR, S_IREAD | S_IWRITE) ) == -1 )
            REPORT_ERROR(ERR_IO_NOTOPEN,"Image Class::ReadData: Error opening the image file.");

        char * map;

        if ( (map = (char*) mmap(0,mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED, mFd, 0)) == (void*) -1 )
            REPORT_ERROR(ERR_MMAP_NOTADDR,"Image Class::ReadData: mmap of image file failed.");
        data.data = reinterpret_cast<T*> (map+mappedOffset);
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

/// @defgroup ImageFormats Image Formats
/// @ingroup Images
// Functions belonging to this topic are commented in rw*.h
//@}
#endif
