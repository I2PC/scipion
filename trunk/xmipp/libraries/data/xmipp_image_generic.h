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

#ifndef IMAGE_GENERIC_H_
#define IMAGE_GENERIC_H_

#include "xmipp_image_base.h"
#include "multidim_array_generic.h"

/** @addtogroup Images
//@{

/**
 *  ImageGeneric class to handle images with independence of data type
 **/
class ImageGeneric
{

public:
    ImageBase* image;
    DataType datatype;
    MultidimArrayGeneric * data;

protected:

    bool    swap; // To store the swap mode of the image when the datatype is read.

public:

    /** Empty constructor.
     *
     * No internal image class is declared.
     */
    ImageGeneric()
    {
        init();
    }
    /** Constructor passing the data type of the image.
     *
     * Defines the data type of the image and then declares the internal image class.
     */
    ImageGeneric(DataType _datatype);

    /** Constructor for reading the image
     *
     * This constructor will be useful for perfom construction and read in a single step.
     */
    ImageGeneric(const FileName &filename);

    /** Copy Constructor */
    ImageGeneric(const ImageGeneric &img);

    /** Destructor.
     */
    ~ImageGeneric();

    /** Initialize the parameters.
     */
    void init();

    /** Clear the parameters and initialize them.
     */
    void clear();

    /** Copy data from other ImageGeneric */
    void copy(const ImageGeneric &img);

    /** Clear the image header
     */
    void clearHeader()
    {
        image->clearHeader();
    }

    /** Init geometry transformation with defaults values
     */
    void initGeometry(const size_t n = 0)
    {
        image->MD[n]=MDL::emptyHeader;
    }

    /** Return geometry row
     */
    MDRow& getGeometry(const size_t n = 0)
    {
        return image->MD[n];
    }

    /** Get Image dimensions
    */
    void getDimensions(int &Xdim, int &Ydim, int &Zdim, size_t &Ndim) const
    {
        image->getDimensions(Xdim, Ydim, Zdim, Ndim);
    }
    void getDimensions(int &Xdim, int &Ydim, int &Zdim) const
    {
        size_t Ndim;
        image->getDimensions(Xdim, Ydim, Zdim, Ndim);
    }

    /** Get number of elements in image
     */
    size_t getSize() const
    {
        return NZYXSIZE(*(data->im));
    }

    /** Get Euler angles from image header
    */
    void getEulerAngles(double &rot, double &tilt, double &psi,
                        size_t n = 0)
    {
        image->getEulerAngles(rot, tilt, psi, n);
    }

    /** Get Tilt angle from image header
    */
    double tilt(const size_t n = 0) const
    {
        return image->tilt(n);
    }

    /** Set the data type for the generic image
     */
    void setDatatype(DataType _datatype);

    /** Set image dataMode */
    void setDataMode(DataMode mode)
    {
        image->setDataMode(mode);
    }

    /** Get the data type
     */
    DataType getDatatype()const
    {
        return datatype;
    }

    /** Get the data type
     */
    int getDatatypeDepth()const
    {
        switch (datatype)
        {
        case Float:
        case UInt:
            return 32;
        case Int:
            return 31;
        case Short:
            return 15;
        case UShort:
            return 16;
        case SChar:
            return 7;
        case UChar:
            return 8;
        }
    }

    /** Check if image is mapped on file
      */
    inline bool isMapped()
    {
        return image->isMapped();
    }

    /** Read image from file.
     */
    int read(const FileName &name, DataMode datamode = DATA, size_t select_img = ALL_IMAGES,
             bool mapData = false);

    /** Read image from file with a header applied.
     */
    int readApplyGeo(const FileName &name, const MDRow &row, bool only_apply_shifts = false,
                     DataMode datamode = DATA, size_t select_img = ALL_IMAGES, bool wrap = WRAP);

    /** Read image from file.
     */
    int readApplyGeo(const FileName &name, const MetaData &md, size_t objId, bool only_apply_shifts = false,
                     DataMode datamode = DATA, size_t select_img = ALL_IMAGES, bool wrap = WRAP);

    /** Read an image from metadata, filename is taken from MDL_IMAGE */
    int readApplyGeo(const MetaData &md, size_t objId, bool only_apply_shifts = false,
                     DataMode datamode = DATA, size_t select_img = ALL_IMAGES, bool wrap = WRAP);

    /** Apply geometry in refering metadata to the image */
    void applyGeo(const MetaData &md, size_t objId, bool only_apply_shifts = false, bool wrap = WRAP);

    /** Read image mapped from file.
     */
    int readMapped(const FileName &name, size_t select_img = ALL_IMAGES, int mode = WRITE_READONLY);

    /* Initially try to read normally, but if there is a memory allocation problem, then
     * try to read from the mapped file.*/
    int readOrReadMapped(const FileName &name, size_t select_img = ALL_IMAGES, int mode = WRITE_READONLY);

    /* Read an image with a lower resolution as a preview image.
    * If Zdim parameter is not passed, then all slices are rescaled.
    * If Ydim is not passed, then Ydim is rescaled same factor as Xdim.
    */
    int readPreview(const FileName &name, int Xdim, int Ydim = -1, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE);

    /** Write image to file.
    */
    inline void write(const FileName &name="", size_t select_img = ALL_IMAGES, bool isStack=false,
                      int mode=WRITE_OVERWRITE, CastWriteMode castMode = CW_CAST, int _swapWrite = 0)
    {
        image->write(name,select_img,isStack,mode,castMode,_swapWrite);
    }

    /* Create an empty image file of format given by filename and map it to memory.
     */
    void mapFile2Write(int Xdim, int Ydim, int Zdim, const FileName &_filename,
                       bool createTempFile=false, size_t select_img = APPEND_IMAGE,
                       bool isStack=false, int mode=WRITE_OVERWRITE, int _swapWrite = 0);

    /* MultidimArrayGeneric data access
     */
    inline MultidimArrayGeneric& operator()()
    {
        return *data;
    }

    inline const MultidimArrayGeneric& operator()() const
    {
        return *data;
    }

    /**
     * assign operator
     */
    ImageGeneric& operator=(const ImageGeneric &img);

    /**
     * equal operator
     */
    bool operator==(const ImageGeneric &i1) const;

    /** Get pixel value
     */
    inline double getPixel(unsigned long n, int k, int i, int j) const
    {
#define GETVALUE(type) return NZYX_ELEM(*(MultidimArray<type>*)data->im,n,k,i,j);
        SWITCHDATATYPE(datatype,GETVALUE)
#undef GETVALUE

    }

    /** Get pixel value
     */
    inline double getPixel(int i, int j) const
    {
#define GETVALUE(type) return (double) dAij(*(MultidimArray<type>*)data->im,i,j);
        SWITCHDATATYPE(datatype,GETVALUE)
#undef GETVALUE

    }

    /** Set pixel value
     */
    inline void setPixel(unsigned long n, int k, int i, int j, double value) const
    {
#define SETVALUE(type) NZYX_ELEM(*(MultidimArray<type>*)data->im,n,k,i,j) = (type) value;
        SWITCHDATATYPE(datatype,SETVALUE)
#undef SETVALUE

    }

    /** Set pixel value
     */
    inline void setPixel(int i, int j, double value) const
    {
#define SETVALUE(type) A2D_ELEM(*(MultidimArray<type>*)data->im,i,j) = (type) value;
        SWITCHDATATYPE(datatype,SETVALUE)
#undef SETVALUE

    }
    void print() const;
    void toString(String &s) const;

    friend std::ostream& operator<<(std::ostream& o, const ImageGeneric& I)
    {
        o << I.image;
    }

    void add(const ImageGeneric &img);
    void subtract(const ImageGeneric &img);

private:
    /* Return the datatype of the image file.
     */
    void getImageType(const FileName &imgName, DataType &datatype);

    /* Return the datatype of the image file and if the bytes must be swapped when read.
     */
    void getImageType(const FileName &imgName, DataType &datatype, bool &swap);
}
;

//@}

/** Create an empty image file
    *
    * An image file, which name and format are given by filename,
    * is created. Only the header info is written, and if image number is given, then disk space
    * is reserved until select_img .
    * Swap the endianess of the image header is also possible.
    */
void createEmptyFile(const FileName &_filename, int Xdim, int Ydim, int Zdim = 1,
                     size_t select_img = APPEND_IMAGE, bool isStack = false,
                     int mode = WRITE_OVERWRITE, int _swapWrite = 0);

#endif /* IMAGE_GENERIC_H_ */
