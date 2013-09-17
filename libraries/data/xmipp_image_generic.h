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

// @addtogroup Images
//@{

#define CHECK_IMG(op)\
 if (image == NULL) \
  REPORT_ERROR(ERR_IMG_UNKNOWN, "ImageGeneric: Image not initialized.");\
 else\
 {\
  op\
 }

/**
 * ImageGeneric class to handle images with independence of data type
 */

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
    void getDimensions(size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim) const;
    void getDimensions(size_t &Xdim, size_t &Ydim, size_t &Zdim) const;
    void getDimensions(ArrayDim &aDim) const;

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

    /** Resizes image dimensions
     */
    void resize(int Xdim, int Ydim, int Zdim, size_t Ndim, bool copy=true)
    {
        image->setDimensions(Xdim, Ydim, Zdim, Ndim);
        data->resize(Ndim, Zdim, Ydim, Xdim, copy);
    }

    /** Set the data type for the generic image
     */
    void setDatatype(DataType _datatype);

    /** Set the data type according to the image file, checking if
     *  file size is correct.
     */
    void setDatatype(const FileName &name, ImageInfo &imgInf);
    void setDatatype(const FileName &name);

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

    /** Get basic information from already read image file
     */
    void getInfo(ImageInfo &imgInfo) const;

    /** Get the data type
     */
    int getDatatypeDepth()const
    {
        switch (datatype)
        {
        case DT_Float:
        case DT_UInt:
            return 32;
        case DT_Int:
            return 31;
        case DT_Short:
            return 15;
        case DT_UShort:
            return 16;
        case DT_SChar:
            return 7;
        case DT_UChar:
            return 8;
        default:
        	REPORT_ERROR(ERR_TYPE_INCORRECT,"Do not know how to handle this type at this point");
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
    int readApplyGeo(const FileName &name, const MDRow &row,
                     const ApplyGeoParams &params = DefaultApplyGeoParams);

    /** Read image from file.
     */
    int readApplyGeo(const FileName &name, const MetaData &md, size_t objId,
                     const ApplyGeoParams &params = DefaultApplyGeoParams);

    /** Read an image from metadata, filename is taken from MDL_IMAGE */
    int readApplyGeo(const MetaData &md, size_t objId,
                     const ApplyGeoParams &params = DefaultApplyGeoParams);

    /** Apply geometry in refering metadata to the image */
    void applyGeo(const MetaData &md, size_t objId,
                  const ApplyGeoParams &params = DefaultApplyGeoParams);
    /* Euler mirror Y ---------------------------------------------------------- */
    void mirrorY(void);

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
    int readPreview(const FileName &name, size_t Xdim, size_t Ydim = 0, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE);

    /** This function allows to read the original image or a preview of it also allowing to select either
     *  a specific image from the stack or a slice from a volume.
     *
     *  In the case of reading images in its real dimensions it is also possible to image map from file.
     */
    int readOrReadPreview(const FileName &name, size_t Xdim, size_t Ydim = 0, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE,
                          bool mapData = false, bool wrap = true);

    /** Returns an image with a lower resolution as a preview image.
     * If Zdim parameter is not passed, then all slices are rescaled.
     * If Ydim is not passed, then Ydim is rescaled same factor as Xdim.
     */
    void getPreview(ImageGeneric &imgOut, int Xdim, int Ydim = -1, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE);

    /* Read an image with a lower resolution as a preview image.
    * Resizing is done in fourier space
    * If Zdim parameter is not passed, then all slices are rescaled.
    * If Ydim is not passed, then Ydim is rescaled same factor as Xdim.
    */
    int readPreviewFourier(const FileName &name, size_t Xdim, size_t Ydim = 0, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE);

    /* Using xvsmooth */
    int readPreviewSmooth(const FileName &name, size_t Xdim, size_t Ydim = 0, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE);


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

    /** It changes the behavior of the internal multidimarray so it points to a specific slice/image
      *  from a stack, volume or stack of volumes. No information is deallocated from memory, so it is
      *  also possible to repoint to the whole stack,volume... (passing select_slice = ALL_SLICES and
      *  selec_img = ALL_IMAGES).
      *
      *  The options for select_slice are:
      *
      *    - a slice number,
      *    - CENTRAL_SLICE, to automatically select the central slice of the volume,
      *    - ALL_SLICES, to recover the whole volume.
      *
      *  The options for selec_img are:
      *
      *    - a image number of the stack,
      *    - ALL_IMAGES, to recover the whole stack.
      *
      *  If a specific slice number is selected, then a specific image from the stack must be
      *  also selected. Otherwise, FIRST_IMAGE is proposed.
      *
      *  If Image Object is read using readPreview method, movePointerTo only works when rescaling
      *  the image in X-Y plane only, but all slices must be read.
      */
    void movePointerTo(int select_slice = ALL_SLICES, size_t select_img = ALL_IMAGES)
    {
        image->movePointerTo(select_slice, select_img);
    }


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

    /** Convert the datatype of the object and cast the image
    */
    void convert2Datatype(DataType datatype, CastWriteMode castMode=CW_CONVERT);

    /* Reslice a volume aligning any X or Y direction with Z axis.
     */
    void reslice(AxisView view);
    /* The resliced volume is returned in out
    */
    void reslice(AxisView view, ImageGeneric &out);

    /**
     * equal for doubles
     */
    bool equal(const ImageGeneric &i1, double accuracy=XMIPP_EQUAL_ACCURACY) const
    {
        return data->equal(MULTIDIM_ARRAY_GENERIC(i1),accuracy);
    }
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
#define GETVALUE(type) return (double) A2D_ELEM(*(MultidimArray<type>*)data->im,i,j);
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
        CHECK_IMG(SWITCHDATATYPE(datatype,SETVALUE))
#undef SETVALUE

    }

    /** Init constant
     */
    inline void initConstant(double value) const
    {
#define INITCONS(type) (*(MultidimArray<type>*)(data->im)).initConstant((type) value);
        CHECK_IMG(SWITCHDATATYPE(datatype,INITCONS))
#undef INITCONS

    }

    /** Init random
     */
    inline void initRandom(double op1, double op2, RandomMode mode = RND_UNIFORM) const
    {
#define INITRND(type) (*(MultidimArray<type>*)(data->im)).initRandom(op1, op2, mode);
        CHECK_IMG(SWITCHDATATYPE(datatype,INITRND))
#undef INITRND

    }

    /** Print image information
     */
    void print() const;

    /** Add to string the image information
     */
    void toString(String &s) const;

    friend std::ostream& operator<<(std::ostream& o, const ImageGeneric& I)
    {
        o << I.image;
        return o;
    }

    /** Addition of the passed image to the internal's
     */
    void add(const ImageGeneric &img);

    /** Subtraction of the passed image to the internal's
     */
    void subtract(const ImageGeneric &img);

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
