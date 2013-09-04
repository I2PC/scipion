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

#ifndef IMAGE_BASE_H_
#define IMAGE_BASE_H_

#include "xmipp_image_macros.h"
#include "multidim_array.h"
#include "transformations.h"
#include "metadata.h"
#include "xmipp_datatype.h"
//
//// Includes for rwTIFF which cannot be inside it
#include "../../external/tiff-3.9.4/libtiff/tiffio.h"
#include "../../external/hdf5-1.8.10/src/hdf5.h"


/* Minimum size of a TIFF file to be mapped to a tempfile in case of mapping from
 * image file is required
 */
const size_t tiff_map_min_size = 0x12000000;

/// @defgroup Images Images
/// @ingroup DataLibrary

//@{
/** Transform type.
 *  This type defines the kind of image.
 */
typedef enum
{
    NoTransform = 0,        // No transform
    Standard = 1,           // Standard transform: origin = (0,0,0)
    Centered = 2,           // Centered transform: origin = (nx/2,ny/2,nz/2)
    Hermitian = 3,          // Hermitian half: origin = (0,0,0)
    CentHerm = 4            // Centered hermitian: origin = (0,ny/2,nz/2)
} TransformType;

/** Write mode
 * This class defines the writing behavior.
 */
typedef enum
{
    WRITE_READONLY,   //only can read the file
    WRITE_OVERWRITE, //forget about the old file and overwrite it
    WRITE_REPLACE,   //replace a particular object by another
    WRITE_APPEND,    //append and object at the end of a stack, so far can not append stacks
    WRITE_LAST_LABEL                       // **** NOTE ****: Do keep this label always at the end
    // it is here for looping purposes
} WriteMode;

/** Data mode
 * This enumerate which data will be read/write from image files
 * We can read/write complete image or only its headers(w/o geometrical information)
 */
typedef enum
{
    _NONE = -2,  // Nothing to do. Used by ImageGeneric to check the right datatype to be used
    HEADER = -1, //Dont read image data, only info from main header(datatype and dimensions)
    _HEADER_ALL = 0, //Read complete header(main and geo), useful for header_extract and header_assign
    DATA = 1, //Read image data and main header, geometrical transformations will be ignored
    _DATA_ALL = 2  //Read data with complete header(the use of this option is not recommended, all Xmipp
    // programs should read and write geo info through metadatas
}DataMode;

/* Cast Write mode
 * This enum defines the cast writing behavior
 */
typedef enum
{
    // prefix needed so extract_image_enums.py script can create the equivalent class for Java
    CW_CAST,       //Only cast the data type
    CW_CONVERT,    //Convert the data from one type to another
    CW_ADJUST,     //Adjust the histogram to fill the gray level range
    CW_LAST_LABEL                       // **** NOTE ****: Do keep this label always at the end
    // it is here for looping purposes
} CastWriteMode;

/** Open File struct
 * This struct is used to share the File handlers with Image Collection class
 */
struct ImageFHandler
{
    FILE*     fimg;       // Image File handler
    FILE*     fhed;       // Image File header handler
    TIFF*     tif;        // TIFF Image file handler
    hid_t     fhdf5;	  // HDF5 File handler
    FileName  fileName;   // Image file name
    FileName  headName;   // Header file name
    FileName  ext_name;   // Filename extension
    bool     exist;       // Shows if the file exists. Equal 0 means file does not exist or not stack.
    int        mode;   // Opening mode behavior
};

struct ImageInfo
{
    size_t    offset;
    DataType  datatype;
    bool      swap;
    ArrayDim  adim;
};

struct ApplyGeoParams
{
    bool only_apply_shifts;
    DataMode datamode;
    size_t select_img;
    bool wrap;

    ApplyGeoParams()
    {
        only_apply_shifts = false;
        datamode = DATA;
        select_img = ALL_IMAGES;
        wrap = WRAP;
    }
};

const ApplyGeoParams DefaultApplyGeoParams;

/// @name Images Speed-up
/// @{

/** Volume Matrix access.
 *
 * This macro does the same as the normal 3D matrix access but in a faster way
 * as no function call is generated.
 *
 * @code
 * VOLMATRIX(V).resize(128, 128, 128);
 *
 * VOLMATRIX(V2) = VOLMATRIX(V1) + VOLMATRIX(V2);
 * @endcode
 */
#define VOLMATRIX(V) ((V).data)

/** Image Matrix access.
 *
 * This macro does the same as the normal 2D matrix access but in a faster way
 * as no function call is generated.
 *
 * @code
 * IMGMATRIX(V).resize(128, 128);
 *
 * IMGMATRIX(V2) = IMGMATRIX(V1) + IMGMATRIX(V2);
 * @endcode
 */
#define IMGMATRIX(I) ((I).data)

/** Pixel access.
 * For fast access to pixel values (and for backwards compatibility of the code)
 */
#define IMGPIXEL(I, i, j) A2D_ELEM(((I).data), (i), (j))

/** Physical pixel access.
 *
 * The physical pixel access gives you access to a pixel by its physical
 * position and not by its logical one. This access shouldn't be used as a
 * custom, use instead the logical access, but there might be cases in which
 * this access might be interesting. Physical positions start at index 0 in C.
 *
 * @code
 * std::cout << "This is the first pixel stored in the Image " <<
 *     DIRECT_IMGPIXEL(V, 0, 0) << std::endl;
 * @endcode
 */
#define DIRECT_IMGPIXEL(I, i, j) DIRECT_A2D_ELEM(((I).data), (i), (j))

/** Voxel access.
 *
 * This macro does the same as the normal voxel access (remember, logical
 * access) but in a faster way as no function call is generated.
 *
 * @code
 * std::cout << "Grey level of voxel (2,-3,-3) of the Volume = " <<
 *     VOLVOXEL(V, 2, -3, -3) << std::endl;
 *
 * VOLVOXEL(I, 2, -3, -3) = VOLVOXEL(I, 2, -3, -2);
 * @endcode
 */
#define VOLVOXEL(V, k, i, j) A3D_ELEM(((V).data), (k), (i), (j))

/** Physical voxel access.
 *
 * The physical voxel access gives you access to a voxel by its physical
 * position and not by its logical one. This access shouldn't be used as a
 * custom, use instead the logical access, but there might be cases in which
 * this access might be interesting. Physical positions start at index 0 in C.
 *
 * @code
 * std::cout << "This is the first voxel stored in the Volume " <<
 *     DIRECT_VOLVOXEL(V, 0, 0, 0) << std::endl;
 * @endcode
 */
#define DIRECT_VOLVOXEL(I, k, i, j) DIRECT_A3D_ELEM(((I).data), (k), (i), (j))
///@}

/** Swapping trigger.
 * Threshold file z size above which bytes are swapped.
 */
#define SWAPTRIG     16776960


/// Image base class
class ImageBase
{
public:
    MultidimArrayBase * mdaBase;    // Pointer to data from Image<template T> casted as MultidimArrayBase
    std::vector<MDRow>  MD;                     // data for each subimage
    MDRow               MDMainHeader;           // data for the file

protected:
    FileName            filename;    // File name
    FileName            tempFilename; // Temporary filename
    FileName            dataFName;   // Data File name without flags
    FILE*               fimg;        // Image File handler
    FILE*               fhed;        // Image File header handler
    TIFF*               tif;         // TIFF Image file hander
    hid_t				fhdf5;       // HDF5 File handler
    ImageFHandler*      hFile;       // Image File handler information structure
    ArrayDim        	aDimFile;   // Image header file information structure (original info from file)
    DataMode            dataMode;    // Flag to force select what will be read/write from image files
    size_t              offset;      // Data offset
    int                 swap;        // Perform byte swapping upon reading
    int                 swapWrite;   // Perform byte swapping upon writing
    TransformType       transform;   // Transform type
    size_t              replaceNsize;// Stack size in the replace case
    bool                _exists;     // does target file exists?  // equal 0 is not exists or not a stack
    bool                mmapOnRead;  // Mapping when reading from file
    bool                mmapOnWrite; // Mapping when writing to file
    int                 mFd;         // Handle the file in reading method and mmap
    size_t              mappedSize;  // Size of the mapped file
    size_t              mappedOffset;// Offset for the mapped file
    size_t        		virtualOffset;// MDA Offset when movePointerTo is used

public:

    /** Init.
     * Initialize everything to 0*/
    void init();

    /** Clear.
     * Initialize everything to 0*/
    virtual void clear()=0;

    /** Clear the header of the image*/
    void clearHeader();

    /** Check whether image is complex based on T*/
    virtual bool isComplexT() const=0;

    /** Check whether image is complex based on transform*/
    bool isComplex() const
    {
        return !(transform==NoTransform);
    }

    /** Destructor.*/
    virtual ~ImageBase()
    {}

    /** Is this file an image
     *
     *  Check whether a real-space image can be read*/
    bool isImage(const FileName &name)
    {
        try
        {
            return !read(name, HEADER);
        }
        catch (XmippError &xe)
        {
            return false;
        }
    }

    /** Check if image is mapped on file
      */
    inline bool isMapped()
    {
        return (mmapOnRead || mmapOnWrite);
    }

    /** Is this file a real-valued image
     *
     *  Check whether a real-space image can be read
     */
    bool isRealImage(const FileName &name)
    {
        return (isImage(name) && !isComplex());
    }

    /** Is this file a complex image
     *
     *  Check whether a fourier-space (complex) image can be read
     */
    bool isComplexImage(const FileName &name)
    {
        return (isImage(name) && isComplex());
    }

    /** Rename the image*/
    void rename (const FileName &name)
    {
        filename = name;
    }

    /** Create a mapped image file
     *
     * An image file, which name and format are given by filename,
     * is created with the given size. Then the image is mapped to this file.
     * The image object must be cleared prior to use this method.
     */
    void mapFile2Write(size_t Xdim, size_t Ydim, size_t Zdim, const FileName &_filename,
                       bool createTempFile = false, size_t select_img = APPEND_IMAGE,
                       bool isStack = false, int mode = WRITE_OVERWRITE);

    /** General read function
     * you can read a single image from a single image file
     * or a single image file from an stack, in the second case
     * the select slide may come in the image name or in the select_img parameter
     * file name takes precedence over select_img
     * If ALL_IMAGES is given the whole object is read.
     *
     * Parameter mapData allows to access to image mapped to disk instead of loaded into
     * memory. In case the mapped image is intended to be modified the parameter
     * mode must be WRITE_REPLACE.
     *
     * This function cannot apply geometrical transformations.
     */
    int read(const FileName &name, DataMode datamode = DATA, size_t select_img = ALL_IMAGES,
             bool mapData = false, int mode = WRITE_READONLY);

    /** General read function
     * you can read a single image from a single image file
     * or a single image file from an stack, in the second case
     * the select slide may come in the image name or in the select_img parameter
     * file name takes precedence over select_img
     * If -1 is given the whole object is read
     *
     * This function can apply geometrical transformations, but cannot map the image in disk
     */
    int readApplyGeo(const FileName &name, const MDRow &row,
                     const ApplyGeoParams &params = DefaultApplyGeoParams);

    /** Read an image from metadata, filename is provided*/
    int readApplyGeo(const FileName &name, const MetaData &md, size_t objId,
                     const ApplyGeoParams &params = DefaultApplyGeoParams);

    /** Read an image from metadata, filename is taken from MDL_IMAGE */
    int readApplyGeo(const MetaData &md, size_t objId,
                     const ApplyGeoParams &params = DefaultApplyGeoParams);

    /** Apply geometry in refering metadata to the image */
    void applyGeo(const MetaData &md, size_t objId,
                  const ApplyGeoParams &params = DefaultApplyGeoParams);

    /* Read an image with a lower resolution as a preview image.
     * If Zdim parameter is not passed, then all slices are rescaled.
     * If Ydim is not passed, then Ydim is rescaled same factor as Xdim.
     */

    /* Read image mapped from file */
    int readMapped(const FileName &name, size_t select_img = ALL_IMAGES, int mode = WRITE_READONLY);

    /* Initially try to read normally, but if there is a memory allocation problem, then
     * try to read from the mapped file.*/
    int readOrReadMapped(const FileName &name, size_t select_img = ALL_IMAGES, int mode = WRITE_READONLY);

    /* Read a thumb/preview image version from file with other dimensions than image originals'.
     * It is also possible to select an specific image from stack or slice from volume.
     * This function is intended for previews of great image files as the image is not copy to memory.
     */
    virtual int readPreview(const FileName &name, size_t Xdim, size_t Ydim=0, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE) = 0;

    /** Returns an image with a lower resolution as a preview image.
      * If Zdim parameter is not passed, then all slices are rescaled.
      * If Ydim is not passed, then Ydim is rescaled same factor as Xdim.
      */
    virtual void getPreview(ImageBase *imgOut, size_t Xdim, size_t Ydim=0, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE) = 0;

    /** This function allows to read the original image or a preview of it also allowing to select either
     *  a specific image from the stack or a slice from a volume.
     *
     *  In the case of reading images in its real dimensions it is also possible to image map from file.
     */
    int readOrReadPreview(const FileName &name, size_t Xdim, size_t Ydim, int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE, bool mapData = false);

    /** General write function
     * select_img= which slice should I replace
     * overwrite = 0, append slice
     * overwrite = 1 overwrite slice
     */
    void write(const FileName &name="", size_t select_img = ALL_IMAGES, bool isStack=false,
               int mode=WRITE_OVERWRITE,CastWriteMode castMode = CW_CAST, int _swapWrite = 0);

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
    virtual void movePointerTo(int select_slice = ALL_SLICES, size_t select_img = ALL_IMAGES) = 0;

    /* Return the datatype of the current image object */
    virtual DataType myT() const = 0;

    /** Check file Datatype is same as T type to use mmap. */
    virtual bool checkMmapT(DataType datatype)=0;

    /** Write an entire page as datatype
     *
     * A page of datasize_n elements T is cast to datatype and written to fimg
     * The memory for the casted page is allocated and freed internally.
     */
    virtual void writePageAsDatatype(FILE * fimg, DataType datatype, size_t datasize_n )=0;

    /** Swap an entire page
      * input pointer char *
      */
    void swapPage(char * page, size_t pageNrElements, DataType datatype, int swap = 1);

    /* Force the swapping of the endianess of the file when writing, except for the case
     * when appending or modifying the file.
     */
    void swapOnWrite()
    {
        swapWrite = true;
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

    /** Get dimensions of the multidimArray inside image.
     *  TODO: This method must be changed to return the size
     *  of the image read from file, i.e. aDimFile, and where this is used
     *  should be used the imageBase::mda->getDimensions instead.
     */
    void getDimensions(size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim) const;
    /** Get Image dimensions
     */
    void getDimensions(ArrayDim &aDim)
    {
        aDim = aDimFile;
    }
    ArrayDim getDimensions()
    {
        return aDimFile;
    }

    /** Get basic information from already read image file
     */
    void getInfo(ImageInfo &imgInfo) const;

    /** Get basic information from image file
     */
    void getInfo(const FileName &name, ImageInfo &imgInfo);

    /** Get whole number of pixels
     */
    virtual size_t getSize() const = 0;

    /** Get Image offset and swap
     */
    void getOffsetAndSwap(size_t &_offset, int &_swap) const
    {
        _offset = offset;
        _swap = swap;
    }

    /** Get Image swap
         */
    int getSwap() const
    {
        return swap;
    }

    /** Return geometry row
     */
    MDRow& getGeometry(const size_t n = 0)
    {
        return MD[n];
    }

    /** Init geometry transformation with defaults values
     */
    void initGeometry(const size_t n = 0)
    {
        MD[n]=MDL::emptyHeader;
    }

    /* Check if the label is in the individual header
    */
    bool individualContainsLabel(MDLabel label) const
    {
        return (!MD.empty() && MD[0].containsLabel(label));
    }

    /* Check if the label is in the main header */
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
    double rot(const size_t n = 0) const;

    /** Get Tilt angle
     *
     * @code
     * std::cout << "Second Euler angle " << I.tilt() << std::endl;
     * @endcode
     */
    double tilt(const size_t n = 0) const;

    /** Get Psi angle
     *
     * @code
     * std::cout << "Third Euler angle " << I.psi() << std::endl;
     * @endcode
     */
    double psi(const size_t n = 0) const;

    /** Get Xoff
     *
     * @code
     * std::cout << "Origin offset in X " << I.Xoff() << std::endl;
     * @endcode
     */
    double Xoff(const size_t n = 0) const;

    /** Get Yoff
     *
     * @code
     * std::cout << "Origin offset in Y " << I.Yoff() << std::endl;
     * @endcode
     */
    double Yoff(const size_t n = 0) const;

    /** Get Zoff
     *
     * @code
     * std::cout << "Origin offset in Z " << I.Zoff() << std::endl;
     * @endcode
     */
    double Zoff(const size_t n = 0) const;

    /** Get Weight
    *
    * @code
    * std::cout << "weight= " << I.weight() << std::endl;
    * @endcode
    */
    double weight(const size_t n = 0) const;

    /** Get Scale factor
    *
    * @code
    * std::cout << "scale= " << I.scale() << std::endl;
    * @endcode
    */
    double scale(const size_t n = 0) const;

    /** Get Flip
    *
    * @code
    * std::cout << "flip= " << flip() << std::endl;
    * @endcode
    */
    bool flip(const size_t n = 0) const;

    /** Get datatype info from MDMainHeader.
     *
     *  In theory, it must be the image file data type.
     *
     * @code
     * std::cout << "datatype= " << datatype() << std::endl;
     * @endcode
     */
    DataType datatype() const;

    /** Sampling RateX
    *
    * @code
    * std::cout << "sampling= " << samplingRateX() << std::endl;
    * @endcode
    */
    double samplingRateX() const;

    /** Set image dimensions */
    virtual void setDimensions(int Xdim, int Ydim, int Zdim, size_t Ndim) = 0;
    virtual void setDimensions(ArrayDim &aDim)
    {
        mdaBase->setDimensions(aDim);
        aDimFile = aDim;
    }


    /** Set the information of the original image dimensions
     */
    void setADimFile(ArrayDim aDim)
    {
        aDimFile = aDim;
    }
    /** Set file name
     */
    void setName(const FileName &_filename)
    {
        filename = _filename;
    }

    /** Set the read/write dataMode of the object
     */
    void setDataMode(DataMode mode)
    {
        dataMode = mode;
    }

    /** Set Euler angles in image header
     */
    void setEulerAngles(double rot, double tilt, double psi,
                        const size_t n = 0);

    /** Get Euler angles from image header
     */
    void getEulerAngles(double &rot, double &tilt, double &psi,
                        const size_t n = 0) const;

    /** Set Rotation angle to image */
    void setRot(double rot, const size_t n = 0)
    {
        MD[n].setValue(MDL_ANGLE_ROT, rot);
    }

    /** Set Tilt angle to image */
    void setTilt(double tilt, const size_t n = 0)
    {
        MD[n].setValue(MDL_ANGLE_TILT, tilt);
    }

    /** Set Rotation angle to image */
    void setPsi(double psi, const size_t n = 0)
    {
        MD[n].setValue(MDL_ANGLE_PSI, psi);
    }

    /** Set origin offsets in image header
     */
    void setShifts(double xoff, double yoff, double zoff = 0.,
                   const size_t n = 0);

    /** Get origin offsets from image header
      */
    void getShifts(double &xoff, double &yoff, double &zoff,
                   const size_t n = 0) const;

    /** Set X offset in image header
     */
    void setXoff(double xoff, const size_t n = 0)
    {
        MD[n].setValue(MDL_SHIFT_X, xoff);
    }

    /** Set Y offset in image header
     */
    void setYoff(double yoff, const size_t n = 0)
    {
        MD[n].setValue(MDL_SHIFT_Y, yoff);
    }

    /** Set Z offset in image header
     */
    void setZoff(double zoff, const size_t n = 0)
    {
        MD[n].setValue(MDL_SHIFT_Z, zoff);
    }

    /** Set scale in image header
     */
    void setScale(double scale, const size_t n = 0)
    {
        MD[n].setValue(MDL_SCALE, scale);
    }

    /** Get scale from image header
     */
    void getScale(double &scale, const size_t n = 0)
    {
        MD[n].getValue(MDL_SCALE, scale);
    }

    /** Set flip in image header
     */
    void setFlip(bool flip, const size_t n = 0)
    {
        MD[n].setValue(MDL_FLIP, flip);
    }

    /** Set Weight in image header
    */
    void setWeight(double weight, const size_t n = 0)
    {
        MD[n].setValue(MDL_WEIGHT, weight);
    }

    /** Get geometric transformation matrix from 2D-image header
      */
    virtual void getTransformationMatrix(Matrix2D<double> &A,
                                         bool only_apply_shifts = false,
                                         const size_t n = 0)=0;

    /** Sum this object with other file and keep in this object
      */
    virtual void sumWithFile(const FileName &fn)=0;

    /** flip multidim array around X axis
     *
     */
    virtual void mirrorY() =0;

protected:

    /**
     *  Specific read functions for different file formats
     */
    /// @defgroup ImageFormats Image Formats
    /// @ingroup Images
    // Functions belonging to this topic are commented in rw*.h
#include "rwDM3.h"
#include "rwIMAGIC.h"
#include "rwMRC.h"
#include "rwINF.h"
#include "rwRAW.h"
#include "rwSPIDER.h"
#include "rwSPE.h"
#include "rwTIA.h"
#include "rwJPEG.h"
#include "rwTIFF.h"
#include "rwEM.h"
#include "rwPIF.h"
#include "rwHDF5.h"

    /// ----------------------------------------------------------

    /** Open file function
      * Open the image file and returns its file hander.
      */
    ImageFHandler* openFile(const FileName &name, int mode = WRITE_READONLY) const;

    /** Close file function.
      * Close the image file according to its name and file handler.
      */
    void closeFile(ImageFHandler* hFile = NULL) const;

    /** Internal apply geometrical transformations */
    virtual void applyGeo(const MDRow &row, bool only_apply_shifts = false, bool wrap = WRAP) = 0;

    /** Internal read image file method.
     */
    int _read(const FileName &name, ImageFHandler* hFile, DataMode datamode = DATA, size_t select_img = ALL_IMAGES,
              bool mapData = false);

    /** Internal write image file method.
     */
    void _write(const FileName &name, ImageFHandler* hFile, size_t select_img = ALL_IMAGES,
                bool isStack=false, int mode=WRITE_OVERWRITE,CastWriteMode castMode = CW_CAST);

    /** Read the raw data
      */
    virtual void readData(FILE* fimg, size_t select_img, DataType datatype, size_t pad) = 0;

    /** Write the raw date after a data type casting.
     */
    virtual void writeData(FILE* fimg, size_t offset, DataType wDType, size_t datasize_n,
                           CastWriteMode castMode=CW_CAST) = 0;

    virtual void setPage2T(size_t offset, char * page, DataType datatype, size_t pageSize ) = 0;
    virtual void getPageFromT(size_t offset, char * page, DataType datatype, size_t pageSize ) = 0;
    virtual void getCastConvertPageFromT(size_t offset, char * page, DataType datatype, size_t pageSize, double min0, double max0, CastWriteMode castMode=CW_CONVERT) const = 0;

    /** Mmap the Image class to an image file.
     */
    virtual void mmapFile() = 0;

    /** Munmap the image file.
     */
    virtual void munmapFile() = 0;

    /** Set datatype info in MDMainHeader
     */
    void setDatatype(DataType datatype);

    /** Show ImageBase */
    friend std::ostream& operator<<(std::ostream& o, const ImageBase& I);

};
//@}
#endif /* IMAGE_BASE_H_ */
