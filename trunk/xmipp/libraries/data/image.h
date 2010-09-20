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

#include <typeinfo>
#include "funcs.h"
#include "memory.h"
#include "multidim_array.h"
#include "transformations.h"
#include "metadata.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

// Includes for rwTIFF which cannot be inside it
#include "../../external/tiff-3.9.4/libtiff/tiffio.h"

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

/** Data type.
 * This class defines the datatype of the data inside this image.
 */
typedef enum
{
    Unknown_Type = 0,       // Undefined data type
    UChar = 1,              // Unsigned character or byte type
    SChar = 2,              // Signed character (for CCP4)
    UShort = 3,             // Unsigned integer (2-byte)
    Short = 4,              // Signed integer (2-byte)
    UInt = 5,               // Unsigned integer (4-byte)
    Int = 6,                // Signed integer (4-byte)
    Long = 7,               // Signed integer (4 or 8 byte, depending on system)
    Float = 8,              // Floating point (4-byte)
    Double = 9,             // Double precision floating point (8-byte)
    ComplexShort = 10,      // Complex two-byte integer (4-byte)
    ComplexInt = 11,        // Complex integer (8-byte)
    ComplexFloat = 12,      // Complex floating point (8-byte)
    ComplexDouble = 13,     // Complex floating point (16-byte)
    Bool = 14               // Boolean (1-byte?)
} DataType;

/** Write mode
 * This class defines the writing behavior.
 */
typedef enum
{
    WRITE_OVERWRITE, //forget about the old file and overwrite it
    WRITE_APPEND,    //append and object at the end of a stack, so far can not append stacks
    WRITE_REPLACE    //replace a particular object by another
} WriteMode;

/// Returns memory size of datatype
unsigned long gettypesize(DataType type);

/// @name ImagesSpeedUp Images Speed-up
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
//@}

/** Swapping trigger.
 * Threshold file z size above which bytes are swapped.
 */
#define SWAPTRIG     65535

/** Template class for images.
 * The image class is the general image handling class.
 */
template<typename T>
class Image
{
public:
    MultidimArray<T>    data;        // The image data array
    std::vector<MDRow>  MD;                     // data for each subimage
    MDRow               MDMainHeader;           // data for the file

private:
    FileName            filename;    // File name
    int                 dataflag;    // Flag to force reading of the data
    unsigned long       i;           // Current image number (may be > NSIZE)
    unsigned long       offset;      // Data offset
    int                 swap;        // Perform byte swapping upon reading
    TransformType       transform;   // Transform type
    int                 replaceNsize;// Stack size in the replace case
    int                 _exists;     // does target file exists?
    // equal 0 is not exists or not a stack
    bool                mmapOn;      // Mapping when loading from file
    int                 mFd;         // Handle the file in reading method and mmap
    size_t              mappedSize;  // Size of the mapped file

    ///Helper functions for formats read/write
    void castTiffTile2T(
        T * ptrDest ,
        char* tif_buf,
        unsigned int x, unsigned int y,
        unsigned int imageWidth, unsigned int imageLength,
        unsigned int tileWidth, unsigned int tileLength,
        unsigned short samplesPerPixel,
        DataType datatype);
    void castTiffLine2T(
        T * ptrDest,
        char* tif_buf,
        unsigned int y,
        unsigned int imageWidth, unsigned int imageLength,
        unsigned short samplesPerPixel,
        DataType datatype);

public:
    /** Empty constructor
     *
     * An empty image is created.
     *
     * @code
     * Image<double> I;
     * @endcode
     */
    Image();

    /** Constructor with size
     *
     * A blank image (0.0 filled) is created with the given size. Pay attention
     * to the dimension order: Y and then X.
     *
     * @code
     * Image I(64,64);
     * @endcode
     */
    Image(int Xdim, int Ydim, int Zdim=1, int Ndim=1);

    /** Clear.
     * Initialize everything to 0
     */
    void clear();

    /** Clear the header of the image
     */
    void clearHeader();

    /** Check whether image is complex based on T
     */
    bool isComplexT() const;

    /** Check whether image is complex based on transform
      */
    bool isComplex() const;

    /** Destructor.
     */
    ~Image();

    /** Specific read functions for different file formats
      */
    int readDM3(int img_select,bool isStack=false);
    int writeDM3(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE);

    int  readIMAGIC(int img_select);
    int  writeIMAGIC(int img_select=-1, int mode=WRITE_OVERWRITE);

    int readMRC(int img_select, bool isStack=false);
    int writeMRC(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE);

    int readINF(int img_select,bool isStack=false);
    int writeINF(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE);

    int readRAW(int img_select,bool isStack=false);

    int  readSPIDER(int img_select,bool isStack=false);
    int  writeSPIDER(int select_img=-1, bool isStack=false, int mode=WRITE_OVERWRITE);

    int readSPE(int img_select,bool isStack=false);
    int writeSPE(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE);

    int readTIA(int img_select,bool isStack=false, double dStddev=5);
    int writeTIA(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE);

    int readTIFF(int img_select, bool isStack=false);
    int writeTIFF(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE, int imParam=NULL);

//#include "rwDM3.h"
//#include "rwIMAGIC.h"
//#include "rwMRC.h"
//#include "rwINF.h"
//#include "rwRAW.h"
//#include "rwSPIDER.h"
//#include "rwSPE.h"
//#include "rwTIA.h"
//#include "rwTIFF.h"

    /** Is this file an image
     *
     *  Check whether a real-space image can be read
     *
     */
    bool isImage(const FileName &name);

    /** Is this file a real-valued image
     *
     *  Check whether a real-space image can be read
     *
     */
    bool isRealImage(const FileName &name);

    /** Is this file a complex image
     *
     *  Check whether a fourier-space (complex) image can be read
     *
     */
    bool isComplexImage(const FileName &name);

    /** Rename the image
      */
    void rename (const FileName &name);

    /** General read function
     * you can read a single image from a single image file
     * or a single image file from an stack, in the second case
     * the select slide may come in the image name or in the select_img parameter
     * file name takes precedence over select_img
     * If -1 is given the whole object is read
     *
     */
    int read(const FileName &name, bool readdata=true, int select_img = 0,
             bool apply_geo = false, bool only_apply_shifts = false,
             MDRow * row = NULL, bool mapData = false);
    /** General write function
     * select_img= which slice should I replace
     * overwrite = 0, append slice
     * overwrite = 1 overwrite slice
     */
    void write(FileName name="",
               int select_img=-1,
               bool isStack=false,
               int mode=WRITE_OVERWRITE);
        /** Cast a page of data from type dataType to type Tdest
     *    input pointer  char *
     */

    void castPage2T(char * page, T * ptrDest, DataType datatype, size_t pageSize );

    /** Cast page from T to datatype
     *  input pointer char *
     */
    void castPage2Datatype(T * srcPtr, char * page, DataType datatype, size_t pageSize );

    /** Check file Datatype is same as T type to use mmap.
     */
    bool checkMmapT(DataType datatype);

    /** Write an entire page as datatype
     *
     * A page of datasize_n elements T is cast to datatype and written to fimg
     * The memory for the casted page is allocated and freed internally.
     */
    void writePageAsDatatype(FILE * fimg, DataType datatype, size_t datasize_n );

    /** Swap an entire page
      * input pointer char *
      */
    void swapPage(char * page, size_t pageNrElements, DataType datatype);

    /** Read the raw data
      */
    void readData(FILE* fimg, int select_img, DataType datatype, unsigned long pad);

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
    MultidimArray<T>& operator()();

    const MultidimArray<T>& operator()() const;

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
    T& operator()(int i, int j) const;
    /** Set pixel
     * (direct access) needed by swig
     */
    void setPixel(int i, int j, T v);
    /** Get pixel
     * (direct acces) needed by swig
     */
    T getPixel(int i, int j) const;
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
    T& operator()(int k, int i, int j) const;
    /** Get file name
     *
     * @code
     * std::cout << "Image name = " << I.name() << std::endl;
     * @endcode
     */
    const FileName & name() const;
    /** Get Image dimensions
     */
    void getDimensions(int &Xdim, int &Ydim, int &Zdim, int &Ndim) const;

    long unsigned int getSize() const;

    /** Get Image offset and swap
     */
    void getOffsetAndSwap(unsigned long &_offset, int &_swap) const;

    /* Is there label in the individual header */
    bool individualContainsLabel(MDLabel label) const;

    /* Is there label in the main header */
    bool mainContainsLabel(MDLabel label) const;

    /** Get Rot angle
    *
    * @code
    * std::cout << "First Euler angle " << I.rot() << std::endl;
    * @endcode
    */
    double rot(const long int n = 0) const;

    /** Get Tilt angle
     *
     * @code
     * std::cout << "Second Euler angle " << I.tilt() << std::endl;
     * @endcode
     */
    double tilt(const long int n = 0) const;

    /** Get Psi angle
     *
     * @code
     * std::cout << "Third Euler angle " << I.psi() << std::endl;
     * @endcode
     */
    double psi(const long int n = 0) const;
    /** Get Xoff
     *
     * @code
     * std::cout << "Origin offset in X " << I.Xoff() << std::endl;
     * @endcode
     */
    double Xoff(const long int n = 0) const;
    /** Get Yoff
     *
     * @code
     * std::cout << "Origin offset in Y " << I.Yoff() << std::endl;
     * @endcode
     */
    double Yoff(const long int n = 0) const;
    /** Get Zoff
     *
     * @code
     * std::cout << "Origin offset in Z " << I.Zoff() << std::endl;
     * @endcode
     */
    double Zoff(const long int n = 0) const;
    /** Get Weight
    *
    * @code
    * std::cout << "weight= " << I.weight() << std::endl;
    * @endcode
    */
    double weight(const long int n = 0) const;
    /** Get Flip
    *
    * @code
    * std::cout << "flip= " << flip() << std::endl;
    * @endcode
    */
    bool flip(const long int n = 0) const;

    /** Data type
        *
        * @code
        * std::cout << "datatype= " << dataType() << std::endl;
        * @endcode
        */
    int dataType() const;

    /** Sampling RateX
    *
    * @code
    * std::cout << "sampling= " << samplingRateX() << std::endl;
    * @endcode
    */
    double samplingRateX() const;

    /** Set file name
     */
    void setName(const FileName &_filename);

    /** Set Euler angles in image header
     */
    void setEulerAngles(double rot, double tilt, double psi,
                        long int n = 0);
    /** Get Euler angles from image header
     */
    void getEulerAngles(double &rot, double &tilt, double &psi,
                        long int n = 0);
    /** Set Rotation angle to image */
    void setRot(double rot, long int n = 0);
    /** Set Tilt angle to image */
    void setTilt(double tilt, long int n = 0);
    /** Set Rotation angle to image */
    void setPsi(double psi, long int n = 0);
    /** Set origin offsets in image header
     */
    void setShifts(double xoff, double yoff, double zoff = 0.,
                   long int n = 0);
    /** Get origin offsets from image header
      */
    void getShifts(double &xoff, double &yoff, double &zoff = 0.,
                   long int n = 0);
    /** Set X offset in image header
     */
    void setXoff(double xoff, long int n = 0);
    /** Set Y offset in image header
     */
    void setYoff(double yoff, long int n = 0);
    /** Set Z offset in image header
     */
    void setZoff(double zoff, long int n = 0);
    /** Set flip in image header
     */
    void setFlip(bool flip, long int n = 0);
    /** Set Weight in image header
    */
    void setWeight(double weight, long int n = 0);
    /** Get geometric transformation matrix from 2D-image header
      */
    Matrix2D< double > getTransformationMatrix(bool only_apply_shifts = false,
            long int n = 0);
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
          o << "Unsigned character or byte type";
          break;
      case SChar:
          o << "Signed character (for CCP4)";
          break;
      case UShort:
          o << "Unsigned integer (2-byte)";
          break;
      case Short:
          o << "Signed integer (2-byte)";
          break;
      case UInt:
          o << "Unsigned integer (4-byte)";
          break;
      case Int:
          o << "Signed integer (4-byte)";
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
      o << "Header size  : " << I.offset << std::endl;
      if (I.individualContainsLabel(MDL_WEIGHT))
          o << "Weight  : " << I.weight() << std::endl;
      if (I.individualContainsLabel(MDL_FLIP))
          o << "Flip    : " << I.flip() << std::endl;
      return o;
    }
    /** Sum this object with other file and keep in this object
      */
    void sumWithFile(const FileName &fn);
};

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

// Explicit template instantiation
template class Image<double>;
template class Image<int>;
template class Image<float>;
template class Image< std::complex<double> >;

/// @defgroup ImageFormats Image Formats
/// @ingroup Images
// Functions belonging to this topic are commented in rw*.h
//@}
#endif
