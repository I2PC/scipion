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

#include <typeinfo>
#include <fcntl.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include "funcs.h"
#include "memory.h"
#include "multidim_array.h"
#include "transformations.h"
#include "metadata.h"

// Includes for rwTIFF which cannot be inside it
#include <cstring>
#include "../../external/tiff-3.9.4/libtiff/tiffio.h"

/* Minimum size of a TIFF file to be mapped to a tempfile in case of mapping from
 * image file is required
 */
const size_t tiff_map_min_size = 200e6;

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
    WRITE_OVERWRITE, //forget about the old file and overwrite it
    WRITE_APPEND,    //append and object at the end of a stack, so far can not append stacks
    WRITE_REPLACE,   //replace a particular object by another
    WRITE_READONLY   //only can read the file
} WriteMode;

/* Cast Write mode
 * This enum defines the cast writing behavior
 */
typedef enum
{
    CAST,       //Only cast the data type
    CONVERT,    //Convert the data from one type to another
    ADJUST      //Adjust the histogram to fill the gray level range
} CastWriteMode;

/** Open File struct
 * This struct is used to share the File handlers with Image Collection class
 */
struct ImageFHandler
{
    FILE*     fimg;       // Image File handler
    FILE*     fhed;       // Image File header handler
    TIFF*     tif;        // TIFF Image file hander
    FileName  fileName;   // Image file name
    FileName  headName;   // Header file name
    FileName  ext_name;   // Filename extension
    bool     exist;       // Shows if the file exists. Equal 0 means file does not exist or not stack.
};


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


// Image base class
class ImageBase
{
public:
    virtual void getDimensions(int &Xdim, int &Ydim, int &Zdim, unsigned long &Ndim) const =0;
    virtual void getEulerAngles(double &rot, double &tilt, double &psi,
                                long int n = 0)=0;
    virtual double tilt(const long int n = 0) const =0;
    virtual int read(const FileName &name, bool readdata=true, int select_img = -1,
                     bool apply_geo = false, bool only_apply_shifts = false,
                     MDRow * row = NULL, bool mapData = false)=0;
    virtual void write(const FileName &name="", int select_img=-1, bool isStack=false,
                       int mode=WRITE_OVERWRITE,bool adjust=false)=0;
    virtual void newMappedFile(int Xdim, int Ydim, int Zdim, int Ndim, FileName _filename)=0;
    virtual void clear()=0;
};

#endif /* IMAGE_BASE_H_ */
