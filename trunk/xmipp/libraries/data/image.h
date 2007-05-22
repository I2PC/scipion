/***************************************************************************
 *
 * Authors: Alberto Pascual Montano (pascual@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef IMAGE_H
#define IMAGE_H

#include "funcs.h"
#include "matrix2d.h"
#include "header.h"
#include "geometry.h"

#include <typeinfo>
#include <complex>
#include <vector>

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
#include <sstream>
#else
#include <strstream.h>
#endif

// This forward declarations are needed for defining operators functions that
// use other clases type

template<typename T>
class ImageImagicT;
template<typename T>
class ImageXmippT;

// String that identifies a selline as an Imagic image
static const char* IMAGIC_TAG = "imagic:";

/// @defgroup Images Images

typedef enum
{
    IBYTE = 1, IFLOAT = 2, I16 = 3
}
Image_Type;

/** Basic image class
 * @ingroup Images
 *
 * The image class is a general class which only contains the image itself and a
 * filename for it. It has got a "float" matrix2D as member, and basically all
 * operations between images are based on that class.
 *
 * This class is the usual one when you want to operate images in memory. Images
 * belonging to this class cannot be saved, instead you could use the class
 * ImageXmipp which inherits from this class, so you have the whole funcionality
 * of the class Image plus the possibility of saving at the end.
 *
 * See "Image logical and physical access" for a detailed information about the
 * pixel accessing, and the conventions used in the image definition.
 */
template<typename T>
class ImageT
{
protected:
    FileName fn_img; // name of the image

public:
    matrix2D< T > img; // matrix with the image

public:
    /// @defgroup ImageConstructors Image constructors
    /// @ingroup Images

    /** Empty constructor
     * @ingroup ImageConstructors
     *
     * An empty image with size 0x0 is created.
     *
     * @code
     * Image I;
     * @endcode
     */
    ImageT()
    {
        fn_img = "";
    }

    /** Constructor with size
     * @ingroup ImageConstructors
     *
     * A blank image (0.0 filled) is created with the given size. Pay attention
     * to the dimension order: Y and then X.
     *
     * @code
     * Image I(64,64);
     * @endcode
     */
    ImageT(int Ydim, int Xdim)
    {
        img.resize(Ydim, Xdim);
        fn_img = "";
    }

    /** Constructor using filename
     * @ingroup ImageConstructors
     *
     * The name for the image is assigned but no load is performed. The image
     * content is empty, ie, with size 0x0.
     *
     * @code
     * Image I("g0ta0001");
     * @endcode
     */
    ImageT(FileName _name)
    {
        fn_img = _name;
    }

    /** Copy constructor
     * @ingroup ImageConstructors
     *
     * @code
     * Image I2(I1);
     * @endcode
     */
    ImageT(const ImageT& I)
    {
        img = I.img;
        fn_img = I.fn_img;
    }

    /// @defgroup ImageOperations Some operations
    /// @ingroup Images

    /** Assignment
     * @ingroup ImageOperations
     */
    ImageT& operator=(const ImageT& I)
    {
        if (this != &I)
        {
            fn_img = I.fn_img;
            img = I.img;
        }
        return *this;
    }

    /** Another function for assignment
     * @ingroup ImageOperations
     */
    void assign(const ImageT& I)
    {
        *this = I;
    }

    /** Assignment from matrix
     * @ingroup ImageOperations
     */
    template<typename T1>
    ImageT& operator=(const matrix2D< T1 >& m)
    {
        if (&img != (matrix2D< T >*) &m)
        {
            fn_img = "";
            type_cast(m, img);
        }
        return *this;
    }

    /** Another function for assignment from matrix
     * @ingroup ImageOperations
     */
    template<typename T1>
    void assign(const matrix2D< T1 >& m)
    {
        *this = m;
    }

    /** Rename image
     * @ingroup ImageOperations
     *
     * Give a new name to the image.
     *
     * @code
     * I.rename("new_name");
     * @endcode
     */
    virtual void rename(FileName newName)
    {
        fn_img = newName;
    }

    /** Empty image
     * @ingroup ImageOperations
     *
     * This function clears the image to a 0x0 image without name.
     *
     * @code
     * I.clear();
     * @endcode
     */
    virtual void clear()
    {
        fn_img = "";
        img.clear();
    }

    /** Sets the origin of the image at its center
     * @ingroup ImageOperations
     *
     * The exact formula for the center of an image is defined in the
     * conventions section, this function modify the indexes of the image in
     * such a way that now the logical origin is at its center and the indexes
     * may take negative values. Use startingX, startingY, finishingX and
     * finishingY to setup loops.
     *
     * @code
     * I.move_origin_to_center();
     * @endcode
     *
     */
    void move_origin_to_center()
    {
        img.startingY() = FIRST_XMIPP_INDEX(img.ydim);
        img.startingX() = FIRST_XMIPP_INDEX(img.xdim);
    }

    /** Fill with 0 and move origin to center
     * @ingroup ImageOperations
     *
     * This function resizes the image to the given size, fill it with 0.0 and
     * ten moves the image logical origin to its center. See previous function.
     *
     * @code
     * I.adapt_to_size(64, 64)
     * @endcode
     */
    void adapt_to_size(int Ydim, int Xdim)
    {
        img.init_zeros(Ydim, Xdim);
        move_origin_to_center();
    }

    /// @defgroup ImageAccess Image access
    /// @ingroup Images

    /** 2D Matrix access
     * @ingroup ImageAccess
     *
     * This operator can be used to access the matrix, and the matrix operations
     * defined in matrix2D. In this way we could resize an image just by
     * resizing its associated matrix or we could add two images by adding their
     * matrices.
     *
     * @code
     * I().resize(128, 128);
     * I2() = I1() + I2();
     * @endcode
     */
    virtual matrix2D< T >&  operator()()
    {
        return img;
    }

    virtual const matrix2D< T >&  operator()() const
    {
        return img;
    }

    /** Another function for 2D Matrix access
     * @ingroup ImageAccess
     */
    virtual void get_matrix2D(matrix2D< T >& m)
    {
        m = img;
    }

    /** Another function for 2D Matrix access
     * @ingroup ImageAccess
     */
    virtual void set_matrix2D(const matrix2D< T >& m)
    {
        img = m;
    }

    /** Pixel access
     * @ingroup ImageAccess
     *
     * This operator is used to access a pixel within the image. This is a
     * logical access, so you could access to negative positions if the image
     * has been defined so (see the general explanation for the class).
     *
     * @code
     * cout << "Grey level of pixel (-3,-3) of the image = " << I(-3, -3)
     * << endl;
     *
     * I(-3, -3) = I(-3, -2);
     * @endcode
     */
    T& operator()(int i, int j) const
    {
        return img(i, j);
    }

    /** Another function for pixel access
     * @ingroup ImageAccess
     */
    T get_pixel(int i, int j) const
    {
        return img(i, j);
    }

    /** Another function for pixel access
     * @ingroup ImageAccess
     */
    void set_pixel(int i, int j, T val)
    {
        img(i, j) = val;
    }

    /** Name access
     * @ingroup ImageAccess
     *
     * This function is used to know the name of the image. It cannot be used to
     * assign a new one.
     *
     * @code
     * cout << "Image name: " << I.name() << endl;
     * @endcode
     */
    const FileName name() const
    {
        return fn_img;
    }

    /** Cout << Image
     * @ingroup ImageAccess
     *
     * Shows name and size
     */
    friend ostream& operator<<(ostream& out, const ImageT& I)
    {
        out << "Image Name   : " << I.fn_img << endl
        << "dimensions   : " << I.img.RowNo() << " x " << I.img.ColNo()
        << "  (rows x columns)" << endl;

        return out;
    }

    /** @defgroup ImageIO I/O functions
     * @ingroup Images
     *
     * All these functions work with the image written in raw floats.
     */


    /** Read image from disk, given the image's dimensions
     * @ingroup ImageIO
     *
     * If the image doesn't exist at the given path then an exception is thrown.
     *
     * The reversed flag indicates if the elements in element_size must be read
     * in a reversed way.
     *
     * Elements are supposed to be in the following order:
     * (y,x)= (0,0)(0,1)(0,2), ..., (0,Xdim-1), (1,0), (1,1), ...
     *
     * The element size can be adjusted so that raw images of bytes (IBYTE),
     * unsigned ints woth 16 bits (I16) and floats (IFLOAT) can be read.
     *
     * The headersize parameter can be used to read raw images
     * preceded by a header. This is the normal image format of many
     * image processing packages.
     *
     * @code
     * I.read(65, 65, "g0ta0001.raw");
     * @endcode
     */
    bool read(const FileName& name, float fIform, int Ydim, int Xdim,
              bool reversed = false, Image_Type image_type = IBYTE,
              int header_size = 0)
    {
        FILE* fh;
        clear();

        if ((fh = fopen(name.c_str(), "rb")) == NULL)
            REPORT_ERROR(1501, "Image::read: File " + name + " not found");
        fseek(fh, header_size, SEEK_SET);

        bool ret;
        if ((ret = ImageT<  T>::read(fh, fIform, Ydim, Xdim, reversed,
                                     image_type)))
            fn_img = name;

        fclose(fh);
        return (ret);
    }

    /** Read image from disk using a file pointer
     * @ingroup ImageIO
     *
     * This is the core routine of the previous one.
     */
    bool read(FILE*& fh, float fIform, int Ydim, int Xdim, bool reversed,
              Image_Type image_type)
    {
        img.resize(Ydim, Xdim);

        FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
        switch (image_type)
        {
        case IBYTE:
            unsigned char u;
            FREAD(&u, sizeof(unsigned char), 1, fh, reversed);
            MULTIDIM_ELEM(img, i) = u;
            break;

        case I16:
            unsigned short us;
            FREAD(&us, sizeof(unsigned short), 1, fh, reversed);
            MULTIDIM_ELEM(img, i) = us;
            break;

        case IFLOAT:
            float f;
            FREAD(&f, sizeof(float), 1, fh, reversed);
            MULTIDIM_ELEM(img, i) = f;
            break;
        }
        return (true);
    }

    /** Load an image
     * @ingroup ImageIO
     *
     * LoadImage is a static function that loads an image file depending on what
     *  kind of image type (e.g., Xmipp, Imagic) it is; it returns a pointer to
     * the base Image class.  Caller is responsible for deleting the memory for
     * the object.  Returns NULL on error.
     */
    static ImageT< T >* LoadImage(const FileName& name, bool apply_geo = false)
    {
        ImageT< T >* ret = NULL;
        if (name.find(IMAGIC_TAG) != string::npos)
        {
            ImageImagicT< T >* i = new ImageImagicT< T >();
            if (i->read(name))
                ret = i;
            else
                delete i;
        }
        else // For now, assume anything else is ImageXmipp type.
        {
            ImageXmippT< T >* i = new ImageXmippT< T >();
            if (i->read(name, false, false, apply_geo, false))
                ret = i;
            else
                delete i;
        }
        return (ret);
    }

    /** Write image to disk
     * @ingroup ImageIO
     *
     * If there is any problem in the writing, an exception is thrown. You can
     * give a name to the written volume different from the one used when it was
     * read. From this point the filename of the volume has changed. This is
     * somehow like the "Save as ..." and "Save".
     *
     * Valid types are IBYTE, I16, IFLOAT.
     *
     * @code
     * I.write(); // Save
     *
     * I.write("g0ta0001.raw") // Save as
     * @endcode
     */
    void write(FileName name = "", bool reversed = false,
               Image_Type image_type = IBYTE)
    {
        FILE* fp;
        if (name != "")
            ImageT< T >::rename(name);

        if ((fp = fopen(fn_img.c_str(), "wb")) == NULL)
        {
            REPORT_ERROR(1503, "Image::write: File " + fn_img +
                         " cannot be saved");
        }

        ImageT< T >::write(fp, reversed, image_type);
        fclose(fp);
    }

    /** Write image to disk using a file pointer
     * @ingroup ImageIO
     *
     * This is the core routine of the previous one.
     */
    void write(FILE*& fh, bool reversed, Image_Type image_type)
    {
        if (XSIZE(img) == 0 || YSIZE(img) == 0)
            return;

        double a, b;
        if (image_type != IFLOAT)
        {
            double min_val, max_val;
            (*this)().compute_double_minmax(min_val, max_val);
            if (image_type == IBYTE)
                a = 255;
            else
                a = 65535;

            a /= (max_val - min_val);
            b = min_val;
        }

        FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
        switch (image_type)
        {
        case IBYTE:
            unsigned char u;
            u = (unsigned char) ROUND(a * (MULTIDIM_ELEM(img, i) - b));
            FWRITE(&u, sizeof(unsigned char), 1, fh, reversed);
            break;

        case I16:
            unsigned short us;
            us = (unsigned short) ROUND(a * (MULTIDIM_ELEM(img, i) - b));
            FWRITE(&us, sizeof(unsigned short), 1, fh, reversed);
            break;

        case IFLOAT:
            float f;
            f = (float) MULTIDIM_ELEM(img, i);
            FWRITE(&f, sizeof(float), 1, fh, reversed);
            break;
        }
    }
};

// Specialized function to read images with complex numbers in them
// Inlined to avoid multiple definitions
template<>
inline bool ImageT< complex< double> >::read(FILE*& fh, float fIform,
        int Ydim, int Xdim, bool reversed, Image_Type image_type)
{
    img.resize(Ydim, Xdim);
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
    {
        float a, b;

        // read real part of a complex number
        FREAD(&a, sizeof(float), 1, fh, reversed);

        // read imaginary part of a complex number
        FREAD(&b, sizeof(float), 1, fh, reversed);

        // Assign the number
        complex< double > c(a, b);
        MULTIDIM_ELEM(img, i) = c;
    }
    return (true);
}

// Specialized function to write images with complex numbers in them
template<>
inline void ImageT< complex< double> >::write(FILE*& fh, bool reversed,
        Image_Type image_type)
{
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
    {
        float a, b;
        a = (float)(MULTIDIM_ELEM(img, i)).real();
        b = (float)(MULTIDIM_ELEM(img, i)).imag();

        FWRITE(&a, sizeof(float), 1, fh, reversed);
        FWRITE(&b, sizeof(float), 1, fh, reversed);
    }
}

/// @defgroup ImageMacros Speed up macros
/// @ingroup Images

/** Matrix access
 * @ingroup ImageMacros
 *
 * This macro does the same as the normal matrix access but in a faster way as
 * no function call is generated. This macro can be used with any of the derived
 * classes from the Image class.
 *
 * @code
 * IMGMATRIX(I).resize(128, 128);
 *
 * IMGMATRIX(I2) = IMGMATRIX(I1) + IMGMATRIX(I2);
 * @endcode
 */
#define IMGMATRIX(I) ((I).img)

/** Array access
 * @ingroup ImageMacros
 *
 * This macro allows you to access to the bidimensional array behind the image
 * (double*).
 */
#define IMGARRAY(I) MAT_ARRAY((I).img)

/** Pixel access
 * @ingroup ImageMacros
 *
 * This macro does the same as the normal pixel access (remember, logical
 * access) but in a faster way as no function call is generated. This macro can
 * be used with any of the derived classes from the Image class, except for the
 * Oversampled Images (see ImageOver) as their logical positions are related in
 * a special way to the physical ones.
 *
 * @code
 * cout << "Grey level of pixel (-3, -3) of the image = "
 *      << IMGPIXEL(I, -3, -3) << endl;
 *
 * IMGPIXEL(I, -3, -3) = IMGPIXEL(I, -3, -2);
 * @endcode
 */
#define IMGPIXEL(I, i, j) MAT_ELEM(((I).img), (i), (j))

/** Physical pixel access
 * @ingroup ImageMacros
 *
 * The physical pixel access gives you access to a pixel by its physical
 * position and not by its logical one. This access shouldn't be used as a
 * custom, use instead the logical access, but there might be cases in which
 * this access might be interesting.
 *
 * Physical positions start at index 0 in C.
 *
 * This macro can be used by any of the derived classes from the Image class.
 *
 * @code
 * cout << "This is the first pixel stored in the image " <<
 *      DIRECT_IMGPIXEL(I, 0, 0) << endl;
 * @endcode
 */
#define DIRECT_IMGPIXEL(I, i, j) DIRECT_MAT_ELEM(((I).img), (i), (j))


/** Xmipp 2D Images
 * @ingroup Images
 *
 * The Xmipp image is a normal image (inherited from Image class) plus a Spider
 * header. This is the appropiate class to work with images in memory which we
 * want to save later.
 *
 * The data in the header is not directly accesible from the programmer and must
 *  be set through object functions, in this way the coherence in the header is
 * assured. See File Formats for more information about the Spider format.
 *
 * In principle, the image starts at pixel (0,0) but this can be modified for
 * any other logical access. See class Image for more information.
 *
 */
template<typename T>
class ImageXmippT: public ImageT< T >
{
protected:
    headerXmipp header; // Declares a header

public:
    /** Empty constructor
     *
     * Creates an empty (0x0) image with no information in the header.
     *
     * @code
     * ImageXmipp IX;
     * @endcode
     */
    ImageXmippT(): ImageT< T >()
    {
        if (typeid(T) == typeid(complex< double >))
            // Sets header of type Image_Fourier(Complex)
            header.headerType() = headerXmipp::IMG_FOURIER;
        else if (typeid(T) == typeid(double))
            // Sets header of type Image_XMipp
            header.headerType() = headerXmipp::IMG_XMIPP;

        // Sets header Slices
        header.Slices() = 1;
    }

    /** Constructor with size
     *
     * Creates a 0.0 filled image of size Ydim x Xdim.
     *
     * @code
     * ImageXmipp IX(64, 64);
     * @endcode
     */
    ImageXmippT(int Ydim, int Xdim) : ImageT< T >(Ydim, Xdim)
    {
        if (typeid(T) == typeid(complex< double >))
            // Sets header of type Image_Fourier(Complex)
            header.headerType() = headerXmipp::IMG_FOURIER;
        else if (typeid(T) == typeid(double))
            // Sets header of type Image_XMipp
            header.headerType() = headerXmipp::IMG_XMIPP;

        header.set_dimension(Ydim, Xdim); // Sets header dimensions
        header.Slices() = 1; // Sets header Slices
        header.set_header(); // Initialize header
        header.set_time(); // Set time and date
        header.set_date();
    }

    /** Constructor with filename, read from disk
     *
     * The filename given must exist, then the file is loaded in the ImageXmipp
     * class structure. You have loaded the image at the declaration time.
     *
     * @code
     * ImageXmipp IX("g1ta0001.spd");
     * @endcode
     */
    ImageXmippT(FileName _name, bool apply_geo = false) : ImageT< T >(_name)
    {
        if (typeid(T) == typeid(complex< double >))
            // Sets header of type Image_Fourier (Complex)
            header.headerType() = headerXmipp::IMG_FOURIER;
        else if (typeid(T) == typeid(double))
            // Sets header of type Image_XMipp
            header.headerType() = headerXmipp::IMG_XMIPP;

        // Read image from file
        read(_name, false, false, apply_geo, false);
    }

    /** Copy constructor
     *
     * @code
     * ImageXmipp IX2(IX1);
     * @endcode
     */
    ImageXmippT(const ImageXmippT< T >& I) : ImageT< T >(I)
    {
        header = I.header;
    }

    /** Empty image
     *
     * All information is cleared.
     *
     * @code
     * IX.clear();
     * @endcode
     */
    void clear()
    {
        clear_header();
        ImageT< T >::clear();
    }

    /** Show the header information of a Xmipp image
     *
     * @code
     * cout << IX;
     * @endcode
     */
    friend ostream& operator<<(ostream& out, const ImageXmippT< T >& I)
    {
        out << (ImageT< T >&) I << I.header;
        return out;
    }

    /** Assignment from another Xmipp image
     *
     * @code
     * ImageXmipp IX1, IX2;
     * IX2 = IX1;
     * @endcode
     */
    ImageXmippT<T>& operator= (const ImageXmippT<T> &op1)
    {
        if (&op1 != this)
        {
            this->ImageT<T>::operator = (op1);
            header = op1.header;
        }
        return *this;
    }

    /** Another function for assignment from another Xmipp image
     */
    void assign(const ImageXmippT< T >& op1)
    {
        *this = op1;
    }

    /** Assignment from a generic image
     *
     * The Euler angles are set to 0.
     *
     * @code
     * Image I;
     * ImageXmipp IX;
     * IX = I;
     * @endcode
     */
    ImageXmippT< T >& operator=(const ImageT< T >& op1)
    {
        if (this != &op1)
        {
            this->ImageT< T >::operator=(op1);
            clear_header();
            adjust_header();
        }

        return *this;
    }

    /** Another function for assignment from a generic image
     */
    void assign(const ImageT< T >& op1)
    {
        *this = op1;
    }

    /** Assignment from a matrix
     *
     * The Euler angles are set to 0 and the image filename is set to "".
     *
     * @code
     * matrix2D< int > m;
     * ImageXmipp IX;
     * IX = m;
     * @endcode
     */
    template<typename T1>
    ImageXmippT< T >& operator=(const matrix2D< T1 >& op1)
    {
        this->ImageT< T >::operator=(op1);
        clear_header();
        adjust_header();
        return *this;
    }

    /** Another function for assignment from a matrix
     */
    template<typename T1>
    void assign(const matrix2D< T1 >& op1)
    {
        *this = op1;
    }

    /** Assignment from any kind of image
     *
     * @code
     * ImageOver IO;
     * ImageXmipp IX;
     * IX.assign_from(&IO);
     * @endcode
     */
    void assign_from(ImageT< T >* v)
    {
        *this = *v;
    }

    /** Read Xmipp image from disk
     *
     * If the image doesn't exist at the given path then an exception is thrown.
     *
     * If apply_geo and/or only_apply_shifts, the image will be mirrored if the
     * img.flip() flag is set to 1.
     *
     * If only_apply_shifts is TRUE, then only the shifts of the header are
     * applied (even if apply_geo is false). If apply_geo is TRUE and
     * only_apply_shifts is false, then shifts and in-plane rotation are
     * applied.
     *
     * If skip_type_check is false, then the routine automatically checks the
     * endianness of the file, and reversed should be set to false.
     *
     * If skip_type_check is TRUE, then the programmer must know whether the
     * endian is reversed or not, and provide this information in reversed.
     *
     * @code
     * IX.read("g1ta0002.spd");
     * @endcode
     */
    virtual bool read(const FileName& name, bool skip_type_check = false,
                      bool reversed = false, bool apply_geo = false,
                      bool only_apply_shifts = false)
    {
        FILE* fp;
        bool ret;

        ImageXmippT< T >::rename(name);
        if ((fp = fopen(ImageT< T >::fn_img.c_str(), "rb")) == NULL)
            REPORT_ERROR(1501, (string) "ImageXmipp::read: File " +
                         ImageT<T>::fn_img + " not found");

        // Read header
        if (!header.read(fp, skip_type_check, reversed))
            REPORT_ERROR(1502, "ImageXmipp::read: File " + ImageT< T >::fn_img +
                         " is not a valid Xmipp file");

        // Read whole image and close file
        if ((ret = ImageT< T >::read(fp, header.fIform(), header.iYdim(),
                                     header.iXdim(), header.reversed(), IFLOAT)))
        {
            if (apply_geo || only_apply_shifts)
            {
                // Apply the geometric transformations in the header to the
                // loaded image. Transform image without wrapping, set new
                // values to first element in the matrix
                matrix2D< double > A =
                    ImageXmippT< T >::get_transformation_matrix(
                        only_apply_shifts);

                if (!A.IsIdent())
                    ImageT< T >::img.self_apply_geom_Bspline(A, 3, IS_INV,
                            WRAP);
            }
//scale value in header is not reliable, do not use it
#undef NEVERDEFINED
#ifdef NEVERDEFINED
            // scale if necessary
            // TODO check this with Carlos
            if ((header.Scale() != 0.) && (header.Scale() != 1.))
            {
                header.set_dimension(header.Ydim() * header.Scale(),
                                     header.Xdim() * header.Scale());

                ImageT< T >::img.self_scale_to_size_Bspline(3, header.iYdim(),
                        header.iXdim());
            }

            header.set_header(); // Set header in a Xmipp consistent state
#endif
        }

        fclose(fp);
        return (ret);
    }

    /** Write Xmipp image to disk
     *
     * If there is any problem in the writing, an exception is thrown. You can
     * give a name to the written image different from the one used when it was
     * read. From this point the filename of the image has changed. This is
     * somehow like the "Save as ..." and "Save".
     *
     * @code
     * IX.write(); // Save
     *
     * IX.write("g1ta0002.spd"); // Save as
     * @endcode
     */
    virtual void write(const FileName& name = "", bool force_reversed = false)
    {
        FILE* fp;
        if (name != "")
            ImageXmippT< T >::rename(name);

        if ((fp = fopen(ImageT<T>::fn_img.c_str(), "wb")) == NULL)
            REPORT_ERROR(1503, (string) "ImageXmipp::write: File " + name +
                         " cannot be written");

        adjust_header();
        bool reversed = (force_reversed) ? !header.reversed() : header.reversed();
        header.write(fp, force_reversed);
        ImageT< T >::write(fp, reversed, IFLOAT);
        fclose(fp);
    }

    /** Resets header
     */
    void clear_header()
    {
        header.clear();
    }

    /** Adjust header
     *
     * Force header to have the dimensions of the image, time, date updated
     */
    void adjust_header()
    {
        if (typeid(T) == typeid(complex< double >))
            // Sets header of type Image_Fourier (Complex)
            header.headerType() = headerXmipp::IMG_FOURIER;
        else if (typeid(T) == typeid(double))
            // Sets header of type Image_XMipp
            header.headerType() = headerXmipp::IMG_XMIPP;

        // Sets header dimensions
        header.set_dimension(YSIZE(ImageT<T>::img), XSIZE(ImageT<T>::img));

        header.Slices() = 1; // Sets header Slices
        header.set_time(); // Set time and date
        header.set_date();
        header.set_title(ImageT< T >::fn_img); // Set title
        header.set_header(); // Initialize header
    }

    /** Change Filename
     *
     * @code
     * IX.rename("newName.spd");
     * @endcode
     */
    void rename(const FileName& newName)
    {
        ImageT< T >::rename(newName);
        header.set_title(newName);
    }

    /** Reversed status
     *
     * This is used for the little/big endian process.
     */
    bool reversed() const
    {
        return header.reversed();
    }

    /** Get geometric transformation matrix from 2D-image header
     */
    matrix2D< double > get_transformation_matrix(bool only_apply_shifts = false)
    {
        matrix2D< double > A(3, 3);
        double psi = realWRAP(header.Psi(), -180, 180);
        double theta = realWRAP(header.Theta(), -180, 180);

        A.init_identity();

        // This is to be compatible with old Xmipp images, that store image
        // translations in another position of the header
        ImageXmippT<T>::check_oldxmipp_header();

        if (only_apply_shifts)
        {
            Euler_angles2matrix(0., 0., 0., A);
            A(0, 2) = -header.fXoff();
            A(1, 2) = -header.fYoff();
        }
        else
        {
            if (theta == 0.)
            {
                // For untilted images: apply Euler matrix
                Euler_angles2matrix(header.Phi(), 0., header.Psi(), A);
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
            A(0, 2) = -header.fXoff();
            A(1, 2) = -header.fYoff();
        }

        // Also for only_apply_shifts: mirror if necessary!
        if (header.Flip() == 1)
        {
            A(0, 0) = -A(0, 0);
            A(0, 1) = -A(0, 1);
        }

        return A;
    }

    /** Check OldXmipp header location for image translation offsets
     *
     * If image appears to be in OldXmipp format, copy offsets to NewXmipp
     * header location
     */
    void check_oldxmipp_header()
    {
        if (header.fXoff() == 0. && header.fYoff() == 0. && header.fZoff() == 0.)
        {
            // Might be oldXmipp image header
            float ang = header.old_rot();
            matrix2D< double > mat = header.fGeo_matrix();
            if (ABS(mat(0, 0) - COSD(ang)) < XMIPP_EQUAL_ACCURACY &&
                ABS(mat(1, 1) - COSD(ang)) < XMIPP_EQUAL_ACCURACY &&
                ABS(mat(0, 1) - SIND(ang)) < XMIPP_EQUAL_ACCURACY &&
                ABS(mat(1, 0) + SIND(ang)) < XMIPP_EQUAL_ACCURACY &&
                mat(2, 0) != 0. && mat(2, 1) != 0.)
            {
                // This indeed seems to be an OldXmipp style header with
                // non-zero offsets
                cerr << "WARNING%% Copying shifts from old to new headers: "
                << (*this).name() << endl;

                header.fXoff() = -(float) mat(2, 0);
                header.fYoff() = -(float) mat(2, 1);

                if (XSIZE(ImageT< T >::img) % 2 == 0)
                    header.fXoff() += 0.5;
                if (YSIZE(ImageT< T >::img) % 2 == 0)
                    header.fYoff() += 0.5;
            }
        }
    }

    /** Set origin offsets
     *
     * @code
     * IX.set_originOffsets(1.51, -3.43);
     * @endcode
     */
    void set_originOffsets(float _Xoff, float _Yoff)
    {
        header.set_originOffsets(_Xoff, _Yoff);
    }

    /** Get origin offsets
     *
     * @code
     * IX.get_originOffsets(Xoff, Yoff);
     * @endcode
     */
    void get_originOffsets(float& _Xoff, float& _Yoff) const
    {
        header.get_originOffsets(_Xoff, _Yoff);
    }

    /** Set Xoff
     *
     * @code
     * IX.Xoff() = 3.50;
     * @endcode
     */
    float& Xoff()
    {
        return header.fXoff();
    }

    /** Another function for set Xoff
     */
    void set_Xoff(float& _Xoff)
    {
        header.fXoff() = _Xoff;
    }

    /** Get Xoff
     *
     * @code
     * cout << "Origin Offset in X-direction: " << IX.Xoff() << endl;
     * @endcode
     */
    float Xoff() const
    {
        return header.fXoff();
    }

    /** Set Yoff
     *
     * @code
     * IX.Yoff() = 3.50;
     * @endcode
     */
    float& Yoff()
    {
        return header.fYoff();
    }

    /** Another function for set Yoff
     */
    void set_Yoff(float& _Yoff)
    {
        header.fYoff() = _Yoff;
    }

    /** Get Yoff
     *
     * @code
     * cout << "Origin Offset in Y-direction: " << IX.Yoff() << endl;
     * @endcode
     */
    float Yoff() const
    {
        return header.fYoff();
    }

    /** Set weight
     *
     * @code
     * IX.weight() = 3.50;
     * @endcode
     */
    float& weight()
    {
        return header.Weight();
    }

    /** Another function for set weight
     */
    void set_weight(float& _Weight)
    {
        header.Weight() = _Weight;
    }

    /** Get weight
     *
     * @code
     * cout << "Weight: " << IX.weight() << endl;
     * @endcode
     */
    float  weight() const
    {
        return header.Weight();
    }

    /** Set flip
     *
     * @code
     * IX.flip() = 1; // flip image
     *
     * IX.flip() = 0; // do NOT flip image
     * @endcode
     */
    float& flip()
    {
        return header.Flip();
    }

    /** Another function for set flip
     */
    void set_flip(float&  _Flip)
    {
        header.Flip() = _Flip;
    }

    /** Get flip
     *
     * @code
     * cout << "Flip: " << IX.flip() << endl;
     * @endcode
     */
    float flip() const
    {
        return header.Flip();
    }

    /** Set angles. (default position)
     *
     * @code
     * IX.set_eulerAngles(30, -10, 350);
     * @endcode
     */
    template<typename T1>
    void set_eulerAngles(T1 _Phi, T1 _Theta, T1 _Psi)
    {
        header.set_eulerAngles((float) _Phi, (float) _Theta, (float) _Psi);
    }

    /** Set angles
     *
     * Spider has three places in which Euler angles can be stored. First
     * alternative position
     *
     * @code
     * IX.set_eulerAngles1(30, -10, 350);
     * @endcode
     */
    template<typename T1>
    void set_eulerAngles1(T1 _Phi1, T1 _Theta1, T1 _Psi1)
    {
        header.set_eulerAngles1((float) _Phi1, (float) _Theta1, (float) _Psi1);
    }

    /** Set angles
     *
     * Spider has three places in which Euler angles can be stored. Second
     * alternative position
     *
     * @code
     * IX.set_eulerAngles2(30, -10, 350);
     * @endcode
     */
    template<typename T1>
    void set_eulerAngles2(T1 _Phi2, T1 _Theta2, T1 _Psi2)
    {
        header.set_eulerAngles2((float) _Phi2, (float) _Theta2, (float) _Psi2);
    }

    /** Clears fFlag flag.
     *
     * The number of triads of Euler angles stored in the header (up to three)
     * is stored here. set_eulerAngles2 makes fFlag=2, set_eulerAngles1 makes
     * fFlag=max(fFlag,1), set_eulerAngles does not change fFlag
     *
     * @code
     * IX.clear_fFlag_flag();
     * @endcode
     */
    void clear_fFlag_flag()
    {
        header.clear_fFlag_flag();
    }

    /** Get angles
     *
     * @code
     * IX.get_eulerAngles(phi, theta, psi);
     * @endcode
     */
    template<typename T1>
    void get_eulerAngles(T1& _Phi, T1& _Theta, T1& _Psi) const
    {
        header.get_eulerAngles(_Phi, _Theta, _Psi);
    }

    /** Get angles
     *
     * Spider has three places in which Euler angles can be stored. First
     * alternative position
     *
     * @code
     * IX.get_eulerAngles1(phi, theta, psi);
     * @encode
     */
    template<typename T1>
    void get_eulerAngles1(T1& _Phi1, T1& _Theta1, T1& _Psi1) const
    {
        header.get_eulerAngles1(_Phi1, _Theta1, _Psi1);
    }

    /** Get angles
     *
     * Spider has three places in which Euler angles can be stored. Second
     * alternative position
     *
     * @code
     * IX.get_eulerAngles2(phi, theta, psi);
     * @endcode
     */
    template<typename T1>
    void get_eulerAngles2(float& _Phi2, float& _Theta2, float& _Psi2) const
    {
        header.get_eulerAngles2(_Phi2, _Theta2, _Psi2);
    }

    /** How many Euler angles are set?
     */
    float Is_flag_set(void)
    {
        return(header.Is_flag_set());
    }

    /** Set old rotational angle
     *
     * @code
     * IX.old_rot() = 30;
     * @endcode
     */
    float& old_rot()
    {
        return header.old_rot();
    }

    /** Another function for set old rotational angle
     */
    void set_old_rot(float& _old_rot)
    {
        header.old_rot() = _old_rot;
    }

    /** Get old rotational angle
     *
     * @code
     * cout << "First Euler angle " << IX.old_rot() << endl;
     * @endcode
     */
    float old_rot() const
    {
        return header.old_rot();
    }

    /** Set Phi
     *
     * @code
     * IX.Phi() = 30;
     * @endcode
     */
    float& Phi()
    {
        return header.Phi();
    }

    /** Another function for set Phi
     */
    void set_Phi(float& _Phi)
    {
        header.Phi() = _Phi;
    }

    /** Get Ph
     *
     * @code
     * cout << "First Euler angle " << IX.Phi() << endl;
     * @endcode
     */
    float Phi() const
    {
        return header.Phi();
    }

    /** Set Rotational angle
     *
     * @code
     * IX.rot() = 30;
     * @endcode
     */
    float& rot()
    {
        return header.Phi();
    }

    /** Another function for set Rotational angle
     */
    void set_rot(float& _rot)
    {
        header.Phi() = _rot;
    }

    /** Get rot
     *
     * @code
     * cout << "First Euler angle " << IX.rot() << endl;
     * @endcode
     */
    float rot() const
    {
        return header.Phi();
    }

    /** Set Phi1. First alternative phi angle
     *
     * @code
     * IX.Phi1() = 30;
     * @endcode
     */
    float& Phi1()
    {
        return header.Phi1();
    }

    /** Another function for set Phi1
     */
    void set_Phi1(float& _Phi1)
    {
        header.Phi1() = _Phi1;
    }

    /** Get Phi1. First alternative phi angle
     *
     * @code
     * cout << "First Euler angle " << IX.Phi1() << endl;
     * @endcode
     */
    float  Phi1() const
    {
        return header.Phi1();
    }

    /** Set 1st Rotational angle. First alternative phi angle
     *
     * @code
     * IX.rot1() = 30;
     * @endcode
     */
    float& rot1()
    {
        return header.Phi1();
    }

    /** Another function for set 1st Rotational angle
     */
    void set_rot1(float& _rot1)
    {
        header.Phi1() = _rot1;
    }

    /** Get rot1. First alternative phi angle
     *
     * @code
     * cout << "First Euler angle " << IX.rot1() << endl;
     * @endcode
     */
    float rot1() const
    {
        return header.Phi1();
    }

    /** Set Phi2. First alternative phi angle
     *
     * @code
     * IX.Phi()=30;
     * @endcode
     */
    float& Phi2()
    {
        return header.Phi2();
    }

    /** Another function for set Phi2
     */
    void set_Phi2(float& _Phi2)
    {
        header.Phi2() = _Phi2;
    }

    /** Get Phi2. Second alternative phi angle
     *
     * @code
     * cout << "First Euler angle " << IX.Phi2() << endl;
     * @endcode
     */
    float Phi2() const
    {
        return header.Phi2();
    }

    /** Set 2sd Rotational angle. First alternative phi angle
     *
     * @endcode
     * IX.rot() = 30;
     * @endcode
     */
    float& rot2()
    {
        return header.Phi2();
    }

    /** Another function for set rot2
     */
    void set_rot2(float& _rot2)
    {
        header.Phi2() = _rot2;
    }

    /** Get rot2. Second alternative phi angle
     *
     * @code
     * cout << "First Euler angle " << IX.rot2() << endl;
     * @endcode
     */
    float rot2() const
    {
        return header.Phi2();
    }

    /** Set Theta
     *
     * @code
     * IX.Theta() = -10;
     * @endcode
     */
    float& Theta()
    {
        return header.Theta();
    }

    /** Another function for set Theta
     */
    void set_Theta(float& _Theta)
    {
        header.Theta() = _Theta;
    }

    /** Get Theta
     *
     * @code
     * cout << "Second Euler angle " << IX.Theta() << endl;
     * @endcode
     */
    float Theta() const
    {
        return header.Theta();
    }

    /** Set Tilting angle
     *
     * @code
     * IX.tilt() = -10;
     * @endcode
     */
    float& tilt()
    {
        return header.Theta();
    }

    /** Another function for set Tilting angle
     */
    void set_tilt(float& _tilt)
    {
        header.Theta() = _tilt;
    }

    /** Get Tilting angle
     *
     * @code
     * cout << "Second Euler angle " << IX.tilt() << endl;
     * @endcode
     */
    float tilt() const
    {
        return header.Theta();
    }

    /** Set Theta1
     *
     * @code
     * IX.Theta1() = -10;
     * @endcode
     */
    float& Theta1()
    {
        return header.Theta1();
    }

    /** Another function for set Theta1.*/
    void set_Theta1(float& _Theta1)
    {
        header.Theta1() = _Theta1;
    }

    /** Get Theta1
     *
     * @code
     * cout << "Second Euler angle " << IX.Theta1() << endl;
     * @endcode
     */
    float Theta1() const
    {
        return header.Theta1();
    }

    /** Set 1st Tilting angle
     *
     * @code
     * IX.tilt1() = -10;
     * @endcode
     */
    float& tilt1()
    {
        return header.Theta1();
    }

    /** Another function for set 1st Tilting angle
     */
    void set_tilt1(float& _tilt1)
    {
        header.Theta1() = _tilt1;
    }

    /** Get 1st Tilting angle
     *
     * @code
     * cout << "Second Euler angle " << IX.tilt1() << endl;
     * @endcode
     */
    float tilt1() const
    {
        return header.Theta1();
    }

    /** Set Theta2
     *
     * @code
     * IX.Theta2() = -10;
     * @endcode
     */
    float& Theta2()
    {
        return header.Theta2();
    }

    /** Another function for set Theta2.*/
    void set_Theta2(float& _Theta2)
    {
        header.Theta2() = _Theta2;
    }

    /** Get Theta2
     *
     * @code
     * cout << "Second Euler angle " << IX.Theta2() << endl;
     * @endcode
     */
    float Theta2() const
    {
        return header.Theta2();
    }

    /** Set 2sd Tilting angle
     *
     * @code
     * IX.tilt2() = -10;
     * @endcode
     */
    float& tilt2()
    {
        return header.Theta2();
    }

    /** Another function for set 2sd Tilting angle
     */
    void set_tilt2(float& _tilt2)
    {
        header.Theta2() = _tilt2;
    }

    /** Get 2sd Tilting angle
     *
     * @code
     * cout << "Second Euler angle " << IX.tilt2() << endl;
     * @endcode
     */
    float tilt2() const
    {
        return header.Theta2();
    }

    /** Set Psi
     *
     * @code
     * IX.Psi() = 350;
     * @endcode
     */
    float& Psi()
    {
        return header.Psi();
    }

    /** Another function for set Psi
     */
    void set_Psi(float& _Psi)
    {
        header.Psi() = _Psi;
    }

    /** Get Psi
     *
     * @code
     * cout << "Third Euler angle " << IX.psi() << endl;
     * @endcode
     */
    float Psi() const
    {
        return header.Psi();
    }

    /** Set psi
     *
     * @code
     * IX.psi() = 350;
     * @endcode
     */
    float& psi()
    {
        return header.Psi();
    }

    /** Another function for set Phi
     */
    void set_psi(float& _psi)
    {
        header.Psi() = _psi;
    }

    /** Get psi
     *
     * @code
     * cout << "Third Euler angle " << IX.psi() << endl;
     * @endcode
     */
    float psi() const
    {
        return header.Psi();
    }

    /** Set Psi1
     *
     * @code
     * IX.Psi1() = 350;
     * @endcode
     */
    float& Psi1()
    {
        return header.Psi1();
    }

    /** Another function for set Psi1
     */
    void set_Psi1(float & _Psi1)
    {
        header.Psi1() = _Psi1;
    }

    /** Get Psi1
     *
     * @code
     * cout << "Third Euler angle " << IX.psi1() << endl;
     * @endcode
     */
    float Psi1() const
    {
        return header.Psi1();
    }

    /** Set psi1
     *
     * @code
     * IX.psi1() = 350;
     * @endcode
     */
    float& psi1()
    {
        return header.Psi1();
    }

    /** Another function for set psi1
     */
    void set_psi1(float& _psi1)
    {
        header.Psi1() = _psi1;
    }

    /** Get psi1
     *
     * @code
     * cout << "Third Euler angle " << IX.psi1() << endl;
     * @endcode
     */
    float psi1() const
    {
        return header.Psi1();
    }

    /** Set Psi2
     *
     * @code
     * IX.Psi2() = 350;
     * @endcode
     */
    float& Psi2()
    {
        return header.Psi2();
    }

    /** Another function for set Psi2
     */
    void set_Psi2(float& _Psi2)
    {
        header.Psi2() = _Psi2;
    }

    /** Get Psi2
     *
     * @code
     * cout << "Third Euler angle " << IX.psi2() << endl;
     * @endcode
     */
    float Psi2() const
    {
        return header.Psi2();
    }

    /** Set psi2
     *
     * @code
     * IX.psi2() = 350;
     * @endcode
     */
    float& psi2()
    {
        return header.Psi2();
    }

    /** Another function for set psi2
     */
    void set_psi2(float& _psi2)
    {
        header.Psi2() = _psi2;
    }

    /** Get psi2
     *
     * @code
     * cout << "Third Euler angle " << IX.psi2() << endl;
     * @endcode
     */
    float psi2() const
    {
        return header.Psi2();
    }
};

/** True if the given volume is an Xmipp image
 */
int Is_ImageXmipp(const FileName &fn, bool skip_type_check = false,
                  bool reversed = false);

/** True if the given volume is a Fourier Xmipp image.
 */
int Is_FourierImageXmipp(const FileName& fn, bool skip_type_check = false,
                         bool reversed = false);

/** Get size of an image
 *
 * It returns -1 if the given image is not an Xmipp file
 */
void GetXmippImageSize(const FileName& fn, int& Ydim, int& Xdim);


typedef ImageT< double > Image;
typedef ImageXmippT< double > ImageXmipp;
typedef ImageT< complex< double > > FourierImage;
typedef ImageXmippT< complex< double > > FourierImageXmipp;

/** Converts a Xmipp Image (real numbers) into  a Fourier Xmipp Image
 * (complex numbers)
 */
void ImageXmipp_to_FourierImageXmipp(ImageXmipp& I, FourierImageXmipp& F);

/** Converts a Fourier Xmipp Image  into a Xmipp Image
 *
 * A simple copy of real parts is done
 */
void FourierImageXmipp_to_ImageXmipp(FourierImageXmipp& F, ImageXmipp& I);


/** Oversampled images
 *
 * The oversampled images are images which allow a more accurate treatment of
 * information by oversampling all pixels. The idea of this class is to have an
 * image with a logical size smaller than its physical one, for this reason you
 * could use non integer logical positions to access to different samples of the
 * "same" pixel. Let's set an example, blob footprints are of size 5x5,
 * for instance, with the center at physical position 3x3. It's convenient to
 * have this footprint defined between -2 and 2 with the origin at the center of
 * the image. But it's also convenient to have more than 25 values of the
 * footprint, this could be done by sampling finer each pixel, say 51 samples
 * per pixel, this would result in an image of 255x255 samples (with the center
 * at [127][127]). But we don't to access to the pixel at position (-120,-120)
 * but to the one at logical position (-2.35,-2.35) which we know is close to
 * the border. The oversampled image class allows you this kind of image access.
 * You have to say from where to where your image is defined for each axis (in
 * this case -2...2 for both) and which are the sampling rates for both axes
 * (51 for both). This is the initialisation process. From now on you can work
 * with your image in the way formerly described.
 *
 * Pay attention to two points:
 *
 * * The sampling rate must be an odd number, this is because so the logical
 * oversampled pixels (-2,-2), (-2, -1) ... are exactly on one cell of the
 * underlying 2D matrix.
 *
 * * As the oversampled pixels are centered with respect to the "superpixels"
 * defined by the 2D matrix, the extent of the oversampled image is a little
 * larger than from (-2,-2) to (2,2), ie, from (-2.5,-2.5) to (2.5,2.5)
 *
 * Oversampled images are normal images except for pixel access which can be
 * done in two fashions: either using the normal Image procedures (in this case,
 * you are restricted to the macro IMGPIXEL), or using the fractional indexes
 * defined for this class. Have a look at the following example of how to
 * compute the blob footprint with this kind of images. Pay special attention to
 * the pixel access at the end of the loop. Notice that x and y moves at values
 * between -2 and 2, which keep logical meaning of what we are doing while all
 * the burden of computing where to put these values at the image is done by the
 * library.
 *
 * @code
 * void footprint_blob
 * (
 *     ImageOver& blobprint, // blob foot_print table
 *     const struct blobtype& blob, // blob description
 *     int istep, // number of foot-print samples per one sample
 *                // on projection plane in u,v directions
 *     int   normalise) // set to 1 if you want normalise. Usually
 *                      // you may omit it and no normalisation is performed
 * {
 *     // Resize output image and redefine the origin of it
 *     int footmax = CEIL(blob.radius);
 *     blobprint.init(-footmax, footmax, istep, -footmax, footmax, istep);
 *
 *     // Run for indexes in the Image class
 *     for (int i = STARTINGY(blobprint()); i <= FINISHINGY(blobprint()); i++)
 *         for (int j = STARTINGX(blobprint()); j <= FINISHINGX(blobprint());
 *              j++)
 *         {
 *             // Compute oversampled index, and then its value
 *             double vi, ui;
 *             IMG2OVER(blobprint, i, j, vi, ui);
 *
 *             double r = sqrt(vi*vi + ui*ui);
 *             IMGPIXEL(blobprint, i, j) = blob_val(r, blob);
 *         }
 *
 *     // Adjust the footprint structure
 *     if (normalise)
 *         blobprint() = blobprint() / blobprint().sum();
 * }
 * @endcode
 *
 * Note: for this class the X axis has been renamed as U, and Y as V.
 */
class ImageOver: public Image
{
public:
    int uistep, vistep; // number of samples per
    // one sample on normal image (50)
    // in u,v directions

    int overumax, overvmax; // table borders in normal units (-2,2)
    int overumin, overvmin; // table borders in normal units (-2,2)
    // They should be an integer number

public:
    /** Prepare the Oversampled Image for work
     *
     * This function defines the oversampled image, it's very important to call
     * it before starting to work with the image. If the sampling rate in any of
     * the directions is an even number, then it is substituted by the next odd
     * number (+1). This is done in order to force the logical oversampled
     * pixels to be exactly in one cell of the underlying 2D matrix.
     *
     * @code
     * IO.init(-2, 2, 51, -2, 2, 51);
     * @endcode
     */
    void init(int _vmin, int _vmax, int _vistep,
              int _umin, int _umax, int _uistep);

    /** Window
     *
     * Set a window in the logical space.
     */
    void window(int _v0, int _u0, int _vF, int _uF);

    /** Empty image
     *
     * A 0x0 image with no information about oversampling is produced. Remember
     * to re-initialise the oversampled image before starting to use it again.
     *
     * @code
     * IO.clear();
     * @endcode
     */
    void clear();

    /** Returns the exact position of a non-integer location
     *
     * The existing indexes (iu,iv) within the overimage are returned according
     * to a certain (u,v) and the Oversampled Image definition. Usually this
     * function is not needed in normal oversampled images operation as you can
     * access to pixels with fractional indexes. In our example of footprints,
     * this funciton would make the conversion between the oversampled position
     * (-1,1) to the image index (-51,51) (what is returned). Don't be mislead
     * by the physical position of the pixel (-51,51) which is (76,178).
     *
     * @code
     * IO.over2img(y, x, iy, ix);
     * @endcode
     */
    void over2img(double v, double u, int& iv, int& iu) const
    {
        if (v < overvmin || v > overvmax)
            REPORT_ERROR(1505, "ImgeOver::over2img: v out of range");

        if (u < overumin || u > overumax)
            REPORT_ERROR(1505, "ImgeOver::over2img: u out of range");

        iu = (int) ROUND((((u) - overumin) * uistep));
        iv = (int) ROUND((((v) - overvmin) * vistep));
    }

    /** Speed up pixel index macro
     *
     * This macro is the same as the function over2img but faster due to no
     * function call is performed.
     *
     * @code
     * OVER2IMG(IO, y, x, iy, ix);
     * @endcode
     */
#define OVER2IMG(IO, v, u, iv, iu) \
    iu = (int) ROUND((((u)-(IO).overumin) * (IO).uistep)); \
    iv = (int) ROUND((((v)-(IO).overvmin) * (IO).vistep));

    /** Returns the logical position of a "physical" location
     *
     * This function makes the conversion from image logical location to
     * oversampled logical position. In our example of footprints the conversion
     * between the oversampled position (-51,51) to the image index (-1,1) (what
     * is returned). Notice that the "physical" position of the pixel (-51,51)
     * is not the true physical one, that would be (76,178).
     *
     * @code
     * IO.img2over(iv, iu, v, u);
     * @endcode
     */
    void img2over(int iv, int iu, double &v, double &u) const
    {
        if (iu < 0 || iu > XSIZE(img))
            REPORT_ERROR(1505, "ImageOver::img2over: iu out of range");

        if (iv < 0 || iv > YSIZE(img))
            REPORT_ERROR(1505, "ImageOver::img2over: iv out of range");

        u = (double) overumin + iu / (double) uistep;
        v = (double) overvmin + iv / (double) vistep;
    }

    /** Speed up pixel index macro
     *
     * This macro is the same as the function img2over but faster due to no
     * function call is performed.
     *
     * @code
     * IMG2OVER(IO, iy, ix, y, x);
     * @endcode
     */
#define IMG2OVER(IO, iv, iu, v, u) \
    u = (double) (IO).overumin + (iu) / (double) ((IO).uistep); \
    v = (double) (IO).overvmin + (iv) / (double) ((IO).vistep);

    /** Constant pixel access with fractional indexes
     *
     * This function returns you a constant pixel access referred with a
     * fractional index. If you want to access to the pixels in the classic
     * Image fashion then you should use the macro IMGPIXEL. An exception is
     * thrown if you try to access an image outside the image.
     *
     * @code
     * cout << IO(1.34,-0.56) << endl;
     * @endcode
     */
    double operator()(double v, double u) const
    {
        if (v < overvmin || v > overvmax)
            REPORT_ERROR(1505, "ImgeOver::over2img: v out of range");

        if (u < overumin || u > overumax)
            REPORT_ERROR(1505, "ImgeOver::over2img: u out of range");

        int iu, iv;
        OVER2IMG(*this, v, u, iv, iu);

        return MAT_ELEM(img, iv, iu);
    }

    /** Pixel access with fractional indexes
     *
     * This function allows you to set pixels referred with a fractional index.
     * If you want to access to the pixels in the classic Image fashion then you
     * should use the macro IMGPIXEL.
     *
     * @code
     * IO(1.34, -0.56) = 1;
     * @endcode
     */
    double& operator()(double v, double u)
    {
        int iu, iv;
        OVER2IMG(*this, v, u, iv, iu);

        return MAT_ELEM(img, iv, iu);
    }

    /** Speed up pixel access macro
     *
     * This macro is a speeded up version of the pixel access routines for
     * oversampled images.
     *
     * @code
     * OVERPIXEL(IO, 1.34, -0.56) = 1;
     * @endcode
     */
#define OVERPIXEL(IO, y, x) IMGPIXEL((IO), \
                                     (int) ROUND(((u) * (IO).uistep)), \
                                     (int) ROUND(((v) * (IO).vistep)))

    // The following two functions have been redefined (they should be
    // inherited from Image) because the compiler doesn't admit this
    // inheritance and the redefinition of the fractional index pixel
    // access.

    /** Matrix access
     *
     * This function allows you to access the matrix2D over which the complex
     * data type is based.
     *
     * @code
     * cout << IO();
     * @endcode
     */
    matrix2D< double >& operator()()
    {
        return img;
    }

    const matrix2D< double >& operator()() const
    {
        return img;
    }

    /** Maximum value in not sampled units on U axis
     *
     * In our example: 2.
     *
     * @code
     * int Umax = IO.umax();
     * @endcode
     */
    int umax() const
    {
        return overumax;
    }

    /** Maximum value in not sampled units on V axis
     *
     * In our example: 2.
     *
     * @code
     * int Vmax = IO.vmax();
     * @endcode
     */
    int vmax() const
    {
        return overvmax;
    }

    /** Minimum value in not sampled units on U axis
     *
     * In our example: -2.
     *
     * @code
     * int Umin = IO.umin();
     * @endcode
     */
    int umin() const
    {
        return overumin;
    }

    /** Minimum value in not sampled units on V axis
     *
     * In our example: -2.
     *
     * @code
     * int Vmin = IO.vmin();
     * @endcode
     */
    int vmin() const
    {
        return overvmin;
    }

    /** Sampling rate in U axis
     *
     * In our example: 51.
     *
     * @code
     * int Usampling = IO.Ustep();
     * @endcode
     */
    int Ustep() const
    {
        return uistep;
    }

    /** Sampling rate in V axis
     *
     * In our example: 51.
     *
     * @code
     * int Vsampling = IO.Vstep();
     * @endcode
     */
    int Vstep() const
    {
        return vistep;
    }

    /** Generate the normal image by averaging
     *
     * This function passes from an oversampled image to a normal one, in our
     * example from the 250x250 image to a 5x5. A pointer to the normal image
     * must be given.
     *
     * @code
     * IO.down_sample(&I);
     * @endcode
     */
    void downsample(Image* I) const;

    /** Generate an oversampled image from a normal one
     *
     * This function passes from a normal image to an oversampled one, in our
     * example from 5x5 to 250x250. First you have to initialise the parameters
     * of the oversampled image and then the oversampling is done by bilinear
     * interpolation.
     *
     * @code
     * oversample(&I);
     * @endcode
     *
     * FIXME NOT IMPLEMENTED!!
     */
    void oversample(Image *I) const;
};

/* Extensions for the Imagic header and image files
 */
static const char* IMAGIC_HEADER_EXT = "hed";
static const char* IMAGIC_IMAGE_EXT = "img";

/* Length of the tag
 */
static const short unsigned IMAGIC_TAG_LEN = 7;

/* Separator character
 */
static const char IMAGIC_TAG_SEP = ':';

/* Constants for the offset into the Imagic header file
 */
static const unsigned IMAGIC_IFOL_OFFSET = 4, IMAGIC_IXLP_OFFSET = 48,
        IMAGIC_IYLP_OFFSET = 52, IMAGIC_TYPE_OFFSET = 56, IMAGIC_WORD_LEN = 4,
                             IMAGIC_RECORD_LEN = 1024;

/* Constants defining the Imagic header and some of its fields
 */
static const unsigned int IMAGIC_HEADER_BLOCK_SIZE = 256;
static const unsigned int IMAGIC_IDX_IMN = 0, IMAGIC_IDX_IFOL = 1,
        IMAGIC_IDX_NHFR = 3, IMAGIC_IDX_NDATE = 5, IMAGIC_IDX_NMONTH = 4,
                                                // These seem reversed in the spec
                                                IMAGIC_IDX_NYEAR = 6, IMAGIC_IDX_NHOUR = 7, IMAGIC_IDX_NMINUT = 8,
                                                                   IMAGIC_IDX_NSEC = 9, IMAGIC_IDX_NPIX2 = 10, IMAGIC_IDX_NPIXEL = 11,
                                                                                                           IMAGIC_IDX_IXLP1 = 12, IMAGIC_IDX_IYLP1 = 13, IMAGIC_IDX_TYPE = 14,
                                                                                                                              IMAGIC_IDX_NAME = 29, IMAGIC_IDX_NAMELEN = 80, IMAGIC_IDX_ARCHTYPE = 68;

/** Types of Imagic files
 *
 * Valid types are IMAGIC_REAL, IMAGIC_INTG, IMAGIC_PACK, IMAGIC_COMP.
 */
enum ImageImagicType
{
    IMAGIC_REAL, IMAGIC_INTG, IMAGIC_PACK, IMAGIC_COMP
};

/* structure that holds Imagic image information */
struct ImageImagicInfo
{
    unsigned int num_img;
    unsigned int xsize, ysize;
    vector< ImageImagicType > img_types;
};

/** Imagic Image class
 */
template<typename T>
class ImageImagicT: public ImageT< T >
{
public:
    /** Empty constructor
     */
    ImageImagicT() : ImageT< T >()
    {
        name_parsed = false;
    }

    /** Constructor with size
     */
    ImageImagicT(int Ydim, int Xdim) :
            ImageT< T >(Ydim, Xdim)
    {
        name_parsed = false;
    }

    /** Constructor with image name
     */
    ImageImagicT(FileName _name) : ImageT< T >(_name)
    {
        name_parsed = false;
    }

    /** Copy constructor
     */
    ImageImagicT(const ImageImagicT& I) : ImageT< T >(I)
    {
        if ((name_parsed = I.name_parsed))
        {
            hedfname = I.hedfname;
            imgfname = I.imgfname;
            imgnum = I.imgnum;
        }
    }

    /** Rename
     */
    virtual void rename(FileName newName)
    {
        if (newName != ImageT< T >::fn_img)
        {
            ImageT< T >::rename(newName);
            name_parsed = false;
        }
    }

    /** Clear
     */
    virtual void clear()
    {
        ImageT< T >::clear();
        name_parsed = false;
    }

    /** Assignment operator
     */
    ImageImagicT< T >& operator=(const ImageImagicT< T >& I)
    {
        if (this != &I)
        {
            this->ImageT< T >::operator=(I);
            if ((name_parsed = I.name_parsed))
            {
                hedfname = I.hedfname;
                imgfname = I.imgfname;
                imgnum = I.imgnum;
            }
        }
        return (*this);
    }

    /** Assignment operator
     */
    ImageImagicT< T >& operator= (const ImageT< T >& I)
    {
        if (this != &I)
        {
            this->ImageT< T >::operator=(I);
            name_parsed = false;
        }
        return (*this);
    }

    /** Assignment operator
     */
    template<typename T1>
    ImageImagicT< T >& operator=(const matrix2D< T1 >& m)
    {
        if (&ImageT< T >::img != (matrix2D< T >*)& m)
        {
            this->ImageT< T >::operator=(m);
            name_parsed = false;
        }
        return (*this);
    }

    /** Show
     */
    friend ostream& operator<<(ostream& out, const ImageImagicT< T >& I)
    {
        out << (ImageT< T >&) I;

        /* COSS: These functions are not defined
           out << "IMAGIC header fname: " << get_hedfname() << endl;
           out << "IMAGIC image fname:  " << get_imgfname() << endl;
           out << "image number:        " << get_imgnum() << endl;
        */

        return (out);
    }

    /** Read file
     */
    bool read(const FileName& name)
    {
        rename(name);
        ImageImagicInfo img_info = ImagicGetImgInfo(getHedFname());

        FileName img_fname = getImgFname();
        if (img_fname == "")
            REPORT_ERROR(1501, "ImageImagic::read: File " + name +
                         " doesn't seem fit Imagic format");

        FILE* img_fh;
        if ((img_fh = fopen(img_fname.c_str(), "rb")) == NULL)
            REPORT_ERROR(1501, "ImageImagic::read: IMAGIC file " + img_fname +
                         " not found");

        const int imgnum = getImgNum();
        const size_t img_offset = img_info.xsize * img_info.ysize;

        // Read the image data
        ImageT<T>::img.resize(img_info.ysize, img_info.xsize);
        const bool reversed = false;
        switch (img_info.img_types[imgnum])
        {
        case IMAGIC_REAL:
        {
            float data;
            const unsigned size = 4;
            fseek(img_fh, imgnum * size * img_offset, SEEK_SET);
            FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(ImageT< T >::img)
            {
                FREAD(&data, size, 1, img_fh, reversed);
                MULTIDIM_ELEM(ImageT< T >::img, i) = data;
            }
            break;
        }
        case IMAGIC_INTG:
        {
            short int data;
            const unsigned size = 2;
            fseek(img_fh, imgnum * size * img_offset, SEEK_SET);
            FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(ImageT< T >::img)
            {
                FREAD(&data, size, 1, img_fh, reversed);
                MULTIDIM_ELEM(ImageT< T >::img, i) = data;
            }
            break;
        }
        case IMAGIC_PACK:
        {
            unsigned char data;
            const unsigned size = 1;
            fseek(img_fh, imgnum * size * img_offset, SEEK_SET);
            FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(ImageT< T >::img)
            {
                FREAD(&data, size, 1, img_fh, reversed);
                MULTIDIM_ELEM(ImageT< T >::img, i) = data;
            }
            break;
        }
        default:
            REPORT_ERROR(1501,
                         "ImageImagicType not supported for this imgtype!");
            break;
        }
        fclose(img_fh);
        return true;
    }

    /** Write file
     *
     * FIXME Not implemented
     */
    void write(const FileName& name = "", bool reversed = false,
               Image_Type image_type = IBYTE)
    {
        REPORT_ERROR(1503, "ImageImagic::write: can't directly save");
    }

    /** Low level write
     *
     * FIXME Not implemented
     */
    void write(FILE*& fh, bool reversed, Image_Type image_type)
    {
        REPORT_ERROR(1503, "ImageImagic::write: can't directly save");
    }

    /** Get Header name
     */
    const FileName& getHedFname()
    {
        parseFname();
        return (hedfname);
    }

    /** Get Image name
     */
    const FileName& getImgFname()
    {
        parseFname();
        return (imgfname);
    }

    /** Get image number
     */
    int getImgNum()
    {
        parseFname();
        return (imgnum);
    }

    virtual void parseFname()
    {
        // Not meant to be called directly, but needs to
        // be kept public
        if (!name_parsed)
        {
            hedfname = "";
            imgfname = "";
            imgnum = -1;

            // Look for special IMAGIC format: 'imagic:<hedfile>:<imgnum>'
            if (ImageT< T >::fn_img.find(IMAGIC_TAG) != string::npos)
            {
                const string::size_type imgnumpos =
                    ImageT< T >::fn_img.rfind(IMAGIC_TAG_SEP);

                if (imgnumpos > IMAGIC_TAG_LEN)
                {
                    hedfname = ImageT< T >::fn_img.substr(IMAGIC_TAG_LEN,
                                                          imgnumpos - IMAGIC_TAG_LEN);
                    imgfname = hedfname.substitute_extension(IMAGIC_HEADER_EXT,
                               IMAGIC_IMAGE_EXT);
                    imgnum = atoi(
                                 (ImageT< T >::fn_img.substr(imgnumpos + 1)).c_str());
                }
            }
            name_parsed = true;
        }
    }

protected:
    bool name_parsed;
    FileName hedfname, imgfname;
    int imgnum;
};

/** Alias for ImageImagic type
 */
typedef ImageImagicT< double > ImageImagic;


/** Creates a string of the format 'imagic:hedfname:num' suitable for use as an
 * image name
 */
inline string ImagicMakeName(const char* hed_fname, unsigned int imgnum)
{
#if GCC_VERSION < 30300
    char aux[15];
    ostrstream ss(aux, sizeof(aux));
#else

    ostringstream ss;
#endif

    ss << IMAGIC_TAG << hed_fname << IMAGIC_TAG_SEP << imgnum;
    return (ss.str());
}

/** Looks at an IMAGIC header file for information about the images it contains
 */
const ImageImagicInfo ImagicGetImgInfo(const FileName& hed_fname);

/** Creates a new Imagic header/image file pair, filling it with the data
 * pointed to by the vector parameter
 */
template<typename T>
bool ImagicWriteImagicFile(const FileName& hed_fname,
                           const vector< ImageT< T >* > & imgs,
                           ImageImagicType img_type = IMAGIC_REAL)
{
    const FileName img_fname = hed_fname.substitute_extension(IMAGIC_HEADER_EXT,
                               IMAGIC_IMAGE_EXT);
    FILE* imagic_hed = fopen(hed_fname.c_str(), "wb");
    FILE* imagic_img = fopen(img_fname.c_str(), "wb");

    if (imagic_hed && imagic_img)
    {
        // Write the header information
        int header_block[IMAGIC_HEADER_BLOCK_SIZE];

        for (unsigned int imgcount = 0; imgcount < imgs.size(); imgcount++)
        {
            const ImageT< T >* image = imgs[imgcount];
            if (!image)
            {
                if (imagic_hed)
                    fclose(imagic_hed);

                if (imagic_img)
                    fclose(imagic_img);

                return (false);
            }

            memset(header_block, 0, sizeof(header_block));
            header_block[IMAGIC_IDX_IMN] = imgcount + 1;
            header_block[IMAGIC_IDX_IFOL] = imgs.size() - (imgcount + 1);
            header_block[IMAGIC_IDX_NHFR] = 1;

            const time_t nowt = time(NULL);
            const struct tm* nowtm = localtime(&nowt);

            header_block[IMAGIC_IDX_NMONTH] = nowtm->tm_mon + 1;
            header_block[IMAGIC_IDX_NDATE] = nowtm->tm_mday;
            header_block[IMAGIC_IDX_NYEAR] = nowtm->tm_year + 1900;
            header_block[IMAGIC_IDX_NHOUR] = nowtm->tm_hour;
            header_block[IMAGIC_IDX_NMINUT] = nowtm->tm_min;
            header_block[IMAGIC_IDX_NSEC] = nowtm->tm_sec;
            header_block[IMAGIC_IDX_NPIX2] = XSIZE((*image)()) *
                                             YSIZE((*image)());
            header_block[IMAGIC_IDX_NPIXEL] = header_block[IMAGIC_IDX_NPIX2];
            header_block[IMAGIC_IDX_IXLP1] = XSIZE((*image)());
            header_block[IMAGIC_IDX_IYLP1] = YSIZE((*image)());

            string formatstr;
            switch (img_type)
            {
            case IMAGIC_REAL:
                formatstr = "REAL";
                break;

            case IMAGIC_INTG:
                formatstr = "INTG";
                break;

            case IMAGIC_PACK:
                formatstr = "PACK";
                break;

            case IMAGIC_COMP:
                formatstr = "COMP";
                break;
            }

            memcpy(&header_block[IMAGIC_IDX_TYPE], formatstr.c_str(), 4);
            strncpy((char*) &header_block[IMAGIC_IDX_NAME],
                    image->name().c_str(), IMAGIC_IDX_NAMELEN);

#if defined(_LINUX)

            static const unsigned int ARCH_VAL = 33686018;
#elif defined(_SUN)

            static const unsigned int ARCH_VAL = 67372036;
#else

            static const unsigned int ARCH_VAL = 33686018;
#endif

            // This next line will generate an error if not using linux or sun!
            header_block[IMAGIC_IDX_ARCHTYPE] = ARCH_VAL;

            fwrite(header_block, sizeof(header_block), 1, imagic_hed);

            // Write the image data to the .img file
            FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY((*image)())
            {
                switch (img_type)
                {
                case IMAGIC_REAL:
                {
                    const float p = (float) MULTIDIM_ELEM((*image)(), i);
                    FWRITE(&p, sizeof(p), 1, imagic_img, false);
                    break;
                }

                case IMAGIC_INTG:
                {
                    const unsigned short p =
                        (unsigned short) MULTIDIM_ELEM((*image)(), i);
                    FWRITE(&p, sizeof(p), 1, imagic_img, false);
                    break;
                }

                case IMAGIC_PACK:
                case IMAGIC_COMP:
                    // NT: TODO: implement these
                    fprintf(stderr, "Unsupported format for writeSelFile!\n");
                    fclose(imagic_hed);
                    fclose(imagic_img);
                    return (false);
                    break;
                }
            }
        }
    }

    if (imagic_hed)
        fclose(imagic_hed);

    if (imagic_img)
        fclose(imagic_img);

    return true;
}

// Specialized function to read images with complex numbers in them
template<>
bool ImageImagicT< complex< double > >::read(const FileName&);

#endif
