/***************************************************************************
 *
 * Authors:      Alberto Pascual Montano (pascual@cnb.uam.es)
 *               Carlos Oscar Sanchez Sorzano
 *               Roberto Marabini
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef VOLUME_H
#define VOLUME_H

#include <iostream>
#include <cstdio>
#include <typeinfo>

#include "funcs.h"
#include "matrix3d.h"
#include "header.h"

/// @defgroup Volumes Volumes

typedef enum
{
    VBYTE = 1,
    VFLOAT = 2,
    VINT = 3,
    VUCHAR = 4,
    V16 = 5
} Volume_Type;

/** Basic volume class.
 * @ingroup Volumes
 *
 * The volume class is a general class which only contains the volume itself
 * and a filename for it. It has got a float matrix3D as member, and basically
 * all operations between volumes are based on that class.
 *
 * This class is the usual one when you want to operate volumes in memory.
 * Volumes belonging to this class cannot be saved, instead you could use the
 * class volumeXmipp which inherits from this class, so you have the whole
 * funcitonality of the class volume plus the possibility of saving at the end.
 *
 * See "Logical and physical access" for a detailed information about the voxel
 * accessing, and the conventions used in the volume definition.
 */
template<typename T>
class VolumeT
{
public:
    FileName fn_img; // name of the image
    matrix3D< T > img; // 3D matrix with the image

    /// @defgroup VolumesConstructors Constructors
    /// @ingroup Volumes

    /** Empty constructor.
     * @ingroup VolumesConstructors
     *
     * An empty Volume with size 0x0x0 is created.
     *
     * @code
     * VolumeT< double > vol;
     * @endcode
     */
    VolumeT()
    {
        fn_img = "";
    }

    /** Constructor with size.
     * @ingroup VolumesConstructors
     *
     * A blank Volume (0.0 filled) is created with the given size. Pay
     * attention to the dimension order: Z, Y and X.
     *
     * @code
     * VolumeT< double > vol(64, 64, 64);
     * @endcode
     */
    VolumeT(int Zdim, int Ydim, int Xdim)
    {
        img.resize(Zdim, Ydim, Xdim);
        fn_img = "";
    }

    /** Constructor using filename.
     * @ingroup VolumesConstructors
     *
     * The name for the volume is assigned but no load is performed. The volume
     * content is empty, ie, with size 0x0x0.
     *
     * @code
     * VolumeT< double > vol("art0001.vol");
     * @endcode
     */
    VolumeT(FileName name)
    {
        fn_img = name;
    }

    /** Copy constructor.
     * @ingroup VolumesConstructors
     *
     * @code
     * VolumeT< double > vol2(vol1);
     * @endcode
     */
    template<typename T2>
    VolumeT(const VolumeT< T2 >& v)
    {
        *this = v;
    }

    /// @defgroup VolumesOperations Some operations.
    /// @ingroup Volumes

    /** Assignment.
     * @ingroup VolumesOperations
     *
     * @code
     * vol2 = vol1;
     * @endcode
     */
    template<typename T2>
    VolumeT& operator=(const VolumeT< T2 >& v)
    {
        if (this != static_cast< VolumeT*>(&v))
        {
            type_cast(v.img, img);
            fn_img = v.fn_img;
        }

        return *this;
    }

    /** Another function for assigment.
     * @ingroup VolumesOperations
     */
    template<typename T2>
    void assign(const VolumeT< T2 >& v)
    {
        *this = v;
    }

    /** Assignment from matrix3D.
     * @ingroup VolumesOperations
     */
    template<typename T2>
    VolumeT& operator=(const matrix3D< T2 >& m)
    {
        // FIXME dangerous cast
        if (&img != (matrix3D< T >*)(&m));
        {
            fn_img = "";
            type_cast(m, img);
        }

        return *this;
    }

    /** Another function for assigment from matrix3D.
     * @ingroup VolumesOperations
     */
    template<typename T2>
    void assign(const matrix3D< T2 >& m)
    {
        *this = m;
    }

    /** Rename Volume.
     * @ingroup VolumesOperations
     *
     * Give a new name to the Volume.
     *
     * @code
     * vol.rename("new_name");
     * @endcode
     */
    virtual void rename(FileName newName)
    {
        fn_img = newName;
    }

    /** Empty Volume.
     * @ingroup VolumesOperations
     *
     * This function clears the Volume to a 0x0x0 Volume without name.
     *
     * @code
     * vol.clear();
     * @endcode
     */
    virtual void clear()
    {
        fn_img = "";
        img.clear();
    }

    /** Sets the origin of the Volume at its center.
     * @ingroup VolumesOperations
     *
     * The exact formula for the center of a Volume is defined in the
     * Conventions Section, this function modify the indexes of the Volume
     * in such a way that now the logical origin is at its center and the
     * indexes may take negative values. Use startingX, finishingX, ... to
     * setup loops.
     *
     * @code
     * vol.moveOriginTo_center();
     * @endcode
     */
    void moveOriginTo_center()
    {
        img.startingZ() = FIRST_XMIPP_INDEX(img.zdim);
        img.startingY() = FIRST_XMIPP_INDEX(img.ydim);
        img.startingX() = FIRST_XMIPP_INDEX(img.xdim);
    }

    /** Fill with 0 and move origin to center.
     * @ingroup VolumesOperations
     *
     * This function resizes the Volume to the given size, fill it with
     * 0.0 and then moves the Volume logical origin to its center. See
     * previous function.
     *
     * @code
     * vol.adapt_to_size(64, 64)
     * @endcode
     */
    void adapt_to_size(int Zdim, int Ydim, int Xdim)
    {
        img.initZeros(Zdim, Ydim, Xdim);
        moveOriginTo_center();
    }

    /// @defgroup VolumesAccess Image access.
    /// @ingroup Volumes

    /** 3D Matrix access.
     * @ingroup VolumesAccess
     *
     * This operator can be used to access the 3D matrix, and the matrix
     * operations defined in matrix3D. In this way we could resize a Volume
     * just by resizing its associated 3D matrix or we could add two Volumes
     * by adding their 3D matrices.
     *
     * @code
     * vol().resize(128, 128, 128);
     *
     * vol2() = vol1() + vol2();
     * @endcode
     */
    matrix3D< T >&  operator()()
    {
        return img;
    }

    const matrix3D< T >&  operator()() const
    {
        return img;
    }

    /** Get matrix3D.
     * @ingroup VolumesAccess
     */
    void get_matrix3D(matrix3D< T >& m)
    {
        m = img;
    }

    /** Set matrix3D.
     * @ingroup VolumesAccess
     */
    void set_matrix3D(const matrix3D< T >& m)
    {
        img = m;
    }

    /** Voxel access.
     * @ingroup VolumesAccess
     *
     * This operator is used to access a voxel within the Volume. This is a
     * logical access, so you could access to negative positions if the
     * Volume has been defined so (see the general explanation for the class).
     *
     * @code
     * std::cout << "Grey level of voxel (2,-3,-3) of the Volume = "
     *     << vol(2, -3, -3) << std::endl;
     *
     * vol(2, -3, -3) = vol(2, -3, -2);
     * @endcode
     */
    T& operator()(int z, int y, int x) const
    {
        return img(z, y, x);
    }

    /** Get voxel.
     * @ingroup VolumesAccess
     */
    T get_voxel(int z, int y, int x) const
    {
        return img(z, y, x);
    }

    /** Set voxel.
     * @ingroup VolumesAccess
     */
    void set_voxel(int z, int y, int x, T val)
    {
        img(z, y, x) = val;
    }

    /** Name access.
     * @ingroup VolumesAccess
     *
     * This function is used to know the name of the Volume. It cannot be used
     * to assign a new one.
     *
     * @code
     * std::cout << "Volume name: " << vol.name() << std::endl;
     * @endcode
     */
    FileName name() const
    {
        return static_cast< FileName >(fn_img);
    }

    /** Cout << Volume.
     * @ingroup VolumesAccess
     *
     * Shows name and size.
     */
    friend std::ostream& operator<<(std::ostream& out, const VolumeT< T >& V)
    {
        out << "Volume Name   : " << V.fn_img << std::endl
        << "dimensions   : " << V.img.SliNo() << " x "
        << V.img.rowNumber() << " x " << V.img.colNumber()
        << "  (slices x rows x columns)" << std::endl;

        return out;
    }

    /** @defgroup VolumesIO I/O functions.
     * @ingroup Volumes
     *
     * All these functions work with the image written in raw floats.
     */

    /** Read Volume from disk, given the Volume's dimensions.
     * @ingroup VolumesIO
     *
     * If the Volume doesn't exist at the given path then an exception is
     * thrown.
     *
     * The reversed flag indicates if the elements in element_size must be read
     * in a reversed way.
     *
     * Elements are supposed to be in the following order (y,x) =
     * (0,0)(0,1)(0,2), ..., (0,Xdim-1), (1,0), (1,1), ...
     *
     * The element size can be adjusted so that raw images of bytes (VBYTE),
     * unsigned ints of 16 bits (V16) and floats (VFLOAT) can be read.
     *
     * The headersize parameter can be used to read raw volumes preceded by a
     * header. This is the normal image format of many image processing
     * packages.
     *
     * @code
     * vol.read(65, 65, 65, "art0001.raw");
     * @endcode
     */
    void read(FileName name,
              int Zdim,
              int Ydim,
              int Xdim,
              bool reversed = false,
              Volume_Type volume_type = VBYTE,
              int header_size = 0)
    {
        FILE* fh;
        clear();
        fn_img = name;

        if ((fh = fopen(fn_img.c_str(), "rb")) == NULL)
            REPORT_ERROR(1501, "Volume::read: File " + fn_img + " not found");

        fseek(fh, header_size, SEEK_SET);
        read(fh, Zdim, Ydim, Xdim, reversed, volume_type);

        fclose(fh);
    }

    /** Read image from disk using a file pointer.
     * @ingroup VolumesIO
     *
     * This is the core routine of the previous one.
     */
    void read(FILE* fh,
              int Zdim,
              int Ydim,
              int Xdim,
              bool reversed,
              Volume_Type volume_type)
    {
        img.resize(Zdim, Ydim, Xdim);

        for (int i = 0; i < img.size; i++)
            switch (volume_type)
            {
            case VBYTE:
                unsigned char u;
                FREAD(&u, sizeof(unsigned char), 1, fh, reversed);
                MULTIDIM_ELEM(img, i) = static_cast< T >(u);
                break;
            case VINT:
                int ii;
                FREAD(&ii, sizeof(int), 1, fh, reversed);
                MULTIDIM_ELEM(img, i) = static_cast< T >(ii);
                break;
                // Integers and floats need to be reversed in identical way
            case V16:
                unsigned short us;
                FREAD(&us, sizeof(unsigned short), 1, fh, reversed);
                MULTIDIM_ELEM(img, i) = static_cast< T >(us);
                break;
            case VFLOAT:
                float f;
                FREAD(&f, sizeof(float), 1, fh, reversed);
                MULTIDIM_ELEM(img, i) = static_cast< T >(f);
                break;
            }
    }

    /** Write Volume to disk.
     * @ingroup VolumesIO
     *
     * If there is any problem in the writing, an exception is thrown. You can
     * give a name to the written Volume different from the one used when it
     * was read. From this point the filename of the Volume has changed. This
     * is somehow like the "Save as ..." and "Save".
     *
     * @code
     * vol.write(); //Save
     * vol.write("art0002.raw"); // Save as
     * @endcode
     */
    void write(FileName name = "",
               bool reversed = false,
               Volume_Type volume_type = VBYTE)
    {
        FILE* fp;
        if (name != "")
            rename(name);

        if ((fp = fopen(fn_img.c_str(), "wb")) == NULL)
            REPORT_ERROR(1503, "Volume::write: File " + fn_img +
                         " cannot be saved");

        write(fp, reversed, volume_type);
        fclose(fp);
    }

    /** Write image to disk using a file pointer.
     * @ingroup VolumesIO
     *
     * This is the core routine of the previous one.
     */
    void write(FILE* fh, bool reversed, Volume_Type volume_type)
    {
        if (XSIZE(img) == 0 || YSIZE(img) == 0 || ZSIZE(img) == 0)
            return;

        double a, b;
        if (volume_type != VFLOAT)
        {
            double min_val, max_val;

            (*this)().computeDoubleMinMax(min_val, max_val);

            if (volume_type == VBYTE)
                a = 255;
            else
                a = 65535;

            a /= (max_val - min_val);
            b = min_val;
        }

        for (int i = 0; i < img.size; i++)
            switch (volume_type)
            {
            case VBYTE:
                unsigned char u;
                u = static_cast< unsigned char >(ROUND(a *
                                                       (MULTIDIM_ELEM(img, i) - b)));
                FWRITE(&u, sizeof(unsigned char), 1, fh, reversed);
                break;
            case V16:
                unsigned short us;
                us = static_cast< unsigned short >(ROUND(a *
                                                   (MULTIDIM_ELEM(img, i) - b)));
                FWRITE(&us, sizeof(unsigned short), 1, fh, reversed);
                break;
            case VFLOAT:
                float f;
                f = static_cast< float >(MULTIDIM_ELEM(img, i));
                FWRITE(&f, sizeof(float), 1, fh, reversed);
                break;
            case VINT:
                int ii;
                ii = static_cast< int >(MULTIDIM_ELEM(img, i));
                FWRITE(&ii, sizeof(int), 1, fh, reversed);
                break;
            }
    }
};


// Specialization for complex numbers
template<>
inline void VolumeT< complex< double > >::read(FILE* fh,
        int Zdim,
        int Ydim,
        int Xdim,
        bool reversed,
        Volume_Type volume_type)
{
    img.resize(Zdim, Ydim, Xdim);

    for (int i = 0; i < img.size; i++)
    {
        float a, b;

        // read real part
        FREAD(&a, sizeof(float), 1, fh, reversed);

        // read imaginary part
        FREAD(&b, sizeof(float), 1, fh, reversed);

        // Assign the number
        std::complex< double > c(a, b);

        MULTIDIM_ELEM(img, i) = c;
    }
}

template<>
inline void VolumeT< complex< double > >::write(FILE *fh,
        bool reversed,
        Volume_Type volume_type)
{
    for (int i = 0; i < img.size; i++)
    {
        float a, b;
        a = static_cast< float >((MULTIDIM_ELEM(img, i)).real());
        b = static_cast< float >((MULTIDIM_ELEM(img, i)).imag());
        FWRITE(&a, sizeof(float), 1, fh, reversed);
        FWRITE(&b, sizeof(float), 1, fh, reversed);
    }
}

/** Write a volume in CCP4 format.
 * @ingroup VolumesIO
 *
 * See CCP4 format at http://www.ccp4.ac.uk/dist/html/maplib.html. Ts is the
 * sampling rate in angstroms per pixel
 */
void write_as_CCP4(VolumeT< double>* V, const FileName& fn, double Ts = 1);

/** @defgroup VolumesSpeedUp Speed up macros
 * @ingroup Volumes
 */

/** 3D Matrix access.
 * @ingroup VolumesSpeedUp
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
#define VOLMATRIX(V) ((V).img)

/** Array access.
 * @ingroup VolumesSpeedUp
 *
 * This macro allows you to access to the tridimensional array behind the image
 * (float***).
 */
#define VOLARRAY(V) (VOL_ARRAY(V).img)

/** Voxel access.
 * @ingroup VolumesSpeedUp
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
#define VOLVOXEL(V, k, i, j) VOL_ELEM(((V).img), (k), (i), (j))

/** Physical voxel access.
 * @ingroup VolumesSpeedUp
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
#define DIRECT_VOLVOXEL(V, k, i, j) DIRECT_VOL_ELEM(((V).img), (k), (i), (j))

/** Xmipp 3D Volumes.
 *
 * The Xmipp volume is a normal volume (inherited from volume class) plus a
 * Spider header. This is the appropiate class to work with volumes in memory
 * which we want to save later. The data in the header is not directly
 * accesible from the programmer and must be set through object functions,
 * in this way the coherence in the header is assured. See File Formats for
 * more information about the Spider format.
 *
 * In principle, the volume starts at voxel (0,0) but this can be modified
 * for any other logical access. See class Volume for more information.
 *
 * The Euler angles are useless in a Xmipp volume, and although the Spider
 * header has got space for them they are not used and cannot be accessed.
 *
 */
template<typename T>
class VolumeXmippT: public VolumeT< T >
{
protected:
    headerXmipp header;

public:

    /** Empty constructor.
     * @ingroup VolumesConstructors
     *
     * Creates an empty (0x0x0) image with no information in the header.
     *
     * @code
     * VolumeXmipp vol;
     * @endcode
     */
    VolumeXmippT(): VolumeT< T >()
    {
        if (typeid(T) == typeid(double))
            header.headerType() = headerXmipp::VOL_XMIPP;
        else if (typeid(T) == typeid(int))
            header.headerType() = headerXmipp::VOL_INT;
        else if (typeid(T) == typeid(complex< double >))
            header.headerType() = headerXmipp::VOL_FOURIER;
        else
        {
            std::cout << "\nError: VolumeXmipp should be" <<
            " complex< double >, double or integer\n";
            exit(0);
        }
    }

    /** Constructor with size.
     * @ingroup VolumesConstructors
     *
     * Creates a 0.0 filled volume of size Zdim x Ydim x Xdim.
     *
     * @code
     * VolumeXmipp< double > vol(64, 64, 64);
     * @endcode
     */
    VolumeXmippT(int Zdim, int Ydim, int Xdim) : VolumeT< T >(Zdim, Ydim, Xdim)
    {
        if (typeid(T) == typeid(double))
            header.headerType() = headerXmipp::VOL_XMIPP;
        else if (typeid(T) == typeid(int))
            header.headerType() = headerXmipp::VOL_INT;
        else if (typeid(T) == typeid(complex< double >))
            header.headerType() = headerXmipp::VOL_FOURIER;
        else
        {
            std::cout << "\nError: VolumeXmipp should be double or integer\n";
            exit(0);
        }

        header.set_dimension(Ydim, Xdim);
        header.Slices() = Zdim;
        header.set_header();
        header.set_time();
        header.set_date();
    }

    /** Constructor with filename, read from disk.
     * @ingroup VolumesConstructors
     *
     * The filename given must exist, then the file is loaded in the VolumeXmipp
     * class structure. You have loaded the volume at the declaration time.
     *
     * @code
     * VolumeXmipp< double > vol("art0001.vol");
     * @endcode
     */
    VolumeXmippT(FileName name) : VolumeT< T >(name)
    {
        if (typeid(T) == typeid(double))
            header.headerType() = headerXmipp::VOL_XMIPP;
        else if (typeid(T) == typeid(int))
            header.headerType() = headerXmipp::VOL_INT;
        else if (typeid(T) == typeid(complex< double >))
            header.headerType() = headerXmipp::VOL_FOURIER;
        else
        {
            std::cout << "\nError: VolumeXmipp should be double or integer\n";
            exit(0);
        }

        read(name);
    }

    /** Copy constructor.
     * @ingroup VolumesConstructors
     *
     * @code
     * VolumeXmipp< double > vol2(vol1);
     * @endcode
     */
    VolumeXmippT(const VolumeXmippT& I): VolumeT< T >(I)
    {
        header = I.header;
    }

    /** Empty image.
     * @ingroup VolumesConstructors
     *
     * All information is cleared.
     *
     * @code
     * vol.clear();
     * @endcode
     */
    void clear()
    {
        clear_header();
        VolumeT< T >::clear();
    }

    /** Show the header information of a Xmipp volume.
     * @ingroup VolumesOperations
     *
     * @code
     * std::cout << vol;
     * @endcode
     */
    friend std::ostream& operator<<(std::ostream& out,
                                    const VolumeXmippT< T >& v)
    {
        // FIXME Previous: old-style cast to reference
        out << static_cast< VolumeT< T > >(v) << v.header;
        return out;
    }

    /** Assignment from another Xmipp volume.
     * @ingroup VolumesOperations
     *
     * @code
     * VolumeXmipp vol1, vol2;
     * vol2 = vol1;
     * @endcode
     */
    VolumeXmippT& operator=(const VolumeXmippT< T >& op1)
    {
        if (&op1 != this)
        {
            this->VolumeT< T >::operator=(op1);
            header = op1.header;
        }

        return *this;
    }

    /** Assignment from a generic image.
     * @ingroup VolumesOperations
     *
     * @code
     * Volume< double > v;
     * VolumeXmipp< double > vol;
     * vol = v;
     * @endcode
     */
    VolumeXmippT& operator=(const VolumeT< T >& op1)
    {
        if (this != &op1)
        {
            this->VolumeT< T >::operator=(op1);
            clear_header();
            adjust_header();
        }

        return *this;
    }

    /** Assignment from a 3D matrix.
     * @ingroup VolumesOperations
     *
     * @code
     * matrix3D< float > m;
     * VolumeXmipp vol;
     * vol = m;
     * @endcode
     */
    template<typename T2>
    VolumeXmippT& operator=(const matrix3D< T2>& op1)
    {
        this->VolumeT< T >::operator=(op1);
        clear_header();
        adjust_header();

        return *this;
    }

    /** Assignment from any kind of volume.
     * @ingroup VolumesOperations
     *
     * @code
     * Volume< double > vol;
     * VolumeXmipp v;
     * v.assign_from(&vol);
     * @endcode
     */
    template<typename T2>
    void assign_from(VolumeT< T2 >* v)
    {
        *this = *v;
    }

    /** Read Xmipp volume from disk.
     * @ingroup VolumesIO
     *
     * If the volume doesn't exist at the given path then an exception is
     * thrown.
     *
     * The type check is a test performed on input image to check if it comes
     * from a big or little endian machine. It is done over the header field
     * fIform. Sometimes, this value is corrupted although the whole image is
     * still valid. You can skip this check and provide the reversed status via
     * force_reversed.
     *
     * @code
     * vol.read("art0001.vol");
     * @endcode
     */
    void read(const FileName& name,
              bool skip_type_check = false,
              bool force_reversed = false)
    {
        FILE* fp;

        rename(name);
        if ((fp = fopen(VolumeT< T >::fn_img.c_str(), "rb")) == NULL)
            REPORT_ERROR(1501,
                         static_cast< std::string >("VolumeXmipp::read: File " +
                                                    VolumeT< T >::fn_img + " not found"));

        // Read header
        if (!header.read(fp, skip_type_check, force_reversed))
            REPORT_ERROR(1502, "VolumeXmipp::read: File " +
                         VolumeT< T >::fn_img + " is not a valid Xmipp file");

        // Read whole image and close file
        VolumeT<T>::read(fp, header.iSlices(), header.iYdim(), header.iXdim(),
                         header.reversed(), VFLOAT);
        fclose(fp);

        header.set_header();  // Set header in a Xmipp consistent state
    }

    /** Write Xmipp volume to disk.
     * @ingroup VolumesIO
     *
     * If there is any problem in the writing, an exception is thrown. You can
     * give a name to the written volume different from the one used when it
     * was read. From this point the filename of the volume has changed. This
     * is somehow like the "Save as ..." and "Save".
     *
     * @code
     * vol.write(); // Save
     * vol.write("art0002.vol"); // Save as
     * @endcode
     *
     * If force_reversed is TRUE then image is saved in reversed mode, if not
     * it is saved in the same mode as it was loaded.
     */
    void write(const FileName& name = "", bool force_reversed = false)
    {
        FILE* fp;
        if (name != "")
            rename(name);

        if ((fp = fopen(VolumeT< T >::fn_img.c_str(), "wb")) == NULL)
            REPORT_ERROR(1503,
                         static_cast< std::string >("VolumeXmipp::write: File " +
                                                    VolumeT< T >::fn_img + " cannot be written"));

        adjust_header();
        header.write(fp, force_reversed);
        VolumeT< T >::write(fp, header.reversed(), VFLOAT);
        fclose(fp);
    }

    /// @defgroup VolumesHeader Header access
    /// @ingroup Volumes

    /** Adjust header.
     * @ingroup VolumesHeader
     *
     * Force header to have the dimensions of the image, time, date updated.
     */
    void adjust_header()
    {
        if (typeid(T) == typeid(double))
            header.headerType() = headerXmipp::VOL_XMIPP;
        else if (typeid(T) == typeid(int))
            header.headerType() = headerXmipp::VOL_INT;
        else if (typeid(T) == typeid(complex< double >))
            header.headerType() = headerXmipp::VOL_FOURIER;

        header.set_dimension(YSIZE(VolumeT< T >::img), XSIZE(VolumeT<T>::img));
        header.Slices() = ZSIZE(VolumeT< T >::img);
        header.set_time();
        header.set_date();
        header.set_title(VolumeT< T >::fn_img);
        header.set_header();
    }

    /** Resets header.
     * @ingroup VolumesHeader
     */
    void clear_header()
    {
        header.clear();
    }

    /** Change Filename.
     * @ingroup VolumesHeader
     *
     * @code
     * vol.rename("newName.spd");
     * @endcode
     */
    void rename(FileName newName)
    {
        VolumeT< T >::rename(newName);
        header.set_title(newName);
    }

    /** Reversed status.
     * @ingroup VolumesHeader
     *
     * This is used for the little/big endian process.
     */
    bool reversed() const
    {
        return header.reversed();
    }
};

// TODO Document
typedef VolumeT< double > Volume;
typedef VolumeXmippT< double > VolumeXmipp;
typedef VolumeT< complex< double > > FourierVolume;
typedef VolumeXmippT< complex< double > > FourierVolumeXmipp;

/** True if the given volume is an Xmipp volume.
 * @ingroup VolumesRelated
 */
int Is_VolumeXmipp(const FileName& fn,
                   bool skip_type_check = false,
                   bool force_reversed = false);

/** True if the given volume is a Fourier Xmipp volume.
 * @ingroup VolumesRelated
 */
int Is_FourierVolumeXmipp(const FileName& fn,
                          bool skip_type_check = false,
                          bool force_reversed = false);

/** Get size of a volume.
 * @ingroup VolumesRelated
 *
 * Returns -1 if the file is not an Xmipp volume.
 */
void GetXmippVolumeSize(const FileName& fn, int &Zdim, int &Ydim, int &Xdim);

#endif
