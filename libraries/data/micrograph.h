/***************************************************************************
 *
 * Authors: Carlos Oscar (coss@cnb.csic.es)
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

#ifndef _XMIPPMICROGRAPH_H
#define _XMIPPMICROGRAPH_H

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */
#include <vector>
#include "filters.h"
#include "xmipp_funcs.h"
#include "multidim_array.h"
#include "xmipp_image.h"
#include "xmipp_fftw.h"
#include "metadata.h"
#include "xmipp_error.h"

/* ************************************************************************* */
/* FORWARD DEFINITIONS                                                       */
/* ************************************************************************* */
// This forward definitions are needed for defining operators functions that
// use other clases type

/* ************************************************************************* */
/* MICROGRAPHY                                                               */
/* ************************************************************************* */
/// @defgroup Micrographs Micrographs
/// @ingroup DataLibrary
//@{
/** Particle coordinates.
    This structure stores the X,Y position of the particle. */
struct Particle_coords
{
    /// Label
    int label;
    /// X position
    int X;
    /// Y position
    int Y;
    /// Valid
    bool valid;
    /// Cost, scaled between 0 and 1
    double cost;
};

/** Micrography class.
    This class manages a large micrograph on disk. The image is not loaded
    into memory, that should avoid memory problems
*/
class Micrograph
{
public:
    /* This image will contain a single particle from the micrograph,
       this is done to avoid asking/freeing memory all time. */
    Image<double>            single_particle;
    std::vector<Particle_coords> coords;
    FileName                 fn_coords;
    FileName                 fn_micrograph;
    FileName                 fn_inf;
    int                      X_window_size;
    int                      Y_window_size;
    size_t                   Xdim,Ydim,Zdim,Ndim;
    int                      datatype;
    int                      swapbyte;
    int                      __offset;
    bool                     compute_transmitance;
    bool                     compute_inverse;
    int                      fh_micrograph;
    std::vector<std::string> labels;
    double                   stdevFilter;
public:
    Image<char>                * auxI;
    Image<unsigned char>       * IUChar;
    Image<short int>           * IShort;
    Image<unsigned short int>  * IUShort;
    Image<int>                 * IInt;
    Image<unsigned int>        * IUInt;
    Image<float>               * IFloat;
    /** Constructor */
    Micrograph();

    /** Destructor */
    ~Micrograph();

    /** Clear */
    void clear();

    /** Get micrograph datatype. */
    int getDatatype() const
    {
        return datatype;
    }
    /** Get micrograph datatype. */
    int getDatatypeDetph() const;

    /** Open micrograph.
        An exception is thrown if the file is not valid. */
    void open_micrograph(const FileName &fn_micrograph);

    /** Close micrograpgh.
        After working with the file, you must close it. */
    void close_micrograph();

    /** Get micrograph filename. */
    const FileName& micrograph_name()
    {
        return(fn_micrograph);
    }

    /** Set micrograph filename. */
    void set_micrograph_name(const FileName& fn)
    {
        fn_micrograph = fn;
    }
    void setStdevFilter(double d)
    {
        stdevFilter=d;
    }

    /** Save coordinates to disk. */
    void write_coordinates(int label, double minCost, const FileName &fn_coords = "");

    /** Read coordinates from disk.
        Coordinates are read into the selected family, the rest of
        families are untouched as well as the coordinates already belonging
        to this family */
    void read_coordinates(int label, const FileName &fn_coords);

    /** Transform all coordinates according to a 3x3 transformation
     * matrix */
    void transform_coordinates(const Matrix2D<double> &M);

    /** Multiply coordinates by a constant */
    void scale_coordinates(const double &c);

    /** Particle number.
        Number of particles in the coordinate list */
    int ParticleNo() const
    {
        return coords.size();
    }

    /** Particle.
        Return the list of particles. */
    std::vector<Particle_coords> & Particles()
    {
        return coords;
    }

    /** Another function for getting the particles.*/
    void get_Particles(std::vector<Particle_coords> & _coords)
    {
        _coords = coords;
    }

    /** Set selfWindow size.
        This selfWindow is set upon each coordinate and is used to cut all
        images. */
    void set_window_size(int _X_window_size, int _Y_window_size)
    {
        X_window_size = _X_window_size;
        Y_window_size = _Y_window_size;
    }

    /** Set Transmitance flag.
        When cutting images, 1/log10 is computed over the pixel values
        if this transmitance_flag=true. This function sets it
        Note: if pixel_value=0, no log is computed */

    void set_transmitance_flag(bool flag_value)
    {
        compute_transmitance =
            flag_value;
    }

    void setDataType(DataType _datatype)
    {
        datatype = _datatype;
    }
    /** Get Transmitance flag.
        When cutting images, 1/log10 is computed over the pixel values
        if transmitance_flag=true. This function reads it
        Note: if pixel_value=0, no log is computed */

    bool read_transmitance_flag(void)
    {
        return compute_transmitance;
    }

    /** Set Log flag.
        When cutting images, the contrast is inverted if inverse flag
        is true.  */

    void set_inverse_flag(bool flag_value)
    {
        compute_inverse =
            flag_value;
    }

    /** Get Log flag.
        When cutting images, the contrast is inverted
        if inverse flag=true.  This function reads it*/

    bool read_inverse_flag(void)
    {
        return compute_inverse;
    }

    /** Templated scissor function.
     *  This is the one actually doing the work
     */
    template <typename T>
    int templateScissor(const Image<T> &I,
                        const Particle_coords &P, MultidimArray<double> &result,
                        double Dmin, double Dmax, double scaleX, double scaleY,
                        bool only_check)
    {
        result.initZeros(Y_window_size, X_window_size);
        int _i0 = ROUND(scaleY * P.Y) + FIRST_XMIPP_INDEX(Y_window_size);
        int _iF = ROUND(scaleY * P.Y) + LAST_XMIPP_INDEX(Y_window_size);
        int _j0 = ROUND(scaleX * P.X) + FIRST_XMIPP_INDEX(X_window_size);
        int _jF = ROUND(scaleX * P.X) + LAST_XMIPP_INDEX(X_window_size);
        int retval = 1;
        double irange=1.0/(Dmax - Dmin);
		size_t i0 = (size_t)_i0;
		size_t iF = (size_t)_iF;
		size_t j0 = (size_t)_j0;
		size_t jF = (size_t)_jF;

        if (_i0 < 0 || iF >= Ydim || _j0 < 0 || jF >= Xdim)
        {
        	retval = 0;
        }
        else
            if (!only_check)
            {

                for (size_t i = i0; i <= iF; i++)
                {
                    int i_i0=i-i0;
                    for (size_t j = j0; j <= jF; j++)
                    {
                        int j_j0=j-j0;
                        double val=IMGPIXEL(I,i,j);
                        if (compute_transmitance)
                        {
                            double temp;
                            if (val < 1)
                                temp = val;
                            else
                                temp = log10(val);
                            if (compute_inverse)
                                A2D_ELEM(result,i_i0, j_j0) = (Dmax - temp) * irange;
                            else
                                A2D_ELEM(result,i_i0, j_j0) = (temp - Dmin) * irange;
                        }
                        else
                        {
                            if (compute_inverse)
                                A2D_ELEM(result, i_i0, j_j0) = (Dmax - val) * irange;
                            else
                                A2D_ELEM(result, i_i0, j_j0) = val;
                        }
                    }
                }
            }
        return retval;
    }

    /** Scissor.
        The single particle is selected by an index within the particle
        coordinate list. If the index is beyond the number of particles
        \ref ParticleNo , or the selfWindow size is not set (\ref set_window_size )
        an exception is thrown.

        Make sure that index n represents a valid particle before cutting it

        The scale affects the particle position, such that the position cut
        is pos*scale, but not the selfWindow size.

        If only check is true then the particle is not scissored, but
        the routine only checks if it can be done.

        Dmax and Dmin are used to invert the image and or compute the
        trnasmitance

        Returns 0 if an error ocurred and 1 if everything is all right*/
    int scissor(const Particle_coords &P, MultidimArray<double> &result,
                double Dmin, double Dmax,
                double scaleX = 1, double scaleY = 1, bool only_check = false);

    /** Access to array of 8 bits. */
    unsigned char * arrayUChar() const
    {
        return (*IUChar)().data;
    }

    /** Another function for access to array of 8 bits.*/
    //This is never used consider delete
    void get_arrayUChar(unsigned char * _m8)
    {
        _m8 = (*IUChar)().data;
    }

    /** Access to array of 16 bits. */
    short int * arrayShort() const
    {
        return (*IShort)().data;
    }

    /** Another function for access to array of 16 bits.*/
    void get_arrayShort(short int * _m16)
    {
        _m16 = (*IShort)().data;
    }

    /** Access to unsigned array of 16 bits. */
    unsigned short int * arrayUShort() const
    {
        return (*IUShort)().data;
    }

    /** Another function for access to unsigned array of 16 bits.*/
    void get_arrayUShort(unsigned short int * _um16)
    {
        _um16 = (*IUShort)().data;
    }
    /** Access to array of 32 bits int. */
    int * arrayInt() const
    {
        return (*IInt)().data;
    }

    /** Another function for access to array of 32 bits int.*/
    void get_arrayInt(int * _m32)
    {
        _m32 = (*IInt)().data;
    }

    /** Access to unsigned array of 32 bits unsig int. */
    unsigned int * arrayUInt() const
    {
        return (*IUInt)().data;
    }

    /** Another function for access to unsigned array of 32 bits unsigned int.*/
    void get_arrayUInt(unsigned int * _um32)
    {
        _um32 = (*IUInt)().data;
    }

    /** Access to array of 32 bits. */
    float * arrayFloat() const
    {
        return (*IFloat)().data;
    }

    /** Another function for access to array of 32 bits.*/
    void get_arrayfloat(float * _mf32)
    {
        _mf32 = (*IFloat)().data;
    }

    /** Pixel access for reading.
        These coordinates follow the physical Xmipp convention
        {../../../Extra_Docs/Conventions.html} for coordinates */
    double operator()(size_t y, size_t x) const
    {
        if (y < 0 || y >= Ydim || x < 0 || x >= Xdim)
            // COSS: REPORT_ERROR(1, "Micrograph::(): index out of range");
            return 0;
        if (datatype == DT_UChar)
        {
            return IMGPIXEL(*IUChar,y,x);
        }
        else if (datatype == DT_UShort)
        {
            return IMGPIXEL(*IUShort,y,x);
        }
        else if (datatype == DT_Short)
        {
            return IMGPIXEL(*IShort,y,x);
        }
        else if (datatype == DT_UInt)
        {
            return IMGPIXEL(*IUInt,y,x);
        }
        else if (datatype == DT_Int)
        {
            return IMGPIXEL(*IInt,y,x);
        }
        else if (datatype == DT_Float)
        {
            return IMGPIXEL(*IFloat,y,x);
        }
        else
            REPORT_ERROR(ERR_TYPE_INCORRECT, "Micrograph::(): unknown datatype");
    }

    /** Micrograph max min*/
    void computeDoubleMinMax(double &Dmin, double &Dmax) const
    {
        if (datatype == DT_UChar)
        {
            return (*IUChar)().computeDoubleMinMax(Dmin,Dmax);
        }
        else if (datatype == DT_UShort)
        {
            return (*IUShort)().computeDoubleMinMax(Dmin,Dmax);
        }
        else if (datatype == DT_Short)
        {
            return (*IShort)().computeDoubleMinMax(Dmin,Dmax);
        }
        else if (datatype == DT_UInt)
        {
            return (*IUInt)().computeDoubleMinMax(Dmin,Dmax);
        }
        else if (datatype == DT_Int)
        {
            return (*IInt)().computeDoubleMinMax(Dmin,Dmax);
        }
        else if (datatype == DT_Float)
        {
            return (*IFloat)().computeDoubleMinMax(Dmin,Dmax);
        }

        else
            REPORT_ERROR(ERR_TYPE_INCORRECT, "Micrograph::computeDoubleMinMax::(): unknown datatype");
    }

    /** Pixel access for writing. */
    //Dangerous function indeed
    //write in default endiam
    void set_val(int y, int x, double new_val)
    {
        if (datatype == DT_UChar)
        {
            IMGPIXEL(*IUChar,y,x) = (unsigned char) new_val;
        }
        else if (datatype == DT_UShort)
        {
            IMGPIXEL(*IUShort,y,x) = (unsigned short) new_val;
        }
        else if (datatype == DT_Short)
        {
            IMGPIXEL(*IShort,y,x) = (short) new_val;
        }
        else if (datatype == DT_UInt)
        {
            IMGPIXEL(*IUInt,y,x) = (unsigned int) new_val;
        }
        else if (datatype == DT_Int)
        {
            IMGPIXEL(*IInt,y,x) = (int) new_val;
        }
        else if (datatype == DT_Float)
        {
            IMGPIXEL(*IFloat,y,x) = (float) new_val;
        }

        else
            REPORT_ERROR(ERR_TYPE_INCORRECT, "Micrograph::set_val::(): unknown datatype");

    }

    /** Produce all single particle images.
        The file fn_micrograph+".sel" is also generated. The angle is the angle
        from the Y axis to the tilt axis, angles are positive clockwise.
        Images are rotated by -ang.
        If this angle is 0 no rotation is applied.*/
    void produce_all_images(int label, double minCost, const FileName &fn_root,
                            const FileName &fn_image = "", double ang = 0,
                            double gamma = 0., double psi = 0., bool rmStack=false);

    /** Search coordinate near a position.
        By default the precission is set to 3 pixels. The index of the coordinate
        within the list is returned. Returns -1 if none. */
    int search_coord_near(int x, int y, int prec = 3) const;

    /** Remove a coordinate from the coordinate list.
        An exception is thrown if the index is out of range within the
        coordinate list */
    void invalidate_coord(int n);

    /** Add coordinate.
        It returns the index of the particle added within the coordinate list. */
    int add_coord(int x, int y, int label, double cost);

    /** Move last coordinate to this position. */
    void move_last_coord_to(int x, int y);

    /** Access to coordinate structure.
        If the index is out of range then an exception is thrown. */
    Particle_coords & coord(int n)
    {
        if (n < 0 || n > ParticleNo())
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "Micrograph::coord(): index out of range");
        return coords[n];
    }

    /** Another function for accessing the coordinate structure.*/
    void get_coord(int n, Particle_coords &_coords)
    {
        _coords = coord(n);
    }

    /** Add label.
        The index assigned to the label is returned */
    int add_label(const std::string &label)
    {
        labels.push_back(label);
        return labels.size() - 1;
    }

    /** Number of labels. */
    int LabelNo()
    {
        return labels.size();
    }

    /** Get a label.
        An exception is thrown if the index is greater than the
        number of labels */
    std::string & get_label(int n)
    {
        if (n < 0 || n > LabelNo())
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "Micrograph::get_label(): index out of range");
        return labels[n];
    }

    /** Another function for get a lebel.*/
    void get_label(int n, std::string &_label)
    {
        _label = get_label(n);
    }

    /** Return micrograph size */
    void size(int &_Xdim, int &_Ydim) const
    {
        _Xdim = Xdim;
        _Ydim = Ydim;
    }

    /// Set micrograph size (when you do not read the file from disk)
    void resize(int Xdim, int Ydim, const FileName &filename="");

    /** Write micrograph.
        Set adjust to true if the values should be scaled within the minimum
        and maximum of the output range.
        Ex: m.write(fn8bits+"%uint8");
    */
    void write(const FileName &fileName,CastWriteMode castMode=CW_CAST);
};

/** Class for aligning two tilted micrographs */
class TiltPairAligner
{
public:
    /// Untilted coordinates
    std::vector<int> coordU;
    /// Tilted coordinates
    std::vector<int> coordT;
public:
    /// Empty constructor
    TiltPairAligner();

    /// Clear set of coordinates
    void clear();

    /// Add coordinates pair
    void addCoordinatePair(int _muX, int _muY, int _mtX, int _mtY);

    /// Adjust passing matrix
    void adjustPassingMatrix(int _muX, int _muY, int _mtX, int _mtY);

    /// Pass to tilted
    void passToTilted(int _muX, int _muY, int &_mtX, int &_mtY);

    /// Pass to untilted
    void passToUntilted(int _mtX, int _mtY, int &_muX, int &_muY);

    /// Compute gamma
    void computeGamma();

    /** Compute alphas.
     * Make sure of calling computeGamma before calling this function.
     */
    void computeAngles(double &ualpha, double &talpha, double &ogamma);
public:
    // For tilted-untilted correspondance
    Matrix2D<double>    Au;     // Untilted "positions"
    Matrix2D<double>    Bt;     // Tilted   "positions"
    Matrix2D<double>    Put;    // Passing matrix from untilted to tilted
    Matrix2D<double>    Ptu;    // Passing matrix from tilted to untilted
    int                 Nu;     // Number of points in matrices
    double              gamma;  // Tilting angle in radians
    double              alpha_u;// Angle of axis with X axis in radians
    double              alpha_t;// Anfle of axis with X axis in radians
    /// Auxiliary vector
    Matrix1D<double> m;
    Matrix2D<double> pair_E;
};

//@}

#endif
