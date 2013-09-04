/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef MASK_H
#define MASK_H

#include "multidim_array.h"
#include "multidim_array_generic.h"
#include "histogram.h"
#include "blobs.h"

#include "xmipp_program.h"
#include "metadata.h"
#include "args.h"
#include "wavelet.h"

#include "xmipp_program.h"
#include "metadata_extension.h"


void apply_geo_binary_2D_mask(MultidimArray< int > &mask,
                              const Matrix2D< double >& A);
void apply_geo_cont_2D_mask(MultidimArray< double >& mask,
                            const Matrix2D< double >& A);

/// @defgroup Masks Masks
/// @ingroup DataLibrary
//@{

#define INNER_MASK 1
#define OUTSIDE_MASK 2
#define NO_ACTIVATE 0
#define ACTIVATE 1

///@name Actual masks
//@{
/** Creates a RaisedCosine mask for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal RaisedCosines, and outside masks are
 * 1 - RaisedCosine. When entering, the mask is initialiazed to 0 and then the
 * mask is created.
 */
void RaisedCosineMask(MultidimArray< double >& mask,
                      double r1, double r2, int mode = INNER_MASK, double x0 = 0,
                      double y0 = 0, double z0 = 0);

/** Creates a RaisedCrown mask for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created within the
 * two the radii indicated with an extra region of (pix_width) pixels. The only
 * two valid modes are INNER_MASK (by default) or OUTSIDE_MASK. Inner mask are
 * normal RaisedCrowns, and outside masks are 1 - RaisedCrowns. When entering,
 * the mask is initialiazed to 0 and then the mask is created.
 */
void RaisedCrownMask(MultidimArray< double >& mask,
                     double r1, double r2, double pix_width,
                     int mode = INNER_MASK,
                     double x0 = 0, double y0 = 0, double z0 = 0);

/** Kaiser selfWindow
 *  The mask is resized.
 *  delta=ripple (in natural units) in the pass band.
 *  Deltaw=transition bandwidth (normalized to 1.0).
 */
void KaiserMask(MultidimArray<double> &mask, double delta = 0.01,
                double Deltaw = 1.0 / 12.0);

/** Creates a sinc mask for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0), by default (0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal sincs, and outside masks are 1 -
 * sinc. When entering the mask is initialiazed to 0 and then the mask is
 * created.
 *
 * Remind that sinc(w*n) is zero at n=1/w;
 */
void SincMask(MultidimArray< double >& mask,
              double omega, int mode = INNER_MASK, double x0 = 0, double y0 = 0, double z0 = 0);

/** Creates a radial-sinc-kaiser mask, the mask is resized.
 *  This function returns a sinc mask windowed by a Kaiser selfWindow.
 *  delta=ripple (in natural units) in the pass band.
 *  Deltaw=transition bandwidth (normalized to 1).
 *  omega=low pass frequency (normalized to 1).
 */
void SincKaiserMask(MultidimArray<double> &mask,
                    double omega, double delta = 0.01, double Deltaw = 1.0 / 12.0);

/** Blackman selfWindow
 *
 * It receives no parameter.
 */
void BlackmanMask(MultidimArray< double >& mask, int mode = INNER_MASK,
                  double x0 = 0, double y0 = 0, double z0 = 0);

/** Creates a sinc-blackman mask, the mask is resized
 *
 * This function returns a sinc mask windowed by a Blackman selfWindow. The selfWindow
 * is designed to cover a certain power of the sinc
 */
void SincBlackmanMask(MultidimArray< double >& mask,
                      double omega, double power_percentage,
                      int mode = INNER_MASK,
                      double x0 = 0, double y0 = 0, double z0 = 0);

/** Creates a circular mask for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. When entering the mask is initialiazed to 0 and then the mask
 * is created.
 */
void BinaryCircularMask(MultidimArray< int >& mask,
                        double radius, int mode = INNER_MASK,
                        double x0 = 0, double y0 = 0, double z0 = 0);


/** Creates a circular mask with blob-shaped edges for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. When entering the mask is initialiazed to 0 and then the mask
 * is created.
 */
void BlobCircularMask(MultidimArray<double> &mask,
                      double r1, blobtype blob, int mode,
                      double x0 = 0, double y0 = 0, double z0 = 0);


/** Creates a crown mask for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created with the two
 * radii indicated. The only two valid modes are INNER_MASK (by default, between
 * the two radii) or OUTSIDE_MASK (the negative of the crown). It is supposed
 * that R1 is smaller than R2.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinaryCrownMask(MultidimArray< int >& mask,
                     double R1, double R2, int mode = INNER_MASK,
                     double x0 = 0, double y0 = 0, double z0 = 0);

/** Creates a crown mask  with blob-shaped edges for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created with the two
 * radii indicated. The only two valid modes are INNER_MASK (by default, between
 * the two radii) or OUTSIDE_MASK (the negative of the crown). It is supposed
 * that R1 is smaller than R2.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BlobCrownMask(MultidimArray<double> &mask,
                   double r1, double r2, blobtype blob, int mode,
                   double x0 = 0, double y0 = 0, double z0 = 0);

/** Creates a gaussian mask for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal sincs, and outside masks are 1 -gaussian.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void GaussianMask(MultidimArray< double >& mask,
                  double sigma, int mode = INNER_MASK,
                  double x0 = 0, double y0 = 0, double z0 = 0);

/** Binary Circular 2D mask in wavelet space */
void BinaryDWTCircularMask2D(MultidimArray< int >& mask,
                             double radius, int smin, int smax,
                             const std::string& quadrant);

/** Creates a 2D separable-sinc-kaiser mask, the mask is resized.
 *  This function returns a sinc mask windowed by a Kaiser selfWindow.
 *  delta=ripple (in natural units) in the pass band.
 *  Deltaw=transition bandwidth (normalized to 1).
 *  omega=low pass frequency (normalized to 1).
 */
void SeparableSincKaiserMask2D(MultidimArray<double> &mask,
                               double omega, double delta = 0.01,
                               double Deltaw = 1.0 / 12.0);

/** Creates a 3x3 mask with value (1 by default) for those 4-neighbours of the
 * central point (0 otherwise).
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask2D_4neig(MultidimArray< int >& mask, int value = 1, int center = NO_ACTIVATE);

/** Creates a 3x3 mask with value1 for those 4-neighbors of the central point
 * and value2 for the 8 neighbours.
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask2D_8neig(MultidimArray< int >& mask, int value1 = 1, int value2 = 1,
                  int center = NO_ACTIVATE);

/** Creates a 3D DWT spherical for already sized masks
 *
 * The mask size must be a power of 2. The radius must be expressed in pixel
 * units corresponding to the size of the image. For instance, a 64x64x64 image
 * might have a radius of 32 pixels to concentrate on the central part only.
 *
 * If quadrant=xxx then 001,010,011,100,101,110 and 111 are generated together
 */
void BinaryDWTSphericalMask2D(MultidimArray< int >& mask,
                              double radius, int smin, int smax,
                              const std::string& quadrant);

/** Creates a 3D Cylinder mask for already sized masks
 *
 * The mask is supposed to be already resized and with its logical origin
 * defined. A cylinder placed logically at (z0,x0,y0), by default (0,0,0), is
 * created with the radius and height indicated. The only two valid modes ar
 * INNER_MASK (by default) or OUTSIDE_MASK. When entering the mask is
 * initialiazed to 0 and then the mask is created.
 */
void BinaryCylinderMask(MultidimArray< int >& mask,
                        double R, double H, int mode = INNER_MASK,
                        double x0 = 0, double y0 = 0, int z0 = 0);

/** Creates a 3D frame mask for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * square placed logically at (x0,y0,z0), by default (0,0,0), is created with
 * the two rectangular dimensions indicated. The only two valid modes are
 * INNER_MASK (by default, between the two radii) or OUTSIDE_MASK (the negative
 * of the crown). It is supposed that R1 is smaller than R2.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinaryFrameMask(MultidimArray< int >& mask,
                     int Xrect, int Yrect, int Zrect, int mode = INNER_MASK,
                     double x0 = 0, double y0 = 0, double z0 = 0);

/** Creates a 3D cone mask for already sized masks
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * cone with angle theta is placed parallel to the Z-axis and centered at
 * (x0,y0,z0). The only two valid modes are INNER_MASK (by default, inside the
 * crown is zero, outside is 1) or OUTSIDE_MASK (the negative of the crown).
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinaryConeMask(MultidimArray< int >& mask,
                    double theta, int mode = INNER_MASK, bool centerOrigin=false);

/** Creates a 3D missing wedge mask
 *
 * The mask is supposed to be resized and with its logical origin already set.
 * theta0 and thetaF are the tilting angles (around the Y-axis) wherin between
 * the data is supposed to be collected. In this region the mask will be one,
 * outside (ie in the missing wedge) it will be zero. The mask is centered at
 * (x0,y0,z0), and rotated with respect to euler angle matrix A.
 */
void BinaryWedgeMask(MultidimArray< int >& mask, double theta0, double thetaF,
                     const Matrix2D< double > &A, bool centerOrigin=false);


/** Creates a 3x3x3 mask with value (1 by default) for those 6-neighbors of the
 * central point (0 otherwise).
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask3D_6neig(MultidimArray< int >& mask, int value = 1, int center = NO_ACTIVATE);

/** Creates a 3x3x3 mask with value1 (1 by default) for those 6-neighbors and
 * value2 for the 18 neighbors of the central point (0 otherwise).
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask3D_18neig(MultidimArray< int >& mask, int value1 = 1, int value2 = 1,
                   int center = NO_ACTIVATE);

/** Creates a 3x3x3 mask with value1 (1 by default) for those 6-neighbors,
 * value2 for the 18 neighbors and value3 for the 26 neighbors of the central
 * point (0 otherwise).
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask3D_26neig(MultidimArray< int >& mask, int value1 = 1, int value2 = 1,
                   int value3 = 1, int center = NO_ACTIVATE);
//@}

/** Parameters for a general Mask.
 *
 * This class contains all parameters needed to generate masks. The class can
 * read parameters from the command line.
 *
 * To read a mask from a file within a program do the following
 *
 * @code
 * Mask_Params Mask;
 * Mask.type = READ_MASK;
 * Mask.fn_mask = "...";
 * Mask.generate_2Dmask();
 *
 * Mask.apply(input_Matrix2D, output_Matrix2D);
 * @endcode
 *
 * To generate a geometric mask within a program do the following:
 *
 * @code
 * Mask mask;
 *
 * // Define an spherical mask of radius 32 (the active part is
 * // within the sphere)
 * mask.type = BINARY_CIRCULAR_MASK;
 * mask.mode = INNER_MASK;
 * mask.R1 = 32;
 *
 * // resize the mask after this pattern
 * mask.resize(input_Matrix2D);
 *
 * // Really generate the mask. It is stored internally
 * mask.generate_2Dmask();
 *
 * // Apply the mask to some image
 * mask.apply_mask(input_Matrix2D, output_Matrix2D);
 * @endcode
 *
 * To read a mask from the command line:
 *
 * @code
 * Mask mask;
 * mask.read(argc, argv);
 * mask.resize(Ydim, Xdim);
 * mask.generate_2Dmask();
 * Mask.apply_mask(input_Matrix2D, output_Matrix2D);
 * @endcode
 */

class Mask
{
public:

#define NO_MASK                   0
#define BINARY_CIRCULAR_MASK      1
#define BINARY_CROWN_MASK         2
#define BINARY_CYLINDER_MASK      3
#define BINARY_FRAME_MASK         4
#define GAUSSIAN_MASK             5
#define RAISED_COSINE_MASK        6
#define BLACKMAN_MASK             7
#define SINC_MASK                 8
#define SINC_BLACKMAN_MASK        9
#define READ_MASK                10
#define RAISED_CROWN_MASK        11
#define BINARY_DWT_CIRCULAR_MASK 12
#define BINARY_DWT_SPHERICAL_MASK 13
#define BINARY_CONE_MASK         14
#define BINARY_WEDGE_MASK        15
#define BLOB_CIRCULAR_MASK       16
#define BLOB_CROWN_MASK          17
#define BINARY_TUBE              18

#define INT_MASK    1
#define DOUBLE_MASK 2
#define ALL_KINDS   INT_MASK | DOUBLE_MASK

    static void defineParams(XmippProgram * program, int allowed_data_types = ALL_KINDS,
                             const char* prefix=NULL, const char* comment=NULL, bool moreOptions=false);
    void readParams(XmippProgram * program);


    /** Mask Type
     *
     * The only valid types are BINARY_CIRCULAR_MASK, BINARY_CROWN_MASK,
     * BINARY_CYLINDER_MASK, BINARY_FRAME_MASK, GAUSSIAN_MASK,
     * RAISED_COSINE_MASK, BLACKMAN_MASK, SINC_MASK, SINC_BLACKMAN_MASK,
     * READ_MASK, RAISED_CROWN_MASK, BINARY_CONE_MASK, BINARY_WEDGE_MASK,
     * BLOB_CIRCULAR_MASK, BLOB_CROWN_MASK, BINARY_TUBE
     */
    int type;

    /** Mode
     * The valid modes are INNER_MASK and OUTSIDE_MASK.
     */
    int mode;

    /** Radius 1
     * Radius for Circular and Cylinder masks and R1 for crowns and raised
     * cosine.
     */
    double R1;

    /** Radius 2
     * R2 for crowns and raised cosine.
     */
    double R2;

    /** Pixel width
     * For raised crowns.
     */
    double pix_width;

    /** Blob parameters
     * For blob_circular and blob_crown masks.
     */
    int blob_order;
    
    /** Blob parameters */
    double blob_radius, blob_alpha;

    /** Height
     * Height for cylinders.
     */
    double H;

    /** Sigma
     * Sigma for gaussians.
     */
    double sigma;

    /** Omega
     * Frequency for sincs
     */
    double omega;

    /** Rectangular X dimension
     */
    int Xrect;

    /** Rectangular Y dimension
     */
    int Yrect;

    /** Rectangular Z dimension
     */
    int Zrect;

    /** Z origin */
    double z0;

    /** Y origin
     */
    double y0;

    /** X origin
     */
    double x0;

    /** Minimum scale for DWT masks
     */
    int smin;

    /** Maximum scale for DWT masks
     */
    int smax;

    /** Quadrant
     * If it is empty then all, except 000, are generated.
     */
    std::string quadrant;

    /** Filename from which the mask is read, if it is the case
     */
    FileName fn_mask;

    /** Geometrix transformation matrix for the mask
     */
    Matrix2D< double > mask_geo;

    /** Allowed data types.
     */
    int allowed_data_types;

    /** integer mask
     */
    MultidimArray< int > imask;

    /** double mask
     */
    MultidimArray< double > dmask;

    /** Mask type
         */
    std::string mask_type;


public:

    /** Constructors
     * Allowed data types are ALL_KINDS, INT_MASK and DOUBLE_MASK used with | .
     */
    Mask(int _allowed_data_type = ALL_KINDS);

    /** Clear
     */
    void clear();

    /** Read from command line
     * An exception is thrown if the read mask is not of an allowed type.
     */
    void read(int argc, const char** argv);

    /** Read parameters
         *
         */

    /** Show
     */
    void show() const;

    /** Usage
     */
    void usage() const;

    /** Save mask as an image
     */
    void write_mask(const FileName& fn);

    /** Return the type of the mask. INT_MASK, DOUBLE_MASK
     */
    int datatype()
    {
        if (type == BINARY_CIRCULAR_MASK || type == BINARY_CROWN_MASK ||
            type == BINARY_CYLINDER_MASK || type == BINARY_FRAME_MASK ||
            type == BINARY_TUBE ||
            type == NO_MASK || type == READ_MASK ||
            type == BINARY_DWT_CIRCULAR_MASK || type == BINARY_CONE_MASK)
            return INT_MASK;

        else if (type == GAUSSIAN_MASK || type == RAISED_COSINE_MASK ||
                 type == SINC_MASK || type == SINC_BLACKMAN_MASK ||
                 type == BLACKMAN_MASK || type == RAISED_CROWN_MASK ||
                 type == BINARY_WEDGE_MASK || type == BLOB_CIRCULAR_MASK ||
                 type == BLOB_CROWN_MASK)
            return DOUBLE_MASK;

        return 0;
    }

    /** Resize and set Xmipp origin
     */
    void resize(size_t Xdim);

    /** Resize and set Xmipp origin
     */
    void resize(size_t Ydim, size_t Xdim);

    /** Resize and set Xmipp origin
     */
    void resize(size_t Zdim, size_t Ydim, size_t Xdim);

    /** Resize after a pattern
     */
    template<typename T>
    void resize(const MultidimArray< T>& m)
    {
        switch (datatype())
        {
        case INT_MASK:
            imask.resizeNoCopy(m);
            break;

        case DOUBLE_MASK:
            dmask.resizeNoCopy(m);
            break;
        }
    }

    /** Generate mask for a resized signal
     * It is supposed that the image is already resized and with its logical
     * origin set.
     */
    void generate_mask(bool apply_geo = false );

    /** Generate mask for an empty signal
     */
    void generate_mask(size_t Zdim, size_t Ydim, size_t Xdim)
    {
        resize(Zdim, Ydim, Xdim);
        generate_mask();
    }

    /** Generate mask for an empty signal
     */
    void generate_mask(int Ydim, int Xdim, bool apply_geo = false)
    {
        resize(Ydim, Xdim);
        generate_mask(apply_geo);
    }

    /** Generate mask for an empty signal
     */
    void generate_mask(int Xdim)
    {
        resize(Xdim);
        generate_mask();
    }

    /** Generate mask for a signal following a pattern
     */
    template<typename T>
    void generate_mask(const MultidimArray< T >& m, bool apply_geo = false)
    {
        resize(m);
        generate_mask(apply_geo);
    }

    /** Apply mask to image
     * subs_val is the substitute value in case of binary masks
     */
    template<typename T>
    void apply_mask(const MultidimArray< T >& I, MultidimArray< T >& result,
                    T subs_val = 0, bool apply_geo = false)
    {
        switch (datatype())
        {
        case INT_MASK:
            if (apply_geo)
                apply_geo_binary_2D_mask(imask, mask_geo);

            apply_binary_mask(imask, I, result, subs_val);
            break;

        case DOUBLE_MASK:
            if (apply_geo)
                apply_geo_cont_2D_mask(dmask, mask_geo);

            apply_cont_mask(dmask, I, result);
            break;
        }
    }


    /** Produce vector from MultidimArray
     *
     * This function returns a 1D vector with all those voxels for which the
     * mask was greater than 0. If the output vector is of size 0, then it is
     * resized to the right size. Otherwise, it is assumed that it has already
     * the right size. The input volume is assumed to be of the same size as the
     * existing mask.
     */
    template<typename T>
    void produce_vector(const MultidimArray< T >& I, MultidimArray< T >& result)
    {
        // Resize the output vector
        if (XSIZE(result) == 0)
        {
            int size = 0;
            switch (datatype())
            {
            case INT_MASK:
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(imask)
                if (DIRECT_A2D_ELEM(imask, i, j) > 0)
                    size++;
                break;

            case DOUBLE_MASK:
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(dmask)
                if (DIRECT_A2D_ELEM(dmask, i, j) > 0)
                    size++;
                break;
            }
            result.initZeros(size);
        }

        int p = 0;
        switch (datatype())
        {
        case INT_MASK:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(imask)
            if (DIRECT_A2D_ELEM(imask, i, j) > 0)
                DIRECT_A1D_ELEM(result, p++) = DIRECT_A3D_ELEM(I, k, i, j);
            break;
        case DOUBLE_MASK:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(dmask)
            if (DIRECT_A2D_ELEM(dmask, i, j) > 0)
                DIRECT_A1D_ELEM(result, p++) = DIRECT_A3D_ELEM(I, k, i, j);
            break;
        }
    }

    /** Get binary mask
     */
    const MultidimArray< int >& get_binary_mask() const
    {
        return imask;
    }

    /** Get binary mask
     */
    MultidimArray< int >& get_binary_mask()
    {
        return imask;
    }

    /** Set binary mask
     */
    void set_binary_mask(MultidimArray< int >& _imask)
    {
        imask = _imask;
    }

    /** Get continuous mask
     */
    const MultidimArray< double >& get_cont_mask() const
    {
        return dmask;
    }

    /** Get continuous mask
     */
    MultidimArray< double >& get_cont_mask()
    {
        return dmask;
    }

    /** Set continuous mask
     */
    void set_cont_mask(MultidimArray< double >& _dmask)
    {
        dmask = _dmask;
    }

    /** Force to be continuous
     *
     * This function is used when you need a binary mask as a double matrix.
     */
    void force_to_be_continuous()
    {
        if (datatype() == INT_MASK)
        {
            typeCast(imask, dmask);
        }
    }

    /** Force to be binary
     *
     * This function is used when you need a double mask as a binary matrix.
     */
    void force_to_be_binary()
    {
        if (datatype() == DOUBLE_MASK)
        {
            typeCast(dmask, imask);
        }
    }

};

/** @name Mask Tools
 *
 * All Mask tools work only in the overlapping area of the given image/volume
 * and the mask in logical coordinates. Ie, if you have a mask defined from -2
 * to 2 and you apply it to an image defined from 0 to 63 then only those values
 * of the mask between 0 and 2 will be applied. The rest of the image will
 * remain untouched. This region where the mask is active within the overlapping
 * area will be called in this documentation: active area.
 */
//@{
/** Apply geometric transformation to a binary (2D) mask
 */
void apply_geo_binary_2D_mask(MultidimArray< int >& mask,
                              const Matrix2D< double >& A);

/** Apply geometric transformation to a continuous (2D) mask
 */
void apply_geo_cont_2D_mask(MultidimArray< double >& mask,
                            const Matrix2D< double >& A);

/** Compute statistics in the active area
 *
 * Only the statistics for values in the overlapping between the mask and the
 * volume for those the mask is not 0 are computed.
 */
template<typename T1, typename T>
void computeStats_within_binary_mask(const MultidimArray< T1 >& mask,
                                     const MultidimArray< T >& m, double& min_val,
                                     double& max_val,
                                     double& avg, double& stddev)
{
	SPEED_UP_tempsInt;
    double sum1 = 0;
    double sum2 = 0;
    int N = 0;

    max_val = min_val = DIRECT_A3D_ELEM(m, 0, 0, 0);

    FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(mask, m)
    {
        if (A3D_ELEM(mask, k, i, j) > 0)
        {
            N++;

            double aux=A3D_ELEM(m, k, i, j);
            // Minimum and maximum
            if (aux < min_val)
                min_val = aux;
            else if (aux > max_val)
                max_val = aux;

            // cumulative sums for average and standard deviation
            sum1 += aux;
            sum2 += aux*aux;
        }
    }

    // average and standard deviation
    avg  = sum1 / (double) N;
    if (N > 1)
        stddev = sqrt(fabs(sum2 / N - avg * avg) * N / (N - 1));
    else
        stddev = 0;
}

inline void computeStats_within_binary_mask(const MultidimArray< int >& mask,
        const MultidimArrayGeneric &m, double& min_val,
        double& max_val,
        double& avg, double& stddev)
{
#define COMPUTESTATS(type) \
computeStats_within_binary_mask(mask,*((MultidimArray<type>*)m.im),min_val,max_val,avg,stddev);

    SWITCHDATATYPE(m.datatype,COMPUTESTATS);

#undef COMPUTESTATS
}

/** Apply mask to a MultidimArray
 *
 * The volume values for which the input mask is 0 are set to 0. The input and
 * output volumes can be the same ones.
 */
template<typename T>
void apply_binary_mask(const MultidimArray< int >& mask, const MultidimArray< T >& m_in,
                       MultidimArray< T >& m_out, T subs_val = (T) 0)
{
    m_out.resize(m_in);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(m_out)
    // If in common with the mask
    if (k >= STARTINGZ(mask) && k <= FINISHINGZ(mask) &&
        i >= STARTINGY(mask) && i <= FINISHINGY(mask) &&
        j >= STARTINGX(mask) && j <= FINISHINGX(mask))
        if (A3D_ELEM(mask, k, i, j) == 0)
            A3D_ELEM(m_out, k, i, j) = subs_val;
        else
            A3D_ELEM(m_out, k, i, j) = A3D_ELEM(m_in, k, i, j);
    // It is not in common, leave the original one
    else
        A3D_ELEM(m_out, k, i, j) = A3D_ELEM(m_in, k, i, j);
}

/** Apply continuous mask to a MultidimArray
 *
 * The image is multiplied by the mask. The input and output matrices can be the
 * same ones. Only the overlapping values are affected by the mask.
 */
template<typename T>
void apply_cont_mask(const MultidimArray< double >& mask, const MultidimArray< T >& m_in,
                     MultidimArray< T >& m_out)
{
    m_out.resize(m_in);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(m_out)
    // If in common with the mask
    if ((k >= STARTINGZ(mask)) && (k <= FINISHINGZ(mask)) &&
        (i >= STARTINGY(mask)) && (i <= FINISHINGY(mask)) &&
        (j >= STARTINGX(mask) && j <= FINISHINGX(mask)))
    {
        A3D_ELEM(m_out, k, i, j) = (T)(A3D_ELEM(m_in, k, i, j) * A3D_ELEM(mask, k, i, j));
    }
    // It is not in common, leave the original one
    else
        A3D_ELEM(m_out, k, i, j) = A3D_ELEM(m_in, k, i, j);
}

/** Compute histogram inside mask within its minimum and maximum value (3D)
 *
 * Given a volume as input, this function returns the histogram of values inside
 * the mask within the minimum and maximum of the volume, in this way all the
 * values in the volume are counted. The volume can be of any numerical type
 * (short int, int, double, ...). The number of steps must always be given.
 */
template<typename T>
void compute_hist_within_binary_mask(const MultidimArray< int >& mask,
                                     MultidimArray< T >& v, Histogram1D &hist,
                                     int no_steps)
{
    T min_val, max_val;
    double avg, stddev;

    computeStats_within_binary_mask(mask, v, min_val, max_val, avg, stddev);
    compute_hist_within_binary_mask(mask, v, hist, min_val, max_val, no_steps);
}

/** Compute histogram inside mask within two values (3D)
 *
 * Given a volume as input, this function returns the histogram of the values
 * inside the mask within two values, the volume values outside this range are
 * not counted. This can be used to avoid the effect of outliers which causes a
 * "compression" in the histogram. The volume can be of any numerical type
 * (short int, int, double, ...). The number of steps must always be given.
 */
template<typename T>
void compute_hist_within_binary_mask(const MultidimArray< int >& mask,
                                     const MultidimArray< T >& v, Histogram1D& hist,
                                     T min, T max, int no_steps)
{
    SPEED_UP_tempsInt;
    hist.init(min, max, no_steps);
    FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(mask, v)
    if (A3D_ELEM(mask, k, i, j) != 0)
    {
    	double value=A3D_ELEM(v, k, i, j);
    	INSERT_VALUE(hist,value);
    }
}

#define COUNT_ABOVE 1
#define COUNT_BELOW 2
#define COUNT_BETWEEN 3

/** Count pixels/voxels with mask and above a threshold
 *
 * Those pixels within the mask with a value greater or equal than a threshold
 * are counted. This function makes a call to count_with_mask
 */
#define count_with_mask_above(mask, m, th) \
    count_with_mask(mask, m, COUNT_ABOVE, th, 0);

/** Count pixels/voxels with mask and below a threshold
 *
 * Those pixels within the mask with a value smaller or equal than a threshold
 * are counted. This function makes a call to count_with_mask
 */
#define count_with_mask_below(mask, m, th) \
    count_with_mask(mask, m, COUNT_BELOW, th, 0);

/** Count pixels/voxels with mask and between two thresholds
 *
 * Those pixels within the mask with a value greater or equal than th1 and
 * smaller or equal than th2 are counted. This function makes a call to
 * count_with_mask
 */
#define count_with_mask_between(mask, m, th1, th2) \
    count_with_mask(mask, m, COUNT_BETWEEN, th1, th2);

/** Count voxels with mask and threshold.
 *
 * This function returns the number of voxels in the ACTIVE area of an volume
 * with a value:
 *
 * COUNT_ABOVE: greater or equal than th1
 * COUNT_BELOW: smaller or equal than th1
 * COUNT_BETWEEN: smaller or equal than th1 and greater or equal than th2
 *
 * For complex matrices the absolute value is compared.
 */
template<typename T>
int count_with_mask(const MultidimArray< int >& mask,
                    const MultidimArray< T >& m, int mode, double th1, double th2)
{
    SPEED_UP_tempsInt;
    int N = 0;
    FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(mask, m)
    if (A3D_ELEM(mask, k, i, j))
        switch (mode)
        {
        case (COUNT_ABOVE):
                        if (A3D_ELEM(m, k, i, j) >= th1)
                            N++;
            break;

        case (COUNT_BELOW):
                        if (A3D_ELEM(m, k, i, j) <= th1)
                            N++;
            break;

        case (COUNT_BETWEEN):
                        if (A3D_ELEM(m, k, i, j) >= th1 && A3D_ELEM(m, k, i, j) <= th2)
                            N++;
            break;
        }
    return N;
}

// Specialization for complex numbers
int count_with_mask(const MultidimArray< int >& mask,
                    const MultidimArray< std::complex< double > > & m, int mode,
                    double th1, double th2);

/** Invert binary mask.
 *
 * 0's are converted in 1's and viceversa
 */
template<typename T>
void invert_binary_mask(MultidimArray< T >& mask)
{
    T* ptr=NULL;
    unsigned long int n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(mask,n,ptr)
    *ptr = 1-(*ptr);
}

/** Range adjust within binary mask
 *
 * Make the grey values of m2 fit, in L2 sense, with those in m1. Only the
 * voxels within the mask are used to compute the linear transformation. If no
 * mask is provided then all voxels are used.
 */
void rangeAdjust_within_mask(const MultidimArray< double >* mask,
                             const MultidimArray< double >& m1,
                             MultidimArray< double >& m2);
//@}

class ProgMask: public XmippMetadataProgram
{
public:

    Mask         mask;
    FileName     fn_mask;
    int          create_mask;
    int          count_above;
    double       th_above;
    int          count_below;
    double       th_below;
    double       subs_val;
    std::string  str_subs_val;
    int          count;
    int          max_length;

    void defineParams();
    void readParams();
    void preProcess();
    void postProcess();

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
};

//@}
#endif
