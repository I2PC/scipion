/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#ifndef MASK_H
#define MASK_H

#include "matrix2d.h"
#include "matrix3d.h"
#include "histogram.h"

void apply_geo_binary_2D_mask(matrix2D< int > &mask,
                              const matrix2D< double >& A);
void apply_geo_cont_2D_mask(matrix2D< double >& mask,
                            const matrix2D< double >& A);

/// @defgroup Masks Masks

/// @defgroup Masks1D 1D masks
/// @ingroup Masks

#define INNER_MASK 1
#define OUTSIDE_MASK 2
#define NO_ACTIVATE 0
#define ACTIVATE 1

/** Creates a 1D RaisedCosine mask for already sized masks
 * @ingroup Masks1D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal RaisedCosines, and outside masks are
 * 1 - RaisedCosine. When entering, the mask is initialiazed to 0 and then the
 * mask is created.
 */
void RaisedCosineMask(matrix1D< double >& mask,
                      double r1, double r2, int mode=INNER_MASK, double x0=0);

/** Creates a 1D RaisedCrown mask for already sized masks
 * @ingroup Masks1D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created within the
 * two the radii indicated with an extra region of <pix_width> pixels. The only
 * two valid modes are INNER_MASK (by default) or OUTSIDE_MASK. Inner mask are
 * normal RaisedCrowns, and outside masks are 1 - RaisedCrowns. When entering,
 * the mask is initialiazed to 0 and then the mask is created.
 */
void RaisedCrownMask(matrix1D< double >& mask,
                     double r1, double r2, double pix_width,
                     int mode=INNER_MASK,
                     double x0=0);

/// @defgroup Masks2D 2D masks
/// @ingroup Masks

/** Creates a 2D circular mask for already sized masks
 * @ingroup Masks2D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. When entering the mask is initialiazed to 0 and then the mask
 * is created.
 */
void BinaryCircularMask(matrix2D< int >& mask,
                        double radius, int mode=INNER_MASK, double x0=0, double y0=0);

/** Creates a 2D DWT circular for already sized masks
 * @ingroup Masks2D
 *
 * The mask size must be a power of 2. The radius must be expressed in pixel
 * units corresponding to the size of the image. For instance, a 64x64 image
 * might have a radius of 32 pixels to concentrate on the central part only. The
 * mask is generated only for the desired masks.
 *
 * If the quadrant="xx" then 01, 10 and 11 are generated together
 */
void BinaryDWTCircularMask(matrix2D< int >& mask,
                           double radius, int smin, int smax,
                           const string& quadrant);

/** Creates a 2D crown mask for already sized masks
 * @ingroup Masks2D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created with the two
 * radii indicated. The only two valid modes are INNER_MASK (by default, between
 * the two radii) or OUTSIDE_MASK (the negative of the crown). It is supposed
 * that R1 is smaller than R2.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinaryCrownMask(matrix2D< int >& mask,
                     double R1, double R2, int mode=INNER_MASK, double x0=0,
                     double y0=0);

/** Creates a 2D frame mask for already sized masks
 * @ingroup Masks2D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * square placed logically at (x0,y0), by default (0,0), is created with the two
 * rectangular dimensions indicated. The only two valid modes are INNER_MASK (by
 * default, between the two radii) or OUTSIDE_MASK (the negative of the crown).
 * It is supposed that R1 is smaller than R2.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinaryFrameMask(matrix2D< int >& mask,
                     int Xrect, int Yrect, int mode=INNER_MASK, double x0=0,
                     double y0=0);

/** Creates a 2D gaussian mask for already sized masks
 * @ingroup Masks2D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal gaussians, and outside masks are 1 -
 * gaussian.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void GaussianMask(matrix2D< double >& mask,
                  double sigma, int mode=INNER_MASK, double x0=0, double y0=0);

/** Creates a 2D RaisedCosine mask for already sized masks
 * @ingroup Masks2D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal RaisedCosines, and outside masks are 1 -
 * RaisedCosine.
 *
 * When entering, the mask is initialiazed to 0 and then the mask is created.
 */
void RaisedCosineMask(matrix2D< double >& mask,
                      double r1, double r2, int mode=INNER_MASK, double x0=0,
                      double y0=0);

/** Creates a 2D RaisedCrown mask for already sized masks
 * @ingroup Masks2D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created within the
 * two the radii indicated with an extra region of <pix_width> pixels. The only
 * two valid modes are INNER_MASK (by default) or OUTSIDE_MASK. Inner mask are
 * normal RaisedCrowns, and outside masks are 1 - RaisedCrowns.
 *
 * When entering, the mask is initialiazed to 0 and then the mask is created.
 */
void RaisedCrownMask(matrix2D< double >& mask,
                     double r1, double r2, double pix_width,
                     int mode=INNER_MASK, double x0=0, double y0=0);

/** 2D blackman window
 * @ingroup Masks2D
 *
 * It receives no parameter.
 */
void BlackmanMask(matrix2D< double >& mask, int mode=INNER_MASK,
                  double x0=0, double y0=0);

/** 2D Kaiser window
 * @ingroup Masks2D
 *  The mask is resized.
 *  delta=ripple (in natural units) in the pass band.
 *  Deltaw=transition bandwidth (normalized to 1.0).
 */
void KaiserMask(matrix2D<double> &mask, double delta=0.01,
                double Deltaw=1.0/12.0);

/** Creates a 2D sinc mask for already sized masks
 * @ingroup Masks2D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal gaussians, and outside masks are 1 -
 * sinc. When entering the mask is initialiazed to 0 and then the mask is
 * created.
 *
 * Remind that sinc(w*n) is zero at n=1/w;
 */
void SincMask(matrix2D< double >& mask,
              double omega, int mode=INNER_MASK, double x0=0, double y0=0);

/** Creates a 2D sinc-blackman mask, the mask is resized
 * @ingroup Masks2D
 *
 * This function returns a sinc mask windowed by a Blackman window. The window
 * is designed to cover a certain power of the sinc
 */
void SincBlackmanMask(matrix2D< double >& mask,
                      double omega, double power_percentage,
                      int mode=INNER_MASK, double x0=0, double y0=0);

/** Creates a 2D radial-sinc-kaiser mask, the mask is resized.
 * @ingroup Masks2D
 *  This function returns a sinc mask windowed by a Kaiser window.
 *  delta=ripple (in natural units) in the pass band.
 *  Deltaw=transition bandwidth (normalized to 1).
 */
void SincKaiserMask(matrix2D<double> &mask,
                    double omega, double delta=0.01, double Deltaw=1.0/12.0);

/** Creates a 2D separable-sinc-kaiser mask, the mask is resized.
 * @ingroup Masks2D
 *  This function returns a sinc mask windowed by a Kaiser window.
 *  delta=ripple (in natural units) in the pass band.
 *  Deltaw=transition bandwidth (normalized to 1).
 */
void SeparableSincKaiserMask(matrix2D<double> &mask,
                             double omega, double delta=0.01,
			     double Deltaw=1.0/12.0);

/** Creates a 3x3 mask with value (1 by default) for those 4-neighbours of the
 * central point (0 otherwise).
 * @ingroup Masks2D
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask2D_4neig(matrix2D< int >& mask, int value=1, int center=NO_ACTIVATE);

/** Creates a 3x3 mask with value1 for those 4-neighbors of the central point
 * and value2 for the 8 neighbours.
 * @ingroup Masks2D
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask2D_8neig(matrix2D< int >& mask, int value1=1, int value2=1,
                  int center=NO_ACTIVATE);

/// @defgroup Masks3D 3D masks
/// @ingroup Masks

/** Creates a 3D spherical mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be already resized and with its logical origin
 * defined. A sphere placed logically at (z0,x0,y0), by default (0,0,0), is
 * created with the radius indicated. The only two valid modes are INNER_MASK
 * (by default)or OUTSIDE_MASK.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinarySphericalMask(matrix3D< int >& mask,
                         double radius, int mode=INNER_MASK, double x0=0,
                         double y0=0, double z0=0);

/** Creates a 3D DWT spherical for already sized masks
 * @ingroup Masks3D
 *
 * The mask size must be a power of 2. The radius must be expressed in pixel
 * units corresponding to the size of the image. For instance, a 64x64x64 image
 * might have a radius of 32 pixels to concentrate on the central part only.
 *
 * If quadrant=xxx then 001,010,011,100,101,110 and 111 are generated together
 */
void BinaryDWTCircularMask(matrix3D< int >& mask,
                           double radius, int smin, int smax,
                           const string& quadrant);

/** Creates a 3D crown mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be already resized and with its logical origin
 * defined. A sphere placed logically at (z0,x0,y0), by default (0,0,0), is
 * created with the radii indicated. The only two valid modes are INNER_MASK
 * (by default)or OUTSIDE_MASK.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinaryCrownMask(matrix3D< int >& mask,
                     double R1, double R2, int mode=INNER_MASK, double x0=0,
                     double y0=0, int z0=0);

/** Creates a 3D Cylinder mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be already resized and with its logical origin
 * defined. A cylinder placed logically at (z0,x0,y0), by default (0,0,0), is
 * created with the radius and height indicated. The only two valid modes ar
 * INNER_MASK (by default) or OUTSIDE_MASK. When entering the mask is
 * initialiazed to 0 and then the mask is created.
 */
void BinaryCylinderMask(matrix3D< int >& mask,
                        double R, double H, int mode=INNER_MASK, double x0=0,
                        double y0=0, int z0=0);

/** Creates a 3D frame mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * square placed logically at (x0,y0,z0), by default (0,0,0), is created with
 * the two rectangular dimensions indicated. The only two valid modes are
 * INNER_MASK (by default, between the two radii) or OUTSIDE_MASK (the negative
 * of the crown). It is supposed that R1 is smaller than R2.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinaryFrameMask(matrix3D< int >& mask,
                     int Xrect, int Yrect, int Zrect, int mode=INNER_MASK,
                     double x0=0, double y0=0, double z0=0);

/** Creates a 3D cone mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * cone with angle theta is placed parallel to the Z-axis and centered at
 * (x0,y0,z0). The only two valid modes are INNER_MASK (by default, inside the
 * crown is zero, outside is 1) or OUTSIDE_MASK (the negative of the crown).
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 */
void BinaryConeMask(matrix3D< int >& mask,
                    double theta, int mode=INNER_MASK);

/** Creates a 3D missing wedge mask
 * @ingroup Masks3D
 *
 * The mask is supposed to be resized and with its logical origin already set.
 * theta0 and thetaF are the tilting angles (around the Y-axis) wherin between
 * the data is supposed to be collected. In this region the mask will be one,
 * outside (ie in the missing wedge) it will be zero. The mask is centered at
 * (x0,y0,z0), and rotated with respect to euler angle matrix A.
 */
void BinaryWedgeMask(matrix3D< double >& mask, double theta0, double thetaF,
                     matrix2D< double > A);

/** Creates a 3D gaussian mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0), by default (0,0), is created with the
 * radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal gaussians, and outside masks are
 *  1 - gaussian. When entering the mask is initialiazed to 0 and then the mask
 * is created.
 */
void GaussianMask(matrix3D< double >& mask,
                  double sigma, int mode=INNER_MASK, double x0=0, double y0=0,
                  double z0=0);

/** Creates a 3D RaisedCosine mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created with
 * the radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal RaisedCosines, and outside masks are 1 -
 * RaisedCosine.
 *
 * When entering, the mask is initialiazed to 0 and then the mask is created.
 */
void RaisedCosineMask(matrix3D< double >& mask,
                      double r1, double r2, int mode=INNER_MASK, double x0=0,
                      double y0=0, double z0=0);

/** Creates a 3D RaisedCrown mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created within
 * the two the radii indicated with an extra region of <pix_width> voxels. The
 * only two valid modes are INNER_MASK (by default) or OUTSIDE_MASK. Inner mask
 * are normal RaisedCrowns, and outside masks are 1 - RaisedCrowns.
 *
 * When entering, the mask is initialiazed to 0 and then the mask is created.
 */
void RaisedCrownMask(matrix3D< double >& mask,
                     double r1, double r2, double pix_width,
                     int mode=INNER_MASK, double x0=0, double y0=0,
                     double z0=0);

/** 3D blackman window
 * @ingroup Masks3D
 *
 * It receives no parameter.
 */
void BlackmanMask(matrix3D< double >& mask, int mode=INNER_MASK,
                  double x0=0, double y0=0, double z0=0);

/** Creates a 3D sinc mask for already sized masks
 * @ingroup Masks3D
 *
 * The mask is supposed to be resized and with its logical origin already set. A
 * circle placed logically at (x0,y0,z0), by default (0,0,0), is created with
 * the radius indicated. The only two valid modes are INNER_MASK (by default) or
 * OUTSIDE_MASK. Inner mask are normal gaussians, and outside masks are 1 -
 * sinc.
 *
 * When entering the mask is initialiazed to 0 and then the mask is created.
 *
 * Remind that sinc(w*t) is zero at t=1/w;
 */
void SincMask(matrix3D< double >& mask,
              double omega, int mode=INNER_MASK, double x0=0, double y0=0,
              double z0=0);

/** Creates a 3D sinc-blackman mask, the mask is resized
 * @ingroup Masks3D
 *
 * This function returns a sinc mask windowed by a Blackman window. The window
 * is designed to cover a certain power of the sinc
 */
void SincBlackmanMask(matrix3D< double >& mask,
                      double omega, double power_percentage,
                      int mode=INNER_MASK, double x0=0, double y0=0,
                      double z0=0);

/** Creates a 3x3x3 mask with value (1 by default) for those 6-neighbors of the
 * central point (0 otherwise).
 * @ingroup Masks3D
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask3D_6neig(matrix3D< int >& mask, int value=1, int center=NO_ACTIVATE);

/** Creates a 3x3x3 mask with value1 (1 by default) for those 6-neighbors and
 * value2 for the 18 neighbors of the central point (0 otherwise).
 * @ingroup Masks3D
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask3D_18neig(matrix3D< int >& mask, int value1=1, int value2=1,
                   int center=NO_ACTIVATE);

/** Creates a 3x3x3 mask with value1 (1 by default) for those 6-neighbors,
 * value2 for the 18 neighbors and value3 for the 26 neighbors of the central
 * point (0 otherwise).
 * @ingroup Masks3D
 *
 * The parameter center controls whether the center pixel is set to 1 or not
 */
void mask3D_26neig(matrix3D< int >& mask, int value1=1, int value2=1,
                   int value3=1, int center=NO_ACTIVATE);

/** Parameters for a general Mask.
 * @ingroup Masks
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
 * Mask.apply(input_matrix2D, output_matrix2D);
 * @endcode
 *
 * To generate a geometric mask within a program do the following:
 *
 * @code
 * Mask_Params Mask;
 *
 * // Define an spherical mask of radius 32 (the active part is
 * // within the sphere)
 * Mask.type = BINARY_CIRCULAR_MASK;
 * Mask.mode = INNER_MASK;
 * Mask.R1 = 32;
 *
 * // resize the mask after this pattern
 * Mask.resize(input_matrix2D);
 *
 * // Really generate the mask. It is stored internally
 * Mask.generate_2Dmask();
 *
 * // Apply the mask to some image
 * Mask.apply_mask(input_matrix2D, output_matrix2D);
 * @endcode
 *
 * To read a mask from the command line:
 *
 * @code
 * Mask_Params mask;
 * mask.read(argc, argv);
 * mask.resize(Ydim, Xdim);
 * mask.generate_2Dmask();
 * Mask.apply_mask(input_matrix2D, output_matrix2D);
 * @endcode
 */
class Mask_Params
{
public:

#define NO_MASK               	  0
#define BINARY_CIRCULAR_MASK  	  1
#define BINARY_CROWN_MASK     	  2
#define BINARY_CYLINDER_MASK  	  3
#define BINARY_FRAME_MASK     	  4
#define GAUSSIAN_MASK         	  5
#define RAISED_COSINE_MASK    	  6
#define BLACKMAN_MASK         	  7
#define SINC_MASK             	  8
#define SINC_BLACKMAN_MASK    	  9
#define READ_MASK             	 10
#define RAISED_CROWN_MASK        11
#define BINARY_DWT_CIRCULAR_MASK 12
#define BINARY_CONE_MASK     	 13
#define BINARY_WEDGE_MASK     	 14

    /** Mask Type
     *
     * The only valid types are BINARY_CIRCULAR_MASK, BINARY_CROWN_MASK,
     * BINARY_CYLINDER_MASK, BINARY_FRAME_MASK, GAUSSIAN_MASK,
     * RAISED_COSINE_MASK, BLACKMAN_MASK, SINC_MASK, SINC_BLACKMAN_MASK,
     * READ_MASK, RAISED_CROWN_MASK, BINARY_CONE_MASK, BINARY_WEDGE_MASK
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
    string quadrant;

    /** Filename from which the mask is read, if it is the case
     */
    FileName fn_mask;

    /** Geometrix transformation matrix for the mask
     */
    matrix2D< double > mask_geo;

    /** Allowed data types.
     */
    int allowed_data_types;

    /** 1D integer mask
     */
    matrix1D< int > imask1D;

    /** 2D integer mask
     */
    matrix2D< int > imask2D;

    /** 3D integer mask
     */
    matrix3D< int > imask3D;

    /** 1D double mask
     */
    matrix1D< double > dmask1D;

    /** 2D double mask
     */
    matrix2D< double > dmask2D;

    /** 3D double mask
     */
    matrix3D< double > dmask3D;

public:

#define INT_MASK    1
#define DOUBLE_MASK 2
#define ALL_KINDS   INT_MASK | DOUBLE_MASK

    /** Constructors
     * Allowed data types are ALL_KINDS, INT_MASK and DOUBLE_MASK used with | .
     */
    Mask_Params(int _allowed_data_type=ALL_KINDS);

    /** Clear
     */
    void clear();

    /** Read from command line
     * An exception is thrown if the read mask is not of an allowed type.
     */
    void read(int argc, char** argv);

    /** Show
     */
    void show() const;

    /** Usage
     */
    void usage() const;

    /** Save 1D mask as a text file
     */
    void write_1Dmask(const FileName& fn);

    /** Save 2D mask as an ImageXmipp
     */
    void write_2Dmask(const FileName& fn);

    /** Save 3D mask as an VolumeXmipp
     */
    void write_3Dmask(const FileName& fn);

    /** Return the type of the mask. INT_MASK, DOUBLE_MASK
     */
    int datatype()
    {
        if (type==BINARY_CIRCULAR_MASK || type==BINARY_CROWN_MASK ||
                type==BINARY_CYLINDER_MASK || type==BINARY_FRAME_MASK ||
                type==NO_MASK || type==READ_MASK ||
                type==BINARY_DWT_CIRCULAR_MASK || type==BINARY_CONE_MASK)
            return INT_MASK;

        else if (type==GAUSSIAN_MASK || type==RAISED_COSINE_MASK ||
                 type==SINC_MASK || type==SINC_BLACKMAN_MASK ||
                 type==BLACKMAN_MASK || type==RAISED_CROWN_MASK ||
                 type==BINARY_WEDGE_MASK)
            return DOUBLE_MASK;

        return 0;
    }

    /** Resize and set Xmipp origin
     */
    void resize(int Xdim);

    /** Resize and set Xmipp origin
     */
    void resize(int Ydim, int Xdim);

    /** Resize and set Xmipp origin
     */
    void resize(int Zdim, int Ydim, int Xdim);

    /** Resize after a pattern
     */
    template<typename T>
    void resize(const matrix1D< T>& m)
    {
        switch (datatype())
        {
        case INT_MASK:
            imask1D.resize(m);
            break;

        case DOUBLE_MASK:
            dmask1D.resize(m);
            break;
        }
    }

    /** Resize after a pattern
     */
    template<typename T>
    void resize(const matrix2D< T >& m)
    {
        switch (datatype())
        {
        case INT_MASK:
            imask2D.resize(m);
            break;

        case DOUBLE_MASK:
            dmask2D.resize(m);
            break;
        }
    }

    /** Resize after a pattern
     */
    template<typename T>
    void resize(const matrix3D< T >& m)
    {
        switch (datatype())
        {
        case INT_MASK:
            imask3D.resize(m);
            break;

        case DOUBLE_MASK:
            dmask3D.resize(m);
            break;
        }
    }

    /** Generate mask for a resized signal
     * It is supposed that the image is already resized and with its logical
     * origin set.
     */
    void generate_1Dmask();

    /** Generate mask for an empty signal
     */
    void generate_1Dmask(int Xdim)
    {
        resize(Xdim);
        generate_1Dmask();
    }

    /** Generate mask for a signal following a pattern
     */
    template<typename T>
    void generate_1Dmask(const matrix1D< T >& m)
    {
        resize(m);
        generate_1Dmask();
    }

    /** Generate mask for a resized image
     * It is supposed that the image is already resized and with its logical
     * origin set
     */
    void generate_2Dmask();

    /** Generate mask for an empty image
     */
    void generate_2Dmask(int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);
        generate_2Dmask();
    }

    /** Generate mask for an image following a pattern
     */
    template<typename T>
    void generate_2Dmask(const matrix2D< T >& m)
    {
        resize(m);
        generate_2Dmask();
    }

    /** Generate mask for a resized volume
     * It is supposed that the image is already resized and with its logical
     * origin set.
     */
    void generate_3Dmask();

    /** Generate mask for an empty volume.
     */
    void generate_3Dmask(int Zdim, int Ydim, int Xdim)
    {
        resize(Zdim, Ydim, Xdim);
        generate_3Dmask();
    }

    /** Generate mask for an image following a pattern
     */
    template<typename T>
    void generate_3Dmask(const matrix3D< T >& m)
    {
        resize(m);
        generate_3Dmask();
    }

    /** Apply mask to signal
     * subs_val is the substitute value in case of binary masks
     */
    template<typename T>
    void apply_mask(const matrix1D< T >& I, matrix1D< T >& result,
                    T subs_val=0)
    {
        switch (datatype())
        {
        case INT_MASK:
            apply_binary_mask(imask1D, I, result, subs_val);
            break;

        case DOUBLE_MASK:
            apply_cont_mask(dmask1D, I, result);
            break;
        }
    }

    /** Apply mask to image
     * subs_val is the substitute value in case of binary masks
     */
    template<typename T>
    void apply_mask(const matrix2D< T >& I, matrix2D< T >& result,
                    T subs_val=0, const bool& apply_geo=false)
    {
        switch (datatype())
        {
        case INT_MASK:
            if (apply_geo)
                apply_geo_binary_2D_mask(imask2D, mask_geo);

            apply_binary_mask(imask2D, I, result, subs_val);
            break;

        case DOUBLE_MASK:
            if (apply_geo)
                apply_geo_cont_2D_mask(dmask2D, mask_geo);

            apply_cont_mask(dmask2D, I, result);
            break;
        }
    }

    /** Apply mask to volume
     * subs_val is the substitute value in case of binary masks
     */
    template<typename T>
    void apply_mask(const matrix3D< T >& I, matrix3D< T >& result, T subs_val=0)
    {
        switch (datatype())
        {
        case INT_MASK:
            apply_binary_mask(imask3D, I, result, subs_val);
            break;

        case DOUBLE_MASK:
            apply_cont_mask(dmask3D, I, result);
            break;
        }
    }

    /** Produce vector from signal
     *
     * This function returns a 1D vector with all those points for which the
     * mask was greater than 0. If the output vector is of size 0, then it is
     * resized to the right size. Otherwise, it is assumed that it has already
     * the right size. The input vector is assumed to be of the same size as the
     * existing mask.
     */
    template<typename T>
    void produce_vector(const matrix1D< T >& I, matrix1D< T >& result)
    {
        // Resize the output vector
        if (XSIZE(result)==0)
        {
            int size=0;
            switch (datatype())
            {
            case INT_MASK:
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(imask1D)
                    if (DIRECT_VEC_ELEM(imask1D, i)>0)
                        size++;
                break;

            case DOUBLE_MASK:
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(dmask1D)
                    if (DIRECT_VEC_ELEM(dmask1D, i)>0)
                        size++;
                break;
            }
            result.init_zeros(size);
        }

        int p=0;
        switch (datatype())
        {
        case INT_MASK:
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(imask1D)
                if (DIRECT_VEC_ELEM(imask1D, i)>0)
                    DIRECT_VEC_ELEM(result, p++) = DIRECT_VEC_ELEM(I, i);
                break;

        case DOUBLE_MASK:
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(dmask1D)
                if (DIRECT_VEC_ELEM(dmask1D, i)>0)
                    DIRECT_VEC_ELEM(result, p++) = DIRECT_VEC_ELEM(I, i);
            break;
        }
    }

    /** Produce vector from image
     *
     * This function returns a 1D vector with all those pixels for which the
     * mask was greater than 0. If the output vector is of size 0, then it is
     * resized to the right size. Otherwise, it is assumed that it has already
     * the right size. The input image is assumed to be of the same size as the
     * existing mask.
     */
    template<typename T>
    void produce_vector(const matrix2D< T >& I, matrix1D< T >& result)
    {
        // Resize the output vector
        if (XSIZE(result)==0)
        {
            int size=0;
            switch (datatype())
            {
            case INT_MASK:
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(imask2D)
                    if (DIRECT_MAT_ELEM(imask2D, i, j)>0)
                        size++;
                break;

            case DOUBLE_MASK:
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(dmask2D)
                    if (DIRECT_MAT_ELEM(dmask2D, i, j)>0)
                        size++;
                break;
            }
            result.init_zeros(size);
        }

        int p=0;
        switch (datatype())
        {
        case INT_MASK:
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(imask2D)
                if (DIRECT_MAT_ELEM(imask2D, i, j)>0)
                    DIRECT_VEC_ELEM(result, p++) = DIRECT_MAT_ELEM(I, i, j);
            break;

        case DOUBLE_MASK:
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(dmask2D)
                if (DIRECT_MAT_ELEM(dmask2D, i, j)>0)
                    DIRECT_VEC_ELEM(result, p++) = DIRECT_MAT_ELEM(I, i, j);
            break;
        }
    }

    /** Produce vector from volume
     *
     * This function returns a 1D vector with all those voxels for which the
     * mask was greater than 0. If the output vector is of size 0, then it is
     * resized to the right size. Otherwise, it is assumed that it has already
     * the right size. The input volume is assumed to be of the same size as the
     * existing mask.
     */
    template<typename T>
    void produce_vector(const matrix3D< T >& I, matrix1D< T >& result)
    {
        // Resize the output vector
        if (XSIZE(result)==0)
        {
            int size=0;
            switch (datatype())
            {
            case INT_MASK:
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(imask3D)
                    if (DIRECT_MAT_ELEM(imask3D, i, j)>0)
                        size++;
                break;

            case DOUBLE_MASK:
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(dmask3D)
                    if (DIRECT_MAT_ELEM(dmask3D, i, j)>0)
                        size++;
                break;
            }
            result.init_zeros(size);
        }

        int p=0;
        switch (datatype())
        {
        case INT_MASK:
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(imask3D)
                if (DIRECT_MAT_ELEM(imask3D, i, j)>0)
                    DIRECT_VEC_ELEM(result, p++) = DIRECT_MAT_ELEM(I, i, j);
            break;
        case DOUBLE_MASK:
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(dmask3D)
                if (DIRECT_MAT_ELEM(dmask3D, i, j)>0)
                    DIRECT_VEC_ELEM(result, p++) = DIRECT_MAT_ELEM(I, i, j);
                break;
        }
    }

    /** Get binary 1D mask
     */
    matrix1D< int >& get_binary_mask1D()
    {
        return imask1D;
    }

    /** Set binary 1D mask
     */
    void set_binary_mask1D(matrix1D< int >& _imask1D)
    {
        imask1D = _imask1D;
    }

    /** Get continuous 1D mask
     */
    matrix1D< double >& get_cont_mask1D()
    {
        return dmask1D;
    }

    /** Set continuous 1D mask
     */
    void set_cont_mask1D(matrix1D< double >& _dmask1D)
    {
        dmask1D = _dmask1D;
    }

    /** Get binary 2D mask
     */
    matrix2D< int >& get_binary_mask2D()
    {
        return imask2D;
    }

    /** Set binary 2D mask
     */
    void set_binary_mask2D(matrix2D< int >& _imask2D)
    {
        imask2D = _imask2D;
    }

    /** Get continuous 2D mask
     */
    matrix2D< double >& get_cont_mask2D()
    {
        return dmask2D;
    }

    /** Set continuous 2D mask
     */
    void set_cont_mask2D(matrix2D< double >& _dmask2D)
    {
        dmask2D = _dmask2D;
    }

    /** Get binary 3D mask
     */
    matrix3D< int >& get_binary_mask3D()
    {
        return imask3D;
    }

    /** Set binary 3D mask
     */
    void set_binary_mask3D(matrix3D< int >& _imask3D)
    {
        imask3D = _imask3D;
    }

    /** Get continuous 3D mask
     */
    matrix3D< double >& get_cont_mask3D()
    {
        return dmask3D;
    }

    /** Set continuous 3D mask
     */
    void set_cont_mask3D(matrix3D< double >& _dmask3D)
    {
        dmask3D = _dmask3D;
    }

    /** Force to be continuous
     *
     * This function is used when you need a binary mask as a double matrix.
     */
    void force_to_be_continuous()
    {
        if (datatype()==INT_MASK)
        {
            type_cast(imask1D, dmask1D);
            type_cast(imask2D, dmask2D);
            type_cast(imask3D, dmask3D);
        }
    }

    /** Force to be binary
     *
     * This function is used when you need a double mask as a binary matrix.
     */
    void force_to_be_binary()
    {
        if (datatype()==DOUBLE_MASK)
        {
            type_cast(dmask1D, imask1D);
            type_cast(dmask2D, imask2D);
            type_cast(dmask3D, imask3D);
        }
    }
};

/** @defgroup MasksTools Tools
 * @ingroup Masks
 *
 * All Mask tools work only in the overlapping area of the given image/volume
 * and the mask in logical coordinates. Ie, if you have a mask defined from -2
 * to 2 and you apply it to an image defined from 0 to 63 then only those values
 * of the mask between 0 and 2 will be applied. The rest of the image will
 * remain untouched. This region where the mask is active within the overlapping
 * area will be called in this documentation: active area.
 */

/** Compute statistics in the active area (2D)
 * @ingroup MasksTools
 *
 * Only the statistics for values in the overlapping between the mask and the
 * image for those the mask is not 0 are computed.
 */
template<typename T>
void compute_stats_within_binary_mask(const matrix2D< int >& mask,
                                      const matrix2D< T >& m, T& min_val,
                                      T& max_val, double& avg, double& stddev)
{
    SPEED_UP_temps;
    double sum1=0;
    double sum2=0;
    int N=0;

    max_val = min_val = DIRECT_MAT_ELEM(m, 0, 0);

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(mask, m)
    {
        if (MAT_ELEM(mask, i, j) != 0)
        {
            N++;

            // Minimum and maximum
            if (MAT_ELEM(m, i, j) < min_val)
                min_val = MAT_ELEM(m, i, j);

            if (MAT_ELEM(m, i, j) > max_val)
                max_val = MAT_ELEM(m, i, j);

            // cumulative sums for average and standard deviation
            sum1 += (double) MAT_ELEM(m, i, j);
            sum2 += ((double) MAT_ELEM(m, i, j)) * ((double) MAT_ELEM(m, i, j));
        }
    }

    // average and standard deviation
    avg = sum1 / (double) N;
    if (N>1)
        stddev = sqrt(ABS(sum2 / N - avg * avg) * N / (N-1));
    else
        stddev = 0;
}

/** Apply geometric transformation to a binary mask
 * @ingroup MasksTools
 */
void apply_geo_binary_2D_mask(matrix2D< int >& mask,
                              const matrix2D< double >& A);

/** Apply geometric transformation to a continuous mask
 * @ingroup MasksTools
 */
void apply_geo_cont_2D_mask(matrix2D< double >& mask,
                            const matrix2D< double >& A);

/** Apply binary mask to an image (1D)
 * @ingroup MasksTools
 *
 * The image values for which the input mask is 0 are set to <subs_val>. The
 * input and output matrices can be the same ones. Only the overlapping values
 * are affected by the mask.
 *
 */
template<typename T>
void apply_binary_mask(const matrix1D< int >& mask, const matrix1D< T >& m_in,
                       matrix1D< T >& m_out, T subs_val=(T) 0)
{
    m_out.resize(m_in);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(m_out)
        // If in common with the mask
        if (i >= STARTINGX(mask) && i <= FINISHINGX(mask))
            if (VEC_ELEM(mask, i) == 0)
                VEC_ELEM(m_out, i) = subs_val;
            else
                VEC_ELEM(m_out, i) = VEC_ELEM(m_in, i);
        // It is not in common, leave the original one
        else
            VEC_ELEM(m_out, i) = VEC_ELEM(m_in, i);
}

/** Apply continuous mask to an image (1D)
 * @ingroup MasksTools
 *
 * The image is multiplied by the mask. The input and output matrices can be the
 * same ones. Only the overlapping values are affected by the mask.
 */
template<typename T>
void apply_cont_mask(const matrix1D< double >& mask, const matrix1D< T >& m_in,
                     matrix1D< T >& m_out)
{
    m_out.resize(m_in);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(m_out)
        // If in common with the mask
        if (i>=STARTINGX(mask) && i<=FINISHINGX(mask))
            VEC_ELEM(m_out, i) = (T) (VEC_ELEM(m_in, i) * VEC_ELEM(mask, i));
        // It is not in common, leave the original one
        else
            VEC_ELEM(m_out, i) = VEC_ELEM(m_in, i);
}

/** Apply binary mask to an image (2D)
 * @ingroup MasksTools
 *
 * The image values for which the input mask is 0 are set to <subs_val>. The
 * input and output matrices can be the same ones. Only the overlapping values
 * are affected by the mask
 */
template<typename T>
void apply_binary_mask(const matrix2D< int >& mask, const matrix2D< T >& m_in,
                       matrix2D< T >& m_out, T subs_val=(T) 0)
{
    m_out.resize(m_in);

    FOR_ALL_ELEMENTS_IN_MATRIX2D(m_out)
        // If in common with the mask
        if (i>=STARTINGY(mask) && i<=FINISHINGY(mask) &&
                j>=STARTINGX(mask) && j<=FINISHINGX(mask))
            if (MAT_ELEM(mask, i, j) == 0)
                MAT_ELEM(m_out, i, j) = subs_val;
            else
                MAT_ELEM(m_out, i, j) = MAT_ELEM(m_in, i, j);
        // It is not in common, leave the original one
        else
            MAT_ELEM(m_out, i, j) = MAT_ELEM(m_in, i, j);
}

/** Apply continuous mask to an image (2D)
 * @ingroup MasksTools
 *
 * The image is multiplied by the mask. The input and output matrices can be the
 * same ones. Only the overlapping values are affected by the mask.
  */
template<typename T>
void apply_cont_mask(const matrix2D< double >& mask, const matrix2D< T >& m_in,
                     matrix2D< T >& m_out)
{
    m_out.resize(m_in);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(m_out)
        // If in common with the mask
        if (i >= STARTINGY(mask) && i <= FINISHINGY(mask) &&
                j >= STARTINGX(mask) && j <= FINISHINGX(mask))
            MAT_ELEM(m_out, i, j) = (T) (MAT_ELEM(m_in, i, j) *
                MAT_ELEM(mask, i, j));
        // It is not in common, leave the original one
        else
            MAT_ELEM(m_out, i, j) = MAT_ELEM(m_in, i, j);
}

/** Compute statistics in the active area (3D)
 * @ingroup MaskTools
 *
 * Only the statistics for values in the overlapping between the mask and the
 * volume for those the mask is not 0 are computed.
 */
template<typename T>
void compute_stats_within_binary_mask(const matrix3D< int >& mask,
                                      const matrix3D< T >& m, T& min_val,
                                      T& max_val,
                                      double& avg, double& stddev)
{
    SPEED_UP_temps;
    double sum1=0;
    double sum2=0;
    int N=0;

    max_val = min_val = DIRECT_VOL_ELEM(m, 0, 0, 0);

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(mask, m)
    {
        if (VOL_ELEM(mask, k, i, j) != 0)
        {
            N++;

            // Minimum and maximum
            if (VOL_ELEM(m, k, i, j) < min_val)
                min_val = VOL_ELEM(m, k, i, j);

            if (VOL_ELEM(m, k, i, j) > max_val)
                max_val = VOL_ELEM(m, k, i, j);

            // cumulative sums for average and standard deviation
            sum1 +=  (double) VOL_ELEM(m, k, i, j);
            sum2 += ((double) VOL_ELEM(m, k, i, j)) *
                ((double) VOL_ELEM(m, k, i, j));
        }
    }

    // average and standard deviation
    avg  = sum1 / (double) N;
    if (N>1)
        stddev = sqrt(ABS(sum2 / N - avg * avg) * N/ (N-1));
    else
        stddev = 0;
}

/** Apply mask to an volume (3D)
 * @ingroup MasksTools
 *
 * The volume values for which the input mask is 0 are set to 0. The input and
 * output volumes can be the same ones.
 */
template<typename T>
void apply_binary_mask(const matrix3D< int >& mask, const matrix3D< T >& m_in,
                       matrix3D< T >& m_out, T subs_val=(T) 0)
{
    m_out.resize(m_in);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(m_out)
        // If in common with the mask
        if (k>=STARTINGZ(mask) && k<=FINISHINGZ(mask) &&
                i>=STARTINGY(mask) && i<=FINISHINGY(mask) &&
                j>=STARTINGX(mask) && j<=FINISHINGX(mask))
            if (VOL_ELEM(mask, k, i, j)==0)
                VOL_ELEM(m_out, k, i, j) = subs_val;
            else
                VOL_ELEM(m_out, k, i, j) = VOL_ELEM(m_in, k, i, j);
        // It is not in common, leave the original one
        else
            VOL_ELEM(m_out, k, i, j) = VOL_ELEM(m_in, k, i, j);
    }

/** Apply continuous mask to an image (3D)
 * @ingroup MasksTools
 *
 * The image is multiplied by the mask. The input and output matrices can be the
 * same ones. Only the overlapping values are affected by the mask.
 */
template<typename T>
void apply_cont_mask(const matrix3D< double >& mask, const matrix3D< T >& m_in,
                     matrix3D< T >& m_out)
{
    m_out.resize(m_in);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(m_out)
        // If in common with the mask
        if (k>=STARTINGZ(mask) && k<=FINISHINGZ(mask) &&
                i>=STARTINGY(mask) && i<=FINISHINGY(mask) &&
                j>=STARTINGX(mask) && j<=FINISHINGX(mask))
            VOL_ELEM(m_out, k, i, j) = (T) (VOL_ELEM(m_in, k, i, j)
                * VOL_ELEM(mask, k, i, j));
        // It is not in common, leave the original one
        else
            VOL_ELEM(m_out, k, i, j) = VOL_ELEM(m_in, k, i, j);
}

/** Compute histogram inside mask within its minimum and maximum value (2D)
 * @ingroup MasksTools
 *
 * Given a matrix as input, this function returns its histogram within the
 * minimum and maximum of the matrix inside the mask, in this way all the values
 * in the matrix are counted. The matrix can be of any numerical type (short
 * int, int, double, ...). The number of steps must always be given.
 */
template<typename T>
void compute_hist_within_binary_mask(const matrix2D< int >& mask,
                                     matrix2D< T >& v, histogram1D& hist,
                                     int no_steps)
{
    T min_val, max_val;
    double avg, stddev;

    compute_stats_within_binary_mask(mask, v, min_val, max_val, avg, stddev);
    compute_hist_within_binary_mask(mask, v, hist, min_val, max_val, no_steps);
}

/** Compute histogram inside mask within two values (2D)
 * @ingroup MasksTools
 *
 * Given a matrix as input, this function returns the histogram of values inside
 * the mask within two values, the matrix values outside this range are not
 * counted. This can be used to avoid the effect of outliers which causes a
 * "compression" in the histogram. The matrix can be of any numerical type
 * (short int, int, double, ...). The number of steps must always be given.
 */
template<typename T>
void compute_hist_within_binary_mask(const matrix2D< int >& mask,
                                     const matrix2D< T >& v, histogram1D &hist,
                                     T min, T max, int no_steps)
{
    SPEED_UP_temps;
    hist.init(min, max, no_steps);

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(mask,v)
        if (MAT_ELEM(mask, i, j) != 0)
            hist.insert_value(MAT_ELEM(v, i, j));
}

/** Compute histogram inside mask within its minimum and maximum value (3D)
 * @ingroup MasksTools
 *
 * Given a volume as input, this function returns the histogram of values inside
 * the mask within the minimum and maximum of the volume, in this way all the
 * values in the volume are counted. The volume can be of any numerical type
 * (short int, int, double, ...). The number of steps must always be given.
 */
template<typename T>
void compute_hist_within_binary_mask(const matrix3D< int >& mask,
                                     matrix3D< T >& v, histogram1D &hist,
                                     int no_steps)
{
    T min_val, max_val;
    double avg, stddev;

    compute_stats_within_binary_mask(mask, v, min_val, max_val, avg, stddev);
    compute_hist_within_binary_mask(mask, v, hist, min_val, max_val, no_steps);
}

/** Compute histogram inside mask within two values (3D)
 * @ingroup MasksTools
 *
 * Given a volume as input, this function returns the histogram of the values
 * inside the mask within two values, the volume values outside this range are
 * not counted. This can be used to avoid the effect of outliers which causes a
 * "compression" in the histogram. The volume can be of any numerical type
 * (short int, int, double, ...). The number of steps must always be given.
 */
template<typename T>
void compute_hist_within_binary_mask(const matrix3D< int >& mask,
                                     const matrix3D< T >& v, histogram1D& hist,
                                     T min, T max, int no_steps)
{
    SPEED_UP_temps;
    hist.init(min, max, no_steps);
    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(mask, v)
        if (VOL_ELEM(mask, k, i, j) != 0)
            hist.insert_value(VOL_ELEM(v, k, i, j));
}

#define COUNT_ABOVE 1
#define COUNT_BELOW 2
#define COUNT_BETWEEN 3

/** Count pixels/voxels with mask and above a threshold
 * @ingroup MasksTools
 *
 * Those pixels within the mask with a value greater or equal than a threshold
 * are counted. This function makes a call to count_with_mask
 */
#define count_with_mask_above(mask, m, th) \
    count_with_mask(mask, m, COUNT_ABOVE, th, 0);

/** Count pixels/voxels with mask and below a threshold
 * @ingroup MasksTools
 *
 * Those pixels within the mask with a value smaller or equal than a threshold
 * are counted. This function makes a call to count_with_mask
 */
#define count_with_mask_below(mask, m, th) \
    count_with_mask(mask, m, COUNT_BELOW, th, 0);

/** Count pixels/voxels with mask and between two thresholds
 * @ingroup MasksTools
 *
 * Those pixels within the mask with a value greater or equal than th1 and
 * smaller or equal than th2 are counted. This function makes a call to
 * count_with_mask
 */
#define count_with_mask_between(mask, m, th1, th2) \
    count_with_mask(mask, m, COUNT_BETWEEN, th1, th2);

/** Count pixels with mask and threshold
 * @ingroup MasksTools
 *
 * This function returns the number of pixels in the ACTIVE area of an image
 * with a value:
 *
 * COUNT_ABOVE: greater or equal than th1
 * COUNT_BELOW: smaller or equal than th1
 * COUNT_BETWEEN: smaller or equal than th1 and greater or equal than th2
 *
 * For complex matrices the absolute value is compared.
 */
template<typename T>
int count_with_mask(const matrix2D< int >& mask,
                    const matrix2D< T >& m, int mode, double th1, double th2)
{
    SPEED_UP_temps;
    int N=0;
    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(mask,m)
    if (MAT_ELEM(mask, i, j))
        switch (mode)
        {
        case (COUNT_ABOVE):
            if (MAT_ELEM(m, i, j) >= th1)
                N++;
            break;

        case (COUNT_BELOW):
            if (MAT_ELEM(m, i, j) <= th1)
                N++;
            break;

        case (COUNT_BETWEEN):
            if (MAT_ELEM(m, i, j) >= th1 && MAT_ELEM(m, i, j)<=th2)
                N++;
            break;
        }
    return N;
}

int count_with_mask(const matrix2D< int >& mask,
                    const matrix2D< complex< double > >& m, int mode,
                    double th1, double th2);

/** Count voxels with mask and threshold
 * @ingroup MasksTools
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
int count_with_mask(const matrix3D< int >& mask,
                    const matrix3D< T >& m, int mode, double th1, double th2)
{
    SPEED_UP_temps;
    int N=0;
    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(mask,m)
        if (VOL_ELEM(mask, k, i, j))
            switch (mode)
            {
            case (COUNT_ABOVE):
                if (VOL_ELEM(m, k, i, j) >= th1)
                    N++;
                break;

            case (COUNT_BELOW):
                if (VOL_ELEM(m, k, i, j) <= th1)
                    N++;
                break;

            case (COUNT_BETWEEN):
                if (VOL_ELEM(m, k, i, j) >= th1 && VOL_ELEM(m, k, i, j) <= th2)
                    N++;
                break;
            }
    return N;
}

int count_with_mask(const matrix3D< int >& mask,
                    const matrix3D< complex< double > >& m, int mode,
                    double th1, double th2);

/** Invert binary mask (2D)
 * @ingroup MasksTools
 *
 * 0's are converted in 1's and viceversa
 */
template<typename T>
void invert_binary_mask(matrix2D< T >& mask)
{
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(mask)
        MULTIDIM_ELEM(mask, i) = (MULTIDIM_ELEM(mask, i) == 1) ? 0 : 1;
}

/** Invert binary mask (3D)
 * @ingroup MasksTools
 *
 * 0's are converted in 1's and viceversa
 */
template<typename T>
void invert_binary_mask(matrix3D< T >& mask)
{
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(mask)
        MULTIDIM_ELEM(mask, i) = (MULTIDIM_ELEM(mask, i) == 1) ? 0 : 1;
}

/** Range adjust within binary mask
 * @ingroup MasksTools
 *
 * Make the grey values of m2 fit, in L2 sense, with those in m1. Only the
 * voxels within the mask are used to compute the linear transformation. If no
 * mask is provided then all voxels are used.
 */
void range_adjust_within_mask(const matrix2D< double >* mask,
                              const matrix2D< double >& m1,
                              matrix2D< double >& m2);

/** Range adjust within binary mask
 * @ingroup MasksTools
 *
 * Make the grey values of m2 fit, in L2 sense, with those in m1. Only the
 * voxels within the mask are used to compute the linear transformation. If no
 * mask is provided then all voxels are used.
 */
void range_adjust_within_mask(const matrix3D< double >* mask,
                              const matrix3D< double >& m1,
                              matrix3D< double >& m2);

#endif
