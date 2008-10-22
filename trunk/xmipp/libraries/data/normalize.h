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

#ifndef NORMALIZE_H
#define NORMALIZE_H

#include "progs.h"
#include "mask.h"

/// @defgroup Normalize Normalization of images and volumes
/// @ingroup DataLibraryPrograms

/**@defgroup NormalizationProcedures Image normalization procedures
   @ingroup Normalize
   This functions implement the normalization of a single image. They should
   be called with all images in the corresponding SelFile. In the
   following documentation m(x) is the mean of x, v(x) is the variance,
   bg(x) is its background.

   The original image is X and is supposed to be related to each single
   projection by I=a*(X+n)+b where a and b are different for every projection.

   Noise is assumed to follow a gaussian distribution N(0,sqrt(v(n)))

   In general the background mask is used only to compute statistics, while
   the mask is the one which is really applied to the image. Supply NULL if
   you don't want any mask to be applied

   When b is used
   it is measured as the mean in the background, while a*sqrt(v(n)) is the
   standard deviation in the same area.
   */
//@{
/** OldXmipp normalization.
   Formula:
   @code
   I'=(I-m(I))/sqrt(v(I))
   @endcode
   Properties:
   @code
   m(I')=0                             m(bg(I'))=-m(x)/sqrt((v(X)+v(n)))
   v(I')=1                             v(bg(I'))=v(n)/(v(X)+v(n))
   @endcode
   Comments: it's not bad but positivity constraints cannot be imposed
*/
void normalize_OldXmipp(Matrix2D<double> &I);

/** Near_OldXmipp normalization.
   Formula:
   @code
   I'=(I-m(I))/a*v(n)
   @endcode
   Properties:
   @code
   m(I')=0                             m(bg(I'))=-m(x)/sqrt(v(n))
   v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
   @endcode
   Comments: it's not bad but positivity constraints cannot be imposed
*/
void normalize_Near_OldXmipp(Matrix2D<double> &I, const Matrix2D<int> &bg_mask);

/** OldXmipp decomposition.
   Formula:
   @code
   I'=(I-b)/a*sqrt(v(n))
   I''=I'*mask
   I'''=(I''-m(I''))/sqrt(v(I''))
   @endcode
   Properties:
   @code
   m(I')=m(X)/sqrt(v(n))               m(bg(I'))=0
   v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
   @endcode
   Comments: it's not bad but positivity constraints cannot be imposed.
      If no mask is applied, then this formula is a beautiful decomposition
      of the OldXmipp method in two steps.
*/
void normalize_OldXmipp_decomposition(Matrix2D<double> &I,
                                      const Matrix2D<int> &bg_mask, const Matrix2D<double> *mask = NULL);

/** Tomography normalization.
   This is similar to the OldXmipp normalization, but the mean and
   standard deviation of the images are computed only within a region
   determined by the tilt angle.
   Formula:
   @code
   I'=(I-m(I))/(sqrt(v(I))*cos(tilt))
   @endcode
*/
void normalize_tomography(Matrix2D<double> &I, double tilt);

/** Michael's normalization.
   Formula:
   @code
   I'=(I-b)/b
   @endcode
   Properties:
   @code
   m(I')=0                             m(bg(I'))=-a*m(x)/b
   v(I')=a^2*(v(X)+v(n))/b^2           v(bg(I'))=a^2*v(n)/b^2
   @endcode
   Comments: it's not bad but positivity constraints cannot be imposed and
      the statistical properties are not so good.
*/
void normalize_Michael(Matrix2D<double> &I, const Matrix2D<int> &bg_mask);

/** NewXmipp's normalization.
   Formula:
   @code
   I'=(I-b)/a*sqrt(v(n))
   // I''=(I'>0)? I':0
   // I'''=I''-a*sqrt(v(n)/2*PI)
   I''''=I'''*mask
   @endcode
   Properties:
   @code
   m(I')=m(X)/sqrt(v(n))               m(bg(I'))=0
   v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
   @endcode
   Comments: In general, we cannot assure that mass projects into positive
      numbers, so the "denoising" capability directly on the images is
disabled. However, a positivity constraint can be applied on the 3D
volume.
*/
void normalize_NewXmipp(Matrix2D<double> &I, const Matrix2D<int> &bg_mask);

/** NewXmipp 2's normalization.
   Formula:
   @code
   I'=(I-m(bg))/(m(I)-m(bg))
   @endcode
   Properties:
   @code
   m(I')=m(X)/sqrt(v(n))               m(bg(I'))=0
   v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
   @endcode
   Comments: In general, we cannot assure that mass projects into positive
      numbers, so the "denoising" capability directly on the images is
disabled. However, a positivity constraint can be applied on the 3D
volume.
*/
void normalize_NewXmipp2(Matrix2D<double> &I, const Matrix2D<int> &bg_mask);

/** Removal of inclined background densities (ramps):
    fitting of a least squares plane throught the pixels in the
    bg_mask, then subtraction of the plane, and division by the
    standard deviation of the pixels in the bg_mask */
void normalize_ramp(Matrix2D<double> &I, const Matrix2D<int> &bg_mask);

/** Removal of neighbouring particles
    .... */
void normalize_remove_neighbours(Matrix2D<double> &I, 
				 const Matrix2D<int> &bg_mask,
                                 const double &threshold);
//@}

/** Normalize parameters
 * @ingroup Normalize
 */
class Normalize_parameters: public Prog_parameters
{
public:

// FIXME Make this an enum or similar
#define NONE 0
#define OLDXMIPP 1
#define NEAR_OLDXMIPP 2
#define NEWXMIPP 3
#define MICHAEL 4
#define NEWXMIPP2 5
#define RANDOM 6
#define RAMP 7
#define NEIGHBOUR 8
#define TOMOGRAPHY 9

    /** Normalizing method.
     * Valid methods are OLDXMIPP, NEAR_OLDXMIPP, NEWXMIPP, NEWXMIPP2, MICHAEL,
     * NONE, RANDOM, RAMP, NEIGHBOUR, TOMOGRAPHY.
     */
    int method;

    /** Nomalization of volumes.
     */
    bool volume;

// TODO Same thing, an enum
#define NONE 0
#define FRAME 1
#define CIRCLE 2

    /** Background mode.
     * Valid modes are NONE, FRAME, CIRCLE.
     */
    int background_mode;

    /** Frame width or circle radius.
     */
    int r;

    /** Upper limit of a in y=ax+b.
     */
    double aF;

    /** Lower limit of a in y=ax+b.
     */
    double a0;

    /** Upper limit of b in y=ax+b.
     */
    double bF;

    /** Lower limit of b in y=ax+b.
     */
    double b0;

    /** Flag for inverting contrast
     */
    bool invert_contrast;

    /** Flags for remving balck/white spots due to dust.
     */
    bool remove_black_dust;
    bool remove_white_dust;

    /** Threshold for removing black/white (dust) spots.
     */
    double thresh_black_dust;
    double thresh_white_dust;

    /** Sigma threshold for neighbour removal.
     */
    double thresh_neigh;

    Matrix2D< int  > bg_mask, bg_mask_bck;
    bool apply_geo;
    bool enable_mask;
    Mask_Params mask_prm;

    /** Read parameters from command line.
     */
    void read(int argc, char** argv);

    /** Produce side information.
     */
    void produce_side_info();

    /** Show parameters.
     */
    void show();

    /** Usage.
     */
    void usage();

    /** Apply inverse geometric transformation.
     * As stored in image header to the mask.
     */
    void apply_geo_mask(ImageXmipp& img);

    /** Apply to an image.
     * The input image is modified.
     */
    void apply(ImageXmipp &img);
};

#endif
