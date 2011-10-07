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

#ifndef NORMALIZE_H
#define NORMALIZE_H

#include "xmipp_program.h"
#include "xmipp_image.h"
#include "mask.h"

/// @defgroup Normalize Normalization of images and volumes
/// @ingroup DataLibrary
//@{
/**@name Image normalization procedures
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
    @ingroup NormalizationProcedures
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
void normalize_OldXmipp(MultidimArray<double> &I);

/** Near_OldXmipp normalization.
    @ingroup NormalizationProcedures
   Formula:
   @code
   I'=(I-m(I))/sqrt(v(bg))
   @endcode
   Properties:
   @code
   m(I')=0                             m(bg(I'))=-m(x)/sqrt(v(n))
   v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
   @endcode
   Comments: it's not bad but positivity constraints cannot be imposed
*/
void normalize_Near_OldXmipp(MultidimArray<double> &I, const MultidimArray<int> &bg_mask);

/** OldXmipp decomposition.
    @ingroup NormalizationProcedures
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
void normalize_OldXmipp_decomposition(MultidimArray<double> &I,
                                      const MultidimArray<int> &bg_mask, const MultidimArray<double> *mask = NULL);

/** Tomography normalization.
    @ingroup NormalizationProcedures
   This is similar to the OldXmipp normalization, but the mean and
   standard deviation of the images are computed only within a region
   determined by the tilt angle.
   Formula for tomography:
   @code
   I'=(I-m(I))/(sqrt(v(I))*cos^2(tilt))
   @endcode

   Formula for tomography0:
   @code
   I'=(I/cos(tilt)-m0)/(sqrt(v(I))*cos(tilt))
   @endcode

   The estimated mean of the image and the local variance are returned in
   sigmai and mui.
*/
void normalize_tomography(MultidimArray<double> &I, double tilt, double &mui,
                          double &sigmai, bool tiltMask,
                          bool tomography0=false, double mu0=0, double sigma0=1);

/** Michael's normalization.
    @ingroup NormalizationProcedures
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
void normalize_Michael(MultidimArray<double> &I, const MultidimArray<int> &bg_mask);

/** NewXmipp's normalization.
    @ingroup NormalizationProcedures
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
void normalize_NewXmipp(MultidimArray<double> &I, const MultidimArray<int> &bg_mask);

/** NewXmipp 2's normalization.
    @ingroup NormalizationProcedures
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
void normalize_NewXmipp2(MultidimArray<double> &I, const MultidimArray<int> &bg_mask);

/** Removal of inclined background densities (ramps).
    @ingroup NormalizationProcedures
    fitting of a least squares plane throught the pixels in the
    bg_mask, then subtraction of the plane, and division by the
    standard deviation of the pixels in the bg_mask */
void normalize_ramp(MultidimArray<double> &I, MultidimArray<int> &bg_mask);

/** Removal of neighbouring particles.
    @ingroup NormalizationProcedures
    .... */
void normalize_remove_neighbours(MultidimArray<double> &I,
                                 const MultidimArray<int> &bg_mask,
                                 const double &threshold);
//@}

/* Normalize program
 */
class ProgNormalize : public XmippMetadataProgram
{
protected:
    /** Possible normalization modes */
    enum NormalizationMode {
        NONE,
        OLDXMIPP,
        NEAR_OLDXMIPP,
        NEWXMIPP,
        MICHAEL,
        NEWXMIPP2,
        RANDOM,
        RAMP,
        NEIGHBOUR,
        TOMOGRAPHY,
        TOMOGRAPHY0
    };

    /** Normalizing method.
     * Valid methods are OLDXMIPP, NEAR_OLDXMIPP, NEWXMIPP, NEWXMIPP2, MICHAEL,
     * NONE, RANDOM, RAMP, NEIGHBOUR, TOMOGRAPHY, TOMOGRAPHY0.
     */
    NormalizationMode method;

    /** Normalization of volumes.
     */
    bool volume;

    /** Possible background modes */
    enum BackgroundMode {
        NOBACKGROUND,
        FRAME,
        CIRCLE
    };

    /** Background mode.
     * Valid modes are NONE, FRAME, CIRCLE.
     */
    BackgroundMode background_mode;

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

    /** Flag for applying a mask depending on the tilt
     */
    bool tiltMask;

    /** Flags for removing balck/white spots due to dust.
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

    MultidimArray<int> bg_mask, bg_mask_bck;
    bool enable_mask;

    /* Mask parameter
     */
    Mask mask_prm;

public:
    // Mean and standard deviation of the image 0. Used for tomography
    double mu0, sigma0;

protected:
    void defineParams();
    void readParams();
    void show();
    void preProcess();
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
};
//@}
#endif
