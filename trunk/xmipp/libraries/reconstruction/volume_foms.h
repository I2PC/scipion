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
#ifndef _volume_FOMs_HH
#define _volume_FOMs_HH

#include <data/volume.h>
#include <data/image.h>

#include "phantom.h"

/**@defgroup FiguresOfMerit Figures of Merit
   @ingroup ReconsLibrary
   Local FOM functions return also a vector with the local error or
   whatever for each feature. Features are numbered from 0 to
   phantom.FeatNo()-1.
*/
//@{
/** Returns the number of voxels in all features and background.
    The volume label must have for each voxel to which feature it belongs
    to. It is supposed that there is no intersection between features and
    the labelling code must be:
    @code
              0         --> background
              1         --> feature 1
              2         --> feature 2
             ...
       phantom.FeatNo() --> last feature
    @endcode
    feat_voxels contains the result with indexes 0 ... FeatNo().
    See also Phantom::label. */
void compute_voxels_in_feat(Volume *vol_label,
                            Matrix1D<double> &feat_voxels);

/** Show voxel values in feat.
    This function shows in two columns the values of the voxels corresponding
    to a feature in the phantom and in the reconstruction. If the selected
    feature is negative then the voxels for that feature is presented, and
    if it is 0 then the background voxels are shown */
void show_voxels_in_feat(const Volume *vol_phantom,
                         const Volume *vol_recons, const Volume *vol_label,
                         int selected_feat, std::ostream &out);

/** Show voxel values in feat using a phantom.
    This function shows in two columns the values of the voxels corresponding
    to a feature in the phantom and in the reconstruction. The advantage
    with respect to the previous function is that now the relative coordinate
    of the point with respecto to the center of the feature is shown */
void show_voxels_in_feat(const Volume *vol_phantom,
                         const Volume *vol_recons, const Volume *volume_label, const Phantom &P,
                         int selected_feat, std::ostream &out);

/** Compute structural consistency for a single feature.
    Under this function several measures are performed over the volume.

    The scL2_FOM is a squared error measure and it is defined as
    @code
                        1   1
    scL2_FOM(f) = 1 - ---- sum 1/2(p(i)-r(i))²
                      N(f) N(f)
    @endcode

    The scL1_FOM is an absolute error measure defined as
    @code
                        1   1
    scL1_FOM(f) = 1 - ---- sum 1/2|p(i)-r(i)|
                      N(f) N(f)
    @endcode

    The scmu_FOM measures the difference between the mean value within
    the feature in the phantom and in the reconstruction.
    @code
                        1
    scmu_FOM(f) = 1 - ---- | p(f).avg() - r(f).avg() |
                        2
    @endcode

    The scdev_FOM measures the difference between the standard deviation
    value within the feature in the phantom and in the reconstruction.
    @code
    scdev_FOM(f) = 1 - | p(f).stddev() - r(f).stddev() |
    @endcode

    The scrange_FOM is the difference between the ranges in the phantom
    and the reconstruction.
    @code
    scrange_FOM(f) = 1 - 1/2(|p(f).max()-r(f).max()| + |p(f).min()-r(f).min()|)
    @endcode

    The sccorr_FOM is the correlation between the phantom and the reconstruction.
    @code
                      1   1                                   1         1
    sccorr_FOM(f) = ---- sum (p(i)-p(f).avg())(r(i)-r(f).avg())/(sum (p(i)-p(f).avg())²*sum (r(i)-r(f).avg())²)
    @endcode  N(f) N(f)      N(f)    N(f)

    The mutual information is defined as
    @code
                        I(x,y)
    sstd::cinf_FOM(x,y) = -----------      I(x) = -sum(p(i)*log(p(i)))
                      I(x)*I(y)
    @endcode

    These calculations will be carried out only on those voxels which are
    labelled as 'f'. If f=-1 the they are carried out over the whole volume.
    If a voxel is outside the mask, then it is not taken into account.

    Activate the tell flag (TRUE or FALSE) if you want the routine to show
    all data known about the structure */
void compute_sc_FOMs(
    const Volume *vol_phantom, const Volume *vol_recons,
    const Volume *vol_label, const Volume *vol_mask, int sel_feat,
    double &scL2_FOM, double &scL1_FOM, double &scmu_FOM, double &scdev_FOM,
    double &scrange_FOM, double &sccorr_FOM, double &scinf_FOM, bool tell = false);

/** Compute histogram based FOMs.
    These FOMs try to measure the separability between the histograms
    of voxels totally inside the features and voxels in a local background
    but totally outside them.

    There are several:
    \\ mean separability:   defined for all features
    \\ border separability: defined for all features
    \\ detectability error: defined for all features
    \\ vertical resolution: only defined for double cylinders

    If a voxel is outside the mask, then it is not taken into account.

    \\========================================================================
    \\The vertical resolution FOM is defined as
    @code
                       r1.avg()+r2.avg()-2*r3.avg()
                   ----------------------------------
                   sqrt(r1.var()+r2.var()+4*r3.var())
    hsvr_FOM = ----------------------------------------
                       p1.avg()+p2.avg()-2*p3.avg()
                   ----------------------------------
                   sqrt(p1.var()+p2.var()+4*p3.var())
    @endcode
    where var() stands for variance, and r1 is a plane in the reconstruction
    in the lower cyilinder, r2 is a plane in the reconstruction in the
    upper cylinder and r3 is in between both cylinders. p1, p2 and p3 are
    the same planes but in the phantom.

    Only those voxels within the planes which intersect with at least
    one corner the phantom are taken into account.

    The returned vector (hsvr_FOMs), as always, starts at 0 (background),
    1 (feature 1),
    ... num_feat, and keeps the hsvr_FOM values for each object. If the
    object is not a double cylinder then the hsvr_FOM is -1.

    The final hsvr_FOM is the average of the hsvr_FOM for all the existing
    double cylinders.
    \\========================================================================
    \\The mean separability FOM measures if the histogram of the foreground
    (feature) has got a different mean from the one of the histogram of
    the local background. The routine studies the confidence in that
    both means are different, and returns the difference of confidence
    between the phantom and the reconstruction.
    @code
    hsmu_FOM(f)=1-|confidence(phantom)-confidence(recons)|

    confidence(recons)=erf(zr/sqrt(2))

    erf(x)=(2 / sqrt(pi)) integral{0 to x} of (e ** -t **2)  dt

    zr=ABS(mrf-mrb)/sqrt(vrf/Nf+vrb/Nb)
    @endcode
    where mrf and mrb are the mean of the reconstruction in the foreground
    and in the background respectively, vrf and vrb are the variances
    in the same places, and Nf and Nb are the number of voxels in
    the foreground and background. Only the inner foreground (defined
    as a feature with the same shape but with a factor scale of 1/3) is
    considered for the foreground statistics.

    I have modified the functional such that what is returned is the
    difference between the two z's.
    @code
    hsmu_FOM(f)=z(phantom)-z(recons)
    @endcode
    In the case that z(phantom)=infinite then hsmu_FOM(f)=z(recons);
    \\The final hsmu_FOM is the average of all the feature hsmu_FOMs
    \\========================================================================
    \\The border separability measures the separability between the
    histograms within and outside the feature near the border. For doing so,
    the mean separability between the background and the outer foreground
    (defined as a feature with the same shape in the region between a scale
    factor of 2/3 and 1).
    @code
    zr=ABS(mro-mrb)/sqrt(vro/No+vrb/Nb)

    hsbr_FOM(f)=z(phantom)-z(recons)
    @endcode
    In the case that z(phantom)=infinite then hsmu_FOM(f)=z(recons);
    \\========================================================================
    \\The detectability error measures which is the intersection area between
    the two histograms in the phantom and the reconstruction. The FOM is
    defined as
    @code
    hsdt_FOM(f) = 1- |dt_err(phantom)-dt_err(recons)|;
    @endcode
    where the dterr is the mentioned interction area.

    The final hsdt_FOM is the average of all the feature hsdt_FOMs
    \\========================================================================
    \\There is an old hot spot detectability FOM which is not computed now.
    @code
               |feat.avg() -  back.avg()|
    z = ------------------------------------
               feat.var()     back.var()
         sqrt(------------ + ------------)
                   Nb             Nf
    @endcode
    he this value was thresholded and all features above the threshold
    gave 1 point and features below the threshold 0 points.

    The value of z is shown in the tell mode.
    \\========================================================================
    \\Set tell flag to force the routine show information about the process.
    If tell!=0 then the information is shown, if furthermore the 4th bit
    (0x8) of the flag is set, then the histograms of the back and foreground
    are shown. If the histograms are to be shown then you may give
    a filename to store them in a GNUplot plottable file or "stdout" to
    show them in the standard output.
*/
void compute_hs_FOMs(Volume *vol_phantom,
                     Volume *vol_recons, const Volume *vol_label, const Volume *vol_mask,
                     int sel_feat, const Phantom &P,
                     int background_mode, double back_param,
                     double &hsmu_FOM, double &hsbr_FOM, double &hsdt_FOM, double &hsvr_FOM,
                     int tell, const FileName &fn_histog = "stdout");

/** Compute directional FOMs.
    The directional FOMs are defined along the line specified by rot and
    psi (1st and 2nd Euler angles respectively).

    The result is stored as a Matrix1D<double>

    The directional FOMs are computed in two steps:
    @code
       1) The volume is rotated after the two Euler angles
       2) With the plane at this position, the FOM is measured.
    @endcode
    The rotation of the volume gives another volume of the same dimensions.

    The two FOMs considered here are the Radon and the Slice histogram FOMs.

    If a voxel is outside the mask, then it is not taken into account.

    ========================================================================
    \\The Radon Transform is defined as the sum of the voxels in the
    Z=constant
    planes of the volume in step 2. The Radon Transform (divided by Xdim*Ydim,
    that is the number of voxels contributing to the sum) is returned
    in a vector which is defined from STARTINGZ(Vol) to FINISHINGZ(vol).

    ========================================================================
    \\The slice histograms are the histograms for the slices of the
    reconstruction at step 2. The histograms of the phantom are not
    calculated as they are assumed to be deltas.

    Histograms are usually dominated by the background value, so you
    can specify a value as a threshold, every value below this threshold
    will be considered to belong to the background and will not be
    taken into account. So the range for these histograms is from
    max_background to the maximum value within the volume divided into
    no_steps intervals.

    The resulting image has got a row of values (the histogram of that
    slice) for each slice, and rows are numbered after the height of the
    plane with respect to the center of the volume (be sure that the center
    of the volume has been calculated (and the volume stores this information)
    before entering this function), so you can find positive and negative
    indexes of rows, and within each row values are numbered as 0, 1, 2, ...,
    no_steps-1. Where the first value corresponds to the histogram interval
    nearest the background threshold and the last one corresponds to the
    interval nearest the maximum value within the volume.

    If the Radon Transforms are to be shown then you may give
    a filename to store them in a GNUplot plottable file or "stdout" to
    show them in the standard output.
*/
void compute_dr_FOMs(const Volume *vol_phantom, const Volume *vol_recons,
                     const Volume *vol_mask,
                     double rot, double tilt,
                     Image *img_histog, int no_steps, double max_background,
                     double &drrt_FOM, int tell, const FileName &fn_radon = "stdout");

/** Distance map for a volume.
    The distance map is another volume which expresses the distance of
    a voxel to the nearest feature (if the voxel belongs completely
    to the background) or to the the background (if the voxel belongs
    completely to a feature). If the voxel is in the border of a feature
    then the distance is 0.
    The distance is measure as the euclidean distance between the two
    voxel centers.

    If a voxel is outside the mask, then it is not taken into account.
    Its distance is -1.
*/
void compute_distance_map(const Volume *vol_label, const Phantom &label,
                          const Volume *vol_mask, Volume *vol_distance);

/** Compute distance based FOMs.
    Based on the distance map you can do the following measures:
    \\- Blurring distance
    \\- Appearance distance

    If a voxel is outside the mask, then it is not taken into account.
    These voxels are labelled with a distance=-1.

    ========================================================================
    \\ The blurring distance somehow reflects up to which distance the
    blurring extends. It is defined as in the following formula:
    @code
                     -1     1    Nb   1    1
    blurring_distance   = ----- sum ----- --- (pi-ri)²
                           Nb   i=1  di    2
    @endcode
    where the sum is extended only to those voxels of the background
    (there are Nb of these voxels). di is the distance value for that
    voxel (see compute_distance_map to know how this distance is
    defined), and pi and ri are the phantom and reconstruction values
    at that voxel respectively.

    Notice that the blurring distance is returned and not its inverse.

    ========================================================================
    \\ The appearance distance try to show up the new masses appearing in the
    reconstruction which were not in the phantom, for this reason
    gives a lot of importance to errors far from the phantom. The formula
    is
    @code
                            1    Nb     1
    appearance_distance = ----- sum di --- (pi-ri)²
                           Nb   i=1     2
    @endcode
*/
void compute_ds_FOMs(const Volume *vol_phantom, const Volume *vol_recons,
                     const Volume *vol_label,
                     const Volume *vol_distance, double &dsbl_FOM, double &dsad_FOM);

/** Show shape.
    The shape is defined as the density as a function of the minimum
    distance to the background or a feature. Previously you MUST have
    calculated the distance map.
    The shape is shown ususally in the stdout but you may give
    a filename to store them in a GNUplot plottable file.

    If a voxel is outside the mask, then it is not taken into account.
    These voxels are labelled with a distance -1.
*/
void show_shape(const Volume *vol_phantom, const Volume *vol_recons,
                const Volume *vol_label,
                const Volume *vol_distance, const FileName &fn_shape);


/* Compute resolution.
   This function is performed via Spider. If spider is not found then -1
   is returned. If the flag SHOW_PROCESS of tell is set the file with
   the resolution analysis is shown.

   An exception is thrown if SPIDER environment variable cannot be found.
*/
/* Compute resolution.
   This function is performed via Bsoft. If there is any problem calling
   the program or reading the output file then -1 is returned.
*/
void compute_resolution(VolumeXmipp &vol_phantom,
                        VolumeXmipp &vol_recons, double &resolution);

/* Compute Fourier Shell Correlation.
   This function is performed via Bsoft. If there is any problem calling
   the program or reading the output file then empty vectors (frequency
   (in 1/Angstroms),FSC) are returned. This function returns the first
   frequency at which the FSC falls below 0.5.
*/
double compute_FSC(VolumeXmipp &vol_phantom,
                   VolumeXmipp &vol_recons, double sampling_rate,
                   Matrix1D<double> &frequency, Matrix1D<double> &FSC);
//@}
#endif
