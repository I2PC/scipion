/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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
#ifndef _PROG_ANGULAR_PREDICT
#define _PROG_ANGULAR_PREDICT

#include <data/xmipp_funcs.h>
#include <data/xmipp_program.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image.h>
#include <classification/pca.h>

#include "angular_distance.h"
#include <data/symmetries.h>

#include <map>
#include <algorithm>

/**@defgroup AngularPredict angular_discrete_assign (Discrete angular assignment)
   @ingroup ReconsLibrary */
//@{
/** Angular Predict parameters. */
class ProgAngularDiscreteAssign: public XmippMetadataProgram
{
public:
    /** Selfile with the reference images */
    FileName fn_ref;
    /** Filename for the symmetry file */
    FileName fn_sym;
    /** Maximum projection change.
        -1 means all are allowed. */
    double max_proj_change;
    /** Maximum psi change.
        -1 means all are allowed */
    double max_psi_change;
    /** Psi step */
    double psi_step;
    /** Maximum shift change */
    double max_shift_change;
    /** Shift step */
    double shift_step;
    /** Threshold for discarding images.
        If it is 20%, only the 20% of the images will be kept each round. */
    double th_discard;
    /** Minimum scale. Finest level */
    int smin;
    /** Maximum scale. Coarsest level */
    int smax;
    /** Way to pick views.
        0 maximum correlation of the first group.
        1 average of the most populated group.
        2 maximum correlation of the most populated group. */
    int pick;
#define TELL_ROT_TILT  1
#define TELL_PSI_SHIFT 2
#define TELL_OPTIONS   4
    /** Show level.*/
    int tell;
    /** Extended usage */
    bool extended_usage;
    /** 5D search, instead of 3D+2D */
    bool search5D;
public:
    // Number of subbands
    int SBNo;
    // Subband size
    Matrix1D<int> SBsize;
    // Selfile with the reference projections
    MetaData SF_ref;
    // Mask disitribution of DWT coefficients.
    // It is created when the training sets
    MultidimArray<int> Mask_no;
    // Vector with all the DWT coefficients of the
    // library
    std::vector<MultidimArray<double> * > library;
    // Vector with all the names of the library images
    std::vector<FileName> library_name;
    // Power of the library images at different
    // subbands
    MultidimArray<double> library_power;
    // Vector with the rotational angles of the library
    std::vector<double> rot;
    // Vector with the tilting angles of the library
    std::vector<double> tilt;
    // Parameters for computing distances
    ProgAngularDistance distance_prm;
public:
    /// Empty constructor
    ProgAngularDiscreteAssign();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /** Produce side info. */
    void preProcess();

    /** Write output metadata */
    void postProcess();

    /** Produce library.*/
    void produce_library();

    /** Build candidate list.
        Build a candidate list with all possible reference projections
        which are not further than the maximum allowed change from
        the given image.

        The list size is the total number of reference images. For each
        image the list is true if it is still a candidate.*/
    void build_ref_candidate_list(const Image<double> &I,
                                  bool * candidate_list, std::vector<double> &cumulative_corr,
                                  std::vector<double> &sumxy);

    /** Refine candidate list via correlation. Given a projection image and the
        list of alive candidates, this function correlates the input image with
        all alive candidates and leave to pass only th% of the images.

        m is the subband being studied*/
    void refine_candidate_list_with_correlation(int m,
            Matrix1D<double> &dwt, bool * candidate_list,
            std::vector<double> &cumulative_corr,
            Matrix1D<double> &x_power,
            std::vector<double> &sumxy, double th = 50);

    /** Evaluate candidates by correlation. The evaluation is returned in
        candidate_rate. Furthermore, this function returns the threshold for
        passing in the "score exam", a 7*/
    double evaluate_candidates(const std::vector<double> &vscore,
                               const std::vector<int> &candidate_idx, std::vector<double> &candidate_rate,
                               double weight);

    /** Group views.
        The input images are supposed to be ordered by rate.
        The groups are also sorted by rate. */
    void group_views(const std::vector<double> &vrot,
                     const std::vector<double> &vtilt, const std::vector<double> &vpsi,
                     const std::vector<int> &best_idx, const std::vector<int> &candidate_idx,
                     std::vector< std::vector<int> > &groups);

    /** Pick the best image from the groups.
        If method == 0 it takes the maximum of the first group (the
        one with best rate). If method==1, it takes the maximum
        of the most populated group. */
    int pick_view(int method,
                  std::vector< std::vector<int> > &groups,
                  std::vector<double> &vscore,
                  std::vector<double> &vrot,
                  std::vector<double> &vtilt,
                  std::vector<double> &vpsi,
                  const std::vector<int> &best_idx,
                  const std::vector<int> &candidate_idx, const std::vector<double> &candidate_rates);

    /** Predict rotational and tilting angles.
        The function returns the two assigned angles and the corresponding
        correlation. The index of the best matching reference image is also
        returned. The function predict shift and psi angle calls this
        one for evaluating each possible combination.*/
    double predict_rot_tilt_angles(Image<double> &I,
                                   double &assigned_rot, double &assigned_tilt, int &best_ref_idx);

    /** Process one image.
        Predict angles and shift.
        This function searches in the shift-psi space and for each combination
        it correlates with the whole reference set. */
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /** Finish processing.
        Close all output files. */
//    void postProcess();
};
//@}
#endif
