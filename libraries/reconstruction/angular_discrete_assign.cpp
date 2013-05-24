/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2002)
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

#include "angular_discrete_assign.h"
#include "fourier_filter.h"

#include <data/args.h>
#include <data/histogram.h>
#include <data/geometry.h>
#include <data/wavelet.h>
#include <data/mask.h>
#include <data/filters.h>
#include <algorithm>

// Empty constructor =======================================================
ProgAngularDiscreteAssign::ProgAngularDiscreteAssign()
{
    produces_a_metadata = true;
    produces_an_output = true;
}

// Read arguments ==========================================================
void ProgAngularDiscreteAssign::readParams()
{
    XmippMetadataProgram::readParams();
    fn_ref = getParam("--ref");
    fn_sym = getParam("--sym");
    max_proj_change = getDoubleParam("--max_proj_change");
    max_psi_change = getDoubleParam("--max_psi_change");
    max_shift_change = getDoubleParam("--max_shift_change");
    psi_step = getDoubleParam("--psi_step");
    shift_step = getDoubleParam("--shift_step");
    th_discard = getDoubleParam("--keep");
    smin = getIntParam("--smin");
    smax = getIntParam("--smax");
    pick = getIntParam("--pick");
    tell = 0;
    if (checkParam("--show_rot_tilt"))
        tell |= TELL_ROT_TILT;
    if (checkParam("--show_psi_shift"))
        tell |= TELL_PSI_SHIFT;
    if (checkParam("--show_options"))
        tell |= TELL_OPTIONS;
    search5D = checkParam("--search5D");
}

// Show ====================================================================
void ProgAngularDiscreteAssign::show()
{
    if (!verbose)
        return;
    XmippMetadataProgram::show();
    std::cout << "Reference images: " << fn_ref << std::endl
    << "Max proj change: " << max_proj_change << std::endl
    << "Max psi change: " << max_psi_change << " step: " << psi_step << std::endl
    << "Max shift change: " << max_shift_change << " step: " << shift_step << std::endl
    << "Keep %: " << th_discard << std::endl
    << "smin: " << smin << std::endl
    << "smax: " << smax << std::endl
    << "Pick: " << pick << std::endl
    << "Show level: " << tell << std::endl
    << "5D search: " << search5D << std::endl
    ;
}

// usage ===================================================================
void ProgAngularDiscreteAssign::defineParams()
{
    addUsageLine("Make a projection assignment using wavelets on a discrete library of projections");
    addUsageLine("+This program assigns Euler angles to experimental projections by matching ");
    addUsageLine("+with ideal projections. This matching is done via a DWT correlation. For ");
    addUsageLine("+every experimental projection, different in-plane rotations and shifts are ");
    addUsageLine("+tried (by exhaustive search). For each possible combination of these two ");
    addUsageLine("+variables, the best correlating ideal projection is sought using a fast ");
    addUsageLine("+multirresolution algorithm.");
    addUsageLine("+The method is fully described at http://www.ncbi.nlm.nih.gov/pubmed/15099579");
    defaultComments["-i"].clear();
    defaultComments["-i"].addComment("List of images to align");
    defaultComments["-i"].addComment("+ Alignment parameters can be provided ");
    defaultComments["-i"].addComment("+ Only the shifts are taken in consideration ");
    defaultComments["-i"].addComment("+ in global searches; in local searches, all ");
    defaultComments["-i"].addComment("+ parameters in the initial docfile are considered.");
    defaultComments["-o"].clear();
    defaultComments["-o"].addComment("Metadata with output alignment");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --ref <selfile>             : Metadata with the reference images and their angles");
    addParamsLine("                               :+Must be created with [[angular_project_library_v3][angular_project_library]]");
    addParamsLine("  [--sym <symmetry_file=\"\">] : Symmetry file if any");
    addParamsLine("                               :+The definition of the symmetry is described at [[transform_symmetrize_v3][transform_symmetrize]]");
    addParamsLine("  [--max_shift_change <r=0>]   : Maximum change allowed in shift");
    addParamsLine("  [--psi_step <ang=5>]         : Step in psi in degrees");
    addParamsLine("  [--shift_step <r=1>]         : Step in shift in pixels");
    addParamsLine("==+Extra parameters==");
    addParamsLine("  [--search5D]                 : Perform a 5D search instead of 3D+2D");
    addParamsLine("  [--max_proj_change <ang=-1>] : Maximum change allowed in rot-tilt");
    addParamsLine("  [--max_psi_change <ang=-1>]  : Maximum change allowed in psi");
    addParamsLine("  [--keep <th=50>]             : How many images are kept each round (%)");
    addParamsLine("  [--smin <s=1>]               : Finest scale to consider (lowest value=0)");
    addParamsLine("  [--smax <s=-1>]              : Coarsest scale to consider (highest value=log2(Xdim))");
    addParamsLine("  [--pick <mth=1>]             : 0 --> maximum of the first group");
    addParamsLine("                               : 1 --> maximum of the most populated");
    addParamsLine("  [--show_rot_tilt]            : Show the rot-tilt process");
    addParamsLine("  [--show_psi_shift]           : Show the psi-shift process");
    addParamsLine("  [--show_options]             : Show final options among which");
    addParamsLine("                               : the angles are selected");
    addExampleLine("Typical use:",false);
    addExampleLine("xmipp_angular_project_library -i referenceVolume.vol -o reference.stk --sampling_rate 5");
    addExampleLine("xmipp_angular_discrete_assign -i projections.sel -o discrete_assignment.xmd --ref reference.doc");
}

// Produce side information ================================================
void ProgAngularDiscreteAssign::preProcess()
{
    // Read input reference image names
    SF_ref.read(fn_ref);
    size_t refYdim, refXdim, refZdim, refNdim;
    getImageSize(SF_ref,refYdim, refXdim, refZdim, refNdim);
    if (refYdim != NEXT_POWER_OF_2(refYdim) || refXdim != NEXT_POWER_OF_2(refXdim))
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     "reference images must be of a size that is power of 2");

    // Produce side info of the angular distance computer
    distance_prm.fn_ang1 = distance_prm.fn_ang2 = "";
    distance_prm.fn_sym = fn_sym;
    distance_prm.produce_side_info();

    // Read the angle file
    rot.resize(SF_ref.size());
    tilt.resize(SF_ref.size());
    int i = 0;
    FOR_ALL_OBJECTS_IN_METADATA(SF_ref)
    {
        SF_ref.getValue(MDL_ANGLE_ROT, rot[i],__iter.objId);
        SF_ref.getValue(MDL_ANGLE_TILT, tilt[i],__iter.objId);
        i++;
    }

    // Build mask for subbands
    Mask_no.resize(refYdim, refXdim);
    Mask_no.initConstant(-1);

    if (smax == -1)
        smax = Get_Max_Scale(refYdim) - 3;
    SBNo = (smax - smin + 1) * 3 + 1;
    SBsize.resize(SBNo);

    Mask Mask(INT_MASK);
    Mask.type = BINARY_DWT_CIRCULAR_MASK;
    Mask.R1 = CEIL((double)refXdim / 2.0);
    Mask.resize(refYdim, refXdim);

    int m = 0, s;
    for (s = smax; s >= smin; s--)
    {
        for (int q = 0; q <= 3; q++)
        {
            if (q == 0 && s != smax)
                continue;
            Mask.smin = s;
            Mask.smax = s;
            Mask.quadrant = Quadrant2D(q);
            Mask.generate_mask();

            const MultidimArray<int> imask2D=Mask.get_binary_mask();
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imask2D)
            if (DIRECT_A2D_ELEM(imask2D, i, j))
            {
                Mask_no(i, j) = m;
                SBsize(m)++;
            }

            m++;
        }
    }

    // Produce library
    produce_library();

    // Save a little space
    SF_ref.clear();
}

// PostProcess ---------------------------------------------------------------
void ProgAngularDiscreteAssign::postProcess()
{
	if (single_image)
		getOutputMd()->write(fn_out);
}

// Produce library -----------------------------------------------------------
void ProgAngularDiscreteAssign::produce_library()
{
    Image<double> I;
    int number_of_imgs = SF_ref.size();
    set_DWT_type(DAUB12);

    // Create space for all the DWT coefficients of the library
    Matrix1D<int> SBidx(SBNo);
    for (int m = 0; m < SBNo; m++)
    {
        MultidimArray<double> *subband = new MultidimArray<double>;
        subband->resize(number_of_imgs, SBsize(m));
        library.push_back(subband);
    }
    library_power.initZeros(number_of_imgs, SBNo);

    if (verbose)
    {
        std::cerr << "Generating reference library ...\n";
        init_progress_bar(number_of_imgs);
    }
    int n = 0, nstep = XMIPP_MAX(1, number_of_imgs / 60); // For progress bar
    FOR_ALL_OBJECTS_IN_METADATA(SF_ref)
    {
        I.readApplyGeo(SF_ref,__iter.objId);
        library_name.push_back(I.name());

        // Make and distribute its DWT coefficients in the different PCA bins
        I().statisticsAdjust(0, 1);
        DWT(I(), I());
        SBidx.initZeros();
        FOR_ALL_ELEMENTS_IN_ARRAY2D(Mask_no)
        {
            int m = Mask_no(i, j);
            if (m != -1)
            {
                double coef = I(i, j), coef2 = coef * coef;
                (*library[m])(n, SBidx(m)++) = coef;
                for (int mp = m; mp < SBNo; mp++)
                    library_power(n, mp) += coef2;
            }
        }

        // Prepare for next iteration
        if (++n % nstep == 0 && verbose)
            progress_bar(n);
    }
    if (verbose)
        progress_bar(SF_ref.size());
}

// Build candidate list ------------------------------------------------------
void ProgAngularDiscreteAssign::build_ref_candidate_list(const Image<double> &I,
        bool *candidate_list, std::vector<double> &cumulative_corr,
        std::vector<double> &sumxy)
{
    int refNo = rot.size();
    cumulative_corr.resize(refNo);
    sumxy.resize(refNo);
    for (int i = 0; i < refNo; i++)
    {
        candidate_list[i] = true;
        sumxy[i] = cumulative_corr[i] = 0;
        if (max_proj_change != -1)
        {
            double dummy_rot = rot[i], dummy_tilt = tilt[i], dummy_psi;
            double ang_distance = distance_prm.SL.computeDistance(
                                      I.rot(), I.tilt(), 0,
                                      dummy_rot, dummy_tilt, dummy_psi,
                                      true, false, false);
            candidate_list[i] = (ang_distance <= max_proj_change);
#ifdef DEBUG

            std::cout << "(" << I.rot() << "," << I.tilt() << ") and ("
            << rot[i] << "," << tilt[i] << ") --> " << ang_distance << std::endl;
#endif

        }
    }
}

// Refine candidate list ---------------------------------------------------
void ProgAngularDiscreteAssign::refine_candidate_list_with_correlation(
    int m,
    Matrix1D<double> &dwt,
    bool *candidate_list, std::vector<double> &cumulative_corr,
    Matrix1D<double> &x_power, std::vector<double> &sumxy,
    double th)
{
    int dimp = SBsize(m);
    int imax = rot.size();
    const MultidimArray<double> &library_m = *(library[m]);
    std::vector<double> sortedCorr;
    sortedCorr.reserve(imax);
    for (int i = 0; i < imax; i++)
    {
        if (candidate_list[i])
        {
            double sumxyp = 0.0;
            // Loop unrolling
            unsigned long jmax=4*(dimp/4);
            for (int j = 0; j < dimp; j+=4)
            {
                int j_1=j+1;
                int j_2=j+2;
                int j_3=j+3;
                sumxyp += VEC_ELEM(dwt,j)   * DIRECT_A2D_ELEM(library_m,i, j) +
                          VEC_ELEM(dwt,j_1) * DIRECT_A2D_ELEM(library_m,i, j_1) +
                          VEC_ELEM(dwt,j_2) * DIRECT_A2D_ELEM(library_m,i, j_2) +
                          VEC_ELEM(dwt,j_3) * DIRECT_A2D_ELEM(library_m,i, j_3);
            }
            for (int j = jmax+1;j<dimp;++j)
                sumxyp += VEC_ELEM(dwt,j)   * DIRECT_A2D_ELEM(library_m,i, j);

            sumxy[i] += sumxyp;

            double corr = sumxy[i] / sqrt(DIRECT_A2D_ELEM(library_power,i, m) *
                                          VEC_ELEM(x_power,m));
            cumulative_corr[i] = corr;
            sortedCorr.push_back(corr);

            if (tell & TELL_ROT_TILT)
            {
                std::cout << "Candidate " << i << " corr= " << cumulative_corr[i]
                << " rot= " << rot[i] << " tilt= " << tilt[i] << std::endl;
            }
        }
    }
    std::sort(sortedCorr.begin(),sortedCorr.end());
    int idx=(int)floor(sortedCorr.size()*(1-th/100.0));

    double corr_th = sortedCorr[idx];

    // Remove all those projections below the threshold
    for (int i = 0; i < imax; i++)
        if (candidate_list[i])
            candidate_list[i] = (cumulative_corr[i] >= corr_th);

    // Show the percentil used
    if (tell & TELL_ROT_TILT)
        std::cout << "# Percentil " << corr_th << std::endl << std::endl;
}

// Predict rot and tilt ----------------------------------------------------
double ProgAngularDiscreteAssign::predict_rot_tilt_angles(Image<double> &I,
        double &assigned_rot, double &assigned_tilt, int &best_ref_idx)
{
    if (XSIZE(I()) != NEXT_POWER_OF_2(XSIZE(I())) ||
        YSIZE(I()) != NEXT_POWER_OF_2(YSIZE(I())))
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     "experimental images must be of a size that is power of 2");

    // Build initial candidate list
    bool* candidate_list=new bool[rot.size()];
    std::vector<double> cumulative_corr;
    std::vector<double> sumxy;
    build_ref_candidate_list(I, candidate_list, cumulative_corr, sumxy);

    // Make DWT of the input image and build vectors for comparison
    std::vector<Matrix1D<double> * > Idwt;
    Matrix1D<double> x_power(SBNo);
    x_power.initZeros();
    Matrix1D<int> SBidx(SBNo);
    SBidx.initZeros();
    for (int m = 0; m < SBNo; m++)
    {
        Matrix1D<double> *subband = new Matrix1D<double>;
        subband->resize(SBsize(m));
        Idwt.push_back(subband);
    }

    I().statisticsAdjust(0, 1);
    DWT(I(), I());
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mask_no)
    {
        int m = DIRECT_A2D_ELEM(Mask_no,i, j);
        if (m != -1)
        {
            double coef = DIRECT_A2D_ELEM(IMGMATRIX(I), i, j), coef2 = coef * coef;
            (*(Idwt[m]))(SBidx(m)++) = coef;
            for (int mp = m; mp < SBNo; mp++)
                VEC_ELEM(x_power,mp) += coef2;
        }
    }

    // Measure correlation for all possible PCAs
    // These variables are used to compute the correlation at a certain
    // scale
    for (int m = 0; m < SBNo; m++)
    {
        // Show image name
        if (tell & TELL_ROT_TILT)
            std::cout << "# " << I.name() << " m=" << m
            << " current rot="  << I.rot()
            << " current tilt=" << I.tilt() << std::endl;
        refine_candidate_list_with_correlation(m, *(Idwt[m]),
                                               candidate_list, cumulative_corr,
                                               x_power, sumxy, th_discard);
    }

    // Select the maximum
    int best_i = -1;
    bool first = true;
    int N_max = 0;
    int imax = rot.size();
    for (int i = 0; i < imax; i++)
        if (candidate_list[i])
        {
            if (first)
            {
                best_i = i;
                first = false;
                N_max = 1;
            }
            else if (cumulative_corr[i] > cumulative_corr[best_i])
                best_i = i;
            else if (cumulative_corr[i] == cumulative_corr[best_i])
                N_max++;
        }

    if (N_max == 0)
    {
        if (verbose)
            std::cerr << "Predict_angles: Empty candidate list for image "
            << I.name() << std::endl;
        assigned_rot = I.rot();
        assigned_tilt = I.tilt();
        return 0;
    }

    // There are several maxima, choose one randomly
    if (N_max != 1)
    {
        int selected = FLOOR(rnd_unif(0, 3));
        for (int i = 0; i < imax && selected >= 0; i++)
            if (candidate_list[i])
                if (cumulative_corr[i] == cumulative_corr[best_i])
                {
                    best_i = i;
                    selected--;
                }
    }

    // Free asked memory
    for (int m = 0; m < SBNo; m++)
        delete Idwt[m];

    assigned_rot    = rot[best_i];
    assigned_tilt   = tilt[best_i];
    best_ref_idx    = best_i;
    delete [] candidate_list;
    return cumulative_corr[best_i];
}

// Evaluate candidates by correlation ----------------------------------------
double ProgAngularDiscreteAssign::evaluate_candidates(
    const std::vector<double> &vscore,
    const std::vector<int> &candidate_idx,
    std::vector<double> &candidate_rate, double weight)
{
    // Compute maximum and minimum of correlations
    int imax = vscore.size();
    double min_score, max_score;
    min_score = max_score = vscore[0];
    for (int i = 1; i < imax; i++)
    {
        if (vscore[i] < min_score)
            min_score = vscore [i];
        if (vscore[i] > max_score)
            max_score = vscore [i];
    }

    // Divide the correlation segment in as many pieces as candidates
    double score_step = (max_score - min_score) / 10;

    int jmax = candidate_idx.size();
    for (int j = 0; j < jmax; j++)
    {
        int i = candidate_idx[j];
        int points;
        if (score_step != 0)
            points = FLOOR((vscore[i] - min_score) / score_step);
        else
            points = 10;
        if (tell & TELL_PSI_SHIFT)
            std::cout << "Candidate (" << i << ") score=" << vscore[i]
            << " points=" << points << std::endl;
        candidate_rate[j] += weight * points;
    }

    if (tell & TELL_PSI_SHIFT)
        std::cout << "Evaluation:" << candidate_rate << std::endl
        << "Threshold for obtaining a 7 in score: "
        << min_score + 7*score_step << std::endl;
    return min_score + 7*score_step;
}

// Group images --------------------------------------------------------------
//#define DEBUG
void ProgAngularDiscreteAssign::group_views(const std::vector<double> &vrot,
        const std::vector<double> &vtilt, const std::vector<double> &vpsi,
        const std::vector<int> &best_idx, const std::vector<int> &candidate_idx,
        std::vector< std::vector<int> > &groups)
{
    for (size_t j = 0; j < best_idx.size(); j++)
    {
        int i = candidate_idx[best_idx[j]];
#ifdef DEBUG

        std::cout << "Looking for a group for image " << best_idx[j] << std::endl;
#endif

        double roti = vrot[i];
        double tilti = vtilt[i];
        double psii = vpsi[i];
        // See if there is any suitable group
        bool assigned = false;
        size_t g;
        for (g = 0; g < groups.size(); g++)
        {
            bool fits = true;
            for (size_t jp = 0; jp < groups[g].size(); jp++)
            {
                int ip = candidate_idx[groups[g][jp]];
                double ang_distance = distance_prm.SL.computeDistance(
                                          vrot[ip], vtilt[ip], vpsi[ip],
                                          roti, tilti, psii, false, false, false);
#ifdef DEBUG

                std::cout << "   comparing with " << groups[g][jp] << " d="
                << ang_distance << std::endl;
#endif

                if (ang_distance > 15)
                {
                    fits = false;
                    break;
                }
            }
            if (fits)
            {
                assigned = true;
                break;
            }
        }

        if (!assigned)
        {
#ifdef DEBUG
            std::cout << "Creating a new group\n";
#endif
            // Create a new group with the first item in the list
            std::vector<int> group;
            group.push_back(best_idx[j]);
            groups.push_back(group);
        }
        else
        {
#ifdef DEBUG
            std::cout << "Assigning to group " << g << std::endl;
#endif
            // Insert the image in the fitting group
            groups[g].push_back(best_idx[j]);
        }
    }

    // Check if there are so many groups as images.
    // If that is the case, then everything is a single group
    // if denoising is not used
    if (groups.size() == best_idx.size())
    {
        groups.clear();
        std::vector<int> group;
        for (size_t j = 0; j < best_idx.size(); j++)
            group.push_back(best_idx[j]);
        groups.push_back(group);
    }
}
#undef DEBUG

// Pick view -----------------------------------------------------------------
int ProgAngularDiscreteAssign::pick_view(int method,
        std::vector< std::vector<int> > &groups,
        std::vector<double> &vscore,
        std::vector<double> &vrot,
        std::vector<double> &vtilt,
        std::vector<double> &vpsi,
        const std::vector<int> &best_idx,
        const std::vector<int> &candidate_idx, const std::vector<double> &candidate_rate)
{
    if (method == 0)
    {
        // This one returns the most scored image of the first group
        double best_rate = -1e38;
        int    best_j, jmax = groups[0].size();
        for (int j = 0; j < jmax; j++)
            // Select the best with the scoreelation
            if (vscore[candidate_idx[groups[0][j]]] > best_rate)
            {
                best_j = j;
                best_rate = vscore[candidate_idx[groups[0][j]]];
            }
        return groups[0][best_j];
    }
    else if (method == 1)
    {
        // Sum the rates in all groups
        std::vector<double> group_rate;
        group_rate.reserve(groups.size());
        int best_g;
        double best_group_rate = -1e38;
        for (size_t g = 0; g < groups.size(); g++)
        {
            double temp_group_rate = 0;
            for (size_t j = 0; j < groups[g].size(); j++)
                temp_group_rate += candidate_rate[groups[g][j]];
            group_rate.push_back(temp_group_rate);
            if (temp_group_rate > best_group_rate)
            {
                best_g = g;
                best_group_rate = group_rate[g];
            }
        }

        // Check that there are not two groups with the same score
        int groups_with_max_rate = 0;
        for (size_t g = 0; g < groups.size(); g++)
            if (group_rate[g] == best_group_rate)
                groups_with_max_rate++;
#ifdef DEBUG

        if (groups_with_max_rate > 1)
        {
            std::cout << "There are two groups with maximum rate\n";
        }
#endif

        // Take the best image within that group
        int best_j;
        double best_rate = -1e38;
        for (size_t j = 0; j < groups[best_g].size(); j++)
        {
#ifdef NEVER_DEFINED
            // Select the best with the rate
            if (candidate_rate[groups[best_g][j]] > best_rate)
            {
                best_j = j;
                best_rate = candidate_rate[groups[best_g][j]];
            }
#endif
            // Select the best with the scoreelation
            if (vscore[candidate_idx[groups[best_g][j]]] > best_rate)
            {
                best_j = j;
                best_rate = vscore[candidate_idx[groups[best_g][j]]];
            }
        }

        // Check that there are not two images with the same rate
        int images_with_max_rate = 0;
        for (size_t j = 0; j < groups[best_g].size(); j++)
#ifdef NEVER_DEFINED
            // Select the best with the rate
            if (candidate_rate[groups[best_g][j]] == best_rate)
                images_with_max_rate++;
#endif
        // Select the best with scoreelation
        if (vscore[candidate_idx[groups[best_g][j]]] == best_rate)
            images_with_max_rate++;
        if (images_with_max_rate > 1)
        {
            // If there are several with the same punctuation take the one
            // with the best scoreelation
            double best_score = -1e38;
            for (size_t j = 0; j < groups[best_g].size(); j++)
            {
                if (vscore[candidate_idx[groups[best_g][j]]] > best_score &&
                    candidate_rate[groups[best_g][j]] == best_rate)
                {
                    best_j = j;
                    best_score = vscore[candidate_idx[groups[best_g][j]]];
                }
            }
        }
        return groups[best_g][best_j];
    }
    REPORT_ERROR(ERR_UNCLASSIFIED,"The code should not have reached this point");
}

// Run ---------------------------------------------------------------------
// Predict shift and psi ---------------------------------------------------
// #define DEBUG
void ProgAngularDiscreteAssign::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    // Read the image and take its angles from the Metadata
    // if they are available. If not, take them from the header.
    // If not, set them to 0.
    Image<double> img;
    img.readApplyGeo(fnImg, rowIn);

    double best_rot, best_tilt, best_psi, best_shiftX, best_shiftY,
    best_score = 0, best_rate;

    Image<double> Ip;
    Ip = img;
    Matrix1D<double> shift(2);

    // Get the 2D alignment shift
    double Xoff = img.Xoff();
    double Yoff = img.Yoff();

    // Establish psi limits
    double psi0, psiF;
    if (max_psi_change == -1)
    {
        psi0 = -180;
        psiF = 180 - psi_step;
    }
    else
    {
        psi0 = img.psi() - max_psi_change;
        psiF = img.psi() + max_psi_change;
    }
    double R2 = max_shift_change * max_shift_change;

    // Search in the psi-shift space
    int N_trials = 0;
    std::vector<double> vshiftX, vshiftY, vpsi, vrot, vtilt, vcorr,
    vproj_error, vproj_compact, vang_jump, vscore;
    std::vector<int>    vref_idx;

    double backup_max_shift_change = max_shift_change;
    if (!search5D)
        max_shift_change = 0;

#ifdef DEBUG

    img.write("PPPoriginal.xmp");
    double bestCorr=-1e38;
#endif

    for (int flip=0; flip<=1; ++flip)
		for (double shiftX = Xoff - max_shift_change; shiftX <= Xoff + max_shift_change; shiftX += shift_step)
			for (double shiftY = Yoff - max_shift_change; shiftY <= Yoff + max_shift_change; shiftY += shift_step)
			{
				if ((shiftX - Xoff)*(shiftX - Xoff) + (shiftY - Yoff)*(shiftY - Yoff) > R2)
					continue;
				for (double psi = psi0; psi <= psiF; psi += psi_step)
				{
					N_trials++;

					// Flip if necessary
					Ip()=img();
					if (flip)
						Ip().selfReverseX();

					// Shift image if necessary
					if (shiftX != 0 || shiftY != 0)
					{
						VECTOR_R2(shift, shiftX, shiftY);
						selfTranslate(LINEAR,Ip(),shift,WRAP);
					}

					// Rotate image if necessary
					// Adding 2 is a trick to avoid that the 0, 90, 180 and 270
					// are treated in a different way
					selfRotate(LINEAR,Ip(),psi + 2, WRAP);
					selfRotate(LINEAR,Ip(),-2, WRAP);
#ifdef DEBUG
					Image<double> Ipsave;
					Ipsave()=Ip();
#endif

					// Project the resulting image onto the visible space
					double proj_error = 0.0, proj_compact = 0.0;

					// Search for the best tilt, rot angles
					double rotp, tiltp;
					int best_ref_idx;
					double corrp =
						predict_rot_tilt_angles(Ip, rotp, tiltp, best_ref_idx);

					double aux_rot = rotp, aux_tilt = tiltp, aux_psi = psi;
					double ang_jump = distance_prm.SL.computeDistance(
										  img.rot(), img.tilt(), img.psi(),
										  aux_rot, aux_tilt, aux_psi,
										  false, false, false);

					double shiftXp=shiftX;
					double shiftYp=shiftY;
					double psip=psi;
					if (flip)
					{
						// std::cout << "       before flipping " << rotp << " " << tiltp << " " << psip << " " << shiftXp << " " << shiftYp << " " << corrp << std::endl;
						shiftXp=-shiftXp;
						double newrot, newtilt, newpsi;
						Euler_mirrorY(rotp,tiltp,psi,newrot,newtilt,newpsi);
						rotp=newrot;
						tiltp=newtilt;
						psip=newpsi;
					}

					vshiftX.push_back(shiftXp);
					vshiftY.push_back(shiftYp);
					vrot.push_back(rotp);
					vtilt.push_back(tiltp);
					vpsi.push_back(psip);
					vcorr.push_back(corrp);
					vproj_error.push_back(proj_error);
					vproj_compact.push_back(proj_compact);
					vang_jump.push_back(ang_jump);
					vref_idx.push_back(best_ref_idx);
					// std::cout << flip << " " << rotp << " " << tiltp << " " << psip << " " << shiftXp << " " << shiftYp << " " << corrp << std::endl;

	#ifdef DEBUG
					if (corrp>bestCorr)
					{
						Ipsave.write("PPPafter_denoising.xmp");
						Image<double> Iref;
						Iref.read(library_name[best_ref_idx]);
						Iref.write("PPPref.xmp");
						std::cerr << "This is index " << vcorr.size()-1 << std::endl;
						std::cerr << "corrp=" << corrp << "\nPress any key\n";
						bestCorr=corrp;
						char c;
						std::cin >> c;
					}
	#endif

				}
			}

    // Compute extrema of all scoring factors
    double max_corr        = vcorr[0],         min_corr        = vcorr[0];
    double max_proj_error  = vproj_error[0],   min_proj_error  = vproj_error[0];
    double max_proj_compact = vproj_compact[0], min_proj_compact = vproj_compact[0];
    for (int i = 1; i < N_trials; i++)
    {
        // Compute extrema of projection error
        if (vproj_error[i] < min_proj_error)
            min_proj_error = vproj_error[i];
        if (vproj_error[i] > max_proj_error)
            max_proj_error = vproj_error[i];

        // Compute extrema of correlation
        if (vcorr[i] < min_corr)
            min_corr = vcorr[i];
        if (vcorr[i] > max_corr)
            max_corr = vcorr[i];

        // Compute extrema of projection error
        if (vproj_compact[i] < min_proj_compact)
            min_proj_compact = vproj_compact[i];
        if (vproj_compact[i] > max_proj_compact)
            max_proj_compact = vproj_compact[i];
    }

    // Score each projection
    vscore.reserve(vcorr.size());
    for (int i = 0; i < N_trials; i++)
    {
        vscore.push_back(vcorr[i]);
        if (tell & TELL_PSI_SHIFT)
            std::cout << "i=" << i
            << " shiftX= " << vshiftX[i] << " shiftY= " << vshiftY[i]
            << " psi= "          << vpsi[i]
            << " rot= "          << vrot[i]
            << " tilt= "         << vtilt[i]
            << " score= "        << vscore[i]
            << " corr= "         << vcorr[i]
            << " proj_error= "   << vproj_error[i]
            << " proj_compact= " << vproj_compact[i]
            << " refidx= "       << vref_idx[i]
            << " ang_jump= "     << vang_jump[i]
            << std::endl;
    }

    // Is the psi range circular?
    bool circular = realWRAP(vpsi[0] - psi_step, -180, 180) ==
                    realWRAP(vpsi[N_trials-1], -180, 180);

    // Compute minimum and maximum of the correlation and projection error
    double avg_score_maxima = 0;
    std::vector<int> local_maxima;
    if (tell & TELL_PSI_SHIFT)
        std::cout << "Local maxima\n";
    for (int i = 0; i < N_trials; i++)
    {
        // Look for the left and right sample
        int il = i - 1, ir = i + 1;
        if (i == 0 && circular)
            il = N_trials - 1;
        else if (i == N_trials - 1)
        {
            if (circular)
                ir = 0;
            else
                ir = -1;
        }

        // Check if the error is a local minimum of the projection error
        // or a local maxima of the correlation
        bool is_local_maxima = true;
        if (il != -1)
            if (vscore[il] >= vscore[i])
                is_local_maxima = false;
        if (ir != -1)
            if (vscore[ir] >= vscore[i])
                is_local_maxima = false;

        // It is a local minimum
        if (is_local_maxima /*|| visible_space*/)
            avg_score_maxima += vscore[i];
        if (is_local_maxima /*|| visible_space*/)
        {
            local_maxima.push_back(i);
            if (tell & TELL_PSI_SHIFT)
                std::cout << "i= " << i
                << " psi= " << vpsi[i] << " rot= " << vrot[i] << " tilt= "
                << vtilt[i] << " score= " << vscore[i] << std::endl;
        }
    }
    avg_score_maxima /= local_maxima.size();
    if (tell & TELL_PSI_SHIFT)
        std::cout << "Avg_maxima=" << avg_score_maxima << std::endl;

    // Remove all local maxima below the average
    int jmax = local_maxima.size();
    std::vector<int> candidate_local_maxima;
    std::vector<double> candidate_rate;
    if (tell & TELL_PSI_SHIFT)
        std::cout << "Keeping ...\n";
    for (int j = 0; j < jmax; j++)
    {
        int i = local_maxima[j];
        if (vscore[i] >= avg_score_maxima)
        {
            candidate_local_maxima.push_back(i);
            candidate_rate.push_back(0);
            if (tell & TELL_PSI_SHIFT)
                std::cout << "i= " << i
                << " psi= " << vpsi[i] << " rot= " << vrot[i] << " tilt= "
                << vtilt[i] << " score= " << vscore[i] << std::endl;
        }
    }
    jmax = candidate_local_maxima.size();

    // Evaluate the local maxima according to their correlations
    evaluate_candidates(vscore, candidate_local_maxima, candidate_rate, 1);

    // Sort the candidates
    if (tell & TELL_PSI_SHIFT)
        std::cout << "\nSelecting image\n";
    MultidimArray<double> score(jmax);
    for (int j = 0; j < jmax; j++)
        score(j) = vscore[candidate_local_maxima[j]];
    MultidimArray<int> idx_score;
    score.indexSort(idx_score);

    if (tell & (TELL_PSI_SHIFT | TELL_OPTIONS))
    {
        std::cout << img.name() << std::endl;
        std::cout.flush();
        for (int j = 0; j < jmax; j++)
        {
            int jp = idx_score(j) - 1;
            int i = candidate_local_maxima[jp];
            std::cout << "i= " << i
            << " psi= " << vpsi[i] << " rot= " << vrot[i] << " tilt= "
            << vtilt[i]
            << " score= " << vscore[i]
            << " corr= " << vcorr[i]
            << " error= " << vproj_error[i]
            << " compact= " << vproj_compact[i]
            << " angjump= " << vang_jump[i]
            << " rate=" << candidate_rate[jp]
            << " reference image #=" << vref_idx[i] + 1 << std::endl;
        }
        std::cout << std::endl;
        std::cout.flush();
    }

    // Consider the top
    int jtop = jmax - 1;
    std::vector<int> best_idx;
    int max_score_diff = 1;
    while (jtop > 0 &&
           candidate_rate[idx_score(jmax-1)-1] -
           candidate_rate[idx_score(jtop-1)-1] <= max_score_diff)
    {
        best_idx.push_back(idx_score(jtop) - 1);
        jtop--;
    }
    best_idx.push_back(idx_score(jtop) - 1);
    if (tell & TELL_PSI_SHIFT)
        std::cout << "Best indices: " << best_idx << std::endl;

    // Pick the best one from the top
    int ibest, jbest;
    if (jtop == jmax - 1)
    {
        // There is only one on the top
        jbest = best_idx[0];
        ibest = candidate_local_maxima[jbest];
    }
    else
    {
        // There are more than one in the top
        // Group the different views
        std::vector< std::vector<int> > groups;
        group_views(vrot, vtilt, vpsi, best_idx, candidate_local_maxima, groups);
        if (tell & TELL_PSI_SHIFT)
            std::cout << "Partition: " << groups << std::endl;

        // Pick the best image from the groups
        jbest = pick_view(pick, groups, vscore, vrot, vtilt, vpsi,
                          best_idx, candidate_local_maxima, candidate_rate);
        ibest = candidate_local_maxima[jbest];
    }

    // Is it a 3D+2D search?
    if (!search5D)
    {
        Image<double> Iref;
        //Iref.readApplyGeo(library_name[vref_idx[ibest]]);
        //TODO: Check if this is correct
        Iref.read(library_name[vref_idx[ibest]]);
        Iref().setXmippOrigin();
        selfRotate(LINEAR,Iref(),-vpsi[ibest]);
        if (Xoff == 0 && Yoff == 0)
            Ip() = img();
        else
        {
            VECTOR_R2(shift, Xoff, Yoff);
            translate(LINEAR,Ip(),img(),shift,WRAP);
        }

        double shiftX, shiftY;
        CorrelationAux aux;
        bestShift(Iref(), Ip(), shiftX, shiftY, aux);
        if (shiftX*shiftX + shiftY*shiftY > R2)
        {
            shiftX = shiftY = 0;
        }
        vshiftX[ibest] = Xoff + shiftX;
        vshiftY[ibest] = Yoff + shiftY;
        max_shift_change = backup_max_shift_change;
    }

    // Take the information of the best image
    best_rot    = vrot[ibest];
    best_tilt   = vtilt[ibest];
    best_psi    = vpsi[ibest];
    best_shiftX = vshiftX[ibest];
    best_shiftY = vshiftY[ibest];
    best_score  = vscore[ibest];
    best_rate   = candidate_rate[jbest];

    if (tell & (TELL_PSI_SHIFT | TELL_OPTIONS))
    {
        std::cout << "Originally it had, psi=" << img.psi() << " rot=" << img.rot()
        << " tilt=" << img.tilt() << std::endl;
        std::cout << "Finally I choose: ";
        if (tell & TELL_PSI_SHIFT)
            std::cout << jbest << "\n";
        std::cout << "psi= " << best_psi << " rot= " << best_rot << " tilt= "
        << best_tilt << " shiftX=" << best_shiftX
        << " shiftY=" << best_shiftY << " score= " << best_score
        << " rate= " << best_rate << std::endl << std::endl;
    }

    // Save results
    rowOut.setValue(MDL_ANGLE_ROT,  best_rot);
    rowOut.setValue(MDL_ANGLE_TILT, best_tilt);
    rowOut.setValue(MDL_ANGLE_PSI,  -best_psi);
    rowOut.setValue(MDL_SHIFT_X,    best_shiftX);
    rowOut.setValue(MDL_SHIFT_Y,    best_shiftY);
    rowOut.setValue(MDL_MAXCC,     best_score);
}
#undef DEBUG

// Finish processing ---------------------------------------------------------
//void ProgAngularDiscreteAssign::postProcess()
//{
//    mdOut.write(fn_out);
//}
