/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.uam.es (2002)
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

#include "angular_discrete_assign.h"
#include "fourier_filter.h"

#include <data/args.h>
#include <data/docfile.h>
#include <data/histogram.h>
#include <data/geometry.h>
#include <data/wavelet.h>
#include <data/mask.h>
#include <data/filters.h>

// Empty constructor =======================================================
Prog_angular_predict_prm::Prog_angular_predict_prm()
{
    MPIversion = false;
    each_image_produces_an_output = true;
}

// Read arguments ==========================================================
void Prog_angular_predict_prm::read(int argc, char **argv)
{
    extended_usage = checkParameter(argc, argv, "-more_help");
    if (extended_usage) REPORT_ERROR(1, "");
    Prog_parameters::read(argc, argv);
    fn_ref = getParameter(argc, argv, "-ref");
    fn_ang = getParameter(argc, argv, "-ang", "");
    fn_out_ang = getParameter(argc, argv, "-oang");
    fn_sym = getParameter(argc, argv, "-sym", "");
    max_proj_change = textToFloat(getParameter(argc, argv, "-max_proj_change", "-1"));
    max_psi_change = textToFloat(getParameter(argc, argv, "-max_psi_change", "-1"));
    max_shift_change = textToFloat(getParameter(argc, argv, "-max_shift_change", "0"));
    psi_step = textToFloat(getParameter(argc, argv, "-psi_step", "5"));
    shift_step = textToFloat(getParameter(argc, argv, "-shift_step", "1"));
    th_discard = textToFloat(getParameter(argc, argv, "-keep", "50"));
    smin = textToInteger(getParameter(argc, argv, "-smin", "1"));
    smax = textToInteger(getParameter(argc, argv, "-smax", "-1"));
    check_mirrors = !checkParameter(argc, argv, "-do_not_check_mirrors");
    pick = textToInteger(getParameter(argc, argv, "-pick", "1"));
    dont_apply_geo = checkParameter(argc, argv, "-dont_apply_geo");
    dont_modify_header = checkParameter(argc, argv, "-dont_modify_header");
    proj_step = textToFloat(getParameter(argc, argv, "-proj_step", "5"));
    tell = 0;
    if (checkParameter(argc, argv, "-show_rot_tilt")) tell |= TELL_ROT_TILT;
    if (checkParameter(argc, argv, "-show_psi_shift")) tell |= TELL_PSI_SHIFT;
    if (checkParameter(argc, argv, "-show_options")) tell |= TELL_OPTIONS;
    search5D = checkParameter(argc, argv, "-5D");
    summaryRootname = getParameter(argc, argv, "-summary", "");
    if (!MPIversion) produce_side_info();
}

// Show ====================================================================
void Prog_angular_predict_prm::show()
{
    Prog_parameters::show();
    std::cout << "Reference images: " << fn_ref << std::endl
    << "Angle file: " << fn_ang << std::endl
    << "Ouput angular file: " << fn_out_ang << std::endl
    << "Max proj change: " << max_proj_change << std::endl
    << "Max psi change: " << max_psi_change << " step: " << psi_step << std::endl
    << "Max shift change: " << max_shift_change << " step: " << shift_step << std::endl
    << "Proj step: " << proj_step << std::endl
    << "Keep %: " << th_discard << std::endl
    << "smin: " << smin << std::endl
    << "smax: " << smax << std::endl
    << "Check mirrors: " << check_mirrors << std::endl
    << "Pick: " << pick << std::endl
    << "Dont apply geo: " << dont_apply_geo << std::endl
    << "Show level: " << tell << std::endl
    << "Modify header:  " << !dont_modify_header << std::endl
    << "5D search: " << search5D << std::endl
    << "Summary: " << summaryRootname << std::endl
    ;
}

// usage ===================================================================
void Prog_angular_predict_prm::usage()
{
    Prog_parameters::usage();
    if (extended_usage) more_usage();
    else
    {
        std::cerr << "   -ref <selfile|volume>    : Selfile with the reference images\n"
        << "                              If a volume is given, supply -proj_step\n"
        << "   -oang <angle file>       : DocFile with output angles\n"
        << "  [-sym <symmetry file>]    : Symmetry file if any\n"
        << "  [-dont_apply_geo]         : do not apply the translations in the header\n"
        << "  [-dont_modify_header]     : Don't save the parameters in the\n"
        << "                              image header\n"
        << "  [-more_help]              : Show all options\n"
        ;
    }
}

void Prog_angular_predict_prm::more_usage()
{
    std::cerr << "   -ref <selfile|volume>    : Selfile with the reference images\n"
    << "                              If a volume is given, supply -proj_step\n"
    << "  [-ang <angle file>]       : DocFile with the angles for the reference\n"
    << "                              produced by xmipp_project\n"
    << "   -oang <angle file>       : DocFile with output angles\n"
    << "  [-sym <symmetry file>]    : Symmetry file if any\n"
    << "  [-max_proj_change <ang=-1>]: Maximum change allowed in rot-tilt\n"
    << "  [-max_psi_change <ang=-1>]: Maximum change allowed in psi\n"
    << "  [-max_shift_change <r=0>] : Maximum change allowed in shift\n"
    << "  [-proj_step <ang=5>]      : Projection (rot-tilt) step\n"
    << "  [-psi_step <ang=5>]       : Step in psi in degrees\n"
    << "  [-shift_step <r=1>]       : Step in shift in pixels\n"
    << "  [-keep <th=50%>]          : How many images are kept each round\n"
    << "  [-smin <s=1>]             : Finest scale to consider (lowest value=0)\n"
    << "  [-smax <s=-1>]            : Coarsest scale to consider (highest value=log2(Xdim))\n"
    << "  [-do_not_check_mirrors]   : Otherwise, mirror versions of the experimental\n"
    << "                              images are also explored\n"
    << "  [-pick <mth=1>]           : 0 --> maximum of the first group\n"
    << "                              1 --> maximum of the most populated\n"
    << "  [-dont_apply_geo]         : do not apply the translations in the header\n"
    << "                              The header rotation is never applied\n"
    << "  [-show_rot_tilt]          : Show the rot-tilt process\n"
    << "  [-show_psi_shift]         : Show the psi-shift process\n"
    << "  [-show_options]           : Show final options among which\n"
    << "                              the angles are selected\n"
    << "  [-dont_modify_header]     : Don't save the parameters in the\n"
    << "                              image header\n"
    << "  [-5D]                     : Perform a 5D search instead of 3D+2D\n"
    << "  [-summary <rootname>]     : Summary rootname\n"
    ;
}

// Produce side information ================================================
void Prog_angular_predict_prm::produce_side_info(int rank)
{
    volume_mode = false;

    // Information for the SF_main
    allow_time_bar = (tell == 0 && rank == 0);

    // Read input reference image names
    if (Is_VolumeXmipp(fn_ref))
    {
        volume_mode = true;
        std::cerr << "Generating reference projections ...\n";

        // Generate the reference projections internally
        randomize_random_generator();
        fn_random = integerToString(ROUND(10000 * rnd_unif()));
        std::string command = (std::string)"-i " + fn_ref + " -o ref" + fn_random +
           "_ -sampling_rate " + integerToString(proj_step);
        if (fn_sym != "") command += (std::string)" -sym " + fn_sym;
        if (MPIversion)
            command=(std::string)"mpirun -np "+integerToString(numberOfProcessors)+
                " `which xmipp_mpi_create_projection_library` "+command;
        else
            command=(std::string)"xmipp_create_projection_library "+command;
        system(command.c_str());
        fn_ang = (std::string) "ref"+fn_random + "__angles.doc";
        fn_ref=(std::string) "ref"+fn_random + "_.sel";
    }
    else
    {
        if (fn_ang == "")
            REPORT_ERROR(1, "Prog_angular_predict_prm::produce_side_info:"
                         " using a reference selfile you must supply -ang option");
    }

    SF_ref.read(fn_ref);
    int refYdim, refXdim;
    SF_ref.ImgSize(refYdim, refXdim);
    if (refYdim != NEXT_POWER_OF_2(refYdim) || refXdim != NEXT_POWER_OF_2(refXdim))
        REPORT_ERROR(1, "Prog_angular_predict_prm::produce_side_info: "
                     "reference images must be of a size that is power of 2");

    // Produce side info of the angular distance computer
    distance_prm.fn_ang1 = distance_prm.fn_ang2 = "";
    distance_prm.fn_sym = fn_sym;
    distance_prm.produce_side_info();

    // Read the angle file
    DocFile DF;
    DF.read(fn_ang);
    DF.go_first_data_line();
    rot.resize(DF.dataLineNo());
    tilt.resize(DF.dataLineNo());
    int i = 0;
    while (!DF.eof())
    {
        rot[i] = DF(0);  // Rotational angle
        tilt[i] = DF(1); // Tilting angle
        i++;
        DF.next_data_line();
    }

    // Resize the predicted vectors
    int number_of_images = get_images_to_process();
    current_img = 0;
    image_name.resize(number_of_images);
    predicted_rot.resize(number_of_images);
    predicted_tilt.resize(number_of_images);
    predicted_psi.resize(number_of_images);
    predicted_shiftX.resize(number_of_images);
    predicted_shiftY.resize(number_of_images);
    predicted_corr.resize(number_of_images);
    predicted_reference.resize(number_of_images);

    // Build mask for subbands
    int Ydim, Xdim;
    SF_ref.ImgSize(Ydim, Xdim);
    Mask_no.resize(Ydim, Xdim);
    Mask_no.init_constant(-1);

    if (smax == -1) smax = Get_Max_Scale(Ydim) - 3;
    SBNo = (smax - smin + 1) * 3 + 1;
    SBsize.resize(SBNo);

    Mask_Params Mask(INT_MASK);
    Mask.type = BINARY_DWT_CIRCULAR_MASK;
    Mask.R1 = CEIL((double)Xdim / 2.0);
    Mask.resize(Ydim, Xdim);

    int m = 0, s;
    for (s = smax; s >= smin; s--)
    {
        for (int q = 0; q <= 3; q++)
        {
            if (q == 0 && s != smax) continue;
            Mask.smin = s;
            Mask.smax = s;
            Mask.quadrant = Quadrant2D(q);
            Mask.generate_2Dmask();

            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mask.imask2D)
            if (DIRECT_MAT_ELEM(Mask.imask2D, i, j))
            {
                Mask_no(i, j) = m;
                SBsize(m)++;
            }

            m++;
        }
    }

    // Produce library
    produce_library(rank);

    // Save a little space
    SF_ref.clear();

    // If dont_modify_header
    if (dont_modify_header)
        each_image_produces_an_output = false;
}

// Produce library -----------------------------------------------------------
void Prog_angular_predict_prm::produce_library(int rank)
{
    ImageXmipp I;
    int number_of_imgs = SF_ref.ImgNo();
    SF_ref.go_first_ACTIVE();
    set_DWT_type(DAUB12);

    // Create space for all the DWT coefficients of the library
    Matrix1D<int> SBidx(SBNo);
    for (int m = 0; m < SBNo; m++)
    {
        Matrix2D<double> *subband = new Matrix2D<double>;
        subband->resize(number_of_imgs, SBsize(m));
        library.push_back(subband);
    }
    library_power.initZeros(number_of_imgs, SBNo);

    if (rank==0) {
       std::cerr << "Generating reference library ...\n";
       init_progress_bar(number_of_imgs);
    }
    int n = 0, nstep = XMIPP_MAX(1, number_of_imgs / 60); // For progress bar
    while (!SF_ref.eof())
    {
        I.read(SF_ref.NextImg(), false, false, true, true);
        library_name.push_back(I.name());

        // Make and distribute its DWT coefficients in the different PCA bins
        I().statistics_adjust(0, 1);
        DWT(I(), I());
        SBidx.initZeros();
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Mask_no)
        {
            int m = Mask_no(i, j);
            if (m != -1)
            {
                double coef = I(i, j), coef2 = coef * coef;
                (*library[m])(n, SBidx(m)++) = coef;
                for (int mp = m; mp < SBNo; mp++) library_power(n, mp) += coef2;
            }
        }

        // Prepare for next iteration
        if (++n % nstep == 0 && rank==0) progress_bar(n);
    }
    if (rank==0) progress_bar(SF_ref.ImgNo());
}

// Build candidate list ------------------------------------------------------
void Prog_angular_predict_prm::build_ref_candidate_list(const ImageXmipp &I,
        std::vector<bool> &candidate_list, std::vector<double> &cumulative_corr,
        std::vector<double> &sumxy)
{
    int refNo = rot.size();
    candidate_list.resize(refNo);
    cumulative_corr.resize(refNo);
    sumxy.resize(refNo);
    for (int i = 0; i < refNo; i++)
    {
        candidate_list[i] = true;
        sumxy[i] = cumulative_corr[i] = 0;
        if (max_proj_change != -1)
        {
            double dummy_rot = rot[i], dummy_tilt = tilt[i], dummy_psi;
            double ang_distance = distance_prm.check_symmetries(
                                      I.rot(), I.tilt(), 0, dummy_rot, dummy_tilt, dummy_psi, true);
            candidate_list[i] = (ang_distance <= max_proj_change);
#ifdef DEBUG
            std::cout << "(" << I.rot() << "," << I.tilt() << ") and ("
            << rot[i] << "," << tilt[i] << ") --> " << ang_distance << std::endl;
#endif
        }
    }
}

// Refine candidate list ---------------------------------------------------
void Prog_angular_predict_prm::refine_candidate_list_with_correlation(
    int m,
    Matrix1D<double> &dwt,
    std::vector<bool> &candidate_list, std::vector<double> &cumulative_corr,
    Matrix1D<double> &x_power, std::vector<double> &sumxy,
    double th)
{
    histogram1D hist;
    hist.init(-1, 1, 201);

    int dimp = SBsize(m);
    int imax = candidate_list.size();
    Matrix2D<double> *library_m = library[m];
    for (int i = 0; i < imax; i++)
    {
        if (candidate_list[i])
        {
            double sumxyp = 0.0;
            for (int j = 0; j < dimp; j++)
                sumxyp += dwt(j) * (*library_m)(i, j);
            sumxy[i] += sumxyp;

            double corr = sumxy[i] / sqrt(library_power(i, m) * x_power(m));
            cumulative_corr[i] = corr;
            hist.insert_value(corr);

            if (tell & TELL_ROT_TILT)
            {
                std::cout << "Candidate " << i << " corr= " << cumulative_corr[i]
                << " rot= " << rot[i] << " tilt= " << tilt[i] << std::endl;
            }
        }
    }
    double corr_th = hist.percentil(100 - th);

    // Remove all those projections below the threshold
    for (int i = 0; i < imax; i++)
        if (candidate_list[i] && cumulative_corr[i] < corr_th)
            candidate_list[i] = false;

    // Show the percentil used
    if (tell & TELL_ROT_TILT)
        std::cout << "# Percentil " << corr_th << std::endl << std::endl;
}

// Predict rot and tilt ----------------------------------------------------
double Prog_angular_predict_prm::predict_rot_tilt_angles(ImageXmipp &I,
        double &assigned_rot, double &assigned_tilt, int &best_ref_idx)
{
    if (XSIZE(I()) != NEXT_POWER_OF_2(XSIZE(I())) ||
        YSIZE(I()) != NEXT_POWER_OF_2(YSIZE(I())))
        REPORT_ERROR(1, "Prog_angular_predict_prm::predict_rot_tilt_angles: "
                     "experimental images must be of a size that is power of 2");

    // Build initial candidate list
    std::vector<bool>   candidate_list;
    std::vector<double> cumulative_corr;
    std::vector<double> sumxy;
    build_ref_candidate_list(I, candidate_list, cumulative_corr, sumxy);
    int imax = candidate_list.size();

    // Make DWT of the input image and build vectors for comparison
    std::vector<Matrix1D<double> * > Idwt;
    Matrix1D<double> x_power(SBNo); x_power.initZeros();
    Matrix1D<int> SBidx(SBNo); SBidx.initZeros();
    for (int m = 0; m < SBNo; m++)
    {
        Matrix1D<double> *subband = new Matrix1D<double>;
        subband->resize(SBsize(m));
        Idwt.push_back(subband);
    }

    I().statistics_adjust(0, 1);
    DWT(I(), I());
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mask_no)
    {
        int m = Mask_no(i, j);
        if (m != -1)
        {
            double coef = DIRECT_IMGPIXEL(I, i, j), coef2 = coef * coef;
            (*(Idwt[m]))(SBidx(m)++) = coef;
            for (int mp = m; mp < SBNo; mp++) x_power(mp) += coef2;
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
    for (int i = 0; i < imax; i++)
        if (candidate_list[i])
            if (first)
            {
                best_i = i;
                first = false;
                N_max = 1;
            }
            else if (cumulative_corr[i] > cumulative_corr[best_i]) best_i = i;
            else if (cumulative_corr[i] == cumulative_corr[best_i])
                N_max++;

    if (N_max == 0)
    {
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
    for (int m = 0; m < SBNo; m++) delete Idwt[m];

    assigned_rot    = rot[best_i];
    assigned_tilt   = tilt[best_i];
    best_ref_idx    = best_i;
    return cumulative_corr[best_i];
}

// Evaluate candidates by correlation ----------------------------------------
double Prog_angular_predict_prm::evaluate_candidates(
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
        if (vscore[i] < min_score) min_score = vscore [i];
        if (vscore[i] > max_score) max_score = vscore [i];
    }

    // Divide the correlation segment in as many pieces as candidates
    double score_step = (max_score - min_score) / 10;

    int jmax = candidate_idx.size();
    for (int j = 0; j < jmax; j++)
    {
        int i = candidate_idx[j];
        int points;
        if (score_step != 0) points = FLOOR((vscore[i] - min_score) / score_step);
        else               points = 10;
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
void Prog_angular_predict_prm::group_views(const std::vector<double> &vrot,
        const std::vector<double> &vtilt, const std::vector<double> &vpsi,
        const std::vector<int> &best_idx, const std::vector<int> &candidate_idx,
        std::vector< std::vector<int> > &groups)
{
    for (int j = 0; j < best_idx.size(); j++)
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
        int g;
        for (g = 0; g < groups.size(); g++)
        {
            bool fits = true;
            for (int jp = 0; jp < groups[g].size(); jp++)
            {
                int ip = candidate_idx[groups[g][jp]];
                double ang_distance = distance_prm.check_symmetries(
                                          vrot[ip], vtilt[ip], vpsi[ip], roti, tilti, psii, false);
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
        for (int j = 0; j < best_idx.size(); j++) group.push_back(best_idx[j]);
        groups.push_back(group);
    }
}
#undef DEBUG

// Pick view -----------------------------------------------------------------
int Prog_angular_predict_prm::pick_view(int method,
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
        for (int g = 0; g < groups.size(); g++)
        {
            double temp_group_rate = 0;
            for (int j = 0; j < groups[g].size(); j++)
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
        for (int g = 0; g < groups.size(); g++)
            if (group_rate[g] == best_group_rate) groups_with_max_rate++;
#ifdef DEBUG
        if (groups_with_max_rate > 1)
        {
            std::cout << "There are two groups with maximum rate\n";
        }
#endif

        // Take the best image within that group
        int best_j;
        double best_rate = -1e38;
        for (int j = 0; j < groups[best_g].size(); j++)
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
        for (int j = 0; j < groups[best_g].size(); j++)
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
            for (int j = 0; j < groups[best_g].size(); j++)
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
}

// Predict shift and psi -----------------------------------------------------
//#define DEBUG
double Prog_angular_predict_prm::predict_angles(ImageXmipp &I,
        double &assigned_shiftX, double &assigned_shiftY,
        double &assigned_rot, double &assigned_tilt, double &assigned_psi)
{
    double best_rot, best_tilt, best_psi, best_shiftX, best_shiftY,
    best_score = 0, best_rate;

    ImageXmipp Ip;
    Ip = I;
    Matrix1D<double> shift(2);

    // Get the 2D alignment shift
    double Xoff, Yoff;
    if (!dont_apply_geo)
    {
        Xoff = I.Xoff();
        Yoff = I.Yoff();
    }
    else
        Xoff = Yoff = 0.0;

    // Establish psi limits
    double psi0, psiF;
    if (max_psi_change == -1)
    {
        psi0 = -180;
        psiF = 180 - psi_step;
    }
    else
    {
        psi0 = I.psi() - max_psi_change;
        psiF = I.psi() + max_psi_change;
    }
    double R2 = max_shift_change * max_shift_change;

    // Search in the psi-shift space
    int N_trials = 0;
    std::vector<double> vshiftX, vshiftY, vpsi, vrot, vtilt, vcorr,
    vproj_error, vproj_compact, vang_jump, vscore;
    std::vector<int>    vref_idx;

    double backup_max_shift_change = max_shift_change;
    if (!search5D) max_shift_change = 0;

    for (int mirror = 0; mirror <= check_mirrors; mirror++)
        for (double shiftX = Xoff - max_shift_change; shiftX <= Xoff + max_shift_change; shiftX += shift_step)
            for (double shiftY = Yoff - max_shift_change; shiftY <= Yoff + max_shift_change; shiftY += shift_step)
            {
                if ((shiftX - Xoff)*(shiftX - Xoff) + (shiftY - Yoff)*(shiftY - Yoff) > R2) continue;
                for (double psi = psi0; psi <= psiF; psi += psi_step)
                {
                    N_trials++;

                    // Shift image if necessary
                    if (shiftX == 0 && shiftY == 0) Ip() = I();
                    else
                    {
                        VECTOR_R2(shift, shiftX, shiftY);
                        I().translate(shift, Ip(), WRAP);
                    }

                    // Rotate image if necessary
                    // Adding 2 is a trick to avoid that the 0, 90, 180 and 270
                    // are treated in a different way
                    Ip().selfRotate(psi + 2, WRAP);
                    Ip().selfRotate(-2, WRAP);

                    // Mirror the image if necessary
                    if (mirror) Ip().selfReverseY();

#ifdef DEBUG
                    Ip.write("PPPbefore_denoising.xmp");
#endif

                    // Project the resulting image onto the visible space
                    double proj_error = 0.0, proj_compact = 0.0;
                    bool   look_for_rot_tilt = true;

                    // Search for the best tilt, rot angles
                    double rotp, tiltp;
                    int best_ref_idx;
                    double corrp =
                        predict_rot_tilt_angles(Ip, rotp, tiltp, best_ref_idx);

                    double aux_rot = rotp, aux_tilt = tiltp, aux_psi = psi;
                    double ang_jump = distance_prm.check_symmetries(
                                          I.rot(), I.tilt(), I.psi(),
                                          aux_rot, aux_tilt, aux_psi, false);

                    // Check if the image was mirrored
                    double rotp_aux = rotp, tiltp_aux = tiltp, psi_aux = psi;
                    if (mirror) Euler_up_down(rotp, tiltp, psi,
                                                  rotp_aux, tiltp_aux, psi_aux);

                    vshiftX.push_back(shiftX);
                    vshiftY.push_back(shiftY);
                    vrot.push_back(rotp_aux);
                    vtilt.push_back(tiltp_aux);
                    vpsi.push_back(psi_aux);
                    vcorr.push_back(corrp);
                    vproj_error.push_back(proj_error);
                    vproj_compact.push_back(proj_compact);
                    vang_jump.push_back(ang_jump);
                    vref_idx.push_back(best_ref_idx);

#ifdef DEBUG
                    Ip.write("PPPafter_denoising.xmp");
                    ImageXmipp Iref((std::string)"g0tb" + integerToString(best_ref_idx + 1, 5) + ".xmp");
                    Iref.write("PPPref.xmp");
                    std::cerr << "Press any key\n";
                    char c;
                    std::cin >> c;
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
        if (vproj_error[i] < min_proj_error) min_proj_error = vproj_error[i];
        if (vproj_error[i] > max_proj_error) max_proj_error = vproj_error[i];

        // Compute extrema of correlation
        if (vcorr[i] < min_corr) min_corr = vcorr[i];
        if (vcorr[i] > max_corr) max_corr = vcorr[i];

        // Compute extrema of projection error
        if (vproj_compact[i] < min_proj_compact) min_proj_compact = vproj_compact[i];
        if (vproj_compact[i] > max_proj_compact) max_proj_compact = vproj_compact[i];
    }

    // Score each projection
    double corr_step = max_corr - min_corr;
    double proj_error_step = max_proj_error - min_proj_error;
    double proj_compact_step = max_proj_compact - min_proj_compact;
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
    if (tell & TELL_PSI_SHIFT) std::cout << "Local maxima\n";
    for (int i = 0; i < N_trials; i++)
    {
        // Look for the left and right sample
        int il = i - 1, ir = i + 1;
        if (i == 0 && circular) il = N_trials - 1;
        else if (i == N_trials - 1)
            if (circular) ir = 0;
            else ir = -1;

        // Check if the error is a local minimum of the projection error
        // or a local maxima of the correlation
        bool is_local_maxima = true;
        if (il != -1)
            if (vscore[il] >= vscore[i]) is_local_maxima = false;
        if (ir != -1)
            if (vscore[ir] >= vscore[i]) is_local_maxima = false;

        // It is a local minimum
        if (is_local_maxima /*|| visible_space*/) avg_score_maxima += vscore[i];
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
    if (tell & TELL_PSI_SHIFT) std::cout << "Keeping ...\n";
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
    if (tell & TELL_PSI_SHIFT) std::cout << "\nSelecting image\n";
    Matrix1D<double> score(jmax);
//   for (int j=0; j<jmax; j++) score(j)=candidate_rate[j];
    for (int j = 0; j < jmax; j++) score(j) = vscore[candidate_local_maxima[j]];
    Matrix1D<int> idx_score = score.indexSort();

    if (tell & (TELL_PSI_SHIFT | TELL_OPTIONS))
    {
        std::cout << I.name() << std::endl;
        std::cout.flush();
        for (int j = 0; j < jmax; j++)
        {
            int jp = idx_score(j) - 1;
            double score = candidate_rate[jp];
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
        ImageXmipp Iref;
        Iref.read(library_name[vref_idx[ibest]]);
        Iref().setXmippOrigin();
        if (Xoff == 0 && Yoff == 0) Ip() = I();
        else
        {
            VECTOR_R2(shift, Xoff, Yoff);
            I().translate(shift, Ip(), WRAP);
        }

        double shiftX, shiftY;
        best_shift(Iref(), Ip(), shiftX, shiftY);
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
        std::cout << "Originally it had, psi=" << I.psi() << " rot=" << I.rot()
        << " tilt=" << I.tilt() << std::endl;
        std::cout << "Finally I choose: ";
        if (tell & TELL_PSI_SHIFT) std::cout << jbest << "\n";
        std::cout << "psi= " << best_psi << " rot= " << best_rot << " tilt= "
        << best_tilt << " shiftX=" << best_shiftX
        << " shiftY=" << best_shiftY << " score= " << best_score
        << " rate= " << best_rate << std::endl << std::endl;
    }

    // Save results
    image_name[current_img]          = I.name();
    assigned_rot    = predicted_rot[current_img]       = best_rot;
    assigned_tilt   = predicted_tilt[current_img]      = best_tilt;
    assigned_psi    = predicted_psi[current_img]       = best_psi;
    assigned_shiftX = predicted_shiftX[current_img]    = best_shiftX;
    assigned_shiftY = predicted_shiftY[current_img]    = best_shiftY;
    predicted_corr[current_img]      = best_score;
    predicted_reference[current_img] = vref_idx[ibest];
    current_img++;
    return best_rate;
}
#undef DEBUG

// Finish processing ---------------------------------------------------------
void Prog_angular_predict_prm::finish_processing()
{
    // Save predicted angles
    int p = predicted_rot.size();
    DocFile DF;
    DF.reserve(p + 1);
    DF.append_comment("Predicted_Rot Predicted_Tilt Predicted_Psi Predicted_ShiftX Predicted_ShiftY Corr");
    Matrix1D<double> v(6);
    for (int i = 0; i < p; i++)
    {
        v(0) = predicted_rot[i];
        v(1) = predicted_tilt[i];
        v(2) = predicted_psi[i];
        v(3) = predicted_shiftX[i];
        v(4) = predicted_shiftY[i];
        v(5) = predicted_corr[i];
        DF.append_comment(image_name[i]);
        DF.append_data_line(v);
    }
    DF.write(fn_out_ang);

    if (summaryRootname != "") produceSummary();

    if (volume_mode)
    {
        FileName fn_root=fn_ref.without_extension();
        system(((std::string)"xmipp_rmsel " + fn_ref + " > /dev/null").c_str());
        system(((std::string)"rm -f " + fn_ang).c_str());
        system(((std::string)"rm -f " + fn_root+"_vectors.doc").c_str());
    }
}

// Produce summary ---------------------------------------------------------
void Prog_angular_predict_prm::produceSummary()
{
    int L = library_name.size();
    int N = predicted_rot.size();
    int Xdim, Ydim;
    SelFile SFin;
    SFin.read(fn_in);
    SFin.ImgSize(Ydim, Xdim);
    SFin.go_first_ACTIVE();

    // Initialize variables for storing the reference weights
    // and the assigned averages
    std::vector< std::vector<FileName> > referenceWeight;
    std::vector< ImageXmipp >  assignedAvg;
    referenceWeight.reserve(L);
    assignedAvg.reserve(L);
    ImageXmipp blankImage(Ydim, Xdim);
    std::vector<FileName> blankList;
    for (int l = 0; l < L; l++)
    {
        referenceWeight.push_back(blankList);
        assignedAvg.push_back(blankImage);
    }

    // Count the images assigned to each reference
    // And compute average
    for (int n = 0; n < N; n++)
    {
        ImageXmipp I;
        I.read(SFin.NextImg(), false, false, true);
        referenceWeight[predicted_reference[n]].push_back(I.name());
        I().selfRotate(I.psi());
        assignedAvg[predicted_reference[n]]() += I();
    }
    for (int l = 0; l < L; l++)
        if (referenceWeight[l].size() != 0)
            assignedAvg[l]() /= referenceWeight[l].size();

    // Construct summary
    DocFile DFsummary;
    SelFile SFsummary;
    SelFile SFcomparison;
    DFsummary.reserve(L + 1);
    DFsummary.append_comment("Rot Tilt Psi X Y Weight");
    Matrix1D<double> v(6);
    for (int l = 0; l < L; l++)
    {
        // Write the reference and assigned reference
        FileName fn_refl = summaryRootname + "_ref" + integerToString(l, 5) + ".xmp";
        FileName fn_avgl = summaryRootname + "_avg" + integerToString(l, 5) + ".xmp";
        ImageXmipp I;
        I.read(library_name[l]);
        I.write(fn_refl);
        assignedAvg[l].write(fn_avgl);
        SFcomparison.insert(fn_refl.remove_directories());
        SFcomparison.insert(fn_avgl.remove_directories());

        // Write corresponding line in the docfile
        v(0) = rot[l];
        v(1) = tilt[l];
        v(2) = 0;
        v(3) = 0;
        v(4) = 0;
        v(5) = referenceWeight[l].size();
        DFsummary.append_comment(fn_refl);
        DFsummary.append_data_line(v);

        // Write the assigned images in a selfile
        SFsummary.insert_comment((std::string)"Images assigned to " + fn_refl);
        SFsummary.insert_comment((std::string)"rot=" + floatToString(rot[l]) + " tilt=" + floatToString(tilt[l]));
        if (referenceWeight[l].size() != 0)
        {
            for (int i = 0; i < referenceWeight[l].size(); i++)
                SFsummary.insert(referenceWeight[l][i]);
            SFsummary.insert_comment("");
        }
    }
    DFsummary.write(summaryRootname + "_summary.doc");
    SFsummary.write(summaryRootname + "_summary.sel");
    SFcomparison.write(summaryRootname + "_comparison.sel");
}
