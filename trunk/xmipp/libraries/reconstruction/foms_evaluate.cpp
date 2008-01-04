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

#include <data/args.h>
#include <data/histogram.h>
#include <data/mask.h>
#include <data/selfile.h>

#include "foms_evaluate.h"
#include "volume_foms.h"

#include <fstream>

/* Compute special statistics ============================================== */
/* Average of a vector without counting the -1, starting at index i0 */
void avg_without__1(Matrix1D<double> &v, double &avg, int i0)
{
    int N = 0;
    avg = 0;
    for (int i = i0; i < XSIZE(v); i++)
        if (v(i) != -1)
        {
            N++;
            avg += v(i);
        }
    if (N != 0) avg /= N;
    else avg = -1;
}

/* Average and standard deviation  without counting the -1 */
void stats__1(const Matrix1D<double> &v, double &avg, double &stddev)
{
    int N = 0;
    avg = stddev = 0;
    for (int i = 0; i < XSIZE(v); i++)
        if (v(i) != -1)
        {
            N++;
            avg += v(i);
            stddev += v(i) * v(i);
        }
    if (N != 0)
    {
        avg /= N;
        stddev = sqrt(XMIPP_MAX(0, stddev / N - avg * avg));
    }
    else
    {
        avg = -1;
        stddev = 0;
    }
}

/* Default values ========================================================== */
void Prog_Evaluate_Parameters::default_values()
{
    back_radius    = 0;
    back_factor    = 1.5;
    back_mode      = ENLARGE_MODE;
    tell           = 0;
    fn_sel         = "";
    fnPhantom     = "";
    fn_recons      = "";
    percent_mass   = 99;
    global_radius  = 0;
    fn_mask        = "";
    fit_gray_scales = true;
    RSrot          = 0;
    RStilt         = 0;
}

/* Read Evaluate parameters from command line ============================== */
void Prog_Evaluate_Parameters::read(int argc, char **argv)
{
    int i;

    // By default the enlarge factor is chosen
    default_values();

    // Read from command line
    fnPhantom       = getParameter(argc, argv, "-p", "");
    if (checkParameter(argc, argv, "-sel"))
        fn_sel = getParameter(argc, argv, "-sel");
    else
    {
        fn_recons        =      getParameter(argc, argv,  "-r");
    }
    percent_mass     = textToFloat(getParameter(argc, argv,  "-mass"       , "99"));
    global_radius    = textToFloat(getParameter(argc, argv,  "-R"          , "0"));
    fn_mask          = getParameter(argc, argv, "-mask", "");
    fit_gray_scales  = checkParameter(argc, argv, "-fit_gray");
    if (checkParameter(argc, argv, "-back_radius"))
    {
        back_radius   = textToFloat(getParameter(argc, argv,  "-back_radius"));
        back_mode = SPHERE_MODE;
    }
    if (checkParameter(argc, argv, "-back_factor"))
    {
        back_factor   = textToFloat(getParameter(argc, argv,  "-back_factor"));
        back_mode = ENLARGE_MODE;
    }
    if ((i = paremeterPosition(argc, argv, "-dir")) != -1)
    {
        if ((++i) < argc)
        {
            if (strcmp(argv[i], "X") == 0)
            {
                RSrot = 0;
                RStilt = 90;
            }
            else if (strcmp(argv[i], "Y") == 0)
            {
                RSrot = 90;
                RStilt = 90;
            }
            else if (strcmp(argv[i], "Z") == 0)
            {
                RSrot = 0;
                RStilt = 0;
            }
            else
            {
                RSrot = textToFloat(argv[i]);
                if ((++i) < argc) RStilt = textToFloat(argv[i]);
            }
        }
    }
    if (checkParameter(argc, argv, "-save_maps"))       tell |= SAVE_MAPS;
    if (checkParameter(argc, argv, "-show_values"))     tell |= SHOW_VALUES;
    if (checkParameter(argc, argv, "-show_process"))    tell |= SHOW_PROCESS;
    if (checkParameter(argc, argv, "-save_histograms")) tell |= SAVE_HISTOGRAMS;
    if (checkParameter(argc, argv, "-only_structural")) tell |= ONLY_STRUCTURAL;
}

/* Evaluate usage ========================================================== */
void Prog_Evaluate_Parameters::usage()
{
    printf("Error in the arguments\n");
    printf("Usage: \n");
    printf("       evaluate <options>\n");
    printf("-p <phantom>             : it can be either a description or a Spider volume\n");
    printf("-r <reconstruction> |    : a Spider volume\n");
    printf("-sel <selfile>           : with all reconstruction names\n");
    printf("[-mass <%%>]             : leave out this percentage of mass in the histograms\n");
    printf("                           [99]\n");
    printf("[-R <r>]                 : The global error will be measured within this radius\n");
    printf("[-mask <surface mask>]   : surface mask applied during reconstruction\n");
    printf("[-mask]                  : only for sel files where the mask name is\n"
           "                           automatically computed\n");
    printf("[-fit_gray]              : Fit gray scales before evaluating\n");
    printf("[-dir <rot> <tilt>]      : to perform the directional FOM and the slice histograms\n");
    printf("                           a direction must be specified by two Euler angles\n");
    printf("                           by default the Z axis is taken [0 0]\n"
           "                           other useful axis are [0 90] --> X\n"
           "                           and [90 90] --> Y\n"
           "[-dir X|Y|Z]\n           : To perform directional FOMs from any of these axis\n"
           "[-back_radius <radius> |]: if this option is not given the background\n"
           " -back_factor <factor> ]   of the features is supposed to be the same\n"
           "                           feature enlarged by 1.25. With this option\n"
           "                           the background will be a sphere of the given\n"
           "                           radius, or the enlarging factor may be changed\n");
    printf("[-save_maps]             : Save different volume maps\n"
           "[-show_values]           : Show values in the features\n"
           "[-show_process]          : Show more information during calculations\n"
           "[-save_histograms]       : Save involved histograms\n"
           "[-only_structural]       : Only compute the structural consistency FOMs\n");
}

/* Show parameters ========================================================= */
std::ostream & operator << (std::ostream &out, const Prog_Evaluate_Parameters &prm)
{
    out << "Evaluating parameters ----------------------\n";
    out << "Phantom         : " << prm.fnPhantom        << std::endl;
    out << "Reconstruction  : " << prm.fn_recons         << std::endl;
    out << "Percent mass    : " << prm.percent_mass      << std::endl;
    out << "RSrot           : " << prm.RSrot             << std::endl;
    out << "RStilt          : " << prm.RStilt            << std::endl;
    out << "Global radius   : " << prm.global_radius     << std::endl;
    out << "Surface mask    : " << prm.fn_mask           << std::endl;
    out << "Fit gray scales : " << prm.fit_gray_scales   << std::endl;
    out << "Back mode       : " << prm.back_mode         << std::endl;
    out << "Back radius     : " << prm.back_radius       << std::endl;
    out << "Back factor     : " << prm.back_factor       << std::endl;
    out << "Tell            : " << prm.tell              << std::endl;
    return out;
}

/* Produce Side information ================================================ */
void EVALUATE_Side_Info::produce_Side_Info(
    const Prog_Evaluate_Parameters &prm)
{

    // Set background mode ..................................................
    if (prm.back_mode == SPHERE_MODE) back_param = prm.back_radius;
    else                            back_param = prm.back_factor;

    // Read reconstruction ..................................................
    fn_root = prm.fn_recons.without_extension();
    vol_recons.read(prm.fn_recons);
    vol_recons.moveOriginTo_center();

    // Read phantom and label ...............................................
    if (Is_VolumeXmipp(prm.fnPhantom))
    {
        descr_mode = XMIPP_PHANTOM;
        vol_phantom.read(prm.fnPhantom);
        vol_phantom.moveOriginTo_center();
        vol_label().resize(vol_phantom());
        vol_label().init_constant(1);
        num_feat = 0;
    }
    else
    {
        std::cerr << "Generating phantom ...\n";
        phantom_descr.read(prm.fnPhantom);
        phantom_descr.draw_in(&vol_phantom);
        phantom_descr.label(&vol_label);
        num_feat = phantom_descr.FeatNo();
        descr_mode = MATH_PHANTOM;
    }

    // Check that both dimensions are equal .................................
    if ((vol_phantom().sliceNumber() != vol_recons().sliceNumber()) ||
        (vol_phantom().rowNumber() != vol_recons().rowNumber()) ||
        (vol_phantom().colNumber() != vol_recons().colNumber()))
    {
        std::cout << "Be careful!!!, volumes with different sizes\n";
        std::cout << "Phantom:        " << vol_phantom().sliceNumber() << " x " <<
        vol_phantom().rowNumber() << " x " << vol_phantom().colNumber() << std::endl;
        std::cout << "Reconstruction: " << vol_recons().sliceNumber() << " x " <<
        vol_recons().rowNumber() << " x " << vol_recons().colNumber() << std::endl;

        cutToCommonSize(vol_phantom(), vol_recons());
        std::cout << "Cutting to common size " << vol_phantom().sliceNumber() << " x " <<
        vol_phantom().rowNumber() << " x " << vol_phantom().colNumber() << std::endl;

        cutToCommonSize(vol_label(), vol_recons());
    }

    // Generate global mask .................................................
    std::cerr << "Generating mask ...\n";
    if (prm.fn_mask != "")
    {
        vol_mask.read(prm.fn_mask);
        invert_binary_mask(vol_mask());
        vol_mask().setXmippOrigin();

        // Find minimum and maximum plane used for the surface
        int ztop = STARTINGZ(vol_mask());
        int zbottom = FINISHINGZ(vol_mask());
        for (int i = STARTINGY(vol_mask()); i <= FINISHINGY(vol_mask()); i++)
            for (int j = STARTINGX(vol_mask()); j <= FINISHINGX(vol_mask()); j++)
            {
                int state;
                state = 0;
                for (int k = STARTINGZ(vol_mask()); k <= FINISHINGZ(vol_mask()); k++)
                {
                    if (state == 0 && VOLVOXEL(vol_mask, k, i, j) == 1)
                    {
                        ztop = XMIPP_MAX(ztop, k);
                        state = 1;
                    }
                    else if (state == 1 && VOLVOXEL(vol_mask, k, i, j) == 0)
                    {
                        zbottom = XMIPP_MIN(zbottom, k);
                        state = 2;
                    }
                }
            }

        // Now replace all mask values outside these two by 0
        FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_mask())
        if (k < ztop || k > zbottom) VOLVOXEL(vol_mask, k, i, j) = 0;
    }
    else
    {
        vol_mask().resize(vol_recons());
        vol_mask().init_constant(1);
    }

    if (prm.global_radius != 0)
    {
        float global_radius2 = prm.global_radius * prm.global_radius;
        FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(vol_mask))
        {
            float r2 = k * k + i * i + j * j;
            if (r2 > global_radius2) VOLVOXEL(vol_mask, k, i, j) = 0;
        }
    }

    // Fit gray values ......................................................
    if (prm.fit_gray_scales)
    {
        std::cerr << "Fitting gray values ...\n";
        range_adjust_within_mask(&(vol_mask()), vol_phantom(), vol_recons());
    }

    // Computing distance map ...............................................
    if (!((prm.tell&ONLY_STRUCTURAL) || descr_mode == XMIPP_PHANTOM))
    {
        std::cerr << "Computing distance map ...\n";
        compute_distance_map(&vol_label, phantom_descr, &vol_mask,
                             &vol_distance);
    }
}

/* Compute FOMs ============================================================ */
void compute_FOMs(const Prog_Evaluate_Parameters &prm,
                  EVALUATE_Side_Info &side, EVALUATE_results &results)
{
    Matrix1D<double> feat_voxels;
    histogram1D hist_recons;

    /* Structural consistency FOMs --------------------------------------------- */
    // Global measures
    std::cerr << "Computing global structural consistency ...\n";
    compute_sc_FOMs(&side.vol_phantom, &side.vol_recons, &side.vol_label,
                    &side.vol_mask, -1, results.scL2_FOM, results.scL1_FOM,
                    results.scmu_FOM, results.scdev_FOM, results.scrange_FOM,
                    results.sccorr_FOM, results.scinf_FOM,
                    prm.tell&SHOW_PROCESS);
#define COMPUTE_THROUGH_SINGLE_VALUE
#ifdef COMPUTE_THROUGH_SINGLE_VALUE
    compute_resolution(side.vol_phantom, side.vol_recons, results.resol_FOM);
#else
    Matrix1D<double> frequency, FSC;
    results.resol_FOM = compute_FSC(side.vol_phantom, side.vol_recons,
                                    1, frequency, FSC);
#endif

    // Local measures
    results.scL2_FOMs.initZeros(side.num_feat + 1); // 0, 1, ..., FeatNo()
    results.scL2_FOMs.init_constant(-1);
    results.scL1_FOMs = results.scmu_FOMs = results.scdev_FOMs =
                                                results.scrange_FOMs = results.sccorr_FOMs = results.scinf_FOMs =
                                                                           results.scL2_FOMs;
    if (side.descr_mode == MATH_PHANTOM)
    {
        std::cerr << "Computing Local structural consistency ...\n";
        for (int i = 0; i <= side.num_feat; i++)
            compute_sc_FOMs(&side.vol_phantom, &side.vol_recons, &side.vol_label,
                            &side.vol_mask, i, results.scL2_FOMs(i), results.scL1_FOMs(i),
                            results.scmu_FOMs(i), results.scdev_FOMs(i),
                            results.scrange_FOMs(i), results.sccorr_FOMs(i),
                            results.scinf_FOMs(i), prm.tell&SHOW_PROCESS);
    }

    // Weighted L2 and L1 measure
    if (side.descr_mode == MATH_PHANTOM)
    {
        compute_voxels_in_feat(&side.vol_label, feat_voxels);
        double voxels_inside_feat = 0;
        for (int i = 0; i < XSIZE(feat_voxels); i++)
            voxels_inside_feat += feat_voxels(i);
        double avg_scL2f, avg_scL1f;
        avg_without__1(results.scL2_FOMs, avg_scL2f, 1);
        avg_without__1(results.scL1_FOMs, avg_scL1f, 1);
        results.scL2w_FOM = 0.5 * (feat_voxels.sum() * (results.scL2_FOMs(0) / feat_voxels(0) +
                                   avg_scL2f / voxels_inside_feat));
        results.scL1w_FOM = 0.5 * (feat_voxels.sum() * (results.scL1_FOMs(0) / feat_voxels(0) +
                                   avg_scL1f / voxels_inside_feat));
    }
    else
    {
        results.scL2w_FOM = results.scL1w_FOM = -1;
    }

    results.hsmu_FOMs.resize(side.num_feat + 1); // 0, 1, ..., FeatNo()
    results.hsmu_FOMs.init_constant(-1);
    results.hsvr_FOMs = results.hsdt_FOMs = results.hsbr_FOMs = results.hsmu_FOMs;
    results.hsmu_FOM = results.hsdt_FOM = results.hsbr_FOM = results.hsvr_FOM = -1;
    results.drrt_FOM = results.dsbl_FOM = results.dsad_FOM = -1;

    if (!(prm.tell&ONLY_STRUCTURAL))
    {
        /* Histogram based FOMs ------------------------------------------------- */
        // Local measures
        if (side.descr_mode == MATH_PHANTOM)
        {
            // FOM for each feature
            std::cerr << "Computing Histogram based FOMs ...\n";
            for (int i = 1; i <= side.num_feat; i++)
                compute_hs_FOMs(&side.vol_phantom, &side.vol_recons, &side.vol_label,
                                &side.vol_mask, i, side.phantom_descr, prm.back_mode,
                                side.back_param, results.hsmu_FOMs(i), results.hsbr_FOMs(i),
                                results.hsdt_FOMs(i), results.hsvr_FOMs(i),
                                prm.tell&(SHOW_PROCESS | SAVE_HISTOGRAMS),
                                side.fn_root + "_eval_histog.plot");

            // Global FOM
            avg_without__1(results.hsmu_FOMs, results.hsmu_FOM, 1);
            avg_without__1(results.hsbr_FOMs, results.hsbr_FOM, 1);
            avg_without__1(results.hsdt_FOMs, results.hsdt_FOM, 1);
            avg_without__1(results.hsvr_FOMs, results.hsvr_FOM, 1);
        }

        /* Directional FOMs ----------------------------------------------------- */
        // Global measures
        std::cerr << "Computing directional FOMs ...\n";
        compute_hist(side.vol_recons(), hist_recons, 200);
        double threshold = hist_recons.percentil(prm.percent_mass);

        compute_dr_FOMs(&side.vol_phantom, &side.vol_recons, &side.vol_mask,
                        prm.RSrot, prm.RStilt,
                        &results.img_histog, 100, threshold, results.drrt_FOM,
                        prm.tell&SHOW_PROCESS, side.fn_root + "_eval_radon.plot");
        if (prm.tell & SAVE_HISTOGRAMS)
            results.img_histog.write(side.fn_root + "_eval_slice_histog.xmp");

        /* Distance map based --------------------------------------------------- */
        if (side.descr_mode == MATH_PHANTOM)
        {
            std::cerr << "Computing distance based FOMs ...\n";
            compute_ds_FOMs(&side.vol_phantom, &side.vol_recons, &side.vol_label,
                            &side.vol_distance, results.dsbl_FOM, results.dsad_FOM);

            if (prm.tell&SHOW_PROCESS)
                show_shape(&side.vol_phantom, &side.vol_recons, &side.vol_label,
                           &side.vol_distance, side.fn_root + "_eval_shape.plot");

        }
    }
}

/* Show FOMs =============================================================== */
void show_FOMs(const Prog_Evaluate_Parameters &prm,
               EVALUATE_Side_Info &side, const EVALUATE_results &results)
{

    // Show Parameters ......................................................
    std::cout << std::endl;
    std::cout << "PHANTOM FILE      : " << prm.fnPhantom << std::endl;
    std::cout << "RECONSTRUCTED FILE: " << prm.fn_recons  << std::endl;
    if (!(prm.tell&ONLY_STRUCTURAL))
    {
        std::cout << "Direction for dFOM: (rot=" << prm.RSrot << "," << "tilt="
        << prm.RStilt << ")" << std::endl;
        std::cout << "Slice histograms  : ignoring initial " << prm.percent_mass <<
        "% mass\n";
    }
    if (prm.global_radius == 0)
        std::cout << "Global FOMs measured over the whole volume\n";
    else
        std::cout << "Global FOMs measured over a sphere of radius "
        << prm.global_radius << std::endl;
    if (prm.back_mode == ENLARGE_MODE)
        std::cout << "Background mode: ENLARGE by " << prm.back_factor << std::endl;
    else
        std::cout << "Background mode: SPHERES of radius " << prm.back_radius << std::endl;

    if (side.descr_mode == MATH_PHANTOM && (prm.tell & SHOW_PROCESS))
    {
        std::cout << "Phantom description ----------------------------------------\n";
        std::cout << side.phantom_descr;
    }

    std::cout << "Structural consistency -------------------------------------\n";

    // Show volume statistics ...............................................
    double avg, stddev, min, max;
    side.vol_phantom().computeStats(avg, stddev, min, max);
    std::cout << "Phantom Stats: \n";
    std::cout << "   ";
    side.vol_phantom().print_stats();
    std::cout << std::endl;
    std::cout << "    range=" << max - min << std::endl;
    side.vol_recons().computeStats(avg, stddev, min, max);
    std::cout << "Recons  Stats: \n";
    std::cout << "   ";
    side.vol_recons().print_stats();
    std::cout << std::endl;
    std::cout << "    range=" << max - min << std::endl;
    std::cout << std::endl;

    // Show Structural consistency ..........................................
    printf("scL2        FOM:%f\n", results.scL2_FOM);
    printf("scL1        FOM:%f\n", results.scL1_FOM);
    printf("scL2w       FOM:%f\n", results.scL2w_FOM);
    printf("scL1w       FOM:%f\n", results.scL1w_FOM);
    printf("scmu        FOM:%f\n", results.scmu_FOM);
    printf("scdev       FOM:%f\n", results.scdev_FOM);
    printf("scrange     FOM:%f\n", results.scrange_FOM);
    printf("sccorr      FOM:%f\n", results.sccorr_FOM);
    printf("scinf       FOM:%f\n", results.scinf_FOM);
    printf("resolution  FOM:%f\n", results.resol_FOM);

    printf("\tFEATURE   scL2     scL1     scmu     scdev   scrange  sccorr    scinf \n");
    printf("\t------- -------- -------- -------- -------- -------- -------- --------\n");
    for (int i = 0; i <= side.num_feat; i++)
    {
        printf("\t  %2d     %1.4f   %1.4f   %1.4f   %1.4f   %1.4f   %1.4f   %1.4f\n",
               i, results.scL2_FOMs(i), results.scL1_FOMs(i), results.scmu_FOMs(i),
               results.scdev_FOMs(i), results.scrange_FOMs(i), results.sccorr_FOMs(i),
               results.scinf_FOMs(i));
    }
    std::cout << std::endl;

    if (!(prm.tell&ONLY_STRUCTURAL))
    {
        // Show histogram based FOMS ............................................
        std::cout << "Histogram based --------------------------------------------\n";
        printf("hsin        FOM:%f\n", results.hsmu_FOM);
        printf("hsbr        FOM:%f\n", results.hsbr_FOM);
        printf("hsdt        FOM:%f\n", results.hsdt_FOM);
        printf("hsvr        FOM:%f\n", results.hsvr_FOM);
        printf("\tFEATURE    hsin       hsbr       hsdt       hsvr\n");
        printf("\t------- ---------- ---------- ---------- ----------\n");
        for (int i = 0; i <= side.num_feat; i++)
        {
            printf("\t  %2d     % 7.2f    % 7.2f    % 7.2f    % 7.2f\n",
                   i, results.hsmu_FOMs(i),  results.hsbr_FOMs(i),
                   results.hsdt_FOMs(i), results.hsvr_FOMs(i));
        }
        std::cout << std::endl;

        // Show Directional FOMs ................................................
        std::cout << "Directional ------------------------------------------------\n";
        printf("scrt        FOM:%f\n", results.drrt_FOM);

        // Show Distance FOMs ...................................................
        std::cout << "Distance based ---------------------------------------------\n";
        printf("scbl        FOM:%f\n", results.dsbl_FOM);
        printf("scad        FOM:%f\n", results.dsad_FOM);
    }

    // Save maps ............................................................
    if (prm.tell&SAVE_MAPS)
    {
        VolumeXmipp save, error;

        // Save generated phantom
        std::cerr << "Saving generated phantom ...\n";
        side.vol_phantom.write(side.fn_root + "_eval_phantom.vol");

        // Save mask
        std::cerr << "Saving evaluation mask ...\n";
        side.vol_mask.write(side.fn_root + "_eval_mask.vol");

        // Save label map
        std::cerr << "Saving label map ...\n";
        side.vol_label.write(side.fn_root + "_eval_label.vol");

        // Save a map of differences
        save() = error() = side.vol_phantom() - side.vol_recons();
        if (side.descr_mode == MATH_PHANTOM) side.phantom_descr.sketch_in(&save);
        std::cerr << "Saving difference map ...\n";
        save.write(side.fn_root + "_eval_difference_map.vol");

        // Save a map of quadratic errors
        save() = error() = error() * error();
        if (side.descr_mode == MATH_PHANTOM) side.phantom_descr.sketch_in(&save);
        std::cerr << "Saving quadratic errors ...\n";
        save.write(side.fn_root + "_evaluateQuadratic_map.vol");

        // Save a absolute difference map
        save() = ABSnD(side.vol_phantom() - side.vol_recons());
        if (side.descr_mode == MATH_PHANTOM) side.phantom_descr.sketch_in(&save);
        std::cerr << "Saving absolute difference map ...\n";
        save.write(side.fn_root + "_eval_absolute_map.vol");

        if (!(prm.tell&ONLY_STRUCTURAL))
        {
            // Save distance map
            std::cerr << "Saving distance map ...\n";
            side.vol_distance.write(side.fn_root + "_eval_distance_map.vol");

            // Save blurring map
            save().resize(side.vol_distance());
            FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(side.vol_distance))
            if (side.vol_distance(k, i, j) == -1) save(k, i, j) = 0;
            else save(k, i, j) = 1 / side.vol_distance(k, i, j) * error(k, i, j);
            std::cerr << "Saving blurring map ...\n";
            save.write(side.fn_root + "_eval_blurring_map.vol");

            // Save appearing map
            FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(side.vol_distance))
            if (side.vol_distance(k, i, j) == -1) save(k, i, j) = 0;
            else save(k, i, j) = side.vol_distance(k, i, j) * error(k, i, j);
            std::cerr << "Saving appearence map ...\n";
            save.write(side.fn_root + "_eval_appearing_map.vol");
        }
    }

    // Show values ..........................................................
    if (prm.tell&SHOW_VALUES)
    {
        int sel_feat = 0;
        std::string fn_out;
        std::ofstream fh_out;

        std::cout << "Name of the filename to dump values: ";
        std::cin >> fn_out;
        fh_out.open(fn_out.c_str(), std::ios::out);
        if (!fh_out)
            REPORT_ERROR(3005, (std::string)"Evaluate show: Could not open " + fn_out
                         + " for output");

        while (sel_feat != 1000)
        {
            std::cout << "What feature do you want to see (" << -side.num_feat << ","
            << side.num_feat << ") (1000=finish): ";
            std::cin  >> sel_feat;
            if (ABS(sel_feat) <= side.num_feat)
                switch (side.descr_mode)
                {
                case XMIPP_PHANTOM:
                    show_voxels_in_feat(&side.vol_phantom, &side.vol_recons,
                                        &side.vol_label, sel_feat, fh_out);
                    break;
                case MATH_PHANTOM:
                    show_voxels_in_feat(&side.vol_phantom, &side.vol_recons,
                                        &side.vol_label, side.phantom_descr, sel_feat, fh_out);
                    break;
                }
        }
        fh_out.close();
    }
}

/* Single step ============================================================= */
void Evaluate_single_step(const Prog_Evaluate_Parameters &prm,
                          EVALUATE_results &results)
{
    std::cout << prm;

// Read volumes, label them and generate global mask
    EVALUATE_Side_Info side;
    side.produce_Side_Info(prm);

// Compute FOMs
    compute_FOMs(prm, side, results);

// Show results
    show_FOMs(prm, side, results);
}

/* Main routine ============================================================ */
void ROUT_Evaluate(Prog_Evaluate_Parameters &prm,
                   EVALUATE_results &results)
{
    if (prm.fn_sel == "") Evaluate_single_step(prm, results);
    else
    {
        SelFile SF(prm.fn_sel);
        SelFile SFmask;
        if (prm.fn_mask != "") SFmask.read(prm.fn_mask);
        FOMs foms(SF.ImgNo()), foms_mean(1), foms_stddev(1);
        int k = 0;
        bool mathematical_phantom;
        while (!SF.eof())
        {
            std::cerr << "Perfoming measure for test number " << k << std::endl;
            prm.fn_recons = SF.NextImg();
            if (prm.fnPhantom == "")
            {
                mathematical_phantom = true;
                prm.fnPhantom = prm.fn_recons.without_extension();
                prm.fnPhantom = prm.fnPhantom.without("_wos");
                int i = prm.fnPhantom.find("_idr");
                if (i != -1) prm.fnPhantom = prm.fnPhantom.replace(i, 6, "");
                prm.fnPhantom += ".descr";
            }
            else mathematical_phantom = false;
            if (!SFmask.eof()) prm.fn_mask = SFmask.NextImg();
            else               prm.fn_mask = "";
            Evaluate_single_step(prm, results);
            foms.set_FOMs(k, results);
            k++;
            if (mathematical_phantom)
            {
                prm.fnPhantom = "";
                mathematical_phantom = false;
            }
        }
        compute_FOMs_stats(foms, 0, foms_mean, foms_stddev);
        std::cout << foms;
        std::cout << "--------------------------------------------------\n";
        std::cout << "After combining " << SF.ImgNo() << " volumes\n";
        show_stats(std::cout, 0, foms_mean, foms_stddev);
    }
}

/* FOMs constructor ======================================================== */
FOMs::FOMs(int n)
{
    scL2.resize(n);
    scL1.resize(n);
    scL2w.resize(n);
    scL1w.resize(n);
    scmu.resize(n);
    scdev.resize(n);
    scrange.resize(n);
    sccorr.resize(n);
    scinf.resize(n);
    scresol.resize(n);
    scL20.resize(n);
    scL10.resize(n);
    scmu0.resize(n);
    scdev0.resize(n);
    scrange0.resize(n);
    scL21.resize(n);
    scL11.resize(n);
    scmu1.resize(n);
    scdev1.resize(n);
    scrange1.resize(n);
    hsvr.resize(n);
    hsmu.resize(n);
    hsbr.resize(n);
    hsdt.resize(n);
    drrt.resize(n);
    dsbl.resize(n);
    dsad.resize(n);
}

/* Set ===================================================================== */
void FOMs::set_FOMs(int k, EVALUATE_results &results)
{
    scL2(k)     = results.scL2_FOM;
    scL1(k)     = results.scL1_FOM;
    scL2w(k)    = results.scL2w_FOM;
    scL1w(k)    = results.scL1w_FOM;
    scmu(k)     = results.scmu_FOM;
    scdev(k)    = results.scdev_FOM;
    scrange(k)  = results.scrange_FOM;
    sccorr(k)   = results.sccorr_FOM;
    scinf(k)    = results.scinf_FOM;
    scresol(k)  = results.resol_FOM;
    scL20(k)    = results.scL2_FOMs(0);
    scL10(k)    = results.scL1_FOMs(0);
    scmu0(k)    = results.scmu_FOMs(0);
    scdev0(k)   = results.scdev_FOMs(0);
    scrange0(k) = results.scrange_FOMs(0);
    if (XSIZE(results.scL2_FOMs) != 1)    scL21(k)    = results.scL2_FOMs(1);
    if (XSIZE(results.scL1_FOMs) != 1)    scL11(k)    = results.scL1_FOMs(1);
    if (XSIZE(results.scmu_FOMs) != 1)    scmu1(k)    = results.scmu_FOMs(1);
    if (XSIZE(results.scdev_FOMs) != 1)   scdev1(k)   = results.scdev_FOMs(1);
    if (XSIZE(results.scrange_FOMs) != 1) scrange1(k) = results.scrange_FOMs(1);
    hsvr(k)     = results.hsvr_FOM;
    hsmu(k)     = results.hsmu_FOM;
    hsbr(k)     = results.hsbr_FOM;
    hsdt(k)     = results.hsdt_FOM;
    drrt(k)     = results.drrt_FOM;
    dsbl(k)     = results.dsbl_FOM;
    dsad(k)     = results.dsad_FOM;
}

/* Compute stats =========================================================== */
void compute_FOMs_stats(const FOMs &foms, int i, FOMs &fmean, FOMs &fstddev)
{
    double avg, stddev;
    stats__1(foms.scL2, avg, stddev);
    fmean.scL2(i)    = avg;
    fstddev.scL2(i)    = stddev;
    stats__1(foms.scL1, avg, stddev);
    fmean.scL1(i)    = avg;
    fstddev.scL1(i)    = stddev;
    stats__1(foms.scL2w, avg, stddev);
    fmean.scL2w(i)   = avg;
    fstddev.scL2w(i)   = stddev;
    stats__1(foms.scL1w, avg, stddev);
    fmean.scL1w(i)   = avg;
    fstddev.scL1w(i)   = stddev;
    stats__1(foms.scmu, avg, stddev);
    fmean.scmu(i)    = avg;
    fstddev.scmu(i)    = stddev;
    stats__1(foms.scdev, avg, stddev);
    fmean.scdev(i)   = avg;
    fstddev.scdev(i)   = stddev;
    stats__1(foms.scrange, avg, stddev);
    fmean.scrange(i) = avg;
    fstddev.scrange(i) = stddev;
    stats__1(foms.sccorr, avg, stddev);
    fmean.sccorr(i)  = avg;
    fstddev.sccorr(i)  = stddev;
    stats__1(foms.scinf, avg, stddev);
    fmean.scinf(i)   = avg;
    fstddev.scinf(i)   = stddev;
    stats__1(foms.scresol, avg, stddev);
    fmean.scresol(i) = avg;
    fstddev.scresol(i) = stddev;
    stats__1(foms.scL20, avg, stddev);
    fmean.scL20(i)   = avg;
    fstddev.scL20(i)   = stddev;
    stats__1(foms.scL10, avg, stddev);
    fmean.scL10(i)   = avg;
    fstddev.scL10(i)   = stddev;
    stats__1(foms.scmu0, avg, stddev);
    fmean.scmu0(i)   = avg;
    fstddev.scmu0(i)   = stddev;
    stats__1(foms.scdev0, avg, stddev);
    fmean.scdev0(i)  = avg;
    fstddev.scdev0(i)  = stddev;
    stats__1(foms.scrange0, avg, stddev);
    fmean.scrange0(i) = avg;
    fstddev.scrange0(i) = stddev;
    stats__1(foms.scL21, avg, stddev);
    fmean.scL21(i)   = avg;
    fstddev.scL21(i)   = stddev;
    stats__1(foms.scL11, avg, stddev);
    fmean.scL11(i)   = avg;
    fstddev.scL11(i)   = stddev;
    stats__1(foms.scmu1, avg, stddev);
    fmean.scmu1(i)   = avg;
    fstddev.scmu1(i)   = stddev;
    stats__1(foms.scdev1, avg, stddev);
    fmean.scdev1(i)  = avg;
    fstddev.scdev1(i)  = stddev;
    stats__1(foms.scrange1, avg, stddev);
    fmean.scrange1(i) = avg;
    fstddev.scrange1(i) = stddev;

    stats__1(foms.hsvr, avg, stddev);
    fmean.hsvr(i)    = avg;
    fstddev.hsvr(i)    = stddev;
    stats__1(foms.hsmu, avg, stddev);
    fmean.hsmu(i)    = avg;
    fstddev.hsmu(i)    = stddev;
    stats__1(foms.hsbr, avg, stddev);
    fmean.hsbr(i)    = avg;
    fstddev.hsbr(i)    = stddev;
    stats__1(foms.hsdt, avg, stddev);
    fmean.hsdt(i)    = avg;
    fstddev.hsdt(i)    = stddev;

    stats__1(foms.drrt, avg, stddev);
    fmean.drrt(i)    = avg;
    fstddev.drrt(i)    = stddev;

    stats__1(foms.dsbl, avg, stddev);
    fmean.dsbl(i)    = avg;
    fstddev.dsbl(i)    = stddev;
    stats__1(foms.dsad, avg, stddev);
    fmean.dsad(i)    = avg;
    fstddev.dsad(i)    = stddev;
}

/* Show ==================================================================== */
std::ostream & operator << (std::ostream &out, const FOMs &foms)
{
    out << "All results in all tests\n";
    out << "--------------------------------------------------\n";
    out << "scL2\n"     << foms.scL2.transpose()     << std::endl;
    out << "scL1\n"     << foms.scL1.transpose()     << std::endl;
    out << "scL2w\n"    << foms.scL2w.transpose()    << std::endl;
    out << "scL1w\n"    << foms.scL1w.transpose()    << std::endl;
    out << "scmu\n"     << foms.scmu.transpose()     << std::endl;
    out << "scdev\n"    << foms.scdev.transpose()    << std::endl;
    out << "scrange\n"  << foms.scrange.transpose()  << std::endl;
    out << "sccorr\n"   << foms.sccorr.transpose()   << std::endl;
    out << "scinf\n"    << foms.scinf.transpose()    << std::endl;
    out << "scresol\n"  << foms.scresol.transpose()  << std::endl;
    out << "scL20\n"    << foms.scL20.transpose()    << std::endl;
    out << "scL10\n"    << foms.scL10.transpose()    << std::endl;
    out << "scmu0\n"    << foms.scmu0.transpose()    << std::endl;
    out << "scdev0\n"   << foms.scdev0.transpose()   << std::endl;
    out << "scrange0\n" << foms.scrange0.transpose() << std::endl;
    out << "scL21\n"    << foms.scL21.transpose()    << std::endl;
    out << "scL11\n"    << foms.scL11.transpose()    << std::endl;
    out << "scmu1\n"    << foms.scmu1.transpose()    << std::endl;
    out << "scdev1\n"   << foms.scdev1.transpose()   << std::endl;
    out << "scrange1\n" << foms.scrange1.transpose() << std::endl;

    out << "hsvr\n"    << foms.hsvr.transpose()    << std::endl;
    out << "hsmu\n"    << foms.hsmu.transpose()    << std::endl;
    out << "hsbr\n"    << foms.hsbr.transpose()    << std::endl;
    out << "hsdt\n"    << foms.hsdt.transpose()    << std::endl;

    out << "drrt\n"    << foms.drrt.transpose()    << std::endl;

    out << "dsbl\n"    << foms.dsbl.transpose()    << std::endl;
    out << "dsad\n"    << foms.dsad.transpose()    << std::endl;
    out << "--------------------------------------------------\n";
    return out;
}

/* Show stats ============================================================== */
void show_stats(std::ostream &out, int i, const FOMs &fmean,
                const FOMs &fstddev)
{
    out << "    scL2:     " << fmean.scL2(i)     << "+-" << fstddev.scL2(i)    << std::endl;
    out << "    scL1:     " << fmean.scL1(i)     << "+-" << fstddev.scL1(i)    << std::endl;
    out << "    scL2w:    " << fmean.scL2w(i)    << "+-" << fstddev.scL2w(i)   << std::endl;
    out << "    scL1w:    " << fmean.scL1w(i)    << "+-" << fstddev.scL1w(i)   << std::endl;
    out << "    scmu:     " << fmean.scmu(i)     << "+-" << fstddev.scmu(i)    << std::endl;
    out << "    scdev:    " << fmean.scdev(i)    << "+-" << fstddev.scdev(i)   << std::endl;
    out << "    scrange:  " << fmean.scrange(i)  << "+-" << fstddev.scrange(i) << std::endl;
    out << "    sccorr:   " << fmean.sccorr(i)   << "+-" << fstddev.sccorr(i)  << std::endl;
    out << "    scinf:    " << fmean.scinf(i)    << "+-" << fstddev.scinf(i)   << std::endl;
    out << "    scresol:  " << fmean.scresol(i)  << "+-" << fstddev.scresol(i) << std::endl;
    out << "    scL20:    " << fmean.scL20(i)    << "+-" << fstddev.scL20(i)   << std::endl;
    out << "    scL10:    " << fmean.scL10(i)    << "+-" << fstddev.scL10(i)   << std::endl;
    out << "    scmu0:    " << fmean.scmu0(i)    << "+-" << fstddev.scmu0(i)   << std::endl;
    out << "    scdev0:   " << fmean.scdev0(i)   << "+-" << fstddev.scdev0(i)  << std::endl;
    out << "    scrange0: " << fmean.scrange0(i) << "+-" << fstddev.scrange0(i) << std::endl;
    out << "    scL21:    " << fmean.scL21(i)    << "+-" << fstddev.scL21(i)   << std::endl;
    out << "    scL11:    " << fmean.scL11(i)    << "+-" << fstddev.scL11(i)   << std::endl;
    out << "    scmu1:    " << fmean.scmu1(i)    << "+-" << fstddev.scmu1(i)   << std::endl;
    out << "    scdev1:   " << fmean.scdev1(i)   << "+-" << fstddev.scdev1(i)  << std::endl;
    out << "    scrange1: " << fmean.scrange1(i) << "+-" << fstddev.scrange1(i) << std::endl;
    out << "    scbl:     " << fmean.dsbl(i)     << "+-" << fstddev.dsbl(i)    << std::endl;
    out << "    scap:     " << fmean.dsad(i)     << "+-" << fstddev.dsad(i)    << std::endl;
    out << "    scrt:     " << fmean.drrt(i)     << "+-" << fstddev.drrt(i)    << std::endl;
    out << "    hsin:     " << fmean.hsmu(i)     << "+-" << fstddev.hsmu(i)    << std::endl;
    out << "    hsbr:     " << fmean.hsbr(i)     << "+-" << fstddev.hsbr(i)    << std::endl;
    out << "    hsdt:     " << fmean.hsdt(i)     << "+-" << fstddev.hsdt(i)    << std::endl;
    out << "    hsvr:     " << fmean.hsvr(i)     << "+-" << fstddev.hsvr(i)    << std::endl;
    out.flush();
}
