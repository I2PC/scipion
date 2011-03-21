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


/* Here are all basic functions which are not extra_parameter dependent for
   the ART process. The extra parameter dependent functions are implemented
   in the Basic_art.inc file and must be included in each specific
   implementation (single particles, crystals, ...) */

#include "basic_art.h"
#include <data/projection.h>
#include "fourier_filter.h"
#include <data/funcs.h>
#include "symmetrize.h"
#include <data/histogram.h>
#include <data/mask.h>
#include <data/wavelet.h>


/* Default values ========================================================== */
void Basic_ART_Parameters::default_values()
{
    fh_hist            = NULL;
    fn_start           = "";
    fn_sym             = "";
    force_sym          = 0;
    do_not_generate_subgroup = false;
    do_not_use_symproj = false;
    fn_surface_mask    = "";
    parallel_mode      = ART;
    block_size        = 1;
    eq_mode            = ARTK;
    random_sort        = false;
    dont_sort          = false;
    sort_last_N        = 2;
    WLS                = false;
    no_it              = 1;
    kappa_list.resize(1);
    kappa_list.initConstant(0.5);
    lambda_list.resize(1);
    lambda_list.initConstant(0.01);
    stop_at            = 0;
    basis.set_default();
    grid_relative_size = 1.41;
    grid_type          = BCC;
    proj_ext           = 0;
    Xoutput_volume_size = 0;
    Youtput_volume_size = 0;
    Zoutput_volume_size = 0;
    R                  = -1;
    print_system_matrix = false;
    tell               = 0;
    save_intermidiate_every = 0;
    is_crystal         = false;
    variability_analysis = false;
    refine             = false;
    noisy_reconstruction = false;

    IMG_Inf            = NULL;
    D                  = NULL;
    Dinv               = NULL;
    GVNeq              = NULL;

    surface_mask       = NULL;
    POCS_freq          = 1;

    known_volume       = -1;
    positivity         = false;
    unmatched          = false;
    ray_length         = -1;
    apply_shifts       = true;

    sampling           = 1.;
    sym_each           = 0;
    ref_trans_after    = -1;
    ref_trans_step     = -1;
    sparseEps          = -1;
    diffusionWeight    = -1;
    max_tilt           = 10.e6;
    grid_relative_size = 1.41;
    fn_control         = "";

    threads            = 1;
}

/* Read ART parameters ===================================================== */
#define GET_PARAM_WITH_DEF(flag,default_value) \
    getParameter(argc,argv,"-"flag,default_value)
#define GET_PARAM(flag) \
    getParameter(argc,argv,"-"flag)
#define CHECK_PARAM(flag) \
    checkParameter(argc,argv,"-"flag)
#define GET_VECTOR_PARAM(flag,length) \
    getVectorParameter(argc, argv, "-"flag,length)

#define GET_ART_PARAMS \
    default_values(); \
    fn_sel             =      GET_PARAM(         "i"                    ); \
    fn_ctf             =      GET_PARAM_WITH_DEF("CTF",     ""          ); \
    unmatched          =      CHECK_PARAM(       "unmatched"            );  \
    if (CHECK_PARAM("o")) \
        fn_root        =      GET_PARAM(         "o"                    );  \
    else fn_root       =      fn_sel.withoutExtension();                   \
    fn_start           =      GET_PARAM_WITH_DEF("start",     ""        );  \
    if      (CHECK_PARAM("pSART"))  parallel_mode=pSART;\
    else if (CHECK_PARAM("pSIRT"))  parallel_mode=pSIRT; \
    else if (CHECK_PARAM("SIRT"))   parallel_mode=SIRT; \
    else if (CHECK_PARAM("pfSIRT")) parallel_mode=pfSIRT; \
    else if (CHECK_PARAM("pBiCAV")) parallel_mode=pBiCAV; \
    else if (CHECK_PARAM("pAVSP"))  parallel_mode=pAVSP; \
    else if (CHECK_PARAM("pCAV"))   parallel_mode=pCAV; \
    else                            parallel_mode=ART; \
    ray_length         = textToInteger(GET_PARAM_WITH_DEF("ray_length","-1"      )); \
    block_size         = textToInteger(GET_PARAM_WITH_DEF("block_size","1" )); \
    fn_sym             =      GET_PARAM_WITH_DEF("sym",       ""        );  \
    fn_control         =      GET_PARAM_WITH_DEF("control",       ""        );  \
    force_sym          = textToInteger(GET_PARAM_WITH_DEF("force_sym","0"        )); \
    do_not_generate_subgroup= CHECK_PARAM(       "no_group"             );  \
    do_not_use_symproj = CHECK_PARAM(       "no_symproj"             );  \
    shiftedTomograms   =      CHECK_PARAM(       "shiftedTomograms"     ); \
    goldmask    = textToFloat(GET_PARAM_WITH_DEF("goldmask","1e6")); \
    fn_surface_mask    =      GET_PARAM_WITH_DEF("surface",   ""        );  \
    random_sort        =      CHECK_PARAM(       "random_sort"          );  \
    dont_sort          =      CHECK_PARAM(       "no_sort"          );  \
    sort_last_N        = textToInteger(GET_PARAM_WITH_DEF("sort_last", "2"       )); \
    WLS                = CHECK_PARAM("WLS"); \
    kappa_list  =  GET_VECTOR_PARAM("k",-1); \
    if (XSIZE(kappa_list)==0)                                              \
    {kappa_list.resize(1); kappa_list.initConstant(0.5);}               \
    no_it              = textToInteger(GET_PARAM_WITH_DEF("n",         "1"       )); \
    stop_at            = textToInteger(GET_PARAM_WITH_DEF("stop_at",   "0"       )); \
    lambda_list        =        GET_VECTOR_PARAM("l",         -1);          \
    if (XSIZE(lambda_list)==0)                                              \
    {lambda_list.resize(1); lambda_list.initConstant(0.01);}               \
    sampling           = textToFloat(GET_PARAM_WITH_DEF("sampling",  "1."        )); \
    sym_each           = textToInteger(GET_PARAM_WITH_DEF("sym_each",  "0"         )); \
    max_tilt           = textToFloat(GET_PARAM_WITH_DEF("max_tilt",  "10E6"      )); \
    ref_trans_after    = textToInteger(GET_PARAM_WITH_DEF("ref_trans_after", "-1"  )); \
    ref_trans_step     = textToFloat(GET_PARAM_WITH_DEF("ref_trans_step", "-1"    )); \
    sparseEps          = textToFloat(GET_PARAM_WITH_DEF("sparse", "-1"  )); \
    diffusionWeight    = textToFloat(GET_PARAM_WITH_DEF("diffusion", "-1"  )); \
    grid_relative_size = textToFloat(GET_PARAM_WITH_DEF("g",         "1.41"      )); \
    R                  = textToFloat(GET_PARAM_WITH_DEF("R",         "-1"        )); \
    POCS_freq          = textToInteger(GET_PARAM_WITH_DEF("POCS_freq", "1"         )); \
    known_volume       = textToFloat(GET_PARAM_WITH_DEF("known_volume","-1"      )); \
    positivity         = CHECK_PARAM("POCS_positivity"); \
    apply_shifts       = !CHECK_PARAM("dont_apply_shifts"); \
    threads            = textToInteger(GET_PARAM_WITH_DEF("thr", "1"         )); \
    if      (grid_relative_size == -1)  grid_relative_size = sqrt (2.0); \
    else if (grid_relative_size == -2)  grid_relative_size = pow (2.0,1.0/3.0); \
    \
    if      (CHECK_PARAM("CAVK"))    eq_mode=CAVK; \
    else if (CHECK_PARAM("CAV"))     eq_mode=CAV; \
    else if (CHECK_PARAM("CAVARTK")) eq_mode=CAVARTK; \
    else                             eq_mode=ARTK; \
    \
    if      (CHECK_PARAM("FCC")) grid_type=FCC; \
    else if (CHECK_PARAM("CC"))  grid_type=CC; \
    else                         grid_type=BCC; \
    proj_ext           = textToInteger(GET_PARAM_WITH_DEF("ext",       "0"    )); \
    \
    if (CHECK_PARAM("small_blobs"))  grid_relative_size=1.41; \
    if (CHECK_PARAM("big_blobs"))    grid_relative_size=2.26; \
    if (CHECK_PARAM("visual_blobs")) grid_relative_size=1.41; \
    if (CHECK_PARAM("voxels")) { \
        grid_relative_size=1; \
        grid_type=CC; \
    } \
    if (CHECK_PARAM("splines")) { \
        grid_relative_size=1; \
        grid_type=CC; \
    } \
    \
    print_system_matrix=CHECK_PARAM("print_system_matrix"); \
    if (CHECK_PARAM("show_error"))        tell |= TELL_SHOW_ERROR; \
    if (CHECK_PARAM("manual_order"))      tell |= TELL_MANUAL_ORDER; \
    if (CHECK_PARAM("only_sym"))          tell |= TELL_ONLY_SYM; \
    if (CHECK_PARAM("save_at_each_step")) tell |= TELL_SAVE_AT_EACH_STEP; \
    if (CHECK_PARAM("save_basis"))        tell |= TELL_SAVE_BASIS; \
    if (CHECK_PARAM("show_stats"))        tell |= TELL_STATS; \
    if (CHECK_PARAM("save_intermidiate")) {\
        tell |= TELL_SAVE_INTERMIDIATE; \
        save_intermidiate_every=textToInteger(GET_PARAM_WITH_DEF("save_intermidiate","0")); \
    } \
    if (CHECK_PARAM("show_iv")) {\
        tell |= TELL_IV; \
        save_intermidiate_every=textToInteger(GET_PARAM_WITH_DEF("show_iv","10")); \
    } \
    if (CHECK_PARAM("variability")) {\
        variability_analysis=true; \
        parallel_mode=SIRT; \
        no_it=1; \
    } \
    refine         = CHECK_PARAM("refine"); \
    if (CHECK_PARAM("noisy_reconstruction")) { \
        if (parallel_mode!=ART) \
            REPORT_ERROR(ERR_ARG_INCORRECT,"Basic_ART_Parameters::read: Noisy reconstructions" \
                         " can only be done for ART"); \
        else noisy_reconstruction=true; \
    }

void Basic_ART_Parameters::read(int argc, char **argv)
{
    GET_ART_PARAMS;
    basis.read(argc, argv);
    //divide by the sampling rate
    if (sampling != 1.)
    {
        basis.set_sampling_rate(sampling);
        grid_relative_size /= sampling;
        if (R != -1.)
            R /= sampling;
        ref_trans_step /= sampling;
    }
    if (CHECK_PARAM("output_size"))
    {
        int i = paremeterPosition(argc, argv, "-output_size");
        if (i + 3 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "Not enough parameters after -output_size");
        Zoutput_volume_size = textToInteger(argv[i+1]);
        Youtput_volume_size = textToInteger(argv[i+2]);
        Xoutput_volume_size = textToInteger(argv[i+3]);
    }
}
#undef GET_PARAM_WITH_DEF
#undef GET_PARAM
#undef CHECK_PARAM
#undef GET_VECTOR_PARAM

#define GET_PARAM_WITH_DEF(flag,default_value) \
    getParameter(fh,flag,0,default_value)
#define GET_PARAM(flag) \
    getParameter(fh,flag,0)
#define CHECK_PARAM(flag) \
    checkParameter(fh,flag)
#define GET_VECTOR_PARAM(flag,length) \
    getVectorParameter(fh,flag,length)
// Read from file
void Basic_ART_Parameters::read(const FileName &fn)
{
    FILE *fh;
    if ((fh = fopen(fn.c_str(), "r")) == NULL)
        REPORT_ERROR(ERR_IO_NOTEXIST,fn);

    GET_ART_PARAMS;
    if (CHECK_PARAM("output_size"))
    {
        int argcp;
        char **argvp = NULL, *copy = NULL;
        generateCommandLine(fh, "output_size", argcp, argvp, copy);
        int i = paremeterPosition(argcp, argvp, "-output_size");
        if (i + 3 >= argcp)
            REPORT_ERROR(ERR_ARG_MISSING, "Not enough parameters after -output_size");
        Zoutput_volume_size = textToInteger(argvp[i+1]);
        Youtput_volume_size = textToInteger(argvp[i+2]);
        Xoutput_volume_size = textToInteger(argvp[i+3]);
    }
    fclose(fh);
    basis.read(fn);
}

/* Usage =================================================================== */
void Basic_ART_Parameters::usage()
{
    std::cerr
    << "Usage: art [Options and Parameters]"
    << "\nOptions:"
    << "\nParameter Values: (note space before value)"
    << "\n    -i selfile           full name of sel file"
    << "\n   [-o name]             name of output files, extensions are added"
    << "\n   [-sym symmfile]       Use a symmetry file"
    << "\n   [-n noit=1]           number of iterations"
    << "\n   [-l lambda=0.01]      relaxation factor (recommended range 0.0 - 0.1)"
    << "\n   [-more_help]          show all parameters"
    << "\n"
    ;
}

void Basic_ART_Parameters::usage_more()
{
    std::cerr
    << "Usage: art [Options and Parameters]"
    << "\nOptions:"
    << "\nParameter Values: (note space before value)"
    << "\nI/O parameters"
    << "\n    -i selfile           full name of sel file"
    << "\n   [-o name]             name of output files, extensions are added"
    << "\n   [-CTF name]           name of a sel file with CTFs"
    << "\n   [-unmatched]          apply unmatched forward/backward projectors"
    << "\n   [-start basisvolume]  Start from basisvolume"
    << "\n   [-sym symmfile]       Use a symmetry file"
    << "\n   [-sym_each n]         Force the reconstruction to be symmetric"
    << "\n                         each n projections"
    << "\n   [-max_tilt n]         Skip projection with absolute tilt angle"
    << "\n                         greater than n\n"
    << "\n   [-ref_trans_after n]  Refine the translation alignment"
    << "\n                         after n projections."
    << "\n   [-ref_trans_step n]   Max displacement in translation alignment"
    << "\n                         This is a double."
    << "\n   [-force_sym <n=0>]    Force the reconstruction to be symmetric"
    << "\n                         n times at each projection"
    << "\n   [-no_group]           Do not generate symmetry subgroup"
    << "\n   [-no_symproj]         Do not use symmetrized projections"
    << "\n   [-sparse <eps=-1>]    Sparsity threshold"
    << "\n   [-diffusion <eps=-1>] Diffusion weight"
    << "\n   [-surface surf_mask]  Use this file as a surface mask"
    << "\n   [-POCS_freq <f=1>]    Impose POCS conditions every <f> projections"
    << "\n   [-known_volume <vol=-1>] Volume of the reconstruction"
    << "\n   [-POCS_positivity]    Apply positivity constraint"
    << "\n   [-goldmask value]     Pixels below this value are not considered for reconstruction"
    << "\n   [-shiftedTomograms]   Remove border pixels created by alignment of tomograms"
    << "\n   [-dont_apply_shifts]  Do not apply shifts as stored in the 2D-image headers"
    << "\n   [-variability]        Perform variability analysis"
    << "\n   [-refine]             Refine input projection before backprojecting"
    << "\n   [-noisy_reconstruction] Perform a companion noisy reconstruction"
    << "\n   [-ray_length <r=-1>]  In basis units\n"
    ;
    std::cerr
    << "\nIteration parameters"
    << "\n   [-n noit=1]           number of iterations"
    << "\n   [-stop_at stop_at=0]  number of images presented"
    << "\n   [-l lambda=0.01 |     relaxation factor (recommended range 0.0 - 0.1)"
    << "\n    -l [lambda0, lambda1, ...]"
    << "\n   [-CAVK|-CAV]          by default, ARTK is applied"
    << "\n   [-sort_last N=2]      Use -1 to sort with all previous projections"
    << "\n   [-random_sort]        by default, perpendicular sort is used for ART"
    << "\n   [-no_sort]            No sort must be applied"
    << "\n   [-WLS]                Perform weighted least squares ART"
    << "\n   [-k kappa=0.5 |       Relaxation factor for WLS residual "
    << "\n    -k [kappa0, kappa1, ...]\n"
    << "\nParallel parameters"
    << "\n                         by default, sequential ART is applied"
    << "\n  [-thr N=1]      Number of threads to use. NOTE: Not available when using MPI."
    << "\n   [-SIRT]               Simultaneous Iterative Reconstruction Technique"
    << "\n   [-pSIRT]              Parallel (MPI) Simultaneous Iterative Reconstruction Technique"
    << "\n   [-pfSIRT]             Parallel (MPI) False Simultaneous Iterative Reconstruction Technique (Faster convergence than pSIRT)"
    << "\n   [-pSART]              Parallel (MPI) Simultaneous ART\n"
    << "\n   [-pAVSP]              Parallel (MPI) Average Strings\n"
    << "\n   [-pBiCAV]             Parallel (MPI) Block Iterative CAV\n"
    << "\n   [-pCAV]               Parallel (MPI) CAV\n"
    << "\n   [-block_size <n=1>]   Number of projections to each block (SART and BiCAV)\n"
    << "\n   [-CAVARTK]            Component Averaging Variant of Block ART\n"
    << "\nGrid parameters"
    << "\n   [-g gridsz=1.41]      relative grid size"
    << "\n                         if gridsz =  -1 => gridsz=2^(1/2)"
    << "\n                                      -2 => gridsz=2^(1/3)"
    << "\n   [-FCC]                 use a FCC grid instead of a BCC"
    << "\n   [-SC]                  use a SC grid instead of a BCC"
    << "\n   [-R interest_sphere=-1] Radius of the interest sphere"
    << "\n   [-ext proj_ext=0]     projection extension"
    << "\n   [-output_size Zsize Ysize Xsize] output volume size in PIXELS\n"
    << "\n   [-sampling=1]         sampling rate,  affects to -r, -g, -R and"
    << "\n                          -ref_trans_step"
    << "\n                         Also to -mod_a and mod_b when processing"
    << "\n                         crystals"
    ;
    basis.usage();
    std::cerr
    << "\nDebugging options"
    << "\n   [-print_system_matrix]print the matrix of the system Ax=b"
    << "\n   [-show_iv <n=10>]     show volumes/images as the reconstruction goes"
    << "\n                         the volume is update every <n> projections\n"
    << "\n   [-show_error]         show error for each projection"
    << "\n   [-show_stats]         give some statistical information during the process"
    << "\n   [-save_at_each_step]  save intermidiate projections"
    << "\n                             PPPtheo, PPPread, PPPcorr, PPPdiff"
    << "\n                             PPPbasis.basis, PPPvol.vol"
    << "\n                             PPPvolPOCS1, PPPvolPOCS2, PPPvolPOCS3"
    << "\n   [-save_intermidiate <n>] save intermidiate volumes (every <n> projections)"
    << "\n                             <fnroot>it<no_it>proj<no_projs>.vol"
    << "\n   [-save_basis]         every time you have to save a volume, save it"
    << "\n                         also in basis"
    << "\n   [-manual_order]       manual selection of projection order"
    << "\n   [-only_sym]           skip all those symmetries different from -1"
    << "\n"
    ;
}

/* ------------------------------------------------------------------------- */
/* Sort_perpendicular                                                        */
/* ------------------------------------------------------------------------- */
void sort_perpendicular(int numIMG, Recons_info *IMG_Inf,
                        MultidimArray<int> &ordered_list, int N)
{
    int   i, j, k;
    MultidimArray<short> chosen(numIMG);     // 1 if that image has been already
    // chosen
    double min_prod;
    int   min_prod_proj;
    Matrix2D<double> v(numIMG, 3);
    Matrix2D<double> euler;
    MultidimArray<double> product(numIMG);

    // Initialisation
    ordered_list.resize(numIMG);
    for (i = 0; i < numIMG; i++)
    {
        Matrix1D<double> z;
        // Initially no image is chosen
        A1D_ELEM(chosen, i) = 0;

        // Compute the Euler matrix for each image and keep only
        // the third row of each one
        //0.f -> double 0. It should be there is the other
        // arguments are doubles because Euler_angles2matrix
        //acepts either all doubles or all doubles
        Euler_angles2matrix(IMG_Inf[i].rot, IMG_Inf[i].tilt, 0.f, euler);
        euler.getRow(2, z);
        v.setRow(i, z);
    }

    // Pick first projection as the first one to be presented
    i = 0;
    A1D_ELEM(chosen, i) = 1;
    A1D_ELEM(ordered_list, 0) = i;

    // Choose the rest of projections
    std::cerr << "Sorting projections ...\n";
    init_progress_bar(numIMG - 1);
    Matrix1D<double> rowj, rowi_1, rowi_N_1;
    for (i = 1; i < numIMG; i++)
    {
        // Compute the product of not already chosen vectors with the just
        // chosen one, and select that which has minimum product
        min_prod = MAXFLOAT;
        v.getRow(A1D_ELEM(ordered_list, i - 1),rowi_1);
        if (N != -1 && i > N)
            v.getRow(A1D_ELEM(ordered_list, i - N - 1),rowi_N_1);
        for (j = 0; j < numIMG; j++)
        {
            if (!A1D_ELEM(chosen, j))
            {
                v.getRow(j,rowj);
                A1D_ELEM(product, j) += ABS(dotProduct(rowi_1,rowj));
                if (N != -1 && i > N)
                    A1D_ELEM(product, j) -= ABS(dotProduct(rowi_N_1,rowj));
                if (A1D_ELEM(product, j) < min_prod)
                {
                    min_prod = A1D_ELEM(product, j);
                    min_prod_proj = j;
                }
            }
        }

        // Store the chosen vector and mark it as chosen
        A1D_ELEM(ordered_list, i) = min_prod_proj;
        A1D_ELEM(chosen, min_prod_proj) = 1;

        // The progress bar is updated only every 10 images
        if (i % 10 == 0)
            progress_bar(i);
    }

    // A final call to progress bar to finish a possible small piece
    progress_bar(numIMG - 1);
    std::cout << std::endl;
}

/* ------------------------------------------------------------------------- */
/* No Sort                                                                   */
/* ------------------------------------------------------------------------- */
void no_sort(int numIMG, MultidimArray<int> &ordered_list)
{
    ordered_list.initLinear(0, numIMG - 1);
}

/* ------------------------------------------------------------------------- */
/* Random Sort                                                               */
/* ------------------------------------------------------------------------- */
void sort_randomly(int numIMG, MultidimArray<int> &ordered_list)
{
    int i;
    MultidimArray<int> chosen;

    // Initialisation
    ordered_list.resize(numIMG);
    chosen.initZeros(numIMG);

    std::cerr << "Randomizing projections ...\n";
    init_progress_bar(numIMG - 1);
    int ptr = 0;
    randomize_random_generator();
    for (int i = numIMG; i > 0; i--)
    {
        // Jump a random number starting at the pointed projection
        int rnd_indx = (int) rnd_unif(0, i) + 1;
        while (rnd_indx > 0)
        {
            // Jump one not chosen image
            ptr = (ptr + 1) % numIMG;
            // Check it is not chosen, if it is, go on skipping
            while (chosen(ptr))
            {
                ptr = (ptr + 1) % numIMG;
            }
            rnd_indx--;
        }

        // Annotate this image
        A1D_ELEM(ordered_list, i - 1) = ptr;
        A1D_ELEM(chosen, ptr) = 1;

        // The progress bar is updated only every 10 images
        if (i % 10 == 0)
            progress_bar(i);
    }

    // A final call to progress bar to finish a possible small piece
    progress_bar(numIMG - 1);
    std::cout << std::endl;
}

/* ------------------------------------------------------------------------- */
/* Produce Side Information                                                  */
/* ------------------------------------------------------------------------- */
//#define DEBUG
void Basic_ART_Parameters::produce_Side_Info(GridVolume &vol_basis0, int level,
        int rank)
{
    MetaData     selfile;

    /* If checking the variability --------------------------------------------- */
    if (variability_analysis)
        parallel_mode = SIRT;

    /* Create history file handler --------------------------------------------- */
    if (level >= FULL)
    {
        fh_hist = new std::ofstream;
        fh_hist->open((fn_root + ".hist").c_str(), std::ios::out);
        if (!fh_hist)
            REPORT_ERROR(ERR_IO_NOWRITE, fn_root + ".hist");
    }

    /* Get True Image number and projection size ------------------------------- */
    if (level >= BASIC)
    {
        //take into account weights here
        if (WLS)
        {
            MetaData SF_aux;
            double weight=0.;
            SF_aux.read(fn_sel);
            SF_aux.removeObjects(MDValueEQ(MDL_ENABLED, -1));
            selfile.clear();
            selfile.importObjects(SF_aux, MDValueRange(MDL_WEIGHT, 1e-9, 99e99));
            if (selfile.size() == 0)
            {
                std::cerr << "there is no input file with weight!=0" << std::endl;
                exit(1);
            }
        }
        else
        {
            selfile.read(fn_sel);
            selfile.removeObjects(MDValueEQ(MDL_ENABLED, -1));
        }
        trueIMG = selfile.size();
        if (trueIMG == 0)
            REPORT_ERROR(ERR_MD_OBJECTNUMBER, "Produce_Basic_ART_Side_Info: No images !!");
        int idum;
        size_t idumLong;
        ImgSize(selfile, projXdim, projYdim, idum, idumLong);
    }

    /* Read symmetry file -------------------------------------------------- */
    if (level >= FULL)
    {
        double accuracy = (do_not_generate_subgroup) ? -1 : 1e-6;
        if (fn_sym != "")
            SL.read_sym_file(fn_sym, accuracy);
        if (!do_not_use_symproj)
            numIMG = trueIMG * (SL.SymsNo() + 1);
        else
            numIMG = trueIMG;
    }

    /* Read surface mask --------------------------------------------------- */
    if (level >= FULL)
    {
        if (fn_surface_mask != "")
        {
            surface_mask = new Image<double>;
            surface_mask->read(fn_surface_mask);
            (*surface_mask)().setXmippOrigin();
        }
    }

    /* Fill ART_sort_info structure and Sort ------------------------------- */
    if (level >= FULL)
    {
        build_recons_info(selfile, fn_ctf, SL, IMG_Inf,
                          do_not_use_symproj);

        if (!(tell&TELL_MANUAL_ORDER))
            if (parallel_mode == SIRT ||
                parallel_mode == pSIRT ||
                parallel_mode == pfSIRT ||
                parallel_mode == pCAV ||
                eq_mode == CAV ||
                rank > 0 || dont_sort)
                no_sort(numIMG, ordered_list);
            else if (random_sort)
                sort_randomly(numIMG, ordered_list);
            else if (sort_last_N != -1)
                sort_perpendicular(numIMG, IMG_Inf, ordered_list,
                                   sort_last_N);
            else
                no_sort(numIMG, ordered_list);
    }

    /* In case of weighted least-squares, find average weight & write residual images ------ */
    if (WLS)
    {
        FileName   fn_resi;
        Projection read_proj;
        double     weight;

        sum_weight = 0.;
        for (int iact_proj = 0; iact_proj < numIMG ; iact_proj++)
        {

            fn_resi = IMG_Inf[iact_proj].fn_proj;
            read_proj.read(fn_resi);
            read_proj().setXmippOrigin();
            weight = read_proj.weight();
            if (weight < 0)
                REPORT_ERROR(ERR_VALUE_INCORRECT,
                		     "BASIC_ART: negative weight not set correctly!");
            sum_weight += weight;
            /*
            read_proj().initZeros();
            fn_resi+="."+fn_root+"_residual";
            if (IMG_Inf[iact_proj].sym>-1)
            fn_resi+=integerToString(IMG_Inf[iact_proj].sym);
            read_proj.write(fn_resi);
            */
        }

        *fh_hist << "WLS-ART% Sum over all weights = " << sum_weight << std::endl;
    }

    /* Setting initial volumes ------------------------------------------------- */
    if (level >= FULL)
    {
        if (fn_start != "")
            vol_basis0.read(fn_start, basis.basisName());
        else
        {
            Grid grid_basis;
            if (R == -1)
            {
                Matrix1D<double> corner;
                if (basis.type == Basis::blobs)
                {
                    if (Zoutput_volume_size == 0)
                        corner = vectorR3((double)projXdim / 2, (double)projXdim / 2,
                                          (double)projXdim / 2);
                    else
                        corner = vectorR3(
                                     (double)Xoutput_volume_size / 2,
                                     (double)Youtput_volume_size / 2,
                                     (double)Zoutput_volume_size / 2);
                }
                else
                {
                    if (Zoutput_volume_size == 0)
                        corner = vectorR3(-(double)FIRST_XMIPP_INDEX(projXdim),
                                          -(double)FIRST_XMIPP_INDEX(projXdim),
                                          -(double)FIRST_XMIPP_INDEX(projXdim));
                    else
                        corner = vectorR3(-(double)FIRST_XMIPP_INDEX(Xoutput_volume_size),
                                          -(double)FIRST_XMIPP_INDEX(Youtput_volume_size),
                                          -(double)FIRST_XMIPP_INDEX(Zoutput_volume_size));
                }
                /* If you substract half the basis radius, you are forcing that the
                   last basis touches slightly the volume border. By not substracting
                   it there is a basis center as near the border as possible. */
                corner = corner + proj_ext/*CO: -blob.radius/2*/;
                switch (grid_type)
                {
                case (CC):
                                grid_basis = Create_CC_grid(grid_relative_size, -corner, corner);
                    break;
                case (FCC):
                                grid_basis = Create_FCC_grid(grid_relative_size, -corner, corner);
                    break;
                case (BCC):
                                grid_basis = Create_BCC_grid(grid_relative_size, -corner, corner);
                    break;
                }
            }
            else
    {
                switch (grid_type)
                {
                case (CC):
                                grid_basis = Create_CC_grid(grid_relative_size, R);
                    break;
                case (FCC):
                                grid_basis = Create_FCC_grid(grid_relative_size, R);
                    break;
                case (BCC):
                                grid_basis = Create_BCC_grid(grid_relative_size, R);
                    break;
                }
            }
            vol_basis0.adapt_to_grid(grid_basis);
        }
    }

    /* Basis side info --------------------------------------------------------- */
    if (level >= BASIC)
{
        basis.set_D(D);
        basis.produce_side_info(vol_basis0.grid());
    }

    /* Express the ray length in basis units ----------------------------------- */
    if (ray_length != -1)
        ray_length *= basis.max_length();
}
#undef DEBUG

/* Count number of equations for CAV --------------------------------------- */
void Basic_ART_Parameters::compute_CAV_weights(GridVolume &vol_basis0,
        int numProjs_node, int debug_level)
{
    if (GVNeq == NULL)
        GVNeq = new GridVolumeT<int>;
    GVNeq->resize(vol_basis0);
    GVNeq->initZeros();

    Projection read_proj;
    if (debug_level > 0)
    {
        std::cerr << "Counting equations ...\n";
        init_progress_bar(numIMG);
    }
    for (int act_proj = 0; act_proj < numProjs_node ; act_proj++)
    {
        read_proj.read(IMG_Inf[ordered_list(act_proj)].fn_proj, apply_shifts);
        read_proj().setXmippOrigin();

        // Projection extension? .........................................
        if (proj_ext != 0)
            read_proj().selfWindow(
                STARTINGY(read_proj()) - proj_ext,
                STARTINGX(read_proj()) - proj_ext,
                FINISHINGY(read_proj()) + proj_ext,
                FINISHINGX(read_proj()) + proj_ext);

        count_eqs_in_projection(*GVNeq, basis, read_proj);

        if (debug_level > 0 &&
            act_proj % XMIPP_MAX(1, numIMG / 60) == 0)
            progress_bar(act_proj);
    }
    if (debug_level > 0)
    {
        progress_bar(numIMG);
        long int Neq = 0, Nunk = 0;
        for (int n = 0; n < GVNeq->VolumesNo(); n++)
            FOR_ALL_ELEMENTS_IN_ARRAY3D((*GVNeq)(n)())
        {
            Neq += (*GVNeq)(n)(k, i, j);
            Nunk++;
        }
        std::cerr << "There are " << Neq << " equations and " << Nunk
        << " unknowns (redundancy=" << 100.0 - 100.0*Nunk / Neq << ")\n";
    }
}

int Basic_ART_Parameters::ProjXdim()
{
    return projXdim;
}

int Basic_ART_Parameters::ProjYdim()
{
    return projYdim;
}

/* Fill Reconstruction info structure -------------------------------------- */
void build_recons_info(MetaData &selfile,
                       const FileName &fn_ctf, const SymList &SL,
                       Recons_info * &IMG_Inf, bool do_not_use_symproj)
{
    Matrix2D<double>  L(4, 4), R(4, 4);  // A matrix from the list
    FileName          fn_proj;
    FileName          fn_ctf1;
    Projection        read_proj;
    bool              is_there_ctf = false;
    bool              is_ctf_unique = false;

    int trueIMG = selfile.size();
    selfile.firstObject();
    int numIMG;
    if (!do_not_use_symproj)
        numIMG = trueIMG * (SL.SymsNo() + 1);
    else
        numIMG = trueIMG;

    // The next two ifs check whether there is a CTF file, and
    // whether it is unique
    if (fn_ctf != "")
    {
        is_there_ctf = true;
        is_ctf_unique = true;
    }
    else if (selfile.containsLabel(MDL_CTFMODEL))
    {
        is_there_ctf = true;
        is_ctf_unique = false;
    }

    if (IMG_Inf != NULL)
        delete [] IMG_Inf;
    if ((IMG_Inf = new Recons_info[numIMG]) == NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Build_Recons_Info: No memory for the sorting");

    int i = 0; // It will account for the number of valid projections processed
    std::cerr << "Reading angle information ...\n";
    init_progress_bar(trueIMG);
    FOR_ALL_OBJECTS_IN_METADATA(selfile)
    {
        selfile.getValue(MDL_IMAGE,fn_proj,__iter.objId);
        if (is_there_ctf && !is_ctf_unique)
            selfile.getValue(MDL_CTFMODEL,fn_ctf1,__iter.objId);
        if (fn_proj != "")
        {
            read_proj.read(fn_proj, false, HEADER);
            // Filling structure
            IMG_Inf[i].fn_proj = fn_proj;
            if (is_ctf_unique)
                IMG_Inf[i].fn_ctf = fn_ctf;
            else if (is_there_ctf)
                IMG_Inf[i].fn_ctf = fn_ctf1;
            IMG_Inf[i].sym     = -1;
            IMG_Inf[i].seed    = ROUND(65535 * rnd_unif());
            read_proj.getEulerAngles(IMG_Inf[i].rot, IMG_Inf[i].tilt, IMG_Inf[i].psi);
            EULER_CLIPPING(IMG_Inf[i].rot, IMG_Inf[i].tilt, IMG_Inf[i].psi);

            // Any symmetry?
            if (SL.SymsNo() > 0 && !do_not_use_symproj)
            {
                for (int j = 0; j < SL.SymsNo(); j++)
                {
                    int sym_index = SYMINDEX(SL, j, i, trueIMG);
                    IMG_Inf[sym_index].fn_proj = IMG_Inf[i].fn_proj;
                    IMG_Inf[sym_index].seed   = IMG_Inf[i].seed;
                    if (is_ctf_unique)
                        IMG_Inf[sym_index].fn_ctf = fn_ctf;
                    else if (is_there_ctf)
                        IMG_Inf[sym_index].fn_ctf = IMG_Inf[i].fn_ctf;
                    IMG_Inf[sym_index].sym = j;
                    SL.get_matrices(j, L, R);
                    L.resize(3, 3); // Erase last row and column
                    R.resize(3, 3); // as only the relative orientation
                    // is useful and not the translation
                    double drot, dtilt, dpsi;
                    Euler_apply_transf(L, R,
                                       IMG_Inf[i].rot, IMG_Inf[i].tilt, IMG_Inf[i].psi,
                                       drot, dtilt, dpsi);
                    IMG_Inf[sym_index].rot = (float)drot;
                    IMG_Inf[sym_index].tilt = (float)dtilt;
                    IMG_Inf[sym_index].psi = (float)dpsi;
                }
            }
        }

        ++i; // I have processed one more image
        if (i % 25 == 0)
            progress_bar(i);
    }
    progress_bar(trueIMG);
}

/* ------------------------------------------------------------------------- */
VariabilityClass::VariabilityClass(Basic_ART_Parameters *_prm,
                                   int _Zoutput_volume_size, int _Youtput_volume_size,
                                   int _Xoutput_volume_size)
{
    prm = _prm;
    Zoutput_volume_size = _Zoutput_volume_size;
    Youtput_volume_size = _Youtput_volume_size;
    Xoutput_volume_size = _Xoutput_volume_size;
    N = 0;
    VAR_state = VAR_none;
}

void VariabilityClass::newIteration()
{}

//#define MODE7
#define MODE8
//#define MODE15
//#define DEBUG
void VariabilityClass::newUpdateVolume(GridVolume *ptr_vol_out,
                                       Projection &read_proj)
{
    Image<double> vol_voxels;

    // Convert from basis to voxels
    prm->basis.changeToVoxels(*ptr_vol_out, &(vol_voxels()),
                              Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
    (*ptr_vol_out).initZeros();
    N++;

    // Make the DWT
    MultidimArray<double> DWTV;
#ifdef MODE7

    int DWT_iterations = 2;
    int keep_from_iteration = 0;
#endif
#ifdef MODE8

    int DWT_iterations = 2;
    int keep_from_iteration = 0;
#endif
#ifdef MODE15

    int DWT_iterations = 2;
    int keep_from_iteration = 0;
#endif

    Bilib_DWT(vol_voxels(), DWTV, DWT_iterations);
#ifdef DEBUG

    vol_voxels.write("PPPVariability.vol");
    char c;
    std::cout << "Press any key\n";
    std::cin >> c;
#endif

    // Select the LLL block and keep it
    int x1, x2, y1, y2, z1, z2;
    SelectDWTBlock(keep_from_iteration, DWTV, "000", x1, x2, y1, y2, z1, z2);
    DWTV.selfWindow(z1, y1, x1, z2, y2, x2);
    STARTINGZ(DWTV) = STARTINGY(DWTV) = STARTINGX(DWTV) = 0;
    VA.push_back(DWTV);
}
#undef DEBUG

#define DEBUG
void VariabilityClass::finishAnalysis()
{
    if (VA.size() == 0)
        return;

    // Coocurrence matrix
    int nmax = VA.size();
    MultidimArray<int> coocurrence(nmax, nmax);

    // Study each voxel
    Image<double> SignificantT2, SignificantMaxRatio, SignificantMinRatio;
    int zsize = ZSIZE(VA[0]) / 2;
    int ysize = YSIZE(VA[0]) / 2;
    int xsize = XSIZE(VA[0]) / 2;
    int zsize2 = ZSIZE(VA[0]) / 4;
    int ysize2 = YSIZE(VA[0]) / 4;
    int xsize2 = XSIZE(VA[0]) / 4;
    SignificantT2().initZeros(zsize, ysize, xsize);
    SignificantMaxRatio().initZeros(zsize, ysize, xsize);
    SignificantMinRatio().initZeros(zsize, ysize, xsize);
    std::cerr << "Classifying voxels ...\n";
    init_progress_bar(MULTIDIM_SIZE(SignificantT2()));
    int counter = 0, nxny = ysize * zsize;
#ifdef MODE7
#define NFEATURES 7
#endif
#ifdef MODE8
#define NFEATURES 8
#endif
#ifdef MODE15
#define NFEATURES 15
#endif

    FOR_ALL_ELEMENTS_IN_ARRAY3D(SignificantT2())
    {
        // Save the data for this voxel
        std::ofstream fh_dat;
        fh_dat.open("PPP.dat");
        if (!fh_dat)
            REPORT_ERROR(ERR_IO_NOTOPEN, "VariabilityClass::finishAnalysis: "
                         "Cannot open PPP.dat for output");
        fh_dat << NFEATURES << " " << nmax << std::endl;
        std::vector< Matrix1D<double> > v;
        v.clear();
        for (int n = 0; n < nmax; n++)
        {
            Matrix1D<double> v_aux(NFEATURES);

#ifdef MODE7

            v_aux(0) = VA[n](k,      i, j + xsize);
            v_aux(1) = VA[n](k, i + ysize,      j);
            v_aux(2) = VA[n](k, i + ysize, j + xsize);
            v_aux(3) = VA[n](k + zsize,      i,      j);
            v_aux(4) = VA[n](k + zsize,      i, j + xsize);
            v_aux(5) = VA[n](k + zsize, i + ysize,      j);
            v_aux(6) = VA[n](k + zsize, i + ysize, j + xsize);
#endif
#ifdef MODE8

            v_aux(0) = VA[n](k,      i,      j);
            v_aux(1) = VA[n](k,      i, j + xsize);
            v_aux(2) = VA[n](k, i + ysize,      j);
            v_aux(3) = VA[n](k, i + ysize, j + xsize);
            v_aux(4) = VA[n](k + zsize,      i,      j);
            v_aux(5) = VA[n](k + zsize,      i, j + xsize);
            v_aux(6) = VA[n](k + zsize, i + ysize,      j);
            v_aux(7) = VA[n](k + zsize, i + ysize, j + xsize);
#endif
#ifdef MODE15

            v_aux(0) = VA[n](k / 2, i / 2,    j / 2);
            v_aux(1) = VA[n](k / 2, i / 2, j / 2 + xsize2);
            v_aux(2) = VA[n](k / 2, i / 2 + ysize2,    j / 2);
            v_aux(3) = VA[n](k / 2, i / 2 + ysize2, j / 2 + xsize2);
            v_aux(4) = VA[n](k / 2 + zsize2,    i / 2,       j / 2);
            v_aux(5) = VA[n](k / 2 + zsize2,    i / 2, j / 2 + xsize2);
            v_aux(6) = VA[n](k / 2 + zsize2, i / 2 + ysize2,       j / 2);
            v_aux(7) = VA[n](k / 2 + zsize2, i / 2 + ysize2, j / 2 + xsize2);
            v_aux(8) = VA[n](k,      i, j + xsize);
            v_aux(9) = VA[n](k, i + ysize,      j);
            v_aux(10) = VA[n](k, i + ysize, j + xsize);
            v_aux(11) = VA[n](k + zsize,      i,      j);
            v_aux(12) = VA[n](k + zsize,      i, j + xsize);
            v_aux(13) = VA[n](k + zsize, i + ysize,      j);
            v_aux(14) = VA[n](k + zsize, i + ysize, j + xsize);
#endif

            v_aux = v_aux * v_aux;
            // COSS: Doesn't work: v_aux=v_aux.sort();

            fh_dat << v_aux.transpose() << std::endl;
            v.push_back(v_aux);
        }
        fh_dat.close();

        // Classify
        system("xmipp_fcmeans -din PPP.dat -std::cout PPP -c 2 -saveclusters > /dev/null");

        // Pick up results
        Matrix2D<double> aux;
        int n_previous;
#define GET_RESULTS(fh,fn,avg,cov,N,idx) \
    fh.open(fn);\
    if (!fh) \
        REPORT_ERROR(ERR_IO_NOTOPEN,(std::string)"VariabilityClass::finishAnalysis: " \
                     "Cannot open "+fn+" for input"); \
    n_previous=-1; \
    while (!fh.eof()) { \
        int n; fh >> n; \
        if (n!=n_previous) { \
            n_previous=n; N++; \
            idx(n)=1; \
            avg+=v[n]; \
            aux.fromVector(v[n]); \
            cov+=aux*aux.transpose(); \
        }\
    } \
    avg/=N; \
    cov/=N; \
    aux.fromVector(avg); \
    cov-=aux*aux.transpose();

        std::ifstream fh_0;
        Matrix1D<double> avg0(NFEATURES);
        Matrix1D<int>    idx0(nmax);
        Matrix2D<double> covariance0(NFEATURES, NFEATURES);
        int N0 = 0;
        GET_RESULTS(fh_0, "PPP.0", avg0, covariance0, N0, idx0);
#ifdef DEBUG

        std::cout << "Class 0 is:\n";
        for (int n = 0; n < idx0.size(); n++)
        {
            if (idx0(n))
            {
                int iact_proj = prm->ordered_list(n);
                std::cout << prm->IMG_Inf[iact_proj].fn_proj << std::endl;
            }
        }
#endif

        std::ifstream fh_1;
        Matrix1D<double> avg1(NFEATURES);
        Matrix1D<int>    idx1(nmax);
        Matrix2D<double> covariance1(NFEATURES, NFEATURES);
        int N1 = 0;
        GET_RESULTS(fh_1, "PPP.1", avg1, covariance1, N1, idx1);
#ifdef DEBUG

        std::cout << "Class 1 is:\n";
        for (int n = 0; n < idx1.size(); n++)
        {
            if (idx1(n))
            {
                int iact_proj = prm->ordered_list(n);
                std::cout << prm->IMG_Inf[iact_proj].fn_proj << std::endl;
            }
        }
#endif

        Matrix2D<double> T2, covariance;
        if (NFEATURES > 1)
        {
            // Perform T2-Hotelling test
            Matrix1D<double> avg_diff = avg1 - avg0;
            covariance = 1.0 / (N0 + N1 - 2) *
                         ((N0 - 1) * covariance0 + (N1 - 1) * covariance1);
            covariance *= (1.0 / N0 + 1.0 / N1);
            aux.fromVector(avg_diff);
            T2 = (double)(N0 + N1 - avg_diff.size() - 1) /
                 ((N0 + N1 - 2) * avg_diff.size()) *
                 aux.transpose() * covariance.inv() * aux;
        }
        else
        {
            // Perform t-test
            double variance = ((N0 - 1) * covariance0(0, 0) + (N1 - 1) * covariance1(0, 0)) /
                              (N0 + N1 - 2);
            double t = (avg1(0) - avg0(0)) / sqrt(variance * (1.0 / N0 + 1.0 / N1));
            T2.initZeros(1, 1);
            T2(0, 0) = t;
        }

        // Analysis of the covariance structure
        Matrix1D<double> eigenvalues;
        if (NFEATURES > 1)
        {
            Matrix2D<double> U, V;
            svdcmp(covariance, U, eigenvalues, V);
        }

        // Analysis of the coocurrences
        for (int n = 0; n < idx0.size(); n++)
            for (int np = n + 1; np < idx0.size(); np++)
                if (idx0(n) && idx0(np))
                    coocurrence(n, np)++;

        for (int n = 0; n < idx1.size(); n++)
            for (int np = n + 1; np < idx1.size(); np++)
                if (idx1(n) && idx1(np))
                    coocurrence(n, np)++;

        // Keep results
        SignificantT2(k, i, j) = T2(0, 0);
        if (NFEATURES > 1)
        {
            SignificantMinRatio(k, i, j) = eigenvalues(1) / eigenvalues(0);
            SignificantMaxRatio(k, i, j) = eigenvalues(NFEATURES - 1) / eigenvalues(0);
        }
#ifdef DEBUG
        std::cout << "T2 for this classification is " << T2(0, 0) << std::endl;
        std::cout << "Eigenvalues are " << eigenvalues.transpose() << std::endl;
#endif

        if (++counter % nxny == 0)
            progress_bar(counter);
    }
    progress_bar(MULTIDIM_SIZE(SignificantT2()));
    SignificantT2.write(prm->fn_root + "_variability.vol");
    SignificantMinRatio.write(prm->fn_root + "_minratio.vol");
    SignificantMaxRatio.write(prm->fn_root + "_maxratio.vol");
    system("rm PPP.dat PPP.cod PPP.vs PPP.err PPP.his PPP.inf PPP.0 PPP.1");

    for (int n = 0; n < nmax; n++)
        for (int np = n + 1; np < nmax; np++)
            coocurrence(np, n) = coocurrence(n, np);
    Image<double> save;
    typeCast(coocurrence, save());
    save.write(prm->fn_root + "_coocurrence.xmp");
}
#undef DEBUG

/* ------------------------------------------------------------------------- */
#define POCS_N_measure  8
#define POCS_N_use      2
/* Constructor ............................................................. */
POCSClass::POCSClass(Basic_ART_Parameters *_prm,
                     int _Zoutput_volume_size, int _Youtput_volume_size,
                     int _Xoutput_volume_size)
{
    prm = _prm;
    POCS_state = POCS_measuring;
    POCS_freq = prm->POCS_freq;
    POCS_i = 0;
    POCS_vec_i = 0;
    POCS_used = 0;
    POCS_N = 0;
    POCS_errors.initZeros(POCS_N_measure);
    Zoutput_volume_size = _Zoutput_volume_size;
    Youtput_volume_size = _Youtput_volume_size;
    Xoutput_volume_size = _Xoutput_volume_size;
    apply_POCS = (prm->surface_mask != NULL ||
                  prm->positivity || (prm->force_sym != 0 && !prm->is_crystal) ||
                  prm->known_volume != -1);
}

void POCSClass::newIteration()
{
    POCS_global_mean_error = 0;
}

void POCSClass::newProjection()
{
    POCS_N = 0;
}

/* Apply ................................................................... */
void POCSClass::apply(GridVolume &vol_basis, int it, int images)
{
    Image<double> vol_POCS, theo_POCS_vol, corr_POCS_vol, vol_voxels;

    if (apply_POCS && POCS_i % POCS_freq == 0)
    {
        Image<double> vol_aux;
        Image<double> *desired_volume = NULL;

        // Compute the corresponding voxel volume
        prm->basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                  Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
        if (prm->tell&TELL_SAVE_AT_EACH_STEP)
        {
            vol_voxels.write("PPPvolPOCS0.vol");
            std::cout << "Stats PPPvolPOCS0.vol: ";
            vol_voxels().printStats();
            std::cout << std::endl;
        }
        // Apply surface restriction
        if (prm->surface_mask != NULL)
        {
            vol_POCS() = (*(prm->surface_mask))();
        }
        else
        {
            vol_POCS().resize(vol_voxels());
            vol_POCS().initZeros();
        }
        if (prm->tell&TELL_SAVE_AT_EACH_STEP)
        {
            vol_POCS.write("PPPvolPOCS1.vol");
            std::cout << "Stats PPPvolPOCS1.vol: ";
            vol_POCS().printStats();
            std::cout << std::endl;
        }
        // Force symmetry
        if (prm->force_sym != 0)
        {
            symmetrizeVolume(prm->SL, vol_voxels(), vol_aux());
            desired_volume = &vol_aux;
            if (prm->tell&TELL_SAVE_AT_EACH_STEP)
            {
                vol_aux.write("PPPvolPOCS2.vol");
                std::cout << "Stats PPPvolPOCS2.vol: ";
                vol_aux().printStats();
                std::cout << std::endl;
            }
        }
        // Apply volume constraint
        if (prm->known_volume != -1)
        {
            Histogram1D hist;
            MultidimArray<int> aux_mask;
            aux_mask.resize(vol_POCS());
            FOR_ALL_ELEMENTS_IN_ARRAY3D(aux_mask)
            aux_mask(k, i, j) = 1 - (int)vol_POCS(k, i, j);
            long mask_voxels = vol_POCS().countThreshold("below", 0.5, 0);
            compute_hist_within_binary_mask(
                aux_mask, vol_voxels(), hist, 300);
            double known_percentage;
            known_percentage = XMIPP_MIN(100, 100 * prm->known_volume / mask_voxels);
            double threshold;
            threshold = hist.percentil(100 - known_percentage);
            FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_voxels())
            if (vol_voxels(k, i, j) < threshold)
                vol_POCS(k, i, j) = 1;
        }
        if (prm->tell&TELL_SAVE_AT_EACH_STEP)
        {
            vol_POCS.write("PPPvolPOCS3.vol");
            std::cout << "Stats PPPvolPOCS3.vol: ";
            vol_POCS().printStats();
            std::cout << std::endl;
        }

        // Do not allow positivity outside interest region
        // and do not allow negativity inside the interest region
        // if positivity restrictions are to be applied
        int bg = (int) vol_POCS().sum();
        int fg = MULTIDIM_SIZE(vol_POCS()) - bg;
        int relax = 0, posi = 0;
        FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_voxels())
        if (vol_POCS(k, i, j) == 1 && vol_voxels(k, i, j) < 0)
        {
            vol_POCS(k, i, j) = 0;
            relax++;
        }
        else if (vol_POCS(k, i, j) == 0 && vol_voxels(k, i, j) < 0 &&
                 prm->positivity)
        {
            vol_POCS(k, i, j) = 1;
            posi++;
        }
        // Debugging messages
        //std::cerr << "Relaxation/Positivity " << (double)relax/(double)bg << " "
        //     << (double)posi/(double)fg << " " << std::endl;

        // Solve volumetric equations
        switch (prm->basis.type)
        {
        case Basis::blobs:
            if (desired_volume == NULL)
                ART_voxels2blobs_single_step(vol_basis, &vol_basis,
                                             prm->basis.blob, prm->D, prm->lambda(it),
                                             &(theo_POCS_vol()), NULL,
                                             &(corr_POCS_vol()),
                                             &(vol_POCS()),
                                             POCS_mean_error, POCS_max_error, VARTK);
            else
            {
                FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_POCS())
                if (vol_POCS(k, i, j) == 1)
                    (*desired_volume)(k, i, j) = 0;
                for (int i = 0; i < prm->force_sym; i++)
                {
                    ART_voxels2blobs_single_step(vol_basis, &vol_basis,
                                                 prm->basis.blob, prm->D, prm->lambda(it),
                                                 &(theo_POCS_vol()), &((*desired_volume)()),
                                                 &(corr_POCS_vol()),
                                                 NULL,
                                                 POCS_mean_error, POCS_max_error, VARTK);
                    if (prm->tell&TELL_SAVE_AT_EACH_STEP)
                        std::cout << "    POCS Iteration " << i
                        << " POCS Error=" <<  POCS_mean_error << std::endl;
                }
            }
            break;
        case Basis::voxels:
            if (desired_volume == NULL)
            {
                FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_POCS())
                if (vol_POCS(k, i, j))
                    vol_basis(0)(k, i, j) = 0;
            }
            else
            {
                vol_basis(0)().initZeros();
                FOR_ALL_ELEMENTS_IN_ARRAY3D((*desired_volume)())
                vol_basis(0)(k, i, j) = (*desired_volume)(k, i, j);
            }
            POCS_mean_error = -1;
            break;
        }
        POCS_i = 1;
        POCS_global_mean_error += POCS_mean_error;
        POCS_N++;

        // Now some control logic
        if (prm->numIMG - images < 100 || images % 100 == 0 ||
            desired_volume != NULL)
        {
            POCS_freq = 1;
            POCS_state = POCS_measuring;
            POCS_vec_i = 0;
        }
        else
        {
            double dummy;
            switch (POCS_state)
            {
            case POCS_measuring:
#ifdef DEBUG_POCS

                std::cout << "M:" << POCS_vec_i << " " << POCS_mean_error << std::endl;
#endif

                POCS_errors(POCS_vec_i++) = POCS_mean_error;
                if (POCS_vec_i == POCS_N_measure)
                {
                    POCS_vec_i = 0;
                    // Change to use state
                    POCS_used = 0;
                    POCS_freq++;
                    POCS_state = POCS_use;
#ifdef DEBUG_POCS

                    std::cerr << "1: Changing to " << POCS_freq << std::endl;
#endif

                }
                break;
            case POCS_use:
                POCS_used++;
                POCS_errors.computeStats(POCS_avg,
                                         POCS_stddev, dummy, POCS_min);
#ifdef DEBUG_POCS

                std::cout << "Reference errors: " << POCS_errors.transpose() << std::endl;
                std::cout << "Checking " << ABS(POCS_mean_error - POCS_avg) << " " << 1.2*1.96*POCS_stddev << std::endl;
#endif

                if (ABS(POCS_mean_error - POCS_avg) < 1.2*1.96*POCS_stddev)
                {
                    if (POCS_mean_error < POCS_avg)
                    {
                        double max_error = POCS_errors(0);
                        POCS_vec_i = 0;
                        for (int i = 1; i < POCS_N_measure; i++)
                            if (POCS_errors(i) > max_error)
                            {
                                max_error = POCS_errors(i);
                                POCS_vec_i = i;
                            }
                        POCS_errors(POCS_vec_i) = POCS_mean_error;
                    }
                    if (POCS_used < POCS_N_use)
                    { // While not enough uses
                    }
                    else if (POCS_freq < 3)
                    { // increase frequency
                        POCS_freq++;
#ifdef DEBUG_POCS

                        std::cerr << "2: Changing to " << POCS_freq << std::endl;
#endif

                        POCS_used = 0;
                    }
                }
                else
                {
                    // It is behaving worse
                    if (POCS_freq > prm->POCS_freq + 1)
                    {
                        POCS_freq = prm->POCS_freq + 1;
                        POCS_used = 0;
#ifdef DEBUG_POCS

                        std::cerr << "3: Changing to " << POCS_freq << std::endl;
#endif

                    }
                    else if (POCS_used > 2)
                    {
                        POCS_freq = prm->POCS_freq;
                        // Change status
                        POCS_used = 0;
                        POCS_state = POCS_lowering;
#ifdef DEBUG_POCS

                        std::cerr << "Lowering\n";
#endif

                    }
                }
                break;
            case POCS_lowering:
                // Lower the POCS error before measuring again
                POCS_errors.computeStats(POCS_avg,
                                         POCS_stddev, POCS_max, dummy);
                POCS_used++;
                if (POCS_mean_error < POCS_max || POCS_used > 2*POCS_N_measure)
                {
                    // Change status
                    POCS_vec_i = 0;
                    POCS_state = POCS_measuring;
                }
                break;
            }
        }
    }
    else
    {
        POCS_i++;
        POCS_mean_error = -1;
    }
}


