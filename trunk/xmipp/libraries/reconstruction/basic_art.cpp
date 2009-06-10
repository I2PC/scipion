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


/* Here are all basic functions which are not extra_parameter dependent for
   the ART process. The extra parameter dependent functions are implemented
   in the Basic_art.inc file and must be included in each specific
   implementation (single particles, crystals, ...) */

#include "basic_art.h"
#include "recons_misc.h"

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
    shiftedTomograms   =      CHECK_PARAM(       "shiftedTomograms"     ); \
    goldmask    = textToFloat(GET_PARAM_WITH_DEF("goldmask","1e6")); \
    fn_ctf             =      GET_PARAM_WITH_DEF("CTF",     ""          ); \
    unmatched          =      CHECK_PARAM(       "unmatched"            );  \
    if (CHECK_PARAM("o")) \
        fn_root        =      GET_PARAM(         "o"                    );  \
    else fn_root       =      fn_sel.without_extension();                   \
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
    if (CHECK_PARAM("noisy_reconstruction")) { \
        if (parallel_mode!=ART) \
            REPORT_ERROR(1,"Basic_ART_Parameters::read: Noisy reconstructions" \
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
        if (R != -1.) R /= sampling;
        ref_trans_step /= sampling;
    }
    if (CHECK_PARAM("output_size"))
    {
        int i = paremeterPosition(argc, argv, "-output_size");
        if (i + 3 >= argc)
            REPORT_ERROR(1, "Not enough parameters after -output_size");
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
        REPORT_ERROR(3005,
                     (std::string)"Basic_ART_Parameters::read: There is a problem "
                     "opening the file " + fn);

    GET_ART_PARAMS;
    if (CHECK_PARAM("output_size"))
    {
        int argcp;
        char **argvp = NULL, *copy = NULL;
        generateCommandLine(fh, "output_size", argcp, argvp, copy);
        int i = paremeterPosition(argcp, argvp, "-output_size");
        if (i + 3 >= argcp)
            REPORT_ERROR(1, "Not enough parameters after -output_size");
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
    << "\n   [-goldmask value]       name of a sel file with masks"
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
    << "\n   [-dont_apply_shifts]  Do not apply shifts as stored in the 2D-image headers\n"
    << "\n   [-variability]        Perform variability analysis"
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
	<< "\n	 [-thr N=1]			   Number of threads to use. NOTE: Not available when using MPI."
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
                        Matrix1D<int> &ordered_list, int N)
{
    int   i, j, k;
    Matrix1D<short> chosen(numIMG);     // 1 if that image has been already
    // chosen
    double min_prod;
    int   min_prod_proj;
    Matrix2D<double> v(numIMG, 3);
    Matrix2D<double> euler;
    Matrix1D<double> product(numIMG);

    // Initialisation
    ordered_list.resize(numIMG);
    for (i = 0; i < numIMG; i++)
    {
        Matrix1D<double> z;
        // Initially no image is chosen
        VEC_ELEM(chosen, i) = 0;

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
    VEC_ELEM(chosen, i) = 1;
    VEC_ELEM(ordered_list, 0) = i;

    // Choose the rest of projections
    std::cerr << "Sorting projections ...\n";
    init_progress_bar(numIMG - 1);
    for (i = 1; i < numIMG; i++)
    {
        // Compute the product of not already chosen vectors with the just
        // chosen one, and select that which has minimum product
        min_prod = MAXFLOAT;
        for (j = 0; j < numIMG; j++)
            if (!VEC_ELEM(chosen, j))
            {
                VEC_ELEM(product, j) +=
                    ABS(dotProduct(v.Row(VEC_ELEM(ordered_list, i - 1)),
		    	v.Row(j)));
                if (N != -1 && i > N)
                    VEC_ELEM(product, j) -=
                        ABS(dotProduct(v.Row(VEC_ELEM(ordered_list, i - N - 1)),
			    v.Row(j)));
                if (VEC_ELEM(product, j) < min_prod)
                {
                    min_prod = VEC_ELEM(product, j);
                    min_prod_proj = j;
                }
            }

        // Store the chosen vector and mark it as chosen
        VEC_ELEM(ordered_list, i) = min_prod_proj;
        VEC_ELEM(chosen, min_prod_proj) = 1;

        // The progress bar is updated only every 10 images
        if (i % 10 == 0) progress_bar(i);
    }

    // A final call to progress bar to finish a possible small piece
    progress_bar(numIMG - 1);
    std::cout << std::endl;
}

/* ------------------------------------------------------------------------- */
/* No Sort                                                                   */
/* ------------------------------------------------------------------------- */
void no_sort(int numIMG, Matrix1D<int> &ordered_list)
{
    ordered_list.initLinear(0, numIMG - 1);
}

/* ------------------------------------------------------------------------- */
/* Random Sort                                                               */
/* ------------------------------------------------------------------------- */
void sort_randomly(int numIMG, Matrix1D<int> &ordered_list)
{
    int i;
    Matrix1D<int> chosen;

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
        VEC_ELEM(ordered_list, i - 1) = ptr;
        VEC_ELEM(chosen, ptr) = 1;

        // The progress bar is updated only every 10 images
        if (i % 10 == 0) progress_bar(i);
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
    SelFile     selfile, selctf;

    /* If checking the variability --------------------------------------------- */
    if (variability_analysis)
        parallel_mode = SIRT;

    /* Create history file handler --------------------------------------------- */
    if (level >= FULL)
    {
        fh_hist = new std::ofstream;
        fh_hist->open((fn_root + ".hist").c_str(), std::ios::out);
        if (!fh_hist)
            REPORT_ERROR(3008, (std::string)"Produce_Basic_ART_Side_Info: Cannot open file "
                         + fn_root + ".hist");
    }

    /* Get True Image number and projection size ------------------------------- */
    if (level >= BASIC)
    {
        //take into account weights here
        if (WLS)
        {
            SelFile SF_aux;
            headerXmipp       head;
            SF_aux.read(fn_sel);
            SF_aux.go_beginning();
            while (!SF_aux.eof())
            {
                head.read(SF_aux.get_current_file());
                if (head.Weight() != 0.)
                {
                    selfile.insert(SF_aux.current());
                }
                SF_aux.NextImg();
            }
            if (selfile.ImgNo() == 0)
            {
                std::cerr << "there is no input file with weight!=0" << std::endl;
                exit(1);
            }
        }
        else
            selfile.read(fn_sel);
        trueIMG = selfile.ImgNo();
        if (trueIMG == 0) REPORT_ERROR(3008, "Produce_Basic_ART_Side_Info: No images !!");
        selfile.ImgSize(projYdim, projXdim);
    }

    /* Get the CTF correction file ----------------------------------------- */
    if (fn_ctf != "")
    {
        selctf.read(fn_ctf);
        if (selctf.ImgNo() != selfile.ImgNo())
            REPORT_ERROR(1, "Basic_ART_Parameters: The number of images in "
                         "the ctf and original selfiles do not match");
    }

    /* Read symmetry file -------------------------------------------------- */
    if (level >= FULL)
    {
        double accuracy = (do_not_generate_subgroup) ? -1 : 1e-6;
        if (fn_sym != "") SL.read_sym_file(fn_sym, accuracy);
        if (!do_not_use_symproj) numIMG = trueIMG * (SL.SymsNo() + 1);
        else                     numIMG = trueIMG;
    }

    /* Read surface mask --------------------------------------------------- */
    if (level >= FULL)
    {
        if (fn_surface_mask != "")
        {
            surface_mask = new VolumeXmipp;
            surface_mask->read(fn_surface_mask);
            (*surface_mask)().setXmippOrigin();
        }
    }

    /* Fill ART_sort_info structure and Sort ------------------------------- */
    if (level >= FULL)
    {
        build_recons_info(selfile, selctf, fn_ctf, SL, IMG_Inf,
            do_not_use_symproj);

        if (!(tell&TELL_MANUAL_ORDER))
            if (parallel_mode == SIRT ||
                parallel_mode == pSIRT ||
                parallel_mode == pfSIRT ||
                parallel_mode == pCAV ||
                eq_mode == CAV ||
                rank > 0 || dont_sort)
                no_sort(numIMG, ordered_list);
            else if (random_sort)     sort_randomly(numIMG, ordered_list);
            else if (sort_last_N != -1) sort_perpendicular(numIMG, IMG_Inf, ordered_list,
                        sort_last_N);
            else                      no_sort(numIMG, ordered_list);
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
            {
                REPORT_ERROR(1, "negative weight not set correctly!");
                exit(1);
            }
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
            vol_basis0.read(fn_start);
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
    if (ray_length != -1) ray_length *= basis.max_length();
}
#undef DEBUG

/* Count number of equations for CAV --------------------------------------- */
void Basic_ART_Parameters::compute_CAV_weights(GridVolume &vol_basis0,
        int numProjs_node, int debug_level)
{
    if (GVNeq == NULL) GVNeq = new GridVolumeT<int>;
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
        read_proj.moveOriginTo_center();

        // Projection extension? .........................................
        if (proj_ext != 0)
            read_proj().window(
                STARTINGY(read_proj()) - proj_ext,
                STARTINGX(read_proj()) - proj_ext,
                FINISHINGY(read_proj()) + proj_ext,
                FINISHINGX(read_proj()) + proj_ext);

        count_eqs_in_projection(*GVNeq, basis, read_proj);

        if (debug_level > 0 &&
            act_proj % XMIPP_MAX(1, numIMG / 60) == 0) progress_bar(act_proj);
    }
    if (debug_level > 0)
    {
        progress_bar(numIMG);
        long int Neq = 0, Nunk = 0;
        for (int n = 0; n < GVNeq->VolumesNo(); n++)
            FOR_ALL_ELEMENTS_IN_MATRIX3D((*GVNeq)(n)())
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
