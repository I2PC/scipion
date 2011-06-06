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
#include "fourier_filter.h"
#include "recons_misc.h"

/* Default values ========================================================== */
void BasicARTParameters::default_values()
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
    kappa_list.resizeNoCopy(1);
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

void BasicARTParameters::defineParams(XmippProgram * program, const char* prefix, const char* comment)
{
    //  char tempLine[256];
    //
    //    if(prefix == NULL)
    //        sprintf(tempLine, "  [--basis <basis_type=blobs>] ");
    //    else
    //        sprintf(tempLine,"%s --basis <basis_type=blobs> ", prefix);
    //    if (comment != NULL)
    //        sprintf(tempLine, "%s : %s", tempLine, comment);
    //
    //    program->addParamsLine(tempLine);
    program->addParamsLine(" == I/O Parameters == ");
    program->addParamsLine("   -i <md_file>                : Metadata file with input projections");
    program->addParamsLine("   [--oroot <rootname>]        : Output rootname. If not supplied, input name is taken without extension.");
    program->addParamsLine("                               :+++The created files are as follows: %BR%");
    program->addParamsLine("                               :+++  =outputname.vol= 3D reconstruction in voxels %BR%");
    program->addParamsLine("                               :+++  =outputname.basis= 3D reconstruction in basis if the =--save_basis= option is enabled). The grid parameters are also stored in the same file %BR%");
    program->addParamsLine("                               :+++  =outputname.hist= History and information about the 3D reconstruction process %BR%");
    program->addParamsLine("   alias -o;");
    program->addParamsLine("   [--ctf <ctf_file=\"\">]     : Metadata file with CTFs");
    program->addParamsLine("   [--unmatched]               : apply unmatched forward/backward projectors");
    program->addParamsLine("   [--start <basisvolume_file=\"\">]  : Start from basisvolume. The reconstruction is performed in the same grid as the one ");
    program->addParamsLine("                               : in which the basis volume was stored (any -FCC or -CC or grid size value are useless)");
    program->addParamsLine("  [--max_tilt <alpha=10.e+6>]  : Skip projection with absolute tilt angle greater than alpha");
    program->addParamsLine("  [--ref_trans_after <n=-1>]   : Refine the translation alignment after n projections");
    program->addParamsLine("  [--ref_trans_step <n=-1>]    : Max displacement in translation alignment. This is a double");
    program->addParamsLine("  [--sparse <eps=-1>]          : Sparsity threshold");
    program->addParamsLine("  [--diffusion <eps=-1>]       : Diffusion weight");
    program->addParamsLine("  [--surface <surf_mask_file=\"\">] : Use this file as a surface mask");
    program->addParamsLine("  [--POCS_freq <f=1>]          : Impose POCS conditions every f projections");
    program->addParamsLine("  [--known_volume <vol=-1>]    : The volume is cut down to this mass, ie, the highest [mass] voxels are kept while ");
    program->addParamsLine("                               : the rest are set to 0");
    program->addParamsLine("  [--POCS_positivity]          : Apply positivity constraint");
    program->addParamsLine("  [--goldmask <value=1.e+6>]   : Pixels below this value are not considered for reconstruction");
    program->addParamsLine("  [--shiftedTomograms]         : Remove border pixels created by alignment of tomograms");
    program->addParamsLine("  [--dont_apply_shifts]        : Do not apply shifts as stored in the 2D-image headers");
    program->addParamsLine("  [--variability]              : Perform variability analysis");
    program->addParamsLine("  [--refine]                   : Refine input projection before backprojecting");
    program->addParamsLine("  [--noisy_reconstruction]     : Perform a companion noisy reconstruction. ");
    program->addParamsLine("                               :+ If given, the algorithm will perform two reconstructions. One with the input ");
    program->addParamsLine("                               :+ data, and another one with pure noise images applying exactly the same procedure ");
    program->addParamsLine("                               :+ as to the signal projections. This reconstruction is further used by ");
    program->addParamsLine("                               :+ the SSNR program in order to calculate the VSSNR or the SSNR itself.");
    program->addParamsLine("                               :+++ The created files are as follows: %BR%");
    program->addParamsLine("                               :+++  =[fn_root]_noise.vol= Reconstruction of the pure noise %BR%");
    program->addParamsLine("                               :+++  =[fn_root]_noise_proj.sel= Selection file with the pure noise images %BR%");
    program->addParamsLine("                               :+++  =[fn_root]_signal_proj.sel= Selection file with the signal images (a reordered version of the input (-i) selfile) %BR%");
    program->addParamsLine("                               :+++  =[fn_root]_noise_proj?????.xmp= Pure noise images used for the reconstruction %BR%");
    program->addParamsLine("  [--ray_length <r=-1>]        : Length of the ray in basis units that will be projected onto the image plane");

    program->addParamsLine(" == Symmetry parameters == ");
    program->addParamsLine("  [--sym <sym_file=\"\">]      : Use a symmetry file. It should give symmetry elements, ie, rotational axis, ");
    program->addParamsLine("                               : symmetry planes or whatever such that new points of view can be obtained");
    program->addParamsLine("  [--sym_each <n=0>]           : Force the reconstruction to be symmetric each n projections");
    program->addParamsLine("  [--force_sym <n=0>]          : Force the reconstruction to be symmetric n times at each projection");
    program->addParamsLine("  [--no_group]                 : Do not generate symmetry subgroup");
    program->addParamsLine("  [--no_symproj]               : Do not use symmetrized projections");

    program->addParamsLine(" == Iteration parameters == ");
    program->addParamsLine("  [-l <...>]                   : Relaxation factor, by default 0.01 (recommended range 0.0 - 0.1). ");
    program->addParamsLine("                               : A list of lambda values is also accepted as \"-l lambda0 lambda1 ...\"");
    program->addParamsLine("  [-n <noit=1>]                : Number of iterations");
    program->addParamsLine("  [--stop_at <it_stop=0>]      : Total number of iterated projections before algorithm stops. ");
    program->addParamsLine("                               :+ For instance, if there are 100 images, with two iterations and we ");
    program->addParamsLine("                               :+ want to stop at the half of the second iteration, then you must set it to 150");
    program->addParamsLine("  [--equation_mode <mode=ARTK> ]: Equation to project onto the hyperplane");
    program->addParamsLine("              where <mode> ");
    program->addParamsLine("        ARTK                   : Block ART");
    program->addParamsLine("        CAV                    : Component Averaging");
    program->addParamsLine("        CAVK                   : Block Component Averaging");
    program->addParamsLine("        CAVARTK                : Component Averaging Variant of Block ART");

    program->addParamsLine("  [--sort_last <N=2>]          : By default the algorithm sorts projections in the most orthogonally possible way. ");
    program->addParamsLine("                               : The most orthogonal way is defined as choosing the projection which maximizes the ");
    program->addParamsLine("                               : dot product with the N previous inserted projections. Use -1 to sort with all ");
    program->addParamsLine("                               : previous projections");
    program->addParamsLine(" or --random_sort              : Instead of orthogonal sort, projections are presented randomly to the algorithm");
    program->addParamsLine(" or --no_sort                  : No sort must be applied");
    program->addParamsLine("  [--WLS]                      : Perform weighted least squares ART");
    program->addParamsLine("  [-k <...> ]                  : Relaxation factor for WLS residual, by default 0.5. ");
    program->addParamsLine("                               : A list of kappa values is also accepted as \"-k kappa0 kappa1 ...\"");

    program->addParamsLine(" == Basis Parameters ==");
    Basis::defineParams(program);

    program->addParamsLine(" == Grid parameters == ");
    program->addParamsLine("   [-g <gridsz=1.41>]          : Relative size of the measuring unit in the grid lattice. By default, a unit in ");
    program->addParamsLine("                               : the grid system equals 1.41 units in the Universal System. This value is optimized ");
    program->addParamsLine("                               : for a BCC structure");
    program->addParamsLine("                               :+++ %BR% <verbatim>");
    program->addParamsLine("                               : if gridsz =  -1 => gridsz=2^(1/2)");
    program->addParamsLine("                               :+++ </verbatim> <verbatim> ");
    program->addParamsLine("                               :              -2 => gridsz=2^(1/3)");
    program->addParamsLine("                               :+++ </verbatim> ");
    program->addParamsLine("   [--grid_type <type=BCC>]    : Shape of the grid structure");
    program->addParamsLine("                where <type>");
    program->addParamsLine("            BCC                : Body Centered Cubic");
    program->addParamsLine("            FCC                : Face Centered  Cubic");
    program->addParamsLine("            SC                 : Simple Cubic");
    program->addParamsLine("   [-R <interest_sphere=-1>]   : Radius of the interest sphere. If provided, ART runs twice faster since only the ");
    program->addParamsLine("                               : sphere with this radius (in pixel units) is reconstructed");
    program->addParamsLine("   [--ext <proj_ext=0>]        : Projection extension. In order to avoid the box effect (those voxels near the volume ");
    program->addParamsLine("                               : borders are brighter than the rest), you can extent your projections resulting in ");
    program->addParamsLine("                               : a slower reconstruction but more accurate. Recommended values to avoid the box ");
    program->addParamsLine("                               : effect go from 5 to 10");
    program->addParamsLine("   [--output_size <Xsize=0> <Ysize=0> <Zsize=0>] : Output volume size in Pixels. Reconstruction size is taken from ");
    program->addParamsLine("                               : the projection size. However, the output size can be different, if the output volume ");
    program->addParamsLine("                               : is bigger, then the volume is zero padded.");
    program->addParamsLine("   [--sampling_rate <Ts=1>]    : Pixel size (Angstrom),  affects to -r, -g, -R and --ref_trans_step");

    program->addParamsLine(" ==+ Parallel parameters == ");
    program->addParamsLine(" : by default, sequential ART is applied");
    program->addParamsLine("   [--thr <N=1>]               : Number of threads to use. NOTE: Not available when using MPI.");
    program->addParamsLine("   [--parallel_mode <mode=ART>]: Parallel mode");
    program->addParamsLine("        where <mode>");
    program->addParamsLine("               ART             : Default");
    program->addParamsLine("               SIRT            : Simultaneous Iterative Reconstruction Technique");
    program->addParamsLine("   pSIRT                       : Parallel (MPI) Simultaneous Iterative Reconstruction Technique");
    program->addParamsLine("   pfSIRT                      : Parallel (MPI) False Simultaneous Iterative Reconstruction Technique (Faster convergence than pSIRT)");
    program->addParamsLine("   pSART                       : Parallel (MPI) Simultaneous ART");
    program->addParamsLine("   pAVSP                       : Parallel (MPI) Average Strings");
    program->addParamsLine("   pBiCAV                      : Parallel (MPI) Block Iterative CAV");
    program->addParamsLine("   pCAV                        : Parallel (MPI) CAV");
    program->addParamsLine("   [--block_size <n=1>]        : Number of projections to each block (SART and BiCAV)");

    program->addParamsLine("==+ Debugging options ==");
    program->addParamsLine("  [--print_system_matrix]      : Print the matrix of the system Ax=b. The format is:");
    program->addParamsLine("                               :+++ %BR% <verbatim>");
    program->addParamsLine("                               : Equation system (Ax=b) ---------------------- ");
    program->addParamsLine("                               :+++ </verbatim> <verbatim> ");
    program->addParamsLine("                               : pixel=<p> --> <b> = <a1> <a2> ... ");
    program->addParamsLine("                               :+++ </verbatim> ");
    program->addParamsLine("                               : I.e., for the pixel p (pixels are numbered lexicographically) with experimental ");
    program->addParamsLine("                               : value b, the equation ax=b is set. a is the corresponding row of matrix A. The ");
    program->addParamsLine("                               : coefficient a_i is equal to the contribution of the basis i to pixel p. x is the ");
    program->addParamsLine("                               : number of basis");
    program->addParamsLine("  [--show_iv <n=10>]           : Show volumes/images as the reconstruction goes. The volume is update every n projections");
    program->addParamsLine("  [--show_error]               : Show error for each projection");
    program->addParamsLine("  [--show_stats]               : Give some statistical information during the process,  they might be useful to see how the process is ");
    program->addParamsLine("                               : going. The mean squared error for each projection is shown by default");
    program->addParamsLine("  [--save_at_each_step]        : Save intermediate projections. This option allows deep debugging as it save all projections and volumes ");
    program->addParamsLine("                               : involved in the reconstruction process. After each step you are asked to press a key, so that you could ");
    program->addParamsLine("                               : inspect carefully the results. The names for these files are:");
    program->addParamsLine("                               : PPPtheo, PPPread, PPPcorr, PPPdiff");
    program->addParamsLine("                               : PPPbasis.basis, PPPvol.vol");
    program->addParamsLine("                               : PPPvolPOCS1, PPPvolPOCS2, PPPvolPOCS3");
    program->addParamsLine("  [--save_intermediate <n=0>]    : Save intermediate volumes (every <n> projections). If not provided, volumes are stored at each iteration ");
    program->addParamsLine("                               : and this parameter must be used at the end of the command to prevent errors. The names for these volumes are:");
    program->addParamsLine("                               :+++ %BR%");
    program->addParamsLine("                               : [filename root]it[it_no].vol Ex: art0001it0.vol ");
    program->addParamsLine("                               :+++ %BR%");
    program->addParamsLine("                               : [filename root]it]it_no].basis If the --save_basis option is enabled");
    program->addParamsLine("                               :+++ %BR%");
    program->addParamsLine("  [--save_basis]               : Save also the 3D reconstruction in basis each time that you have to save the reconstructed volume");
    program->addParamsLine("  [--manual_order]             : You are prompted to give the number of the following projection to be presented to the algorithm");
    program->addParamsLine("  [--only_sym]                 : Skip all those projections generated by symmetry (symmetries different from -1)");

}

void BasicARTParameters::readParams(XmippProgram * program)
{
    default_values();

    fn_sel = program->getParam("-i");

    if (program->checkParam("-o"))
        fn_root = program->getParam("-o");
    else
        fn_root = fn_sel.withoutExtension();

    fn_ctf = program->getParam("--ctf");
    unmatched = program->checkParam("--unmatched");
    fn_start = program->getParam("--start");
    max_tilt = program->getDoubleParam("--max_tilt");
    ref_trans_after = program->getIntParam("--ref_trans_after");
    ref_trans_step  = program->getIntParam("--ref_trans_step");
    sparseEps = program->getDoubleParam("--sparse");
    diffusionWeight = program->getDoubleParam("--diffusion");
    fn_surface_mask = program->getParam("--surface");
    POCS_freq = program->getIntParam("--POCS_freq");
    known_volume = program->getDoubleParam("--known_volume");
    positivity = program->checkParam("--POCS_positivity");
    goldmask = program->getDoubleParam("--goldmask");
    shiftedTomograms = program->checkParam("--shiftedTomograms");
    apply_shifts = !program->checkParam("--dont_apply_shifts");


    ray_length = program->getDoubleParam("--ray_length");

    // Symmetry parameters
    fn_sym = program->getParam("--sym");
    sym_each = program->getIntParam("--sym_each");
    force_sym = program->getIntParam("--force_sym");
    do_not_generate_subgroup = program->checkParam("--no_group");
    do_not_use_symproj = program->checkParam("--no_symproj");

    // Iteration parameters
    StringVector list;
    program->getListParam("-l", list);
    size_t listSize = list.size();

    if (listSize != 0)
    {
        lambda_list.resizeNoCopy(listSize);

        for (size_t k = 0; k < listSize; k++)
            VEC_ELEM(lambda_list, k) = textToFloat(list[k]);
    }

    no_it = program->getIntParam("-n");
    stop_at = program->getIntParam("--stop_at");

    String tempString = program->getParam("--equation_mode");
    if (tempString == "CAVK")
        eq_mode = CAVK;
    else if (tempString == "CAV")
        eq_mode = CAV;
    else if (tempString == "CAVARTK")
        eq_mode = CAVARTK;
    else
        eq_mode = ARTK;

    sort_last_N = program->getIntParam("--sort_last");
    random_sort = program->checkParam("--random_sort");
    dont_sort   = program->checkParam("--no_sort");
    WLS         = program->checkParam("--WLS");

    list.clear();
    program->getListParam("-k", list);
    listSize = list.size();
    if (listSize != 0)
    {
        kappa_list.resizeNoCopy(listSize);

        for (size_t k = 0; k < listSize; k++)
            VEC_ELEM(kappa_list, k) = textToFloat(list[k]);
    }

    // Basis parameters
    basis.readParams(program);

    // Grid parameters

    if (basis.type == Basis::voxels || basis.type == Basis::splines)
        grid_type = CC;

    if (program->checkParam("-g"))
    {
        grid_relative_size = program->getDoubleParam("-g");
        if (grid_relative_size == -1)
            grid_relative_size = sqrt (2.0);
        else if (grid_relative_size == -2)
            grid_relative_size = pow (2.0,1.0/3.0);
    }
    else
        grid_relative_size = basis.grid_relative_size;

    tempString = program->getParam("--grid_type");

    if (tempString == "BCC")
        grid_type = BCC;
    if (tempString == "FCC")
        grid_type = FCC;
    else if (tempString == "SC")
        grid_type = CC;
    else
        grid_type = BCC;

    R = program->getDoubleParam("-R");
    proj_ext = program->getIntParam("--ext");
    Xoutput_volume_size = program->getDoubleParam("--output_size", 0);
    Youtput_volume_size = program->getDoubleParam("--output_size", 1);
    Zoutput_volume_size = program->getDoubleParam("--output_size", 2);
    sampling = program->getDoubleParam("--sampling_rate");

    // Parallel parameters
    threads = program->getIntParam("--thr");

    tempString = program->getParam("--parallel_mode");

    if (tempString == "pSART")
        parallel_mode=pSART;
    else if (tempString == "pSIRT")
        parallel_mode=pSIRT;
    else if (tempString == "SIRT")
        parallel_mode=SIRT;
    else if (tempString == "pfSIRT")
        parallel_mode=pfSIRT;
    else if (tempString == "pBiCAV")
        parallel_mode=pBiCAV;
    else if (tempString == "pAVSP")
        parallel_mode=pAVSP;
    else if (tempString == "pCAV")
        parallel_mode=pCAV;
    else
        parallel_mode=ART;

    block_size = program->getIntParam("--block_size");
    //    fn_control = program->getParam("--control");


    // Debugging parameters

    print_system_matrix = program->checkParam("--print_system_matrix");
    if (program->checkParam("--show_error"))
        tell |= TELL_SHOW_ERROR;
    if (program->checkParam("--manual_order"))
        tell |= TELL_MANUAL_ORDER;
    if (program->checkParam("--only_sym"))
        tell |= TELL_ONLY_SYM;
    if (program->checkParam("--save_at_each_step"))
        tell |= TELL_SAVE_AT_EACH_STEP;
    if (program->checkParam("--save_basis"))
        tell |= TELL_SAVE_BASIS;
    if (program->checkParam("--show_stats"))
        tell |= TELL_STATS;
    if (program->checkParam("--save_intermediate"))
    {
        tell |= TELL_SAVE_INTERMIDIATE;
        save_intermidiate_every = program->getIntParam("save_intermediate");
    }
    if (program->checkParam("--show_iv"))
    {
        tell |= TELL_IV;
        save_intermidiate_every = program->getIntParam("save_intermediate");
    }



    if (program->checkParam("--variability"))
    {
        variability_analysis = true;
        parallel_mode = SIRT;
        no_it = 1;
    }

    refine = program->checkParam("--refine");

    if (program->checkParam("--noisy_reconstruction"))
    {
        if (parallel_mode != ART)
            REPORT_ERROR(ERR_ARG_INCORRECT,"BasicARTParameters::read: Noisy reconstructions" \
                         " can only be done for ART");
        else
            noisy_reconstruction = true;
    }


    //divide by the sampling rate
    if (sampling != 1.)
    {
        basis.set_sampling_rate(sampling);
        grid_relative_size /= sampling;
        if (R != -1.)
            R /= sampling;
        ref_trans_step /= sampling;
    }

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
            REPORT_ERROR(ERR_ARG_INCORRECT,"BasicARTParameters::read: Noisy reconstructions" \
                         " can only be done for ART"); \
        else noisy_reconstruction=true; \
    }

void BasicARTParameters::read(int argc, char **argv)
{
    //    GET_ART_PARAMS;
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
void BasicARTParameters::read(const FileName &fn)
{
    FILE *fh;
    if ((fh = fopen(fn.c_str(), "r")) == NULL)
        REPORT_ERROR(ERR_IO_NOTEXIST,fn);

    //    GET_ART_PARAMS;
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
void BasicARTParameters::usage()
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

void BasicARTParameters::usage_more()
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
    << "\n   [-save_at_each_step]  save intermediate projections"
    << "\n                             PPPtheo, PPPread, PPPcorr, PPPdiff"
    << "\n                             PPPbasis.basis, PPPvol.vol"
    << "\n                             PPPvolPOCS1, PPPvolPOCS2, PPPvolPOCS3"
    << "\n   [-save_intermediate <n>] save intermediate volumes (every <n> projections)"
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
void sortPerpendicular(int numIMG, ReconsInfo *IMG_Inf,
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

    // Initialization
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
void noSort(int numIMG, MultidimArray<int> &ordered_list)
{
    ordered_list.initLinear(0, numIMG - 1);
}

/* ------------------------------------------------------------------------- */
/* Random Sort                                                               */
/* ------------------------------------------------------------------------- */
void sortRandomly(int numIMG, MultidimArray<int> &ordered_list)
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
void BasicARTParameters::produceSideInfo(GridVolume &vol_basis0, int level,
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
        buildReconsInfo(selfile, fn_ctf, SL, IMG_Inf,
                          do_not_use_symproj);

        if (!(tell&TELL_MANUAL_ORDER))
            if (parallel_mode == SIRT ||
                parallel_mode == pSIRT ||
                parallel_mode == pfSIRT ||
                parallel_mode == pCAV ||
                eq_mode == CAV ||
                rank > 0 || dont_sort)
                noSort(numIMG, ordered_list);
            else if (random_sort)
                sortRandomly(numIMG, ordered_list);
            else if (sort_last_N != -1)
                sortPerpendicular(numIMG, IMG_Inf, ordered_list,
                                   sort_last_N);
            else
                noSort(numIMG, ordered_list);
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
void BasicARTParameters::compute_CAV_weights(GridVolume &vol_basis0,
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

int BasicARTParameters::ProjXdim()
{
    return projXdim;
}

int BasicARTParameters::ProjYdim()
{
    return projYdim;
}

/* ------------------------------------------------------------------------- */
/* Update residual vector for WLS                                            */
/* ------------------------------------------------------------------------- */
void updateResidualVector(BasicARTParameters &prm, GridVolume &vol_basis,
                          double &kappa, double &pow_residual_vol, double &pow_residual_imgs)
{
    GridVolume       residual_vol;
    Projection       read_proj, dummy_proj, new_proj;
    FileName         fn_resi, fn_tmp;
    double           sqrtweight, dim2, norma, normb, apply_kappa;
    ImageOver        *footprint = (ImageOver *) & prm.basis.blobprint;
    ImageOver        *footprint2 = (ImageOver *) & prm.basis.blobprint2;
    Matrix2D<double> *A = NULL;
    std::vector<MultidimArray<double> > newres_imgs;
    MultidimArray<int>    mask;

    residual_vol.resize(vol_basis);
    residual_vol.initZeros();

    // Calculate volume from all backprojected residual images
    std::cerr << "Backprojection of residual images " << std::endl;
    if (!(prm.tell&TELL_SHOW_ERROR))
        init_progress_bar(prm.numIMG);

    for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++)
    {
        // backprojection of the weighted residual image
        sqrtweight = sqrt(prm.residual_imgs[iact_proj].weight() / prm.sum_weight);
        read_proj = prm.residual_imgs[iact_proj];
        read_proj() *= sqrtweight;
        dummy_proj().resize(read_proj());

        dummy_proj.set_angles(prm.IMG_Inf[iact_proj].rot,
                              prm.IMG_Inf[iact_proj].tilt,
                              prm.IMG_Inf[iact_proj].psi);

        project_GridVolume(residual_vol, prm.basis, dummy_proj,
                           read_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                           prm.IMG_Inf[iact_proj].rot,
                           prm.IMG_Inf[iact_proj].tilt,
                           prm.IMG_Inf[iact_proj].psi, BACKWARD, prm.eq_mode,
                           prm.GVNeq, NULL, NULL, prm.ray_length, prm.threads);

        if (!(prm.tell&TELL_SHOW_ERROR))
            if (iact_proj % XMIPP_MAX(1, prm.numIMG / 60) == 0)
                progress_bar(iact_proj);
    }
    if (!(prm.tell&TELL_SHOW_ERROR))
        progress_bar(prm.numIMG);

    // Convert to voxels: solely for output of power of residual volume
    Image<double>      residual_vox;
    int Xoutput_volume_size = (prm.Xoutput_volume_size == 0) ?
                              prm.projXdim : prm.Xoutput_volume_size;
    int Youtput_volume_size = (prm.Youtput_volume_size == 0) ?
                              prm.projYdim : prm.Youtput_volume_size;
    int Zoutput_volume_size = (prm.Zoutput_volume_size == 0) ?
                              prm.projXdim : prm.Zoutput_volume_size;
    prm.basis.changeToVoxels(residual_vol, &(residual_vox()),
                             Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
    pow_residual_vol = residual_vox().sum2() / MULTIDIM_SIZE(residual_vox());
    residual_vox.clear();

    std::cerr << "Projection of residual volume; kappa = " << kappa << std::endl;
    if (!(prm.tell&TELL_SHOW_ERROR))
        init_progress_bar(prm.numIMG);

    // Now that we have the residual volume: project in all directions
    pow_residual_imgs = 0.;
    new_proj().resize(read_proj());
    mask.resize(read_proj());
    BinaryCircularMask(mask, YSIZE(read_proj()) / 2, INNER_MASK);

    dim2 = (double)YSIZE(read_proj()) * XSIZE(read_proj());
    for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++)
    {
        project_GridVolume(residual_vol, prm.basis, new_proj,
                           dummy_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                           prm.IMG_Inf[iact_proj].rot,
                           prm.IMG_Inf[iact_proj].tilt,
                           prm.IMG_Inf[iact_proj].psi, FORWARD, prm.eq_mode,
                           prm.GVNeq, A, NULL, prm.ray_length, prm.threads);

        sqrtweight = sqrt(prm.residual_imgs[iact_proj].weight() / prm.sum_weight);

        // Next lines like normalization in [EHL] (2.18)?
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(new_proj())
        {
            dAij(dummy_proj(), i, j) = XMIPP_MAX(1., dAij(dummy_proj(), i, j)); // to avoid division by zero
            dAij(new_proj(), i, j) /= dAij(dummy_proj(), i, j);
        }
        new_proj() *= sqrtweight * kappa;

        /*
        fn_tmp="residual_"+integerToString(iact_proj);
        dummy_proj()=1000*prm.residual_imgs[iact_proj]();
        dummy_proj.write(fn_tmp+".old");
        */

        prm.residual_imgs[iact_proj]() -= new_proj();
        pow_residual_imgs += prm.residual_imgs[iact_proj]().sum2();

        // Mask out edges of the images
        apply_binary_mask(mask, prm.residual_imgs[iact_proj](), prm.residual_imgs[iact_proj](), 0.);

        /*
        dummy_proj()=1000*new_proj();
        dummy_proj.write(fn_tmp+".change");
        dummy_proj()=1000*prm.residual_imgs[iact_proj]();
        dummy_proj.write(fn_tmp+".new");
        */

        if (!(prm.tell&TELL_SHOW_ERROR))
            if (iact_proj % XMIPP_MAX(1, prm.numIMG / 60) == 0)
                progress_bar(iact_proj);
    }

    pow_residual_imgs /= dim2;
    newres_imgs.clear();

    if (!(prm.tell&TELL_SHOW_ERROR))
        progress_bar(prm.numIMG);

}


