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

/* Desctructor */
BasicARTParameters::~BasicARTParameters()
{
    delete fh_hist;
    delete[] IMG_Inf;
    delete D;
    delete Dinv;
    delete GVNeq;
    delete surface_mask;
}

/* Default values ========================================================== */
void BasicARTParameters::defaultValues()
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
    basis.setDefault();
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

void BasicARTParameters::defineParams(XmippProgram * program, bool mpiMode)
{

    program->addParamsLine(" == I/O Parameters == ");
    program->addParamsLine("   -i <md_file>                : Metadata file with input projections");
    program->addParamsLine("   [-o <volume_file=\"rec_art.vol\">]  : Filename for output volume.");
    program->addParamsLine("             : Rootname for rest of output files is taken from volume filename");
    program->addParamsLine("                               :+++The created files are as follows: %BR%");
    program->addParamsLine("                               :+++  =outputname.vol= 3D reconstruction in voxels %BR%");
    program->addParamsLine("                               :+++  =outputname.basis= 3D reconstruction in basis if the =--save_basis= option is enabled). The grid parameters are also stored in the same file %BR%");
    program->addParamsLine("                               :+++  =outputname.hist= History and information about the 3D reconstruction process %BR%");
    program->addParamsLine("   [--ctf <ctf_file=\"\">]     : Metadata file with CTFs");
    program->addParamsLine("   [--unmatched]               : Apply unmatched forward/backward projectors");
    program->addParamsLine("   [--start <basisvolume_file=\"\">]  : Start from this basis volume. The reconstruction is performed in the same grid as the one ");
    program->addParamsLine("                               : in which the basis volume was stored (any -FCC or -CC or grid size value are useless)");
    program->addParamsLine("  [--max_tilt <alpha=10.e+6>]  : Skip projections with absolute tilt angle greater than alpha.");
    program->addParamsLine("                               : It means that if alpha=40, then only images with tilt angle ");
    program->addParamsLine("                               : within the ranges 0+/-40 and 180+/-40 will be processed. (Projection");
    program->addParamsLine("                               : tilt angles are forced to be in the range 0-360)");
    program->addParamsLine("  [--ref_trans_after <n=-1>]   : Refine the translation alignment after n projections. (Integer type)");
    program->addParamsLine("  [--ref_trans_step <v=-1>]    : Maximum displacement in translation alignment. (Double type)");
    program->addParamsLine("  [--sparse <eps=-1>]          : Sparsity threshold");
    program->addParamsLine("  [--diffusion <eps=-1>]       : Diffusion weight");
    program->addParamsLine("  [--surface <surf_mask_file=\"\">] : Mask for the surface constraint. It says where the volume is known to be 0");
    program->addParamsLine("  [--POCS_freq <f=1>]          : Impose POCS conditions every f projections");
    program->addParamsLine("  [--known_volume <value=-1>]  : The volume is cut down to this mass, ie, the highest [value] voxels are kept while ");
    program->addParamsLine("                               : the rest are set to 0");
    program->addParamsLine("  [--POCS_positivity]          : Force the resulting volume to be positive");
    program->addParamsLine("  [--goldmask <value=1.e+6>]   : Pixels below this value are considered to come frome gold beads and are not used for reconstruction");
    program->addParamsLine("  [--shiftedTomograms]         : Remove external zero-valued border pixels created by alignment of tomograms");
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
    program->addParamsLine("                               :+++  =[fn_root]_noise_proj.stk= Pure noise images used for the reconstruction %BR%");
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
    program->addParamsLine("   [-g <gridsz=1.41>]          : Relative size of the measuring unit in the grid lattice in pixels. By default, a unit in ");
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
    program->addParamsLine("   [--sampling_rate <Ts=1>]    : Pixel size (Angstroms)");

    program->addParamsLine(" ==+ Parallel parameters == ");
    program->addParamsLine(" : by default, sequential ART is applied");
    program->addParamsLine("   [--thr <N=1>]               : Number of threads to use. NOTE: Not available when using MPI.");
    program->addParamsLine("   [--parallel_mode <mode=ART>]: Parallelization algorithm to use with threads or MPI program version");
    program->addParamsLine("        where <mode>");
    program->addParamsLine("               ART             : Default");
    program->addParamsLine("               SIRT            : Simultaneous Iterative Reconstruction Technique");
    program->addParamsLine("   pSIRT                       : Parallel (MPI) Simultaneous Iterative Reconstruction Technique");
    program->addParamsLine("   pfSIRT                      : Parallel (MPI) False Simultaneous Iterative Reconstruction Technique (Faster convergence than pSIRT)");
    program->addParamsLine("   pSART                       : Parallel (MPI) Simultaneous ART");
    program->addParamsLine("   pAVSP                       : Parallel (MPI) Average Strings");
    program->addParamsLine("   pBiCAV                      : Parallel (MPI) Block Iterative CAV");
    program->addParamsLine("   pCAV                        : Parallel (MPI) CAV");
    program->addParamsLine("   [--block_size <n=1>]        : Number of projections for each block (SART and BiCAV)");

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
    defaultValues();

    fn_sel = program->getParam("-i");
    fn_out = program->getParam("-o");
    fn_root = fn_out.withoutExtension();

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
    Xoutput_volume_size = program->getIntParam("--output_size", 0);
    Youtput_volume_size = program->getIntParam("--output_size", 1);
    Zoutput_volume_size = program->getIntParam("--output_size", 2);
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

    verbose = program->getIntParam("--verbose");

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

    // Measures are given in pixels, independent of pixel size
    //    //divide by the sampling rate
    //    if (sampling != 1.)
    //    {
    //        basis.setSamplingRate(sampling);
    //        grid_relative_size /= sampling;
    //        if (R != -1.)
    //            R /= sampling;
    //        ref_trans_step /= sampling;
    //    }


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
            SF_aux.read(fn_sel);
            if (SF_aux.containsLabel(MDL_ENABLED))
                SF_aux.removeObjects(MDValueEQ(MDL_ENABLED, -1));
            selfile.clear();
            selfile.importObjects(SF_aux, MDValueRange(MDL_WEIGHT, 1e-9, 99e99));
            if (selfile.size() == 0)
                REPORT_ERROR(ERR_MD_OBJECTNUMBER, "There is no input file with weight!=0");
        }
        else
        {
            selfile.read(fn_sel);
            if (selfile.containsLabel(MDL_ENABLED))
                selfile.removeObjects(MDValueEQ(MDL_ENABLED, -1));
        }
        trueIMG = selfile.size();
        if (trueIMG == 0)
            REPORT_ERROR(ERR_MD_OBJECTNUMBER, "Produce_Basic_ART_Side_Info: No images !!");
        size_t idum, idumLong;
        getImageSize(selfile, projXdim, projYdim, idum, idumLong);
    }

    /* Read symmetry file -------------------------------------------------- */
    if (level >= FULL)
    {
        double accuracy = (do_not_generate_subgroup) ? -1 : 1e-6;
        if (fn_sym != "")
            SL.readSymmetryFile(fn_sym, accuracy);
        if (!do_not_use_symproj)
            numIMG = trueIMG * (SL.symsNo() + 1);
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
        {
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
        if ( !fn_start.empty() )
        {
            if (fn_start.contains("basis")) // A basis file
                vol_basis0.read(fn_start, basis.basisName());
            else // If it is a volume of voxels
            {
                Image<double> imTemp;
                imTemp.read(fn_start);
                basis.changeFromVoxels(imTemp(), vol_basis0, grid_type, grid_relative_size,
                                       NULL, NULL, R, threads);
            }
        }
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
                /* If you subtract half the basis radius, you are forcing that the
                last basis touches slightly the volume border. By not subtracting
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
        basis.setD(D);
        basis.produceSideInfo(vol_basis0.grid());
    }

    /* Express the ray length in basis units ----------------------------------- */
    if (ray_length != -1)
        ray_length *= basis.maxLength();

    /* With CAV equalization mode weights must be calculated, but for the parallel cases
       where weights are calculated in a parallel manner.*/
    if(eq_mode == CAV && parallel_mode != pCAV && parallel_mode != pBiCAV)
        computeCAVWeights(vol_basis0, numIMG, verbose-1);
}
#undef DEBUG

/* Count number of equations for CAV --------------------------------------- */
void BasicARTParameters::computeCAVWeights(GridVolume &vol_basis0,
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
        ReconsInfo &imgInfo = IMG_Inf[ordered_list(act_proj)];

        read_proj.read(imgInfo.fn_proj, apply_shifts, DATA, &imgInfo.row);
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
        for (size_t n = 0; n < GVNeq->VolumesNo(); n++)
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




