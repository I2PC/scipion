/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2007)
 *
 * Unidad de Bioinformatica del Centro Nacional de Biotecnologia , CSIC
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
#include "ml_tomo.h"
#include <data/metadata_extension.h>

//#define JM_DEBUG

// For blocking of threads
pthread_mutex_t mltomo_weightedsum_update_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mltomo_selfile_access_mutex = PTHREAD_MUTEX_INITIALIZER;

// Usage ===================================================================
void
ProgMLTomo::defineParams()
{
    addUsageLine(
        "Align and classify 3D images with missing data regions in Fourier space,");
    addUsageLine(
        "e.g. subtomograms or RCT reconstructions, by a 3D multi-reference refinement");
    addUsageLine("based on a maximum-likelihood (ML) target function.");
    addUsageLine(
        "+++For several cases, this method has been shown to be able to both align and "
        "classify in a completely *reference-free* manner, by starting from random assignments of "
        "the orientations and classes. The mathematical details behind this approach are explained "
        "in detail in");
    addUsageLine("+++Scheres et al. (2009) Structure, 17, 1563-1572", true);
    addUsageLine(
        "+++*Please cite this paper if this program is of use to you!* %BR% "
        "There also exists a standardized python script [[http://newxmipp.svn.sourceforge.net/viewvc/"
        "newxmipp/trunk/xmipp/applications/scripts/protocols/protocol__mltomo.py]"
        "[xmipp_protocol__mltomo.py]] for this program. "
        "Thereby, rather than executing the command line options explained below, the user"
        " can submit his jobs through a convenient GUI in the [[GettingStartedWithProtocols][Xmipp_protocols]], "
        "although we still recommend reading this page carefully in order "
        "to fully understand the options given in the protocol. Note that this protocol is available "
        "from the main xmipp_protocols setup window by pressing the _Additional protocols_ button.)");

    addParamsLine(
        "   -i <metadata>                : Metadata file with input images (and angles) ");
    addParamsLine(
        "   --nref <int=0>               : Number of references to generate automatically (recommended)");
    addParamsLine(
        "   OR --ref <file=\"\">         : or metadata file with initial references/single reference image ");
    addParamsLine(" [ --oroot <rootname=mltomo> ]  : Output rootname ");
    addParamsLine("        alias -o;");
    //docfile not implemented, it can go in the -i option
    //addParamsLine(" [ --doc <metadata=\"\"> ]      : MetaData with orientations  for all particles ");
    addParamsLine(" [ --sym <symgroup=c1> ]        : Symmetry group ");
    addParamsLine(
        " [ --missing <metadata=\"\"> ]  : Metadata file with missing data region definitions");

    addParamsLine("== Angular sampling ==");
    addParamsLine(
        " [ --ang <float=10.> ]           : Angular sampling rate (in degrees)");
    addParamsLine(
        " [ --ang_search <float=-1.> ]    : Angular search range around orientations from Metadata ");
    addParamsLine(
        "                                 : (by default, exhaustive searches are performed)");
    addParamsLine(
        " [ --dont_limit_psirange ]       : Exhaustive psi searches when using -ang_search (only for c1 symmetry)");
    addParamsLine(
        " [ --limit_trans <float=-1.> ]   : Maximum allowed shifts (negative value means no restriction)");
    addParamsLine(
        " [ --tilt0+ <float=-91.> ]       : Limit tilt angle search from tilt0 to tiltF (in degrees) ");
    addParamsLine(
        " [ --tiltF+ <float=91.> ]        : Limit tilt angle search from tilt0 to tiltF (in degrees) ");
    addParamsLine(
        " [ --psi_sampling+ <float=-1.> ] : Angular sampling rate for the in-plane rotations(in degrees)");

    addParamsLine("== Regularization ==");
    addParamsLine(
        " [ --reg0 <float=0.> ]          : Initial regularization parameters (in N/K^2) ");
    addParamsLine(
        " [ --regF <float=0.> ]          : Final regularization parameters (in N/K^2) ");
    addParamsLine(
        " [ --reg_steps <int=5> ]        : Number of iterations in which the regularization is changed from reg0 to regF");

    addParamsLine("== Others ==");
    addParamsLine(
        " [ --dont_rotate ]              : Keep orientations from Metadata fixed, only translate and classify ");
    addParamsLine(
        " [ --dont_align ]               : Keep orientations and tran Metadata (otherwise start from random)");
    addParamsLine(
        " [ --only_average ]             : Keep orientations and classes from docfile, only output weighted averages ");
    addParamsLine(
        " [ --perturb ]                  : Apply random perturbations to angular sampling in each iteration");
    addParamsLine(
        " [ --dim <int=-1> ]             : Use downscaled (in fourier space) images of this size ");
    addParamsLine(
        " [ --maxres <float=0.5> ]       : Maximum resolution (in pixel^-1) to use ");
    addParamsLine(
        " [ --thr <int=1> ]              : Number of shared-memory threads to use in parallel ");

    addParamsLine("==+ Additional options: ==");
    addParamsLine(
        " [ --impute_iter <int=1> ]      : Number of iterations for inner imputation loop ");
    addParamsLine(
        " [ --iter <int=25> ]           : Maximum number of iterations to perform ");
    addParamsLine(
        " [ --istart <int=1> ]             : number of initial iteration ");
    addParamsLine(
        " [ --noise <float=1.> ]          : Expected standard deviation for pixel noise ");
    addParamsLine(
        " [ --offset <float=3.> ]         : Expected standard deviation for origin offset [pix]");
    addParamsLine(
        " [ --frac <metadata=\"\"> ]     : Metadata with expected model fractions (default: even distr.)");
    addParamsLine(
        " [ --restart <logfile> ]        : restart a run with all parameters as in the logfile ");
    addParamsLine(
        " [ --fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise ");
    addParamsLine(
        " [ --fix_sigma_offset]          : Do not re-estimate the standard deviation in the origin offsets ");
    addParamsLine(
        " [ --fix_fractions]             : Do not re-estimate the model fractions ");
    addParamsLine(" [ --eps <float=5e-5> ]         : Stopping criterium ");
    addParamsLine(
        " [ --pixel_size <float=1> ]     : Pixel size (in Anstrom) for resolution in FSC plots ");

    addParamsLine(
        " [ --mask <maskfile=\"\"> ]          : Mask particles; only valid in combination with -dont_align ");
    addParamsLine(
        " [ --maxCC ]                    : Use constrained cross-correlation and weighted averaging instead of ML ");
    addParamsLine(
        " [ --dont_impute ]              : Use weighted averaging, rather than imputation ");
    addParamsLine(
        " [ --noimp_threshold <float=1.>] : Threshold to avoid division by zero for weighted averaging ");

    addParamsLine("==+++++ Hidden arguments ==");
    addParamsLine(" [--trymindiff_factor <float=0.9>]");
    addParamsLine(" [ --no_SMALLANGLE <float=5e-5> ]: ");
    addParamsLine(" [ --keep_angles]: ");
    addParamsLine(" [ --filter]: ");
    addParamsLine(" [ --ini_filter]: ");

    addExampleLine("+++---+++Input images file", false);
    addExampleLine(
        "+++ The input metadata should contain the =image= column = and =missingRegionNumber=, "
        "indicating the subtomogran filename and the missing region number, respectively. It can"
        "also contains columns with angles and shift information. The output will be a metadata"
        "with the same format. Follow is an example:", false);
    addExampleLine("+++ # XMIPP_STAR_1 *", true);
    addExampleLine("+++ #             ", true);
    addExampleLine("+++ data_          ", true);
    addExampleLine("+++ loop_          ", true);
    addExampleLine("+++  _image        ", true);
    addExampleLine("+++  _missingRegionNumber  ", true);
    addExampleLine("+++  _angleRot     ", true);
    addExampleLine("+++  _angleTilt    ", true);
    addExampleLine("+++  _anglePsi     ", true);
    addExampleLine("+++  _shiftX       ", true);
    addExampleLine("+++  _shiftY       ", true);
    addExampleLine("+++  _shiftZ       ", true);
    addExampleLine("+++  _ref          ", true);
    addExampleLine("+++  _logLikelihood        ", true);
    addExampleLine(
        "+++   32_000001.scl  0.00 0.000000 0.000  0.000  0.000000     0.000000  0    0.000000 1",
        true);
    addExampleLine(
        "+++   32_000002.scl  0.00 0.000000 0.000  0.000  0.000000     0.000000  0    0.000000 1",
        true);
    //  where:
    //  =ref= is the number of the corresponding reference (from 1 to NUMBER_OF_REFERENCES)
    //  =missingRegionNumber= is the number of the corresponding missing region number (starting at 1) from the missing regionn Metadata.
    //  =Pmax= is the height of the maximum of the (normalized) probability distribution. (Pmax=1 means delta functions, small Pmax values mean broad distributions)
    //  =dLL= is the contribution to the log-likelihood function for this image. The larger this value the better the image fits to the model.
    //
    //
    //  addExampleLine("+++ ---+++MissingDocFileFormat",false);
    //  The format of the docfile to define missing data regions is as follows:
    addExampleLine("+++---+++Input missing regions file", false);
    addExampleLine("+++ # XMIPP_STAR_1 *", true);
    addExampleLine("+++ # Wedgeinfo", true);
    addExampleLine("+++ data_", true);
    addExampleLine("+++ loop_          ", true);
    addExampleLine("+++  _missingRegionNumber ", true);
    addExampleLine("+++  _missingRegionType   ", true);
    addExampleLine("+++  _missingRegionThetaY0", true);
    addExampleLine("+++  _missingRegionThetaYF", true);
    addExampleLine("+++   1  wedge_y  -64  64", true);

    addExampleLine(
        "+++ The first column =missingRegionNumber= (starting at 1) is required for each type "
        "of missing region, this number should appears in the input images metadata %BR% ",
        false);
    addExampleLine(
        "+++ Here =missingRegionType= can be one of the following: %BR% ", false);
    addExampleLine(
        "+++   * =wedge_y= for a missing wedge where the tilt axis is along Y,",
        false);
    addExampleLine(
        "+++    colums =missingRegionThetaY0= and =missingRegionThetaYF= are used",
        false);
    addExampleLine(
        "+++   * =wedge_x= for a missing wedge where the tilt axis is along X,",
        false);
    addExampleLine(
        "+++    colums =missingRegionThetaX0= and =missingRegionThetaXF= are used",
        false);
    addExampleLine(
        "+++   * =pyramid= for a missing pyramid where the tilt axes are along Y and X,",
        false);
    addExampleLine("+++    same columns as =wedge_y= and =wedge_x= are used",
                   false);
    addExampleLine(
        "+++   * =cone= for a missing cone (pointing along Z)                    ",
        false);
    addExampleLine("+++    column =missingRegionThetaY0= is used", false);

}

// Read arguments ==========================================================
void
ProgMLTomo::readParams()
{
    // Generate new command line for restart procedure
    cline = "";
    int argc2 = 0;
    char ** argv2 = NULL;

    if (checkParameter(argc, argv, "-restart"))
    {
        REPORT_ERROR(ERR_DEBUG_TEST,
                     "restart procedure temporarily de-activated.... sorry!");
        /*
         std::string comment;
         FileName fn_sel;
         MetaData DFi;
         DFi.read(getParameter(argc, argv, "-restart"));
         DFi.go_beginning();
         comment = (DFi.get_current_line()).get_text();
         if (strstr(comment.c_str(), "ml_tomo-logfile") == NULL)
         {
         std::cerr << "Error!! MetaData is not of ml_tomo-logfile type. " << std::endl;
         exit(1);
         }
         else
         {
         char *copy;
         int n = 0;
         int nmax = DFi.dataLineNo();
         SFr.reserve(nmax);
         copy = NULL;
         DFi.next();
         comment = " -frac " + DFi.name();
         fn_sel = DFi.name();
         fn_sel = fn_sel.without_extension() + "_restart.sel";
         comment += " -ref " + fn_sel;
         comment += (DFi.get_current_line()).get_text();
         DFi.next();
         cline = (DFi.get_current_line()).get_text();
         // Also read additional options upon restart
         for (int i=1; i <argc; i+=2)
         {
         std::cerr<<"i= "<<i<<" argc= "<<argc<<" argv[i]"<<argv[i]<<std::endl;
         if ((strcmp("-restart", argv[i]) != 0))
         {
         comment += " ";
         comment += argv[i];
         comment += " ";
         comment += argv[i+1];
         comment += " ";
         }
         }
         comment = comment + cline;
         // regenerate command line
         generateCommandLine(comment, argc2, argv2, copy);
         // Read images names from restart file
         DFi.next();
         while (n < nmax)
         {
         n++;
         DFi.next();
         if (DFi.get_current_line().Is_comment()) fn_sel = ((DFi.get_current_line()).get_text()).erase(0, 3);
         SFr.insert(fn_sel, SelLine::ACTIVE);
         DFi.adjust_to_data_line();
         }
         fn_sel = DFi.name();
         fn_sel = fn_sel.without_extension() + "_restart.sel";
         SFr.write(fn_sel);
         SFr.clear();
         }
         */
    }
    else
    {
        // no restart, just copy argc to argc2 and argv to argv2
        argc2 = argc;
        argv2 = (char **)argv;
        for (int i = 1; i < argc2; i++)
        {
            cline = cline + (std::string) argv2[i] + " ";
        }
    }

    fn_sel = getParam("-i");
    nr_ref = getIntParam("--nref");
    fn_ref = getParam("--ref");
    if (!fn_ref.empty() && (!compareImageSize(fn_sel, fn_ref)))
        REPORT_ERROR(ERR_GRID_SIZE, "Reference and images aren't of same size");
    //fn_doc = getParam("--doc");

    fn_root = getParam("--oroot");
    fn_sym = getParam("--sym");
    Niter = getIntParam("--iter");
    Niter2 = getIntParam("--impute_iter");
    istart = getIntParam("--istart");
    sigma_noise = getDoubleParam("--noise");
    sigma_offset = getDoubleParam("--offset");
    fn_frac = getParam("--frac");
    fix_fractions = checkParam("--fix_fractions");
    fix_sigma_offset = checkParam("--fix_sigma_offset");
    fix_sigma_noise = checkParam("--fix_sigma_noise");
    eps = getDoubleParam("--eps");
    no_SMALLANGLE = checkParam("--no_SMALLANGLE");
    do_keep_angles = checkParam("--keep_angles");
    do_perturb = checkParam("--perturb");
    pixel_size = getDoubleParam("--pixel_size");

    // Low-pass filter
    do_filter = checkParam("--filter");
    do_ini_filter = checkParam("--ini_filter");

    // regularization
    reg0 = getDoubleParam("--reg0");
    regF = getDoubleParam("--regF");
    reg_steps = getIntParam("--reg_steps");

    // ML (with/without) imputation, or maxCC
    do_ml = !checkParam("--maxCC");
    do_impute = !checkParam("--dont_impute");
    noimp_threshold = getDoubleParam("--noimp_threshold");

    // Angular sampling
    angular_sampling = getDoubleParam("--ang");
    psi_sampling = getDoubleParam("--psi_sampling");
    tilt_range0 = getDoubleParam("--tilt0");
    tilt_rangeF = getDoubleParam("--tiltF");
    ang_search = getDoubleParam("--ang_search");
    do_limit_psirange = !checkParam("--dont_limit_psirange");
    limit_trans = getDoubleParam("--limit_trans");

    // Skip rotation, only translate and classify
    dont_rotate = checkParam("--dont_rotate");
    // Skip rotation and translation, only classify
    dont_align = checkParam("--dont_align");
    // Skip rotation and translation and classification
    do_only_average = checkParam("--only_average");

    // For focussed classification (only in combination with dont_align)
    fn_mask = getParam("--mask");

    // Missing data structures
    fn_missing = getParam("--missing");
    dim = getIntParam("--dim");
    maxres = getDoubleParam("--maxres");

    // Hidden arguments
    trymindiff_factor = getDoubleParam("--trymindiff_factor");

    // Number of threads
    threads = getIntParam("--thr");

}

void
ProgMLTomo::run()
{
    double LL, sumw_allrefs, sumcorr;
    bool converged;
    std::vector<double> conv;
    double wsum_sigma_noise, wsum_sigma_offset;
    std::vector<MultidimArray<double> > wsumimgs; //3D
    std::vector<MultidimArray<double> > wsumweds; //3D
    std::vector<MultidimArray<double> > fsc;
    MultidimArray<double> sumw; //1D
    MultidimArray<double> P_phi, Mr2, Maux; //3D
    FileName fn_img, fn_tmp;
    MultidimArray<double> oneline(0); //1D

    produceSideInfo();
    show();
    produceSideInfo2();
    Maux.resize(dim, dim, dim);
    Maux.setXmippOrigin();

    // Loop over all iterations
    for (int iter = istart; iter <= Niter; iter++)
    {

        if (verbose)
            std::cout << "  Multi-reference refinement:  iteration " << iter << " of "
            << Niter << std::endl;
        // Save old reference images
        for (int refno = 0; refno < nr_ref; refno++)
            Iold[refno]() = Iref[refno]();
        //#define DEBUG5
#ifdef DEBUG5

        try
        {
            //wrong
            MDimg.write("0.xmd");
        }
        catch (XmippError XE)
        {
            std::cout << XE;
            exit(0);
        }
#endif
        // Integrate over all images
        expectation(MDimg, Iref, iter, LL, sumcorr, wsumimgs, wsumweds,
                    wsum_sigma_noise, wsum_sigma_offset, sumw);
#ifdef DEBUG5

        try
        {
            //wrong
            MDimg.write("1.xmd");
        }
        catch (XmippError XE)
        {
            std::cout << XE;
            exit(0);
        }
#endif
        // Update model parameters
        maximization(wsumimgs, wsumweds, wsum_sigma_noise, wsum_sigma_offset, sumw,
                     sumcorr, sumw_allrefs, fsc, iter);
#ifdef DEBUG5

        try
        {
            //wrong
            MDimg.write("2.xmd");
        }
        catch (XmippError XE)
        {
            std::cout << XE;
            exit(0);
        }
#endif
        // Check convergence
        converged = checkConvergence(conv);
#ifdef DEBUG5

        try
        {
            //wrong
            MDimg.write("3.xmd");
        }
        catch (XmippError XE)
        {
            std::cout << XE;
            exit(0);
        }
#endif
        addPartialDocfileData(docfiledata, myFirstImg, myLastImg);

#ifdef DEBUG5

        try
        {
            //wrong
            MDimg.write("4.xmd");
        }
        catch (XmippError XE)
        {
            std::cout << XE;
            exit(0);
        }
#endif
        writeOutputFiles(iter, wsumweds, sumw_allrefs, LL, sumcorr, conv, fsc);
        if (converged)
        {
            if (verbose)
                std::cout << " Optimization converged!" << std::endl;
            break;
        }

    } // end loop iterations

    writeOutputFiles(-1, wsumweds, sumw_allrefs, LL, sumcorr, conv, fsc);

}

// Show ====================================================================
void
ProgMLTomo::show()
{

    if (verbose)
    {
        // To screen
        std::cout
        << " -----------------------------------------------------------------"
        << std::endl;
        std::cout
        << " | Read more about this program in the following publication:    |"
        << std::endl;
        std::cout
        << " |  Scheres ea.  (2009) Structure, 17, 1563-1572                 |"
        << std::endl;
        std::cout
        << " |                                                               |"
        << std::endl;
        std::cout
        << " |   *** Please cite it if this program is of use to you! ***    |"
        << std::endl;
        std::cout
        << " -----------------------------------------------------------------"
        << std::endl;
        std::cout << "--> Maximum-likelihood multi-reference refinement "
        << std::endl;
        std::cout << "  Input images            : " << fn_sel << " ("
        << nr_images_global << ")" << std::endl;
        if (!fn_ref.empty())
            std::cout << "  Reference image(s)      : " << fn_ref << std::endl;
        else
            std::cout << "  Number of references:   : " << nr_ref << std::endl;
        std::cout << "  Output rootname         : " << fn_root << std::endl;
        if (!(dont_align || do_only_average || dont_rotate))
        {
            std::cout << "  Angular sampling rate   : " << angular_sampling
            << " degrees" << std::endl;
            if (ang_search > 0.)
            {
                std::cout << "  Local angular searches  : " << ang_search << " degrees"
                << std::endl;
                if (!do_limit_psirange)
                    std::cout
                    << "                          : but with complete psi searches"
                    << std::endl;
            }
        }
        if (limit_trans >= 0.)
            std::cout << "  Maximum allowed shifts  : " << limit_trans << " pixels"
            << std::endl;
        std::cout << "  Symmetry group          : " << fn_sym << std::endl;
        std::cout << "  Stopping criterium      : " << eps << std::endl;
        std::cout << "  initial sigma noise     : " << sigma_noise << std::endl;
        std::cout << "  initial sigma offset    : " << sigma_offset << std::endl;
        std::cout << "  Maximum resolution      : " << maxres << " pix^-1"
        << std::endl;
        std::cout << "  Use images of size      : " << dim << std::endl;
        if (reg0 > 0.)
        {
            std::cout << "  Regularization from     : " << reg0 << " to " << regF
            << " in " << reg_steps << " steps" << std::endl;
        }
        if (!fn_missing.empty())
        {
            std::cout << "  Missing data info       : " << fn_missing << std::endl;
            if (do_impute && do_ml)
                std::cout << "  Missing data treatment  : imputation " << std::endl;
            else
                std::cout << "  Missing data treatment  : conventional division "
                << std::endl;
        }
        if (!fn_frac.empty())
            std::cout << "  Initial model fractions : " << fn_frac << std::endl;

        if (dont_rotate)
        {
            std::cout << "  -> Skip rotational searches, only translate and classify "
            << std::endl;
        }
        if (dont_align)
        {
            std::cout << "  -> Skip alignment, only classify " << std::endl;
            if (do_mask)
                std::cout << "  -> Classify within mask " << fn_mask << std::endl;
        }
        if (do_only_average)
        {
            std::cout
            << "  -> Skip alignment and classification, only calculate weighted average "
            << std::endl;
        }
        if (!do_ml)
        {
            std::cout
            << "  -> Use constrained correlation coefficient instead of ML-imputation approach."
            << std::endl;
        }
        if (do_perturb)
        {
            std::cout << "  -> Perturb angular sampling." << std::endl;
        }
        if (fix_fractions)
        {
            std::cout << "  -> Do not update estimates of model fractions."
            << std::endl;
        }
        if (fix_sigma_offset)
        {
            std::cout << "  -> Do not update sigma-estimate of origin offsets."
            << std::endl;
        }
        if (fix_sigma_noise)
        {
            std::cout << "  -> Do not update sigma-estimate of noise." << std::endl;
        }
        if (threads > 1)
        {
            std::cout << "  -> Using " << threads << " parallel threads" << std::endl;
        }

        std::cout
        << " -----------------------------------------------------------------"
        << std::endl;

    }

}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produce_Side_info2)
void
ProgMLTomo::produceSideInfo()
{

    FileName fn_img, fn_tmp, fn_base, fn_tmp2;
    Image<double> vol;
    Matrix1D<double> offsets(3), dum;
    MultidimArray<double> Maux, Maux2;
    MultidimArray<std::complex<double> > Faux;
    Matrix1D<int> center(3);
    MultidimArray<int> radial_count;
    size_t xdim, ydim, zdim;
    size_t ndim;

#ifdef  DEBUG

    std::cerr<<"Start produceSideInfo"<<std::endl;
#endif

    // For several random operations
    randomize_random_generator();

    // Read metadatafile with experimental images
    MDimg.read(fn_sel);
    MDimg.removeDisabled();
    nr_images_global = MDimg.size();
    //Get all images id for threads work on it
    MDimg.findObjects(imgs_id);

    // Check whether MDimg contains orientations
    //  mdimg_contains_angles = (MDimg.containsLabel(MDL_ANGLE_ROT)
    //      && MDimg.containsLabel(MDL_ANGLE_TILT) && MDimg.containsLabel(MDL_ANGLE_PSI)
    //      && MDimg.containsLabel(MDL_SHIFT_X) && MDimg.containsLabel(MDL_SHIFT_Y)
    //      && MDimg.containsLabel(MDL_SHIFT_Z));
    // Check angles and shifts are not present on input metadata, fill it with 0 value
    MDLabel labels[6] = { MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_SHIFT_X, MDL_SHIFT_Y, MDL_SHIFT_Z };
    String labelStr = "WARNING: Labels ";
    mdimg_contains_angles = true;
    for (int i = 0; i < 6; ++i)
    {
        if (!MDimg.containsLabel(labels[i]))
        {
            mdimg_contains_angles = false;
            MDimg.fillConstant(labels[i], "0.0");
            labelStr += MDL::label2Str(labels[i]) + " ";
            std::cerr << "missing label: "  << labels[i] << " " << labelStr<<std::endl;
        }
    }
    if (!mdimg_contains_angles)
        std::cerr << labelStr + "are mmissing in input metadata and are filled with value 0.0" <<std::endl;

    if (!MDimg.containsLabel(MDL_MISSINGREGION_NR))
    {
        if (MDmissing.size() > 1)
            REPORT_ERROR(ERR_MD_OBJECTNUMBER, "missingRegionNumber is missing from input images metadata"
                         "and your missing region metadata contains more than one missing region");
        int missno=0;
        MDmissing.getValue(MDL_MISSINGREGION_NR, missno, MDmissing.firstObject());
        MDimg.fillConstant(MDL_MISSINGREGION_NR, formatString("%d", missno));
        std::cerr << "WARNING: missingRegionNumber is missing from input images metadata,"
        "column filled with unique missing region provided" <<std::endl;
    }

    //WHY ROB <<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //MDimg.fillConstant(MDL_ANGLE_PSI, "0");

    // Get original dimension
    getImageSize(MDimg, xdim, ydim, zdim, ndim);
    if (xdim != ydim || xdim != zdim)
        REPORT_ERROR(ERR_MULTIDIM_SIZE, "Only cubic volumes are allowed");
    oridim = (int) xdim;

    // Downscaled dimension
    if (dim < 0)
        dim = oridim;
    if (dim > oridim)
        REPORT_ERROR(ERR_MULTIDIM_DIM,
                     "dim should be smaller than the size of the images");
    // Keep both dim and oridim even or uneven
    if (oridim % 2 != dim % 2)
        dim++;
    hdim = dim / 2;
    dim3 = dim * dim * dim;
    ddim3 = (double) dim3;
    scale_factor = (double) dim / (double) oridim;
    sigma_offset *= scale_factor;

    if (regF > reg0)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "regF should be smaller than reg0!");
    reg_current = reg0;

    // Make real-space masks
    MultidimArray<int> int_mask(dim, dim, dim);
    int_mask.setXmippOrigin();
    real_mask.resize(dim, dim, dim);
    real_mask.setXmippOrigin();
    BinaryCircularMask(int_mask, hdim, INNER_MASK);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(int_mask)
    {
        DIRECT_MULTIDIM_ELEM(real_mask,n) =
            (double) DIRECT_MULTIDIM_ELEM(int_mask,n);
    }
    real_omask.resize(real_mask);
    real_omask = 1. - real_mask;

    if (fn_mask.empty())
    {
        do_mask = false;
    }
    else
    {
        if (!dont_align)
            REPORT_ERROR(ERR_ARG_DEPENDENCE,
                         "option -mask is only valid in combination with -dont_align");
        Imask.read(fn_mask);
        if (Imask().computeMin() < 0. || Imask().computeMax() > 1.)
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         "mask should have values within the range [0,1]");
        Imask().setXmippOrigin();
        reScaleVolume(Imask(), true);
        // Remove any borders from the mask (to prevent problems rotating it later on)
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(int_mask)
        {
            DIRECT_MULTIDIM_ELEM(Imask(),n) *=
                (double) DIRECT_MULTIDIM_ELEM(int_mask,n);
        }
        nr_mask_pixels = Imask().sum();
    }

    // Set-up fourier-space mask
    MultidimArray<double> cosine_mask(dim, dim, dim);
    cosine_mask.setXmippOrigin();
    fourier_mask.resize(dim, dim, hdim + 1);
    fourier_imask.resize(dim, dim, hdim + 1);
    int r_maxres = XMIPP_MIN(hdim, FLOOR(maxres * dim));
    RaisedCosineMask(cosine_mask, r_maxres - 2, r_maxres, INNER_MASK);
    int yy, zz;
    int dimb = (dim - 1) / 2;
    for (int z = 0, ii = 0; z < dim; z++)
    {
        if (z > dimb)
            zz = z - dim;
        else
            zz = z;
        for (int y = 0; y < dim; y++)
        {
            if (y > dimb)
                yy = y - dim;
            else
                yy = y;
            for (int xx = 0; xx < hdim + 1; xx++, ii++)
            {
                if (xx <= FINISHINGX(cosine_mask))
                {
                    DIRECT_MULTIDIM_ELEM(fourier_mask,ii) = A3D_ELEM(cosine_mask,xx,yy,zz);
                    if (A3D_ELEM(cosine_mask,xx,yy,zz) > 0.)
                        DIRECT_MULTIDIM_ELEM(fourier_imask,ii) = 1;
                    else
                        DIRECT_MULTIDIM_ELEM(fourier_imask,ii) = 0;
                }
            }
        }
    }
    // exclude origin voxel from fourier masks
    DIRECT_MULTIDIM_ELEM(fourier_mask,0) = 0.;
    DIRECT_MULTIDIM_ELEM(fourier_imask,0) = 0;

    // Precalculate sampling
    do_sym = false;
    if (dont_align || do_only_average || dont_rotate)
    {
        if (!mdimg_contains_angles)
            REPORT_ERROR(
                ERR_ARG_MISSING,
                "Options --dont_align, --dont_rotate and --only_average require that angle information is present in input images metadata");
        ang_search = -1.;
    }
    else
    {
        mysampling.setSampling(angular_sampling);
        if (!mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         (std::string)"Invalid symmetry " + fn_sym);
        mysampling.SL.readSymmetryFile(fn_sym);
        // Check whether we are using symmetry
        if (mysampling.SL.symsNo() > 0)
        {
            do_sym = true;
            if (verbose)
            {
                std::cout << "  --> Using symmetry version of the code: " << std::endl;
                std::cout
                << "      [-psi, -tilt, -rot] is applied to the references, thereby "
                << std::endl;
                std::cout << "      tilt and psi are sampled on an hexagonal lattice, "
                << std::endl;
                std::cout << "      and rot is sampled linearly " << std::endl;
                std::cout
                << "      Note that --dont_limit_psirange option is not allowed.... "
                << std::endl;
            }
            if (!do_limit_psirange && ang_search > 0.)
                REPORT_ERROR(ERR_ARG_INCORRECT,
                             "exhaustive psi-angle search only allowed for C1 symmetry");
        }
        mysampling.fillLRRepository();
        // by default max_tilt= +91., min_tilt= -91.
        mysampling.computeSamplingPoints(false, // half sphere?
                                         tilt_rangeF, tilt_range0);
        mysampling.removeRedundantPointsExhaustive(symmetry, sym_order, false, // half sphere?
                0.75 * angular_sampling);
        if (psi_sampling < 0)
            psi_sampling = angular_sampling;

    }
    readMissingInfo();////////////////////////////
    // Get number of references
    if (do_only_average)
    {
        int refno;
        nr_ref = 0;
        FOR_ALL_OBJECTS_IN_METADATA(MDimg)
        {
            if (MDimg.getValue(MDL_REF, refno, __iter.objId))
            {
                nr_ref = XMIPP_MAX(refno, nr_ref);
            }
            else
            {
                nr_ref = 1;
                break;
            }
        }
        fn_ref = "tt";
        Niter = 1;
        do_impute = false;
        do_ml = false;
    }
    else
    {
    	do_generate_refs = fn_ref.empty();
        if (do_generate_refs) //assign bool value and asking at same time
            generateInitialReferences();
    }



} //end of function produceSideInfo

// Generate initial references =============================================
void
ProgMLTomo::generateInitialReferences()
{

    //MetaData SFtmp;
    Image<double> Iave1, Iave2, Itmp, Iout;
    MultidimArray<double> Msumwedge1, Msumwedge2;
    MultidimArray<unsigned char> Mmissing;
    MultidimArray<std::complex<double> > Fave;
    std::vector<MultidimArray<double> > fsc; //1D
    Matrix2D<double> my_A; //2D
    Matrix1D<double> my_offsets(3); //1D
    FileName fn_tmp;
    //SelLine line;
    MetaData MDrand;
    std::vector<MetaData> MDtmp(nr_ref);
    int Nsub = ROUND((double)nr_images_global / nr_ref);
    double my_rot, my_tilt, my_psi, resolution;
    //DocLine DL;
    int missno, c = 0, cc = XMIPP_MAX(1, nr_images_global / 60);

#ifdef  DEBUG

    std::cerr<<"Start generateInitialReferences"<<std::endl;
#endif
#ifdef  DEBUG_JM

    std::cerr << "DEBUG_JM: entering ProgMLTomo::generateInitialReferences" <<std::endl;
#endif

    if (verbose)
    {
        std::cout
        << "  Generating initial references by averaging over random subsets"
        << std::endl;
        init_progress_bar(nr_images_global);
    }

    fsc.clear();
    fsc.resize(nr_ref + 1);
    Iave1().resize(dim, dim, dim);
    Iave1().setXmippOrigin();
    Iave2().resize(dim, dim, dim);
    Iave2().setXmippOrigin();
    Msumwedge1.resize(dim, dim, hdim + 1);
    Msumwedge2.resize(dim, dim, hdim + 1);
    // Make random subsets and calculate average images in random orientations
    // and correct using conventional division by the sum of the wedges
    MDrand.randomize(MDimg);
    MDrand.split(nr_ref, MDtmp);
    MDref.clear(); //Just to be sure

    for (int refno = 0; refno < nr_ref; ++refno)
    {
        Iave1().initZeros();
        Iave2().initZeros();
        Msumwedge1.initZeros();
        Msumwedge2.initZeros();

        MetaData &md = MDtmp[refno];
        FOR_ALL_OBJECTS_IN_METADATA(md)
        {
            md.getValue(MDL_IMAGE, fn_tmp, __iter.objId);
            Itmp.read(fn_tmp);
            Itmp().setXmippOrigin();
            reScaleVolume(Itmp(), true);
            if (do_keep_angles || dont_align || dont_rotate)
            {
                // angles from MetaData
                md.getValue(MDL_ANGLE_ROT, my_rot, __iter.objId);
                md.getValue(MDL_ANGLE_TILT, my_tilt, __iter.objId);
                md.getValue(MDL_ANGLE_PSI, my_psi, __iter.objId);

                md.getValue(MDL_SHIFT_X, my_offsets(0), __iter.objId);
                md.getValue(MDL_SHIFT_Y, my_offsets(1), __iter.objId);
                md.getValue(MDL_SHIFT_Z, my_offsets(2), __iter.objId);
                my_offsets *= scale_factor;

                selfTranslate(LINEAR, Itmp(), my_offsets, DONT_WRAP);
            }
            else
            {
                // reset to random angles
                my_rot = 360. * rnd_unif(0., 1.);
                my_tilt = ACOSD((2.*rnd_unif(0., 1.) - 1));
                my_psi = 360. * rnd_unif(0., 1.);
            }

            Euler_angles2matrix(my_rot, my_tilt, my_psi, my_A, true);
            selfApplyGeometry(LINEAR, Itmp(), my_A, IS_NOT_INV, DONT_WRAP, 0.);

            int iran_fsc = ROUND(rnd_unif());
            if (iran_fsc == 0)
                Iave1() += Itmp();
            else
                Iave2() += Itmp();
            if (do_missing)
            {
                md.getValue(MDL_MISSINGREGION_NR, missno, __iter.objId);
                --missno;
                getMissingRegion(Mmissing, my_A, missno);
                if (iran_fsc == 0)
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Msumwedge1)
                    DIRECT_MULTIDIM_ELEM(Msumwedge1,n)+=DIRECT_MULTIDIM_ELEM(Mmissing,n);
                else
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Msumwedge2)
                    DIRECT_MULTIDIM_ELEM(Msumwedge2,n)+=DIRECT_MULTIDIM_ELEM(Mmissing,n);
            }
            c++;
            if (verbose && c % cc == 0)
                progress_bar(c);
        }

        // Calculate resolution
        calculateFsc(Iave1(), Iave2(), Msumwedge1, Msumwedge2, fsc[0],
                     fsc[refno + 1], resolution);
        Iave1() += Iave2();
        Msumwedge1 += Msumwedge2;

        if (do_missing)
        {
            // 1. Correct for missing wedge by division of sumwedge
            transformer.FourierTransform(Iave1(), Fave, true);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fave)
            {
                if (DIRECT_MULTIDIM_ELEM(Msumwedge1,n) > noimp_threshold)
                {
                    DIRECT_MULTIDIM_ELEM(Fave,n) /= DIRECT_MULTIDIM_ELEM(Msumwedge1,n);
                }
                else
                {
                    DIRECT_MULTIDIM_ELEM(Fave,n) = 0.;
                }
            }
            transformer.inverseFourierTransform(Fave, Iave1());
        }
        else
        {
            Iave1() /= (double) Nsub;
        }

        // Enforce fourier_mask, symmetry and omask
        if (do_ini_filter)
            postProcessVolume(Iave1, resolution);
        else
            postProcessVolume(Iave1);

        fn_tmp = FN_ITER_VOL(0, "ref", refno);
        Iout = Iave1;
        reScaleVolume(Iout(), false);
        Iout.write(fn_tmp);
        MDref.setValue(MDL_IMAGE, fn_tmp, MDref.addObject());
        // Also write out average wedge of this reference
        fn_tmp = FN_ITER_VOL(0, "wedge", refno);
        Iout() = Msumwedge1;
        reScaleVolume(Iout(), false);
        Iout.write(fn_tmp);
    }
    if (verbose)
        progress_bar(nr_images_global);
    fn_ref = FN_ITER_MD(0);
    MDref.write(fn_ref);

#ifdef  DEBUG

    std::cerr<<"Finished generateInitialReferences"<<std::endl;
#endif
}

/// Set the number of images, this function is useful only for MPI
void
ProgMLTomo::setNumberOfLocalImages()
{
    nr_images_local = nr_images_global;
    // By default set myFirst and myLast equal to 0 and N - 1
    // respectively, this should be changed when using MPI
    // by calling setWorkingImages before produceSideInfo2
    myFirstImg = 0;
    myLastImg = nr_images_global - 1;
}
// Read reference images to memory and initialize offset vectors
// This side info is NOT general, i.e. in parallel mode it is NOT the
// same for all processors! (in contrast to produce_Side_info)
void
ProgMLTomo::produceSideInfo2(int nr_vols)
{

#ifdef  JM_DEBUG

    std::cerr<<"Start produceSideInfo2"<<std::endl;
#endif

    MetaData DF, DFsub;
    FileName fn_tmp;
    Image<double> img, Vaux;
    std::vector<Matrix1D<double> > Vdm;

    setNumberOfLocalImages();
    // Store tomogram angles, offset vectors and missing wedge parameters
    imgs_optrefno.assign(nr_images_local, 0.);
    imgs_optangno.assign(nr_images_local, 0.);
    imgs_optpsi.assign(nr_images_local, 0.);
    imgs_trymindiff.assign(nr_images_local, -1.);
    Matrix1D<double> dum(3);
    if (do_missing)
    {
        if (!MDimg.containsLabel(MDL_MISSINGREGION_NR))
            REPORT_ERROR(
                ERR_MD_BADLABEL,
                "Expecting missing region label 'missingRegionNumber' on input metadata");
        imgs_missno.assign(nr_images_local, -1);
    }
    if (dont_align || dont_rotate || do_only_average)
        imgs_optoffsets.assign(nr_images_local, dum);
    //--------Setup for Docfile -----------
    docfiledata.initZeros(nr_images_local, MLTOMO_DATALINELENGTH);
    // Read in MetaData with wedge info, optimal angles, etc.
    size_t count = 0, imgno;
    int refno;
    MetaData MDsub;
    MDRow row;
    double rot, tilt, psi;

    //FOR_ALL_OBJECTS_IN_METADATA(MDimg)
    for (size_t img = myFirstImg; img <= myLastImg; ++img)
    {

        imgno = img - myFirstImg;
        MDimg.getRow(row, imgs_id[imgno]);
        // Get missing wedge type
        if (do_missing)
        {
            row.getValue(MDL_MISSINGREGION_NR, imgs_missno[imgno]);
            //The label MDL_MISSINGREGION_NR should be the index of missing regions, starting at 1
            --imgs_missno[imgno];
        }
        // Generate a MetaData from the (possible subset of) images in SF
        if (ang_search > 0. || dont_align || do_only_average || dont_rotate)
        {
            row.getValue(MDL_ANGLE_ROT, rot);
            row.getValue(MDL_ANGLE_PSI, psi);
            row.getValue(MDL_ANGLE_TILT, tilt);

            if (do_sym)
            {
                imgs_optpsi[imgno] = rot; // for limited psi (now rot) searches
                // Reverse rotation applied to the references!!
                row.setValue(MDL_ANGLE_ROT, -psi);
                row.setValue(MDL_ANGLE_TILT, -tilt);
                row.setValue(MDL_ANGLE_PSI, -rot);
            }
            else
                imgs_optpsi[imgno] = psi; // for limited psi searches
            std::cerr << "r t p " << rot << " " << tilt << " " << psi <<std::endl;
            MDsub.addRow(row);
            imgs_optangno[imgno] = count++;

            row.getValue(MDL_REF, refno);
            --refno;
            imgs_optrefno[imgno] = (refno < 0 || refno >= nr_ref) ? 0 : refno;

            if (dont_align || dont_rotate || do_only_average)
            {
                row.getValue(MDL_SHIFT_X, imgs_optoffsets[imgno](0));
                row.getValue(MDL_SHIFT_Y, imgs_optoffsets[imgno](1));
                row.getValue(MDL_SHIFT_Z, imgs_optoffsets[imgno](2));
            }
        }
    }////////////////////<<<<<<<<<<<<<<<<<<<<<<<<
    // Set up local searches
    if (ang_search > 0.)
    {
        std::cerr << "myLastImg" << myLastImg <<std::endl;
        std::cerr << "ang_search " << ang_search << std::endl;
        mysampling.setNeighborhoodRadius(ang_search);
        std::cerr << "MDsub.size " << MDsub.size() << std::endl;

        mysampling.fillExpDataProjectionDirectionByLR(MDsub);
        std::cerr << "exp_data_projection_direction_by_L_R "
        << mysampling.exp_data_projection_direction_by_L_R.size()
        <<std::endl;
        mysampling.removePointsFarAwayFromExperimentalData();
        std::cerr << "no_redundant_sampling_points_vector"
        << mysampling.no_redundant_sampling_points_vector.size() <<std::endl;
        mysampling.computeNeighbors();
        //mysampling.saveSamplingFile("NEW");
        //MDsub.write("NEWexp.xmd");
        convert_refno_to_stack_position.resize(mysampling.numberSamplesAsymmetricUnit, -1);
        for (size_t i = 0; i < mysampling.no_redundant_sampling_points_index.size(); i++)
            convert_refno_to_stack_position[i] =
                mysampling.no_redundant_sampling_points_index[i];

        //exit(1);
    }


    // Calculate angular sampling
    // A. from MetaData entries (MDsub) only
    if (dont_align || dont_rotate || do_only_average)
    {

        AnglesInfo myinfo;
        all_angle_info.clear();
        nr_ang = 0;
        /*DFsub.go_first_data_line();
         while (!DFsub.eof())*/
        FOR_ALL_OBJECTS_IN_METADATA(MDsub)
        {
            MDsub.getRow(row, __iter.objId);
            row.getValue(MDL_ANGLE_ROT, myinfo.rot);
            row.getValue(MDL_ANGLE_TILT, myinfo.tilt);
            row.getValue(MDL_ANGLE_PSI, myinfo.psi);
            if (do_sym)
            {
                rot = myinfo.rot;
                myinfo.rot = -myinfo.psi;
                myinfo.tilt *= -1;
                myinfo.psi = -rot;
            }
            Euler_angles2matrix(myinfo.rot, myinfo.tilt, myinfo.psi, myinfo.A, true);
            myinfo.direction = nr_ang++;
            all_angle_info.push_back(myinfo);
        }
    }
    // B. from mysampling
    else
    {

        int nr_psi = CEIL(360. / psi_sampling);
        AnglesInfo myinfo;
        all_angle_info.clear();
        nr_ang = 0;
        for (size_t i = 0; i < mysampling.no_redundant_sampling_points_angles.size(); i++)
        {
            rot = XX(mysampling.no_redundant_sampling_points_angles[i]);
            tilt = YY(mysampling.no_redundant_sampling_points_angles[i]);
            for (int ipsi = 0; ipsi < nr_psi; ipsi++)
            {
                psi = (double) (ipsi * 360. / nr_psi);
                if (!no_SMALLANGLE)
                    psi += SMALLANGLE;
                if (do_sym)
                {
                    // Inverse rotation because A_rot is being applied to the reference
                    // and rot, tilt and psi apply to the experimental subtomogram!
                    myinfo.rot  = -psi;
                    myinfo.tilt = -tilt;
                    myinfo.psi  = -rot;
                }
                else
                {
                    // Only in C1 we get away with reverse order rotations...
                    // This is ugly, but allows -dont_limit_psirange....
                    myinfo.rot  = rot;
                    myinfo.tilt = tilt;
                    myinfo.psi  = psi;
                }
                Euler_angles2matrix(myinfo.rot, myinfo.tilt, myinfo.psi, myinfo.A,
                                    true);
                if (ang_search > 0.)
                    myinfo.direction = convert_refno_to_stack_position[i];
                else
                    myinfo.direction = i;////convert_refno_to_stack_position
                all_angle_info.push_back(myinfo);
                nr_ang++;
            }
        }
    }

    // Copy all rot, tilr, psi and A into _ori equivalents
    for (int angno = 0; angno < nr_ang; angno++)
    {
        all_angle_info[angno].rot_ori  = all_angle_info[angno].rot;
        all_angle_info[angno].tilt_ori = all_angle_info[angno].tilt;
        all_angle_info[angno].psi_ori  = all_angle_info[angno].psi;
        all_angle_info[angno].A_ori    = all_angle_info[angno].A;
    }

#ifdef JM_DEBUG  //////////////////////////////////////////////////////////////////

    MetaData DFt;
    size_t _id;
    for (int angno = 0; angno < nr_ang; angno++)
    {
        _id = DFt.addObject();

        double rot=all_angle_info[angno].rot;
        double tilt=all_angle_info[angno].tilt;
        double psi=all_angle_info[angno].psi;
        DFt.setValue(MDL_ANGLE_ROT,rot,_id);
        DFt.setValue(MDL_ANGLE_TILT,tilt,_id);
        DFt.setValue(MDL_ANGLE_PSI,psi,_id);
        DFt.setValue(MDL_ANGLE_PSI,psi,_id);
        //        DFt.append_angles(rot,tilt,psi,"rot","tilt","psi");
    }
    FileName fnt;
    fnt = fn_root + "_angles1.doc";
    DFt.write(fnt);
    mysampling.saveSamplingFile(fn_root + "_sampling.doc");
#endif

    // Prepare reference images
    if (do_only_average)
    {
        img().initZeros(dim, dim, dim);
        img().setXmippOrigin();
        for (int refno = 0; refno < nr_ref; refno++)
        {
            Iref.push_back(img);
            Iold.push_back(img);
        }
    }
    else
    {
        // Read in all reference images in memory
        MDref.read(fn_ref);
        nr_ref = MDref.size();
        FileName fn_img;
        FOR_ALL_OBJECTS_IN_METADATA(MDref)
        {
            MDref.getValue(MDL_IMAGE, fn_img, __iter.objId);
            img.read(fn_img);
            img().setXmippOrigin();

            if (do_mask)
            {
                img() *= Imask();
            }
            reScaleVolume(img(), true);
            // From now on we will assume that all references are omasked, so enforce this here
            maskSphericalAverageOutside(img());
            // Rotate some arbitrary (off-axis) angle and rotate back again to remove high frequencies
            // that will be affected by the interpolation due to rotation
            // This makes that the A2 values of the rotated references are much less sensitive to rotation
            Matrix2D<double> E_rot;
            Euler_angles2matrix(32., 61., 53., E_rot, true);
            selfApplyGeometry(LINEAR, img(), E_rot, IS_NOT_INV, DONT_WRAP,
                              DIRECT_MULTIDIM_ELEM(img(), 0));
            selfApplyGeometry(LINEAR, img(), E_rot, IS_INV, DONT_WRAP,
                              DIRECT_MULTIDIM_ELEM(img(), 0));
            Iref.push_back(img);
            Iold.push_back(img);

        }
    }

    // Prepare prior alpha_k
    alpha_k.resize(nr_ref);
    if (MDref.containsLabel(MDL_WEIGHT))
        //if (fn_frac != "")
    {
        // read in model fractions if given on command line
        double sumfrac = 0.;
        //MetaData DF;
        //DocLine DL;
        //DF.read(fn_frac);
        //DF.go_first_data_line();
        //for (int refno = 0; refno < nr_ref; refno++)
        int refno = 0;
        FOR_ALL_OBJECTS_IN_METADATA(MDref)
        {
            MDref.getValue(MDL_WEIGHT, alpha_k(refno), __iter.objId);
            //DL = DF.get_current_line();
            //alpha_k(refno) = DL[0];
            sumfrac += alpha_k(refno);
            ++refno;
            //DF.next_data_line();
        }
        if (ABS(sumfrac - 1.) > 1e-3)
            if (verbose)
                std::cerr << " ->WARNING: Sum of all expected model fractions ("
                << sumfrac << ") is not one!" << std::endl;
        for (int refno = 0; refno < nr_ref; refno++)
        {
            alpha_k(refno) /= sumfrac;
        }
    }
    else
    {
        // Even distribution
        alpha_k.initConstant(1. / (double) nr_ref);
    }

    //exit(1);

    // Regularization (do not regularize during restarts!)
    if (istart == 1)
        regularize(istart - 1);

    //#define DEBUG_GENERAL
#ifdef JM_DEBUG

    //    std::cerr<<"nr images= "<<SF.ImgNo()<<std::endl;
    std::cerr<<"do_generate_refs ="<<do_generate_refs<<std::endl;
    std::cerr<<"nr_ref= "<<nr_ref<<std::endl;
    std::cerr<<"nr_miss= "<<nr_miss<<std::endl;
    std::cerr<<"dim= "<<dim<<std::endl;
    std::cerr<<"Finished produceSideInfo"<<std::endl;
    std::cerr<<"nr_ang= "<<nr_ang<<std::endl;
    //    std::cerr<<"nr_psi= "<<nr_psi<<std::endl;
#endif

#ifdef JM_DEBUG

    std::cerr<<"Finished produceSideInfo2"<<std::endl;
#endif

}

void
ProgMLTomo::perturbAngularSampling()
{

    Matrix2D<double> R, I(3, 3);
    double ran1, ran2, ran3;
    I.initIdentity();

    // Note that each mpi process will have a different random perturbation!!
    ran1 = rnd_gaus(0., angular_sampling / 3.);
    ran2 = rnd_gaus(0., angular_sampling / 3.);
    ran3 = rnd_gaus(0., angular_sampling / 3.);
    Euler_angles2matrix(ran1, ran2, ran3, R, true);
    R.resize(3, 3);

    //#define DEBUG_PERTURB
#ifdef DEBUG_PERTURB

    std::cerr<<" random perturbation= "<<ran1<<" "<<ran2<<" "<<ran3<<std::endl;
#endif

    for (int angno = 0; angno < nr_ang; angno++)
    {
        Euler_apply_transf(I, R, all_angle_info[angno].rot_ori,
                           all_angle_info[angno].tilt_ori, all_angle_info[angno].psi_ori,
                           all_angle_info[angno].rot, all_angle_info[angno].tilt,
                           all_angle_info[angno].psi);
        Euler_angles2matrix(all_angle_info[angno].rot, all_angle_info[angno].tilt,
                            all_angle_info[angno].psi, all_angle_info[angno].A, true);
    }

}

void
ProgMLTomo::readMissingInfo()
{
    // Read in MetaData with information about the missing wedges
    nr_miss = 0;
    do_missing = !fn_missing.empty();
    if (do_missing) //if provided missing wedges metadata
    {
        MDmissing.read(fn_missing);
        nr_miss = MDmissing.size();
        if (nr_miss <= 0)
            REPORT_ERROR(ERR_ARG_INCORRECT, "Empty metadata with missing wedges info");

        MultidimArray<double> Mcomplete(dim, dim, dim);
        MultidimArray<unsigned char> Mmissing(dim, dim, hdim + 1);
        Matrix2D<double> I(4, 4);
        I.initIdentity();
        //Create a clean missing info
        MissingInfo emptyMissingInfo;
        emptyMissingInfo.do_limit_x = emptyMissingInfo.do_limit_y =
                                          emptyMissingInfo.do_cone = false;
        emptyMissingInfo.thx0 = emptyMissingInfo.thxF = emptyMissingInfo.thy0 =
                                    emptyMissingInfo.thyF = 0.;
        emptyMissingInfo.tg0_x = emptyMissingInfo.tgF_x = emptyMissingInfo.tg0_y =
                                     emptyMissingInfo.tgF_y = 0.;
        emptyMissingInfo.nr_pixels = 0.;
        all_missing_info.assign(nr_miss, emptyMissingInfo);

        int missno = 0;
        String missingType;
        MDRow row;

        FOR_ALL_OBJECTS_IN_METADATA(MDmissing)
        {
            //Note: Now we are assuming that all missing regions will be numerate
            //from 0 to n - 1 (missno), and images belonging to missno 0 will
            //have the value 1 in MDL_MISSINGREGION_NR and so on...
            MissingInfo &myinfo = all_missing_info[missno];
            MDmissing.getRow(row, __iter.objId);
            //MDmissing.getValue(MDL_MISSINGREGION_NR, missno, __iter.objId);
            row.getValue(MDL_MISSINGREGION_TYPE, missingType);
            if (missingType == "wedge_y")
                myinfo.type = MISSING_WEDGE_Y;
            else if (missingType == "wedge_x")
                myinfo.type = MISSING_WEDGE_X;
            else if (missingType == "pyramid")
                myinfo.type = MISSING_PYRAMID;
            else if (missingType == "cone")
                myinfo.type = MISSING_CONE;
            else
                REPORT_ERROR(
                    ERR_ARG_INCORRECT,
                    formatString("Unrecognized type of missing region: '%s'", missingType.c_str()));

            if (myinfo.type == MISSING_WEDGE_Y || myinfo.type == MISSING_PYRAMID)
            {
                myinfo.do_limit_x = true;
                row.getValue(MDL_MISSINGREGION_THY0, myinfo.thy0);
                row.getValue(MDL_MISSINGREGION_THYF, myinfo.thyF);
                myinfo.tg0_y = -tan(PI * (-90. - myinfo.thyF) / 180.);
                myinfo.tgF_y = -tan(PI * (90. - myinfo.thy0) / 180.);
            }
            if (myinfo.type == MISSING_WEDGE_X || myinfo.type == MISSING_PYRAMID)
            {
                myinfo.do_limit_y = true;
                row.getValue(MDL_MISSINGREGION_THX0, myinfo.thx0);
                row.getValue(MDL_MISSINGREGION_THXF, myinfo.thxF);
                myinfo.tg0_x = -tan(PI * (-90. - myinfo.thxF) / 180.);
                myinfo.tgF_x = -tan(PI * (90. - myinfo.thx0) / 180.);
            }
            if (myinfo.type == MISSING_CONE)
            {
                myinfo.do_cone = true;
                row.getValue(MDL_MISSINGREGION_THY0, myinfo.thy0);
                myinfo.tg0_y = -tan(PI * (-90. - myinfo.thy0) / 180.);
            }

            getMissingRegion(Mmissing, I, missno);

            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Mcomplete)
            {
                if (j < XSIZE(Mmissing)
                   )
                    myinfo.nr_pixels += DIRECT_A3D_ELEM(Mmissing,k,i,j);
                else
                    myinfo.nr_pixels +=
                        DIRECT_A3D_ELEM( Mmissing,(dim-k)%dim,(dim-i)%dim,dim-j);
            }

            ++missno;
        }
    }
}

void
ProgMLTomo::getMissingRegion(MultidimArray<unsigned char> &Mmissing,
                             const Matrix2D<double> &A, const int missno)
{
#ifdef DEBUG_JM
    std::cerr << "   DEBUG_JM: entering ProgMLTomo::getMissingRegion" <<std::endl;
    std::cerr << "   DEBUG_JM:    missno: " << missno << std::endl;
#endif

    if (missno < 0 || missno >= nr_miss)
        REPORT_ERROR(
            ERR_ARG_INCORRECT,
            "ProgMLTomo::getMissingRegion: missno should be between 0 and nr_miss - 1");

    Matrix2D<double> Ainv = A.inv();
    double xp, yp, zp;
    double limx0 = 0., limxF = 0., limy0 = 0., limyF = 0., lim = 0.;
    bool is_observed;
    Mmissing.resizeNoCopy(dim, dim, hdim + 1);

    MissingInfo &myinfo = all_missing_info[missno];

    //#define DEBUG_WEDGE
#ifdef DEBUG_WEDGE

    std::cerr<<"do_limit_x= "<<do_limit_x<<std::endl;
    std::cerr<<"do_limit_y= "<<do_limit_y<<std::endl;
    std::cerr<<"do_cone= "<<do_cone<<std::endl;
    std::cerr<<"tg0_y= "<<tg0_y<<std::endl;
    std::cerr<<"tgF_y= "<<tgF_y<<std::endl;
    std::cerr<<"tg0_x= "<<tg0_x<<std::endl;
    std::cerr<<"tgF_x= "<<tgF_x<<std::endl;
    std::cerr<<"XMIPP_EQUAL_ACCURACY= "<<XMIPP_EQUAL_ACCURACY<<std::endl;
    std::cerr<<"Ainv= "<<Ainv<<" A= "<<A<<std::endl;
#endif

    int zz, yy;
    int dimb = (dim - 1) / 2;
    double Ainv_zz_x, Ainv_zz_y, Ainv_zz_z;
    double Ainv_yy_x, Ainv_yy_y, Ainv_yy_z;

    for (int z = 0, ii = 0; z < dim; ++z)
    {
        zz = (z > dimb ) ? z - dim : z;
        Ainv_zz_x = dMij(Ainv, 0, 2) * zz;
        Ainv_zz_y = dMij(Ainv, 1, 2) * zz;
        Ainv_zz_z = dMij(Ainv, 2, 2) * zz;

        for (int y = 0; y < dim; ++y)
        {
            yy = (y > dimb ) ? y - dim : y;
            Ainv_yy_x = Ainv_zz_x + dMij(Ainv, 0, 1) * yy;
            Ainv_yy_y = Ainv_zz_y + dMij(Ainv, 1, 1) * yy;
            Ainv_yy_z = Ainv_zz_z + dMij(Ainv, 2, 1) * yy;

            for (int xx = 0; xx < hdim + 1; ++xx, ++ii)
            {

                unsigned char maskvalue = DIRECT_MULTIDIM_ELEM(fourier_imask,ii);
                if (!maskvalue)
                {
                    DIRECT_MULTIDIM_ELEM(Mmissing,ii) = 0;
                }
                else
                {
                    // Rotate the wedge
                    xp = dMij(Ainv, 0, 0) * xx + Ainv_yy_x;
                    yp = dMij(Ainv, 1, 0) * xx + Ainv_yy_y;
                    zp = dMij(Ainv, 2, 0) * xx + Ainv_yy_z;

                    // Calculate the limits
                    if (myinfo.do_cone)
                    {
                        lim = myinfo.tg0_y * zp;
                        lim *= lim;
                        is_observed = (xp * xp + yp * yp >= lim);
                    }
                    else
                    {
                        is_observed = false; // for pyramid
                        if (myinfo.do_limit_x)
                        {
                            limx0 = myinfo.tg0_y * zp;
                            limxF = myinfo.tgF_y * zp;
                            if (zp >= 0)
                                is_observed = (xp <= limx0 || xp >= limxF);
                            else
                                is_observed = (xp <= limxF || xp >= limx0);
                        }
                        if (myinfo.do_limit_y && !is_observed)
                        {
                            limy0 = myinfo.tg0_x * zp;
                            limyF = myinfo.tgF_x * zp;
                            if (zp >= 0)
                                is_observed = (yp <= limy0 || yp >= limyF);
                            else
                                is_observed = (yp <= limyF || yp >= limy0);
                        }
                    }

                    if (is_observed)
                        DIRECT_MULTIDIM_ELEM(Mmissing,ii) = maskvalue;
                    else
                        DIRECT_MULTIDIM_ELEM(Mmissing,ii) = 0;
                }
            }
        }
    }

#ifdef DEBUG_WEDGE
    Image<double> test(dim,dim,dim), ori;
    ori()=Mmissing;
    ori.write("oriwedge.fft");
    test().initZeros();
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(test())
    if (j<hdim+1)
        DIRECT_A3D_ELEM(test(),k,i,j)=
            DIRECT_A3D_ELEM(ori(),k,i,j);
    else
        DIRECT_A3D_ELEM(test(),k,i,j)=
            DIRECT_A3D_ELEM(ori(),
                            (dim-k)%dim,
                            (dim-i)%dim,
                            dim-j);
    test.write("wedge.ftt");
    CenterFFT(test(),true);
    test.write("Fwedge.vol");
    test().setXmippOrigin();
    MultidimArray<int> ress(dim,dim,dim);
    ress.setXmippOrigin();
    //ROB
    //    BinaryWedgeMask(test(),all_missing_info[missno].thy0, all_missing_info[missno].thyF, A);
    BinaryWedgeMask(test(),all_missing_info[missno].thy0, all_missing_info[missno].thyF, A.inv());
    test.write("Mwedge.vol");
#endif
#ifdef DEBUG_JM

    std::cerr << "   DEBUG_JM: leaving ProgMLTomo::getMissingRegion" <<std::endl;
#endif
}

// Calculate probability density function of all in-plane transformations phi
void
ProgMLTomo::calculatePdfTranslations()
{

#ifdef  DEBUG
    std::cerr<<"start calculatePdfTranslations"<<std::endl;
#endif

    double r2, pdfpix;
    P_phi.resize(dim, dim, dim);
    P_phi.setXmippOrigin();
    Mr2.resize(dim, dim, dim);
    Mr2.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY3D(P_phi)
    {
        r2 = (double) (j * j + i * i + k * k);

        if (limit_trans >= 0. && r2 > limit_trans * limit_trans)
            A3D_ELEM(P_phi, k, i, j) = 0.;
        else
        {
            if (sigma_offset > XMIPP_EQUAL_ACCURACY)
            {
                pdfpix = exp(-r2 / (2. * sigma_offset * sigma_offset));
                pdfpix *= pow(2. * PI * sigma_offset * sigma_offset, -3. / 2.);
            }
            else
            {
                if (k == 0 && i == 0 && j == 0)
                    pdfpix = 1.;
                else
                    pdfpix = 0.;
            }
            A3D_ELEM(P_phi, k, i, j) = pdfpix;
            A3D_ELEM(Mr2, k, i, j) = (float) r2;
        }
    }

    // Re-normalize for limit_trans
    if (limit_trans >= 0.)
    {
        double sum = P_phi.sum();
        P_phi /= sum;
    }

    //#define  DEBUG_PDF_SHIFT
#ifdef DEBUG_PDF_SHIFT
    std::cerr<<" Sum of translation pdfs (should be one) = "<<P_phi.sum()<<std::endl;
    Image<double> Vt;
    Vt()=P_phi;
    Vt.write("pdf.vol");
#endif

#ifdef  DEBUG

    std::cerr<<"finished calculatePdfTranslations"<<std::endl;
#endif

}

void
ProgMLTomo::maskSphericalAverageOutside(MultidimArray<double> &Min)
{
    double outside_density = 0., sumdd = 0.;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(real_omask)
    {
        outside_density += DIRECT_MULTIDIM_ELEM(Min,n)
                           * DIRECT_MULTIDIM_ELEM(real_omask,n);
        sumdd += DIRECT_MULTIDIM_ELEM(real_omask,n);
    }
    outside_density /= sumdd;

    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(real_omask)
    {
        DIRECT_MULTIDIM_ELEM(Min,n) *= DIRECT_MULTIDIM_ELEM(real_mask,n);
        DIRECT_MULTIDIM_ELEM(Min,n) += (outside_density)
                                       * DIRECT_MULTIDIM_ELEM(real_omask,n);
    }
}

void
ProgMLTomo::reScaleVolume(MultidimArray<double> &Min, bool down_scale)
{
    MultidimArray<std::complex<double> > Fin;
    MultidimArray<double> Mout;
    FourierTransformer local_transformer_in, local_transformer_out;

    if (oridim == dim)
        return;

    int newdim = down_scale ? dim : oridim;

    local_transformer_in.setReal(Min);
    local_transformer_in.FourierTransform();
    local_transformer_in.getCompleteFourier(Fin);
    Mout.resize(newdim, newdim, newdim);
    Mout.setXmippOrigin();
    local_transformer_out.setReal(Mout);

    CenterFFT(Fin, true);
    Fin.setXmippOrigin();
    Fin.selfWindow(STARTINGZ(Mout), STARTINGY(Mout), STARTINGX(Mout),
                   FINISHINGZ(Mout), FINISHINGY(Mout), FINISHINGX(Mout));
    CenterFFT(Fin, false);
    local_transformer_out.setFromCompleteFourier(Fin);
    local_transformer_out.inverseFourierTransform();

    //#define DEBUG_RESCALE_VOLUME
#ifdef DEBUG_RESCALE_VOLUME

    Image<double> Vt;
    Vt()=Min;
    Vt.write("Min.vol");
    Vt()=Mout;
    Vt.write("Mout.vol");
#endif

    Min = Mout;

}

void
ProgMLTomo::postProcessVolume(Image<double> &Vin, double resolution)
{

    MultidimArray<std::complex<double> > Faux;
#ifdef DEBUG

    std::cerr<<"start postProcessVolume"<<std::endl;
#endif

    // Fourier transform
    transformer.FourierTransform(Vin(), Faux, true);

    // Apply fourier mask
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
    {
        DIRECT_MULTIDIM_ELEM(Faux,n) *= DIRECT_MULTIDIM_ELEM(fourier_mask,n);
    }

    // Low-pass filter for given resolution
    if (resolution > 0.)
    {
        Matrix1D<double> w(3);
        double raised_w = 0.02;
        double w1 = resolution;
        double w2 = w1 + raised_w;

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Faux)
        {
            FFT_IDX2DIGFREQ(k, ZSIZE(Vin()), ZZ(w));
            FFT_IDX2DIGFREQ(i, YSIZE(Vin()), YY(w));
            FFT_IDX2DIGFREQ(j, XSIZE(Vin()), XX(w));
            double absw = w.module();
            if (absw > w2)
                DIRECT_A3D_ELEM(Faux,k,i,j) *= 0.;
            else if (absw > w1 && absw < w2)
                DIRECT_A3D_ELEM(Faux,k,i,j) *= (1 + cos(PI / raised_w * (absw - w1)))
                                               / 2;
        }
    }

    // Inverse Fourier Transform
    transformer.inverseFourierTransform(Faux, Vin());

    // Spherical mask with average outside density
    maskSphericalAverageOutside(Vin());

    // Apply symmetry
    if (mysampling.SL.symsNo() > 0)
    {
        // Note that for no-imputation this is not correct!
        // One would have to symmetrize the missing wedges and the sum of the images separately
        if (do_missing && !do_impute)
            std::cerr
            << " WARNING: Symmetrization and dont_impute together are not implemented correctly...\n Proceed at your own risk"
            << std::endl;

        Image<double> Vaux, Vsym = Vin;
        Matrix2D<double> L(4, 4), R(4, 4);
        Matrix1D<double> sh(3);
        Vaux().initZeros(dim, dim, dim);
        Vaux().setXmippOrigin();
        for (int isym = 0; isym < mysampling.SL.symsNo(); isym++)
        {
            mysampling.SL.getMatrices(isym, L, R);
            mysampling.SL.getShift(isym, sh);
            R(3, 0) = sh(0) * dim;
            R(3, 1) = sh(1) * dim;
            R(3, 2) = sh(2) * dim;
            applyGeometry(LINEAR, Vaux(), Vin(), R.transpose(), IS_NOT_INV, DONT_WRAP,
                          DIRECT_MULTIDIM_ELEM(Vin(), 0));
            Vsym() += Vaux();
        }
        Vsym() /= mysampling.SL.symsNo() + 1.;
        Vin = Vsym;
    }

#ifdef DEBUG
    std::cerr<<"finished postProcessVolume"<<std::endl;
#endif

}

// Calculate FT of each reference and calculate A2 =============
void
ProgMLTomo::precalculateA2(std::vector<Image<double> > &Iref)
{

#ifdef DEBUG
    std::cerr<<"start precalculateA2"<<std::endl;
    TimeStamp t0;
    time_config();
    annotate_time(&t0);
#endif
#ifdef DEBUG_JM

    std::cerr << "DEBUG_JM: entering ProgMLTomo::precalculateA2" <<std::endl;
#endif

    double AA, stdAA, corr;
    Matrix2D<double> A_rot_inv(4, 4), I(4, 4);
    MultidimArray<double> Maux(dim, dim, dim);
    MultidimArray<unsigned char> Mmissing;
    MultidimArray<std::complex<double> > Faux, Faux2;

    A2.clear();
    corrA2.clear();
    I.initIdentity();
    Maux.setXmippOrigin();
    for (int refno = 0; refno < nr_ref; refno++)
    {
        MultidimArray<double> &Iref_refno=Iref[refno]();
        // Calculate A2 for all different orientations
        for (int angno = 0; angno < nr_ang; angno++)
        {
            A_rot_inv = ((all_angle_info[angno]).A).inv();
            // use DONT_WRAP and put density of first element outside
            // i.e. assume volume has been processed with omask
            applyGeometry(LINEAR, Maux, Iref_refno, A_rot_inv, IS_NOT_INV,
                          DONT_WRAP, DIRECT_MULTIDIM_ELEM(Iref_refno,0));
            //#define DEBUG_PRECALC_A2_ROTATE
#ifdef DEBUG_PRECALC_A2_ROTATE

            Image<double> Vt;
            Vt()=Maux;
            Vt.write("rot.vol");
            Vt()=Iref[refno]();
            Vt.write("ref.vol");
            std::cerr<<" angles= "<<(all_angle_info[angno]).rot<<" "<<(all_angle_info[angno]).tilt<<" "<<(all_angle_info[angno]).psi<<std::endl;
            std::cerr<<" Written volume rot.vol and ref.vol, press any key to continue ..."<<std::endl;
            char c;
            std::cin >> c;
#endif

            AA = Maux.sum2();
            if (angno == 0)
                stdAA = AA;
            if (AA > 0)
            {
                corr = sqrt(stdAA / AA);
                Maux *= corr;
            }
            else
                corr = 1.;
            corrA2.push_back(corr);
            if (do_missing)
            {
                transformer.FourierTransform(Maux, Faux, false);
                // Save original copy of Faux in Faux2
                Faux2=Faux;

                for (int missno = 0; missno < nr_miss; ++missno)
                {
                    getMissingRegion(Mmissing, I, missno);
                    Faux.initZeros();
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                    if (DIRECT_MULTIDIM_ELEM(Mmissing,n))
                        DIRECT_MULTIDIM_ELEM(Faux,n) = DIRECT_MULTIDIM_ELEM(Faux2,n);
                    transformer.inverseFourierTransform();
                    A2.push_back(Maux.sum2());
                    //#define DEBUG_PRECALC_A2
#ifdef DEBUG_PRECALC_A2

                    std::cerr<<"rot= "<<all_angle_info[angno].rot<<" tilt= "<<all_angle_info[angno].tilt<<" psi= "<<all_angle_info[angno].psi<<std::endl;
                    std::cerr<<"refno= "<<refno<<" angno= "<<angno<<" missno= "<<missno<<" A2= "<<Maux.sum2()<<" corrA2= "<<corr<<std::endl;
                    //#define DEBUG_PRECALC_A2_b
#ifdef DEBUG_PRECALC_A2_b

                    Image<double> tt;
                    tt()=Maux;
                    tt.write("refrotwedge.vol");
                    std::cerr<<"press any key"<<std::endl;
                    char c;
                    std::cin >> c;
#endif
#endif

                }
            }
            else
            {
                A2.push_back(Maux.sum2());
#ifdef DEBUG_PRECALC_A2

                std::cerr<<"refno= "<<refno<<" angno= "<<angno<<" A2= "<<Maux.sum2()<<std::endl;
#endif

            }
        }
    }

#ifdef DEBUG
    std::cerr<<"finished precalculateA2"<<std::endl;
    print_elapsed_time(t0);
#endif
}

// Maximum Likelihood calculation for one image ============================================
// Integration over all translation, given  model and in-plane rotation
void
ProgMLTomo::expectationSingleImage(MultidimArray<double> &Mimg, int imgno,
                                   const int missno, double old_psi, std::vector<Image<double> > &Iref,
                                   std::vector<MultidimArray<double> > &wsumimgs,
                                   std::vector<MultidimArray<double> > &wsumweds, double &wsum_sigma_noise,
                                   double &wsum_sigma_offset, MultidimArray<double> &sumw, double &LL,
                                   double &dLL, double &fracweight, double &sumfracweight, double &trymindiff,
                                   int &opt_refno, int &opt_angno, Matrix1D<double> &opt_offsets)
{

#ifdef DEBUG
    std::cerr<<"start expectationSingleImage"<<std::endl;
    TimeStamp t0, t1;
    time_config();
    annotate_time(&t0);
    annotate_time(&t1);
#endif

    //#define DEBUG_JM
#ifdef DEBUG_JM

    std::cerr << "DEBUG_JM: entering ProgMLTomo::expectationSingleImage" <<std::endl;
    std::cerr << "   DEBUG_JM: imgno: " << imgno << std::endl;
    std::cerr << "DEBUG_JM: missno: " << missno << std::endl;
#endif

    MultidimArray<double> Maux, Maux2, Mweight, Mzero(dim, dim,
            hdim + 1), Mzero2(dim, dim, dim);
    MultidimArray<unsigned char> Mmissing;
    MultidimArray<std::complex<double> > Faux, Faux2(dim, dim, hdim + 1), Fimg,
    Fimg0, Fimg_rot;
    std::complex<double> complex_zero=0;
    std::vector<MultidimArray<double> > mysumimgs;
    std::vector<MultidimArray<double> > mysumweds;
    MultidimArray<double> refw;
    double sigma_noise2, aux, pdf, fracpdf, myA2, mycorrAA, myXi2, A2_plus_Xi2,
    myXA;
    double diff, mindiff, my_mindiff;
    double my_sumweight, weight;
    double wsum_sc, wsum_sc2, wsum_offset;
    double wsum_corr, sum_refw, maxweight, my_maxweight;
    int sigdim;
    int ioptx = 0, iopty = 0, ioptz = 0;
    bool is_ok_trymindiff = false;
    int my_nr_ang, old_optangno = opt_angno, old_optrefno = opt_refno;
    std::vector<double> all_Xi2;
    Matrix2D<double> A_rot(4, 4), I(4, 4), A_rot_inv(4, 4);
    bool is_a_neighbor, is_within_psirange = true;
    FourierTransformer local_transformer;

    my_nr_ang = (dont_align || dont_rotate) ? 1 : nr_ang;

    // Only translations smaller than 6 sigma_offset are considered!
    // TODO: perhaps 3 sigma??
    I.initIdentity();
    sigdim = 2 * CEIL(sigma_offset * 6);
    sigdim++; // (to get uneven number)
    sigdim = XMIPP_MIN(dim, sigdim);
    // Setup matrices and constants
    Maux.resize(dim, dim, dim);
    Maux2.resize(dim, dim, dim);
    Maux.setXmippOrigin();
    Maux2.setXmippOrigin();
    if (dont_align)
        Mweight.initZeros(1, 1, 1);
    else
        Mweight.initZeros(sigdim, sigdim, sigdim);
    Mweight.setXmippOrigin();
    Mzero.initZeros();
    Mzero2.initZeros();
    Mzero2.setXmippOrigin();

    sigma_noise2 = sigma_noise * sigma_noise;

    Maux = Mimg;
    // Apply inverse rotation to the mask and apply
    if (do_mask)
    {
        if (!dont_align)
            REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,
                         "BUG: !dont_align and do_mask cannot coincide at this stage...");
        MultidimArray<double> Mmask;
        A_rot = (all_angle_info[opt_angno]).A;
        A_rot_inv = A_rot.inv();
        applyGeometry(LINEAR, Mmask, Imask(), A_rot_inv, IS_NOT_INV, DONT_WRAP, 0.);
        Maux *= Mmask;
    }
    // Calculate the unrotated Fourier transform with enforced wedge of Mimg (store in Fimg0)
    // also calculate myXi2;
    // Note that from here on local_transformer will act on Maux <-> Faux
    local_transformer.FourierTransform(Maux, Faux, false);
    if (do_missing)
    {
        // Enforce missing wedge
        getMissingRegion(Mmissing, I, missno);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
        if (!DIRECT_MULTIDIM_ELEM(Mmissing,n))
            DIRECT_MULTIDIM_ELEM(Faux,n) = complex_zero;
    }
    Fimg0 = Faux;
    // BE CAREFUL: inverseFourierTransform messes up Faux!
    local_transformer.inverseFourierTransform();
    myXi2 = Maux.sum2();

    // To avoid numerical problems, subtract smallest difference from all differences.
    // That way: Pmax will be one and all other probabilities will be [0,1>
    // But to find mindiff I first have to loop over all hidden variables...
    // Fortunately, there is some flexibility as the reliable domain of the exp-functions goes from -700 to 700.
    // Therefore, make a first guess for mindiff (trymindiff) and check whether the difference
    // with the real mindiff is not larger than 500 (to be on the save side).
    // If that is the case: OK; if not: do the entire loop again.
    // The efficiency of this will depend on trymindiff.
    // First cycle use trymindiff = trymindiff_factor * 0.5 * Xi2
    // From then on: use mindiff from the previous iteration
    if (trymindiff < 0.)
        // 90% of Xi2 may be a good idea (factor half because 0.5*diff is calculated)
        trymindiff = trymindiff_factor * 0.5 * myXi2;
    int redo_counter = 0;
    while (!is_ok_trymindiff)
    {
        // Initialize mindiff, weighted sums and maxweights
        mindiff = 99.e99;
        wsum_corr = wsum_offset = wsum_sc = wsum_sc2 = 0.;
        maxweight = sum_refw = 0.;
        mysumimgs.assign(nr_ref, Mzero2);
        if (do_missing)
            mysumweds.assign(nr_ref, Mzero);
        refw.initZeros(nr_ref);
        // The real stuff: now loop over all orientations, references and translations
        // Start the loop over all refno at old_optangno (=opt_angno from the previous iteration).
        // This will speed-up things because we will find Pmax probably right away,
        // and this will make the if-statement that checks SIGNIFICANT_WEIGHT_LOW
        // effective right from the start
        for (int aa = old_optangno; aa < old_optangno + my_nr_ang; aa++)
        {
            int angno = aa;
            if (angno >= nr_ang)
                angno -= nr_ang;

            // See whether this image is in the neighborhood for this imgno
            if (ang_search > 0.)
            {
                is_a_neighbor = false;
                if (do_limit_psirange)
                {
                    // Only restrict psi range for directions away from the poles...
                    // Otherwise the is_a_neighbour construction may give errors:
                    // (-150, 0, 150) and (-100, 0, 150) would be considered as equal, since
                    // the distance between (-150,0) and (-100,0) is zero AND
                    // the distance between 150 and 150 is also zero...
                    if (ABS(realWRAP(all_angle_info[angno].tilt, -90., 90.))
                        < 1.1 * ang_search)
                        is_within_psirange = true;
                    else if ((do_sym
                              && ABS(realWRAP(old_psi - (all_angle_info[angno]).rot,-180.,180.))
                              <= ang_search)
                             || (!do_sym
                                 && ABS(realWRAP(old_psi - (all_angle_info[angno]).psi,-180.,180.))
                                 <= ang_search))
                        is_within_psirange = true;
                    else
                        is_within_psirange = false;
                }

                if (!do_limit_psirange || is_within_psirange)
                {
                    for (size_t i = 0; i < mysampling.my_neighbors[imgno].size(); i++)
                    {
                        if (mysampling.my_neighbors[imgno][i]
                            == (all_angle_info[angno]).direction)
                        {
                            is_a_neighbor = true;
                            break;
                        }
                    }
                }
            }
            else
            {
                is_a_neighbor = true;
            }

            // If it is in the neighborhood: proceed
            if (is_a_neighbor)
            {
                A_rot = (all_angle_info[angno]).A;
                A_rot_inv = A_rot.inv();

                // Start the loop over all refno at old_optrefno (=opt_refno from the previous iteration).
                // This will speed-up things because we will find Pmax probably right away,
                // and this will make the if-statement that checks SIGNIFICANT_WEIGHT_LOW
                // effective right from the start
                for (int rr = old_optrefno; rr < old_optrefno + nr_ref; rr++)
                {
                    int refno = rr;
                    if (refno >= nr_ref)
                        refno -= nr_ref;

                    fracpdf = alpha_k(refno) * (1. / nr_ang);
                    // Now (inverse) rotate the reference and calculate its Fourier transform
                    // Use DONT_WRAP and assume map has been omasked
                    applyGeometry(LINEAR, Maux2, Iref[refno](), A_rot_inv, IS_NOT_INV,
                                  DONT_WRAP, DIRECT_MULTIDIM_ELEM(Iref[refno](),0));
                    mycorrAA = corrA2[refno * nr_ang + angno];
                    Maux = Maux2 * mycorrAA;
                    local_transformer.FourierTransform();
                    if (do_missing)
                        myA2 = A2[refno * nr_ang * nr_miss + angno * nr_miss + missno];
                    else
                        myA2 = A2[refno * nr_ang + angno];
                    A2_plus_Xi2 = 0.5 * (myA2 + myXi2);
                    // A. Backward FFT to calculate weights in real-space
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Faux)
                    {
                        dAkij(Faux,k,i,j) = dAkij(Fimg0,k,i,j) * conj(dAkij(Faux,k,i,j));
                    }
                    local_transformer.inverseFourierTransform();
                    CenterFFT(Maux, true);

                    // B. Calculate weights for each pixel within sigdim (Mweight)
                    my_sumweight = my_maxweight = 0.;
                    FOR_ALL_ELEMENTS_IN_ARRAY3D(Mweight)
                    {

                        pdf = fracpdf * A3D_ELEM(P_phi, k, i, j);
                        // Sjors 5 aug 2010
                        // With -limit_trans pdf can actually be zero, and mindiff is irrelevant for those.
                        if (pdf > 0.)
                        {
                            myXA = A3D_ELEM(Maux, k, i, j) * ddim3;
                            diff = A2_plus_Xi2 - myXA;
                            mindiff = XMIPP_MIN(mindiff,diff);
                            //#define DEBUG_MINDIFF
#ifdef  DEBUG_MINDIFF

                            if (mindiff < 0)
                            {
                                std::cerr<<"k= "<<k<<" j= "<<j<<" i= "<<i<<std::endl;
                                std::cerr<<"xaux="<<STARTINGX(Maux)<<" xw= "<<STARTINGX(Mweight)<<std::endl;
                                std::cerr<<"yaux="<<STARTINGY(Maux)<<" yw= "<<STARTINGY(Mweight)<<std::endl;
                                std::cerr<<"zaux="<<STARTINGZ(Maux)<<" zw= "<<STARTINGZ(Mweight)<<std::endl;
                                std::cerr<<"diff= "<<diff<<" A2_plus_Xi"<<std::endl;
                                std::cerr<<" mycorrAA= "<<mycorrAA<<" "<<std::endl;
                                std::cerr<<"debug mindiff= " <<mindiff<<" trymindiff= "<<trymindiff<< std::endl;
                                std::cerr<<"A2_plus_Xi2= "<<A2_plus_Xi2<<" myA2= "<<myA2<<" myXi2= "<<myXi2<<std::endl;
                                std::cerr.flush();
                                Image<double> tt;
                                tt()=Maux;
                                tt.write("Maux.vol");
                                tt()=Mweight;
                                tt.write("Mweight.vol");
                                std::cerr<<"Ainv= "<<A_rot_inv<<std::endl;
                                std::cerr<<"A= "<<A_rot_inv.inv()<<std::endl;
                                exit(1);
                            }
#endif
                            // Normal distribution
                            aux = (diff - trymindiff) / sigma_noise2;
                            // next line because of numerical precision of exp-function
                            if (aux > 1000.)
                                weight = 0.;
                            else
                                weight = exp(-aux) * pdf;
                            A3D_ELEM(Mweight, k, i, j) = weight;
                            // Accumulate sum weights for this (my) matrix
                            my_sumweight += weight;
                            // calculate weighted sum of (X-A)^2 for sigma_noise update
                            wsum_corr += weight * diff;
                            // calculated weighted sum of offsets as well
                            wsum_offset += weight * A3D_ELEM(Mr2, k, i, j);
                            // keep track of optimal parameters
                            my_maxweight = XMIPP_MAX(my_maxweight, weight);
                            if (weight > maxweight)
                            {
                                maxweight = weight;
                                ioptz = k;
                                iopty = i;
                                ioptx = j;
                                opt_angno = angno;
                                opt_refno = refno;
                            }
                        } // close if (pdf > 0.)

                    } // close for over all elements in Mweight

                    // C. only for significant settings, store weighted sums
                    if (my_maxweight > SIGNIFICANT_WEIGHT_LOW * maxweight)
                    {
                        //#define DEBUG_EXP_A2
#ifdef DEBUG_EXP_A2
                        std::cout<<" imgno= "<<imgno<<" refno= "<<refno<<" angno= "<<angno<<" A2= "<<myA2<<" <<Xi2= "<<myXi2<<" my_maxweight= "<<my_maxweight<<std::endl;
#endif

                        sum_refw += my_sumweight;
                        refw(refno) += my_sumweight;
                        // Back from smaller Mweight to original size of Maux
                        Maux.initZeros();
                        FOR_ALL_ELEMENTS_IN_ARRAY3D(Mweight)
                        {
                            A3D_ELEM(Maux, k, i, j) = A3D_ELEM(Mweight, k, i, j);
                        }
                        // Use forward FFT in convolution theorem again
                        CenterFFT(Maux, false);
                        local_transformer.FourierTransform();
                        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Faux)
                        {
                            dAkij(Faux, k, i, j) = conj(dAkij(Faux,k,i,j))
                                                   * dAkij(Fimg0,k,i,j);
                        }
                        local_transformer.inverseFourierTransform();
                        maskSphericalAverageOutside(Maux);
                        selfApplyGeometry(LINEAR, Maux, A_rot, IS_NOT_INV, DONT_WRAP,
                                          DIRECT_MULTIDIM_ELEM(Maux,0));
                        if (do_missing)
                        {
                            // Store sum of wedges!
                            getMissingRegion(Mmissing, A_rot, missno);
                            MultidimArray<double> &mysumweds_refno=mysumweds[refno];
                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mmissing)
                            if (DIRECT_MULTIDIM_ELEM(Mmissing,n))
                                DIRECT_MULTIDIM_ELEM(mysumweds_refno,n) += my_sumweight;

                            // Again enforce missing region to avoid filling it with artifacts from the rotation
                            local_transformer.FourierTransform();
                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mmissing)
                            if (!DIRECT_MULTIDIM_ELEM(Mmissing,n))
                                DIRECT_MULTIDIM_ELEM(Faux,n)=complex_zero;
                            local_transformer.inverseFourierTransform();
                        }
                        // Store sum of rotated images
                        mysumimgs[refno] += Maux * ddim3;
                    } // close if SIGNIFICANT_WEIGHT_LOW
                } // close for refno
            } // close if is_a_neighbor
        } // close for angno

        // Now check whether our trymindiff was OK.
        // The limit of the exp-function lies around
        // exp(700)=1.01423e+304, exp(800)=inf; exp(-700) = 9.85968e-305; exp(-800) = 0
        // Use 500 to be on the save side?
        if (ABS((mindiff - trymindiff) / sigma_noise2) > 500.)
        {

            //#define DEBUG_REDOCOUNTER
#ifdef DEBUG_REDOCOUNTER
            std::cerr<<"repeating mindiff "<<redo_counter<<"th time"<<std::endl;
            std::cerr<<"trymindiff= "<<trymindiff<<" mindiff= "<<mindiff<<std::endl;
            std::cerr<<"diff= "<<ABS((mindiff - trymindiff) / sigma_noise2)<<std::endl;
#endif
            // Re-do whole calculation now with the real mindiff
            trymindiff = mindiff;
            redo_counter++;
            // Never re-do more than once!
            if (redo_counter > 1)
                REPORT_ERROR(ERR_DEBUG_IMPOSIBLE, "ml_tomo BUG, redo_counter > 1");
        }
        else
        {
            is_ok_trymindiff = true;
            my_mindiff = trymindiff;
            trymindiff = mindiff;
        }
    }

    fracweight = maxweight / sum_refw;
    wsum_sc /= sum_refw;
    wsum_sc2 /= sum_refw;

    // Calculate remaining optimal parameters

    XX(opt_offsets) = -(double) ioptx;
    YY(opt_offsets) = -(double) iopty;
    ZZ(opt_offsets) = -(double) ioptz;

    // From here on lock threads
    pthread_mutex_lock(&mltomo_weightedsum_update_mutex);

    // Update all global weighted sums after division by sum_refw
    wsum_sigma_noise += (2 * wsum_corr / sum_refw);
    wsum_sigma_offset += (wsum_offset / sum_refw);
    sumfracweight += fracweight;
    // Randomly choose 0 or 1 for FSC calculation
    int iran_fsc = ROUND(rnd_unif());
    for (int refno = 0; refno < nr_ref; refno++)
    {
        sumw(refno) += refw(refno) / sum_refw;
        // Sum mysumimgs to the global weighted sum
        wsumimgs[iran_fsc * nr_ref + refno] += (mysumimgs[refno]) / sum_refw;
        if (do_missing)
        {
            wsumweds[iran_fsc * nr_ref + refno] += mysumweds[refno] / sum_refw;
        }
    }

    // 1st term: log(refw_i)
    // 2nd term: for subtracting mindiff
    // 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
    // TODO: check this!!
    if (do_missing)

        dLL =
            log(sum_refw) - my_mindiff / sigma_noise2
            - all_missing_info[missno].nr_pixels
            * log(sqrt(2. * PI * sigma_noise2));
    else
        dLL = log(sum_refw) - my_mindiff / sigma_noise2
              - ddim3 * log(sqrt(2. * PI * sigma_noise2));
    LL += dLL;

    pthread_mutex_unlock(&mltomo_weightedsum_update_mutex);

#ifdef DEBUG_JM

    std::cerr << "   DEBUG_JM: my_mindiff: " << my_mindiff << std::endl;
    std::cerr << "   DEBUG_JM: sigma_noise2: " << sigma_noise2 << std::endl;
    //std::cerr << "   DEBUG_JM: miss_nr_pixels[missno]: " << miss_nr_pixels[missno] << std::endl;

    std::cerr << "   DEBUG_JM: sum_refw: " << sum_refw << std::endl;
    std::cerr << "   DEBUG_JM: wsum_sigma_noise: " << wsum_sigma_noise << std::endl;
    std::cerr << "   DEBUG_JM: wsum_sigma_offset: " << wsum_sigma_offset << std::endl;
    std::cerr << "   DEBUG_JM: sumfracweight: " << sumfracweight << std::endl;
    std::cerr << "   DEBUG_JM: dLL: " << dLL << std::endl;
    std::cerr << "   DEBUG_JM: LL: " << LL << std::endl;
    std::cerr << "DEBUG_JM: leaving ProgMLTomo::expectationSingleImage" <<std::endl;
#endif
#undef DEBUG_JM

#ifdef DEBUG

    std::cerr<<"finished expectationSingleImage"<<std::endl;
    print_elapsed_time(t0);
#endif
}

void
ProgMLTomo::maxConstrainedCorrSingleImage(MultidimArray<double> &Mimg,
        int imgno, int missno, double old_psi, std::vector<Image<double> > &Iref,
        std::vector<MultidimArray<double> > &wsumimgs,
        std::vector<MultidimArray<double> > &wsumweds, MultidimArray<double> &sumw,
        double &maxCC, double &sumCC, int &opt_refno, int &opt_angno,
        Matrix1D<double> &opt_offsets)
{
    //#define DEBUG
#ifdef DEBUG
    std::cerr<<"start maxConstrainedCorrSingleImage"<<std::endl;
    TimeStamp t0, t1;
    time_config();
    annotate_time(&t0);
    annotate_time(&t1);
#endif
#ifdef DEBUG_JM

    std::cerr << "DEBUG_JM: entering ProgMLTomo::maxConstrainedCorrSingleImage" <<std::endl;
#endif

    MultidimArray<double> Mimg0, Maux, Mref;
    MultidimArray<unsigned char> Mmissing;
    MultidimArray<std::complex<double> > Faux, Fimg0, Fref;
    FourierTransformer local_transformer;
    Matrix2D<double> A_rot(4, 4), I(4, 4), A_rot_inv(4, 4);
    std::complex<double> complex_zero=0;
    bool is_a_neighbor;
    double img_stddev, ref_stddev, corr, maxcorr = -9999.;
    int ioptx, iopty, ioptz;
    int my_nr_ang, old_optangno = opt_angno, old_optrefno = opt_refno;
    bool is_within_psirange = true;

    my_nr_ang = (dont_align || dont_rotate) ? 1 : nr_ang;
    I.initIdentity();
    Maux.resize(dim, dim, dim);
    Maux.setXmippOrigin();

    Maux = Mimg;
    // Apply inverse rotation to the mask and apply
    if (do_mask)
    {
        if (!dont_align)
            REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,
                         "BUG: !dont_align and do_mask cannot coincide at this stage...");
        MultidimArray<double> Mmask;
        A_rot = (all_angle_info[opt_angno]).A;
        A_rot_inv = A_rot.inv();
        applyGeometry(LINEAR, Mmask, Imask(), A_rot_inv, IS_NOT_INV, DONT_WRAP, 0.);
        Maux *= Mmask;
    }
    // Calculate the unrotated Fourier transform with enforced wedge of Mimg (store in Fimg0)
    // Note that from here on local_transformer will act on Maux <-> Faux
    local_transformer.FourierTransform(Maux, Faux, false);
    if (do_missing)
    {
        // Enforce missing wedge
        getMissingRegion(Mmissing, I, missno);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mmissing)
        if (!DIRECT_MULTIDIM_ELEM(Mmissing,n))
            DIRECT_MULTIDIM_ELEM(Faux,n)=complex_zero;
        // BE CAREFUL: inverseFourierTransform messes up Faux!!
        Fimg0 = Faux;
        local_transformer.inverseFourierTransform();
    }
    else
        Fimg0 = Faux;
    Mimg0 = Maux;

    if (do_only_average)
    {
        maxcorr = 1.;
        ioptz = 0;
        iopty = 0;
        ioptx = 0;
        opt_angno = old_optangno;
        opt_refno = old_optrefno;
    }
    else
    {
        // Calculate stddev of (wedge-inforced) image
        img_stddev = Mimg0.computeStddev();
        // Loop over all orientations
        for (int aa = old_optangno; aa < old_optangno + my_nr_ang; aa++)
        {
            int angno = aa;
            if (angno >= nr_ang)
                angno -= nr_ang;

            // See whether this image is in the neighborhoood for this imgno
            if (ang_search > 0.)
            {
                is_a_neighbor = false;
                if (do_limit_psirange)
                {
                    // Only restrict psi range for directions away from the poles...
                    // Otherwise the is_a_neighbour construction may give errors:
                    // (-150, 0, 150) and (-100, 0, 150) would be considered as equal, since
                    // the distance between (-150,0) and (-100,0) is zero AND
                    // the distance between 150 and 150 is also zero...
                    if (ABS(realWRAP(all_angle_info[angno].tilt, -90., 90.))
                        < 1.1 * ang_search)
                        is_within_psirange = true;
                    else if ((do_sym
                              && ABS(realWRAP(old_psi - (all_angle_info[angno]).rot,-180.,180.))
                              <= ang_search)
                             || (!do_sym
                                 && ABS(realWRAP(old_psi - (all_angle_info[angno]).psi,-180.,180.))
                                 <= ang_search))
                        is_within_psirange = true;
                    else
                        is_within_psirange = false;
                }

                if (!do_limit_psirange || is_within_psirange)
                {
                    for (size_t i = 0; i < mysampling.my_neighbors[imgno].size(); i++)
                    {
                        if (mysampling.my_neighbors[imgno][i]
                            == (all_angle_info[angno]).direction)
                        {
                            is_a_neighbor = true;
                            break;
                        }
                    }
                }
            }
            else
            {
                is_a_neighbor = true;
            }
            // If it is in the neighborhoood: proceed
            if (is_a_neighbor)
            {
                A_rot = (all_angle_info[angno]).A;
                A_rot_inv = A_rot.inv();
                //        std::cerr << "A_rot" << A_rot << std::endl;
                //        std::cerr << "all_angle_info[angno].rot" << all_angle_info[angno].rot << std::endl;
                //        std::cerr << "all_angle_info[angno].tilt" << all_angle_info[angno].tilt << std::endl;
                //        std::cerr << "all_angle_info[angno].psi" << all_angle_info[angno].psi << std::endl;

                // Loop over all references
                for (int rr = old_optrefno; rr < old_optrefno + nr_ref; rr++)
                {
                    int refno = rr;
                    if (refno >= nr_ref)
                        refno -= nr_ref;

                    // Now (inverse) rotate the reference and calculate its Fourier transform
                    // Use DONT_WRAP because the reference has been omasked
                    applyGeometry(LINEAR, Maux, Iref[refno](), A_rot_inv, IS_NOT_INV,
                                  DONT_WRAP, DIRECT_MULTIDIM_ELEM(Iref[refno](),0));
                    local_transformer.FourierTransform();
                    if (do_missing)
                    {
                        // Enforce wedge on the reference
                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                        {
                            DIRECT_MULTIDIM_ELEM(Faux,n) *= DIRECT_MULTIDIM_ELEM(Mmissing,n);
                        }
                        // BE CAREFUL! inverseFourierTransform messes up Faux
                        Fref = Faux;
                        local_transformer.inverseFourierTransform();
                    }
                    else
                        Fref = Faux;
                    // Calculate stddev of (wedge-inforced) reference
                    Mref = Maux;
                    ref_stddev = Mref.computeStddev();

                    // Calculate correlation matrix via backward FFT
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Faux)
                    {
                        dAkij(Faux,k,i,j) = dAkij(Fimg0,k,i,j) * conj(dAkij(Fref,k,i,j));
                    }
                    local_transformer.inverseFourierTransform();
                    CenterFFT(Maux, true);
                    //#define DEBUG_IMG
#ifdef DEBUG_IMG

                    static size_t jjjj=1;
                    Maux.write("kk.spi");
                    Image<double> tmpImg;
                    tmpImg() = Maux;
                    String s;
                    formatStringFast(s, "%d@%s.%s", jjjj++,"Maux","spi");
                    std::cerr << s <<std::endl;
                    tmpImg.write(s);

#endif

                    if (dont_align)
                    {
                        corr = A3D_ELEM(Maux,0,0,0) / (img_stddev * ref_stddev);
                        if (corr > maxcorr)
                        {
                            maxcorr = corr;
                            ioptz = 0;
                            iopty = 0;
                            ioptx = 0;
                            opt_angno = angno;
                            opt_refno = refno;
                        }
                    }
                    else
                    {
                        FOR_ALL_ELEMENTS_IN_ARRAY3D(Maux)
                        {
                            corr = A3D_ELEM(Maux,k,i,j) / (img_stddev * ref_stddev);
                            if (corr > maxcorr)
                            {
                                maxcorr = corr;
                                ioptz = k;
                                iopty = i;
                                ioptx = j;
                                opt_angno = angno;
                                opt_refno = refno;
                            }
                        }
                    }
                }
            }
        }
    }

    // Now that we have optimal hidden parameters, update output

    XX(opt_offsets) = -(double) ioptx;
    YY(opt_offsets) = -(double) iopty;
    ZZ(opt_offsets) = -(double) ioptz;
    selfTranslate(LINEAR, Mimg0, opt_offsets, DONT_WRAP);
    A_rot = (all_angle_info[opt_angno]).A;
    maskSphericalAverageOutside(Mimg0);
    selfApplyGeometry(LINEAR, Mimg0, A_rot, IS_NOT_INV, DONT_WRAP,
                      DIRECT_MULTIDIM_ELEM(Mimg0,0));
    maxCC = maxcorr;

    if (do_missing)
    {
        // Store sum of wedges
        getMissingRegion(Mmissing, A_rot, missno);
        Maux = Mimg0;
        // Again enforce missing region to avoid filling it with artifacts from the rotation
        local_transformer.FourierTransform();
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
        {
            DIRECT_MULTIDIM_ELEM(Faux,n) *= DIRECT_MULTIDIM_ELEM(Mmissing,n);
        }
        local_transformer.inverseFourierTransform();
        Mimg0 = Maux;
    }

    // Randomly choose 0 or 1 for FSC calculation
    int iran_fsc = ROUND(rnd_unif());

    // From here on lock threads to add to sums
    pthread_mutex_lock(&mltomo_weightedsum_update_mutex);

    sumCC += maxCC;
    sumw(opt_refno) += 1.;
    wsumimgs[iran_fsc * nr_ref + opt_refno] += Mimg0;
    if (do_missing)
    {
        MultidimArray<double> &wsumweds_i=wsumweds[iran_fsc * nr_ref + opt_refno];
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mmissing)
        DIRECT_MULTIDIM_ELEM(wsumweds_i,n)=DIRECT_MULTIDIM_ELEM(Mmissing,n);
    }

    pthread_mutex_unlock(&mltomo_weightedsum_update_mutex);

#ifdef DEBUG

    std::cerr<<"finished maxConstrainedCorrSingleImage"<<std::endl;
    print_elapsed_time(t0);
#endif

}

void *
threadMLTomoExpectationSingleImage(void * data)
{
    structThreadExpectationSingleImage * thread_data =
        (structThreadExpectationSingleImage *) data;

    // Variables from above
    int thread_id = thread_data->thread_id;
    ProgMLTomo *prm = thread_data->prm;
    MetaData *MDimg = thread_data->MDimg;
    double *wsum_sigma_noise = thread_data->wsum_sigma_noise;
    double *wsum_sigma_offset = thread_data->wsum_sigma_offset;
    double *sumfracweight = thread_data->sumfracweight;
    double *LL = thread_data->LL;
    std::vector<MultidimArray<double> > *wsumimgs = thread_data->wsumimgs;
    std::vector<MultidimArray<double> > *wsumweds = thread_data->wsumweds;
    std::vector<Image<double> > *Iref = thread_data->Iref;

    MultidimArray<double> *sumw = thread_data->sumw;
    std::vector<size_t> * imgs_id = thread_data->imgs_id;
    ThreadTaskDistributor * distributor = thread_data->distributor;

    //#define DEBUG_THREAD
#ifdef DEBUG_THREAD

    std::cerr<<"start threadMLTomoExpectationSingleImage"<<std::endl;
#endif

    // Local variables
    Image<double> img;
    FileName fn_img, fn_trans;
    std::vector<MultidimArray<double> > allref_offsets;
    Matrix1D<double> opt_offsets(3);
    float old_psi = -999.;
    double fracweight, trymindiff, dLL;
    int opt_refno, opt_angno, missno;

    //try
    {
        //Aproximate number of images
        size_t myNum = MDimg->size() / prm->threads;
        //First and last image to process in each block
        size_t firstIndex, lastIndex;

        if (prm->verbose && thread_id == 0)
            init_progress_bar(myNum);

        // Loop over all images
        size_t cc = 0;
        MultidimArray<double> &docfiledata = prm->docfiledata;

        //Work while there are tasks to do
        while (distributor->getTasks(firstIndex, lastIndex))
        {
            for (size_t imgno = firstIndex; imgno <= lastIndex; ++imgno)
            {
                //TODO: Check if really needed the mutexes
                //only for read from MetaData
                pthread_mutex_lock(&mltomo_selfile_access_mutex);

                MDimg->getValue(MDL_IMAGE, fn_img, (*imgs_id)[imgno + prm->myFirstImg]);

                pthread_mutex_unlock(&mltomo_selfile_access_mutex);

                img.read(fn_img);
                img().setXmippOrigin();

                if (prm->dont_align || prm->do_only_average)
                {
                    selfTranslate(LINEAR, img(), prm->imgs_optoffsets[imgno], DONT_WRAP);
                }
                prm->reScaleVolume(img(), true);
                missno = (prm->do_missing) ? prm->imgs_missno[imgno] : -1;
                // These three parameters speed up expectationSingleImage
                trymindiff = prm->imgs_trymindiff[imgno];
                opt_refno = prm->imgs_optrefno[imgno];
                opt_angno = prm->imgs_optangno[imgno];
                old_psi = prm->imgs_optpsi[imgno];

                if (prm->do_ml)
                {
                    // A. Use maximum likelihood approach
                    prm->expectationSingleImage(img(), imgno, missno, old_psi, *Iref,
                                                *wsumimgs, *wsumweds, *wsum_sigma_noise, *wsum_sigma_offset,
                                                *sumw, *LL, dLL, fracweight, *sumfracweight, trymindiff,
                                                opt_refno, opt_angno, opt_offsets);
                }
                else
                {
                    // B. Use constrained correlation coefficient approach
                    prm->maxConstrainedCorrSingleImage(img(), imgno, missno, old_psi,
                                                       *Iref, *wsumimgs, *wsumweds, *sumw, fracweight, *sumfracweight,
                                                       opt_refno, opt_angno, opt_offsets);
                }
                // Store for next iteration
                prm->imgs_trymindiff[imgno] = trymindiff;
                prm->imgs_optrefno[imgno] = opt_refno;
                prm->imgs_optangno[imgno] = opt_angno;
                // Output MetaData
                //FIXME: THIS ALSO LIKE JoseMiguel did ML2D
                dAij(docfiledata, imgno, 0) = (prm->all_angle_info[opt_angno]).rot; // rot
                dAij(docfiledata, imgno, 1) = (prm->all_angle_info[opt_angno]).tilt; // tilt
                dAij(docfiledata, imgno, 2) = (prm->all_angle_info[opt_angno]).psi; // psi
                if (prm->dont_align || prm->do_only_average)
                {
                    dAij(docfiledata, imgno, 3) = prm->imgs_optoffsets[imgno](0); // Xoff
                    dAij(docfiledata, imgno, 4) = prm->imgs_optoffsets[imgno](1); // Yoff
                    dAij(docfiledata, imgno, 5) = prm->imgs_optoffsets[imgno](2); // zoff
                }
                else
                {
                    dAij(docfiledata, imgno, 3) = opt_offsets(0) / prm->scale_factor; // Xoff
                    dAij(docfiledata, imgno, 4) = opt_offsets(1) / prm->scale_factor; // Yoff
                    dAij(docfiledata, imgno, 5) = opt_offsets(2) / prm->scale_factor; // Zoff
                }
                dAij(docfiledata, imgno, 6) = (double) (opt_refno + 1); // Ref
                dAij(docfiledata, imgno, 7) = (double) (missno + 1); // Wedge number

                dAij(docfiledata, imgno, 8) = fracweight; // P_max/P_tot
                dAij(docfiledata, imgno, 9) = dLL; // log-likelihood

                if (prm->verbose && thread_id == 0)
                    progress_bar(cc);
                cc++;
            }
        }
        if (prm->verbose && thread_id == 0)
            progress_bar(myNum);
    }

#ifdef DEBUG_THREAD

    std::cerr<<"finished threadMLTomoExpectationSingleImage"<<std::endl;
#endif
    return NULL;
}

void
ProgMLTomo::expectation(MetaData &MDimg, std::vector<Image<double> > &Iref,
                        int iter, double &LL, double &sumfracweight,
                        std::vector<MultidimArray<double> > &wsumimgs,
                        std::vector<MultidimArray<double> > &wsumweds, double &wsum_sigma_noise,
                        double &wsum_sigma_offset, MultidimArray<double> &sumw)
{
    //#define DEBUG
#ifdef DEBUG
    std::cerr<<"start expectation"<<std::endl;
#endif

    MultidimArray<double> Mzero(dim, dim, hdim + 1), Mzero2(dim, dim, dim);

    // Perturb all angles
    // (Note that each mpi process will have a different random perturbation)
    if (do_perturb)
        perturbAngularSampling();

    if (do_ml)
    {
        // Precalculate A2-values for all references
        precalculateA2(Iref);

        // Pre-calculate pdf of all in-plane transformations
        calculatePdfTranslations();
    }

    // Initialize weighted sums
    LL = 0.;
    wsum_sigma_noise = 0.;
    wsum_sigma_offset = 0.;
    sumfracweight = 0.;
    sumw.initZeros(nr_ref);
    //Create a task distributor to distribute images to process
    //each thread will process images one by one
    ThreadTaskDistributor * distributor = new ThreadTaskDistributor(
                                              nr_images_local, 1);

    Mzero.initZeros();
    Mzero2.initZeros();
    Mzero2.setXmippOrigin();
    wsumimgs.assign(2 * nr_ref, Mzero2);
    wsumweds.assign(2 * nr_ref, Mzero);
    // Call threads to calculate the expectation of each image in the selfile
    pthread_t * th_ids = (pthread_t *) malloc(threads * sizeof(pthread_t));
    structThreadExpectationSingleImage * threads_d =
        (structThreadExpectationSingleImage *) malloc(
            threads * sizeof(structThreadExpectationSingleImage));
    for (int c = 0; c < threads; c++)
    {
        threads_d[c].thread_id = c;
        threads_d[c].thread_num = threads;
        threads_d[c].prm = this;
        threads_d[c].MDimg = &MDimg;
        threads_d[c].iter = &iter;
        threads_d[c].wsum_sigma_noise = &wsum_sigma_noise;
        threads_d[c].wsum_sigma_offset = &wsum_sigma_offset;
        threads_d[c].sumfracweight = &sumfracweight;
        threads_d[c].LL = &LL;
        threads_d[c].wsumimgs = &wsumimgs;
        threads_d[c].wsumweds = &wsumweds;
        threads_d[c].Iref = &Iref;
        threads_d[c].sumw = &sumw;
        threads_d[c].imgs_id = &imgs_id;
        threads_d[c].distributor = distributor;
        pthread_create(
            (th_ids + c), NULL, threadMLTomoExpectationSingleImage, (void *)(threads_d+c) );
    }
    // Wait for threads to finish and get joined MetaData
    for (int c = 0; c < threads; c++)
    {
        pthread_join(*(th_ids + c), NULL);
    }
    //Free some memory
    delete distributor;

    //FIXME
    // Send back output in the form of a MetaData
    //use docfiledata in writingOutput files
    /*    SF.go_beginning();
     for (int imgno = 0; imgno < SF.ImgNo(); imgno++)
     {
     DFo.append_comment(SF.NextImg());
     DFo.append_data_line(docfiledata[imgno]);
     }*/

#ifdef DEBUG

    std::cerr<<"finished expectation"<<std::endl;
#endif
#undef DEBUG
}

// Update all model parameters
void
ProgMLTomo::maximization(std::vector<MultidimArray<double> > &wsumimgs,
                         std::vector<MultidimArray<double> > &wsumweds, double &wsum_sigma_noise,
                         double &wsum_sigma_offset, MultidimArray<double> &sumw,
                         double &sumfracweight, double &sumw_allrefs,
                         std::vector<MultidimArray<double> > &fsc, int iter)
{
#ifdef DEBUG
    std::cerr<<"started maximization"<<std::endl;
#endif

    MultidimArray<double> rmean_sigma2, rmean_signal2;
    Matrix1D<int> center(3);
    MultidimArray<int> radial_count;
    MultidimArray<std::complex<double> > Faux(dim, dim, hdim + 1), Fwsumimgs;
    MultidimArray<double> Maux, Msumallwedges(dim, dim, hdim + 1);
    std::vector<double> refs_resol(nr_ref);
    Image<double> Vaux;
    FileName fn_tmp;

    // Calculate resolutions
    fsc.clear();
    fsc.resize(nr_ref + 1);
    for (int refno = 0; refno < nr_ref; refno++)
    {
        calculateFsc(wsumimgs[refno], wsumimgs[nr_ref + refno], wsumweds[refno],
                     wsumweds[nr_ref + refno], fsc[0], fsc[refno + 1], refs_resol[refno]);
        // Restore original wsumimgs and wsumweds
        wsumimgs[refno] += wsumimgs[nr_ref + refno];
        wsumweds[refno] += wsumweds[nr_ref + refno];
    }

    // Update the reference images
    sumw_allrefs = 0.;
    Msumallwedges.initZeros();
    for (int refno = 0; refno < nr_ref; refno++)
    {
        sumw_allrefs += sumw(refno);
        transformer.FourierTransform(wsumimgs[refno], Fwsumimgs, true);
        if (do_missing)
        {
            Msumallwedges += wsumweds[refno];
            // Only impute for ML-approach (maxCC approach uses division by sumwedge)
            if (do_ml && do_impute)
            {
                transformer.FourierTransform(Iref[refno](), Faux, true);
                if (sumw(refno) > 0.)
                {
                    // inner iteration: marginalize over missing data only
                    for (int iter2 = 0; iter2 < Niter2; iter2++)
                    {
                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                        {
                            // Impute old reference for missing pixels
                            DIRECT_MULTIDIM_ELEM(Faux,n) *= (1.
                                                             - DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) / sumw(refno));
                            // And sum the weighted sum for observed pixels
                            DIRECT_MULTIDIM_ELEM(Faux,n) += (DIRECT_MULTIDIM_ELEM(Fwsumimgs,n)
                                                             / sumw(refno));
                        }
                    }
                }
                // else do nothing (i.e. impute old reference completely
            }
            else
            {
                // no imputation: divide by number of times a pixel has been observed
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                {
                    if (DIRECT_MULTIDIM_ELEM(wsumweds[refno],n) > noimp_threshold)
                    {
                        DIRECT_MULTIDIM_ELEM(Faux,n) = DIRECT_MULTIDIM_ELEM(Fwsumimgs,n)
                                                       / DIRECT_MULTIDIM_ELEM(wsumweds[refno],n);
                        // TODO:  CHECK THE FOLLOWING LINE!!
                        DIRECT_MULTIDIM_ELEM(Faux,n) *= sumw(refno) / sumw(refno);
                    }
                    else
                    {
                        DIRECT_MULTIDIM_ELEM(Faux,n) = 0.;
                    }
                }
            }
        }
        else
        {
            // No missing data
            if (sumw(refno) > 0.)
            {
                // if no missing data: calculate straightforward average,
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
                {
                    DIRECT_MULTIDIM_ELEM(Faux,n) = DIRECT_MULTIDIM_ELEM(Fwsumimgs,n)
                                                   / sumw(refno);
                }
            }
        }
        // Back to real space
        transformer.inverseFourierTransform(Faux, Iref[refno]());
    }

    // Average fracweight
    sumfracweight /= sumw_allrefs;

    // Update the model fractions
    if (!fix_fractions)
    {
        for (int refno = 0; refno < nr_ref; refno++)
        {
            if (sumw(refno) > 0.)
                alpha_k(refno) = sumw(refno) / sumw_allrefs;
            else
                alpha_k(refno) = 0.;
        }
    }

    // Update sigma of the origin offsets
    if (!fix_sigma_offset)
    {
        sigma_offset = sqrt(wsum_sigma_offset / (2. * sumw_allrefs));
    }

    // Update the noise parameters
    if (!fix_sigma_noise)
    {
        if (do_missing)
        {
            double sum_complete_wedge = 0.;
            double sum_complete_fourier = 0.;
            MultidimArray<double> Mcomplete(dim, dim, dim);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Mcomplete)
            {
                if (j < XSIZE(Msumallwedges))
                {
                    sum_complete_wedge += DIRECT_A3D_ELEM(Msumallwedges,k,i,j);
                    sum_complete_fourier += DIRECT_A3D_ELEM(fourier_imask,k,i,j);
                }
                else
                {
                    sum_complete_wedge += DIRECT_A3D_ELEM(Msumallwedges,
                                                          (dim-k)%dim,(dim-i)%dim,dim-j);
                    sum_complete_fourier += DIRECT_A3D_ELEM(fourier_imask,
                                                            (dim-k)%dim,(dim-i)%dim,dim-j);
                }
            }
            //#define DEBUG_UPDATE_SIGMA
#ifdef DEBUG_UPDATE_SIGMA
            std::cerr<<" sum_complete_wedge= "<<sum_complete_wedge<<" = "<<100*sum_complete_wedge/(sumw_allrefs*sum_complete_fourier)<<"%"<<std::endl;
            std::cerr<<" ddim3= "<<ddim3<<std::endl;
            std::cerr<<" sum_complete_fourier= "<<sum_complete_fourier<<std::endl;
            std::cerr<<" sumw_allrefs= "<<sumw_allrefs<<std::endl;
            std::cerr<<" wsum_sigma_noise= "<<wsum_sigma_noise<<std::endl;
            std::cerr<<" sigma_noise_old= "<<sigma_noise<<std::endl;
            std::cerr<<" sigma_new_impute= "<<sqrt((sigma_noise*sigma_noise*(sumw_allrefs*sum_complete_fourier - sum_complete_wedge ) + wsum_sigma_noise)/(sumw_allrefs * sum_complete_fourier))<<std::endl;
            std::cerr<<" sigma_new_noimpute= "<<sqrt(wsum_sigma_noise / (sum_complete_wedge))<<std::endl;
            std::cerr<<" sigma_new_nomissing= "<<sqrt(wsum_sigma_noise / (sumw_allrefs * sum_complete_fourier))<<std::endl;
#endif

            if (do_impute)
            {
                for (int iter2 = 0; iter2 < Niter2; iter2++)
                {
                    sigma_noise *= sigma_noise;
                    sigma_noise *= sumw_allrefs * sum_complete_fourier
                                   - sum_complete_wedge;
                    sigma_noise += wsum_sigma_noise;
                    sigma_noise = sqrt(
                                      sigma_noise / (sumw_allrefs * sum_complete_fourier));
                }
            }
            else
            {
                sigma_noise = sqrt(wsum_sigma_noise / (sum_complete_wedge));
            }
        }
        else
        {
            sigma_noise = sqrt(wsum_sigma_noise / (sumw_allrefs * ddim3));
        }
        // Correct sigma_noise for number of pixels within the mask
        if (do_mask)
            sigma_noise *= ddim3 / nr_mask_pixels;
    }

    // post-process reference volumes
    for (int refno = 0; refno < nr_ref; refno++)
    {
        if (do_filter)
            postProcessVolume(Iref[refno], refs_resol[refno]);
        else
            postProcessVolume(Iref[refno]);
    }

    // Regularize
    regularize(iter);

#ifdef DEBUG

    std::cerr<<"finished maximization"<<std::endl;
#endif
}

void
ProgMLTomo::calculateFsc(MultidimArray<double> &M1, MultidimArray<double> &M2,
                         MultidimArray<double> &W1, MultidimArray<double> &W2,
                         MultidimArray<double> &freq, MultidimArray<double> &fsc, double &resolution)
{

    MultidimArray<std::complex<double> > FT1, FT2;
    FourierTransformer transformer1, transformer2;
    transformer1.FourierTransform(M1, FT1, false);
    transformer2.FourierTransform(M2, FT2, false);

    int vsize = hdim + 1;
    Matrix1D<double> f(3);
    MultidimArray<double> num, den1, den2;
    double w1, w2, R;
    freq.initZeros(vsize);
    fsc.initZeros(vsize);
    num.initZeros(vsize);
    den1.initZeros(vsize);
    den2.initZeros(vsize);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FT1)
    {
        FFT_IDX2DIGFREQ(j, XSIZE(M1), XX(f));
        FFT_IDX2DIGFREQ(i, YSIZE(M1), YY(f));
        FFT_IDX2DIGFREQ(k, ZSIZE(M1), ZZ(f));
        R = f.module();
        if (R > 0.5)
            continue;
        int idx = ROUND(R*XSIZE(M1));
        if (do_missing)
        {
            w1 = dAkij(W1, k, i, j);
            w2 = dAkij(W2, k, i, j);
        }
        else
        {
            w1 = w2 = 1.;
        }
        std::complex<double> z1 = w2 * dAkij(FT1, k, i, j);
        std::complex<double> z2 = w1 * dAkij(FT2, k, i, j);
        double absz1 = abs(z1);
        double absz2 = abs(z2);
        num(idx) += real(conj(z1) * (z2));
        den1(idx) += absz1 * absz1;
        den2(idx) += absz2 * absz2;
    }

    FOR_ALL_ELEMENTS_IN_ARRAY1D(freq)
    {
        freq(i) = (double) i / (XSIZE(M1) * pixel_size);
        if (num(i) != 0)
            fsc(i) = num(i) / sqrt(den1(i) * den2(i));
        else
            fsc(i) = 0;
    }

    // Find resolution limit acc. to FSC=0.5 criterion
    int idx;
    for (idx = 1; idx < XSIZE(freq); idx++)
        if (fsc(idx) < 0.5)
            break;
    idx = XMIPP_MAX(idx - 1, 1);
    resolution = freq(idx);

    //#define DEBUG_FSC
#ifdef DEBUG_FSC

    FOR_ALL_ELEMENTS_IN_ARRAY1D(freq)
    {
        std::cerr<<freq(i)<<" "<<fsc(i)<<std::endl;
    }
    std::cerr<<"resolution= "<<resolution<<std::endl;
#endif

}

// Apply regularization
void
ProgMLTomo::regularize(int iter)
{

    // Update regularization constant in a linear manner
    reg_current = reg0 - (double) iter * (reg0 - regF) / (double) reg_steps;
    reg_current = XMIPP_MAX(reg_current, regF);

    if (reg_current > 0.)
    {
#ifdef DEBUG
        std::cerr<<"start regularize"<<std::endl;
#endif
        // Write out original volumes (before regularization)
        FileName fnt;
        for (int refno = 0; refno < nr_ref; refno++)
        {
            fnt = formatString("%s_it%06d_oriref%06d.vol", fn_root.c_str(), iter,
                               refno + 1);
            Iref[refno].write(fnt);
        }
        // Normalized regularization (in N/K)
        double reg_norm = reg_current * (double) nr_images_global / (double) nr_ref;
        MultidimArray<std::complex<double> > Fref, Favg, Fzero(dim, dim, hdim + 1);
        MultidimArray<double> Mavg, Mdiff;
        double sum_diff2 = 0.;

        //#define DEBUG_REGULARISE
#ifdef DEBUG_REGULARISE

        std::cerr<<"written out oriref volumes"<<std::endl;
#endif
        // Calculate  FT of average of all references
        // and sum of squared differences between all references
        for (int refno = 0; refno < nr_ref; refno++)
        {
            if (refno == 0)
                Mavg = Iref[refno]();
            else
                Mavg += Iref[refno]();
            for (int refno2 = 0; refno2 < nr_ref; refno2++)
            {
                Mdiff = Iref[refno]() - Iref[refno2]();
                sum_diff2 += Mdiff.sum2();
            }
        }
        transformer.FourierTransform(Mavg, Favg, true);
        Favg *= std::complex<double>(reg_norm, 0.);
#ifdef DEBUG_REGULARISE

        std::cerr<<"calculated Favg"<<std::endl;
#endif

        // Update the regularized references
        for (int refno = 0; refno < nr_ref; refno++)
        {
            transformer.FourierTransform(Iref[refno](), Fref, true);
            double sumw = alpha_k(refno) * (double) nr_images_global;
            // Fref = (sumw*Fref + reg_norm*Favg) /  (sumw + nr_ref * reg_norm)
#ifdef DEBUG_REGULARISE

            if (verb>0)
                std::cerr<<"refno= "<<refno<<" sumw = "<<sumw<<" reg_current= "<<reg_current<<" reg_norm= "<<reg_norm<<" Fref1= "<<DIRECT_MULTIDIM_ELEM(Fref,1) <<" Favg1= "<<DIRECT_MULTIDIM_ELEM(Favg,1)<<" (sumw + nr_ref * reg_norm)= "<<(sumw + nr_ref * reg_norm)<<std::endl;
#endif

            Fref *= std::complex<double>(sumw, 0.);
            Fref += Favg;
            Fref /= std::complex<double>((sumw + nr_ref * reg_norm), 0.);
#ifdef DEBUG_REGULARISE

            if (verb>0)
                std::cerr<<" newFref1= "<<DIRECT_MULTIDIM_ELEM(Fref,1) <<std::endl;
#endif

            transformer.inverseFourierTransform(Fref, Iref[refno]());
        }

        // Update the regularized sigma_noise estimate
        if (!fix_sigma_noise)
        {
            double reg_sigma_noise2 = sigma_noise * sigma_noise * ddim3
                                      * (double) nr_images_global;
#ifdef DEBUG_REGULARISE

            if (verb>0)
                std::cerr<<"reg_sigma_noise2= "<<reg_sigma_noise2<<" nr_exp_images="<<nr_images_global<<" ddim3= "<<ddim3<<" sigma_noise= "<<sigma_noise<<" sum_diff2= "<<sum_diff2<<" reg_norm= "<<reg_norm<<std::endl;
#endif

            reg_sigma_noise2 += reg_norm * sum_diff2;
            sigma_noise = sqrt(
                              reg_sigma_noise2 / ((double) nr_images_global * ddim3));
#ifdef DEBUG_REGULARISE

            if (verb>0)
                std::cerr<<"new sigma_noise= "<<sigma_noise<<std::endl;
#endif

        }

#ifdef DEBUG
        std::cerr<<"finished regularize"<<std::endl;
#endif

    }
}
// Check convergence
bool
ProgMLTomo::checkConvergence(std::vector<double> &conv)
{
#ifdef DEBUG
    std::cerr<<"started checkConvergence"<<std::endl;
#endif

    bool converged = true;
    double convv;
    MultidimArray<double> Maux;

    Maux.resize(dim, dim, dim);
    Maux.setXmippOrigin();

    conv.clear();
    for (int refno = 0; refno < nr_ref; refno++)
    {
        if (alpha_k(refno) > 0.)
        {
            Maux = Iold[refno]() * Iold[refno]();
            convv = 1. / (Maux.computeAvg());
            Maux = Iold[refno]() - Iref[refno]();
            Maux = Maux * Maux;
            convv *= Maux.computeAvg();
            conv.push_back(convv);
            if (convv > eps)
                converged = false;
        }
        else
        {
            conv.push_back(-1.);
        }
    }

#ifdef DEBUG
    std::cerr<<"finished checkConvergence"<<std::endl;
#endif
#undef DEBUG

    return converged;
}

void
ProgMLTomo::addPartialDocfileData(const MultidimArray<double> &data,
                                  size_t first, size_t last)
{
    size_t index;
    MDRow row;

    for (size_t imgno = first; imgno <= last; ++imgno)
    {
        index = imgno - first ;
        MDimg.getRow(row, imgs_id[imgno]);
        row.setValue(MDL_ANGLE_ROT, dAij(data, index, 0));
        row.setValue(MDL_ANGLE_TILT, dAij(data, index, 1));
        row.setValue(MDL_ANGLE_PSI, dAij(data, index, 2));
        row.setValue(MDL_SHIFT_X, dAij(data, index, 3));
        row.setValue(MDL_SHIFT_Y, dAij(data, index, 4));
        row.setValue(MDL_SHIFT_Z, dAij(data, index, 5));
        row.setValue(MDL_REF, ROUND(dAij(data, index, 6)));
        if (do_missing)
            row.setValue(MDL_REF, ROUND(dAij(data, index, 7)));
        MDimg.setRow(row, imgs_id[imgno]);
    }
}

void
ProgMLTomo::writeOutputFiles(const int iter,
                             std::vector<MultidimArray<double> > &wsumweds, double &sumw_allrefs,
                             double &LL, double &avefracweight, std::vector<double> &conv,
                             std::vector<MultidimArray<double> > &fsc)
{

    FileName fn_tmp, fn_base, fn_tmp2;
    MultidimArray<double> fracline(2);
    MetaData SFo, SFc;
    MetaData DFl;
    MetaData MDo, MDSel;
    std::string comment;
    Image<double> Vt;

    DFl.clear();
    SFo.clear();
    SFc.clear();

    fn_base = fn_root;
    if (iter >= 0)
    {
        fn_base = formatString("%s_it%06d", fn_root.c_str(), iter);
    }

    // Write out current reference images and fill sel & log-file
    //Re-use the MDref that were read in produceSideInfo2
    int refno = 0;
    //for (int refno = 0; refno < nr_ref; refno++)
    FOR_ALL_OBJECTS_IN_METADATA(MDref)
    {
        fn_tmp = formatString("%s_ref%06d.vol", fn_base.c_str(), refno + 1);
        Vt = Iref[refno];
        reScaleVolume(Vt(), false);
        Vt.write(fn_tmp);

        MDref.setValue(MDL_IMAGE, fn_tmp, __iter.objId);
        MDref.setValue(MDL_ENABLED, 1, __iter.objId);
        //SFo.insert(fn_tmp, SelLine::ACTIVE);
        //fracline(0) = alpha_k(refno);
        //fracline(1) = 1000 * conv[refno]; // Output 1000x the change for precision
        MDref.setValue(MDL_WEIGHT, alpha_k(refno), __iter.objId);
        MDref.setValue(MDL_SIGNALCHANGE, conv[refno] * 1000, __iter.objId);
        //DFl.insert_comment(fn_tmp);
        //DFl.insert_data_line(fracline);

        if (iter >= 1 && do_missing)
        {
            //double sumw = alpha_k(refno)*sumw_allrefs;
            Vt().resize(dim, dim, dim);
            Vt().setXmippOrigin();
            Vt().initZeros();

            size_t shdim=hdim;
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Vt())
            if (j < shdim + 1)
                DIRECT_A3D_ELEM(Vt(),k,i,j) = DIRECT_A3D_ELEM(wsumweds[refno],k,i,j)
                                              / sumw_allrefs;
            else
                DIRECT_A3D_ELEM(Vt(),k,i,j) = DIRECT_A3D_ELEM(wsumweds[refno],
                                              (dim-k)%dim,
                                              (dim-i)%dim,
                                              dim-j) / sumw_allrefs;

            CenterFFT(Vt(), true);
            reScaleVolume(Vt(), false);
            fn_tmp = formatString("%s_wedge%06d.vol", fn_base.c_str(), refno + 1);
            Vt.write(fn_tmp);
        }
    }
    fn_tmp = fn_base + "_ref.xmd";
    MDref.write(fn_tmp);

    // Write out FSC curves
    std::ofstream fh;
    fn_tmp = fn_base + ".fsc";
    fh.open((fn_tmp).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR(ERR_IO_NOTOPEN, (std::string)"Cannot write file: " + fn_tmp);

    fh << "# freq. FSC refno 1-n... \n";
    for (int idx = 1; idx < hdim + 1; idx++)
    {
        fh << fsc[0](idx) << " ";
        for (int refno = 0; refno < nr_ref; refno++)
        {
            fh << fsc[refno + 1](idx) << " ";
        }
        fh << "\n";
    }
    fh.close();

    // Write out images with new docfile data
    fn_tmp = fn_base + "_img.xmd";
    MDimg.write(fn_tmp);

    if (iter >= 1)
    {
        //Write out log-file
        MDo.clear();
        MDo.setColumnFormat(false);
        MDo.setComment(cline);
        size_t id = MDo.addObject();
        MDo.setValue(MDL_LL, LL, id);
        MDo.setValue(MDL_PMAX, avefracweight, id);
        MDo.setValue(MDL_SIGMANOISE, sigma_noise, id);
        MDo.setValue(MDL_SIGMAOFFSET, sigma_offset, id);
        MDo.setValue(MDL_ITER, iter, id);
        fn_tmp = fn_base + "_img.xmd";
        MDo.setValue(MDL_IMGMD, fn_tmp, id);
        fn_tmp = fn_base + "_ref.xmd";
        MDo.setValue(MDL_REFMD, fn_tmp, id);

        /*        fn_tmp = fn_base + "_log.xmd";
         MDLog.write(fn_tmp);

         DFl.go_beginning();
         comment = "ml_tomo-logfile: Number of images= " + floatToString(sumw_allrefs);
         comment += " LL= " + floatToString(LL, 15, 10) + " <Pmax/sumP>= " + floatToString(avefracweight, 10, 5);
         DFl.insert_comment(comment);
         comment = "-noise " + floatToString(sigma_noise, 15, 12) + " -offset " + floatToString(sigma_offset/scale_factor, 15, 12) + " -istart " + integerToString(iter + 1);
         DFl.insert_comment(comment);
         DFl.insert_comment(cline);
         DFl.insert_comment("columns: model fraction (1); 1000x signal change (2); resolution (3)");
         fn_tmp = fn_base + ".log";
         DFl.write(fn_tmp);*/

        // Write out MetaData with optimal transformation & references
        //fn_tmp = fn_base + ".doc";
        //DFo.write(fn_tmp);
        // Also write out selfiles of all experimental images,
        // classified according to optimal reference image
        fn_tmp = fn_root + "_ref";
        for (int refno = 0; refno < nr_ref; refno++)
        {
            /*DFo.go_beginning();
             SFo.clear();
             for (int n = 0; n < DFo.dataLineNo(); n++)
             {
             DFo.next();
             fn_tmp = ((DFo.get_current_line()).get_text()).erase(0, 3);
             DFo.adjust_to_data_line();
             if ((refno + 1) == (int)DFo(6))
             SFo.insert(fn_tmp, SelLine::ACTIVE);
             }
             fn_tmp.compose(fn_base + "_ref", refno + 1, "sel");
             SFo.write(fn_tmp);
             */
            MDo.clear();
            MDo.importObjects(MDimg, MDValueEQ(MDL_REF, refno + 1));
            fn_tmp.compose(fn_tmp, refno + 1, "");
            fn_tmp += "_img.xmd";
            MDo.write(fn_tmp);
        }
    }

}

