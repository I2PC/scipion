/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
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

#include "ml_refine3d.h"

#include <data/xmipp_fft.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/morphology.h>
#include <data/grids.h>
#include <data/blobs.h>
#include <data/symmetries.h>
#include <data/projection.h>
#include "directions.h"
#include "reconstruct_art.h"
#include "reconstruct_fourier.h"
#include "fourier_filter.h"
#include "symmetrize.h"
#include "volume_segment.h"


//#define DEBUG
//Macro to obtain the iteration base name
#define FN_ITER(iter, suffix) formatString("%sextra/iter%03d/%s",fn_root.c_str(), (iter), (suffix))
#undef FN_ITER_BASE
#define FN_ITER_BASE(iter) FN_ITER(iter, "vol")
#define FN_INITIAL_BASE FN_ITER_BASE(0)
#define FN_PROJECTIONS_MD FN_EXTRA("projections.xmd")
#define FN_PROJECTIONS    FN_EXTRA("projections.stk")

#define FN_NOISE_VOLBASE  FN_EXTRA("noise_vol")
#define FN_CREF_VOLBASE   FN_EXTRA("cref_vol")
// Filename generation for a volume from a base
#define COMPOSE_VOL_FN(fn, volno, base) fn.compose(base, volno, "vol")


ProgMLRefine3D::ProgMLRefine3D(bool fourier)
{
    fourier_mode = fourier;
    if (fourier)
        ml2d = new ProgMLF2D();
    else
        ml2d = new ProgML2D();
    rank = 0;
    size = 1;
}

ProgMLRefine3D::~ProgMLRefine3D()
{
    CLOSE_LOG();
    delete ml2d;
}

// Usage ===================================================================
void ProgMLRefine3D::defineParams()
{
    addUsageLine("Separate structurally heterogenous data sets into homogeneous classes by a");
    addUsageLine("multi-reference 3D-angular refinement using a maximum-likelihood(ML) target function.");
    //Add some params from 2D
    //Setting some flags for the params definitions
    ml2d->defaultNiter = 25;
    ml2d->referenceExclusive = false;
    ml2d->allowFastOption = false;
    ml2d->allowRestart = false;
    ml2d->allowIEM = true;

    //basic params
    ml2d->defineBasicParams(this);
    addParamsLine(" [ --ang <float=10> ]           : Angular sampling (degrees) ");
    //extra params
    ml2d->defineAdditionalParams(this, "==+ ML additional options: ==");
    addParamsLine("==+ Additional options: ==");
    addParamsLine("--recons <recons_type=wlsART>       : Reconstruction method to be used");
    addParamsLine("       where <recons_type>           ");
    addParamsLine("           wlsART  <params=\"\">     : wlsART parameters");
    addParamsLine("           fourier <params=\"\">     : fourier parameters");
    addParamsLine(" [ --nostart ]                      : Start wlsART reconstructions from all-zero volumes ");
    addParamsLine(" [ --sym <symfile=c1> ]             : Symmetry group ");
    addParamsLine(" [ --low_pass <freq=-1> ]           : Low-pass filter volume every iteration ");
    addParamsLine(" [ --sym_mask <maskfile=\"\"> ]     : Local symmetry (only inside mask) ");
    addParamsLine(" [ --tilt <min=-91.> <max=91.> ]    : Minimum and maximum values for restriction tilt angle search ");
    addParamsLine(" [ --perturb ]                      : Randomly perturb reference projection directions ");
    //hidden params
    ml2d->defineHiddenParams(this);
    addParamsLine(" [--solvent <filename=\"\">]");
    addParamsLine(" [--prob_solvent]");
    addParamsLine(" [--threshold_solvent <float=999.>]");
    addParamsLine(" [--deblob_solvent]");
    addParamsLine(" [--dilate_solvent <int=0>]");
    addParamsLine(" [--skip_reconstruction]");
}

// Read ===================================================================
void ProgMLRefine3D::readParams()
{
    bool do_restart = false;

    // Generate new command line for restart procedure
    /// FIXME: restart has to be re-thought
    /*
    if (checkParam( "-restart"))
{
        String   comment, cline = "";
        DocFile  DFi;
        FileName fn_tmp;

        do_restart = true;
        DFi.read(getParameter(argc, argv, "-restart"));
        DFi.go_beginning();
        comment = (DFi.get_current_line()).get_text();
        if (fourier_mode && strstr(comment.c_str(), "MLFalign2D-logfile") == NULL)
        {
            std::cerr << "Error!! Docfile is not of MLFalign2D-logfile type. " << std::endl;
            exit(1);
        }
        else if (!fourier_mode && strstr(comment.c_str(), "MLalign2D-logfile") == NULL)
        {
            std::cerr << "Error!! Docfile is not of MLalign2D-logfile type. " << std::endl;
            exit(1);
        }
        else
        {
            char *copy;
            copy = NULL;
            DFi.next();
            // get new part of command line (with -istart)
            comment = (DFi.get_current_line()).get_text();
            DFi.next();
            // get original command line
            cline = (DFi.get_current_line()).get_text();
            comment = comment + cline;
            // regenerate command line
            argv2 = NULL;
            argc2 = 0;
            generateCommandLine(comment, argc2, argv2, copy);
            // Get number of volumes and names to generate SFvol
            if (fourier_mode)
                fn_root = getParam( "-o", "mlf3d");
            else
                fn_root = getParam( "-o", "ml3d");
            fn_vol = getParam( "-vol");
            istart = getIntParam( "-istart"));
            if (Is_VolumeXmipp(fn_vol))
            {
                SFvol.clear();
                SFvol.addObject();
                SFvol.setValue(MDL_IMAGE, fn_vol);
                SFvol.setValue(MDL_ENABLED, 1);
            }
            else
            {
                SFvol.read(fn_vol);
            }
            SFvol.removeObjects(MDValueEQ(MDL_ENABLED, -1));
            Nvols = SFvol.size();

            SFvol.clear();
            for (int ivol = 0; ivol < Nvols; ivol++)
            {
                fn_tmp = fn_root + "_it";
                fn_tmp.compose(fn_tmp, istart - 1, "");
                if (Nvols > 1)
                {
                    fn_tmp += "_vol";
                    fn_tmp.compose(fn_tmp, ivol + 1, "");
                }
                fn_tmp += ".vol";
                SFvol.addObject();
                SFvol.setValue(MDL_IMAGE, fn_tmp);
                SFvol.setValue(MDL_ENABLED, 1):
                    }
            fn_vol = fn_root + "_it";
            fn_vol.compose(fn_vol, istart - 1, "");
            fn_vol += "_restart.sel";
            SFvol.write(fn_vol);
        }
}
    else
{
    */
    // no restart, just copy argc to argc2 and argv to argv2
    //argc2 = argc;
    //argv2 = argv;
    //   }


    //Read Refine3d parameters
    fn_sel = getParam( "-i");
    fn_root = getParam("--oroot");


    //    if (fourier_mode)
    //        fn_root = getParam( "-o", "mlf3d");
    //    else
    //        fn_root = getParam( "-o", "ml3d");

    if (!do_restart)
    {
        // Fill volume selfile
        fn_ref = getParam( "--ref");
        mdVol.read(fn_ref);
        Nvols = mdVol.size();
    }

    angular = getDoubleParam( "--ang");
    fn_sym = getParam( "--sym");
    eps = getDoubleParam( "--eps");
    Niter = getIntParam( "--iter");
    istart = 1;//getIntParam( "-istart");
    tilt_range0 = getDoubleParam( "--tilt", 0);
    tilt_rangeF = getDoubleParam( "--tilt", 1);
    fn_symmask = getParam( "--sym_mask");
    lowpass = getDoubleParam( "--low_pass");

    wlsart_no_start = checkParam( "--nostart");
    //    if (checkParam("--wlsart"))
    //    {
    //        wlsart_lambda = getDoubleParam("--wlsart", 0);
    //        wlsart_kappa = getDoubleParam("--wlsart", 1);
    //        wlsart_Niter = getIntParam("--wlsart", 2);
    //    }
    do_perturb = checkParam( "--perturb");

    // Hidden for now
    fn_solv = getParam( "--solvent");

    if (STR_EQUAL(getParam("--recons"), "wlsART"))
        recons_type = RECONS_ART;
    else if (STR_EQUAL(getParam("--recons"), "fourier"))
        recons_type = RECONS_FOURIER;

    do_prob_solvent = checkParam( "--prob_solvent");
    threshold_solvent = getDoubleParam( "--threshold_solvent");
    do_deblob_solvent = checkParam( "--deblob_solvent");
    dilate_solvent = getIntParam( "--dilate_solvent");
    skip_reconstruction = checkParam( "--skip_reconstruction");

    // Checks
    if (lowpass > 0.5)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "Digital frequency for low-pass filter should be smaller than 0.5");

    //Read ml2d params
    ml2d->do_ML3D = true;
    ml2d->verbose = verbose; // 2d inherits verbosity from 3d
    ml2d->read(argc, argv, false);

    if (!checkParam("--psi_step"))
        ml2d->psi_step = angular;

    ml2d->fn_img = fn_sel;
    ml2d->fn_root = fn_root + ml2d->defaultRoot;
    ml2d->fn_ref = FN_PROJECTIONS_MD;
    ml2d->do_mirror = true;

    //add empty string string now, this will be later
    //overwrite with the iteration base name
    reconsOutFnBase.push_back("");
    reconsMdFn.push_back("");

    if (fourier_mode)
    {
        //For fourier case add also the _noise and _cref
        reconsOutFnBase.push_back(FN_CREF_VOLBASE);
        reconsOutFnBase.push_back(FN_NOISE_VOLBASE);
        {//make fn_root local scope
            String fn_root = this->fn_root + ml2d->defaultRoot; // this is need for the following macro
            reconsMdFn.push_back(FN_CREF_IMG_MD);
        }
        reconsMdFn.push_back(FN_NOISE_IMG_MD);

    }
    CREATE_LOG(LOG_FN(fn_root));
}

void ProgMLRefine3D::show()
{
    LOG_FUNCTION();

    if (verbose == 0)
        return;

    std::cout << " -----------------------------------------------------------------" << std::endl;
    std::cout << " | Read more about this program in the following publication:    |" << std::endl;
    if (fourier_mode)
        std::cout << " |  Scheres ea. (2007)  Structure, 15, 1167-1177                 |" << std::endl;
    else
        std::cout << " |  Scheres ea. (2007)  Nature Methods, 4, 27-29                 |" << std::endl;
    std::cout << " |                                                               |" << std::endl;
    std::cout << " |    *** Please cite it if this program is of use to you! ***   |" << std::endl;
    std::cout << " -----------------------------------------------------------------" << std::endl;
    std::cout << "--> Maximum-likelihood multi-reference 3D-refinement" << std::endl;
    if (Nvols == 1)
        std::cout << "  Initial reference volume : " << fn_ref << std::endl;
    else
    {
        std::cout << "  Selfile with references  : " << fn_ref << std::endl;
        std::cout << "    with # of volumes      : " << Nvols << std::endl;
    }
    std::cout << "  Experimental images:     : " << fn_sel << std::endl;
    std::cout << "  Angular sampling rate    : " << angular << std::endl;
    std::cout << "  Symmetry group:          : " << fn_sym << std::endl;
    if (fn_symmask != "")
        std::cout << "  Local symmetry mask      : " << fn_symmask << std::endl;
    std::cout << "  Output rootname          : " << fn_root << std::endl;
    std::cout << "  Convergence criterion    : " << eps << std::endl;
    if (lowpass > 0)
        std::cout << "  Low-pass filter          : " << lowpass << std::endl;
    if (tilt_range0 > -91. || tilt_rangeF < 91.)
        std::cout << "  Limited tilt range       : " << tilt_range0 << "  " << tilt_rangeF << std::endl;
    if (wlsart_no_start)
        std::cout << "  -> Start wlsART reconstructions from all-zero volumes " << std::endl;
    if (recons_type == RECONS_FOURIER)
        std::cout << "  -> Use fourier-interpolation instead of wlsART for reconstruction" << std::endl;
    if (do_prob_solvent)
        std::cout << "  -> Perform probabilistic solvent flattening" << std::endl;
    std::cout << " -----------------------------------------------------------------" << std::endl;
}

void ProgMLRefine3D::createSampling()
{
    LOG_FUNCTION();
    FileName fn_sym_loc;
    // Precalculate sampling
    mysampling.setSampling(angular);
    fn_sym_loc = fn_symmask.empty() ? fn_sym : "c1";

    if (!mysampling.SL.isSymmetryGroup(fn_sym_loc, symmetry, sym_order))
        REPORT_ERROR(ERR_NUMERICAL, (String)"ml_refine3d::run Invalid symmetry" +  fn_sym_loc);
    mysampling.SL.readSymmetryFile(fn_sym_loc);
    mysampling.computeSamplingPoints(true, tilt_rangeF, tilt_range0);
    mysampling.removeRedundantPointsExhaustive(symmetry, sym_order, true, 0.75 * angular);
    nr_projections = mysampling.no_redundant_sampling_points_angles.size();
}

// Fill sampling and create DFlib
void ProgMLRefine3D::produceSideInfo()
{
    LOG_FUNCTION();
    //Create sampling
    createSampling();
    show();
    // Write starting volume(s) to disc with correct name for iteration loop
    copyVolumes();
    // Project volumes and store projections in a metadata
    projectVolumes(ml2d->MDref);
//    //FIXME: this is for concurrency problem...remove after that
//    FileName myImg = fn_root + formatString("images_node%02d.xmd", rank);
//    MetaData(fn_sel).write(myImg);
//    ml2d->fn_img = myImg;
    //2d initialization
    LOG("before ml2d->produceSideInfo");
    ml2d->produceSideInfo();
}

void ProgMLRefine3D::produceSideInfo2()
{
    LOG_FUNCTION();
    ml2d->produceSideInfo2();
    ml2d->refs_per_class = nr_projections;
    ml2d->show();
    Nvols *= ml2d->factor_nref;
    ml2d->Iold.clear(); // To save memory
}

void ProgMLRefine3D::run()
{
    LOG_FUNCTION();
    bool converged = false;


    // Get input parameters
    produceSideInfo();
    produceSideInfo2();
    ml2d->createThreads();

    //Local image to read data
    Image<double> img;
    FileName fn;
    bool doProject = false;

    // Loop over all iterations
    for (ml2d->iter = ml2d->istart; !converged && ml2d->iter <= ml2d->Niter; ml2d->iter++)
    {
        iter = ml2d->iter; //keep updated the iter class variable

        //Make path for iterations result files
        getIterExtraPath(fn_root, iter);

        if (verbose)
            std::cout << formatString("--> 3D-EM volume refinement:  iteration %d of %d", iter, Niter) << std::endl;

        LOG("==============================================");
        LOG(formatString("ML3D: Iteration %d of %d", iter, Niter));
        LOG_LEVEL(Iteration);

        for (ml2d->current_block = 0; ml2d->current_block < ml2d->blocks; ml2d->current_block++)
        {
            LOG(formatString("ML3D: BEGIN BLOCK %d of %d", ml2d->current_block, ml2d->blocks));
            LOG_LEVEL(Block);
            // Project volumes, already done for first iteration, first block
            if (doProject)// || ml2d->current_block > 0)
            {
                projectVolumes(ml2d->MDref);
                int refno = 0;

                // Read new references from disc (I could just as well keep them in memory, maybe...)
                FOR_ALL_OBJECTS_IN_METADATA(ml2d->MDref)
                {
                    ml2d->MDref.getValue(MDL_IMAGE, fn, __iter.objId);
                    img.read(fn);
                    img().setXmippOrigin();
                    ml2d->model.Iref[refno]() = img();
                    if (++refno == ml2d->model.n_ref) //avoid reading noise and c_ref projections
                        break;
                }
            }
            LOG("Calling ML2D Expectation");
            // Integrate over all images
            ml2d->expectation();
            LOG("Calling ML2D Maximization");
            ml2d->maximization();

            //do not reconstruction on special iteration 0 until last block
            if (iter > SPECIAL_ITER || ml2d->current_block == ml2d->blocks - 1)//last block
            {

                LOG("Writing ML2D references");
                // Write out 2D reference images (to be used in reconstruction)
                ml2d->writeOutputFiles(ml2d->model, OUT_REFS);

                // Jump out before 3D reconstruction
                // (Useful for some parallelization protocols)
                if (skip_reconstruction)
                    exit(1);

                if (fourier_mode)
                    makeNoiseImages();
                // Reconstruct new volumes from the reference images
                //update the base name for current iteration
                reconsOutFnBase[0] = FN_ITER_BASE(iter);
                reconsMdFn[0] = ml2d->outRefsMd;
                reconstructVolumes();
                // Update the reference volume selection file
                updateVolumesMetadata();
                // post-process the volumes
                postProcessVolumes();
                doProject = true;
            }
        } // end loop blocks

        if (fourier_mode)
            calculate3DSSNR(ml2d->spectral_signal);

        // Check convergence
        converged = checkConvergence();

        // End 2D iteration
        LOG("Calling ml2d->endIteration");
        ml2d->endIteration();
    } // end loop iterations

    if (verbose)
    {
        std::cout << (converged ?
                      "--> Optimization converged!" :
                      "--> Optimization was stopped before convergence was reached!")
        << std::endl;
    }

    // Write converged output ML2D files
    LOG("Calling ml2d->writeOutputFiles, OUT_FINAL");
    ml2d->writeOutputFiles(ml2d->model, OUT_FINAL);
    ml2d->destroyThreads();

}//end of function run

void ProgMLRefine3D::createEmptyFiles(int type)
{
    size_t dim, idum, idumLong;
    getImageSizeFromFilename(fn_sel, dim, idum, idum, idumLong);
    Image<double> img;

    if (type == EMPTY_PROJECTIONS)
        createEmptyFile(FN_PROJECTIONS, dim, dim, 1, Nvols*nr_projections, true);
    //    {
    //        img().initZeros(dim, dim);
    //        img.write(FN_PROJECTIONS, Nvols * nr_projections, true, WRITE_OVERWRITE);
    //    }
    else if (type == EMPTY_VOLUMES)
    {
        img().initZeros(dim, dim, dim);
        for (size_t i = 0; i < reconsOutFnBase.size(); ++i)
            createEmptyFile(reconsOutFnBase[i], dim, dim, dim, Nvols, true);
        //img.write(reconsOutFnBase[i], Nvols, true, WRITE_OVERWRITE);
    }
}

// Projection of the reference (blob) volume =================================
void ProgMLRefine3D::projectVolumes(MetaData &mdProj)
{

    LOG_FUNCTION();

    Image<double>                vol;
    FileName                      fn_base = FN_PROJECTIONS, fn_tmp;
    Projection                    proj;
    double                       rot, tilt, psi = 0.;
    size_t                        nl, nr_dir, id, bar_step;
    int                           volno;

    // Here all nodes fill SFlib and DFlib, but each node actually projects
    // only a part of the projections. In this way parallellization is obtained
    // Total number of projections
    nl = Nvols * nr_projections;
    bar_step = XMIPP_MAX(1, nl / 60);

    // Initialize projections output metadata
    mdProj.clear();

    if (verbose)
    {
        std::cout << formatString("--> projecting %d volumes x %d projections...", Nvols, nr_projections) << std::endl;
        //init_progress_bar(nl);
    }

    createEmptyFiles(EMPTY_PROJECTIONS);

    // Loop over all reference volumes
    volno = nr_dir = 0;

    //std::cerr << "DEBUG_JM: ProgMLRefine3D::projectVolumes" <<std::endl;
    MDIterator iter(mdVol);
    for (size_t i = 0; i < Nvols; ++i)
    //FOR_ALL_OBJECTS_IN_METADATA(mdVol)
    {
        mdVol.getValue(MDL_IMAGE, fn_tmp, iter.objId);
        //std::cerr << "DEBUG_JM: fn_tmp: " << fn_tmp << std::endl;
        vol.read(fn_tmp);
        vol().setXmippOrigin();
        ++volno;

        for (int ilib = 0; ilib < nr_projections; ++ilib)
        {
            ++nr_dir;
            fn_tmp.compose(nr_dir, fn_base);
            rot = XX(mysampling.no_redundant_sampling_points_angles[ilib]);
            tilt = YY(mysampling.no_redundant_sampling_points_angles[ilib]);

            // Parallelization: each rank projects and writes a different direction
            if (nr_dir % size == rank)
            {
                projectVolume(vol(), proj, vol().rowNumber(), vol().colNumber(), rot, tilt, psi);
                //proj.setEulerAngles(rot, tilt, psi);
                //std::cerr << formatString("DEBUG_JM: Proyecting vol: %s, rot: %f, tilt: %f, psi: %f", fn_tmp.c_str(), rot, tilt, psi) <<std::endl;
                proj.write(fn_tmp);
                //proj.write(formatString("%s_iter%d.stk", fn_tmp.c_str(), iter));
            }

            id = mdProj.addObject();
            mdProj.setValue(MDL_IMAGE, fn_tmp, id);
            mdProj.setValue(MDL_ENABLED, 1, id);
            mdProj.setValue(MDL_ANGLE_ROT, rot, id);
            mdProj.setValue(MDL_ANGLE_TILT, tilt, id);
            mdProj.setValue(MDL_ANGLE_PSI, psi, id);
            mdProj.setValue(MDL_REF3D, volno, id);

            if (verbose && (nr_dir % bar_step == 0))
                progress_bar(nr_dir);
        }
        iter.moveNext();
    }

    if (verbose)
    {
        //progress_bar(nl);
        std::cout << " -----------------------------------------------------------------" << std::endl;
    }

    // Only the master write the complete SFlib
    if (rank == 0)
    {
        fn_tmp = FN_PROJECTIONS_MD;
        mdProj.write(fn_tmp);
        //Just copying projections md for debugging
        //fn_tmp.copyFile(formatString("%s_iter%d_projections.xmd", fn_root.c_str(), iter));
        //Also copying projections stack
        //FileName(FN_PROJECTIONS).copyFile(formatString("%s_iter%d_projections.stk", fn_root.c_str(), iter));
    }

}

// Make noise images for 3D SSNR calculation ===================================
void ProgMLRefine3D::makeNoiseImages()
{

    Image<double> img;
    std::vector<Image<double> > & Iref = ml2d->model.Iref;
    FileName   fn_noise(FN_NOISE_IMG), fn_img;
    MetaData    mdNoise(ml2d->MDref);
    int refno = 0;
    MDRow row;

    FOR_ALL_OBJECTS_IN_METADATA(mdNoise)
    {
        img = Iref[refno];
        img().initZeros();
        img().addNoise(0, 1, "gaussian");

        if (Iref[refno].weight() > 1.)
            img() /= sqrt(Iref[refno].weight());
        fn_img.compose(++refno, fn_noise);
        img.write(fn_img);
        mdNoise.setValue(MDL_IMAGE, fn_img, __iter.objId);
        if (refno == ml2d->model.n_ref)
            break;
    }
    fn_noise = FN_NOISE_IMG_MD;
    mdNoise.write(fn_noise);
}

ProgReconsBase * ProgMLRefine3D::createReconsProgram(FileName &input, FileName &output)
{
    //get reconstruction extra params
    String arguments = getParam("--recons", 1) +
        formatString(" -v 0 --thr %d -i %s -o %s", ml2d->threads, input.c_str(), output.c_str());
    ProgReconsBase * program;
//    std::cerr << "DEBUG_JM: ProgMLRefine3D::createReconsProgram" <<std::endl;
//    std::cerr << "DEBUG_JM: arguments: " << arguments << std::endl;
//    std::cerr << "DEBUG_JM: input: " << input << std::endl;

    if (recons_type == RECONS_FOURIER)
    {
        program = new ProgRecFourier();
        //force use of weights and the verbosity will be the same of this program
        //-i and -o options are passed for avoiding errors, this should be changed
        //when reconstructing
        arguments += " --weight";
        program->read(arguments);
        return program;
    }
    else if (recons_type == RECONS_ART)//use of wlsArt
    {
        //REPORT_ERROR(ERR_NOT_IMPLEMENTED,"not implemented reconstruction throught wlsArt");
        program = new ProgReconsART();
        FileName fn_tmp(arguments);
        arguments += " --WLS";
        if (fn_symmask.empty() && checkParam("--sym"))
          arguments += " --sym " + fn_sym;
        if (!fn_tmp.contains("-n "))
          arguments += " -n 10";
        if (!fn_tmp.contains("-l "))
          arguments += " -l 0.2";
        if (!fn_tmp.contains("-k "))
          arguments += " -k 0.5";

        bool noise_vols = input.contains("_cref_") || input.contains("_noise_");
        if (!noise_vols && !wlsart_no_start)
        {
             arguments += " --save_basis";
             if (iter > 1)
             {
               String currIter = FN_ITER_BASE(iter);
               String prevIter = FN_ITER_BASE(iter - 1);
               arguments += " --start ";
               arguments += output.replaceSubstring(currIter, prevIter).replaceExtension("basis");
             }
        }
        program->read(arguments);
        return program;
        //        if (fn_symmask != "")
        //            art_prm.fn_sym = "";


        //        BasicARTParameters   art_prm;
        //        Plain_ART_Parameters   dummy;
        //        GridVolume             new_blobs;
        //        GridVolume             start_blobs;
        //        if (verbose)
        //            std::cerr << "--> weighted least-squares ART reconstruction " << std::endl;
        //
        //        // Read ART parameters from command line & I/O with outer loop of Refine3d
        //        art_prm.read(argc, argv);
        //        art_prm.WLS = true;
        //        if (fn_symmask != "")
        //            art_prm.fn_sym = "";
        //        if (!checkParam( "-n"))
        //            art_prm.no_it = 10;
        //        if (!checkParam( "-l"))
        //        {
        //            art_prm.lambda_list.resize(1);
        //            art_prm.lambda_list.initConstant(0.2);
        //        }
        //        if (!checkParam( "-k"))
        //        {
        //            art_prm.kappa_list.resize(1);
        //            art_prm.kappa_list.initConstant(0.5);
        //        }
        //        art_prm.fn_sel = "xxxx"; //fixme
        //        art_prm.fn_root = "xxxx";
        //        if (noise == 1 || noise == 2)
        //        {
        //            art_prm.fn_start = "";
        //            art_prm.tell = false;
        //        }
        //        else if (!wlsart_no_start)
        //        {
        //            art_prm.tell = TELL_SAVE_BASIS;
        //            art_prm.fn_start = fn_blob;
        //        }
        //        // Reconstruct using weighted least-squares ART
        //        Basic_ROUT_Art(art_prm, dummy, new_vol, new_blobs);
    }
    return NULL;
}

// Reconstruction using the ML-weights ==========================================
void ProgMLRefine3D::reconstructVolumes()
{
  LOG_FUNCTION();

    FileName fn_vol, fn_vol_prev, fn_one;
    MetaData mdOne, mdProj, mdOutVols;
    size_t id;


    ProgReconsBase * reconsProgram;
    int volno_index  = 0;

    createEmptyFiles(EMPTY_VOLUMES);

    for (size_t i = 0; i < reconsOutFnBase.size(); ++i)
    {
        for (int volno = 1; volno <= (int)Nvols; ++volno)
        {
            volno_index = Nvols * i + volno - 1;
            String &fn_base = reconsOutFnBase[i];
            COMPOSE_VOL_FN(fn_vol, volno, fn_base);
            //for now each node reconstruct one volume
            if (volno_index % size == rank)
            {
                //fn_vol.compose(volno, fn_base);
				mdProj.read(reconsMdFn[i]);
                fn_one.compose(fn_base, volno, "projections.xmd");
                // Select only relevant projections to reconstruct
                mdOne.importObjects(mdProj, MDValueEQ(MDL_REF3D, volno));
                mdOne.write(fn_one);
                // Set input/output for the reconstruction algorithm
				reconsProgram = createReconsProgram(fn_one, fn_vol);
                reconsProgram->run();
				delete reconsProgram;
            }
            //Store output volumes, avoid noise and cref volumes
            //Only need to be done by master node
            if (i == 0 && rank == 0)
            {
               id = mdOutVols.addObject();
               mdOutVols.setValue(MDL_IMAGE, fn_vol, id);
               mdOutVols.setValue(MDL_ENABLED, 1, id);
            }
        }//for volno
    }//for reconsOutFnBase

    if (rank == 0)
    {
      FileName fn = FN_ITER_VOLMD();
      mdOutVols.write(fn);
    }

}

void ProgMLRefine3D::calculate3DSSNR(MultidimArray<double> &spectral_signal)
{

  LOG_FUNCTION();

    MetaData                    mdNoiseAll, mdNoiseOne;
    MultidimArray<std::complex<double> >  Faux;
    Image<double>              vol, nvol;
    FileName                    fn_tmp, fn_tmp2;
    MultidimArray<double>      alpha_signal, alpha_noise, input_signal, avg_alphaS, avg_alphaN;
    MultidimArray<double>      alpha_T, alpha_N, Msignal, Maux, Mone, mask;
    Projection                  proj;
    size_t                      c, dim, idum;
    size_t                      idumLong;
    double                      volweight, weight, rot, tilt, psi = 0.;
    Matrix1D<int>               center(2);
    MultidimArray<int>          radial_count;

    // Read in noise reconstruction and calculate alpha's
    mdNoiseAll.read(FN_NOISE_IMG_MD);
    getImageSize(mdNoiseAll, dim, idum, idum, idumLong);

    center.initZeros();
    proj().resize(dim, dim);
    proj().setXmippOrigin();
    mask.resize(dim, dim);
    mask.setXmippOrigin();
    RaisedCosineMask(mask, dim / 2 - 2, dim / 2);

    if (verbose)
    {
        std::cout << "--> calculating 3D-SSNR ..." << std::endl;
        initProgress(mdNoiseAll.size());
    }

    FileName fn_noise_base = FN_NOISE_VOLBASE;
    FileName fn_cref_base = FN_CREF_VOLBASE;
    double inv_dim2 = 1. / (double)(dim * dim);

    for (int volno = 1; volno <= (int)Nvols; ++volno)
    {
    	COMPOSE_VOL_FN(fn_tmp, volno, fn_noise_base);
    	COMPOSE_VOL_FN(fn_tmp2, volno, fn_cref_base);
//        fn_tmp.compose(volno, fn_noise_base);
//        fn_tmp2.compose(volno, fn_cref_base);

        nvol.read(fn_tmp);
        vol.read(fn_tmp2);
        nvol().setXmippOrigin();
        vol().setXmippOrigin();
        Mone.resize(dim, dim);
        Mone.initConstant(inv_dim2);
        Mone.setXmippOrigin();

        mdNoiseOne.clear();
        mdNoiseOne.importObjects(mdNoiseAll, MDValueEQ(MDL_REF3D, volno));
        c = 0;
        volweight = 0.;
        bool first_time = true;

        FOR_ALL_OBJECTS_IN_METADATA(mdNoiseOne)
        {
            mdNoiseOne.getValue(MDL_WEIGHT, weight, __iter.objId);
            mdNoiseOne.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
            mdNoiseOne.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);

            // accumulate alpha denominator
            SUM_INIT(alpha_N, Mone * weight);
            // alpha nominator
            projectVolume(nvol(), proj, dim, dim, rot, tilt, psi);
            apply_cont_mask(mask, proj(), proj());
            FourierTransform(proj(), Faux);
            FFT_magnitude(Faux, Maux);
            CenterFFT(Maux, true);
            Maux.setXmippOrigin();
            Maux *= Maux;
            SUM_INIT(alpha_T, Maux * weight);
            // input signal
            projectVolume(vol(), proj, dim, dim, rot, tilt, psi);
            apply_cont_mask(mask, proj(), proj());
            FourierTransform(proj(), Faux);
            FFT_magnitude(Faux, Maux);
            CenterFFT(Maux, true);
            Maux.setXmippOrigin();
            Maux *= Maux;
            SUM_INIT(Msignal, Maux * weight);
            volweight += weight;
            setProgress(++c);
            first_time = false;
        }

        alpha_T.setXmippOrigin();
        alpha_N.setXmippOrigin();
        Msignal.setXmippOrigin();
        alpha_signal.initZeros();
        alpha_noise.initZeros();
        input_signal.initZeros();
        radialAverage(alpha_T, center, alpha_signal, radial_count, true);
        radialAverage(alpha_N, center, alpha_noise, radial_count, true);
        radialAverage(Msignal, center, input_signal, radial_count, true);

//        std::cerr << "DEBUG_JM: alpha_T: " << alpha_T << std::endl;
//        std::cerr << "DEBUG_JM: alpha_N: " << alpha_N << std::endl;
//        std::cerr << "DEBUG_JM: Msignal: " << Msignal << std::endl;

        input_signal *= 1./volweight;

        // Calculate spectral_signal =input_signal/alpha!!
        // Also store averages of alphaN and alphaS for output
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(input_signal)
        {
          dAi(input_signal, i) = (dAi(alpha_signal, i) > 0.) ? dAi(input_signal, i) * dAi(alpha_noise, i) / dAi(alpha_signal, i) : 0.;

        }

        if (volno == 1)
        {
            spectral_signal = input_signal;
            avg_alphaN = alpha_noise;
            avg_alphaS = alpha_signal;
        }
        else
        {
            spectral_signal += input_signal;
            avg_alphaN += alpha_noise;
            avg_alphaS += alpha_signal;
        }
    }
    endProgress();

    double inv_Nvols = 1. / (double)Nvols;
    spectral_signal *= inv_Nvols;
    avg_alphaN *= inv_Nvols;
    avg_alphaS *= inv_Nvols;

    if (verbose)
    {
        fn_tmp = getIterExtraPath(fn_root, iter) + "3dssnr.log";
        std::ofstream out(fn_tmp.c_str(), std::ios::out);
        out  << "#        signal    1/alpha    alpha-S    alpha-N" << std::endl;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(spectral_signal)
        {
            if (i > 0 && i < dim / 2)
            {
                out.width(5);
                out  << integerToString(i);
                out.width(10);
                out <<  floatToString(DIRECT_A1D_ELEM(spectral_signal, i));
                out.width(10);
                out <<  floatToString(DIRECT_A1D_ELEM(avg_alphaN, i) / DIRECT_A1D_ELEM(avg_alphaS, i));
                out.width(10);
                out <<  floatToString(DIRECT_A1D_ELEM(avg_alphaS, i));
                out.width(10);
                out <<  floatToString(DIRECT_A1D_ELEM(avg_alphaN, i));
                out << std::endl;
            }
        }
        out.close();
    }
}

void ProgMLRefine3D::copyVolumes()
{
  LOG_FUNCTION();
    ImageGeneric img;
    getIterExtraPath(fn_root, 0); //Create folder to store volume
    FileName fn_vol, fn_base = FN_INITIAL_BASE;
    size_t volno = 0;

    FOR_ALL_OBJECTS_IN_METADATA(mdVol)
    {
        mdVol.getValue(MDL_IMAGE, fn_vol, __iter.objId);
        img.read(fn_vol);
        //fn_vol.compose(++volno, fn_base);
        COMPOSE_VOL_FN(fn_vol, ++volno, fn_base);
        img.write(fn_vol);
        mdVol.setValue(MDL_IMAGE, fn_vol, __iter.objId);
        mdVol.setValue(MDL_ENABLED, 1, __iter.objId);
    }
}

void ProgMLRefine3D::updateVolumesMetadata()
{
  LOG_FUNCTION();
    FileName fn_vol, fn_base;
    mdVol.clear();

    for (size_t i = 0; i < reconsOutFnBase.size(); ++i)
    {
        fn_base = reconsOutFnBase[i];
        for (size_t volno = 1; volno <= Nvols; ++volno)
        {
        	COMPOSE_VOL_FN(fn_vol, volno, fn_base);
            //fn_vol.compose(volno, fn_base);
            mdVol.setValue(MDL_IMAGE, fn_vol, mdVol.addObject());
            mdVol.setValue(MDL_ENABLED, 1);
        }
    }
}

// Modify reference volume ======================================================
void ProgMLRefine3D::postProcessVolumes()
{
  LOG_FUNCTION();

  ProgVolumeSegment            segm_prm;
    FileName               fn_vol, fn_tmp;
    Image<double>          vol, Vaux, Vsymmask, Vsolv;
    MultidimArray<int>     mask3D;
    double                 in, out;
    Sampling               locsampling;

    // Use local sampling because of symmask
    if (!locsampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
        REPORT_ERROR(ERR_NUMERICAL, (String)"ml_refine3d::run Invalid symmetry" +  fn_sym);
    locsampling.SL.readSymmetryFile(fn_sym);

    if ( !(fn_sym == "c1" || fn_sym == "C1" ) || (lowpass > 0) ||
         (fn_solv != "") || (do_prob_solvent) || (threshold_solvent != 999))
    {
    	LOG_LEVEL(postProcessVolumes_IF);

        FOR_ALL_OBJECTS_IN_METADATA(mdVol)
        {
            mdVol.getValue(MDL_IMAGE, fn_vol, __iter.objId);
            // Read corresponding volume from disc
            LOG("   ProgMLRefine3D::postProcessVolumes READING vol");
            vol.read(fn_vol);
            vol().setXmippOrigin();
            // Store the original volume on disc
            fn_tmp = fn_vol;
            fn_tmp.insertBeforeExtension(".original");
            //fn_tmp = fn_vol + ".original";
            LOG("   ProgMLRefine3D::postProcessVolumes writing ORIGINAL vol");
            vol.write(fn_tmp);

            // Symmetrize if requested
            if (!fn_sym.empty())
            {
            	 LOG("   ProgMLRefine3D::postProcessVolumes applying SYMMETRY");
                Vaux().resize(vol());
                symmetrizeVolume(locsampling.SL, vol(), Vaux());
                // Read local symmetry mask if requested
                if (!fn_symmask.empty())
                {
                    Vsymmask.read(fn_symmask);
                    Vsymmask().setXmippOrigin();
                    if (Vsymmask().computeMax() > 1. || Vsymmask().computeMin() < 0.)
                        REPORT_ERROR(ERR_VALUE_INCORRECT, "ERROR: sym_mask should have values between 0 and 1!");
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Vsymmask())
                    {
                        in = dAkij(Vsymmask(), k, i, j);
                        out = 1. - in;
                        dAkij(vol(), k, i, j) = out * dAkij(vol(), k, i, j) + in * dAkij(Vaux(), k, i, j);
                    }
                    Vsymmask.clear();
                }
                else
                {
                    vol = Vaux;
                }
                Vaux.clear();
            }

            // Filtering the volume
            if (lowpass > 0)
            {
            	LOG("   ProgMLRefine3D::postProcessVolumes applying LOWPASS");
                FourierFilter fmask;
                fmask.raised_w = 0.02;
                fmask.FilterShape = RAISED_COSINE;
                fmask.FilterBand = LOWPASS;
                fmask.w1 = lowpass;
                fmask.applyMaskSpace(vol());
            }

            // Different types of solvent flattening
            if (do_prob_solvent || !fn_solv.empty() || (threshold_solvent != 999))
            {
            	LOG("   ProgMLRefine3D::postProcessVolumes applying SOLVENT");
                if (do_prob_solvent)
                {
                    // A. Probabilistic solvent flattening
                    // Write already processed volume to disc (for segment program)
                    vol.write(fn_vol);
                    segm_prm.read(argc, argv);
                    segm_prm.fn_vol = fn_vol;
                    segm_prm.fn_mask = fn_vol + ".solv";
                    segm_prm.do_prob = true;
                    segm_prm.show();
                    segm_prm.produce_side_info();
                    segm_prm.segment(Vsolv);
                }
                else if (threshold_solvent != 999)
                {
                    // B. Perform flooding and separate_objects-like solvent mask
                    Vsolv = vol;
                    Vsolv().threshold("below", threshold_solvent, 0.);
                    // The following is because binarize() seems buggy
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Vsolv())
                    {
                        if (dAkij(Vsolv(), k, i, j) != 0.)
                            dAkij(Vsolv(), k, i, j) = 1.;
                    }
                }
                else if (fn_solv != "")
                {
                    // C. Read user-provided solvent mask from disc
                    Vsolv.read(fn_solv);
                    if (Vsolv().computeMax() > 1. || Vsolv().computeMin() < 0.)
                        REPORT_ERROR(ERR_VALUE_INCORRECT, "ERROR: solvent mask should have values between 0 and 1!");
                }
                // Binarize Vsolv, avoiding buggy Vsolv().binarize()
                if (do_deblob_solvent || dilate_solvent > 0)
                {
                    Vsolv().threshold("below", 0.5, 0.);
                    // The following is because binarize() seems buggy
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Vsolv())
                    {
                        if (dAkij(Vsolv(), k, i, j) != 0.)
                            dAkij(Vsolv(), k, i, j) = 1.;
                    }
                }
                if (do_deblob_solvent)
                {
                    int object_no;
                    double nr_vox, max_vox = 0.;
                    Image<double> label;
                    object_no = labelImage3D(Vsolv(), label());
                    max_vox = 0;
                    for (int o = 0; o <= object_no; o++)
                    {
                        Vaux() = label();
                        FOR_ALL_ELEMENTS_IN_ARRAY3D(Vaux())
                        {
                            Vaux(k, i, j) = Vaux(k, i, j) == o;
                        }
                        nr_vox = Vaux().sum();
                        if (o != 0 && (nr_vox > max_vox))
                        {
                            max_vox = nr_vox;
                            Vsolv() = Vaux();
                        }
                    }
                    label.clear();
                }
                // Dilate solvent mask (only for binary masks)
                // Dilate several times, result is summed iteratively
                if (dilate_solvent > 0)
                {
                    Image<double> Vsum;
                    Vsum() = Vsolv();
                    for (int i = 0; i < dilate_solvent; i++)
                    {
                        dilate3D(Vsolv(), Vaux(), 18, 0, 1);
                        Vsum() = Vsum() + Vaux();
                        Vsolv() = Vaux();
                    }
                    Vsum() /= (double)(dilate_solvent + 1);
                    Vsolv() = Vsum();
                    Vsum.clear();
                }
                // Apply solvent mask
                Vsolv() = 1. - Vsolv();
                double solvavg = 0., sumsolv = 0.;
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(vol())
                {
                    solvavg += dAkij(Vsolv(), k, i, j) * dAkij(vol(), k, i, j);
                    sumsolv += dAkij(Vsolv(), k, i, j);
                }
                solvavg /= sumsolv;
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Vsolv())
                {
                    dAkij(vol(), k, i, j) -= dAkij(Vsolv(), k, i, j) * (dAkij(vol(), k, i, j) - solvavg);
                }
            }

            // (Re-) write post-processed volume to disc
            LOG("Before WRITING vol");
            vol.write(fn_vol);
            LOG("After WRITING vol");

        }
        if (verbose)
            std::cout << " -----------------------------------------------------------------" << std::endl;
    }
}

// Convergence check ===============================================================
bool ProgMLRefine3D::checkConvergence()
{
    LOG_FUNCTION();

    Image<double>        vol, old_vol, diff_vol;
    FileName               fn_base, fn_base_old, fn_vol;
    Mask            mask_prm;
    MultidimArray<int>    mask3D;
    double                 signal, change;
    int                    dim;
    bool                   converged = true;

    if (iter == 0)
        return false;

    if (verbose)
        std::cout << "--> Checking convergence " << std::endl;

    fn_base = FN_ITER_BASE(iter);
    fn_base_old = FN_ITER_BASE(iter - 1);

    for (size_t volno = 1; volno <= Nvols; ++volno)
    {
        // Read corresponding volume from disc
    	COMPOSE_VOL_FN(fn_vol, volno, fn_base);
        //fn_vol.compose(volno, fn_base);
    	//std::cerr << "DEBUG_JM: fn_vol: " << fn_vol << std::endl;
        vol.read(fn_vol);
        vol().setXmippOrigin();
        dim = vol().rowNumber();
        old_vol().initZeros(vol());
        diff_vol().initZeros(vol());

        // Only consider voxels within the spherical mask
        mask_prm.R1 = dim / 2;
        mask_prm.type = BINARY_CIRCULAR_MASK;
        mask_prm.mode = INNER_MASK;
        mask_prm.generate_mask(vol());
        //fn_vol.compose(volno, fn_base_old);
        COMPOSE_VOL_FN(fn_vol, volno, fn_base_old);
        old_vol.read(fn_vol);
        //std::cerr << "DEBUG_JM: oldVol: " << fn_vol << std::endl;
        diff_vol() = vol() - old_vol();
        mask_prm.apply_mask(old_vol(), old_vol());
        mask_prm.apply_mask(diff_vol(), diff_vol());
        change = diff_vol().sum2();
        signal = old_vol().sum2();
        //std::cerr << "DEBUG_JM: change: " << change << std::endl;
        //std::cerr << "DEBUG_JM: signal: " << signal << std::endl;
        if (change / signal > eps)
            converged = false;
        if (verbose)
            std::cout << formatString("--> Relative signal change volume %d = %f", volno, change / signal) << std::endl;
    }

    return converged;

}
