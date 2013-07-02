/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include "idr_xray_tomo.h"
#include "data/xmipp_image_convert.h"
#include "reconstruct_fourier.h"
#include "reconstruct_art.h"

void ProgIDRXrayTomo::defineParams()
{
    addUsageLine("Iterative Data Refinement applied to the reconstruction of X-ray tomography projections.");
    addSeeAlsoLine("xray_psf_create, xray_import, xray_project");
    //Params
    //    projParam.defineParams(this); // Projection parameters
    addParamsLine("-i <md_file>  : Metadata file with input projections");
    addParamsLine("alias --input;");
    addParamsLine(" [--oroot <rootname=\"idr_xray\">]  : Rootname used for refined projections");
    addParamsLine(" [-o <Vol_file=\"idr_recon_vol.mrc\">]  : Filename for refined output volume");
    addParamsLine("alias --output;");
    addParamsLine(" [--start <start_vol_file=\"\">]  : Start from this basis volume. The reconstruction is performed in the same grid as the one ");
    addParamsLine(" [-s <Ts>]          : Sampling rate (nm)");
    addParamsLine("alias --sampling_rate;");
    addParamsLine("[--thr <threads=1>]           : Number of concurrent threads.");

    addParamsLine("== Iterations options == ");
    addParamsLine("  [-l <...>]                  : Relaxation factor, by default 1.8 (recommended range 0.0 - 2.0). ");
    addParamsLine("                              : A list of lambda values is also accepted as \"-l lambda0 lambda1 ...\"");
    addParamsLine("  [-n <noit=1>]               : Number of iterations");

    addParamsLine("== Reconstruction options == ");
    addParamsLine("[--recons <recons_type=ART>]        : Reconstruction method to be used");
    addParamsLine("       where <recons_type>           ");
    addParamsLine("           ART  <params=\"\">        : wlsART parameters");
    addParamsLine("           fourier <params=\"\">     : fourier parameters");
    addParamsLine("           tomo3d  <params=\"\">     : tomo3d parameters");


    addParamsLine("== Xray options == ");
    addParamsLine("[--psf <psf_param_file=\"\">] : XRay-Microscope parameters file. If not set, then default parameters are chosen.");
    addParamsLine("[--threshold <thr=0.0>]       :+ Normalized threshold relative to maximum of PSF to reduce the volume into slabs");
    //Examples
    //    addExampleLine("Generating a set of projections using a projection parameter file and two threads:", false);
    //    addExampleLine("xmipp_xray_project -i volume.vol --oroot images --params projParams.xmd -s 10 --psf psf_560.xmd --thr 2");
    //    addExampleLine("Generating a single projection at 45 degrees around X axis:", false);
    //    addExampleLine("xmipp_xray_project -i volume.vol -o image.spi --angles 45 0 90 -s 10 --psf psf_560.xmd");
    //    addExampleLine("Generating a single projection at 45 degrees around Y axis with default PSF:", false);
    //    addExampleLine("xmipp_xray_project -i volume.vol -o image.spi --angles 45 90 90 -s 10");
    //Example projection file
    //    addExampleLine("In the following link you can find an example of projection parameter file:",false);
    //    addExampleLine(" ",false);
    //    addExampleLine("http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/tomoProjection.param",false);
    //    addExampleLine(" ",false);
    //    addExampleLine("In the following link you can find an example of X-ray microscope parameters file:",false);
    //    addExampleLine(" ",false);
    //    addExampleLine("http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/xray_psf.xmd",false);
}

/* Read from command line ================================================== */
void ProgIDRXrayTomo::readParams()
{
    fnInputProjMD = getParam("-i");
    fnOutVol = getParam("-o");
    fnRootInter = getParam("--oroot");
    fnStart = getParam("--start");
    //    projParam.readParams(this);

    fnPSF = getParam("--psf");
    sampling  = (checkParam("-s"))? getDoubleParam("-s")*1e-9 : -1 ;
    psfThr = getDoubleParam("--threshold");
    nThr = getIntParam("--thr");

    // Iteration parameters
    StringVector list;
    getListParam("-l", list);
    size_t listSize = list.size();

    if (listSize != 0)
    {
        lambda_list.resizeNoCopy(listSize);

        for (size_t k = 0; k < listSize; k++)
            VEC_ELEM(lambda_list, k) = textToFloat(list[k]);
    }
    else
    {
        lambda_list.resizeNoCopy(1);
        VEC_ELEM(lambda_list, 0) = 1.8;
    }

    itNum = getIntParam("-n");

    // Reconstruction
    if (STR_EQUAL(getParam("--recons"), "ART"))
        reconsMethod = RECONS_ART;
    else if (STR_EQUAL(getParam("--recons"), "fourier"))
        reconsMethod = RECONS_FOURIER;
    else if (STR_EQUAL(getParam("--recons"), "tomo3d"))
        reconsMethod = RECONS_TOMO3D;



}//readParams



void ProgIDRXrayTomo::preRun()
{
    psf.read(fnPSF);
    psf.verbose = verbose;
    psf.nThr = nThr;

    psf.calculateParams(sampling, sampling, psfThr);

    // Threads stuff
    barrier = new Barrier(nThr);
    //Create threads to start working
    thMgr = new ThreadManager(nThr);

    // Reading the input projections. Tilt angle values are needed
    projMD.read(fnInputProjMD);

    // Setting the intermediate filenames
    fnInterProjs = fnRootInter + "_inter_projs";
    fnInterProjs += (reconsMethod == RECONS_TOMO3D)? ".mrc" : ".stk";

    fnInterProjsMD = fnInterProjs.replaceExtension("xmd");
    fnInterAngles = fnRootInter + "_angles.txt";

    // We make sure to delete intermediate files
    fnInterProjs.deleteFile();
    fnInterProjsMD.deleteFile();

    FileName fnProjs;
    projMD.getValue(MDL_IMAGE, fnProjs, projMD.firstObject());
    // We generate the MD of temporary projections to be used for further reconstructions
    interProjMD = projMD;
    interProjMD.replace(MDL_IMAGE, fnProjs.removePrefixNumber(), fnInterProjs);
    interProjMD.write(fnInterProjsMD);

    String arguments = formatString("-i %s -o %s -s",
                                    fnInputProjMD.c_str(), fnInterProjs.c_str());
    ProgConvImg conv;
    conv.read(arguments);
    conv.run();

    // Saving angles in plain text to pass to tomo3D
    if ( reconsMethod == RECONS_TOMO3D )
    {
        std::vector<MDLabel> desiredLabels;
        desiredLabels.push_back(MDL_ANGLE_TILT);
        projMD.writeText(fnInterAngles, &desiredLabels);
    }

    if (fnStart.empty()) // if initial volume is not passed,then we must calculate it
    {
        // Projections must be in an MRC stack to be passed
        if ( reconsMethod == RECONS_TOMO3D )
            reconstruct(fnInterProjsMD, fnOutVol);
        else
        {
            reconstruct(fnInputProjMD, fnOutVol);
            phantom.read(fnOutVol);
        }
    }
    else
        phantom.read(fnStart);
}

void ProgIDRXrayTomo::run()
{
    preRun();

    Projection proj, stdProj;
    Image<double> fixedProj, prevFProj;
    MultidimArray<double> mFixedProj, mPrevFProj;


    double rot, tilt, psi;

    double lambda = VEC_ELEM(lambda_list, 0);

    FileName fnProj;

    for (int iter = 0; iter < itNum; iter++)
    {
        std::cout << "== Iteration " << iter << std::endl;
        std::cout << "Projecting ..." <<std::endl;

        // Fix the reconstructed volume to physical density values
        double factor = 1/sampling;
        MULTIDIM_ARRAY(phantom.iniVol) *= factor;

        initProgress(projMD.size());
        size_t imgNo = 1; // Image number
        double meanError = 0.;

        FOR_ALL_OBJECTS_IN_METADATA(projMD)
        {
            projMD.getValue(MDL_IMAGE , fnProj, __iter.objId);
            projMD.getValue(MDL_ANGLE_ROT , rot, __iter.objId);
            projMD.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
            projMD.getValue(MDL_ANGLE_PSI , psi, __iter.objId);
            proj.setAngles(rot, tilt, psi);

            XrayRotateAndProjectVolumeOffCentered(phantom, psf, proj, stdProj, YSIZE(MULTIDIM_ARRAY(phantom.iniVol)),
                                                  XSIZE(MULTIDIM_ARRAY(phantom.iniVol)));

            /* Image normalization to optimize for reconstruction.
             * To correct for self-attenuation (without considering 3DPSF)
             * I = -ln(Iab)
             */
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(MULTIDIM_ARRAY(proj))
            dAi(MULTIDIM_ARRAY(proj),n) = -log(1. - dAi(MULTIDIM_ARRAY(proj),n));


            fixedProj.read(fnProj);
            mFixedProj.alias(MULTIDIM_ARRAY(fixedProj));

            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mFixedProj)
            dAi(mFixedProj, n) = (dAi(mFixedProj, n) - dAi(MULTIDIM_ARRAY(proj),n))*lambda +  dAi(MULTIDIM_ARRAY(stdProj),n);

            prevFProj.read(fnInterProjs, DATA, imgNo); // To calculate meanError
            mPrevFProj.alias(MULTIDIM_ARRAY(prevFProj));

            // Write the refined projection
            fixedProj.write(fnInterProjs, imgNo, true, WRITE_REPLACE);

            // debug stuff //
            if (verbose > 5)
        {
            proj.write(fnRootInter + "_debug_proj.stk", imgNo , true, WRITE_REPLACE)
                ;
                stdProj.write(fnRootInter + "_debug_std_proj.stk", imgNo , true, WRITE_REPLACE);

                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mFixedProj)
                dAi(mFixedProj, n) = dAi(MULTIDIM_ARRAY(stdProj),n) - dAi(MULTIDIM_ARRAY(proj),n)*lambda ;

                fixedProj.write(fnRootInter + "_debug_diff_proj.stk", imgNo , true, WRITE_REPLACE);
            }

            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mFixedProj)
            meanError += fabs(dAi(mPrevFProj, n) - dAi(mFixedProj, n));

            ++imgNo;
            // Update progress bar
            setProgress();
        }
        // Once calculated all fixed projection, let reconstruct again
        endProgress();

        meanError /= (NZYXSIZE(mFixedProj)*projMD.size());

        std::cout << "Mean error = " << meanError <<std::endl;

        std::cout << "Reconstructing ..." << std::endl;
        reconstruct(fnInterProjsMD, fnOutVol);

        if ( (iter < itNum-1) && (reconsMethod != RECONS_TOMO3D) ) // Since we process the output of tomo3d, this is already read
            phantom.read(fnOutVol);

    }

    // Finish progress bar

    if (reconsMethod == RECONS_TOMO3D) // Trick to save the processed output of tomo3d
        phantom.iniVol.write(fnOutVol);

}

void ProgIDRXrayTomo::reconstruct(const FileName &fnProjsMD, const FileName &fnVol)
{
    if (reconsMethod == RECONS_TOMO3D)
    {
        MetaData MD(fnProjsMD);
        FileName fnProjs;
        MD.getValue(MDL_IMAGE, fnProjs, MD.firstObject());

        reconsTomo3D(fnInterAngles, fnProjs.removePrefixNumber(), fnVol);
        phantom.read(fnVol);
        MULTIDIM_ARRAY(phantom.iniVol).resetOrigin();
        MULTIDIM_ARRAY(phantom.iniVol).reslice(VIEW_Y_NEG);
        //FIXME: this is a temporary fixing to match the reconstructed volume with phantom
        double factor = 1/13.734;
        MULTIDIM_ARRAY(phantom.iniVol) *= factor;
    }
    else
    {
        reconsProgram = createReconsProgram(fnProjsMD, fnVol);
        reconsProgram->run();
        delete reconsProgram;
    }
}

int reconsTomo3D(const MetaData& MD, const FileName& fnOut, const String& params)
{
    // Saving angles in plain text to pass to tomo3D
    FileName fnAngles("idr_xray_tomo_to_tomo3d_angles.txt");
    std::vector<MDLabel> desiredLabels;
    desiredLabels.push_back(MDL_ANGLE_TILT);
    MD.writeText(fnAngles, &desiredLabels);

    FileName fnProjs;

    MD.getValue(MDL_IMAGE, fnProjs, MD.firstObject());

    // Projections must be in an MRC stack to be passed
    if ( !(fnProjs.isInStack() && (fnProjs.contains("mrc") || fnProjs.contains("st"))))
    {
        //FIXME: To be implemented the conversion to an MRC stack

        REPORT_ERROR(ERR_IMG_NOREAD, "reconsTromo3D: Image format cannot be read by tomo3D.");
    }

    return reconsTomo3D(fnAngles, fnProjs.removePrefixNumber(), fnOut, params);
}

int reconsTomo3D(const FileName& fnAngles, const FileName& fnProjs,
                 const FileName& fnOut, const String& params)
{
    String exec("tomo3d -n -f -v 0 -a "+fnAngles+" -i "+fnProjs.removeBlockNameOrSliceNumber()+" -o "+fnOut+" "+params);
    return system(exec.c_str());
}


ProgReconsBase * ProgIDRXrayTomo::createReconsProgram(const FileName &input, const FileName &output)
{
    //get reconstruction extra params
    String arguments = getParam("--recons", 1) +
                       formatString(" -v 0 --thr %d -i %s -o %s", nThr, input.c_str(), output.c_str());
    ProgReconsBase * program;
    //    std::cerr << "DEBUG_JM: ProgMLRefine3D::createReconsProgram" <<std::endl;
    //    std::cerr << "DEBUG_JM: arguments: " << arguments << std::endl;
    //    std::cerr << "DEBUG_JM: input: " << input << std::endl;

    if (reconsMethod == RECONS_FOURIER)
    {
        program = new ProgRecFourier();
        //force use of weights and the verbosity will be the same of this program
        //-i and -o options are passed for avoiding errors, this should be changed
        //when reconstructing
        //        arguments += " --weight";
        program->read(arguments);
        return program;
    }
    else if (reconsMethod == RECONS_ART)//use of Art
    {
        //REPORT_ERROR(ERR_NOT_IMPLEMENTED,"not implemented reconstruction throught wlsArt");
        program = new ProgReconsART();
        FileName fn_tmp(arguments);
        //        arguments += " --WLS";
        //        if (fn_symmask.empty() && checkParam("--sym"))
        //            arguments += " --sym " + fn_sym;
        if (!fn_tmp.contains("-n "))
            arguments += " -n 10";
        if (!fn_tmp.contains("-l "))
            arguments += " -l 0.2";
        //        if (!fn_tmp.contains("-k "))
        //            arguments += " -k 0.5";

        program->read(arguments);
        return program;
    }
    return NULL;
}



