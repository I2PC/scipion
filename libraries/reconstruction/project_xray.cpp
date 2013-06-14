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
#include "project_xray.h"

void ProgXrayProject::defineParams()
{
    addUsageLine("Generate projections as in a X-ray microscope from a 3D Xmipp volume.");
    addSeeAlsoLine("xray_psf_create, xray_import, phantom_create");
    //Params
    projParam.defineParams(this); // Projection parameters
    addParamsLine("== Xray options == ");
    addParamsLine(" [--sampling_rate <Ts>]       : Sampling rate of the volume to be projected (nm).");
    addParamsLine("                              : If empty, same value as X-Y plane sampling from PSF.");
    addParamsLine("alias -s;");
    addParamsLine("[--psf <psf_param_file=\"\">] : XRay-Microscope parameters file. If not set, then default parameters are chosen.");
    addParamsLine("[--focal_shift+ <delta_z=0.0>] : Shift along optical axis between the center of the phantom and the center of the psf in microns");
    addParamsLine("alias -f;");
    addParamsLine("[--threshold+ <thr=0.0>]     : Normalized threshold relative to maximum of PSF to reduce the volume into slabs");
    addParamsLine("[--std_proj]                  : Save also standard projections, adding _std suffix to filenames");
    addParamsLine("[--thr <threads=1>]           : Number of concurrent threads");
    //Examples
    addExampleLine("Generating a set of projections using a projection parameter file and two threads:", false);
    addExampleLine("xmipp_xray_project -i volume.vol --oroot images --params projParams.xmd -s 10 --psf psf_560.xmd --thr 2");
    addExampleLine("Generating a single projection at 45 degrees around X axis:", false);
    addExampleLine("xmipp_xray_project -i volume.vol -o image.spi --angles 45 0 90 -s 10 --psf psf_560.xmd");
    addExampleLine("Generating a single projection at 45 degrees around Y axis with default PSF:", false);
    addExampleLine("xmipp_xray_project -i volume.vol -o image.spi --angles 45 90 90 -s 10");
    //Example projection file
    addExampleLine("In the following link you can find an example of projection parameter file:",false);
    addExampleLine(" ",false);
    addExampleLine("http://sourceforge.net/p/testxmipp/code/ci/3.0/tree/input/tomoProjection.param",false);
    addExampleLine(" ",false);
    addExampleLine("In the following link you can find an example of X-ray microscope parameters file:",false);
    addExampleLine(" ",false);
    addExampleLine("http://sourceforge.net/p/testxmipp/code/ci/3.0/tree/input/xray_psf.xmd",false);
}

/* Read from command line ================================================== */
void ProgXrayProject::readParams()
{
    //    fn_proj_param = getParam("-i");
    //    fn_sel_file   = getParam("-o");
    projParam.readParams(this);

    fn_psf_xr = getParam("--psf");
    dxo  = (checkParam("-s"))? getDoubleParam("-s")*1e-9 : -1 ;
    psfThr = getDoubleParam("--threshold");
    nThr = getIntParam("--thr");
    save_std_projs = checkParam("--std_proj");

    psf.read(fn_psf_xr);
    psf.setFocalShift(getDoubleParam("--focal_shift")*1e-6);
    psf.verbose = verbose;
    psf.nThr = nThr;
}

//Some global threads management variables
Mutex mutex;
Barrier * barrier;
ThreadManager * thMgr;


void ProgXrayProject::preRun()
{
    randomize_random_generator();

    psf.calculateParams(dxo, -1, psfThr);
    phantom.read(projParam);
    //Correct the rotation axis displacement in projectionParams from pixels to meters
    projParam.raxis *= psf.dxo;

    // Threads stuff

    barrier = new Barrier(nThr);

    //Create threads to start working
    thMgr = new ThreadManager(nThr);

}

void ProgXrayProject::postRun()
{
    //Terminate threads and free memory
    //    if (td != NULL)
    //        delete td;
    delete thMgr;
    delete barrier;
}

void ProgXrayProject::run()
{
    preRun();

    int expectedNumProjs = FLOOR((projParam.tiltF-projParam.tilt0)/projParam.tiltStep);
    int numProjs=0;

    if (verbose > 0)
    {
        std::cerr << "Projecting ...\n";
        if (!(projParam.show_angles))
            init_progress_bar(expectedNumProjs);
    }
    projMD.setComment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
    double tRot,tTilt,tPsi,rot,tilt,psi;
    FileName fn_proj;  // Projection name
    size_t idx = FIRST_IMAGE;

    // Project
    for (double angle=projParam.tilt0; angle<=projParam.tiltF; angle+=projParam.tiltStep)
    {
        if (projParam.singleProjection)
            fn_proj = projParam.fnOut;
        else
            fn_proj.compose(idx, projParam.fnRoot + ".stk");

        // Choose Center displacement ........................................
        double shiftX = rnd_gaus(projParam.Ncenter_avg, projParam.Ncenter_dev);
        double shiftY = rnd_gaus(projParam.Ncenter_avg, projParam.Ncenter_dev);
        Matrix1D<double> inPlaneShift(3);
        VECTOR_R3(inPlaneShift,shiftX,shiftY,0);

        projParam.calculateProjectionAngles(proj,angle, 0,inPlaneShift);

        //Reset thread task distributor
        //        td->clear();

        // Really project ....................................................
        if (!projParam.only_create_angles)
            XrayRotateAndProjectVolumeOffCentered(phantom, psf, proj, stdProj, projParam.proj_Ydim, projParam.proj_Xdim);

        // Add noise in angles and voxels ....................................
        proj.getEulerAngles(tRot, tTilt,tPsi);

        rot  = tRot  + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);
        tilt = tTilt + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);
        psi  = tPsi  + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);

        proj.setEulerAngles(rot,tilt,psi);

        size_t objId = projMD.addObject();
        if (!projParam.only_create_angles)
        {
            // Here we subtract the differences to the background illumination (at this moment normalized to 1) //TODO
            MULTIDIM_ARRAY(proj) = 1.0 - MULTIDIM_ARRAY(proj);
            proj.write(fn_proj, ALL_IMAGES, !projParam.singleProjection, WRITE_REPLACE);

            if (save_std_projs)
                stdProj.write(fn_proj.insertBeforeExtension("_std"), ALL_IMAGES, !projParam.singleProjection, WRITE_REPLACE);

            projMD.setValue(MDL_IMAGE,fn_proj,objId);
        }

        projMD.setValue(MDL_ANGLE_ROT,rot,objId);
        projMD.setValue(MDL_ANGLE_TILT,tilt,objId);
        projMD.setValue(MDL_ANGLE_PSI,psi,objId);
        projMD.setValue(MDL_ANGLE_ROT2,tRot,objId);
        projMD.setValue(MDL_ANGLE_TILT2,tTilt,objId);
        projMD.setValue(MDL_ANGLE_PSI2,tPsi,objId);
        projMD.setValue(MDL_SHIFT_X,shiftX,objId);
        projMD.setValue(MDL_SHIFT_Y,shiftY,objId);

        IMGMATRIX(proj).addNoise(projParam.Npixel_avg, projParam.Npixel_dev, "gaussian");

        // Save ..............................................................
        if (projParam.show_angles)
        {
            std::cout << "N      Rot     Tilt     Psi" <<std::endl;
            std::cout << idx << "\t" << proj.rot() << "\t"
            << proj.tilt() << "\t" << proj.psi() << std::endl;
        }
        else if ((expectedNumProjs % XMIPP_MAX(1, numProjs / 60))  == 0 && verbose > 0)
            progress_bar(numProjs);

        numProjs++;
        idx++;
    }
    if (!(projParam.show_angles))
        progress_bar(expectedNumProjs);

    // Save metadata file with angles and shift info
    if (!projParam.singleProjection)
    {
        projMD.setComment("Angles rot,tilt and psi contain noisy projection angles and rot2,tilt2 and psi2 contain actual projection angles");
        projMD.write(projParam.fnRoot + ".xmd");
    }

    postRun();

    return;
}

void XrayRotateAndProjectVolumeOffCentered(XrayProjPhantom &phantom, XRayPSF &psf, Projection &P, Projection &standardP,
        int Ydim, int Xdim)
{

    int iniXdim, iniYdim, iniZdim, newXdim, newYdim, newZdim, rotXdim, rotZdim;
    int yinit, xinit;

    iniXdim = XSIZE(MULTIDIM_ARRAY(phantom.iniVol));
    iniYdim = YSIZE(MULTIDIM_ARRAY(phantom.iniVol));
    iniZdim = ZSIZE(MULTIDIM_ARRAY(phantom.iniVol));

    // Projection offset in pixels
    Matrix1D<double> offsetNV(3);
    P.getShifts(XX(offsetNV), YY(offsetNV), ZZ(offsetNV));

    //    xOffsetN = P.Xoff()/psf.dxo;
    //    yOffsetN = P.Yoff()/psf.dxo;


    /* If Xdim is greater than Zdim, then as we rotate the volume the rotated Zdim must be
     * great enough to cover all the volume.
     */
    /* By now, we are going to keep it unused */
    //    if (true)
    //    {
    //        double radi, theta0, theta, tilt;
    //
    //        radi = sqrt(iniZdim*iniZdim + iniXdim*iniXdim);
    //        theta0 = atan2(iniZdim, iniXdim);
    //        tilt = ((P.tilt() < 90)? P.tilt(): 180 - P.tilt())* PI / 180;
    //
    //        rotZdim = XMIPP_MAX(ROUND(radi * sin(ABS(tilt) + ABS(theta0))), iniZdim);
    //        rotXdim = XMIPP_MAX(ROUND(radi * cos(ABS(tilt) - ABS(theta0))), iniXdim);
    //        rotYdim = iniYdim;
    //
    //        //    rotZdim = iniZdim;
    //        //    rotXdim = iniXdim;
    //
    //        newXdim = rotXdim + 2*ABS(XX(offsetNV));
    //        newYdim = iniYdim + 2*ABS(YY(offsetNV));
    //        newZdim = rotZdim;
    //    }
    //    else
    {
        rotXdim = iniXdim;
        rotZdim = iniZdim;
    }

    newXdim = (int)(rotXdim + 2*fabs(XX(offsetNV)));
    newYdim = (int)(iniYdim + 2*fabs(YY(offsetNV)));
    newZdim = rotZdim;

    // We set the dimensions only to obtain the values of starting X,Y,Z
    //    phantom.rotVol.setDimensions(rotXdim, rotYdim, rotZdim, 1);
    phantom.rotVol.setDimensions(newXdim, newYdim, newZdim, 1);

    phantom.rotVol.setXmippOrigin();

    yinit = STARTINGY((phantom.rotVol));
    if (YY(offsetNV) > 0)
        yinit -= (int)(2 * YY(offsetNV));

    xinit = STARTINGX((phantom.rotVol));
    if (XX(offsetNV) > 0)
        xinit -= (int)(2 * XX(offsetNV));

    if (psf.verbose > 1)
    {
        if (newXdim != iniXdim || newYdim != iniYdim)
        {
            std::cout << std::endl;
            std::cout << "--------------------------------" << std::endl;
            std::cout << "XrayProject::Volume_offCentered:" << std::endl;
            std::cout << "--------------------------------" << std::endl;
            std::cout << "(X,Y,Z) shifts = (" << XX(offsetNV) << "," << YY(offsetNV) << ","
            << ZZ(offsetNV) << ") um" << std::endl;
            std::cout << "Image resize (Nx,Ny): (" << iniXdim << "," << iniYdim << ") --> ("
            << newXdim << "," << newYdim << ") " << std::endl;
        }
#ifdef DEBUG

        std::cout <<"yoffsetN "<< YY(offsetNV) <<std::endl;
        std::cout <<"xoffsetN "<< XX(offsetNV) <<std::endl;
        std::cout <<"yinit    " << yinit  <<std::endl;
        std::cout <<"yend     "    << yend  <<std::endl;
        std::cout <<"xinit    "   << xinit  <<std::endl;
        std::cout <<"xend     "    << xend  <<std::endl;
        std::cout <<"zinit    "   << zinit  <<std::endl;
        std::cout <<"zend     "    << zend  <<std::endl;
#endif

    }



    //    phantom.rotVol.setMmap(true);
    phantom.rotVol.resizeNoCopy(newZdim, newYdim, newXdim);
    phantom.rotVol.setXmippOrigin();

    // Rotate volume ....................................................

    Matrix2D<double> R, T;

    ZZ(offsetNV) = 0; // We are not interested in apply this Zshift
    translation3DMatrix(offsetNV, T);
    Euler_angles2matrix(P.rot(), P.tilt(), P.psi(), R, true);

    double outside = 0; //phantom.iniVol.getPixel(0,0,0,0);
    MULTIDIM_ARRAY(phantom.iniVol).setXmippOrigin();

    applyGeometry(1, phantom.rotVol, MULTIDIM_ARRAY(phantom.iniVol), T*R, IS_NOT_INV, DONT_WRAP, outside);

    psf.adjustParam(phantom.rotVol);


    /// Calculation of IgeoVol to optimize process
    MultidimArray<double> IgeoVol, IgeoZb;
    IgeoZb.resize(1, 1, YSIZE(phantom.rotVol), XSIZE(phantom.rotVol),false);
    IgeoZb.initConstant(1.);

    calculateIgeo(phantom.rotVol, psf.dzo, IgeoVol, MULTIDIM_ARRAY(standardP), IgeoZb, psf.nThr, thMgr);



    projectXrayVolume(phantom.rotVol, IgeoVol, psf, P, NULL, thMgr);

    int outXDim = XMIPP_MIN(Xdim,iniXdim);
    int outYDim = XMIPP_MIN(Ydim,iniYdim);

    P().selfWindow(-ROUND(outYDim/2),
                   -ROUND(outXDim/2),
                   -ROUND(outYDim/2) + outYDim -1,
                   -ROUND(outXDim/2) + outXDim -1);

}//XrayRotateAndProjectVolumeOffCentered

void projectXrayVolume(MultidimArray<double> &muVol,
                       MultidimArray<double> &IgeoVol,
                       XRayPSF &psf, Projection &P, MultidimArray<double> * projNorm, ThreadManager * thMgr)
{

    /* Now we have to split the phantom rotated volume into slabs, taking into
     * account the possible different resolution and Z size.
     *
     * PhantomSlabIdx are the indexes of rotVol defining the slabs according to
     * the splitting done to PSFVol.
     *
     * psfSlicesIdx are the indexes of the z planes of rotVol whose psf are
     * used for the slabs.
     */
    std::vector<int> phantomSlabIdx, psfSlicesIdx;

    // Search for the PSFslab of the beginning of the volume
    int firstSlab = (int)(STARTINGZ(muVol)*psf.dzo/psf.dzoPSF);

    if (!XMIPP_EQUAL_ZERO(psf.slabThr))
    {
        for (size_t kk = 0; kk < psf.slabIndex.size(); ++kk)
        {
            if (firstSlab < psf.slabIndex[kk])
            {
                firstSlab = kk-1;
                break;
            }
        }

        // Searching the equivalent index in rotvol for the indexes for the Slabs of PSFVol
        phantomSlabIdx.push_back(STARTINGZ(muVol));

        for (size_t kk = firstSlab+1; kk < psf.slabIndex.size(); ++kk)
        {
            int tempK = (int)(psf.slabIndex[kk] * psf.dzoPSF / psf.dzo);

            if (tempK <= FINISHINGZ(muVol))
            {
                phantomSlabIdx.push_back(tempK);
                int tempKK = psf.slabIndex[kk-1];
                int psfMeanSlice = (int)((tempK + tempKK)*0.5 * psf.dzoPSF / psf.dzo);
                psfSlicesIdx.push_back(psfMeanSlice);
            }
            else
            {
                phantomSlabIdx.push_back(FINISHINGZ(muVol));
                int psfMeanSlice = (tempK + phantomSlabIdx[kk-1])/2;
                psfSlicesIdx.push_back(psfMeanSlice);
                continue;
            }
        }
    }
    else
    {
        for (int k = STARTINGZ(muVol); k <= FINISHINGZ(muVol); ++k)
        {
            phantomSlabIdx.push_back(k);
            psfSlicesIdx.push_back(k);
        }
        // To define the last range
        phantomSlabIdx.push_back(FINISHINGZ(muVol)+1);
    }

    XrayThreadArgument projThrData;
    Barrier barrier(thMgr->threads);

    projThrData.psf = &psf;
    projThrData.muVol = &muVol;
    projThrData.IgeoVol = &IgeoVol;
    projThrData.projOut = &P;
    projThrData.projNorm = projNorm;
    projThrData.phantomSlabIdx = & phantomSlabIdx;
    projThrData.psfSlicesIdx = & psfSlicesIdx;

    //Create the job handler to distribute thread jobs
    size_t blockSize, numberOfJobs = psfSlicesIdx.size() ;
    blockSize = 1;
    ThreadTaskDistributor td(numberOfJobs, blockSize);
    projThrData.td = &td;
    projThrData.barrier = &barrier;

    //the really really final project routine, I swear by Snoopy.
    //    project_xr(psf,side.rotPhantomVol,P, idxSlice);

    thMgr->run(threadXrayProject,(void*) &projThrData);
    //    P.write("projection_after_threads.spi");
}

/* Produce Side Information ================================================ */
void XrayProjPhantom::read(
    const ParametersProjectionTomography &prm)
{
    read(prm.fnPhantom);
}

void XrayProjPhantom::read(const FileName &fnVol)
{
    //    iniVol.readMapped(prm.fnPhantom);
    iniVol.read(fnVol);
    //    MULTIDIM_ARRAY(iniVol).setXmippOrigin();
    //    rotVol.resizeNoCopy(MULTIDIM_ARRAY(iniVol));

}

/// Generate an X-ray microscope projection for volume vol using the microscope configuration psf
void threadXrayProject(ThreadArgument &thArg)
{

    int thread_id = thArg.thread_id;

    XrayThreadArgument *dataThread = (XrayThreadArgument*) thArg.data;
    const XRayPSF &psf = *(dataThread->psf);
    MultidimArray<double> &muVol =  *(dataThread->muVol);
    MultidimArray<double> &IgeoVol =  *(dataThread->IgeoVol);
    Image<double> &imOutGlobal = *(dataThread->projOut);
    MultidimArray<double> * projNorm = dataThread->projNorm;
    std::vector<int> &phantomSlabIdx = *(dataThread->phantomSlabIdx);
    std::vector<int> &psfSlicesIdx = *(dataThread->psfSlicesIdx);
    ParallelTaskDistributor * td = dataThread->td;
    Barrier * barrier = dataThread->barrier;

    size_t first = 0, last = 0;

    if (thread_id == 0)
    {
        muVol.setXmippOrigin();
        MULTIDIM_ARRAY(imOutGlobal).resizeNoCopy(psf.Niy, psf.Nix);
        MULTIDIM_ARRAY(imOutGlobal).initZeros();
        MULTIDIM_ARRAY(imOutGlobal).setXmippOrigin();
    }

    barrier->wait();

    MultidimArray<double> imOut(psf.Niy, psf.Nix), projNormTemp(psf.Niy, psf.Nix);
    imOut.setXmippOrigin();
    projNormTemp.setXmippOrigin();

    MultidimArray<double> imTemp(psf.Noy, psf.Nox),intExp(psf.Noy, psf.Nox),imTempSc(imOut),*imTempP;
    intExp.setXmippOrigin();
    imTemp.setXmippOrigin();
    imTempSc.setXmippOrigin();

    // Parallel thread jobs
    while (td->getTasks(first, last))
    {
        //#define DEBUG
#ifdef DEBUG
        Image<double> _Im;
#endif

        int kini = phantomSlabIdx[first];
        int kend = phantomSlabIdx[last + 1] - 1 ;

        //        double deltaZSlab =  psf.dzo*(kend-kini+1);

        intExp.initZeros();

        FOR_ALL_ELEMENTS_IN_ARRAY2D(intExp)
        {
            for (int k = kini; k <= kend; k++)
            {
                double tempValue = A3D_ELEM(muVol,k,i,j);
                A2D_ELEM(intExp,i,j) += tempValue;
            }
            A2D_ELEM(imTemp,i,j) = A3D_ELEM(IgeoVol,kini,i,j)*(1 - exp( -A2D_ELEM(intExp,i,j) * psf.dzo));
        }
#ifdef DEBUG
        _Im().alias(intExp);
        _Im.write("projectXR-intExp.spi");

        _Im().alias(imTemp);
        _Im.write("projectXR-imTemp.spi");
#endif

        switch (psf.AdjustType)
        {
        case PSFXR_INT:
            imTempP = &imTempSc;
            scaleToSize(LINEAR,*imTempP,imTemp,psf.Nix,psf.Niy);
            break;

        case PSFXR_STD:
            imTempP = &imTemp;
            break;

        case PSFXR_ZPAD:
            //    (*imTempSc).resize(imTemp);
            imTempP = &imTempSc;
            imTemp.window(*imTempP,-ROUND(psf.Niy/2)+1,-ROUND(psf.Nix/2)+1,ROUND(psf.Niy/2)-1,ROUND(psf.Nix/2)-1);
            break;
        }

#ifdef DEBUG
        _Im().alias(*imTempP);
        _Im.write("projectXR-imTempEsc_before.spi");
#endif

        psf.applyOTF(*imTempP, imOutGlobal.Zoff()- psfSlicesIdx[first]*psf.dzo);

#ifdef DEBUG

        _Im.write("projectXR-imTempEsc_after.spi");
#endif

        // Adding to the thread temporal imOut
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
        dAij(imOut,i,j) += dAij(*imTempP,i,j);

#ifdef DEBUG

        _Im().alias(imOut);
        _Im.write("projectXR-imout.spi");
#endif

        /// Calculate projNorm

        if (projNorm != NULL)
        {

            FOR_ALL_ELEMENTS_IN_ARRAY2D(imTemp)
            A2D_ELEM(imTemp,i,j) = A3D_ELEM(IgeoVol,kini,i,j) * psf.dzo;

            //            IgeoVol.getSlice(kini, imTemp);
            psf.applyOTF(imTemp, imOutGlobal.Zoff()- psfSlicesIdx[first]*psf.dzo);

            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imTemp)
            dAij(projNormTemp,i,j) += dAij(imTemp,i,j)*dAij(imTemp,i,j);
        }


    }

    //Lock to update the total addition and multiply for the pixel width (to avoid multiply several times)
    mutex.lock();
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imOut)
    dAij(MULTIDIM_ARRAY(imOutGlobal),i,j) += dAij(imOut,i,j);

    if (projNorm != NULL)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(projNormTemp)
        dAij(*projNorm,i,j) += dAij(projNormTemp,i,j);
    }
    mutex.unlock();
#ifdef DEBUG

    imOutGlobal.write("projectXR-imoutGlobal.spi");
#endif

    barrier->wait();

    if (thread_id==0)
    {
        // The output is the addition of the  accumulated differences, the real projection is background illumination minus this
        //        MULTIDIM_ARRAY(imOutGlobal) = 1.0- MULTIDIM_ARRAY(imOutGlobal);

        switch (psf.AdjustType)
        {
        case PSFXR_STD:
            // Do nothing;
            break;
        case PSFXR_INT:
            selfScaleToSize(LINEAR,imOutGlobal(), psf.Nox, psf.Noy);
            break;
        case PSFXR_ZPAD:
            MULTIDIM_ARRAY(imOutGlobal).selfWindow(-ROUND(psf.Noy/2)+1,-ROUND(psf.Nox/2)+1,ROUND(psf.Noy/2)-1,ROUND(psf.Nox/2)-1);
            break;
        }

        MULTIDIM_ARRAY(imOutGlobal).setXmippOrigin();

#ifdef DEBUG

        imOutGlobal.write("projectXR-imoutGlobal_negative.spi");
#endif

    }
    //    std::cerr << "th" << thread_id << ": Thread Finished" <<std::endl;

#undef DEBUG
}

void calculateIgeo(MultidimArray<double> &muVol, double sampling,
                   MultidimArray<double> &IgeoVol, MultidimArray<double> &cumMu,
                   MultidimArray<double> &IgeoZb, int nThreads, ThreadManager * ThrMgr)
{
    // Checking sizes
    if ( (XSIZE(muVol) != XSIZE(IgeoZb)) || (YSIZE(muVol) != YSIZE(IgeoZb)) )
        REPORT_ERROR(ERR_MULTIDIM_SIZE, "CalculateIgeo: IgeoZb size does not match muVol size.");

    IgeoVol.initZeros(muVol);
    cumMu.initZeros(YSIZE(muVol), XSIZE(muVol));

    ThreadManager * myThMgr;

    if (ThrMgr == NULL)
        myThMgr = new ThreadManager(nThreads);
    else
        myThMgr = ThrMgr;

    CIGTArgument iGeoArgs;
    iGeoArgs.samplingZ = sampling;
    iGeoArgs.muVol = &muVol;
    iGeoArgs.IgeoVol = &IgeoVol;
    iGeoArgs.cumMu = &cumMu;
    iGeoArgs.IgeoZb = & IgeoZb;
    iGeoArgs.td = new ThreadTaskDistributor(YSIZE(muVol),YSIZE(muVol)/nThreads/2);

    myThMgr->run(calculateIgeoThread,&iGeoArgs);

    delete iGeoArgs.td;

    if (ThrMgr == NULL)
        delete myThMgr;
}

void calculateIgeoThread(ThreadArgument &thArg)
{
    CIGTArgument *dataThread = (CIGTArgument*) thArg.data;

    double sampling = dataThread->samplingZ;
    MultidimArray<double> &muVol = *(dataThread->muVol);
    MultidimArray<double> &IgeoVol = *(dataThread->IgeoVol);
    MultidimArray<double> &cumMu = *(dataThread->cumMu);
    MultidimArray<double> &IgeoZb = *(dataThread->IgeoZb);
    ParallelTaskDistributor * td = dataThread->td;

    // Resizing
    //    IgeoVol.resize(muVol, false);

    size_t first, last;
    while (td->getTasks(first, last))
    {
        //        first--;
        //        last--;
        for (size_t i = first; i <= last; ++i)
        {
            //            for (int j=0; j<XSIZE(muVol); ++j)
            //                dAkij(IgeoVol,0,i,j) = dAij(IgeoZb,i,j)*exp(-dAkij(muVol,0,i,j)*sampling);
            //
            //            for (int k = 1; k < ZSIZE(muVol); ++k)
            //                for (int j=0; j<XSIZE(muVol); ++j)
            //                    dAkij(IgeoVol,k,i,j) = dAkij(IgeoVol,k-1,i,j)*exp(-dAkij(muVol,k,i,j)*sampling);

            for (size_t k = 0; k < ZSIZE(muVol); ++k)
                for (size_t j=0; j<XSIZE(muVol); ++j)
                {
                    dAij(cumMu,i,j) += dAkij(muVol,k,i,j)*sampling;
                    dAkij(IgeoVol,k,i,j) = dAij(IgeoZb,i,j)*exp(-dAij(cumMu,i,j));
                }
        }
    }
}

void projectXrayGridVolume(
    MultidimArray<double> &muVol,                  // Volume
    XRayPSF &psf,                   // Basis
    MultidimArray<double> &IgeoVol,    /// Vol with the Igeometrical distribution along specimen volume
    Projection       &proj,               // Projection
    MultidimArray<double> *projNorm,     // Projection of a unitary volume
    int              FORW,                // 1 if we are projecting a volume
    //   norm_proj is calculated
    // 0 if we are backprojecting
    //   norm_proj must be valid
    ThreadManager * thMgr)
{
    XrayThreadArgument projThrData;
    Barrier barrier(thMgr->threads);

    projThrData.psf = &psf;
    projThrData.muVol = &muVol;
    projThrData.IgeoVol = &IgeoVol;
    projThrData.projOut = &proj;
    projThrData.projNorm = projNorm;
    projThrData.forw = FORW;
    projThrData.phantomSlabIdx = NULL;
    projThrData.psfSlicesIdx = NULL;

    //    //Create the job handler to distribute thread jobs
    //    size_t blockSize, numberOfJobs = psfSlicesIdx.size() ;
    //    blockSize = 1;
    //    ThreadTaskDistributor td(numberOfJobs, blockSize);
    //    projThrData.td = &td;
    //    projThrData.barrier = &barrier;

    //the really really final project routine, I swear by Snoopy.
    //    project_xr(psf,side.rotPhantomVol,P, idxSlice);

    thMgr->run(projectXraySimpleGridThread,(void*) &projThrData);

}

void projectXraySimpleGridThread(ThreadArgument &thArg)
{

    int threadId = thArg.thread_id;

    XrayThreadArgument *dataThread = (XrayThreadArgument*) thArg.data;
    const XRayPSF &psf = *(dataThread->psf);
    MultidimArray<double> *muVol =  (dataThread->muVol);
    MultidimArray<double> *IgeoVol =  (dataThread->IgeoVol);
    Projection *proj = (dataThread->projOut);
    MultidimArray<double> &projNorm = *(dataThread->projNorm);
    int forw = dataThread->forw;

    projectXraySimpleGrid(muVol, psf, IgeoVol, proj, projNorm , forw,
                          threadId, thArg.threads);




}

void projectXraySimpleGrid(MultidimArray<double> *vol, const XRayPSF &psf,
                           MultidimArray<double> *IgeoVol,
                           Projection *proj, MultidimArray<double> &projNorm , int FORW,
                           int threadId, int numthreads)
{

    MultidimArray<double> psfVol;
    psfVol.alias(psf.psfVol());


    Matrix1D<double> zero(3);                // Origin (0,0,0)
    Matrix1D<double> prjPix(3);       // Position of the pixel within the projection
    Matrix1D<double> prjX(3);         // Coordinate: Projection of the
    Matrix1D<double> prjY(3);         // 3 grid vectors
    Matrix1D<double> prjZ(3);
    Matrix1D<double> prjOrigin(3);    // Coordinate: Where in the 2D
    // projection plane the origin of
    // the grid projects
    Matrix1D<double> prjDir(3);       // Projection direction

    Matrix1D<double> actprj(3);       // Coord: Actual position inside
    // the projection plane
    Matrix1D<double> beginZ(3);       // Coord: Plane coordinates of the
    // projection of the 3D point
    // (z0,YY(lowest),XX(lowest))
    Matrix1D<double> univ_beginY(3);  // Coord: coordinates of the
    // grid point
    // (z0,y0,XX(lowest))
    Matrix1D<double> univ_beginZ(3);  // Coord: coordinates of the
    // grid point
    // (z0,YY(lowest),XX(lowest))
    Matrix1D<double> beginY(3);       // Coord: Plane coordinates of the
    // projection of the 3D point
    // (z0,y0,XX(lowest))
    double XX_footprint_size;                // The footprint is supposed
    double YY_footprint_size;                // to be defined between
    // (-vmax,+vmax) in the Y axis,
    // and (-umax,+umax) in the X axis
    // This footprint size is the
    // R2 vector (umax,vmax).

    int XX_corner2, XX_corner1;              // Coord: Corners of the
    int YY_corner2, YY_corner1;              // footprint when it is projected
    // onto the projection plane
    // directions respectively
    // inside the blobprint
    double        vol_corr;                  // Correction for a volume element
    int           N_eq;                      // Number of equations in which
    // a blob is involved
    int           i, j, k;                   // volume element indexes

    // Project grid axis ....................................................
    // These vectors ((1,0,0),(0,1,0),...) are referred to the grid
    // coord. system and the result is a 2D vector in the image plane
    // The unit size within the image plane is 1, ie, the same as for
    // the universal grid.
    // actprj is reused with a different purpose
    //    VECTOR_R3(actprj, 1, 0, 0);
    //    Uproject_to_plane(actprj, proj->euler, prjX);
    //    VECTOR_R3(actprj, 0, 1, 0);
    //    Uproject_to_plane(actprj, proj->euler, prjY);
    //    VECTOR_R3(actprj, 0, 0, 1);
    //    Uproject_to_plane(actprj, proj->euler, prjZ);
    proj->euler.getCol(0, prjX);
    proj->euler.getCol(1, prjY);
    proj->euler.getCol(2, prjZ);

    // This time the origin of the grid is in the universal coord system
    //    Uproject_to_plane(grid->origin, proj->euler, prjOrigin);
    prjOrigin.initZeros();

    // Get the projection direction .........................................
    proj->euler.getRow(2, prjDir);

    // Footprint size .......................................................
    // The following vectors are integer valued vectors, but they are
    // stored as real ones to make easier operations with other vectors
    XX_footprint_size = FINISHINGX(psfVol) + XMIPP_EQUAL_ACCURACY;
    YY_footprint_size = FINISHINGY(psfVol) + XMIPP_EQUAL_ACCURACY;

    // Project the whole grid ...............................................
    // Corner of the plane defined by Z. These coordinates try to be within
    // the valid indexes of the projection (defined between STARTING and
    // FINISHING values, but it could be that they may lie outside.
    // These coordinates need not to be integer, in general, they are
    // real vectors.
    // The vectors returned by the projection routines are R3 but we
    // are only interested in their first 2 components, ie, in the
    // in-plane components

    // This type conversion gives more speed

    int ZZ_lowest = STARTINGZ(*vol);
    if( threadId != -1 )
        ZZ_lowest += threadId;
    int YY_lowest = STARTINGY(*vol);
    int XX_lowest = STARTINGX(*vol);
    int ZZ_highest = FINISHINGZ(*vol);
    int YY_highest = FINISHINGY(*vol);
    int XX_highest = FINISHINGX(*vol);

    beginZ = (double)XX_lowest * prjX + (double)YY_lowest * prjY + (double)ZZ_lowest * prjZ + prjOrigin;


#ifdef DEBUG_LITTLE

    int condition;
    condition = threadId==1;
    if (condition || numthreads==1)
    {
        std::cout << "Equation mode " << eq_mode << std::endl;
        std::cout << "Footprint size " << YY_footprint_size << "x"
        << XX_footprint_size << std::endl;
        std::cout << "rot=" << proj->rot() << " tilt=" << proj->tilt()
        << " psi=" << proj->psi() << std::endl;
        std::cout << "Euler matrix " << proj->euler;
        std::cout << "Projection direction " << prjDir << std::endl;
        std::cout << *grid;
        std::cout << "image limits (" << x0 << "," << y0 << ") (" << xF << ","
        << yF << ")\n";
        std::cout << "prjX           " << prjX.transpose()      << std::endl;
        std::cout << "prjY           " << prjY.transpose()      << std::endl;
        std::cout << "prjZ           " << prjZ.transpose()      << std::endl;
        std::cout << "prjOrigin      " << prjOrigin.transpose() << std::endl;
        std::cout << "beginZ(coord)  " << beginZ.transpose()    << std::endl;
        std::cout << "lowest         " << XX_lowest  << " " << YY_lowest
        << " " << XX_lowest  << std::endl;
        std::cout << "highest        " << XX_highest << " " << YY_highest
        << " " << XX_highest << std::endl;
        std::cout << "Stats of input basis volume ";
        (*vol)().printStats();
        std::cout << std::endl;
        std::cout.flush();
    }
#endif

    for (k = ZZ_lowest; k <= ZZ_highest; k += numthreads)
    {
        // Corner of the row defined by Y
        beginY = beginZ;
        for (i = YY_lowest; i <= YY_highest; i++)
        {
            // First point in the row
            actprj = beginY;
            for (j = XX_lowest; j <= XX_highest; j++)
            {
                // Ray length interesting
                bool ray_length_interesting = true;

                //                // Points out of 3DPSF
                //                ray_length_interesting = (ABS(zCenterDist) <= ZZ_footprint_size);
                //                // There is still missing the shift of the volume from focal plane

                if (ray_length_interesting)
                {
                    // Be careful that you cannot skip any basis, although its
                    // value be 0, because it is useful for norm_proj
#ifdef DEBUG
                    condition = threadId == 1;
                    if (condition)
                    {
                        printf("\nProjecting grid coord (%d,%d,%d) ", j, i, k);
                        std::cout << "Vol there = " << A3D_ELEM(*vol, k, i, j) << std::endl;
                        printf(" 3D universal position (%f,%f,%f) \n",
                               XX(univ_position), YY(univ_position), ZZ(univ_position));
                        std::cout << " Center of the basis proj (2D) " << XX(actprj) << "," << YY(actprj) << std::endl;
                        Matrix1D<double> aux;
                        Uproject_to_plane(univ_position, proj->euler, aux);
                        std::cout << " Center of the basis proj (more accurate) " << aux.transpose() << std::endl;
                    }
#endif

                    // Search for integer corners for this basis
                    XX_corner1 = CEIL(XMIPP_MAX(STARTINGX(IMGMATRIX(*proj)), XX(actprj) - XX_footprint_size));
                    YY_corner1 = CEIL(XMIPP_MAX(STARTINGY(IMGMATRIX(*proj)), YY(actprj) - YY_footprint_size));
                    XX_corner2 = FLOOR(XMIPP_MIN(FINISHINGX(IMGMATRIX(*proj)), XX(actprj) + XX_footprint_size));
                    YY_corner2 = FLOOR(XMIPP_MIN(FINISHINGY(IMGMATRIX(*proj)), YY(actprj) + YY_footprint_size));

#ifdef DEBUG

                    if (condition)
                    {
                        std::cout << "Clipped and rounded Corner 1 " << XX_corner1
                        << " " << YY_corner1 << " " << std::endl;
                        std::cout << "Clipped and rounded Corner 2 " << XX_corner2
                        << " " << YY_corner2 << " " << std::endl;
                    }
#endif

                    // Check if the voxel falls outside the projection plane
                    // COSS: I have changed here
                    if (XX_corner1 <= XX_corner2 && YY_corner1 <= YY_corner2
                        && INSIDEZ(*IgeoVol,ZZ(actprj)))
                    {

                        double IgeoValue = A3D_ELEM(*IgeoVol, (int)ZZ(actprj), (int)YY(actprj), (int)XX(actprj));

                        if (!FORW)
                            vol_corr = 0;

                        N_eq = 0;
                        for (int y = YY_corner1; y <= YY_corner2; y++)
                        {
                            for (int x = XX_corner1; x <= XX_corner2; x++)
                            {
                                //                                if (!((mask != NULL) && A2D_ELEM(*mask,y,x)<0.5))
                                {
#ifdef DEBUG
                                    if (condition)
                                    {
                                        std::cout << "Position in projection (" << x << ","
                                        << y << ") ";
                                        double y, x;
                                        if (basis->type == Basis::blobs)
                                        {
                                            std::cout << "in footprint ("
                                            << foot_U << "," << foot_V << ")";
                                            IMG2OVER(basis->blobprint, foot_V, foot_U, y, x);
                                            std::cout << " (d= " << sqrt(y*y + x*x) << ") ";
                                            fflush(stdout);
                                        }
                                    }
#endif

                                    double a, a2;

                                    a = A3D_ELEM(psfVol, (int)ZZ(actprj), y, x) * IgeoValue;
                                    a2 = a * a;

#ifdef DEBUG

                                    if (condition)
                                        std::cout << "=" << a << " , " << a2;
#endif

                                    if (FORW)
                                    {
                                        /// ARTK Equation mode
                                        IMGPIXEL(*proj, y, x) += A3D_ELEM(*vol, k, i, j) * a;
                                        A2D_ELEM(projNorm , y, x) += a2;

#ifdef DEBUG

                                        if (condition)
                                        {
                                            std::cout << " proj= " << IMGPIXEL(*proj, y, x)
                                            << " norm_proj=" << A2D_ELEM(projNorm , y, x) << std::endl;
                                            std::cout.flush();
                                        }
#endif

                                    }
                                    else
                                    {
                                        vol_corr += A2D_ELEM(projNorm , y, x) * a;
                                        if (a != 0)
                                            N_eq++;
#ifdef DEBUG

                                        if (condition)
                                        {
                                            std::cout << " corr_img= " << A2D_ELEM(projNorm , y, x)
                                            << " correction=" << vol_corr << std::endl;
                                            std::cout.flush();
                                        }
#endif

                                    }
                                }// close possible mask constraint
                            }
                        } // Project this basis

                        if (!FORW)
                        {
                            A3D_ELEM(*vol, k, i, j) += vol_corr;

#ifdef DEBUG

                            if (condition)
                            {
                                printf("\nFinal value at (%d,%d,%d) ", j, i, k);
                                std::cout << " = " << A3D_ELEM(*vol, k, i, j) << std::endl;
                                std::cout.flush();
                            }
#endif

                        }
                    } // If not collapsed
                    //                    number_of_basis++;
                } // If interesting

                // Prepare for next iteration
                V3_PLUS_V3(actprj, actprj, prjX);
            }
            V3_PLUS_V3(beginY, beginY, prjY);
        }
        V3_PLUS_V3(beginZ, beginZ, prjZ * numthreads );
    }
}


