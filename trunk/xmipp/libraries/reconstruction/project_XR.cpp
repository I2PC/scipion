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
#include "project_XR.h"

void ProgXrayProject::defineParams()
{
    addUsageLine("Generate projections as in a X-ray microscope from a 3D Xmipp volume.");
    addSeeAlsoLine("xray_psf_create, xray_import, phantom_create");
    //Params
    projParam.defineParams(this); // Projection parameters
    addParamsLine("== Xray options == ");
    addParamsLine(" [-s <Ts>]     : Sampling rate of the volume to be projected (nm).");
    addParamsLine("                              : If empty, same value as X-Y plane sampling from PSF.");
    addParamsLine("alias --sampling_rate;");
    addParamsLine("[--psf <psf_param_file=\"\">] : XRay-Microscope parameters file. If not set, then default parameters are chosen.");
    addParamsLine("[--thr <threads=1>]           : Number of concurrent threads.");
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
    addExampleLine("http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/tomoProjection.param",false);
    addExampleLine(" ",false);
    addExampleLine("In the following link you can find an example of X-ray microscope parameters file:",false);
    addExampleLine(" ",false);
    addExampleLine("http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/xray_psf.xmd",false);
}

/* Read from command line ================================================== */
void ProgXrayProject::readParams()
{
    //    fn_proj_param = getParam("-i");
    //    fn_sel_file   = getParam("-o");
    projParam.readParams(this);

    fn_psf_xr = getParam("--psf");
    dxo  = (checkParam("-s"))? getDoubleParam("-s")*1e-9 : -1 ;
    nThr = getIntParam("--thr");

    psf.read(fn_psf_xr);
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

    psf.calculateParams(dxo);
    phantom.read(projParam);
    //Correct the rotation axis displacement in projectionParams from pixels to meters
    projParam.raxis *= psf.dxo;

    // Threads stuff
    mainDataThread = new XrayThread;

    mainDataThread->psf= &psf;
    mainDataThread->vol = &phantom.rotVol;
    mainDataThread->imOut = &proj;

    //Create the job handler to distribute thread jobs

    size_t blockSize, numberOfJobs = ZSIZE(MULTIDIM_ARRAY_BASE(phantom.iniVol));
    blockSize = (nThr == 1) ? numberOfJobs : numberOfJobs/nThr/2;
    td = new ThreadTaskDistributor(numberOfJobs, blockSize);
    mainDataThread->td = td;

    barrier = new Barrier(nThr-1);

    //Create threads to start working
    thMgr = new ThreadManager(nThr,(void*) mainDataThread);

}

void ProgXrayProject::postRun()
{
    //Terminate threads and free memory
    delete td;
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
        td->clear();

        // Really project ....................................................
        if (!projParam.only_create_angles)
            XrayProjectVolumeOffCentered(phantom, psf, proj,projParam.proj_Ydim, projParam.proj_Xdim, idx);

        // Add noise in angles and voxels ....................................
        proj.getEulerAngles(tRot, tTilt,tPsi);

        rot  = tRot  + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);
        tilt = tTilt + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);
        psi  = tPsi  + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);

        proj.setEulerAngles(rot,tilt,psi);

        size_t objId = projMD.addObject();
        if (!projParam.only_create_angles)
        {
            proj.write(fn_proj, ALL_IMAGES, !projParam.singleProjection, WRITE_REPLACE);
            projMD.setValue(MDL_IMAGE,fn_proj,objId);
        }

        projMD.setValue(MDL_ANGLEROT,rot,objId);
        projMD.setValue(MDL_ANGLETILT,tilt,objId);
        projMD.setValue(MDL_ANGLEPSI,psi,objId);
        projMD.setValue(MDL_ANGLEROT2,tRot,objId);
        projMD.setValue(MDL_ANGLETILT2,tTilt,objId);
        projMD.setValue(MDL_ANGLEPSI2,tPsi,objId);
        projMD.setValue(MDL_SHIFTX,shiftX,objId);
        projMD.setValue(MDL_SHIFTY,shiftY,objId);

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
        projMD.write(projParam.fnRoot + ".sel");
    }

    postRun();

    return;
}

void XrayProjectVolumeOffCentered(XrayProjPhantom &phantom, XRayPSF &psf, Projection &P,
                                  int Ydim, int Xdim, int idxSlice)
{

    int iniXdim, iniYdim, iniZdim, newXdim, newYdim, newZdim, rotXdim, rotYdim, rotZdim;
    int zinit, zend, yinit, yend, xinit, xend;

    iniXdim = XSIZE(MULTIDIM_ARRAY_BASE(phantom.iniVol));
    iniYdim = YSIZE(MULTIDIM_ARRAY_BASE(phantom.iniVol));
    iniZdim = ZSIZE(MULTIDIM_ARRAY_BASE(phantom.iniVol));

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
        rotYdim = iniYdim;
        rotZdim = iniZdim;
    }

    newXdim = rotXdim + 2*ABS(XX(offsetNV));
    newYdim = iniYdim + 2*ABS(YY(offsetNV));
    newZdim = rotZdim;

    // We set the dimensions only to obtains the values of starting X,Y,Z
    //    phantom.rotVol.setDimensions(rotXdim, rotYdim, rotZdim, 1);
    phantom.rotVol.setDimensions(newXdim, newYdim, newZdim, 1);

    phantom.rotVol.setXmippOrigin();

    zinit = STARTINGZ((phantom.rotVol));
    zend = zinit + rotZdim - 1;

    yinit = STARTINGY((phantom.rotVol));
    if (YY(offsetNV) > 0)
        yinit -= 2 * YY(offsetNV);
    yend = yinit + newYdim -1;

    xinit = STARTINGX((phantom.rotVol));
    if (XX(offsetNV) > 0)
        xinit -= 2 * XX(offsetNV);
    xend = xinit + newXdim - 1;

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
    phantom.iniVol.data->setXmippOrigin();

    applyGeometry(1, phantom.rotVol, MULTIDIM_ARRAY_GENERIC(phantom.iniVol), T*R, IS_NOT_INV, DONT_WRAP, outside);
//
//    Image<double> tempvol;
//
//    tempvol.data.alias(phantom.rotVol);
//
//    tempvol.write("rotvol_2.vol");
//
//    exit(0);

    //    Euler_rotate(MULTIDIM_ARRAY_GENERIC(phantom.iniVol), P.rot(), P.tilt(), P.psi(),phantom.rotVol);

    // Correct the shift position due to tilt axis is out of optical axis
    //    phantom.rotVol.selfWindow(zinit, yinit, xinit, zend, yend, xend);

    psf.adjustParam(phantom.rotVol);

    //the really really final project routine, I swear by Snoopy.
    //    project_xr(psf,side.rotPhantomVol,P, idxSlice);
    thMgr->run(threadXrayProject);

    int outXDim = XMIPP_MIN(Xdim,iniXdim);
    int outYDim = XMIPP_MIN(Ydim,iniYdim);

    P().selfWindow(-ROUND(outYDim/2),
                   -ROUND(outXDim/2),
                   -ROUND(outYDim/2) + outYDim -1,
                   -ROUND(outXDim/2) + outXDim -1);

}

/* Produce Side Information ================================================ */
void XrayProjPhantom::read(
    const ParametersProjectionTomography &prm)
{
    //    iniVol.readMapped(prm.fnPhantom);
    iniVol.read(prm.fnPhantom);
    MULTIDIM_ARRAY_GENERIC(iniVol).setXmippOrigin();
    //    rotVol.resizeNoCopy(MULTIDIM_ARRAY(iniVol));
}

/// Generate an X-ray microscope projection for volume vol using the microscope configuration psf
void threadXrayProject(ThreadArgument &thArg)
{

    int thread_id = thArg.thread_id;

    XrayThread *dataThread = (XrayThread*) thArg.workClass;
    XRayPSF psf = *(dataThread->psf);
    MultidimArray<double> &vol =  *(dataThread->vol);
    Image<double> &imOutGlobal = *(dataThread->imOut);
    ParallelTaskDistributor * td = dataThread->td;

    size_t first = 0, last = 0, priorLast = 0;

    if (thread_id == 0)
    {
        vol.setXmippOrigin();
        MULTIDIM_ARRAY(imOutGlobal).resizeNoCopy(psf.Niy, psf.Nix);
        MULTIDIM_ARRAY(imOutGlobal).initZeros();
        MULTIDIM_ARRAY(imOutGlobal).setXmippOrigin();
    }

    barrier->wait();

    MultidimArray<double> imOut(psf.Niy, psf.Nix);
    imOut.setXmippOrigin();

    MultidimArray<double> imTemp(psf.Noy, psf.Nox),intExp(psf.Noy, psf.Nox),imTempSc(imOut),*imTempP;
    intExp.setXmippOrigin();
    imTemp.setXmippOrigin();
    imTempSc.setXmippOrigin();

    // Parallel thread jobs
    while (td->getTasks(first, last))
    {
        // Previous addition for threads calculating intermediate slices
        for (int k = STARTINGZ(vol) + (int)priorLast; k <= STARTINGZ(vol) + (int)first - 1 ; k++)
        {
            FOR_ALL_ELEMENTS_IN_ARRAY2D(intExp)
            A2D_ELEM(intExp,i,j) += A3D_ELEM(vol,k,i,j);
        }

        //#define DEBUG
#ifdef DEBUG
        Image<double> _Im;
#endif

        int kini = STARTINGZ(vol) + first;
        int kend = STARTINGZ(vol) + last;

        for (int k = kini; k <= kend; k++)
        {
            FOR_ALL_ELEMENTS_IN_ARRAY2D(intExp)
            {
                double tempValue = A3D_ELEM(vol,k,i,j);
                A2D_ELEM(intExp,i,j) += tempValue;
                A2D_ELEM(imTemp,i,j) = (exp(-A2D_ELEM(intExp,i,j)*psf.dzo))*tempValue*psf.dzo;
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

            psf.applyOTF(*imTempP, imOutGlobal.Zoff()- k*psf.dzo);

#ifdef DEBUG

            _Im.write("projectXR-imTempEsc_after.spi");
#endif

            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
            dAij(imOut,i,j) += dAij(*imTempP,i,j);

#ifdef DEBUG

            _Im().alias(imOut);
            _Im.write("projectXR-imout.spi");
#endif

        }
        priorLast = last + 1;
        //        std::cerr << "th" << thread_id << ": Finished work from " << first << " to " << last <<std::endl;
    }

    //Lock to update the total addition
    mutex.lock();
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imOut)
    dAij(MULTIDIM_ARRAY(imOutGlobal),i,j) += dAij(imOut,i,j);
    mutex.unlock();

    barrier->wait();

    if (thread_id==0)
    {
        MULTIDIM_ARRAY(imOutGlobal) = 1- MULTIDIM_ARRAY(imOutGlobal);

        switch (psf.AdjustType)
        {
        case PSFXR_INT:
            selfScaleToSize(LINEAR,imOutGlobal(), psf.Nox, psf.Noy);
            break;
        case PSFXR_ZPAD:
            MULTIDIM_ARRAY(imOutGlobal).selfWindow(-ROUND(psf.Noy/2)+1,-ROUND(psf.Nox/2)+1,ROUND(psf.Noy/2)-1,ROUND(psf.Nox/2)-1);
            break;
        }

        MULTIDIM_ARRAY(imOutGlobal).setXmippOrigin();
    }
    //    std::cerr << "th" << thread_id << ": Thread Finished" <<std::endl;
}
