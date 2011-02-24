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


void ProgXrayProject::run()
{
    Projection   proj;
    MetaData     projMD;

    randomize_random_generator();

    psf.calculateParams(dxo);

    XrayProjPhantom phantom;
    phantom.read(projParam);


    // Project

    // Threads stuff
    XrayThread *dataThread = new XrayThread;

    dataThread->psf= &psf;
    dataThread->vol = &phantom.rotVol;
    dataThread->imOut = &proj;

    longint blockSize, numberOfJobs = ZSIZE(MULTIDIM_ARRAY(phantom.iniVol));
    numberOfThreads = psf.nThr;

    blockSize = (numberOfThreads == 1) ? numberOfJobs : numberOfJobs/numberOfThreads/2;

    //Create the job handler to distribute thread jobs
    td = new ThreadTaskDistributor(numberOfJobs, blockSize);
    barrier = new Barrier(numberOfThreads-1);

    //Create threads to start working
    thMgr = new ThreadManager(numberOfThreads,(void*) dataThread);

    int expectedNumProjs = FLOOR((projParam.tiltF-projParam.tilt0)/projParam.tiltStep);
    int numProjs=0;

    std::cerr << "Projecting ...\n";
    if (!(projParam.show_angles))
        init_progress_bar(expectedNumProjs);

    projMD.setComment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
    double tRot,tTilt,tPsi,rot,tilt,psi;
    FileName fn_proj;  // Projection name
    size_t idx = FIRST_IMAGE;

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
        projMD.setValue(MDL_ANGLEROT,tRot,objId);
        projMD.setValue(MDL_ANGLETILT,tTilt,objId);
        projMD.setValue(MDL_ANGLEPSI,tPsi,objId);
        projMD.setValue(MDL_ANGLEROT2,rot,objId);
        projMD.setValue(MDL_ANGLETILT2,tilt,objId);
        projMD.setValue(MDL_ANGLEPSI2,psi,objId);
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
        else if ((expectedNumProjs % XMIPP_MAX(1, numProjs / 60))  == 0)
            progress_bar(numProjs);

        numProjs++;
        idx++;
    }
    if (!(projParam.show_angles))
        progress_bar(expectedNumProjs);

    // Save metadata file with angles and shift info
    if (!projParam.singleProjection)
        projMD.write(projParam.fnRoot + ".sel");

    //Terminate threads and free memory
    delete td;
    delete thMgr;
    delete barrier;

    return;
}

/* Produce Side Information ================================================ */
void XrayProjPhantom::read(
    const ParametersProjectionTomography &prm)
{
    iniVol.read(prm.fnPhantom);
    MULTIDIM_ARRAY(iniVol).setXmippOrigin();
    rotVol.resizeNoCopy(MULTIDIM_ARRAY(iniVol));
}

void XrayProjectVolumeOffCentered(XrayProjPhantom &phantom, XRayPSF &psf, Projection &P,
                                  int Ydim, int Xdim, int idxSlice)
{

    int iniXdim, iniYdim, iniZdim, newXdim, newYdim;
    int xOffsetN, yOffsetN, zinit, zend, yinit, yend, xinit, xend;

    iniXdim = XSIZE(MULTIDIM_ARRAY(phantom.iniVol));
    iniYdim = YSIZE(MULTIDIM_ARRAY(phantom.iniVol));
    iniZdim = ZSIZE(MULTIDIM_ARRAY(phantom.iniVol));

    // Projection offset in pixels
    xOffsetN = P.Xoff()/psf.dxo;
    yOffsetN = P.Yoff()/psf.dxo;

    newXdim = iniXdim + 2*ABS(xOffsetN);
    newYdim = iniYdim + 2*ABS(yOffsetN);

    zinit = STARTINGZ(MULTIDIM_ARRAY(phantom.iniVol));
    zend = zinit + iniZdim - 1;

    if (yOffsetN<=0)
        yinit = STARTINGY(MULTIDIM_ARRAY(phantom.iniVol));
    else
        yinit = STARTINGY(MULTIDIM_ARRAY(phantom.iniVol)) - 2 * yOffsetN;

    yend = yinit + newYdim -1;

    if (xOffsetN<=0)
        xinit = STARTINGX(MULTIDIM_ARRAY(phantom.iniVol));
    else
        xinit = STARTINGX(MULTIDIM_ARRAY(phantom.iniVol)) - 2 * xOffsetN;

    xend = xinit + newXdim - 1;


    if (psf.verbose > 1)
    {
        int finalXdim = xend - xinit + 1;
        int finalYdim = yend - yinit + 1;

        if (finalXdim != iniXdim || finalYdim != iniYdim)
        {
            std::cout << std::endl;
            std::cout << "--------------------------------" << std::endl;
            std::cout << "XrayProject::Volume_offCentered:" << std::endl;
            std::cout << "--------------------------------" << std::endl;
            std::cout << "(X,Y,Z) shifts = (" << P.Xoff()*1e6 << "," << P.Yoff()*1e6 << ","
            << P.Zoff()*1e6 << ") um" << std::endl;
            std::cout << "Image resize (Nx,Ny): (" << iniXdim << "," << iniYdim << ") --> ("
            << finalXdim << "," << finalYdim << ") " << std::endl;
        }

#ifdef DEBUG

        std::cout <<"yoffsetN "<< yOffsetN <<std::endl;
        std::cout <<"xoffsetN "<< xOffsetN <<std::endl;
        std::cout <<"yinit    " << yinit  <<std::endl;
        std::cout <<"yend     "    << yend  <<std::endl;
        std::cout <<"xinit    "   << xinit  <<std::endl;
        std::cout <<"xend     "    << xend  <<std::endl;
        std::cout <<"zinit    "   << zinit  <<std::endl;
        std::cout <<"zend     "    << zend  <<std::endl;
#endif

    }

    // Rotate volume ....................................................
    //    applyGeometry(LINEAR,volTemp(), V, Euler_rotation3DMatrix(rot, tilt, psi), IS_NOT_INV, DONT_WRAP);

    phantom.rotVol.resizeNoCopy(iniZdim,iniYdim,iniXdim);
    //    side.rotPhantomVol().zinit = zinit;
    //    side.rotPhantomVol().yinit = yinit;
    //    side.rotPhantomVol().xinit = xinit;
    phantom.rotVol.setXmippOrigin();

    Euler_rotate(MULTIDIM_ARRAY(phantom.iniVol), P.rot(), P.tilt(), P.psi(),phantom.rotVol);

    // Correct the shift position due to tilt axis is out of optical axis
    phantom.rotVol.selfWindow(zinit, yinit, xinit, zend, yend, xend);

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

/// Generate an X-ray microscope projection for volume vol using the microscope configuration psf
void threadXrayProject(ThreadArgument &thArg)
{
    int thread_id = thArg.thread_id;

    XrayThread *dataThread = (XrayThread*) thArg.workClass;
    XRayPSF psf = *(dataThread->psf);
    MultidimArray<double> &vol =  *(dataThread->vol);
    Image<double> &imOutGlobal = *(dataThread->imOut);

    long long int first = -1, last = -1, priorLast = -1;

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
        //        std::cerr << "th" << thread_id << ": working from " << first << " to " << last <<std::endl;

        //        if (numberOfThreads == 1)
        //            first = last -50;
        //        if (first>300)

        {
            // Previous addition for threads calculating intermediate slices
            for (int k = STARTINGZ(vol) + priorLast + 1; k <= STARTINGZ(vol) + first - 1 ; k++)
            {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(intExp)
                A2D_ELEM(intExp,i,j) += A3D_ELEM(vol,k,i,j);
            }

            //#define DEBUG
#ifdef DEBUG
            Image<double> _Im;
#endif

            for (int k = STARTINGZ(vol) + first; k <= STARTINGZ(vol) + last; k++)
            {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(intExp)
                {
                    double tempValue = A3D_ELEM(vol,k,i,j);
                    A2D_ELEM(intExp,i,j) += tempValue;
                    A2D_ELEM(imTemp,i,j) = (exp(-A2D_ELEM(intExp,i,j)*psf.dzo))*tempValue*psf.dzo;
                }
#ifdef DEBUG
                //                _Im().resizeNoCopy(intExp);
                //                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(intExp)
                //                dAij(_Im(),i,j) = dAij(intExp,i,j);
                _Im().alias(intExp);
                _Im.write("psfxr-intExp.spi");
                //                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imTemp)
                //                dAij(_Im(),i,j) = dAij(imTemp,i,j);
                _Im().alias(imTemp);
                _Im.write("psfxr-imTemp.spi");
#endif

                //                psf.Z = psf.Zo + psf.DeltaZo - k*psf.dzo + imOutGlobal.Zoff();

                switch (psf.AdjustType)
                {
                case PSFXR_INT:
                    imTempP = &imTempSc;
                    scaleToSize(LINEAR,*imTempP,imTemp,psf.Nix,psf.Niy);
                    //          imTemp.scaleToSize(psf.Niy, psf.Nix, *imTempP);
                    break;

                case PSFXR_STD:
                    imTempP = &imTemp;
                    break;

                case PSFXR_ZPAD:
                    //    (*imTempSc).resize(imTemp);
                    imTempSc = imTemp;
                    imTempSc.selfWindow(-ROUND(psf.Niy/2)+1,-ROUND(psf.Nix/2)+1,ROUND(psf.Niy/2)-1,ROUND(psf.Nix/2)-1);
                    imTempP = &imTempSc;
                    break;
                }

#ifdef DEBUG
                //                _Im().resize(*imTempP);
                //                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
                //                dAij(_Im(),i,j) = dAij(*imTempP,i,j);
                _Im().alias(*imTempP);
                _Im.write("psfxr-imTempEsc.spi");
#endif

                psf.applyOTF(*imTempP, imOutGlobal.Zoff()- k*psf.dzo);

#ifdef DEBUG

                //                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
                //                dAij(_Im(),i,j) = dAij(*imTempP,i,j);
                _Im.write("psfxr-imTempEsc2.spi");
#endif

                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
                dAij(imOut,i,j) += dAij(*imTempP,i,j);

#ifdef DEBUG

                _Im().alias(imOut);
                _Im.write("psfxr-imout.spi");
#endif

            }
        }
        priorLast = last;
        //        std::cerr << "th" << thread_id << ": Finished work from " << first << " to " << last <<std::endl;
    }

    //Lock to update the total summatory
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
