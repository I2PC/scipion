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
#include <data/psf_xr.h>

void ProgProjectXR::defineParams()
{
    addUsageLine("Generate projections as in a X-ray microscope from a 3D Xmipp volume.");
    addParamsLine(" -i <Proj_param_file>   :MetaData file with projection parameters.");
    addParamsLine("                         :Check the manual for a description of the parameters.");
    addParamsLine(" alias --input;");
    addParamsLine("[-psf <psf_param_file=\"\">] : XRay-Microscope parameters file. If not set, then default parameters are chosen.");
    addParamsLine("[-o <sel_file=\"\">]         : Output Metadata file with all the generated projections.");
    addParamsLine("[-show_angles]          : Print angles value for each projection.");
    addParamsLine("[-only_create_angles]   : Projections are not calculated, only the angles values.");
    addParamsLine("[-thr <threads=1>]      : Number of concurrent threads.");
}

/* Read from command line ================================================== */
void ProgProjectXR::readParams()
{
    fn_proj_param = getParam("-i");
    fn_sel_file   = getParam("-o");
    fn_psf_xr   = getParam("-psf");

    only_create_angles = checkParam("-only_create_angles");
    verbose = getParam("-v");
    nThr = getIntParam("-thr");
    tell=0;
    if (checkParam("-show_angles"))
        tell |= TELL_SHOW_ANGLES;
}


void ProgProjectXR::run()
{
    Projection         proj;
    MetaData           SF;

    randomize_random_generator();

    // Read Microscope optics parameters and produce side information
    XRayPSF psf;
    psf.verbose = verbose;
    psf.nThr = nThr;
    psf.read(fn_psf_xr);
    psf.produceSideInfo();
    if(psf.verbose)
        psf.show();


    // Read projection parameters and produce side information
    ParametersProjectionXR proj_prm;
    proj_prm.read(fn_proj_param);
    proj_prm.tell=tell;
    PROJECT_XR_Side_Info side;
    side.produce_Side_Info(proj_prm);


    //    psf.adjustParam(side.phantomVol);

    // Project
    int ProjNo = 0;
    if (!only_create_angles)
    {
        // Really project
        ProjNo = PROJECT_XR_Effectively_project(proj_prm, side,
                                                proj, psf, SF);
        // Save SelFile
        if (fn_sel_file != "")
            SF.write(fn_sel_file);
    }
    else
    {
        side.DF.write("/dev/stdout");
    }
    return;
}

void ParametersProjectionXR::calculateProjectionAngles(Projection &P, double angle, double inplaneRot,
        const Matrix1D<double> &rinplane)
{
    // double axisRot, double axisTilt,
    //                          const Matrix1D<double> &raxis,

    // Find Euler rotation matrix
    Matrix1D<double> axis;
    Euler_direction(axisRot,axisTilt,0,axis);
    Matrix2D<double> Raxis, Rinplane;
    rotation3DMatrix(angle,axis,Raxis,false);
    rotation3DMatrix(inplaneRot,'Z',Rinplane,false);
    double rot, tilt, psi;
    Euler_matrix2angles(Rinplane*Raxis, rot, tilt, psi);
    P.set_angles(rot, tilt, psi);

    // Find displacement because of axis offset and inplane shift
    Matrix1D<double> roffset=Rinplane*(raxis-Raxis*raxis)+rinplane;

    P.setShifts(XX(roffset), YY(roffset), ZZ(roffset));
}

/* Produce Side Information ================================================ */
void PROJECT_XR_Side_Info::produce_Side_Info(
    const ParametersProjectionXR &prm)
{
    phantomVol.read(prm.fnPhantom);
    phantomVol().setXmippOrigin();
    rotPhantomVol = phantomVol;
}

/* Effectively project ===================================================== */
int PROJECT_XR_Effectively_project(ParametersProjectionXR &prm,
                                   PROJECT_XR_Side_Info &side,
                                   Projection &proj,
                                   XRayPSF &psf,
                                   MetaData &SF)
{

    // Threads stuff

    XrayThread *dataThread = new XrayThread;

    dataThread->psf= &psf;
    dataThread->vol = &side.rotPhantomVol;
    dataThread->imOut = &proj;

    longint blockSize, numberOfJobs= side.phantomVol().zdim;
    numberOfThreads = psf.nThr;

    blockSize = (numberOfThreads == 1) ? numberOfJobs : numberOfJobs/numberOfThreads/2;

    //Create the job handler to distribute thread jobs
    td = new ThreadTaskDistributor(numberOfJobs, blockSize);
    barrier = new Barrier(numberOfThreads-1);

    //Create threads to start working
    thMgr = new ThreadManager(numberOfThreads,(void*) dataThread);

    int expectedNumProjs = FLOOR((prm.tiltF-prm.tilt0)/prm.tiltStep);
    int numProjs=0;

    std::cerr << "Projecting ...\n";
    if (!(prm.tell&TELL_SHOW_ANGLES))
        init_progress_bar(expectedNumProjs);

    SF.clear();
    MetaData DF_movements;
    DF_movements.setComment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
    double tRot,tTilt,tPsi,rot,tilt,psi;

    int idx=prm.starting;
    for (double angle=prm.tilt0; angle<=prm.tiltF; angle+=prm.tiltStep)
    {
        FileName fn_proj;              // Projection name
        fn_proj.compose(prm.fnProjectionSeed, idx,
                        prm.fn_projection_extension);

        // Choose Center displacement ........................................
        double shiftX     = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        double shiftY    = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        Matrix1D<double> inPlaneShift(3);
        VECTOR_R3(inPlaneShift,shiftX,shiftY,0);

        prm.calculateProjectionAngles(proj,angle, 0,inPlaneShift);

        //Reset thread task distributor
        td->clear();

        // Really project ....................................................
        project_xr_Volume_offCentered(side, psf, proj,prm.proj_Ydim, prm.proj_Xdim, idx);

        // Add noise in angles and voxels ....................................
        proj.getEulerAngles(tRot, tTilt,tPsi);

        rot  = tRot  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        tilt = tTilt + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        psi  = tPsi  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);

        proj.setEulerAngles(rot,tilt,psi);

        size_t objId=DF_movements.addObject();
        DF_movements.setValue(MDL_ANGLEROT,tRot,objId);
        DF_movements.setValue(MDL_ANGLETILT,tTilt,objId);
        DF_movements.setValue(MDL_ANGLEPSI,tPsi,objId);
        DF_movements.setValue(MDL_ANGLEROT2,rot,objId);
        DF_movements.setValue(MDL_ANGLETILT2,tilt,objId);
        DF_movements.setValue(MDL_ANGLEPSI2,psi,objId);
        DF_movements.setValue(MDL_SHIFTX,shiftX,objId);
        DF_movements.setValue(MDL_SHIFTY,shiftY,objId);

        IMGMATRIX(proj).addNoise(prm.Npixel_avg, prm.Npixel_dev, "gaussian");

        // Save ..............................................................
        if (prm.tell&TELL_SHOW_ANGLES)
            std::cout << idx << "\t" << proj.rot() << "\t"
            << proj.tilt() << "\t" << proj.psi() << std::endl;
        else if ((expectedNumProjs % XMIPP_MAX(1, numProjs / 60))  == 0)
            progress_bar(numProjs);

        proj.write(fn_proj);
        numProjs++;
        idx++;
        objId=SF.addObject();
        SF.setValue(MDL_IMAGE,fn_proj,objId);
        SF.setValue(MDL_ENABLED,1,objId);
    }
    if (!(prm.tell&TELL_SHOW_ANGLES))
        progress_bar(expectedNumProjs);

    DF_movements.write(prm.fnProjectionSeed + "_movements.txt");

    //Terminate threads and free memory
    delete td;
    delete thMgr;
    delete barrier;

    return numProjs;
}

void project_xr_Volume_offCentered(PROJECT_XR_Side_Info &side, XRayPSF &psf, Projection &P,
                                   int Ydim, int Xdim, int idxSlice)
{

    int iniXdim, iniYdim, iniZdim, newXdim, newYdim;
    int xOffsetN, yOffsetN, zinit, zend, yinit, yend, xinit, xend;

    iniXdim = side.phantomVol().xdim;
    iniYdim = side.phantomVol().ydim;
    iniZdim = side.phantomVol().zdim;

    // Projection offset in pixels
    xOffsetN = P.Xoff()/psf.dxo;
    yOffsetN = P.Yoff()/psf.dxo;

    newXdim = iniXdim + 2*ABS(xOffsetN);
    newYdim = iniYdim + 2*ABS(yOffsetN);

    zinit = side.phantomVol().zinit;
    zend = zinit + iniZdim - 1;

    if (yOffsetN<=0)
        yinit = side.phantomVol().yinit ;
    else
        yinit = side.phantomVol().yinit - 2 * yOffsetN;

    yend = yinit + newYdim -1;

    if (xOffsetN<=0)
        xinit = side.phantomVol().xinit;
    else
        xinit = side.phantomVol().xinit - 2 * xOffsetN;

    xend = xinit + newXdim - 1;


    if (psf.verbose)
    {
        std::cout << std::endl;
        std::cout << "project_XR::Volume_offCentered:" << std::endl;
        std::cout << "(X,Y,Z) shifts = (" << P.Xoff()*1e6 << "," << P.Yoff()*1e6 << ","
        << P.Zoff()*1e6 << ") um" << std::endl;
        std::cout << "Image resize (Nx,Ny): (" << iniXdim << "," << iniYdim << ") --> ("
        << xend - xinit +1<< "," << yend - yinit +1 << ") " << std::endl;

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

    side.rotPhantomVol().resize(iniZdim,iniYdim,iniXdim);
    //    side.rotPhantomVol().zinit = zinit;
    //    side.rotPhantomVol().yinit = yinit;
    //    side.rotPhantomVol().xinit = xinit;
    side.rotPhantomVol().setXmippOrigin();

    Euler_rotate(side.phantomVol(), P.rot(), P.tilt(), P.psi(),side.rotPhantomVol());

    // Correct the shift position due to tilt axis is out of optical axis
    side.rotPhantomVol().window(zinit, yinit, xinit, zend, yend, xend);

    psf.adjustParam(side.rotPhantomVol);

    //the really really final project routine, I swear by Snoopy.
    //    project_xr(psf,side.rotPhantomVol,P, idxSlice);
    thMgr->run(thread_project_xr);

    int outXDim = XMIPP_MIN(Xdim,iniXdim);
    int outYDim = XMIPP_MIN(Ydim,iniYdim);

    P().window(-ROUND(outYDim/2),
               -ROUND(outXDim/2),
               -ROUND(outYDim/2) + outYDim -1,
               -ROUND(outXDim/2) + outXDim -1);

}

/// Generate an X-ray microscope projection for volume vol using the microscope configuration psf
void project_xr(XRayPSF &psf, Image<double> &vol, Image<double> &imOut, int idxSlice)
{

    XrayThread *dataThread = new XrayThread;

    dataThread->psf= &psf;
    dataThread->vol = &vol;
    dataThread->imOut = &imOut;

    longint blockSize, numberOfJobs= vol().zdim;
    numberOfThreads = psf.nThr;

    blockSize = (numberOfThreads == 1) ? numberOfJobs : numberOfJobs/numberOfThreads;

    //Create the job handler to distribute jobs
    td = new ThreadTaskDistributor(numberOfJobs, blockSize);
    //    td = new FileTaskDistributor(numberOfJobs,blockSize);
    barrier = new Barrier(numberOfThreads-1);

    //Create threads to start working

    ThreadManager * thMgr = new ThreadManager(numberOfThreads,(void*) dataThread);


    //    if (numberOfThreads==1)
    //    {
    //        ThreadArgument thArg;
    //        thArg.thread_id = 0;
    //        thArg.workClass = dataThread;
    //        thread_project_xr(thArg);
    //    }
    //    else
    thMgr->run(thread_project_xr);


    //Terminate threads and free memory
    delete td;
    delete thMgr;


}


void thread_project_xr(ThreadArgument &thArg)
{

    int thread_id = thArg.thread_id;

    XrayThread *dataThread = (XrayThread*) thArg.workClass;
    XRayPSF psf = *(dataThread->psf);
    Image<double> &vol =  *(dataThread->vol);
    Image<double> &imOutGlobal = *(dataThread->imOut);

    long long int first = -1, last = -1, priorLast = -1;

    if (thread_id == 0)
    {
        vol().setXmippOrigin();
        imOutGlobal().resizeNoCopy(psf.Niy, psf.Nix);
        imOutGlobal().initZeros();
        imOutGlobal().setXmippOrigin();
    }

    barrier->wait();

    Image<double> imOut(psf.Nix,psf.Niy);
    imOut().setXmippOrigin();

    MultidimArray<double> imTemp(psf.Noy, psf.Nox),intExp(psf.Noy, psf.Nox),imTempSc(imOut()),*imTempP;
    intExp.setXmippOrigin();
    imTemp.setXmippOrigin();
    imTempSc.setXmippOrigin();

    // Parallel thread jobs
    while (td->getTasks(first, last))
    {
        std::cerr << "th" << thread_id << ": working from " << first << " to " << last <<std::endl;

        //        if (numberOfThreads == 1)
        //            first = last -50;
        //        if (first>300)

        {
            for (int k=(vol()).zinit + priorLast + 1; k<=(vol()).zinit + first - 1 ; k++)
            {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(intExp)
                intExp(i, j) += vol(k, i, j);
            }


            //#define DEBUG
#ifdef DEBUG
            Image<double> _Im(imOut);
#endif

            for (int k=((vol()).zinit + first); k<=((vol()).zinit + last); k++)
            {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(intExp)
                {
                    intExp(i, j) += vol(k, i, j);
                    imTemp(i, j) = (exp(-intExp(i,j)*psf.dzo))*vol(k,i,j)*psf.dzo;
                    //            imTemp(i, j) = 1./(exp(intExp(i,j)))*vol(k,i,j);
                }
#ifdef DEBUG
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imTemp)
                dAij(_Im(),i,j) = dAij(imTemp,i,j);
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
                    imTempSc.window(-ROUND(psf.Niy/2)+1,-ROUND(psf.Nix/2)+1,ROUND(psf.Niy/2)-1,ROUND(psf.Nix/2)-1);
                    imTempP = &imTempSc;
                    break;
                }

#ifdef DEBUG

                _Im().resize(intExp);
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(intExp)
                dAij(_Im(),i,j) = dAij(intExp,i,j);
                _Im.write("psfxr-intExp.spi");
                _Im().resize(*imTempP);
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
                dAij(_Im(),i,j) = dAij(*imTempP,i,j);
                _Im.write("psfxr-imTempEsc.spi");
#endif

                psf.applyOTF(*imTempP, imOutGlobal.Zoff()- k*psf.dzo);

#ifdef DEBUG

                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
                dAij(_Im(),i,j) = dAij(*imTempP,i,j);
                _Im.write("psfxr-imTempEsc2.spi");
#endif

                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
                dAij(imOut(),i,j) += dAij(*imTempP,i,j);

                //        imOut.write("psfxr-imout.spi");
            }
        }
        priorLast = last;

        std::cerr << "th" << thread_id << ": Finished work from " << first << " to " << last <<std::endl;

    }

    //Lock to update the total summatory
    mutex.lock();
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
    dAij(imOutGlobal(),i,j) += dAij(imOut(),i,j);
    mutex.unlock();

    barrier->wait();

    if (thread_id==0)
    {
        imOutGlobal() = 1-imOutGlobal();

        switch (psf.AdjustType)
        {
        case PSFXR_INT:
            //    imOut().selfScaleToSize(psf.Noy, psf.Nox);
            selfScaleToSize(LINEAR,imOutGlobal(), psf.Nox, psf.Noy);
            break;
        case PSFXR_ZPAD:
            imOutGlobal().window(-ROUND(psf.Noy/2)+1,-ROUND(psf.Nox/2)+1,ROUND(psf.Noy/2)-1,ROUND(psf.Nox/2)-1);
            break;
        }

        imOutGlobal().setXmippOrigin();
    }
    std::cerr << "th" << thread_id << ": Thread Finished" <<std::endl;
}
