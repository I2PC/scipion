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


#include "mpi_project_XR.h"

/* Structure with the information about the different
 * projections to be distributed among nodes */
struct mpiProjData
{
    FileName fn_proj;
    double rot, tRot;
    double tilt, tTilt;
    double psi, tPsi;
    double xShift;
    double yShift;
    double zShift;
};

ProgMPIXrayProject::~ProgMPIXrayProject()
{
    delete node;
}

void ProgMPIXrayProject::defineParams()
{
    ProgXrayProject::defineParams();
    clearUsage();
    addUsageLine("MPI Generate projections as in a X-ray microscope from a 3D Xmipp volume.");

}
void ProgMPIXrayProject::read(int argc, char** argv)
{
    node = new MpiNode(argc, argv);
    if (!node->isMaster())
        verbose = 0;
    ProgXrayProject::read(argc, (const char **)argv);
}

void ProgMPIXrayProject::run()
{
    preRun();

    // Project
    projMD.setComment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");

    //    // Assign mpi node information
    //    MpiNode & node = *node;

    // Calculation of data to be distributed in nodes -----------------------------------------
    std::vector<mpiProjData> mpiData;
    mpiProjData data;

    size_t idx = FIRST_IMAGE;
    size_t id;

    for (double angle = projParam.tilt0; angle <= projParam.tiltF; angle += projParam.tiltStep)
    {
        if (projParam.singleProjection)
            data.fn_proj = projParam.fnOut;
        else
            data.fn_proj.compose(idx, projParam.fnRoot + ".stk");

        // Choose Center displacement ........................................
        double shiftX = rnd_gaus(projParam.Ncenter_avg, projParam.Ncenter_dev);
        double shiftY = rnd_gaus(projParam.Ncenter_avg, projParam.Ncenter_dev);
        Matrix1D<double> inPlaneShift(3);
        VECTOR_R3(inPlaneShift,shiftX,shiftY,0);

        projParam.calculateProjectionAngles(proj,angle, 0,inPlaneShift);
        proj.getEulerAngles(data.tRot, data.tTilt,data.tPsi);
        proj.getShifts(data.xShift, data.yShift, data.zShift);

        // Add noise in angles and voxels ....................................

        data.rot  = data.tRot  + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);
        data.tilt = data.tTilt + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);
        data.psi  = data.tPsi  + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);

        //MPI            proj.setEulerAngles(rot,tilt,psi);

        if (node->isMaster())
        {
            id = projMD.addObject();
            projMD.setValue(MDL_IMAGE,data.fn_proj, id);

            projMD.setValue(MDL_ANGLE_ROT,data.tRot, id);
            projMD.setValue(MDL_ANGLE_TILT,data.tTilt, id);
            projMD.setValue(MDL_ANGLE_PSI,data.tPsi, id);
            projMD.setValue(MDL_ANGLE_ROT2,data.rot, id);
            projMD.setValue(MDL_ANGLE_TILT2,data.tilt, id);
            projMD.setValue(MDL_ANGLE_PSI2,data.psi, id);
            projMD.setValue(MDL_SHIFT_X,shiftX, id);
            projMD.setValue(MDL_SHIFT_Y,shiftY, id);
        }

        mpiData.push_back(data);
        idx++;
    }

    // Creation of MPI Job Handler file

    MpiTaskDistributor *jobHandler;
    long long int nodeBlockSize = 1;
    jobHandler = new MpiTaskDistributor(mpiData.size(), nodeBlockSize, node);

    size_t first = 0, last = 0;

    if (node->isMaster())
    {
        // Save metadata file with angles and shift info
        if (!projParam.singleProjection)
            projMD.write(projParam.fnRoot + ".sel");

        if (!(projParam.show_angles))
            init_progress_bar(mpiData.size());

        // Creation of output file to reserve space
        createEmptyFile(mpiData[mpiData.size()-1].fn_proj,
                        XMIPP_MIN(XSIZE(MULTIDIM_ARRAY(phantom.iniVol)),(size_t)projParam.proj_Xdim),
                        XMIPP_MIN(YSIZE(MULTIDIM_ARRAY(phantom.iniVol)),(size_t)projParam.proj_Ydim));

        std::cout << "Projecting ...\n";
    }

    node->barrierWait();

    MPI_Barrier(MPI_COMM_WORLD);

    // Parallel node jobs
    while (jobHandler->getTasks(first, last))
    {
        for (size_t k = first; k <= last; k++)
        {
            std::cout << "Node: " << node->rank << " - Task: " << k <<std::endl;
            proj.setEulerAngles(mpiData[k].tRot, mpiData[k].tTilt,mpiData[k].tPsi);
            proj.setShifts(mpiData[k].xShift, mpiData[k].yShift, mpiData[k].zShift);

            //Reset thread task distributor
            td->clear();

            // Really project ....................................................
            if (!projParam.only_create_angles)
            {
                XrayRotateAndProjectVolumeOffCentered(phantom, psf, proj, stdProj, projParam.proj_Ydim, projParam.proj_Xdim);
                proj.write(mpiData[k].fn_proj);
            }

            // Add noise in angles and voxels ....................................
            //            proj.setEulerAngles(mpiData[k].rot,mpiData[k].tilt,mpiData[k].psi); // Geo info is not stored in header file
            IMGMATRIX(proj).addNoise(projParam.Npixel_avg, projParam.Npixel_dev, "gaussian");

            if (projParam.show_angles)
                std::cout << "Node: " << node->rank << "\t" << proj.rot() << "\t"
                << proj.tilt() << "\t" << proj.psi() << std::endl;

            if (node->isMaster() && !(projParam.show_angles))
                progress_bar(k+1);
        }
    }
    jobHandler->wait();

    delete jobHandler;
    postRun();

    return;
}

///* Effectively project ===================================================== */
//int PROJECT_mpi_XR_Effectively_project(
//    ParametersXrayProjectMPI &prm,
//    XrayProjPhantom &side,
//    Projection &proj,
//    XRayPSF &psf,
//    MetaData &SF)
//{
//    // Threads stuff
//
//    XrayThread *dataThread = new XrayThread;
//
//    dataThread->psf= &psf;
//    dataThread->vol = &side.rotVol;
//    dataThread->imOut = &proj;
//
//    size_t threadBlockSize, numberOfJobs= side.iniVol().zdim;
//    numberOfThreads = psf.nThr;
//
//    threadBlockSize = (numberOfThreads == 1) ? numberOfJobs : numberOfJobs/numberOfThreads/2;
//
//    //Create the job handler to distribute thread jobs
//    td = new ThreadTaskDistributor(numberOfJobs, threadBlockSize);
//    barrier = new Barrier(numberOfThreads-1);
//
//    //Create threads to start working
//    thMgr = new ThreadManager(numberOfThreads,(void*) dataThread);
//
//    int numProjs=0;
//    SF.clear();
//    MpiNode & node = *prm.node;
//    MetaData DF_movements;
//    DF_movements.setComment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
//    double tRot,tTilt,tPsi,rot,tilt,psi;
//
//    // Calculation of data to be distributed in nodes
//    std::vector<mpiProjData> mpiData;
//    mpiProjData data;
//
//    int idx=prm.starting;
//    size_t id;
//    for (double angle=prm.tilt0; angle<=prm.tiltF; angle+=prm.tiltStep)
//    {
//        data.fn_proj.compose(prm.fnOut, idx,
//                             prm.fn_projection_extension);
//
//        // Choose Center displacement ........................................
//        double shiftX     = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
//        double shiftY    = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
//        Matrix1D<double> inPlaneShift(3);
//        VECTOR_R3(inPlaneShift,shiftX,shiftY,0);
//
//        prm.calculateProjectionAngles(proj,angle, 0,inPlaneShift);
//        proj.getEulerAngles(data.tRot, data.tTilt,data.tPsi);
//        proj.getShifts(data.xShift, data.yShift, data.zShift);
//
//        // Add noise in angles and voxels ....................................
//
//        data.rot  = data.tRot  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
//        data.tilt = data.tTilt + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
//        data.psi  = data.tPsi  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
//
//        //MPI            proj.setEulerAngles(rot,tilt,psi);
//
//        if (node.isMaster())
//        {
//            id = DF_movements.addObject();
//            DF_movements.setValue(MDL_ANGLE_ROT,data.tRot, id);
//            DF_movements.setValue(MDL_ANGLE_TILT,data.tTilt, id);
//            DF_movements.setValue(MDL_ANGLE_PSI,data.tPsi, id);
//            DF_movements.setValue(MDL_ANGLE_ROT2,data.rot, id);
//            DF_movements.setValue(MDL_ANGLE_TILT2,data.tilt, id);
//            DF_movements.setValue(MDL_ANGLE_PSI2,data.psi, id);
//            DF_movements.setValue(MDL_SHIFT_X,shiftX, id);
//            DF_movements.setValue(MDL_SHIFT_Y,shiftY, id);
//
//            id = SF.addObject();
//            SF.setValue(MDL_IMAGE,data.fn_proj, id);
//            SF.setValue(MDL_ENABLED,1, id);
//        }
//        mpiData.push_back(data);
//        idx++;
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    // Creation of Job Handler file
//
//    FileTaskDistributor *jobHandler;
//    long long int nodeBlockSize = 1;
//    jobHandler = new FileTaskDistributor(mpiData.size(), nodeBlockSize, &node);
//
//
//    if (node.isMaster())
//    {
//        DF_movements.write(prm.fnOut + "_movements.txt");
//        std::cerr << "Projecting ...\n";
//    }
//
//
//    long long int first = -1, last = -1;
//
//    if (!(prm.show_angles))
//        init_progress_bar(mpiData.size());
//
//    // Parallel node jobs
//    while (jobHandler->getTasks(first, last))
//    {
//        for (long long int k = first; k <= last; k++)
//        {
//            std::cout << "Node: " << node.rank << " - Task: " << k <<std::endl;
//            proj.setEulerAngles(mpiData[k].tRot, mpiData[k].tTilt,mpiData[k].tPsi);
//            proj.setShifts(mpiData[k].xShift, mpiData[k].yShift, mpiData[k].zShift);
//
//            //Reset thread task distributor
//            td->clear();
//            // Really project ....................................................
//            XrayProjectVolumeOffCentered(side, psf, proj,prm.proj_Ydim, prm.proj_Xdim);
//
//
//            // Add noise in angles and voxels ....................................
//            proj.setEulerAngles(mpiData[k].rot,mpiData[k].tilt,mpiData[k].psi);
//            IMGMATRIX(proj).addNoise(prm.Npixel_avg, prm.Npixel_dev, "gaussian");
//
//            // Save ..............................................................
//            if (prm.show_angles)
//                std::cout << "Node: " << node.rank << "\t" << proj.rot() << "\t"
//                << proj.tilt() << "\t" << proj.psi() << std::endl;
//
//            proj.write(mpiData[k].fn_proj);
//
//            progress_bar(k+1);
//            numProjs++;
//        }
//    }
//    delete jobHandler;
//    //Terminate threads and free memory
//    delete td;
//    delete thMgr;
//    delete barrier;
//
//    return numProjs;
//}
