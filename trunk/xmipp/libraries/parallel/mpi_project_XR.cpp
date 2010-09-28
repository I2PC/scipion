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

MPIProgProjectXR::~MPIProgProjectXR()
{
    delete node;
}

void MPIProgProjectXR::defineParams()
{
    ProgProjectXR::defineParams();
    clearUsage();
    addUsageLine("MPI Generate projections as in a X-ray microscope from a 3D Xmipp volume.");

}
void MPIProgProjectXR::read(int argc, char** argv)
{
    node = new MpiNode(argc, argv);
    ProgProjectXR::read(argc, argv);
}

void MPIProgProjectXR::run()
{
    randomize_random_generator();

    Projection  proj;
    MetaData    SF;

    // Read Microscope optics parameters and produce side information
    XRayPSF psf;
    psf.verbose = verbose;
    psf.nThr = nThr;
    psf.read(fn_psf_xr);
    psf.produceSideInfo();

    // Read projection parameters and produce side information
    Projection_mpi_XR_Parameters mpi_proj_prm;
    mpi_proj_prm.node = node;
    mpi_proj_prm.read(fn_proj_param);
    mpi_proj_prm.tell=tell;
    PROJECT_XR_Side_Info side;

    side.produce_Side_Info(mpi_proj_prm);
    //    psf.adjustParam(side.phantomVol);

    // Project
    int ProjNo = 0;
    if (!only_create_angles)
    {
        // Really project
        ProjNo = PROJECT_mpi_XR_Effectively_project(mpi_proj_prm, side,
                 proj, psf, SF);
        // Save SelFile
        if (node->isMaster() && fn_sel_file != "")
            SF.write(fn_sel_file);
    }
    else
    {
        psf.adjustParam(side.phantomVol);
        if (node->isMaster())
            side.DF.write("/dev/stdout");
    }
    return;
}

/* Read parameters --------------------------------------------------------- */
void Projection_mpi_XR_Parameters::read(const FileName &fn_proj_param)
{
    Projection_XR_Parameters::read(fn_proj_param);
}


//Some global variables
Mutex mutex;
Barrier * barrier;
ThreadManager * thMgr;
ParallelTaskDistributor * td;
//FileTaskDistributor     * td;
int numberOfThreads;


/* Effectively project ===================================================== */
int PROJECT_mpi_XR_Effectively_project(
    Projection_mpi_XR_Parameters &prm,
    PROJECT_XR_Side_Info &side,
    Projection &proj,
    XRayPSF &psf,
    MetaData &SF)
{

    XrayThread *dataThread = new XrayThread;

    dataThread->psf= &psf;
    dataThread->vol = &side.rotPhantomVol;
    dataThread->imOut = &proj;

    longint threadBlockSize, numberOfJobs= side.phantomVol().zdim;
    numberOfThreads = psf.nThr;

    threadBlockSize = (numberOfThreads == 1) ? numberOfJobs : numberOfJobs/numberOfThreads/2;

    //Create the job handler to distribute jobs
    td = new ThreadTaskDistributor(numberOfJobs, threadBlockSize);
    //    td = new FileTaskDistributor(numberOfJobs,threadBlockSize);
    barrier = new Barrier(numberOfThreads-1);

    //Create threads to start working

    thMgr = new ThreadManager(numberOfThreads,(void*) dataThread);





    int numProjs=0;
    SF.clear();
    MpiNode & node = *prm.node;
    MetaData DF_movements;
    DF_movements.setComment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
    double tRot,tTilt,tPsi,rot,tilt,psi;

    // Calculation of data to be distributed in nodes
    std::vector<mpiProjData> mpiData;
    mpiProjData data;

    int idx=prm.starting;
    for (double angle=prm.tilt0; angle<=prm.tiltF; angle+=prm.tiltStep)
    {
        data.fn_proj.compose(prm.fnProjectionSeed, idx,
                             prm.fn_projection_extension);

        // Choose Center displacement ........................................
        double shiftX     = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        double shiftY    = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        Matrix1D<double> inPlaneShift(3);
        VECTOR_R3(inPlaneShift,shiftX,shiftY,0);

        prm.calculateProjectionAngles(proj,angle, 0,inPlaneShift);
        proj.getEulerAngles(data.tRot, data.tTilt,data.tPsi);
        proj.getShifts(data.xShift, data.yShift, data.zShift);

        // Add noise in angles and voxels ....................................

        data.rot  = data.tRot  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        data.tilt = data.tTilt + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        data.psi  = data.tPsi  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);

        //MPI            proj.setEulerAngles(rot,tilt,psi);

        if (node.isMaster())
        {
            DF_movements.addObject();
            DF_movements.setValue(MDL_ANGLEROT,data.tRot);
            DF_movements.setValue(MDL_ANGLETILT,data.tTilt);
            DF_movements.setValue(MDL_ANGLEPSI,data.tPsi);
            DF_movements.setValue(MDL_ANGLEROT2,data.rot);
            DF_movements.setValue(MDL_ANGLETILT2,data.tilt);
            DF_movements.setValue(MDL_ANGLEPSI2,data.psi);
            DF_movements.setValue(MDL_SHIFTX,shiftX);
            DF_movements.setValue(MDL_SHIFTY,shiftY);

            SF.addObject();
            SF.setValue(MDL_IMAGE,data.fn_proj);
            SF.setValue(MDL_ENABLED,1);
        }
        //MPI            IMGMATRIX(proj).addNoise(prm.Npixel_avg, prm.Npixel_dev, "gaussian");
        mpiData.push_back(data);
        idx++;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Creation of Job Handler file

    FileTaskDistributor *jobHandler;
    char fName[L_tmpnam];
    long long int nodeBlockSize = 1;


    if (node.isMaster())
    {
        DF_movements.write(prm.fnProjectionSeed + "_movements.txt");
        std::cerr << "Projecting ...\n";
    }

    jobHandler = new FileTaskDistributor(mpiData.size(), nodeBlockSize, &node);

    long long int first = -1, last = -1;

    if (!(prm.tell&TELL_SHOW_ANGLES))
        init_progress_bar(mpiData.size());

    // Parallel jobs
    while (jobHandler->getTasks(first, last))
    {
        for (long long int k = first; k <= last; k++)
        {
            std::cout << "Node: " << node.rank << " - Task: " << k <<std::endl;
            proj.setEulerAngles(mpiData[k].tRot, mpiData[k].tTilt,mpiData[k].tPsi);
            proj.setShifts(mpiData[k].xShift, mpiData[k].yShift, mpiData[k].zShift);

            //Reset thread task distributor
            td->clear();
            // Really project ....................................................
            project_xr_Volume_offCentered(side, psf, proj,prm.proj_Ydim, prm.proj_Xdim);


            // Add noise in angles and voxels ....................................
            proj.setEulerAngles(mpiData[k].rot,mpiData[k].tilt,mpiData[k].psi);
            IMGMATRIX(proj).addNoise(prm.Npixel_avg, prm.Npixel_dev, "gaussian");

            // Save ..............................................................
            if (prm.tell&TELL_SHOW_ANGLES)
                std::cout << "Node: " << node.rank << "\t" << proj.rot() << "\t"
                << proj.tilt() << "\t" << proj.psi() << std::endl;

            proj.write(mpiData[k].fn_proj);

            progress_bar(k+1);
            numProjs++;
        }
    }
    delete jobHandler;
    //Terminate threads and free memory
    delete td;
    delete thMgr;
    delete barrier;

    return numProjs;
}
