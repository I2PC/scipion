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

#include "projectTomography.h"
#include "directions.h"

#include <data/args.h>

void ProgProjectTomography::defineParams()
{
    addUsageLine("Generate projections as if it were a single axis tomographic series.");
    addParamsLine(" -i <Proj_param_file>   :MetaData file with projection parameters.");
    addParamsLine("                        :Check the manual for a description of the parameters.");
    addParamsLine(" alias --input;");
    addParamsLine("[-o <metadata_file=\"\">]    : Output Metadata file with all the generated projections.");
    addParamsLine("[-show_angles]          : Print angles value for each projection.");
    addParamsLine("[-only_create_angles]   : Projections are not calculated, only the angles values.");
}

void ProgProjectTomography::readParams()
{
    fn_proj_param = getParam("-i");
    fn_sel_file   = getParam("-o");

    only_create_angles = checkParam("-only_create_angles");
    tell=0;
    if (checkParam("-show_angles"))
        tell |= TELL_SHOW_ANGLES;
}

void ProgProjectTomography::run()
{

    Projection         proj;
    MetaData           SF;

    randomize_random_generator();

    // Read projection parameters and produce side information
    ParametersProjectionTomography proj_prm;
    PROJECT_Tomography_Side_Info side;
    proj_prm.read(fn_proj_param);
    proj_prm.tell = tell;
    side.produce_Side_Info(proj_prm);

    // Project
    int ProjNo = 0;
    if (!only_create_angles)
    {
        // Really project
        ProjNo = PROJECT_Tomography_Effectively_project(proj_prm, side,
                 proj, SF);
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

/* Produce Side Information ================================================ */
void PROJECT_Tomography_Side_Info::produce_Side_Info(
    const ParametersProjectionTomography &prm)
{
    phantomVol.read(prm.fnPhantom);
    phantomVol().setXmippOrigin();
}

/* Effectively project ===================================================== */
int PROJECT_Tomography_Effectively_project(
    const ParametersProjectionTomography &prm,
    PROJECT_Tomography_Side_Info &side, Projection &proj, MetaData &SF)
{
    int expectedNumProjs = FLOOR((prm.tiltF-prm.tilt0)/prm.tiltStep);
    int numProjs=0;
    SF.clear();
    std::cerr << "Projecting ...\n";
    if (!(prm.tell&TELL_SHOW_ANGLES))
        init_progress_bar(expectedNumProjs);
    MetaData DF_movements;
    DF_movements.setComment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
    double tRot,tTilt,tPsi,rot,tilt,psi;


    int idx=prm.starting;
    size_t id;
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

        // Really project ....................................................
        project_Volume_offCentered(side.phantomVol(), proj,
                                   prm.proj_Ydim, prm.proj_Xdim, prm.axisRot, prm.axisTilt,
                                   prm.raxis, angle, 0, inPlaneShift);

        // Add noise in angles and voxels ....................................
        proj.getEulerAngles(tRot, tTilt,tPsi);

        rot  = tRot  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        tilt = tTilt + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        psi  = tPsi  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);

        proj.setEulerAngles(rot,tilt,psi);

        id=DF_movements.addObject();
        DF_movements.setValue(MDL_ANGLEROT,tRot,id);
        DF_movements.setValue(MDL_ANGLETILT,tTilt,id);
        DF_movements.setValue(MDL_ANGLEPSI,tPsi,id);
        DF_movements.setValue(MDL_ANGLEROT2,rot,id);
        DF_movements.setValue(MDL_ANGLETILT2,tilt,id);
        DF_movements.setValue(MDL_ANGLEPSI2,psi,id);
        DF_movements.setValue(MDL_SHIFTX,shiftX,id);
        DF_movements.setValue(MDL_SHIFTY,shiftY,id);

        IMGMATRIX(proj).addNoise(prm.Npixel_avg, prm.Npixel_dev, "gaussian");

        // Save ..............................................................
        if (prm.tell&TELL_SHOW_ANGLES)
            std::cout << idx << "\t" << tRot << "\t"
            << tTilt << "\t" << tPsi << std::endl;
        else if ((expectedNumProjs % XMIPP_MAX(1, numProjs / 60))
                 == 0)
            progress_bar(numProjs);
        proj.write(fn_proj);
        numProjs++;
        idx++;
        size_t objId = SF.addObject();
        SF.setValue(MDL_IMAGE,fn_proj,objId);
        SF.setValue(MDL_ENABLED,1,objId);
    }
    if (!(prm.tell&TELL_SHOW_ANGLES))
        progress_bar(expectedNumProjs);

    DF_movements.write(prm.fnProjectionSeed + "_movements.txt");
    return numProjs;
}
