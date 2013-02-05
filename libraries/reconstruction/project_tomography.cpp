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

#include "project_tomography.h"
#include "directions.h"

#include <data/args.h>

void ProgProjectTomography::defineParams()
{
    addUsageLine("Generate projections as if it were a single axis tomographic series.");
    addSeeAlsoLine("phantom_create, phantom_project, xray_project");
    //Params
    projParam.defineParams(this);
    //Examples
    addExampleLine("Generating a set of projections using a projection parameter:", false);
    addExampleLine("xmipp_xray_project -i volume.vol --oroot images --params projParams.xmd");
    addExampleLine("Generating a single projection at 45 degrees around X axis:", false);
    addExampleLine("xmipp_xray_project -i volume.vol -o image.spi --angles 45 0 90");
    addExampleLine("Generating a single projection at 45 degrees around Y axis:", false);
    addExampleLine("xmipp_xray_project -i volume.vol -o image.spi --angles 45 90 90");
    //Example projection file
    addExampleLine("In the following link you can find an example of projection parameter file:",false);
    addExampleLine(" ",false);
    addExampleLine("http://sourceforge.net/p/testxmipp/code/ci/3.0/tree/input/tomoProjection.param",false);
}

void ProgProjectTomography::readParams()
{
    projParam.readParams(this);
}

void ProgProjectTomography::run()
{

    Projection         proj;
    MetaData           projMD;

    randomize_random_generator();

    TomoProjectSideInfo side;
    side.produceSideInfo(projParam);

    // Project
    int expectedNumProjs = FLOOR((projParam.tiltF-projParam.tilt0)/projParam.tiltStep);
    int numProjs=0;

    std::cerr << "Projecting ...\n";
    if (!(projParam.show_angles))
        init_progress_bar(expectedNumProjs);

    projMD.setComment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
    double tRot,tTilt,tPsi,rot,tilt,psi;
    FileName fn_proj;              // Projection name
    int idx = 1;
    size_t objId;

    for (double angle=projParam.tilt0; angle<=projParam.tiltF; angle+=projParam.tiltStep)
    {
        if (projParam.singleProjection)
            fn_proj = projParam.fnOut;
        else
            fn_proj.compose(idx, projParam.fnRoot + ".stk");

        // Choose Center displacement ........................................
        double shiftX     = rnd_gaus(projParam.Ncenter_avg, projParam.Ncenter_dev);
        double shiftY    = rnd_gaus(projParam.Ncenter_avg, projParam.Ncenter_dev);
        Matrix1D<double> inPlaneShift(3);
        VECTOR_R3(inPlaneShift,shiftX,shiftY,0);

        projParam.calculateProjectionAngles(proj,angle, 0,inPlaneShift);

        // Really project ....................................................
        if (!projParam.only_create_angles)
            projectVolumeOffCentered(side.phantomVol(), proj,
                                     projParam.proj_Ydim, projParam.proj_Xdim);

        // Add noise in angles and voxels ....................................
        proj.getEulerAngles(tRot, tTilt,tPsi);

        rot  = tRot  + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);
        tilt = tTilt + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);
        psi  = tPsi  + rnd_gaus(projParam.Nangle_avg,  projParam.Nangle_dev);

        proj.setEulerAngles(rot,tilt,psi);

        objId = projMD.addObject();
        if (!projParam.only_create_angles)
        {
            proj.write(fn_proj, ALL_IMAGES, !projParam.singleProjection, WRITE_REPLACE);
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
            std::cout << idx << "\t" << tRot << "\t"
            << tTilt << "\t" << tPsi << std::endl;
        else if ((expectedNumProjs % XMIPP_MAX(1, numProjs / 60))
                 == 0)
            progress_bar(numProjs);

        numProjs++;
        idx++;
    }
    if (!projParam.show_angles)
        progress_bar(expectedNumProjs);

    if (!projParam.singleProjection)
    {
        projMD.setComment("Angles rot,tilt and psi contain noisy projection angles and rot2,tilt2 and psi2 contain actual projection angles");
        projMD.write(projParam.fnRoot + ".sel");
    }

    return;
}

/* Produce Side Information ================================================ */
void TomoProjectSideInfo::produceSideInfo(
    const ParametersProjectionTomography &prm)
{
    phantomVol.read(prm.fnPhantom);
    phantomVol().setXmippOrigin();
}
