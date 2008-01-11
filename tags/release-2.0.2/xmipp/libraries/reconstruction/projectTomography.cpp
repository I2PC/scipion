/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "projectTomography.h"
#include "directions.h"

#include <data/args.h>

/* Read from command line ================================================== */
void Prog_Project_Tomography_Parameters::read(int argc, char **argv)
{
    fn_proj_param = getParameter(argc, argv, "-i");
    fn_sel_file   = getParameter(argc, argv, "-o", "");
    only_create_angles = checkParameter(argc, argv, "-only_create_angles");
    tell=0;
    if (checkParameter(argc, argv, "-show_angles"))
        tell |= TELL_SHOW_ANGLES;
}

/* Usage =================================================================== */
void Prog_Project_Tomography_Parameters::usage()
{
    printf("\nUsage:\n\n");
    printf("project -i <Parameters File> \n"
           "       [-o <sel_file>]\n"
           "       [-show_angles]\n"
           "       [-only_create_angles]\n");
    printf(
        "\tWhere:\n"
        "\t<Parameters File>:  File containing projection parameters\n"
        "\t                    check the manual for a description of the parameters\n"
        "\t<sel_file>:         This is a selection file with all the generated\n"
        "\t                    projections\n");
}

/* Read Projection Parameters ============================================== */
void Projection_Tomography_Parameters::read(const FileName &fn_proj_param)
{
    FILE    *fh_param;
    char    line[201];
    int     lineNo = 0;
    char    *auxstr;

    if ((fh_param = fopen(fn_proj_param.c_str(), "r")) == NULL)
        REPORT_ERROR(3005,
                     (string)"Projection_Tomography_Parameters::read: There is a problem "
                     "opening the file " + fn_proj_param);
    while (fgets(line, 200, fh_param) != NULL)
    {
        if (line[0] == 0)    continue;
        if (line[0] == '#')  continue;
        if (line[0] == '\n') continue;
        switch (lineNo)
        {
        case 0:
            fnPhantom = firstWord(line, 3007,
                                    "Projection_Tomography_Parameters::read: Phantom name not found");
            lineNo = 1;
            break;
        case 1:
            fnProjectionSeed =
                firstWord(line, 3007,
                           "Projection_Tomography_Parameters::read: Error in Projection seed");
            // Next two parameters are optional
            auxstr = nextToken();
            if (auxstr != NULL) starting =
                    textToInteger(auxstr, 3007,
                         "Projection_Tomography_Parameters::read: Error in First "
                         "projection number");
            fn_projection_extension = nextToken();
            lineNo = 2;
            break;
        case 2:
            proj_Xdim = textToInteger(firstToken(line), 3007,
                             "Projection_Tomography_Parameters::read: Error in X dimension");
            proj_Ydim = textToInteger(nextToken(), 3007,
                             "Projection_Tomography_Parameters::read: Error in Y dimension");
            lineNo = 3;
            break;
        case 3:
            axisRot = textToDouble(firstToken(line), 3007,
                             "Projection_Tomography_Parameters::read: Error in axisRot");
            axisTilt = textToDouble(nextToken(), 3007,
                             "Projection_Tomography_Parameters::read: Error in axisTilt");
            lineNo = 4;
            break;
        case 4:
	    raxis.resize(3);
	    XX(raxis) = textToDouble(firstToken(line), 3007,
                             "Projection_Tomography_Parameters::read: Error in X component of raxis");
	    YY(raxis) = textToDouble(nextToken(), 3007,
                             "Projection_Tomography_Parameters::read: Error in Y component of raxis");
	    ZZ(raxis) = textToDouble(nextToken(), 3007,
                             "Projection_Tomography_Parameters::read: Error in Z component of raxis");
	    lineNo = 5;
            break;
        case 5:
	    tilt0 = textToDouble(firstToken(line), 3007,
                             "Projection_Tomography_Parameters::read: Error in tilt0");
	    tiltF = textToDouble(nextToken(), 3007,
                             "Projection_Tomography_Parameters::read: Error in tiltF");
	    tiltStep = textToDouble(nextToken(), 3007,
                             "Projection_Tomography_Parameters::read: Error in tiltStep");
	    lineNo = 6;
	    break;
        case 6:
            Nangle_dev = textToFloat(firstWord(line), 3007,
                                  "Projection_Tomography_Parameters::read: Error in angular noise");
            auxstr = nextToken();
            if (auxstr != NULL)
                Nangle_avg = textToFloat(auxstr, 3007,
                                      "Projection_Tomography_Parameters::read: Error in angular bias");
            else Nangle_avg = 0;
            lineNo = 7;
            break;
        case 7:
            Npixel_dev = textToFloat(firstWord(line), 3007,
                              "Projection_Tomography_Parameters::read: Error in pixel noise");
            auxstr = nextToken();
            if (auxstr != NULL)
                Npixel_avg = textToFloat(auxstr, 3007,
                                  "Projection_Tomography_Parameters::read: Error in pixel bias");
            else Npixel_avg = 0;
            lineNo = 8;
            break;
        case 8:
            Ncenter_dev = textToFloat(firstWord(line), 3007,
                               "Projection_Tomography_Parameters::read: Error in center noise");
            auxstr = nextToken();
            if (auxstr != NULL)
                Ncenter_avg = textToFloat(auxstr, 3007,
                                   "Projection_Tomography_Parameters::read: Error in center bias");
            else Ncenter_avg = 0;
            lineNo = 9;
            break;
        } /* switch end */
    } /* while end */
    if (lineNo != 9)
        REPORT_ERROR(3007, (string)"Projection_Tomography_Parameters::read: I "
                     "couldn't read all parameters from file " + fn_proj_param);
    fclose(fh_param);
}

/* Produce Side Information ================================================ */
void PROJECT_Tomography_Side_Info::produce_Side_Info(
    const Projection_Tomography_Parameters &prm,
    const Prog_Project_Tomography_Parameters &prog_prm)
{
    phantomVol.read(prm.fnPhantom);
    phantomVol().setXmippOrigin();
}

/* Effectively project ===================================================== */
int PROJECT_Tomography_Effectively_project(
    const Projection_Tomography_Parameters &prm,
    PROJECT_Tomography_Side_Info &side, Projection &proj, SelFile &SF)
{
    int expectedNumProjs = FLOOR((prm.tiltF-prm.tilt0)/prm.tiltStep);
    int numProjs=0;
    SF.clear();
    std::cerr << "Projecting ...\n";
    if (!(prm.tell&TELL_SHOW_ANGLES)) init_progress_bar(expectedNumProjs);
    SF.reserve(expectedNumProjs);
    DocFile DF_movements;
    DF_movements.append_comment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
    Matrix1D<double> movements(8);

    int idx=prm.starting;
    for (double angle=prm.tilt0; angle<=prm.tiltF; angle+=prm.tiltStep)
    {
        FileName fn_proj;              // Projection name
        fn_proj.compose(prm.fnProjectionSeed, idx,
                        prm.fn_projection_extension);

        // Choose Center displacement ........................................
        double shiftX     = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        double shiftY 	  = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        movements(6) = shiftX;
        movements(7) = shiftY;
	Matrix1D<double> inPlaneShift(3);
	VECTOR_R3(inPlaneShift,shiftX,shiftY,0);

        // Really project ....................................................
        project_Volume_offCentered(side.phantomVol(), proj,
	    prm.proj_Ydim, prm.proj_Xdim, prm.axisRot, prm.axisTilt,
	    prm.raxis, angle, 0, inPlaneShift);

        // Add noise in angles and voxels ....................................
	movements(0)=proj.rot();
	movements(1)=proj.tilt();
	movements(2)=proj.psi();
        movements(3)=rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        movements(4)=rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        movements(5)=rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        proj.set_eulerAngles(movements(0)+movements(3),
	    movements(1)+movements(4),movements(2)+movements(5));
        IMGMATRIX(proj).add_noise(prm.Npixel_avg, prm.Npixel_dev, "gaussian");

        // Save ..............................................................
        if (prm.tell&TELL_SHOW_ANGLES)
	    std::cout << idx << "\t" << movements(0) << "\t"
	    	      << movements(1) << "\t" << movements(2) << std::endl;
        else if ((expectedNumProjs % XMIPP_MAX(1, numProjs / 60)) == 0)
            progress_bar(numProjs);
        proj.write(fn_proj);
        numProjs++;
	idx++;
        SF.insert(fn_proj, SelLine::ACTIVE);
        DF_movements.append_data_line(movements);
    }
    if (!(prm.tell&TELL_SHOW_ANGLES)) progress_bar(expectedNumProjs);

    DF_movements.write(prm.fnProjectionSeed + "_movements.txt");
    return numProjs;
}

/* ROUT_project ============================================================ */
int ROUT_Tomography_project(Prog_Project_Tomography_Parameters &prm,
    Projection &proj, SelFile &SF)
{
    randomize_random_generator();

    // Read projection parameters and produce side information
    Projection_Tomography_Parameters proj_prm;
    PROJECT_Tomography_Side_Info side;
    proj_prm.read(prm.fn_proj_param);
    proj_prm.tell=prm.tell;
    side.produce_Side_Info(proj_prm, prm);

    // Project
    int ProjNo = 0;
    if (!prm.only_create_angles)
    {
        // Really project
        ProjNo = PROJECT_Tomography_Effectively_project(proj_prm, side,
    	    proj, SF);
        // Save SelFile
        if (prm.fn_sel_file != "") SF.write(prm.fn_sel_file);
    }
    else
    {
        cout << side.DF;
    }
    return ProjNo;
}
