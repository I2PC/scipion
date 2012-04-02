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

#include "project.h"
#include "directions.h"
#include "project_real_shears.h"

#include <data/args.h>

/* Read from command line ================================================== */
void ProgProject::readParams()
{
	fnPhantom = getParam("-i");
    fnOut = getParam("-o");
    samplingRate  = getDoubleParam("--sampling_rate");
    singleProjection = false;
    shears  = checkParam("--shears");

    bool doParams = checkParam("--params");
    bool doAngles = checkParam("--angles");

    if (doParams && doAngles)
        REPORT_ERROR(ERR_ARG_BADCMDLINE, "--params and --angles are mutually exclusive");
    else if (!doParams && !doAngles)
        REPORT_ERROR(ERR_ARG_BADCMDLINE, "You should provide --params or --angles");

    if (doParams)
    {
        fn_proj_param = getParam("--params");
        if (checkParam("--sym"))
            fn_sym        = getParam("--sym");
        only_create_angles = checkParam("--only_create_angles");
    }
    else //doAngles = true
    {
        singleProjection = true;
        only_create_angles = false;
        projSize = getIntParam("--xdim");
        rotSingle  = getDoubleParam("--angles",0);
        tiltSingle = getDoubleParam("--angles",1);
        psiSingle  = getDoubleParam("--angles",2);
    }
}

/* Usage =================================================================== */
void ProgProject::defineParams()
{
    addUsageLine("Generates projections from a volume");
    addSeeAlsoLine("tomo_project, xray_project, phantom_create");
    addExampleLine("Generating a set of projections",false);
    addExampleLine("xmipp_phantom_project -i volume.vol -o images.stk --params projectionParams.xmd");
    addExampleLine("Generating a top view from Z",false);
    addExampleLine("xmipp_phantom_project -i volume.vol -o image.xmp --angles 0 0 0");
    addExampleLine("Generating a side view from Y",false);
    addExampleLine("xmipp_phantom_project -i volume.vol -o image.xmp --angles 90 90 0");
    addExampleLine("Generating a side view from X",false);
    addExampleLine("xmipp_phantom_project -i volume.vol -o image.xmp --angles 0 90 0");
    addExampleLine("+In the following links you can find some examples of projection parameter files",false);
    addExampleLine("+ ",false);
    addExampleLine("+http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/phantomProject.param",false);
    addExampleLine("+ ",false);
    addExampleLine("+http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/uniformProjection_xmd.param",false);
    addExampleLine("+ ",false);
    addExampleLine("+http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/clusterProjection_xmd.param",false);
    addExampleLine("Creating a 2D crystal",false);
    addExampleLine("In order to create a 2D crystal, you can pass --params as a projection file with a second block for crystal projection.: ",false);
    addExampleLine(" xmipp_phantom_project   -i cylinder_with_axis.descr --oroot MRCproj --params MRCCrystalProj_xmd.param");
    addExampleLine("+In the following links you can find some examples of projection parameter files",false);
    addExampleLine("+ ",false);
    addExampleLine("+http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/Crystal/MRCCrystalProj_xmd.param",false);
    addExampleLine("+ ",false);
    addExampleLine("+http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/Crystal/cylinder_with_axis.descr",false);
    addExampleLine("+ ",false);
    addExampleLine("+http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/Crystal/MRC_crystal_projection_xmd.param",false);
    addParamsLine("   -i <volume_file>                           : Voxel volume, PDB or description file");
    addParamsLine("   -o <image_file>                            : Output stack or image");
    addParamsLine("  [--sampling_rate <Ts=1>]                    : It is only used for PDB phantoms");
    addParamsLine("  [--shears]                                  : Use real-shears algorithm");
    addParamsLine("                                              :+This algorithm is slower but more accurate. For a full description see");
    addParamsLine("                                              :+this Ph.D. [[http://biblion.epfl.ch/EPFL/theses/2003/2901/EPFL_TH2901.pdf][thesis]] (Chapter 3).");
    addParamsLine("== Generating a set of projections == ");
    addParamsLine("  [--params <parameters_file>]           : File containing projection parameters");
    addParamsLine("                                         : Check the manual for a description of the parameters");
    addParamsLine("  [--sym <sym_file>]                     : It is used for computing the assymetric unit");
    addParamsLine("  [--only_create_angles]                 : Do not create projections");
    addParamsLine("== Generating a single projection == ");
    addParamsLine("  [--angles <rot> <tilt> <psi>]          : Angles for a single projection");
    addParamsLine("  [--xdim <size=-1>]                     : Size of the projection");
    addParamsLine("                                         : For geometric descriptions and voxel volumes");
    addParamsLine("                                         : this parameter is not necessary");
}

/* Run ===================================================================== */
void ProgProject::run()
{
    Projection proj;
    MetaData SF;
    ROUT_project(*this, proj, SF);
}

/* Projection parameters from program parameters =========================== */
void ParametersProjection::from_prog_params(const ProgProject &prog_prm)
{
    read(prog_prm.fn_proj_param);
}

/* Read Projection Parameters ============================================== */
int translate_randomness(char * str)
{
    if (str == NULL)
        return ANGLE_RANGE_DETERMINISTIC;
    if (strcmp(str, "random_group") == 0)
        return ANGLE_RANGE_RANDOM_GROUPS;
    if (strcmp(str, "random") == 0)
        return ANGLE_RANGE_RANDOM;
    if (strcmp(str, "even") == 0)
        return ANGLE_EVENLY;
    REPORT_ERROR(ERR_PARAM_INCORRECT,
                 (String)"Prog_Project_Parameters::read: Not recognized randomness: "
                 + str);
}

void ParametersProjection::read(const FileName &fn_proj_param)
{
	if (fn_proj_param.isMetaData())
	{
		MetaData MD;
		MD.read((std::string)"block1@"+fn_proj_param.c_str());
		if (MD.isEmpty())
			REPORT_ERROR(ERR_IO_NOTOPEN,
						 (String)"Prog_Project_Parameters::read: There is a problem "
						 "opening the file " + fn_proj_param);

		std::vector <double> ParamVec;
		std::string RandStr;
		char * RandChar;
		size_t objId;

		// Read data from the MetaData
		objId = MD.firstObject();
		MD.getValue(MDL_DIMENSIONS_2D, ParamVec, objId);
		proj_Xdim = (int)ParamVec[0];proj_Ydim = (int)ParamVec[1];
		if (!MD.getValue(MDL_PRJ_ANGFILE, fn_angle, objId))
			fn_angle = "NULL";
		else
		if (!fn_angle.exists())
			REPORT_ERROR(ERR_IO_NOTEXIST, (String)"Prog_Project_Parameters::read: "
			"file " + fn_angle + " doesn't exist");
		if (MD.getValue(MDL_PRJ_ROT_RANGE,ParamVec, objId))
		{
			enable_angle_range = true;
			rot_range.ang0 = ParamVec[0];
			if (ParamVec.size() == 1)
			{
				rot_range.angF = rot_range.ang0;
				rot_range.samples = 1;
			}
			else
			{
				rot_range.angF = ParamVec[1];
				rot_range.samples = ParamVec[2];
				if (rot_range.ang0 == rot_range.angF)
					rot_range.samples = 1;
			}
			if (!MD.getValue(MDL_PRJ_ROT_RANDSTR,RandStr,objId))
				rot_range.randomness = ANGLE_RANGE_DETERMINISTIC;
			else
			{
				RandChar = new char[RandStr.length() + 1];
				strcpy(RandChar, RandStr.c_str());
				rot_range.randomness = translate_randomness(RandChar);
			}
			MD.getValue(MDL_PRJ_TILT_RANGE,ParamVec, objId);
			tilt_range.ang0 = ParamVec[0];
			if (ParamVec.size() == 1)
			{
				tilt_range.angF = tilt_range.ang0;
				tilt_range.samples = 1;
			}
			else
			{
				tilt_range.angF = ParamVec[1];
				tilt_range.samples = ParamVec[2];
				if (tilt_range.ang0 == tilt_range.angF)
					tilt_range.samples = 1;
			}
			if (!MD.getValue(MDL_PRJ_TILT_RANDSTR,RandStr,objId))
				tilt_range.randomness = ANGLE_RANGE_DETERMINISTIC;
			else
			{
				RandChar = new char[RandStr.length() + 1];
				strcpy(RandChar, RandStr.c_str());
				tilt_range.randomness = translate_randomness(RandChar);
			}
			MD.getValue(MDL_PRJ_PSI_RANGE,ParamVec, objId);
			psi_range.ang0 = ParamVec[0];
			if (ParamVec.size() == 1)
			{
				psi_range.angF = tilt_range.ang0;
				psi_range.samples = 1;
			}
			else
			{
				psi_range.angF = ParamVec[1];
				psi_range.samples = ParamVec[2];
				if (psi_range.ang0 == psi_range.angF)
					psi_range.samples = 1;
			}
			if (!MD.getValue(MDL_PRJ_PSI_RANDSTR,RandStr,objId))
				psi_range.randomness = ANGLE_RANGE_DETERMINISTIC;
			else
			{
				RandChar = new char[RandStr.length() + 1];
				strcpy(RandChar, RandStr.c_str());
			    psi_range.randomness = translate_randomness(RandChar);
			}
		}
		else
			enable_angle_range = false;
		MD.getValue(MDL_PRJ_ROT_Noise,ParamVec, objId);
		rot_range.Ndev = ParamVec[0];
		if (ParamVec.size()<2)
			rot_range.Navg = 0;
		else
			rot_range.Navg = ParamVec[1];
		MD.getValue(MDL_PRJ_TILT_Noise,ParamVec, objId);
		tilt_range.Ndev = ParamVec[0];
		if (ParamVec.size()<2)
			tilt_range.Navg = 0;
		else
			tilt_range.Navg = ParamVec[1];
		MD.getValue(MDL_PRJ_PSI_Noise,ParamVec, objId);
		psi_range.Ndev = ParamVec[0];
		if (ParamVec.size()<2)
			psi_range.Navg = 0;
		else
			psi_range.Navg = ParamVec[1];
		MD.getValue(MDL_NOISE_PIXEL_LEVEL,ParamVec, objId);
		Npixel_dev = ParamVec[0];
		if (ParamVec.size() < 2)
			Npixel_avg = 0;
		else
			Npixel_avg = ParamVec[1];
		MD.getValue(MDL_NOISE_COORD,ParamVec, objId);
		Ncenter_dev = ParamVec[0];
		if (ParamVec.size() < 2)
			Ncenter_avg = 0;
		else
			Ncenter_avg = ParamVec[1];
	}
	else
	{
		FILE    *fh_param;
		char    line[201];
		int     lineNo = 0;
		char    *auxstr;

		if ((fh_param = fopen(fn_proj_param.c_str(), "r")) == NULL)
			REPORT_ERROR(ERR_IO_NOTOPEN,
						 (String)"Prog_Project_Parameters::read: There is a problem "
						 "opening the file " + fn_proj_param);
		while (fgets(line, 200, fh_param) != NULL)
		{
			if (line[0] == 0)
				continue;
			if (line[0] == '#')
				continue;
			if (line[0] == '\n')
				continue;
			switch (lineNo)
			{
			case 0:
				proj_Xdim = textToInteger(firstToken(line));
				proj_Ydim = textToInteger(nextToken());
				lineNo = 1;
				break;
			case 1:
				// Angle file
				fn_angle = firstWord(line);
				if (fn_angle == "NULL")
					;
				else if (!fn_angle.exists())
					REPORT_ERROR(ERR_IO_NOTEXIST, (String)"Prog_Project_Parameters::read: "
								 "file " + fn_angle + " doesn't exist");
				lineNo = 2;
				break;
			case 2:
				// theta init
				auxstr = firstWord(line);
				if (strcmp(auxstr, "NULL") != 0)
				{
					enable_angle_range = true;
					rot_range.ang0 = textToFloat(auxstr);
					auxstr = nextToken();
					if (auxstr == NULL)
					{
						// Fixed mode
						rot_range.randomness = ANGLE_RANGE_DETERMINISTIC;
						rot_range.angF = rot_range.ang0;
						rot_range.samples = 1;
					}
					else
					{
						rot_range.angF = textToFloat(auxstr);
						rot_range.samples = textToInteger(nextToken());
						if (rot_range.ang0 == rot_range.angF)
							rot_range.samples = 1;
						rot_range.randomness = translate_randomness(nextToken());
					}
					lineNo = 3;
				}
				else
				{
					enable_angle_range = false;
					lineNo = 5;
				}
				break;
			case 3:
				tilt_range.ang0 = textToFloat(firstToken(line));
				auxstr = nextToken();
				if (auxstr == NULL)
				{
					// Fixed mode
					tilt_range.randomness = ANGLE_RANGE_DETERMINISTIC;
					tilt_range.angF = tilt_range.ang0;
					tilt_range.samples = 1;
				}
				else
				{
					tilt_range.angF = textToFloat(auxstr);
					tilt_range.samples = textToInteger(nextToken());
					if (tilt_range.ang0 == tilt_range.angF)
						tilt_range.samples = 1;
					tilt_range.randomness = translate_randomness(nextToken());
				}
				lineNo = 4;
				break;
			case 4:
				psi_range.ang0 = textToFloat(firstToken(line));
				auxstr = nextToken();
				if (auxstr == NULL)
				{
					// Fixed mode
					psi_range.randomness = ANGLE_RANGE_DETERMINISTIC;
					psi_range.angF = psi_range.ang0;
					psi_range.samples = 1;
				}
				else
				{
					psi_range.angF = textToFloat(auxstr);
					psi_range.samples = textToInteger(nextToken());
					if (psi_range.ang0 == psi_range.angF)
						psi_range.samples = 1;
					psi_range.randomness = translate_randomness(nextToken());
				}
				lineNo = 5;
				break;
			case 5:
				rot_range.Ndev = textToFloat(firstWord(line));
				auxstr = nextToken();
				if (auxstr != NULL)
					rot_range.Navg = textToFloat(auxstr);
				else
					rot_range.Navg = 0;
				lineNo = 6;
				break;
			case 6:
				tilt_range.Ndev = textToFloat(firstWord(line));
				auxstr = nextToken();
				if (auxstr != NULL)
					tilt_range.Navg = textToFloat(auxstr);
				else
					tilt_range.Navg = 0;
				lineNo = 7;
				break;
			case 7:
				psi_range.Ndev = textToFloat(firstWord(line));
				auxstr = nextToken();
				if (auxstr != NULL)
					psi_range.Navg = textToFloat(auxstr);
				else
					psi_range.Navg = 0;
				lineNo = 8;
				break;
			case 8:
				Npixel_dev = textToFloat(firstWord(line));
				auxstr = nextToken();
				if (auxstr != NULL)
					Npixel_avg = textToFloat(auxstr);
				else
					Npixel_avg = 0;
				lineNo = 9;
				break;
			case 9:
				Ncenter_dev = textToFloat(firstWord(line));
				auxstr = nextToken();
				if (auxstr != NULL)
					Ncenter_avg = textToFloat(auxstr);
				else
					Ncenter_avg = 0;
				lineNo = 10;
				break;
			} /* switch end */
		} /* while end */
		if (lineNo != 10)
			REPORT_ERROR(ERR_PARAM_MISSING, formatString("Prog_Project_Parameters::read: I "
						 "couldn't read all parameters from file %s, only read %d lines", fn_proj_param.c_str(), lineNo));
		fclose(fh_param);
	}
}

/* Generate angles ========================================================= */
// This function generates the angles for a given angle ("rot", "tilt"
// or "psi") according to the projection parameters. The output document
// file is supposed to be large enough to hold all angles
// Some aliases
#define Nrot  prm.rot_range.samples
#define Ntilt prm.tilt_range.samples
#define Npsi  prm.psi_range.samples
#define proj_number(base,irot,itilt,ipsi) base+irot*Ntilt*Npsi+itilt*Npsi+ipsi
void generate_angles(int ExtProjs, const Angle_range &range,
                     MetaData &DF, char ang_name, const ParametersProjection &prm)
{
    double ang;
    int   N1, N2;
    int   i, j, k;
    size_t iproj, idx;
    int   limit;
    MetaData DFaux=DF;

    // Select loop limit ....................................................
    switch (range.randomness)
    {
    case ANGLE_RANGE_DETERMINISTIC:
        limit = range.samples;
        break;
    case ANGLE_RANGE_RANDOM_GROUPS:
        limit = range.samples;
        break;
    case ANGLE_RANGE_RANDOM       :
        limit = Nrot * Ntilt * Npsi;
        break;
    }

    // Which column to write in the document file ...........................
    switch (ang_name)
    {
    case 'r':
        idx = 0;
        break;
    case 't':
        idx = 1;
        break;
    case 'p':
        idx = 2;
        break;
    }

    double unif_min = cos(DEG2RAD(range.angF));
    double unif_max = cos(DEG2RAD(range.ang0));
    for (i = 0; i < limit; i++)
    {
        // Select angle .....................................................
        if (range.randomness == ANGLE_RANGE_DETERMINISTIC)
        {
            if (range.samples > 1)
                ang = range.ang0 +
                      (range.angF - range.ang0) / (double)(range.samples - 1) * i;
            else
                ang = range.ang0;
        }
        else
        {
            switch (ang_name)
            {
            case 'r':
            case 'p':
                ang = rnd_unif(range.ang0, range.angF);
                break;
            case 't':
                ang = RAD2DEG(acos(rnd_unif(unif_min, unif_max)));
                break;
            }
        }

        // Copy this angle to those projections belonging to this group .....
        // If there is any group
        if (range.randomness != ANGLE_RANGE_RANDOM)
        {
            switch (ang_name)
            {
            case 'r':
                N1 = Ntilt;
                N2 = Npsi;
                break;
            case 't':
                N1 = Nrot;
                N2 = Npsi;
                break;
            case 'p':
                N1 = Nrot;
                N2 = Ntilt;
                break;
            }
            for (j = 0; j < N1; j++)
                for (k = 0; k < N2; k++)
                {
                    switch (ang_name)
                    {
                    case 'r':
                        iproj = proj_number(ExtProjs, i, j, k);
                        break;
                    case 't':
                        iproj = proj_number(ExtProjs, j, i, k);
                        break;
                    case 'p':
                        iproj = proj_number(ExtProjs, j, k, i);
                        break;
                    }
                    size_t idx_tmp=DFaux.firstObject(MDValueEQ(MDL_ORDER,iproj));
                    if (idx_tmp==BAD_OBJID)
                    {
                        idx_tmp=DFaux.addObject();
                        DFaux.setValue(MDL_ORDER,iproj,idx_tmp);
                    }
                    switch (idx)
                    {
                    case 0:
                        DFaux.setValue(MDL_ANGLEROT,ang,idx_tmp);
                        break;
                    case 1:
                        DFaux.setValue(MDL_ANGLETILT,ang,idx_tmp);
                        break;
                    case 2:
                        DFaux.setValue(MDL_ANGLEPSI,ang,idx_tmp);
                        break;
                    }
                }
        }
        else
        {
            size_t iproj=ExtProjs + i;
            size_t dfidx=DFaux.firstObject(MDValueEQ(MDL_ORDER,iproj));
            if (dfidx==BAD_OBJID)
            {
                dfidx=DFaux.addObject();
                DFaux.setValue(MDL_ORDER,iproj,dfidx);
            }
            switch (idx)
            {
            case 0:
                DFaux.setValue(MDL_ANGLEROT,ang,dfidx);
                break;
            case 1:
                DFaux.setValue(MDL_ANGLETILT,ang,dfidx);
                break;
            case 2:
                DFaux.setValue(MDL_ANGLEPSI,ang,dfidx);
                break;
            }
        }
    }
    DF.sort(DFaux,MDL_ORDER);
}

/* Generate evenly distributed angles ====================================== */
void generate_even_angles(int ExtProjs, int Nrottilt, MetaData &DF,
                          const ParametersProjection &prm)
{
    // We will run over the tilt angle in a deterministic way
    // then for every tilt angle, a rot_step is computed so that
    // it keeps the same distance in the circle generated by tilt
    // as the sample distance at the equator (tilt=90).
    MetaData DFaux;
    int N = 0;
    int limit = prm.tilt_range.samples;
    double rot_step_at_equator = (prm.rot_range.angF - prm.rot_range.ang0) /
                                 (double)(Nrot - 1);
    for (int i = 0; i < limit; i++)
    {
        // Compute the corresponding deterministic tilt
        double tilt = prm.tilt_range.ang0 +
                      (prm.tilt_range.angF - prm.tilt_range.ang0) /
                      (double)(Ntilt - 1) * i;
        // Now compute the corresponding rotational angles
        double rot_step;
        if (tilt != 0 && tilt != 180)
            rot_step = rot_step_at_equator / sin(DEG2RAD(tilt));
        else
            rot_step = prm.rot_range.angF - prm.rot_range.ang0 + 1;
        for (double rot = prm.rot_range.ang0; rot <= prm.rot_range.angF; rot += rot_step)
        {
            // Copy this angle to those projections belonging to this group .....
            // If there is any group
            for (int k = 0; k < Npsi; k++)
            {
                // Select psi
                double psi;
                if (prm.psi_range.randomness == ANGLE_RANGE_DETERMINISTIC)
                {
                    if (prm.psi_range.samples > 1)
                        psi = prm.psi_range.ang0 +
                              (prm.psi_range.angF - prm.psi_range.ang0) /
                              (double)(prm.psi_range.samples - 1) * k;
                    else
                        psi = prm.psi_range.ang0;
                }
                else
                    psi = rnd_unif(prm.psi_range.ang0, prm.psi_range.angF);

                size_t iproj = ExtProjs + N + Nrottilt * k;
                size_t idx_tmp=DFaux.firstObject(MDValueEQ(MDL_ORDER,iproj));
                if (idx_tmp==BAD_OBJID)
                {
                    idx_tmp=DFaux.addObject();
                    DFaux.setValue(MDL_ORDER,iproj,idx_tmp);
                }
                DFaux.setValue(MDL_ANGLEROT,rot,idx_tmp);
                DFaux.setValue(MDL_ANGLETILT,tilt,idx_tmp);
                DFaux.setValue(MDL_ANGLEPSI,psi,idx_tmp);
            }
            N++;
        }
    }
    DF.sort(DFaux,MDL_ORDER);
    DF.removeLabel(MDL_ORDER);
}

// See generate_even_angles for comments
int count_even_angles(const ParametersProjection &prm)
{
    int N = 0;
    int limit = prm.tilt_range.samples;
    double rot_step_at_equator = (prm.rot_range.angF - prm.rot_range.ang0) /
                                 (double)(Nrot - 1);
    for (int i = 0; i < limit; i++)
    {
        double tilt = prm.tilt_range.ang0 +
                      (prm.tilt_range.angF - prm.tilt_range.ang0) /
                      (double)(Ntilt - 1) * i;
        double rot_step;
        if (tilt != 0 && tilt != 180)
            rot_step = rot_step_at_equator / sin(DEG2RAD(tilt));
        else
            rot_step = prm.rot_range.angF - prm.rot_range.ang0 + 1;
        for (double rot = prm.rot_range.ang0; rot <= prm.rot_range.angF; rot += rot_step)
            N++;
    }
    N++; // This shouldn't be necessary but some GCC optimization
    // sometimes doesn't do well its work. For instance if we
    // add std::cout << N after N++ in the loop, then it works perfectly
    return N;
}

/* Assign angles =========================================================== */
int Assign_angles(MetaData &DF, const ParametersProjection &prm,
                  const FileName &fn_sym)
{
    int ExtProjs = 0, IntProjs = 0;    // External and internal projections
    int Nrottilt;                      // Number of evenly distributed

    // External generation mode
    if (prm.fn_angle != "NULL")
    {
        DF.read(prm.fn_angle);
        ExtProjs = DF.size();
    }

    // Internal generation mode
    if (prm.enable_angle_range)
    {
        randomize_random_generator();
        if (prm.rot_range.randomness != ANGLE_EVENLY)
            IntProjs = Nrot * Ntilt * Npsi;
        else
        {
            if (fn_sym == "")
            {
                Nrottilt = count_even_angles(prm);
                IntProjs = Nrottilt * Npsi;
            }
            else
                IntProjs = 0;
        }
        if (prm.rot_range.randomness != ANGLE_EVENLY)
        {
            generate_angles(ExtProjs, prm.rot_range,  DF, 'r', prm);
            generate_angles(ExtProjs, prm.tilt_range, DF, 't', prm);
            generate_angles(ExtProjs, prm.psi_range,  DF, 'p', prm);
            DF.removeLabel(MDL_ORDER);
        }
        else
        {
            if (fn_sym == "")
                generate_even_angles(ExtProjs, Nrottilt, DF, prm);
            else
            {
                SymList SL;
                if (fn_sym != "")
                    SL.readSymmetryFile(fn_sym);
                std::vector<double> rotList, tiltList;
                double rot_step_at_equator = (prm.rot_range.angF - prm.rot_range.ang0) /
                                             (double)(Nrot - 1);
                make_even_distribution(rotList, tiltList, rot_step_at_equator, SL, true);
                size_t NN=rotList.size();
                for (size_t n=0; n<NN; ++n)
                {
                    size_t DFid=DF.addObject();
                    DF.setValue(MDL_ANGLEROT,rotList[n],DFid);
                    DF.setValue(MDL_ANGLETILT,tiltList[n],DFid);
                    DF.setValue(MDL_ANGLEPSI,0.0,DFid);
                }
            }
        }
    }

    // Exit
    return ExtProjs + IntProjs;
}

/* Produce Side Information ================================================ */
void PROJECT_Side_Info::produce_Side_Info(ParametersProjection &prm,
        ProgProject &prog_prm)
{
    // Generate Projection angles
    if (!prog_prm.singleProjection)
        Assign_angles(DF, prm, prog_prm.fn_sym);
    else
    {
        size_t DFid=DF.addObject();
        DF.setValue(MDL_ANGLEROT,prog_prm.rotSingle,DFid);
        DF.setValue(MDL_ANGLETILT,prog_prm.tiltSingle,DFid);
        DF.setValue(MDL_ANGLEPSI,prog_prm.psiSingle,DFid);
    }

    // Load Phantom and set working mode
    if (prog_prm.fnPhantom.getExtension()=="descr")
    {
        phantomDescr.read(prog_prm.fnPhantom);
        phantomMode = XMIPP;
        if (prog_prm.singleProjection)
            if (prog_prm.projSize==-1)
                prm.proj_Xdim=prm.proj_Ydim=phantomDescr.xdim;
            else
                prm.proj_Xdim=prm.proj_Ydim=prog_prm.projSize;
    }
    else if (prog_prm.fnPhantom.getExtension()=="pdb")
    {
        phantomPDB.read(prog_prm.fnPhantom);
        const double highTs=1.0/12.0;
        int M=ROUND(prog_prm.samplingRate/highTs);
        interpolator.setup(M,prog_prm.samplingRate/M,true);
        phantomMode = PDB;
        if (prog_prm.singleProjection)
            if (prog_prm.projSize==-1)
                REPORT_ERROR(ERR_ARG_MISSING,"--xdim");
            else
                prm.proj_Xdim=prm.proj_Ydim=prog_prm.projSize;
    }
    else
    {
        phantomVol.read(prog_prm.fnPhantom);
        phantomVol().setXmippOrigin();
        phantomMode = VOXEL;
        if (prog_prm.singleProjection)
            if (prog_prm.projSize==-1)
                prm.proj_Xdim=prm.proj_Ydim=XSIZE(phantomVol());
            else
                prm.proj_Xdim=prm.proj_Ydim=prog_prm.projSize;
    }
}

/* Effectively project ===================================================== */
int PROJECT_Effectively_project(const String &fnOut,
                                bool singleProjection,
                                bool shears,
                                const ParametersProjection &prm,
                                PROJECT_Side_Info &side,
                                const Crystal_Projection_Parameters &prm_crystal,
                                Projection &proj, MetaData &SF)
{
    int NumProjs = 0;
    SF.clear();
    FileName(fnOut).deleteFile();
    std::cerr << "Projecting ...\n";
    init_progress_bar(side.DF.size());
    SF.setComment("First set of angles=actual angles; Second set of angles=noisy angles");

//#define DEBUG
#ifdef DEBUG

    MetaData mdShifts;
    MetaData mdRotations;
    /**
     * Here we create two auxiliary metadata files to check that the
     * Alignment information that is created by project is correct
     * We can not just use the metadata file created by project because this is
     * designed to be applied to the volume (during the reconstruction) but
     * not to the images.
     *
     * Helpful batch file:
     * xmipp_phantom_project  -i Phantom/asy.vol -o Images/proj_shifts.stk --params Phantom/asy2.param
    cp Images/proj_shifts.stk kk.stk
    xmipp_show Images/proj_shifts.stk
    xmipp_transform_geometry -i shifts.xmd  --apply_transform -o Images/kk.stk
    mv Images/kk.stk  Images/proj_shifts.stk
    xmipp_show Images/proj_shifts.stk
    xmipp_transform_geometry -i rotations.xmd  --apply_transform -o Images/kk.stk
    mv Images/kk.stk  Images/proj_shifts.stk
    xmipp_show Images/proj_shifts.stk
     *
     */
#endif

    int projIdx=FIRST_IMAGE;
    FileName fn_proj;              // Projection name
    RealShearsInfo *Vshears=NULL;
    if (shears && side.phantomMode==PROJECT_Side_Info::VOXEL)
        Vshears=new RealShearsInfo(side.phantomVol());
    FOR_ALL_OBJECTS_IN_METADATA(side.DF)
    {
        size_t DFmov_objId=SF.addObject();
        if (singleProjection)
            fn_proj=fnOut;
        else
            fn_proj.compose(projIdx,fnOut);
        SF.setValue(MDL_IMAGE,fn_proj,DFmov_objId);
        SF.setValue(MDL_ENABLED,1,DFmov_objId);

        // Choose angles .....................................................
        double rot, tilt, psi;         // Actual projecting angles
        side.DF.getValue(MDL_ANGLEROT,rot,__iter.objId);
        side.DF.getValue(MDL_ANGLETILT,tilt,__iter.objId);
        side.DF.getValue(MDL_ANGLEPSI,psi,__iter.objId);
        SF.setValue(MDL_ANGLEROT,rot,DFmov_objId);
        SF.setValue(MDL_ANGLETILT,tilt,DFmov_objId);
        SF.setValue(MDL_ANGLEPSI,psi,DFmov_objId);
        if ((NumProjs % XMIPP_MAX(1, side.DF.size() / 60)) == 0)
            progress_bar(NumProjs);

        // Choose Center displacement ........................................
        double shiftX = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        double shiftY = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        SF.setValue(MDL_SHIFTX,-shiftX,DFmov_objId);
        SF.setValue(MDL_SHIFTY,-shiftY,DFmov_objId);

#ifdef DEBUG

        mdShifts.addObject();
        mdShifts.setValue(MDL_IMAGE,fn_proj,DFmov_objId);
        mdShifts.setValue(MDL_SHIFTX,-shiftX,DFmov_objId);
        mdShifts.setValue(MDL_SHIFTY,-shiftY,DFmov_objId);
        mdShifts.setValue(MDL_ENABLED,1,DFmov_objId);


        mdRotations.addObject();
        mdRotations.setValue(MDL_IMAGE,fn_proj,DFmov_objId);
        mdRotations.setValue(MDL_ANGLEROT,-rot,DFmov_objId);
        mdRotations.setValue(MDL_ANGLETILT,-tilt,DFmov_objId);
        mdRotations.setValue(MDL_ANGLEPSI,-psi,DFmov_objId);


#endif
        // Really project ....................................................
        if (side.phantomMode==PROJECT_Side_Info::VOXEL)
        {
            if (shears)
                projectVolume(*Vshears, proj, prm.proj_Ydim, prm.proj_Xdim,
                              rot, tilt, psi);
            else
                projectVolume(side.phantomVol(), proj, prm.proj_Ydim, prm.proj_Xdim,
                              rot, tilt, psi);
            Matrix1D<double> shifts(2);
            XX(shifts) = shiftX;
            YY(shifts) = shiftY;
            selfTranslate(LINEAR,IMGMATRIX(proj), shifts);
        }
        else if (side.phantomMode==PROJECT_Side_Info::PDB)
        {
            PDBPhantom aux;
            aux=side.phantomPDB;
            aux.shift(shiftX,shiftY,0);
            projectPDB(side.phantomPDB, side.interpolator,
                       proj, prm.proj_Ydim, prm.proj_Xdim,
                       rot, tilt, psi);
        }
        else if (side.phantomMode==PROJECT_Side_Info::XMIPP)
        {
            Phantom aux;
            aux = side.phantomDescr;
            aux.shift(shiftX, shiftY, 0);
            if (prm_crystal.crystal_Xdim == 0)
            {
                // Project a single mathematical volume
                aux.project_to(proj,
                               prm.proj_Ydim, prm.proj_Xdim,
                               rot, tilt, psi);
            }
            else
            { // Project mathematical volume as a crystal
                project_crystal(aux, proj, prm, side, prm_crystal, rot, tilt, psi);
            }
        }

        // Add noise in angles and voxels ....................................
        SF.setValue(MDL_ANGLEROT2,rot,DFmov_objId);
        SF.setValue(MDL_ANGLETILT2,tilt,DFmov_objId);
        SF.setValue(MDL_ANGLEPSI2,psi,DFmov_objId);
        rot  += rnd_gaus(prm.rot_range.Navg,  prm.rot_range.Ndev);
        tilt += rnd_gaus(prm.tilt_range.Navg, prm.tilt_range.Ndev);
        psi  += rnd_gaus(prm.psi_range.Navg,  prm.psi_range.Ndev);
        SF.setValue(MDL_ANGLEROT,rot,DFmov_objId);
        SF.setValue(MDL_ANGLETILT,tilt,DFmov_objId);
        SF.setValue(MDL_ANGLEPSI,psi,DFmov_objId);
        IMGMATRIX(proj).addNoise(prm.Npixel_avg, prm.Npixel_dev, "gaussian");

        // Save ..............................................................
        if (singleProjection)
            proj.write(fn_proj);
        else
            proj.write(fn_proj,projIdx,true,WRITE_APPEND);
        projIdx++;
        NumProjs++;
    }
    progress_bar(side.DF.size());

#ifdef DEBUG

    mdShifts.write("shifts.xmd");
    mdRotations.write("rotations.xmd");
#endif
#undef DEBUG

    return NumProjs;
}

/* ROUT_project ============================================================ */
int ROUT_project(ProgProject &prm, Projection &proj, MetaData &SF)
{
    randomize_random_generator();
    // Read projection parameters and produce side information
    ParametersProjection proj_prm;
    PROJECT_Side_Info side;
    bool doCrystal = false;
    if (!prm.singleProjection)
        proj_prm.from_prog_params(prm);
    side.produce_Side_Info(proj_prm, prm);
    Crystal_Projection_Parameters crystal_proj_prm;
	MetaData MD;
	size_t objId;

	MD.read((std::string)"block1@"+prm.fn_proj_param);
	objId = MD.firstObject();
	doCrystal = MD.getValue(MDL_CRYSTAL_PROJ,doCrystal,objId);
	if (doCrystal)
    {
		crystal_proj_prm.read(prm.fn_proj_param,
    		                  (side.phantomDescr).phantom_scale);
        // if not null read doc file with unitcell shift
        // format h, k, shift_X shift_Y shift_Z
        if (crystal_proj_prm.DF_shift_bool == true)
            crystal_proj_prm.DF_shift.read(crystal_proj_prm.fn_shift);
        double my_scale = (side.phantomDescr).phantom_scale;
        FOR_ALL_OBJECTS_IN_METADATA(crystal_proj_prm.DF_shift)
        {
            double xcell, ycell;
            crystal_proj_prm.DF_shift.getValue(MDL_CELLX,xcell,__iter.objId);
            crystal_proj_prm.DF_shift.getValue(MDL_CELLY,ycell,__iter.objId);
            crystal_proj_prm.DF_shift.setValue(MDL_CELLX,xcell*my_scale,__iter.objId);
            crystal_proj_prm.DF_shift.setValue(MDL_CELLY,ycell*my_scale,__iter.objId);

            double x,y,z;
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFTX,x,__iter.objId);
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFTY,y,__iter.objId);
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFTZ,z,__iter.objId);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFTX,x*my_scale,__iter.objId);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFTY,y*my_scale,__iter.objId);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFTZ,z*my_scale,__iter.objId);

            crystal_proj_prm.DF_shift.getValue(MDL_SHIFT_CRYSTALX,x,__iter.objId);
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFT_CRYSTALY,y,__iter.objId);
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFT_CRYSTALZ,z,__iter.objId);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFT_CRYSTALX,x*my_scale,__iter.objId);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFT_CRYSTALY,y*my_scale,__iter.objId);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFT_CRYSTALZ,z*my_scale,__iter.objId);
        }
    }

    //#define DEBUG1
#ifdef DEBUG1
    if (crystal_proj_prm.DF_shift_bool == true)
        crystal_proj_prm.DF_shift.write("DEBUG1_shifts");
#endif
#undef DEBUG1

    int ProjNo = 0;
    if (!prm.only_create_angles)
    {
        // Really project
        if (prm.singleProjection)
            ProjNo = PROJECT_Effectively_project(prm.fnOut, prm.singleProjection, prm.shears,
                                                 proj_prm, side, crystal_proj_prm, proj, SF);
        else
        {
            FileName stackName;
            if (prm.fnOut.hasStackExtension())
            	stackName = prm.fnOut;
            else
            	stackName = prm.fnOut.removeAllExtensions() + ".stk";
            FileName mdName = prm.fnOut.removeAllExtensions() + ".xmd";
            ProjNo = PROJECT_Effectively_project(stackName, prm.singleProjection, prm.shears,
                                                 proj_prm, side, crystal_proj_prm, proj, SF);
            SF.setComment("Angles rot,tilt and psi contain noisy projection angles and rot2,tilt2 and psi2 contain actual projection angles");
            SF.write(mdName);
        }
    }
    else
        if (!prm.singleProjection)
            side.DF.write(prm.fnOut.removeAllExtensions()+".xmd");
    return ProjNo;
}
