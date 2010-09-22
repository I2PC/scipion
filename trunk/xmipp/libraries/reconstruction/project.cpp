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

#include <data/args.h>

/* Read from command line ================================================== */
void Prog_Project_Parameters::read(int argc, char **argv)
{
    fn_proj_param = getParameter(argc, argv, "-i");
    fn_sel_file   = getParameter(argc, argv, "-o", "");
    fn_crystal    = getParameter(argc, argv, "-crystal", "");
    fn_sym        = getParameter(argc, argv, "-sym", "");
    samplingRate  = textToFloat(
                        getParameter(argc, argv, "-sampling_rate", "1"));
    only_create_angles = checkParameter(argc, argv, "-only_create_angles");
    if (checkParameter(argc, argv, "-show_angles"))
        tell |= TELL_SHOW_ANGLES;
}

/* Usage =================================================================== */
void Prog_Project_Parameters::usage()
{
    printf("\nUsage:\n\n");
    printf("project -i <Parameters File> \n"
           "       [-o <sel_file>]\n"
           "       [-show_angles]\n"
           "       [-sym <sym_file>]\n"
           "       [-sampling_rate <Ts=1>\n"
           "       [-only_create_angles]\n"
           "       [-crystal <crystal_parameters_file>]\n");
    printf(
        "\tWhere:\n"
        "\t<Parameters File>:  File containing projection parameters\n"
        "\t                    check the manual for a description of the parameters\n"
        "\t<sel_file>:         This is a selection file with all the generated\n"
        "\t                    projections\n");
    printf(
        "\t<sym_file>:         This is a symmetry description file, used\n"
        "\t                    for computing the assymetric projection unit\n");
    printf(
        "\t<Ts>:               It is only used for PDB projections\n");
}

/* Projection parameters from program parameters =========================== */
void Projection_Parameters::from_prog_params(
    const Prog_Project_Parameters &prog_prm)
{
    read(prog_prm.fn_proj_param);
    tell = prog_prm.tell;
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
                 (std::string)"Prog_Project_Parameters::read: Not recognized randomness: "
                 + str);
}

void Projection_Parameters::read(const FileName &fn_proj_param)
{
    FILE    *fh_param;
    char    line[201];
    int     lineNo = 0;
    char    *auxstr;

    if ((fh_param = fopen(fn_proj_param.c_str(), "r")) == NULL)
        REPORT_ERROR(ERR_IO_NOTOPEN,
                     (std::string)"Prog_Project_Parameters::read: There is a problem "
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
            fnPhantom = firstWord(line);
            lineNo = 1;
            break;
        case 1:
            fnProjectionSeed =
                firstWord(line);
            // Next two parameters are optional
            auxstr = nextToken();
            if (auxstr != NULL)
                starting =
                    textToInteger(auxstr);
            fn_projection_extension = nextToken();
            lineNo = 2;
            break;
        case 2:
            proj_Xdim = textToInteger(firstToken(line));
            proj_Ydim = textToInteger(nextToken());
            lineNo = 3;
            break;
        case 3:
            // Angle file
            fn_angle = firstWord(line);
            if (fn_angle == "NULL")
                ;
            else if (!exists(fn_angle))
                REPORT_ERROR(ERR_IO_NOTEXIST, (std::string)"Prog_Project_Parameters::read: "
                             "file " + fn_angle + " doesn't exist");
            lineNo = 4;
            break;
        case 4:
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
                lineNo = 5;
            }
            else
            {
                enable_angle_range = false;
                lineNo = 7;
            }
            break;
        case 5:
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
            lineNo = 6;
            break;
        case 6:
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
            lineNo = 7;
            break;
        case 7:
            rot_range.Ndev = textToFloat(firstWord(line));
            auxstr = nextToken();
            if (auxstr != NULL)
                rot_range.Navg = textToFloat(auxstr);
            else
                rot_range.Navg = 0;
            lineNo = 8;
            break;
        case 8:
            tilt_range.Ndev = textToFloat(firstWord(line));
            auxstr = nextToken();
            if (auxstr != NULL)
                tilt_range.Navg = textToFloat(auxstr);
            else
                tilt_range.Navg = 0;
            lineNo = 9;
            break;
        case 9:
            psi_range.Ndev = textToFloat(firstWord(line));
            auxstr = nextToken();
            if (auxstr != NULL)
                psi_range.Navg = textToFloat(auxstr);
            else
                psi_range.Navg = 0;
            lineNo = 10;
            break;
        case 10:
            Npixel_dev = textToFloat(firstWord(line));
            auxstr = nextToken();
            if (auxstr != NULL)
                Npixel_avg = textToFloat(auxstr);
            else
                Npixel_avg = 0;
            lineNo = 11;
            break;
        case 11:
            Ncenter_dev = textToFloat(firstWord(line));
            auxstr = nextToken();
            if (auxstr != NULL)
                Ncenter_avg = textToFloat(auxstr);
            else
                Ncenter_avg = 0;
            lineNo = 12;
            break;
        } /* switch end */
    } /* while end */
    if (lineNo != 12)
        REPORT_ERROR(ERR_PARAM_MISSING, (std::string)"Prog_Project_Parameters::read: I "
                     "couldn't read all parameters from file " + fn_proj_param);
    fclose(fh_param);
}

/* Write =================================================================== */
void Projection_Parameters::write(const FileName &fn_proj_param) const
{
    FILE *fh_param;

    if ((fh_param = fopen(fn_proj_param.c_str(), "w")) == NULL)
        REPORT_ERROR(ERR_IO_NOTOPEN,
                     (std::string)"Prog_Project_Parameters::write: There is a problem "
                     "opening the file " + fn_proj_param + " for output");

    fprintf(fh_param,
            "# Volume and projection files -----------------------------------\n");
    fprintf(fh_param,
            "# volume description file or volume file\n");
    fprintf(fh_param,
            "%s\n", fnPhantom.c_str());
    fprintf(fh_param,
            "# projection seed, first projection number (by default, 1) and extension\n");
    fprintf(fh_param,
            "%s %d %s\n", fnProjectionSeed.c_str(), starting,
            fn_projection_extension.c_str());
    fprintf(fh_param,
            "# Y and X projection dimensions\n");
    fprintf(fh_param,
            "%d %d\n", proj_Ydim, proj_Xdim);
    fprintf(fh_param,
            "#\n");

    fprintf(fh_param,
            "# Angle Definitions ---------------------------------------------\n");
    if (fn_angle == "")
    {
        fprintf(fh_param, "%d ", rot_range.ang0);
        if (rot_range.angF != rot_range.ang0)
        {
            fprintf(fh_param, "%d %d ", rot_range.angF, rot_range.samples);
            switch (rot_range.randomness)
            {
            case (ANGLE_RANGE_RANDOM_GROUPS):
                            fprintf(fh_param, "random_group\n");
                break;
            case (ANGLE_RANGE_RANDOM):
                            fprintf(fh_param, "random\n");
                break;
            case (ANGLE_EVENLY):
                            fprintf(fh_param, "even\n");
                break;
            default:
                fprintf(fh_param, "\n");
                break;
            }
        }

        fprintf(fh_param, "%d ", tilt_range.ang0);
        if (tilt_range.angF != tilt_range.ang0)
{
            fprintf(fh_param, "%d %d ", tilt_range.angF, tilt_range.samples);
            switch (tilt_range.randomness)
            {
            case (ANGLE_RANGE_RANDOM_GROUPS):
                            fprintf(fh_param, "random_group\n");
                break;
            case (ANGLE_RANGE_RANDOM):
                            fprintf(fh_param, "random\n");
                break;
            case (ANGLE_EVENLY):
                            fprintf(fh_param, "even\n");
                break;
            default:
                fprintf(fh_param, "\n");
                break;
            }
        }

        fprintf(fh_param, "%d ", psi_range.ang0);
        if (psi_range.angF != psi_range.ang0)
{
            fprintf(fh_param, "%d %d ", psi_range.angF, psi_range.samples);
            switch (psi_range.randomness)
            {
            case (ANGLE_RANGE_RANDOM_GROUPS):
                            fprintf(fh_param, "random_group\n");
                break;
            case (ANGLE_RANGE_RANDOM):
                            fprintf(fh_param, "random\n");
                break;
            case (ANGLE_EVENLY):
                            fprintf(fh_param, "even\n");
                break;
            default:
                fprintf(fh_param, "\n");
                break;
            }
        }
    }
    fprintf(fh_param,
            "# Noise description ----------------------------------------------\n");
    fprintf(fh_param,
            "#     noise (and bias) applied to rotational angle\n");
    fprintf(fh_param, "%f ", rot_range.Ndev);
    if (rot_range.Navg != 0)
        fprintf(fh_param, "%f \n", rot_range.Navg);
    else
        fprintf(fh_param, "\n");

    fprintf(fh_param,
            "#     noise (and bias) applied to tilting angle\n");
    fprintf(fh_param, "%f ", tilt_range.Ndev);
    if (tilt_range.Navg != 0)
        fprintf(fh_param, "%f \n", tilt_range.Navg);
    else
        fprintf(fh_param, "\n");

    fprintf(fh_param,
            "#     noise (and bias) applied to psi angle\n");
    fprintf(fh_param, "%f ", psi_range.Ndev);
    if (psi_range.Navg != 0)
        fprintf(fh_param, "%f \n", psi_range.Navg);
    else
        fprintf(fh_param, "\n");

    fprintf(fh_param,
            "#     Noise (and bias) applied to pixels\n");
    fprintf(fh_param, "%f ", Npixel_dev);
    if (Npixel_avg != 0)
        fprintf(fh_param, "%f \n", Npixel_avg);
    else
        fprintf(fh_param, "\n");

    fprintf(fh_param,
            "#     Noise (and bias) applied to particle center coordenates\n");
    fprintf(fh_param, "%f ", Ncenter_dev);
    if (Ncenter_avg != 0)
        fprintf(fh_param, "%f \n", Ncenter_avg);
    else
        fprintf(fh_param, "\n");

    fclose(fh_param);
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
                     MetaData &DF, char ang_name, const Projection_Parameters &prm)
{
    double ang;
    int   N1, N2;
    int   i, j, k;
    long int   iproj, idx;
    int   limit;

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
                        iproj = proj_number(ExtProjs, i, j, k)+prm.starting;
                        break;
                    case 't':
                        iproj = proj_number(ExtProjs, j, i, k)+prm.starting;
                        break;
                    case 'p':
                        iproj = proj_number(ExtProjs, j, k, i)+prm.starting;
                        break;
                    }
                    long int idx_tmp=DF.gotoFirstObject(MDValueEQ(MDL_OBJID,iproj));
                    if (idx_tmp==NO_OBJECTS_STORED || idx_tmp==NO_OBJECT_FOUND)
                    {
                        DF.addObject(iproj);
                        //DF.setValue(MDL_OBJID,iproj);
                    }
                    switch (idx)
                    {
                    case 0:
                        DF.setValue(MDL_ANGLEROT,ang);
                        break;
                    case 1:
                        DF.setValue(MDL_ANGLETILT,ang);
                        break;
                    case 2:
                        DF.setValue(MDL_ANGLEPSI,ang);
                        break;
                    }
                }
        }
        else
        {
            long int iproj=ExtProjs + i + prm.starting;
            unsigned long int dfidx=DF.gotoFirstObject(MDValueEQ(MDL_OBJID,iproj));
            if (dfidx==NO_OBJECTS_STORED || dfidx==NO_OBJECT_FOUND)
            {
                DF.addObject(dfidx);
                //DF.setValue(MDL_OBJID,iproj);
            }
            switch (idx)
            {
            case 0:
                DF.setValue(MDL_ANGLEROT,ang);
                break;
            case 1:
                DF.setValue(MDL_ANGLETILT,ang);
                break;
            case 2:
                DF.setValue(MDL_ANGLEPSI,ang);
                break;
            }
        }
    }
}

/* Generate evenly distributed angles ====================================== */
void generate_even_angles(int ExtProjs, int Nrottilt, MetaData &DF,
                          const Projection_Parameters &prm)
{
    // We will run over the tilt angle in a deterministic way
    // then for every tilt angle, a rot_step is computed so that
    // it keeps the same distance in the circle generated by tilt
    // as the sample distance at the equator (tilt=90).
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

                long int iproj = ExtProjs + N + Nrottilt * k;
                unsigned long int idx_tmp=DF.gotoFirstObject(MDValueEQ(MDL_OBJID,iproj));
                if (idx_tmp==NO_OBJECTS_STORED || idx_tmp==NO_OBJECT_FOUND)
                {
                    DF.addObject(iproj);
                    //DF.setValue(MDL_OBJID,iproj);
                }
                DF.setValue(MDL_ANGLEROT,rot);
                DF.setValue(MDL_ANGLETILT,tilt);
                DF.setValue(MDL_ANGLEPSI,psi);
            }
            N++;
        }
    }
}

// See generate_even_angles for comments
int count_even_angles(const Projection_Parameters &prm)
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
int Assign_angles(MetaData &DF, const Projection_Parameters &prm,
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
            //DF.write("a.db");
            //MDSql::dumpToFile()
        }
        else
        {
            if (fn_sym == "")
                generate_even_angles(ExtProjs, Nrottilt, DF, prm);
            else
            {
                SymList SL;
                if (fn_sym != "")
                    SL.read_sym_file(fn_sym);
                MetaData DF_aux;
                double rot_step_at_equator = (prm.rot_range.angF - prm.rot_range.ang0) /
                                             (double)(Nrot - 1);
                make_even_distribution(DF_aux, rot_step_at_equator, SL, true);
                FOR_ALL_OBJECTS_IN_METADATA(DF_aux)
                {
                    double rot,tilt,psi;
                    DF_aux.getValue(MDL_ANGLEROT,rot);
                    DF_aux.getValue(MDL_ANGLETILT,tilt);
                    DF_aux.getValue(MDL_ANGLEPSI,psi);
                    DF.addObject();
                    DF.setValue(MDL_ANGLEROT,rot);
                    DF.setValue(MDL_ANGLETILT,tilt);
                    DF.setValue(MDL_ANGLEPSI,psi);
                }
            }
        }
    }

    // Exit
    return ExtProjs + IntProjs;
}

/* Produce Side Information ================================================ */
void PROJECT_Side_Info::produce_Side_Info(const Projection_Parameters &prm,
        const Prog_Project_Parameters &prog_prm)
{
    // Generate Projection angles
    Assign_angles(DF, prm, prog_prm.fn_sym);

    // Load Phantom and set working mode
    if (prm.fnPhantom.getExtension()=="descr")
    {
        phantomDescr.read(prm.fnPhantom);
        phantomMode = XMIPP;
    }
    else if (prm.fnPhantom.getExtension()=="pdb")
    {
        phantomPDB.read(prm.fnPhantom);
        const double highTs=1.0/12.0;
        int M=ROUND(prog_prm.samplingRate/highTs);
        interpolator.setup(M,prog_prm.samplingRate/M,true);
        phantomMode = PDB;
    }
    else
    {
        phantomVol.read(prm.fnPhantom);
        phantomVol().setXmippOrigin();
        phantomMode = VOXEL;
    }
}

/* Effectively project ===================================================== */
int PROJECT_Effectively_project(const Projection_Parameters &prm,
                                PROJECT_Side_Info &side,
                                const Crystal_Projection_Parameters &prm_crystal,
                                Projection &proj, MetaData &SF)
{
    int NumProjs = 0;
    SF.clear();
    std::cerr << "Projecting ...\n";
    if (!(prm.tell&TELL_SHOW_ANGLES))
        init_progress_bar(side.DF.size());
    MetaData DF_movements;
    DF_movements.setComment("First set of angles=actual angles; Second set of angles=noisy angles");

    int projIdx=prm.starting;
    FOR_ALL_OBJECTS_IN_METADATA(side.DF)
    {
        double rot, tilt, psi;         // Actual projecting angles
        FileName fn_proj;              // Projection name
        fn_proj.compose(prm.fnProjectionSeed, projIdx++, prm.fn_projection_extension);
        DF_movements.addObject(side.DF.getActiveObject());
        DF_movements.setValue(MDL_IMAGE,fn_proj);
        DF_movements.setValue(MDL_ENABLED,1);

        // Choose angles .....................................................
        side.DF.getValue(MDL_ANGLEROT,rot);
        side.DF.getValue(MDL_ANGLETILT,tilt);
        side.DF.getValue(MDL_ANGLEPSI,psi);
        DF_movements.setValue(MDL_ANGLEROT,rot);
        DF_movements.setValue(MDL_ANGLETILT,tilt);
        DF_movements.setValue(MDL_ANGLEPSI,psi);
        if (prm.tell&TELL_SHOW_ANGLES)
        {
            std::cout << rot << "\t" << tilt << "\t" << psi << std::endl;
        }
        else if ((NumProjs % XMIPP_MAX(1, side.DF.size() / 60)) == 0)
            progress_bar(NumProjs);

        // Choose Center displacement ........................................
        double shiftX = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        double shiftY = rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
        DF_movements.setValue(MDL_SHIFTX,-shiftX);
        DF_movements.setValue(MDL_SHIFTY,-shiftY);

        // Really project ....................................................
        if (side.phantomMode==PROJECT_Side_Info::VOXEL)
        {
            project_Volume(side.phantomVol(), proj, prm.proj_Ydim, prm.proj_Xdim,
                           rot, tilt, psi);
            selfTranslate(LINEAR,IMGMATRIX(proj),vectorR2(shiftX, shiftY));
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
        rot  += rnd_gaus(prm.rot_range.Navg,  prm.rot_range.Ndev);
        tilt += rnd_gaus(prm.tilt_range.Navg, prm.tilt_range.Ndev);
        psi  += rnd_gaus(prm.psi_range.Navg,  prm.psi_range.Ndev);
        DF_movements.setValue(MDL_ANGLEROT2,rot);
        DF_movements.setValue(MDL_ANGLETILT2,tilt);
        DF_movements.setValue(MDL_ANGLEPSI2,psi);
        proj.setEulerAngles(rot, tilt, psi);
        proj.setShifts(-shiftX,-shiftY);
        IMGMATRIX(proj).addNoise(prm.Npixel_avg, prm.Npixel_dev, "gaussian");

        // Save ..............................................................
        proj.write(fn_proj);
        NumProjs++;
        SF.addObject();
        SF.setValue(MDL_IMAGE,fn_proj);
        SF.setValue(MDL_ENABLED,1);
    }
    if (!(prm.tell&TELL_SHOW_ANGLES))
        progress_bar(side.DF.size());

    DF_movements.write(prm.fnProjectionSeed + "_movements.txt");
    return NumProjs;
}

/* ROUT_project ============================================================ */
int ROUT_project(Prog_Project_Parameters &prm, Projection &proj, MetaData &SF)
{
    randomize_random_generator();
    // Read projection parameters and produce side information
    Projection_Parameters proj_prm;
    PROJECT_Side_Info side;
    proj_prm.from_prog_params(prm);
    side.produce_Side_Info(proj_prm, prm);
    Crystal_Projection_Parameters crystal_proj_prm;

    if (prm.fn_crystal != "")
    {
        crystal_proj_prm.read(prm.fn_crystal,
                              (side.phantomDescr).phantom_scale);

        // if not null read doc file with unitcell shift
        // format h, k, shift_X shift_Y shift_Z
        if (crystal_proj_prm.DF_shift_bool == true)
            crystal_proj_prm.DF_shift.read(crystal_proj_prm.fn_shift);
        double my_scale = (side.phantomDescr).phantom_scale;
        FOR_ALL_OBJECTS_IN_METADATA(crystal_proj_prm.DF_shift)
        {
            double xcell, ycell;
            crystal_proj_prm.DF_shift.getValue(MDL_CELLX,xcell);
            crystal_proj_prm.DF_shift.getValue(MDL_CELLY,ycell);
            crystal_proj_prm.DF_shift.setValue(MDL_CELLX,xcell*my_scale);
            crystal_proj_prm.DF_shift.setValue(MDL_CELLY,ycell*my_scale);

            double x,y,z;
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFTX,x);
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFTY,y);
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFTZ,z);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFTX,x*my_scale);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFTY,y*my_scale);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFTZ,z*my_scale);

            crystal_proj_prm.DF_shift.getValue(MDL_SHIFT_CRYSTALX,x);
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFT_CRYSTALY,y);
            crystal_proj_prm.DF_shift.getValue(MDL_SHIFT_CRYSTALZ,z);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFT_CRYSTALX,x*my_scale);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFT_CRYSTALY,y*my_scale);
            crystal_proj_prm.DF_shift.setValue(MDL_SHIFT_CRYSTALZ,z*my_scale);
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
        ProjNo = PROJECT_Effectively_project(proj_prm, side, crystal_proj_prm,
                                             proj, SF);
        // Save SelFile
        if (prm.fn_sel_file != "")
            SF.write(prm.fn_sel_file);
    }
    else
    {
        side.DF.write("/dev/stdout");
    }
    return ProjNo;
}
