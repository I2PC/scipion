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

#include <data/args.h>
#include <data/psf_xr.h>

/* Read from command line ================================================== */
void Prog_Project_XR_Parameters::read(int argc, char **argv)
{
    if (argc==0)
    {
        usage();
        exit(0);
    }

    verbose = checkParameter(argc, argv, "-v");
    if (verbose)
        fn_proj_param = getParameter(argc, argv, "-i", "");
    else
        fn_proj_param = getParameter(argc, argv, "-i");
    fn_sel_file   = getParameter(argc, argv, "-o", "");
    fn_psf_xr   = getParameter(argc, argv, "-psf", "");
    only_create_angles = checkParameter(argc, argv, "-only_create_angles");
    tell=0;
    if (checkParameter(argc, argv, "-show_angles"))
        tell |= TELL_SHOW_ANGLES;
}

/* Usage =================================================================== */
void Prog_Project_XR_Parameters::usage()
{
    printf("\nUsage:\n\n");
    printf("project_xr -i <Parameters File> \n"
           "       [-psf <psf_param_file>]\n"
           "       [-o <sel_file>]\n"
           "       [-show_angles]\n"
           "       [-v]\n"
           "       [-only_create_angles]\n");
    printf(
        "\tWhere:\n"
        "\t<Parameters File>:  File containing projection parameters\n"
        "\t                    check the manual for a description of the parameters\n"
        "\t<psf_param_file>:  File containing X-ray microscope parameters\n"
        "\t<sel_file>:         This is a selection file with all the generated\n"
        "\t                    projections\n");
}

/* Read Projection Parameters ============================================== */
void Projection_XR_Parameters::read(const FileName &fn_proj_param)
{
    FILE    *fh_param;
    char    line[201];
    int     lineNo = 0;
    char    *auxstr;

    if (fn_proj_param=="")
        REPORT_ERROR(ERR_ARG_MISSING,
                     (std::string)"Projection_XR_Parameters::read: There is no parameters file.");


    if ((fh_param = fopen(fn_proj_param.c_str(), "r")) == NULL)
        REPORT_ERROR(ERR_IO_NOTOPEN,
                     (std::string)"Projection_XR_Parameters::read: There is a problem "
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
            axisRot = textToFloat(firstToken(line));
            axisTilt = textToFloat(nextToken());
            lineNo = 4;
            break;
        case 4:
            raxis.resize(3);
            XX(raxis) = textToFloat(firstToken(line));
            YY(raxis) = textToFloat(nextToken());
            ZZ(raxis) = textToFloat(nextToken());
            lineNo = 5;
            break;
        case 5:
            tilt0 = textToFloat(firstToken(line));
            tiltF = textToFloat(nextToken());
            tiltStep = textToFloat(nextToken());
            lineNo = 6;
            break;
        case 6:
            Nangle_dev = textToFloat(firstWord(line));
            auxstr = nextToken();
            if (auxstr != NULL)
                Nangle_avg = textToFloat(auxstr);
            else
                Nangle_avg = 0;
            lineNo = 7;
            break;
        case 7:
            Npixel_dev = textToFloat(firstWord(line));
            auxstr = nextToken();
            if (auxstr != NULL)
                Npixel_avg = textToFloat(auxstr);
            else
                Npixel_avg = 0;
            lineNo = 8;
            break;
        case 8:
            Ncenter_dev = textToFloat(firstWord(line));
            auxstr = nextToken();
            if (auxstr != NULL)
                Ncenter_avg = textToFloat(auxstr);
            else
                Ncenter_avg = 0;
            lineNo = 9;
            break;
        } /* switch end */
    } /* while end */
    if (lineNo != 9)
        REPORT_ERROR(ERR_PARAM_MISSING, (std::string)"Projection_XR_Parameters::read: I "
                     "couldn't read all parameters from file " + fn_proj_param);
    fclose(fh_param);
}


void Projection_XR_Parameters::calculateProjectionAngles(Projection &P, double angle, double inplaneRot,
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
    const Projection_XR_Parameters &prm,
    const Prog_Project_XR_Parameters &prog_prm)
{
    phantomVol.read(prm.fnPhantom);
    phantomVol().setXmippOrigin();
    rotPhantomVol = phantomVol;
}

/* Effectively project ===================================================== */
int PROJECT_XR_Effectively_project(Projection_XR_Parameters &prm,
                                   PROJECT_XR_Side_Info &side,
                                   Projection &proj,
                                   XRayPSF &psf,
                                   MetaData &SF)
{
    int expectedNumProjs = FLOOR((prm.tiltF-prm.tilt0)/prm.tiltStep);
    int numProjs=0;
    SF.clear();
    std::cerr << "Projecting ...\n";
    if (!(prm.tell&TELL_SHOW_ANGLES))
        init_progress_bar(expectedNumProjs*side.phantomVol().zdim);
    //    SF.reserve(expectedNumProjs);
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

        // Really project ....................................................
        project_xr_Volume_offCentered(side, psf, proj,prm.proj_Ydim, prm.proj_Xdim, idx);

        // Add noise in angles and voxels ....................................
        proj.getEulerAngles(tRot, tTilt,tPsi);

        rot  = tRot  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        tilt = tTilt + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);
        psi  = tPsi  + rnd_gaus(prm.Nangle_avg,  prm.Nangle_dev);

        proj.setEulerAngles(rot,tilt,psi);

        DF_movements.addObject();
        DF_movements.setValue(MDL_ANGLEROT,tRot);
        DF_movements.setValue(MDL_ANGLETILT,tTilt);
        DF_movements.setValue(MDL_ANGLEPSI,tPsi);
        DF_movements.setValue(MDL_ANGLEROT2,rot);
        DF_movements.setValue(MDL_ANGLETILT2,tilt);
        DF_movements.setValue(MDL_ANGLEPSI2,psi);
        DF_movements.setValue(MDL_SHIFTX,shiftX);
        DF_movements.setValue(MDL_SHIFTY,shiftY);

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
        SF.addObject();
        SF.setValue(MDL_IMAGE,fn_proj);
        SF.setValue(MDL_ENABLED,1);
    }
    if (!(prm.tell&TELL_SHOW_ANGLES))
        progress_bar(expectedNumProjs);

    DF_movements.write(prm.fnProjectionSeed + "_movements.txt");
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

    //    xOffsetN = XX(roffset)/psf.dxo;
    //    yOffsetN = YY(roffset)/psf.dxo;

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

#if DEBUG

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

    std::vector<MDLabel> labels = P.MD.getActiveLabels() ;


    //the really really final project routine, I swear by Snoopy.
    project_xr(psf,side.rotPhantomVol,P, idxSlice);


    int outXDim = XMIPP_MIN(Xdim,iniXdim);
    int outYDim = XMIPP_MIN(Ydim,iniYdim);

    P().window(-ROUND(outYDim/2),
               -ROUND(outXDim/2),
               -ROUND(outYDim/2) + outYDim -1,
               -ROUND(outXDim/2) + outXDim -1);

}




/* ROUT_project ============================================================ */
int ROUT_XR_project(Prog_Project_XR_Parameters &prm,
                    Projection &proj, MetaData &SF)
{
    randomize_random_generator();

    // Read Microscope optics parameters and produce side information
    XRayPSF psf;
    psf.verbose = prm.verbose;
    psf.read(prm.fn_psf_xr);
    psf.produceSideInfo();
    if(psf.verbose)
        psf.show();


    // Read projection parameters and produce side information
    Projection_XR_Parameters proj_prm;
    proj_prm.read(prm.fn_proj_param);
    proj_prm.tell=prm.tell;
    PROJECT_XR_Side_Info side;
    side.produce_Side_Info(proj_prm, prm);


    //    psf.adjustParam(side.phantomVol);

    // Project
    int ProjNo = 0;
    if (!prm.only_create_angles)
    {
        // Really project
        ProjNo = PROJECT_XR_Effectively_project(proj_prm, side,
                                                proj, psf, SF);
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
