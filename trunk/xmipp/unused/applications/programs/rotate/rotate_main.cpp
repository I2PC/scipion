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

#include <data/progs.h>
#include <data/args.h>
#include <data/geometry.h>
#include <data/metadata.h>

class Rotate_parameters: public Prog_parameters
{
public:
    bool Euler_mode;
    double rot, tilt, psi;
    bool Align_mode;
    bool Axis_mode;
    Matrix1D<double> axis;
    double ang;
    bool wrap;
    bool linear;
    bool write_matrix;

    Matrix2D<double> A3D, A2D;
    // Also allow rotation of docfiles
    FileName fn_DFin, fn_DFout;
    MetaData DF;
    int col_rot, col_tilt, col_psi;

    void read(int argc, char **argv)
    {
        // Do not read Prog params for -doc options
        if (!checkParameter(argc, argv, "-doc"))
            Prog_parameters::read(argc, argv);

        Euler_mode = Align_mode = Axis_mode = false;
        if (checkParameter(argc, argv, "-euler"))
        {
            Euler_mode = true;
            int i = paremeterPosition(argc, argv, "-euler");
            if (i + 3 >= argc)
                REPORT_ERROR(ERR_ARG_MISSING, "Not enough parameters after -euler");
            rot  = textToFloat(argv[i+1]);
            tilt = textToFloat(argv[i+2]);
            psi  = textToFloat(argv[i+3]);
            Euler_angles2matrix(rot, tilt, psi, A3D, true);
        }
        else if (checkParameter(argc, argv, "-alignWithZ"))
        {
            Align_mode = true;
            axis = getVectorParameter(argc, argv, "-alignWithZ", 3);
            alignWithZ(axis, A3D);
        }
        else
        {
            Axis_mode = true;
            if (checkParameter(argc, argv, "-axis"))
                axis = getVectorParameter(argc, argv, "-axis", 3);
            else
                axis = vectorR3(0., 0., 1.);
            ang = textToFloat(getParameter(argc, argv, "-ang"));
            rotation3DMatrix(ang, axis, A3D);

            A2D.resize(3,3);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(A2D)
            {
            	dMij(A2D,i,j) = dMij(A3D,i,j);
            }
        }
        wrap = !checkParameter(argc, argv, "-dont_wrap");
        linear = checkParameter(argc, argv, "-linear");
        write_matrix = checkParameter(argc, argv, "-write_matrix");

        if (checkParameter(argc, argv, "-inverse"))
        {
            A3D=A3D.inv();
            A2D=A2D.inv();
        }

        fn_DFin = getParameter(argc, argv, "-doc","");
        fn_DFout = getParameter(argc, argv, "-o","");
    }

    void show()
    {
        Prog_parameters::show();
        if (Euler_mode)
            std::cout << "Euler angles (rot, tilt, psi): " << rot << " " << tilt
            << " " << psi << std::endl;
        else if (Align_mode)
            std::cout << "Aligning " << axis.transpose() << " with Z\n";
        else if (Axis_mode)
            std::cout << "Rotating " << ang << " degrees around " << axis.transpose()
            << std::endl;
        if (!wrap)
            std::cout << "Do not wrap."<<std::endl;
        if (linear)
            std::cout << "Use linear for interpolation."<<std::endl;
        if (write_matrix)
            std::cout << "Write out transformation matrix to the screen."<<std::endl;

    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-euler <rot> <tilt> <psi>        : Rotate with these Euler angles\n"
        << "  [-alignWithZ \"[<x>,<y>,<z>]\"]     : Align (x,y,z) with Z\n"
        << "                                      Notice that brackets for the\n"
        << "                                      vector must be written and do not\n"
        << "                                      represent optional parameters\n"
        << "  [-axis \"[<x>,<y>,<z>]\" -ang <ang>]: Rotate <ang> degrees around (x,y,z),\n"
        << "                                      by default (0,0,1)\n"
        << "  [-inverse ]                       : Use the inverse rotation \n"
        << "  [-dont_wrap]                      : By default, the image/volume is wrapped\n"
        << "  [-write_matrix]                   : Print transformation matrix to screen\n"
        << "  [-linear]                         : Use linear interpolation (instead of B-splines)\n"
        << "\n"
        << " OR rather than rotating image/volume(s), rotate all angles in a docfile\n"
        << "  [-doc <docfile>]                  : Input docfile \n"
        << "  [-o <output docfile>]             : Output docfile \n";

    }

    // Rotate all angles in a docfile
    void rotateAnglesInDocFile()
    {

        Matrix2D< double > I(3,3);
        I.initIdentity();
        A3D.resize(3,3);

        double rot, tilt, psi, newrot, newtilt, newpsi;

        DF.read(fn_DFin);
        DF.firstObject();
        FOR_ALL_OBJECTS_IN_METADATA(DF)
        {
            DF.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
            DF.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
            DF.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
            Euler_apply_transf(A3D, I, rot, tilt, psi, newrot, newtilt, newpsi);
            DF.setValue(MDL_ANGLE_ROT, newrot, __iter.objId);
            DF.setValue(MDL_ANGLE_TILT, newtilt, __iter.objId);
            DF.setValue(MDL_ANGLE_PSI, newpsi, __iter.objId);
        }
        std::cerr<<" Written output docfile "<<fn_DFout<<std::endl;
    }


};

bool process_img(Image<double> &img, const Prog_parameters *prm)
{
    Rotate_parameters *eprm = (Rotate_parameters *) prm;
    Matrix2D<double> AA;
    if (img().getDim()==2)
        AA=eprm->A2D;
    else if (img().getDim()==3)
        AA=eprm->A3D;
    else
        REPORT_ERROR(ERR_MULTIDIM_DIM,"Cannot rotate this image, wrong dimension...");

    if (eprm->write_matrix)
        std::cerr<<"Transformation matrix = "<<AA<<std::endl;
    if (eprm->linear)
        selfApplyGeometry(LINEAR, img(), AA,  IS_NOT_INV, eprm->wrap);
    else
        selfApplyGeometry(BSPLINE3, img(), AA,  IS_NOT_INV, eprm->wrap);

    return true;
}

int main(int argc, char **argv)
{

    Rotate_parameters prm;

    // Also allow rotation of the angles in a docfile.
    if (checkParameter(argc, argv, "-doc"))
    {
        prm.read(argc, argv);
        prm.rotateAnglesInDocFile();
    }
    else
    {
        // Normal rotation of images and volumes
        SF_main(argc, argv, &prm, (void*)&process_img);
    }
}
