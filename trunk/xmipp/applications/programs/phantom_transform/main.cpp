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

#include <data/progs.h>
#include <data/args.h>
#include <data/geometry.h>
#include <reconstruction/phantom.h>
#include <interface/pdb.h>

class Phantom_transform_parameters
{
public:
    FileName fn_in, fn_out;
    bool Euler_mode, PDBmode, centerPDB;
    double rot, tilt, psi;
    bool Align_mode;
    bool Axis_mode;
    Matrix1D<double> axis;
    double ang;
    Matrix1D<double> shift;
    Matrix1D<double> scale;

    Matrix2D<double> A3D;

    void read(int argc, char **argv)
    {
        fn_in = getParameter(argc, argv, "-i");
        fn_out = getParameter(argc, argv, "-o");
	PDBmode=checkParameter(argc, argv, "-pdb");
	centerPDB=checkParameter(argc,argv, "-centerPDB");
        Euler_mode = Align_mode = Axis_mode = false;
        if (checkParameter(argc, argv, "-euler"))
        {
            Euler_mode = true;
            int i = paremeterPosition(argc, argv, "-euler");
            if (i + 3 >= argc)
                REPORT_ERROR(1, "Not enough parameters after -euler");
            rot  = textToFloat(argv[i+1]);
            tilt = textToFloat(argv[i+2]);
            psi  = textToFloat(argv[i+3]);
            A3D = Euler_rotation3DMatrix(rot, tilt, psi);
        }
        else if (checkParameter(argc, argv, "-alignWithZ"))
        {
            Align_mode = true;
            axis = getVectorParameter(argc, argv, "-alignWithZ", 3);
            A3D = alignWithZ(axis);
        }
        else if (checkParameter(argc, argv, "-axis"))
        {
            Axis_mode = true;
            axis = getVectorParameter(argc, argv, "-axis", 3);
            ang = textToFloat(getParameter(argc, argv, "-ang"));
            A3D = rotation3DMatrix(ang, axis);
        }
        else
            A3D.initIdentity(4);
        if (checkParameter(argc, argv, "-shift"))
            shift = getVectorParameter(argc, argv, "-shift", 3);
        else shift.initZeros(3);
        if (checkParameter(argc, argv, "-scale"))
            scale = getVectorParameter(argc, argv, "-scale", 3);
        else
        {
            scale.resize(3);
            scale.init_constant(1);
        }

        // Apply shift
        A3D(0, 3) = XX(shift);
        A3D(1, 3) = YY(shift);
        A3D(2, 3) = ZZ(shift);

        // Apply scale
        A3D(0, 0) *= XX(scale);
        A3D(0, 1) *= YY(scale);
        A3D(0, 2) *= ZZ(scale);
        A3D(1, 0) *= XX(scale);
        A3D(1, 1) *= YY(scale);
        A3D(1, 2) *= ZZ(scale);
        A3D(2, 0) *= XX(scale);
        A3D(2, 1) *= YY(scale);
        A3D(2, 2) *= ZZ(scale);
    }

    void show()
    {
        std::cout << "Input file : " << fn_in  << std::endl
                  << "Output file: " << fn_out << std::endl;
	if (PDBmode) std::cout << "Input file is PDB\n";
        if (Euler_mode)
            std::cout << "Euler angles (rot, tilt, psi): " << rot << " " << tilt
                      << " " << psi << std::endl;
        else if (Align_mode)
            std::cout << "Aligning " << axis.transpose() << " with Z\n";
        else if (Axis_mode)
            std::cout << "Rotating " << ang << " degrees around "
	              << axis.transpose() << std::endl;
    }

    void usage()
    {
        std::cerr << "Usage: phantom_transform [Options]\n"
        	  << "   -i <input filename>              : Phantom description file\n"
        	  << "   -o <output filename>             : Phantom description file\n"
		  << "  [-pdb]                            : Use PDBs instead of phantom descriptions\n"
		  << "  [-center_pdb]                     : Substract the center of mass from coordinates\n"
        	  << "  [-euler <rot> <tilt> <psi>        : Rotate with these Euler angles\n"
        	  << "  [-alignWithZ [<x>,<y>,<z>]]     : Align (x,y,z) with Z\n"
        	  << "                                      Notice that brackets for the\n"
        	  << "                                      vector must be written and do not\n"
        	  << "                                      represent optional parameters\n"
        	  << "  [-axis [<x>,<y>,<z>] -ang <ang>]  : Rotate <ang> degrees around (x,y,z),\n"
        	  << "                                      by default (0,0,1)\n"
        	  << "  [-shift [<x>,<y>,<z>]             : Shift vector\n"
         	  << "  [-scale [<x>,<y>,<z>]             : Scale vector\n"
        ;
    }
};

bool process_phantom(const FileName &fn_in, const FileName &fn_out,
                     const Phantom_transform_parameters *prm)
{
    Phantom P;
    P.read(fn_in, false); // Read phantom without applying scale
    P.selfApplyGeometry(prm->A3D, IS_NOT_INV);
    P.write(fn_out);
    return true;
}

int main(int argc, char **argv)
{
    Phantom_transform_parameters prm;
    try
    {
        prm.read(argc, argv);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(1);
    }
    try
    {
        prm.show();
	if (prm.PDBmode)
	    applyGeometry(prm.fn_in, prm.fn_out, prm.A3D, prm.centerPDB);
	else
            process_phantom(prm.fn_in, prm.fn_out, &prm);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE << std::endl;
    }
    return 0;
}
