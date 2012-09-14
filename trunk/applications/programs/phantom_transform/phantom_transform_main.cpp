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

#include <data/xmipp_program.h>
#include <data/geometry.h>
#include <data/phantom.h>
#include <data/pdb.h>

class Phantom_transform_parameters: public XmippProgram
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

    void defineParams()
    {
        addUsageLine("Apply a geometrical transformation to a phantom description or PDB.");
        addSeeAlsoLine("phantom_create, volume_from_pdb");
        addParamsLine(" -i <file>       : Phantom description file (.descr) or PDB (.pdb)");
        addParamsLine("[-o <file=\"\">] : Phantom description file (.descr) or PDB (.pdb)");
        addParamsLine("                 : For PDB files, you must explicitly give an output file, different from input");
        addParamsLine(" --operation <op>: Operation to perform");
        addParamsLine("      where <op>");
        addParamsLine("            shift <x> <y> <z> : Shift vector");
        addParamsLine("            scale <x> <y> <z> : Scale vector");
        addParamsLine("            rotate_euler <rot> <tilt> <psi> :Rotate with these Euler angles ");
        addParamsLine("            rotate_align_with_z <x> <y> <z> : Align (x,y,z) with Z");
        addParamsLine("            rotate_axis <x> <y> <z> <ang>: Rotate <ang> degrees around (x,y,z)");
        addParamsLine(" [--center_pdb]  : Substract the center of mass from coordinates.");
        addParamsLine("                 :+Only valid for PDB files");
        addExampleLine("xmipp_phantom_transform -i model.pdb -o shifted.pdb --operation shift 1 2 3");
    }

    void readParams()
    {
        fn_in = getParam("-i");
        fn_out = getParam("-o");
        if (fn_out=="")
        	fn_out=fn_in;
        PDBmode=fn_in.contains(".pdb");
        String operation;
        shift.initZeros(3);
        scale.resize(3);
        scale.initConstant(1);
        A3D.initIdentity(4);
        Euler_mode = Align_mode = Axis_mode = false;
        centerPDB = checkParam("--center_pdb");
        operation=getParam("--operation");
        if (operation=="shift")
        {
            double x=getDoubleParam("--operation",1);
            double y=getDoubleParam("--operation",2);
            double z=getDoubleParam("--operation",3);
            shift=vectorR3(x,y,z);
        }
        else if (operation=="scale")
        {
            double x=getDoubleParam("--operation",1);
            double y=getDoubleParam("--operation",2);
            double z=getDoubleParam("--operation",3);
            scale=vectorR3(x,y,z);
        }
        else if (operation=="rotate_euler")
        {
            Euler_mode = true;
            rot  = getDoubleParam("--operation",1);
            tilt = getDoubleParam("--operation",2);
            psi  = getDoubleParam("--operation",3);
            Euler_angles2matrix(rot, tilt, psi, A3D, true);
        }
        else if (operation=="rotate_align_with_z")
        {
            Align_mode = true;
            double x=getDoubleParam("--operation",1);
            double y=getDoubleParam("--operation",2);
            double z=getDoubleParam("--operation",3);
            axis=vectorR3(x,y,z);
            alignWithZ(axis, A3D);
        }
        else if (operation=="rotate_axis")
        {
            Axis_mode = true;
            double x=getDoubleParam("--operation",1);
            double y=getDoubleParam("--operation",2);
            double z=getDoubleParam("--operation",3);
            ang=getDoubleParam("--operation",4);
            axis=vectorR3(x,y,z);
            rotation3DMatrix(ang, axis, A3D);
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
        if (verbose==0)
            return;
        std::cout << "Input file : " << fn_in  << std::endl
        << "Output file: " << fn_out << std::endl;
        if (PDBmode)
            std::cout << "Input file is PDB\n";
        if (centerPDB)
            std::cout << "Centering PDB\n";
        if (Euler_mode)
            std::cout << "Euler angles (rot, tilt, psi): " << rot << " " << tilt
            << " " << psi << std::endl;
        else if (Align_mode)
            std::cout << "Aligning " << axis.transpose() << " with Z\n";
        else if (Axis_mode)
            std::cout << "Rotating " << ang << " degrees around "
            << axis.transpose() << std::endl;
        std::cout << "Transformation matrix\n" << A3D << std::endl;
    }

    void run()
    {
        show();
        if (PDBmode)
            applyGeometryToPDBFile(fn_in, fn_out, A3D, centerPDB);
        else
        {
            Phantom P;
            P.read(fn_in, false); // Read phantom without applying scale
            P.selfApplyGeometry(A3D, IS_NOT_INV);
            P.write(fn_out);
        }
    }
};

int main(int argc, char **argv)
{
    Phantom_transform_parameters prm;
    prm.read(argc, argv);
    return prm.tryRun();
}
