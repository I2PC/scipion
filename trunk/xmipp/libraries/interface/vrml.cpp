/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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

#include "vrml.h"

VrmlFile::VrmlFile(FileName vrmlfilename)
{
    VrmlFile::fh_out.open(vrmlfilename.c_str(), ios::out); //trunc file
    VrmlFile::need_footer = 0;
    VrmlFile::fh_out << "#VRML V2.0 utf8" << endl;
    VrmlFile::fh_out << "# (XMIPP libraries) by R. Marabini" << endl ;
}/* VrmlFile */

void VrmlFile::Sphere(matrix1D<double> XYZ,
                      matrix1D<double> RGB,
                      double radius)
{
    double X, Y, Z;
    double R, G, B;

    X = XX(XYZ);
    Y = YY(XYZ);
    Z = ZZ(XYZ);
    R = XX(RGB);
    G = YY(RGB);
    B = ZZ(RGB);

    if ((VrmlFile::need_footer) == 1)
        VrmlFile::fh_out << "    ]\n}\n";

    VrmlFile::fh_out << "Transform {" << endl
    << "    children [" << endl
    << "       Transform { translation " << X << " " << Y << " " << Z << endl
    << "        children[" << endl
    << "        DEF BallRow Group {" << endl
    << "            children [" << endl
    << "                Shape {" << endl
    << "                   appearance Appearance {" << endl
    << "                        material Material { emissiveColor " << R << " " << G << " " << B << endl
    << "                                            diffuseColor  " << R << " " << G << " " << B << endl
    << "                                            specularColor 1 1 1"               << endl
    << "                                            shininess     0.1" <<          " }" << endl
    << "                    }" << endl
    << "                    geometry Sphere { radius " << radius << " }" << endl
    << "                }" << endl
    << "            ]" << endl
    << "        }" << endl
    << "        ]" << endl
    << "        }" << endl;
    VrmlFile::need_footer = 1;

}// sphere end

void VrmlFile::Add_sphere(const matrix1D<double> XYZ)
{
    double X, Y, Z;
    X = XX(XYZ);
    Y = YY(XYZ);
    Z = ZZ(XYZ);

    VrmlFile::fh_out << "       Transform { translation  " << X << " " << Y << " " << Z << " children [USE BallRow] }     " << endl;
}//add sphere


void VrmlFile::Trans_Sphere(const matrix1D<double> XYZ,
                            const matrix1D<double> RGB,
                            const double radius)
{
    double X, Y, Z;
    double R, G, B;

    X = XX(XYZ);
    Y = YY(XYZ);
    Z = ZZ(XYZ);
    R = XX(RGB);
    G = YY(RGB);
    B = ZZ(RGB);

    if (VrmlFile::need_footer == 1)
        VrmlFile::fh_out << "    ]\n}\n";
    VrmlFile::fh_out
    << "# Trans_Sphere " << endl
    << "Transform {" << endl
    << "    children [" << endl
    << "       Transform { translation " << X << " " << Y << " " << Z << endl
    << "        children[" << endl
    << "        DEF BallRow Group {" << endl
    << "            children [" << endl
    << "                Shape {" << endl
    << "                   appearance Appearance {" << endl
    << "                       material Material { emissiveColor " << R << " " << G << " " << B << endl
    << "                                           diffuseColor  " << R << " " << G << " " << B << endl
    << "                                           transparency  0.5}" << endl
    << "                }" << endl
    << "                    geometry Sphere { radius " << radius << " }" << endl
    << "            }" << endl
    << "        ]" << endl
    << "    }" << endl
    << "    ]" << endl
    << "    }" << endl;
    VrmlFile::need_footer = 1;
//      << "    ]" << endl
//      << "}" << endl;
}//transparent sphere end

/*-----------------------------------------------
Draw coordinate system, x axis is red, y axis is green, z axis is blue.
 Parameters are the scale of the coordinate system and the relative size
 between the cone and the cylinder. Defaults for both parameters are 1.
-------------------------------------------------*/
void VrmlFile::Axis(double scale, double Cone_scale)
{
    if (VrmlFile::need_footer == 1)
        VrmlFile::fh_out << "    ]\n}\n";
//#ifdef NEFERER
    VrmlFile::fh_out
    << "#X\n"
    << "Transform {\n"
    << " rotation 0 0 1 -1.57 \n"
    << " scale " << scale << " " << scale << " " << scale << "\n"
    << " children [\n"
    << "  Transform{\n"
    << "   translation 0.0  0.5 0.0\n"
    << "   children [ \n"
    << "    DEF line Shape{\n"
    << "        appearance Appearance {\n"
    << "         material Material { \n"
    << "       diffuseColor   1.0 0.0 0.0\n"
    << "       emissiveColor  1.0 0.0 0.0\n"
    << "      }\n"
    << "     }\n"
    << "     geometry Cylinder {\n"
    << "      height 1.00\n"
    << "      radius 0.01\n"
    << "     }\n"
    << "    }\n"
    << "   ]\n"
    << "  } \n"
    << "  Transform{\n"
    << "   translation 0.0 1.0 0.0\n"
    << "   children [\n"
    << "    DEF head Shape{\n"
    << "        appearance Appearance {\n"
    << "         material Material { \n"
    << "       diffuseColor  0.0 0.0 1.0\n"
    << "       emissiveColor 1.0 0.0 0.0\n"
    << "      }\n"
    << "     }\n"
    << "     geometry Cone {\n"
    << "      height       " << 0.50 * Cone_scale << "\n"
    << "      bottomRadius " << 0.05 * Cone_scale << "\n"
    << "     }\n"
    << "    }\n"
    << "   ]\n"
    << "  } \n"
    << " ]\n"
    << "}\n"
    << "\n"
    << "#Y\n"
    << "Transform {\n"
    << " rotation 1 0 0 0. \n"
    << " scale " << scale << " " << scale << " " << scale << "\n"
    << " children [\n"
    << "  Transform{\n"
    << "   translation 0.0  0.5 0.0\n"
    << "   children [ \n"
    << "    DEF line Shape{\n"
    << "        appearance Appearance {\n"
    << "         material Material { \n"
    << "       diffuseColor  .0 1.0 0.0\n"
    << "       emissiveColor .0 1.0 0.0\n"
    << "      }\n"
    << "     }\n"
    << "     geometry Cylinder {\n"
    << "      height 1.00\n"
    << "      radius 0.01\n"
    << "     }\n"
    << "    }\n"
    << "   ]\n"
    << "  } \n"
    << "  Transform{\n"
    << "   translation 0.0 1.0 0.0\n"
    << "   children [\n"
    << "    DEF head Shape{\n"
    << "        appearance Appearance {\n"
    << "         material Material { \n"
    << "       diffuseColor  0.0 0.0 1.0\n"
    << "       emissiveColor 1.0 0.0 0.0\n"
    << "      }\n"
    << "     }\n"
    << "     geometry Cone {\n"
    << "      height       " << 0.50 * Cone_scale << "\n"
    << "      bottomRadius " << 0.05 * Cone_scale << "\n"
    << "     }\n"
    << "    }\n"
    << "   ]\n"
    << "  } \n"
    << " ]\n"
    << "}\n"
    << "\n"
    << "\n"
    << "#Z\n"
    << "Transform {\n"
    << "rotation 1 0 0 1.57 \n"
    << " scale " << scale << " " << scale << " " << scale << "\n"
    << " children [\n"
    << "  Transform{\n"
    << "   translation 0.0  0.5 0.0\n"
    << "   children [ \n"
    << "    DEF line Shape{\n"
    << "        appearance Appearance {\n"
    << "         material Material { \n"
    << "       diffuseColor  0.0 0.0 1.0\n"
    << "       emissiveColor 0.0 0.0 1.0\n"
    << "      }\n"
    << "     }\n"
    << "     geometry Cylinder {\n"
    << "      height 1.00\n"
    << "      radius 0.01\n"
    << "     }\n"
    << "    }\n"
    << "   ]\n"
    << "  }\n"
    << "  Transform{\n"
    << "   translation 0.0 1.0 0.0\n"
    << "   children [\n"
    << "    DEF head Shape{\n"
    << "        appearance Appearance {\n"
    << "         material Material { \n"
    << "       diffuseColor  0.0 0.0 1.0\n"
    << "       emissiveColor 1.0 0.0 0.0\n"
    << "      }\n"
    << "     }\n"
    << "     geometry Cone {\n"
    << "      height       " << 0.50 * Cone_scale << "\n"
    << "      bottomRadius " << 0.05 * Cone_scale << "\n"
    << "     }\n"
    << "    }\n"
    << "   ]\n"
    << "  }\n";
    VrmlFile::need_footer = 1;
// ]
//}
//#endif
}//end Axis
