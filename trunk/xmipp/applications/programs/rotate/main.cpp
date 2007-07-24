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
#include <data/gridding.h>

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
    bool gridding;

    Matrix2D<double> A3D, A2D;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
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
        else
        {
            Axis_mode = true;
            if (checkParameter(argc, argv, "-axis"))
                axis = getVectorParameter(argc, argv, "-axis", 3);
            else
                axis = vectorR3(0., 0., 1.);
            ang = textToFloat(getParameter(argc, argv, "-ang"));
            A3D = rotation3DMatrix(ang, axis);
            A2D = A3D;
            A2D.window(0, 0, 2, 2);
        }
        wrap = !checkParameter(argc, argv, "-dont_wrap");
        gridding = checkParameter(argc, argv, "-gridding");

    }

    void show()
    {
        Prog_parameters::show();
        if (Euler_mode)
            cout << "Euler angles (rot, tilt, psi): " << rot << " " << tilt
            << " " << psi << endl;
        else if (Align_mode)
            cout << "Aligning " << axis.transpose() << " with Z\n";
        else if (Axis_mode)
            cout << "Rotating " << ang << " degrees around " << axis.transpose()
		 << endl;
       if (gridding)
	   cout << "Use reverse gridding for interpolation."<<endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        cerr << "  [-euler <rot> <tilt> <psi>        : Rotate with these Euler angles\n"
        << "  [-alignWithZ \"[<x>,<y>,<z>]\"]   : Align (x,y,z) with Z\n"
        << "                                      Notice that brackets for the\n"
        << "                                      vector must be written and do not\n"
        << "                                      represent optional parameters\n"
        << "  [-axis \"[<x>,<y>,<z>]\" -ang <ang>]: Rotate <ang> degrees around (x,y,z),\n"
        << "                                      by default (0,0,1)\n"
        << "  [-dont_wrap]                      : By default, the image/volume is wrapped\n"
	<< "  [-gridding]                       : Use reverse gridding for interpolation\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Rotate_parameters *eprm = (Rotate_parameters *) prm;
    Image img_out;
    if (XSIZE(eprm->A2D) != 0)
    {
	if (eprm->gridding)
	{
	    KaiserBessel kb;
	    produceGriddingImage(img(),img_out(),kb);
	    applyGeometryGridding(img(), eprm->A2D, img_out(), kb, IS_NOT_INV, eprm->wrap);
	}
	else
	{
	    applyGeometryBSpline(img_out(), eprm->A2D, img(), 3, IS_NOT_INV, eprm->wrap);
	    img() = img_out();
	}
    }
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Rotate_parameters *eprm = (Rotate_parameters *) prm;
    Volume vol_out;
    if (eprm->gridding)
    {
	KaiserBessel kb;
	produceGriddingVolume(vol(),vol_out(),kb);
	applyGeometryGridding(vol(), eprm->A3D, vol_out(), kb, IS_NOT_INV, eprm->wrap);
    }
    else
    {
	applyGeometryBSpline(vol_out(), eprm->A3D, vol(), 3, IS_NOT_INV, eprm->wrap);
	vol() = vol_out();
    }
    return true;
}

int main(int argc, char **argv)
{
    Rotate_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Rotate {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Rotate/Help/rotate.html";
      help="Rotate volumes and images";
      OPEN MENU menu_rotate;
      COMMAND LINES {
 + usual: xmipp_rotate
               #include "prog_line.mnu"
               $ROTATION_METHOD
               [-euler $ROT $TILT $PSI]
               [-alignWithZ] [-axis]["["$X","$Y","$Z"]"]
               [-ang $ANG]
               [-dont_wrap]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        $ROTATION_METHOD {
           label="Rotation action";
           type=list {
              "Euler rotation" {OPT(-euler)=1; OPT(-alignWithZ)=0;
                                OPT(-axis)=0; OPT($X)=0;}
              "Align with Z"   {OPT(-euler)=0; OPT(-alignWithZ)=1;
                                OPT(-axis)=0; OPT($X)=1;}
              "Around an axis" {OPT(-euler)=0; OPT(-alignWithZ)=0;
                                OPT(-axis)=1; OPT($X)=1;}
           };
        }
        OPT(-euler) {label="Euler rotation";}
           $ROT  {type=float; label="Rotational angle";}
           $TILT {type=float; label="Tilting angle";}
           $PSI  {type=float; label="In-plane rotation";}
        OPT(-alignWithZ) {label="Align with Z";}
        OPT($X) {label="Axis";}
           $Z  {type=float; label="Z ";}
           $Y  {type=float; label="Y ";}
           $X  {type=float; label="X ";}
        OPT(-axis) {label="Rotate around an axis";}
        $ANG {type=float; label="Angle";}
        OPT(-dont_wrap) {label="Do not wrap";}
      }
   }

   MENU menu_rotate {
      #include "prog_menu.mnu"
      "Rotation parameters"
      $ROTATION_METHOD
      OPT(-euler)
      OPT(-alignWithZ)
      OPT(-axis)
      OPT($X)
      OPT($ANG)
      OPT(-dont_wrap)
   }
*/
