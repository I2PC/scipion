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

#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippGeometry.hh>

class Rotate_parameters: public Prog_parameters {
public:
   bool Euler_mode;
      double rot, tilt, psi;
   bool Align_mode;
   bool Axis_mode;
      matrix1D<double> axis;
      double ang;
   bool wrap;
   
   matrix2D<double> A3D, A2D;

   void read(int argc, char **argv) _THROW {
      Prog_parameters::read(argc,argv);
      Euler_mode=Align_mode=Axis_mode=FALSE;
      if (check_param(argc,argv,"-euler")) {
         Euler_mode=TRUE;
         int i=position_param(argc,argv,"-euler");
         if (i+3>=argc)
            REPORT_ERROR(1,"Not enough parameters after -euler");
         rot  = AtoF(argv[i+1]);
         tilt = AtoF(argv[i+2]);
         psi  = AtoF(argv[i+3]);
         A3D=Euler_rot3D_matrix(rot,tilt,psi);
      } else if (check_param(argc,argv,"-align_with_Z")) {
         Align_mode=TRUE;
         axis=get_vector_param(argc,argv,"-align_with_Z",3);
         A3D=align_with_Z(axis);
      } else {
         Axis_mode=TRUE;
         if (check_param(argc, argv, "-axis"))
            axis=get_vector_param(argc, argv, "-axis", 3);
         else 
            axis=vector_R3(0.,0.,1.);
         ang=AtoF(get_param(argc, argv, "-ang"));
         A3D=rot3D_matrix(ang,axis);
         A2D=A3D; A2D.window(0,0,2,2);
      }
      wrap=!check_param(argc,argv,"-dont_wrap");
   }
   
   void show() {
      Prog_parameters::show();
      if (Euler_mode)
         cout << "Euler angles (rot, tilt, psi): " << rot << " " << tilt
              << " " << psi << endl;
      else if (Align_mode)
         cout << "Aligning " << axis.transpose() << " with Z\n";
      else if (Axis_mode)
         cout << "Rotating " << ang << " degrees around " << axis.transpose()
              << endl;
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "  [-euler <rot> <tilt> <psi>        : Rotate with these Euler angles\n"
           << "  [-align_with_Z [<x>,<y>,<z>]]     : Align (x,y,z) with Z\n"
           << "                                      Notice that brackets for the\n"
           << "                                      vector must be written and do not\n"
           << "                                      represent optional parameters\n"
           << "  [[-axis [<x>,<y>,<z>]] -ang <ang>]: Rotate <ang> degrees around (x,y,z),\n"
           << "                                      by default (0,0,1)\n"
           << "  [-dont_wrap]                      : By default, the volume is wrapped\n";
   }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Rotate_parameters *eprm=(Rotate_parameters *) prm;
   Image img_out;
   if (XSIZE(eprm->A2D)!=0) {
      apply_geom(img_out(),eprm->A2D,img(),IS_NOT_INV,eprm->wrap);
      img()=img_out();
   }
   return TRUE;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Rotate_parameters *eprm=(Rotate_parameters *) prm;
   Volume vol_out;
   apply_geom(vol_out(),eprm->A3D,vol(),IS_NOT_INV,eprm->wrap);
   vol()=vol_out();
   return TRUE;
}

int main (int argc, char **argv) {
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
               [-align_with_Z] [-axis]["["$X","$Y","$Z"]"]
               [-ang $ANG]
               [-dont_wrap]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        $ROTATION_METHOD {
           label="Rotation action";
           type=list {
              "Euler rotation" {OPT(-euler)=1; OPT(-align_with_Z)=0;
                                OPT(-axis)=0; OPT($X)=0;}
              "Align with Z"   {OPT(-euler)=0; OPT(-align_with_Z)=1;
                                OPT(-axis)=0; OPT($X)=1;}
              "Around an axis" {OPT(-euler)=0; OPT(-align_with_Z)=0;
                                OPT(-axis)=1; OPT($X)=1;}
           };
        }
        OPT(-euler) {label="Euler rotation";}
           $ROT  {type=float; label="Rotational angle";}
           $TILT {type=float; label="Tilting angle";}
           $PSI  {type=float; label="In-plane rotation";}
        OPT(-align_with_Z) {label="Align with Z";}
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
      OPT(-align_with_Z)
      OPT(-axis)
      OPT($X)
      OPT($ANG)
      OPT(-dont_wrap)
   }
*/
