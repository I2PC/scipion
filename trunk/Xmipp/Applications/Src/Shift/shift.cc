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

class Shift_parameters: public Prog_parameters {
public:
   matrix1D<double> shift;
   bool wrap;

   void read(int argc, char **argv) _THROW {
      Prog_parameters::read(argc,argv);
      shift=get_vector_param(argc,argv,"-shift",-1);
      wrap=!check_param(argc,argv,"-dont_wrap");
   }
   
   void show() {
      Prog_parameters::show();
      if (wrap) cout << "Wrapping image/volume\n";
      else      cout << "Not wrapping image/volume\n";
      cout << "Shift: " << shift.transpose() << endl;
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "   -shift [<x>,<y>[,<z>]]   : Shift by (x,y,z)\n"
           << "  [-dont_wrap]              : By default, the volume is wrapped\n";
   }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Shift_parameters *eprm=(Shift_parameters *) prm;
   img().translate(eprm->shift,eprm->wrap);
   return TRUE;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Shift_parameters *eprm=(Shift_parameters *) prm;
   vol().translate(eprm->shift,eprm->wrap);
   return TRUE;
}

int main (int argc, char **argv) {
   Shift_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Shift {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Shift/Help/shift.html";
      help="Shift volumes and images";
      OPEN MENU menu_shift;
      COMMAND LINES {
	+ usual: xmipp_shift
               #include "prog_line.mnu"
                -shift "["$X","$Y[","$Z]"]"
               [-dont_wrap]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
            $X  {type=float; label="Shift X ";}
            $Y  {type=float; label="Shift Y ";}
        OPT($Z) {type=float; label="Shift Z ";}
        OPT(-dont_wrap) {label="Do not wrap";}
      }
   }

   MENU menu_shift {
      #include "prog_menu.mnu"
      "Shift parameters"
      $X
      $Y
      OPT($Z)
      OPT(-dont_wrap)
   }
*/
