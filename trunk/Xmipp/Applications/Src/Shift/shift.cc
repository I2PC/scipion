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
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippGeometry.hh>

class Shift_parameters: public Prog_parameters {
public:
   matrix1D<double> shift;
   bool             wrap;
   DocFile          DF_shifts;
   int              colX;

   void read(int argc, char **argv) _THROW {
      Prog_parameters::read(argc,argv);
      int i=position_param(argc,argv,"-shift");
      if (i==-1)
         REPORT_ERROR(1,"Shift:: Cannot find -shift");
      else if (i+1>=argc)
         REPORT_ERROR(1,"Shift: Not enough parameters after -shift");
      if (exists(argv[i+1])) {
         DF_shifts.read(argv[i+1]);
	 colX=AtoI(get_param(argc,argv,"-colX","3"));
      }
      else
         shift=get_vector_param(argc,argv,"-shift",-1);
      wrap=!check_param(argc,argv,"-dont_wrap");
   }
   
   void show() {
      Prog_parameters::show();
      if (wrap) cout << "Wrapping image/volume\n";
      else      cout << "Not wrapping image/volume\n";
      if (DF_shifts.name()=="")
         cout << "Shift: " << shift.transpose() << endl;
      else {
         cout << "Shift: " << DF_shifts.name() << endl;
	 cout << "ColX:  " << colX << endl;
      }
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "   -shift [<x>,<y>[,<z>]]   : Shift by (x,y,z)\n"
           << "   -shift <DocFile>         : Shifts are stored in a file\n"
	   << "  [-colX <col=3>]           : Column of the X shift in the DocFile\n"
	   << "                              First column is number 1.\n"
           << "  [-dont_wrap]              : By default, the volume is wrapped\n";
   }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Shift_parameters *eprm=(Shift_parameters *) prm;
   if (eprm->DF_shifts.name()=="")
      img().self_translate(eprm->shift,eprm->wrap);
   else {
      matrix1D<double> shift(2);
      XX(shift)=eprm->DF_shifts(eprm->colX);
      YY(shift)=eprm->DF_shifts(eprm->colX+1);
      img().self_translate(shift,eprm->wrap);
      eprm->DF_shifts.next_data_line();
   }
   return TRUE;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Shift_parameters *eprm=(Shift_parameters *) prm;
   if (eprm->DF_shifts.name()=="")
      vol().self_translate(eprm->shift,eprm->wrap);
   else {
      matrix1D<double> shift(3);
      XX(shift)=eprm->DF_shifts(eprm->colX);
      YY(shift)=eprm->DF_shifts(eprm->colX+1);
      ZZ(shift)=eprm->DF_shifts(eprm->colX+2);
      vol().self_translate(shift,eprm->wrap);
      eprm->DF_shifts.next_data_line();
   }
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
