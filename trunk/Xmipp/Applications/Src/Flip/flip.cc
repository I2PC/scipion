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

class Flip_parameters: public Prog_parameters {
public:
   bool flipX;
   bool flipY;
   bool flipZ;

   void read(int argc, char **argv) {
      Prog_parameters::read(argc,argv);
      flipX=check_param(argc,argv,"-flipX");
      flipY=check_param(argc,argv,"-flipY");
      flipZ=check_param(argc,argv,"-flipZ");
   }
   
   void show() {
      Prog_parameters::show();
      cout << "FlipX = " << flipX << endl
           << "FlipY = " << flipY << endl
           << "FlipZ = " << flipZ << endl;
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "  [-flipX]                  : Flip along X\n"
           << "  [-flipY]                  : Flip along Y\n"
           << "  [-flipZ]                  : Flip along Z\n";
   }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Flip_parameters *eprm=(Flip_parameters *) prm;
   if (eprm->flipX) img().self_reverseX();
   if (eprm->flipY) img().self_reverseY();
   return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Flip_parameters *eprm=(Flip_parameters *) prm;
   if (eprm->flipX) vol().self_reverseX();
   if (eprm->flipY) vol().self_reverseY();
   if (eprm->flipZ) vol().self_reverseZ();
   return true;
}

int main (int argc, char **argv) {
   Flip_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Flip {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Flip/Help/flip.html";
      help="Flip (mirror) volumes and images";
      OPEN MENU menu_flip;
      COMMAND LINES {
	+ usual: xmipp_flip
               #include "prog_line.mnu"
               [-flipX]
               [-flipY]
               [-flipZ]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        OPT(-flipX) {label="Flip X";}
        OPT(-flipY) {label="Flip Y";}
        OPT(-flipZ) {label="Flip Z";}
      }
   }

   MENU menu_flip {
      #include "prog_menu.mnu"
      "Flipping parameters"
      OPT(-flipX)
      OPT(-flipY)
      OPT(-flipZ)
   }
*/
