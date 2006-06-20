/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
#include <XmippData/xmippFFT.hh>

bool process_img(ImageXmipp &img, FourierImageXmipp &IMG,
   const Prog_parameters *prm) {
   FourierTransform(img(), IMG());
   //copy angle information
   double phi,theta,psi;
   img.get_eulerAngles(phi,theta,psi);
   IMG.set_eulerAngles(phi,theta,psi);
   return true;
}

bool process_vol(VolumeXmipp &vol, FourierVolumeXmipp &VOL,
   const Prog_parameters *prm) {
   FourierTransform(vol(), VOL());
   return true;
}

int main (int argc, char **argv) {
   Prog_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol,
      IMAGE2FOURIER);
}

/* Menus =================================================================== */
/*Colimate:
   PROGRAM FFT {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/FFT/Help/FFT.html";
      help="Fourier Transform of volumes and images";
      OPEN MENU menu_FFT;
      COMMAND LINES {
	+ usual: xmipp_FFT
               #include "prog_line.mnu"
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
      }
   }

   MENU menu_FFT {
      #include "prog_menu.mnu"
   }
*/
