/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es
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

/* INCLUDES ---------------------------------------------------------------- */
#include <Reconstruction/Programs/Prog_adjust_surface.hh>

/* MAIN -------------------------------------------------------------------- */
int main (int argc,char *argv[]) {
   Prog_Adjust_Surface_Parameters      prog_prm;

// Check the command line
   try {prog_prm.read(argc,argv);}
   catch (Xmipp_error &XE) {cout << XE; prog_prm.usage(); exit(1);}

try {
// Really adjust surface
  prog_prm.produce_Side_Info();
  ROUT_adjust_surface(prog_prm);
} catch (Xmipp_error XE) {cout << XE;}

}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Adjust_surface {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Adjust_surface/Help/adjust_surface.html";
      help="Find the best position and height of a surface within
            a volume";
      OPEN MENU menu_adjust_surface;
      COMMAND LINES {
         + usual: xmipp_adjust_surface -i $SURFACE [-bottom_surface]
                     [-phantom]
                      -vol $INPUT_VOL
                     [-o $OUTPUT_SURF]
                     [-ztop $ZTOP0 $ZTOPF]
                     [-zbottom $ZBOTTOM0 $ZBOTTOMF]
                     $CORRELATION_METHOD
                     [-manual_order]
      }
      PARAMETER DEFINITIONS {
         $INPUT_VOL {
            label="Input density volume";
            help="On top of which the surface will be fitted";
            type=file existing;
         }
         $SURFACE {
            label="Surface image";
            help="Image with the surface mask profile";
            type=file existing;
         }
         OPT(-bottom_surface) {label="Suface is a bottom surface";}
         OPT(-phantom) {label="Surface comes from a phantom";}
         $OUTPUT_SURF {
            label="Output surface";
            help="The input surface is rewritten if this file is not
                  provided";
            type=file;
         }
         OPT(-ztop) {
            label="Highest position range";
            help="The surface is forced to fit with its highest position
                  within this range";
         }
            $ZTOP0 {
               label="Initial";
               help="In logical coordinates";
               type=float;
            }
            $ZTOPF {COPY $ZTOP0; label="Final";}
         OPT(-zbottom) {
            label="Lowest position range";
            help="The surface is forced to fit with its lowest position
                  within this range";
         }
            $ZBOTTOM0 {COPY $ZTOP0;}
            $ZBOTTOMF {COPY $ZTOPF;}
         $CORRELATION_METHOD {
            label="Correlation method";
            type=Exclusion {
               "3D" {}
               "2D" {-corr_2D}
               "Gradient" {-corr_grad}
            };
         }
         OPT(-manual_order) {label="Supply positions manually";}
      }
   }
   
   MENU menu_adjust_surface {
      "I/O parameters"
      $INPUT_VOL
      $SURFACE
      OPT(-bottom_surface)
      OPT(-phantom)
      OPT(-o)
      "Adjusting options"
      OPT(-ztop)
      OPT(-zbottom)
      $CORRELATION_METHOD
      "Debugging options"
      OPT(-manual_order)
   }
*/
