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
#include <Reconstruction/Programs/Prog_draw_surface.hh>

/* MAIN -------------------------------------------------------------------- */
int main (int argc,char *argv[]) {
   Prog_Draw_Surface_Parameters prm;

// Check the command line
   try {prm.read(argc,argv);}
   catch (Xmipp_error &XE) {cout << XE; prm.usage(); exit(1);}

// Execute program
   try {
      ROUT_draw_surface(prm);
   } catch (Xmipp_error XE) {cout << XE;}
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Draw_surface {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Draw_surface/Help/draw_surface.html";
      help="Draw a surface image onto a density volume";
      OPEN MENU menu_draw_surface;
      COMMAND LINES {
         + usual: xmipp_draw_surface -i $INPUT_VOL -s $SURFACE
                     [-o $OUTPUT_VOL] [-ztop $ZTOP -zbottom $ZBOTTOM]
                     [-color $COLOR]
      }
      PARAMETER DEFINITIONS {
         $INPUT_VOL {
            label="Input density volume";
            help="On top of which the surface will be drawn";
            type=file existing;
         }
         $SURFACE {
            label="Surface image";
            help="Image with the surface mask profile";
            type=file existing;
         }
         $OUTPUT_VOL {
            label="Output volume";
            help="The input volume is rewritten if this file is not
                  provided";
            type=file;
         }
         OPT(-ztop) {
            label="Supply surface position";
            help="The surface is forced to fit within this position";
         }
         $ZTOP {
            label="Highest position";
            help="In logical coordinates";
            type=float;
         }
         $ZBOTTOM {COPY $ZTOP; label="Lowest position";}
         $COLOR {label="Density for surface"; type=float;}
      }
   }
   
   MENU menu_draw_surface {
      "I/O parameters"
      $INPUT_VOL
      OPT(-o)
      $SURFACE
      "Drawing options"
      OPT(-ztop)
      OPT(-color)
   
   }
*/
