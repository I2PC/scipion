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
#include <Reconstruction/Programs/Prog_surface.hh>

/* MAIN -------------------------------------------------------------------- */
int main (int argc,char *argv[]) {
   Prog_Surface_Parameters      prog_prm;

// Check the command line
   try {prog_prm.read(argc,argv);}
   catch (Xmipp_error &XE) {cout << XE; prog_prm.usage(); exit(1);}

try {
// Really create surface
  prog_prm.produce_Side_Info();
  ROUT_surface(prog_prm);
} catch (Xmipp_error XE) {cout << XE;}

}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Surface {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Surface/Help/surface.html";
      help="Generate surfaces for phantoms or combine surfaces into a 3D mask";
      OPEN MENU menu_surface;
      COMMAND LINES {
	+ create_mask: surface -i $FILE_IN -o $FILE_OUT
           [-bottom $BOTTOM_SURF] [-zbottom $ZBOTTOM]
           [-top    $TOP_SURF]    [-ztop    $ZTOP]
        + create_top:    surface -i $FILE_IN OPT(-top)     OPT(-ztop)
        + create_bottom: surface -i $FILE_IN OPT(-bottom)  OPT(-zbottom)
        + combine_surfaces: surface -o $FILE_OUT OPT(-top) OPT(-bottom)
      }
      PARAMETER DEFINITIONS {
        $FILE_IN {
	   label="Input Phantom";
	   help="Phantom description";
	   type=file existing;
	}
        $FILE_OUT {
	   label="Output Mask";
	   type=file;
	}
        $BOTTOM_SURF {
           label="Bottom surface file";
           help="You have to activate it in the create_bottom mode";
           type=file;
        }
        $ZBOTTOM {
           label="Maximum bottom plane";
           help="Limit bottom ray.
                  The bottom surface probe cannot go beyond this point";
           type=FLOAT;
        }
        $TOP_SURF {
           label="top surface file";
           help="You have to activate it in the create_top mode";
           type=file;
        }
        $ZTOP {
           label="Maximum top plane";
           help="Limit top ray.
                  The top surface probe cannot go beyond this point";
           type=FLOAT;
        }
      }
   }

   MENU menu_surface {
      "I/O Parameters"
      $FILE_IN
      $FILE_OUT
      "Surface restrictions"
      OPT(-bottom)
      OPT(-zbottom)
      OPT(-top)
      OPT(-ztop)
   }
*/
