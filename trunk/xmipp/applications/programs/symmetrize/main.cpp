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

#include <reconstruction/symmetrize.h>

int main (int argc, char *argv[]) {
   Symmetrize_Parameters prm;

   try {
      prm.read(argc, argv);
   }
   catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(1);}

   try {
      ROUT_symmetrize(prm);
   } catch (Xmipp_error XE) {cout << XE;}
   exit(0);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Symmetrize {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Symmetrize/Help/symmetrize.html";
      help="Symmetrize a volume using a symmetry description file";
      OPEN MENU menu_symmetrize;
      COMMAND LINES {
	+ usual: symmetrize -i $INPUT_VOL [-o $OUTPUT_VOL] -sym $SYMFILE
      }
      PARAMETER DEFINITIONS {
        $INPUT_VOL {
	   label="Input volume";
           help="Xmipp volume, if an output is not given this is
                 the output, too";
	   type=file existing;
	}
        $OUTPUT_VOL {
	   label="Output volume";
	   help="Xmipp format";
	   type=file;
	}
        $SYMFILE {
           label="Symmetry volume";
	   help="This file has got a complex structure, better see
                 the Web help";
           type=file existing;
        }
      }
   }

   MENU menu_symmetrize {
      "I/O variables"
      $INPUT_VOL
      OPT($OUTPUT_VOL)
      $SYMFILE
   }
*/
