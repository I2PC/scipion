/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (1998)
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

#include <reconstruction/project.h>

int main(int argc, char *argv[])
{
    Prog_Project_Parameters      prog_prm;
    Projection                   proj;
    SelFile                      SF;

// Check the command line
    try
    {
        prog_prm.read(argc, argv);
    }
    catch (Xmipp_error &XE)
    {
        cout << XE;
        prog_prm.usage();
        exit(1);
    }

    try
    {
// Really project
        ROUT_project(prog_prm, proj, SF);
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }

}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Project {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Project/Help/project.html";
      help="Project a volume along the desired directions";
      OPEN MENU menu_project;
      COMMAND LINES {
 + usual: project -i $FILE_IN [-o $SEL_OUT] [-show_angles]
                         [-crystal $CRYSTAL_PRM]
      }
      PARAMETER DEFINITIONS {
        $FILE_IN {
    label="Projection Parameters file";
    help="This file has got a complex structure, better see
                 the Web help";
    type=file existing;
 }
        $SEL_OUT {
    label="Output sel file";
    help="If not given, there is no output selfile";
    type=file;
 }
        OPT(-show_angles) {
           label="Show angles in screen";
           help="The angle order is rot, tilt, psi";
        }
        $CRYSTAL_PRM {
    label="Crystal parameters";
           help="This file has got a complex structure, better see
                 the Web help";
    type=file;
 }
      }
   }

   MENU menu_project {
      "I/O parameters"
      $FILE_IN
      OPT(-o)
      "Options"
      OPT(-show_angles)
      OPT(-crystal)
   }
*/
