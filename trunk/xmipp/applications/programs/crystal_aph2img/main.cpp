/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu )
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

#include <data/args.h>
#include <reconstruction/crystal_aph2img.h>

void Usage(char *argv[]);

int main(int argc, char **argv)
{
    FileName                  fn_prm;
    Spot2RealSpace2D_Parameters prm;
    Projection                  prj;

    // Get command line parameters ------------------------------------------
    try
    {
        fn_prm = getParameter(argc, argv, "-i");
        prm.read_from_file(fn_prm);
    }
    catch (Xmipp_error XE)
    {
        Usage(argv);
    }

    // Main program ---------------------------------------------------------
    try
    {
        ROUT_Spots2RealSpace(prm, prj);
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }
    exit(0);
}

void Usage(char *argv[])
{
    cout << "Purpose:\n"
    << "    Perform a DISCRETE but not FAST Fourier Transform\n"
    << "    The input is an aph (MRC) file\n"

    << "Usage:" << argv[0] << " -i filename" << endl << endl
    << "\t-i           : Parameters file" << endl
    ;
    exit(1);

}

/* Menu -------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Spots2RealSpace2D {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Spots2RealSpace2D/Help/spots2realspace2D.html";
      help="Generate an image suitable for ART+crystals starting from
            a set of amplitudes and phases provided by MRC";
      OPEN MENU Spots2RealSpace2D;
      COMMAND LINES {
         + usual: spots2realspace2D $FILE_IN
      }
      PARAMETER DEFINITIONS {
         $FILE_IN {
            label="Input File";
            help="File with a complex structure, better see Web help";
            type=FILE EXISTING;
         }
      }
   }
   MENU Spots2RealSpace2D {
      "I/O Parameters"
      $FILE_IN
   }
*/
