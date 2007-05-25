/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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
#include <data/matrix3d.h>
#include <data/matrix2d.h>
#include <data/geometry.h>
#include <data/volume.h>

void Usage(char *argv[]);

int main(int argc, char **argv)
{
    FileName       fn_in;     // input file
    FileName       fn_out;    // output file
    double d00, d01, d10, d11;

    VolumeXmipp    V_in, V_out;
    matrix2D<double> D(4, 4);

    D.initZeros();
    D(0, 0) = D(1, 1) = D(2, 2) = D(3, 3) = 1.0f;

    // Get command line parameters ------------------------------------------
    try
    {
        D(0, 0) = AtoF(get_param(argc, argv, "-d00", "1.0"));
        D(0, 1) = AtoF(get_param(argc, argv, "-d01", "0.0"));
        D(1, 0) = AtoF(get_param(argc, argv, "-d10", "0.0"));
        D(1, 1) = AtoF(get_param(argc, argv, "-d11", "1.0"));
        fn_in  = get_param(argc, argv, "-i");
        fn_out = get_param(argc, argv, "-o");
    }
    catch (Xmipp_error XE)
    {
        Usage(argv);
    }

    cout << "Matrix D: \n" << D;
    // Main program ---------------------------------------------------------


    if (!Is_VolumeXmipp(fn_in))
    {
        cout << "File: " <<  fn_in << " is not an Spider volume, bye."
        << endl;
        exit(1);
    }
    V_in.read(fn_in);

    apply_geom(V_out.img, D, V_in.img, IS_INV, DONT_WRAP);
    V_out.write(fn_out);
    exit(0);
}

void Usage(char *argv[])
{
    cout << "Purpose:\n";
    cout << "    Applies a lineal transformation defined by a 2x2 matrix";
    cout << " to all the planes perpendicular to Z that form a volume\n\n";
    cout << "Usage:" << argv[0] ;
    cout << argv[0] << " -i input file -o outout_file"
    << " -d00 d00 -d01 d01 -d10 d10 -d11 d11\n"
    << "where d00, d01, d10 and d11"
    << "are the element of the 2x2 matrix " << endl << endl;
    exit(1);

}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Deform {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Deform/Help/deform.html";
      help="Applies a linear transformation to all Z planes in a volume";
      OPEN MENU menu_deform;
      COMMAND LINES {
         +usual: xmipp_deform -i $FILE_IN -o $FILE_OUT
                    -d00 $D00 -d01 $D01 -d10 $D10 -d11 $D11
      }
      PARAMETER DEFINITIONS {
         $FILE_IN  {label="Input volume"; type=file existing;}
         $FILE_OUT {label="Output volume"; type=file; by default=$FILE_IN;}
         $D00 {label="D00"; type=float;}
         $D01 {label="D01"; type=float;}
         $D10 {label="D10"; type=float;}
         $D11 {label="D11"; type=float;}
      }
   }

   MENU menu_deform {
      "I/O Parameters"
      $FILE_IN
      $FILE_OUT
      "Deformation matrix"
      {$D00 $D01}
      {$D10 $D11}
   }
*/
