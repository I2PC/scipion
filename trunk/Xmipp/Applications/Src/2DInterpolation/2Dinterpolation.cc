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

#include <XmippData/xmippImages.hh>
#include <XmippData/xmippArgs.hh>

void Usage(char *argv[]);

int main (int argc, char **argv) {
FileName       fn_in;    // input file
double x,y;
ImageXmipp    M;
   // Get command line parameters ------------------------------------------
   try {
      x  = AtoF(get_param(argc, argv, "-u" ));
      y  = AtoF(get_param(argc, argv, "-v" ));
      fn_in = get_param(argc, argv, "-i");		
   } catch (Xmipp_error XE) {Usage(argv);}

   // Main program ---------------------------------------------------------
  
   try {
     M.read(fn_in);
   } catch (Xmipp_error XE) 
   {
   cerr << XE;
   exit(1);
   }
  M().set_Xmipp_origin();
          
  cout << "The value at (" << x << ", ";
  cout << y << ")= ";
  cout << M().interpolated_elem(x,  y) << endl;
  
  exit(0);    
}

void Usage (char *argv[]) {
    cout << "Purpose:\n";
    cout << "    Calculation of the value of the image stored in the file";
    cout << " input_file at the point (u,v) where u and v are any float";
    cout << " number\n\n";
    cout << "Usage:" << argv[0] 
         << " -i input_file"
         << " -u u"
         << " -v v"
         << endl << endl;
    exit(1);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Pixel_value {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/2DInterpolation/Help/2Dinterpolation.html";
      help="Compute the value of an image at a non integer position by
            interpolation";
      OPEN MENU menu_2Dinterpolation;
      COMMAND LINES {
	+ usual: xmipp_2Dinterpolation
                   -i $FILE_IN -u $X -v $Y
      }
      PARAMETER DEFINITIONS {
        $FILE_IN {label="Input image"; type=file existing;}
        $X {type=float; label="X position";}
        $Y {type=float; label="Y position";}
      }
   }

   MENU menu_2Dinterpolation {
      $FILE_IN
      $X
      $Y
   }
*/
