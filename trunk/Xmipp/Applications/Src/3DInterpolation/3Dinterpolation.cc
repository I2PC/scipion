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

#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippArgs.hh>

void Usage(char *argv[]);

int main (int argc, char **argv) {
FileName       fn_in;    // input file
double x,y,z;
VolumeXmipp    V;
   // Get command line parameters ------------------------------------------
   try {
      x  = AtoF(get_param(argc, argv, "-x" ));
      y  = AtoF(get_param(argc, argv, "-y" ));
      z  = AtoF(get_param(argc, argv, "-z" ));
      fn_in = get_param(argc, argv, "-i");		
   } catch (Xmipp_error XE) {Usage(argv);}

   // Main program ---------------------------------------------------------
  
 if (!Is_VolumeXmipp(fn_in))
    {cout << "File: " <<  fn_in<<" is not an Spider volume, bye." 
	  << endl; exit(1);}
  V.read(fn_in);
  V().set_Xmipp_origin();
          
  cout << "The value at (" << x << ", ";
  cout << y << ", ";
  cout << z << ")=";
  cout <<   V().interpolated_elem(x,  y,  z, 0) << endl;
  
  exit(0);   
}

void Usage (char *argv[]) {
    cout << "Purpose:\n";
    cout << "    Calculation of the value of the volume stored in the file";
    cout << " input_file at the point (x,y,z) where x,y and z are any float";
    cout << " number\n\n";
    cout << "Usage:" << argv[0] 
         << " -i input_file"
         << " -x x"
         << " -y y"
         << " -z z"<< endl << endl;
    exit(1);

}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Voxel_value {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/3DInterpolation/Help/3Dinterpolation.html";
      help="Compute the value of a volume at a non integer position by
            interpolation";
      OPEN MENU menu_3Dinterpolation;
      COMMAND LINES {
	+ usual: xmipp_3Dinterpolation
                   -i $FILE_IN -x $X -y $Y -z $Z
      }
      PARAMETER DEFINITIONS {
        $FILE_IN {label="Input volume"; type=file existing;}
        $X {type=float; label="X position";}
        $Y {type=float; label="Y position";}
        $Z {type=float; label="Z position";}
      }
   }

   MENU menu_3Dinterpolation {
      $FILE_IN
      $X
      $Y
      $Z
   }
*/
