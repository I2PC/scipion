/***************************************************************************
 *
 * Authors:     Roberto Marabini 
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

#include <XmippData/xmippArgs.hh>
#include <Reconstruction/Programs/Prog_euler.hh>

void Usage(char *argv[]);

int main (int argc, char **argv) {
double rot, tilt, psi ;
double d00, d01, d10, d11; 
   // Get command line parameters ------------------------------------------
   try {
      if (argc<2) REPORT_ERROR(1,"");
      rot  = AtoF(get_param(argc, argv, "-rot" , "30.0"));
      tilt = AtoF(get_param(argc, argv, "-tilt", "20.0"));
      psi  = AtoF(get_param(argc, argv, "-psi" , "0.0"));
      d00  = AtoF(get_param(argc, argv, "-d00" ,
        				  "1.0087"));
      d01  = AtoF(get_param(argc, argv, "-d01" , "0.0"));
      d10  = AtoF(get_param(argc, argv, "-d10" , "0.0"));
      d11  = AtoF(get_param(argc, argv, "-d11" ,
                                                "1.0087"));
//      if(check_param(argc, argv, "-help")) Usage(argv);  						
   } catch (Xmipp_error XE) {Usage(argv);}

   // Main program ---------------------------------------------------------
   try {
       
       ROUT_EULER(rot,tilt,psi,d00,d01,d10,d11);
   } catch (Xmipp_error XE) {cout<<XE;}
   exit(0);
}

void Usage (char *argv[]) {
    cout << "Purpose:\n";
    cout << "    Calculate new  projection angle";
    cout << " after deformation  transform\n";
    cout << "Usage:" << argv[0] 
         <<" -rot rot" 
         << " -tilt tilt"
         << " -psi psi"
         << " -d00 d00"
         << " -d01 d01"
         << " -d10 d10"
         << " -d11 d11"<< endl << endl;
    exit(1);

}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Euler {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Euler/Help/euler.html";
      help="Compute new rot, tilt and psi after a deformation of the volume";
      OPEN MENU menu_euler;
      COMMAND LINES {
         +usual: xmipp_euler -rot $ROT -tilt $TILT -psi $PSI
                    -d00 $D00 -d01 $D01 -d10 $D10 -d11 $D11 
      }
      PARAMETER DEFINITIONS {
         $ROT  {label="Rotational angle"; type=float;}
         $TILT {label="Tilting angle"; type=float;}
         $PSI  {label="In-plane rotational angle"; type=float;}
         $D00 {label="D00"; type=float;}
         $D01 {label="D01"; type=float;}
         $D10 {label="D10"; type=float;}
         $D11 {label="D11"; type=float;}
      }
   }

   MENU menu_euler {
      "Projection angles"
      $ROT
      $TILT
      $PSI
      "Deformation matrix"
      {$D00 $D01}
      {$D10 $D11}
   }
*/
