/***************************************************************************
 *
 * Authors:    Carlos Oscar Sánchez Sorzano      coss@cnb.uam.es
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
//Sun Nov 14 22:07:48 EST 1999: added discrete option (R. Marabini)

#include <Reconstruction/Programs/Prog_random_phantom.hh>

/* Program ================================================================= */
int main (int argc,char **argv ) {
// Parameters
   Prog_Random_Phantom_Parameters prog_prm;
   Phantom                        Realization;

// Read parameters from command line
   try {prog_prm.read(argc,argv);}
   catch (Xmipp_error XE) {cout << XE; prog_prm.usage(); exit(1);}

// Generate Realization and write to disk
   try{ROUT_random_phantom(prog_prm, Realization);}
   catch (Xmipp_error XE) {cout << XE;}
   exit(0);
}
/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Random_phantom {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Random_phantom/Help/random_phantom.html";
      help="Produce a realizaion of a random phantom family";
      OPEN MENU menu_random_phantom;
      COMMAND LINES {
	+ single: random_phantom -i $PHANTOM_FAMILY -o $PHANTOM
                  [-min_volume $MINVOLUME] [-discrete]
      }
      PARAMETER DEFINITIONS {
        $PHANTOM_FAMILY {
	   label="Phantom Family Description";
	   help="Input: This file has got a complex structure, better see 
                 the Web help";
	   type=file existing;
	}
        $PHANTOM {
	   label="Phantom realisation";
	   help="Output";
	   type=file;
	}
        $MINVOLUME {
	   label="Minimum volume";
           help="All features have got at least this volume";
	   type=float;
           by default=0;
	}
        OPT(-discrete) {
           label="Discrete values";
           help="Density values are forced to be integer";
        }
      }
   }

   MENU menu_random_phantom {
      "I/O variables"
      $PHANTOM_FAMILY
      $PHANTOM
      "Options"
      OPT(-min_volume)
      OPT(-discrete)
   }
*/
