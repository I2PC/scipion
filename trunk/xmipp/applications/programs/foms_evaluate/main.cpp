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

#include <reconstruction/foms_evaluate.h>

int main (int argc, char *argv[]) {
   Prog_Evaluate_Parameters prog_prm;
   EVALUATE_results         results;

// Read parameters from command line
   try {prog_prm.read(argc,argv);}
   catch (Xmipp_error XE) {cout << XE; prog_prm.usage(); exit(1);}

// Evaluate
   try {ROUT_Evaluate(prog_prm, results);}
   catch (Xmipp_error XE) {cout << XE;}
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Evaluate {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Evaluate/Help/evaluate.html";
      help="Evaluate FOMs using a phantom and a reconstruction";
      OPEN MENU menu_evaluate;
      COMMAND LINES {
	+ single: evaluate -p $PHANTOM -r $RECONS [-mass $MASS] [-R $R]
                  [-dir $ROT $TILT $AXIS] [-back_radius $r]
                  [-back_factor $FACTOR] [-save_maps]
                  [-show_values] [-show_process] [-save_histograms]
                  [-only_structural]
        + multiple: evaluate -sel $SELFILE OPT(-mass) OPT(-R)
                  OPT(-dir) OPT(-back_radius)
                  OPT(-back_factor) OPT(-save_maps)
                  OPT(-show_values) OPT(-show_process) OPT(-save_histograms)
                  OPT(-only_structural)
      }
      PARAMETER DEFINITIONS {
        $PHANTOM {
	   label="Phantom Description File";
	   help="This file has got a complex structure, better see
                 the Web help";
	   type=file existing;
	}
        $RECONS {
	   label="Reconstructed file";
	   help="Xmipp format";
	   type=file existing;
	}
        $MASS {
           label="Percentage of mass out slice histograms";
	   type=float;
           by default="99";
	}
        $R {
	   label="Global Radius Mask";
	   help="The volume is only evaluated in this sphere";
	   type=FLOAT;
	}
        OPT(-dir) {
           label="Directional FOMs";
           help="For the Radon Transform";
        }
           $ROT  {type=FLOAT; label="Rotational angle";}
           $TILT {type=FLOAT; label="Tilting angle";}
           $AXIS {
               label="System axis";
               type=LIST {
                  "X" {$ROT=0;  $TILT=90;}
                  "Y" {$ROT=90; $TILT=90;}
                  "Z" {$ROT=0;  $TILT=0;}
               };
               by default=2;
           }
        $r {
	   label="Surrounding ackground radius";
           help="For all histogram measures";
	   type=FLOAT;
	}
        $FACTOR {
	   label="Background enlarging factor";
	   type=FLOAT;
           help="For all histogram measures";
           by default=1.25;
	}
        OPT(-save_maps) {
           label="Save Intermidiate maps";
           help="Squared errors, distances, ...";
        }
        OPT(-show_values) {
           label="Show value inside features";
           help="Interactive";
        }
        OPT(-show_process) {
           label="Show process";
           help="Show a lot of information during the process";
        }
        OPT(-save_histograms) {
           label="Save histograms";
           help="Slice and feature histograms";
        }
        OPT(-only_structural) {
           label="Only structural";
           help="Compute only structural FOMs";
        }
        $SELFILE {
            label="Selection file";
            help="Only the reconstructions in this file";
            type=file;
        }
      }
   }

   MENU menu_evaluate {
      "I/O variables"
      $PHANTOM
      $RECONS
      $SELFILE
      "FOMs definition"
      OPT(-mass)
      OPT(-R)
      OPT(-dir)
      OPT(-back_radius)
      OPT(-back_factor)
      "Debugging options"
      OPT(-save_maps)
      OPT(-show_values)
      OPT(-show_process)
      OPT(-save_histograms)
      OPT(-only_structural)
   }
*/
