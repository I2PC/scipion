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

#include <reconstruction/reconstruct_art.h>
#include <reconstruction/art_crystal.h>

int main (int argc, char *argv[]) {
// Variables
   Basic_ART_Parameters   art_prm;
   Plain_ART_Parameters   dummy;
   Crystal_ART_Parameters crystal_art_prm;
   VolumeXmipp            vol_voxels;
   GridVolume             vol_blobs;
   GridVolume             *vol_blobs_var=NULL;
   int                    crystal_mode;

// Read Art Parameters
   try {
      art_prm.read(argc, argv);
      // Crystal
      crystal_mode = check_param(argc,argv,"-crystal");
      if (crystal_mode) crystal_art_prm.read(argc, argv, art_prm);
   } catch (Xmipp_error &XE) {
      cout << XE;
      bool usage_more=check_param(argc,argv,"-more_help");
      if (usage_more) {
	 art_prm.usage_more();
	 crystal_art_prm.usage_more();
      } else
         art_prm.usage();
      exit(1);
   }

// Call main ART routine
   try {
      if (!crystal_mode)
         Basic_ROUT_Art(art_prm,dummy,vol_voxels,vol_blobs);
      else
         Basic_ROUT_Art(art_prm,crystal_art_prm,vol_voxels,vol_blobs);
   cerr.flush();
   } catch (Xmipp_error XE) {cout << XE; exit(1);}
   exit(0);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Art {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Art/Help/art.html";
      help="3D reconstruction algorithm from projections";
      OPEN MENU menu_art;
      COMMAND LINES {
	+ usual: art -i $SELFILE -o $FNROOT [-start $STARTVOL] [-sym $SYMFILE]
                 [-l $LAMBDA $STANDARD_PRMS -n $NO_IT] [-stop_at $STOP_AT]
                 [-SIRT] [-random_sort]
                 [-r $BLOB_RADIUS -m $BLOB_ORDER -a $BLOB_ALPHA]
                 [-g $REL_SIZE] [-ext $EXT] [$OTHER_GRIDS]
                 [-show_error] [-show_stats] [-save_at_each_step]
                 [-save_intermidiate] [-save_blobs] [-manual_order]
                 [-crystal -lattice_a "[" $AX "," $AY "]"
                    -lattice_b "[" $BX "," $BY "]"]
                 [-surface $SURFACE_MASK] [-surface_freq $SURFACE_FREQ]
      }
      PARAMETER DEFINITIONS {
        $SELFILE {
	   label="Sel File";
	   help="File with all input projections";
	   type=file existing;
	}
        $FNROOT {
	   label="FileName root";
	   help="A .hist and .vol files are generated";
	   type=file;
	}
        $STARTVOL {
	   label="Initial volume";
	   help="Blob Xmipp format";
	   type=file existing;
	}
        $SYMFILE {
	   label="Symmetry file";
	   help="This file has got a complex structure, better see Web help";
	   type=file existing;
	}
        OPT(-l) {label="Iteration parameters";}
           $LAMBDA {
	      label="Relaxation parameter";
              help="Controls the convergence speed";
	      type=FLOAT [-2...2];
              by default="0.0047";
	   }
           $NO_IT {
              label="Number of iterations";
              help="Usually no more than 1 or 2";
              type=integer [1 ...];
              by default="1";
           }
           $STANDARD_PRMS {
              label="Standard Parameters";
              help="Some useful combinations";
              type=list {
                  "User defined"
                  "Negative staining" {$LAMBDA=0.047;}
                  "Cryo microscopy"   {$LAMBDA=0.0047;}
              };
           }
        $STOP_AT {
           label="Stop after this number of images";
           type=integer;
           by default=0;
        }
        OPT(-SIRT) {label="Run in SIRT mode";}
        OPT(-random_sort) {
           label="Random sort of projections";
           help="Only valid in ART";
        }
        OPT(-r) {label="Blob parameters";}
           $BLOB_RADIUS {
	      label="Blob radius";
              help="in Voxels";
	      type=FLOAT [0...];
              by default=2;
	   }
           $BLOB_ORDER {
              label="Blob order";
              help="Derivative order of the Bessel function, must be a natural";
              type=NATURAL;
              by default=2;
           }
           $BLOB_ALPHA {
              label="Blob smoothness";
              help="The higher the value, the sharper the blob";
              type=FLOAT [0...];
              by default=3.6;
           }
        $REL_SIZE {
           label="Grid spacing: relative size";
           help="in voxels";
           type=float;
           by default=2.26;
        }
        $EXT {
           label="Grid extension: frame size";
           help="in voxels, must be a natural";
           type=NATURAL;
           by default=6;
        }
        $OTHER_GRIDS {
           label="Use other grids than BCC";
           type=exclusion {
              "FCC" {-FCC}
              "CC"  {-CC}
           };
        }
        OPT(-show_error)        {label="Show error on each projection";}
        OPT(-show_stats)        {label="Show statistics on each projection";}
        OPT(-save_at_each_step) {label="Save volumes at each projection";}
        OPT(-save_intermidiate) {label="Save volumes at each iteration";}
        OPT(-save_blobs)        {label="Save blob volumes each time you have to save";}
        OPT(-manual_order)      {label="Projection order is given manually";}
        OPT(-crystal) {label="Crystal reconstruction";}
           $AX {label="Lattice A vector (X component)"; type=FLOAT;}
           $AY {label="Lattice A vector (Y component)"; type=FLOAT;}
           $BX {label="Lattice B vector (X component)"; type=FLOAT;}
           $BY {label="Lattice B vector (Y component)"; type=FLOAT;}
        $SURFACE_MASK {
           label="Surface Volume Mask";
           help="This mask should be generated by the Surface program";
           type=file existing;
        }
        $SURFACE_FREQ {
           label="Insert surface restrictions every";
           help="Number of projections after which the surface restriction
                is applied";
           type=natural;
        }
      }
   }

   MENU menu_art {
      "I/O variables"
      $SELFILE
      $FNROOT
      OPT(-start)
      OPT(-sym)
      "Iteration parameters"
      OPT(-l)
      OPT(-stop_at)
      OPT(-SIRT)
      OPT(-random_sort)
      "Blob parameters"
      OPT(-r)
      "Grid parameters"
      $Art_00007
      OPT(-g)
      OPT(-ext)
      OPT($OTHER_GRIDS)
      "Crystal parameters"
      OPT(-crystal)
      "Surface restrictions"
      OPT(-surface)
      OPT(-surface_freq)
      "Debugging options"
      OPT(-show_error)
      OPT(-show_stats)
      OPT(-save_at_each_step)
      OPT(-save_intermidiate)
      OPT(-save_blobs)
      OPT(-manual_order)
   }
*/
