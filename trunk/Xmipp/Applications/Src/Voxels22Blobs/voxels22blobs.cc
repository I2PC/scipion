/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2000)
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

/* INCLUDES ---------------------------------------------------------------- */
#include <Reconstruction/blobs.hh>
#include <Reconstruction/grids.hh>

/* PROTOTYPES -------------------------------------------------------------- */
void Usage();

/* MAIN -------------------------------------------------------------------- */
int main (int argc,char *argv[]) {
   FileName fn_in, fn_out;
   bool   voxels_to_blobs;
   struct blobtype blob;
   float lambda;
   float final_error;
   float grid_relative_size;
   #define  CC 0
   #define FCC 1
   #define BCC 2
   int grid_type;
   double R;

// Read the command line ---------------------------------------------------
   try {
      if (check_param(argc,argv,"-voxels")) {
         voxels_to_blobs=true;
         fn_in=get_param(argc,argv,"-voxels");
         grid_relative_size = AtoF(get_param(argc, argv, "-g", "1.41"));
         if      (check_param(argc, argv, "-FCC")) grid_type=FCC;
         else if (check_param(argc, argv, "-CC"))  grid_type=CC;
         else                                      grid_type=BCC;
      } else if (check_param(argc,argv,"-blobs")) {
         voxels_to_blobs=false;
         fn_in=get_param(argc,argv,"-blobs");
      } else
         REPORT_ERROR(1,"Voxels22blobs: Not recognised input file type");
      fn_out=get_param(argc,argv,"-o");
      lambda             = AtoF(get_param(argc, argv, "-l",    "0.05" ));
      final_error        = AtoF(get_param(argc, argv, "-final_error","0.01"));
      blob.radius        = AtoF(get_param(argc, argv, "-r",    "2"    ));
      blob.order         = AtoI(get_param(argc, argv, "-m",    "2"    ));
      blob.alpha         = AtoF(get_param(argc, argv, "-a",    "10.4"  ));
      R                  = AtoF(get_param(argc, argv, "-R",    "-1"   ));
   } catch (Xmipp_error &XE) {cout << XE; Usage(); exit(1);}
try {
// Really convert ----------------------------------------------------------
   GridVolume  vol_blobs;
   VolumeXmipp vol_voxels;
   if (voxels_to_blobs) {
      vol_voxels.read(fn_in);
      vol_voxels().set_Xmipp_origin();
      voxels2blobs(&(vol_voxels()), blob, vol_blobs, grid_type,
         grid_relative_size, lambda, NULL, NULL, final_error, false, R);
      vol_blobs.write(fn_out);
   } else {
      vol_blobs.read(fn_in);
      blobs2voxels(vol_blobs, blob, &(vol_voxels()));
      vol_voxels.write(fn_out);
   }
} catch (Xmipp_error XE) {cout << XE;}

}

/* Usage ------------------------------------------------------------------- */
void Usage() {
   cout << "Usage: Voxels22blobs [Parameters]\n"
        << "   (-voxels | -blobs) <file_in>       : Input file\n"
        << "    -o <file_out>                     : of the opposite type\n"
        << "   [-r <blob radius=2>]               : blob radius\n"
        << "   [-m <blob order=2>]                : blob derivative order\n"
        << "   [-a <blob alpha=10.4>]             : controls smoothness\n"
        << "Only if voxels:\n"
        << "   [-g <grid_relative_size=1.41>]     : size between grid samples\n"
        << "   [-FCC | -CC]                       : by default, BCC grid\n"
        << "   [-l <lambda=0.05>]                 : convergence rate\n"
        << "   [-final_error <error=0.01>]        : minimum change percentage\n"
        << "   [-R <R=-1>]                        : interest radius\n";
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Voxels22Blobs {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Window/Help/window.html";
      help="Convert a blob volume to voxels and viceversa";
      OPEN MENU menu_voxels22blobs;
      COMMAND LINES {
	+ to_blobs: voxels22blobs -voxels $FILE_IN -o $FILE_OUT
                       [-r $BLOB_RADIUS -m $BLOB_ORDER -a $BLOB_ALPHA]
                       [-g $REL_SIZE] [$OTHER_GRIDS]
                       [-l $LAMBDA -final_error $ERROR]
        + to_voxels: voxels22blobs -blobs $FILE_IN -o $FILE_OUT
                       OPT(-r)
      }
      PARAMETER DEFINITIONS {
        $FILE_IN {
	   label="Input file";
	   help="Either voxel or blob volume";
	   type=file existing;
	}
        $FILE_OUT {
	   label="Output file";
	   help="Either blob or voxel volume";
	   type=file;
	}
        OPT(-r) {label="Blob shape";}
           $BLOB_RADIUS {
	      label="Blob radius";
              help="in Voxels";
	      type=FLOAT [0...];
              by default="2";
	   }
           $BLOB_ORDER {
              label="Blob order";
              help="Derivative order of the Bessel function, must be a natural";
              type=NATURAL;
              by default="2";
           }
           $BLOB_ALPHA {
              label="Blob smoothness";
              help="The higher the value, the sharper the blob";
              type=FLOAT [0...];
              by default="3.6";
           }
        OPT(-g) {label="Grid Spacing";}
           $REL_SIZE {
              label="Relative size";
              help="in voxels";
              type=float [0...];
              by default="2.26";
           }
        $OTHER_GRIDS {
           label="Use other grids than BCC";
           type=exclusion {
              "FCC" {-FCC}
              "CC"  {-CC}
           };
        }
        OPT(-l) {label="Iteration parameters";}
           $LAMBDA {
	      label="Relaxation parameter";
              help="Controls the convergence speed";
	      type=FLOAT [-2...2];
              by default="0.05";
	   }
           $ERROR {
              label="Final error";
              help="If the mean change in the volume is less than this value
                    as a percentage, the conversion is over";
              type=float [0...];
              by default="0.01";
           }
      }
   }

   MENU menu_voxels22blobs {
      "I/O parameters"
      $FILE_IN
      $FILE_OUT
      "Only Voxels --> Blobs"
      OPT(-g)
      OPT($OTHER_GRIDS)
      OPT(-l)
      "Blob shape"
      OPT(-r)
   }
*/
