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

#include <XmippData/xmippMicrograph.hh>
#include <XmippData/xmippArgs.hh>

void Usage();

int main (int argc, char **argv) {
   FileName fn_orig, fn_micrograph, fn_pos, fn_root;
   int      Ydim, Xdim;
   int      startN;
   bool     reverse_endian;
   bool     compute_transmitance=FALSE;
   bool     compute_inverse=FALSE;
   double   alpha;
   try {
      fn_micrograph = get_param(argc,argv,"-i");
      fn_pos        = get_param(argc,argv,"-pos");
      fn_root       = get_param(argc,argv,"-root");
      fn_orig       = get_param(argc,argv,"-orig","");
      Xdim          = AtoI(get_param(argc,argv,"-Xdim"));
      alpha         = AtoF(get_param(argc,argv,"-alpha","0"));
      if (check_param(argc,argv,"-Ydim"))
         Ydim       = AtoI(get_param(argc,argv,"-Ydim"));
      else Ydim=Xdim;
      startN        = AtoI(get_param(argc,argv,"-start","1"));
      reverse_endian= check_param(argc,argv,"-reverse_endian");
      compute_inverse= check_param(argc,argv,"-invert");
      compute_transmitance= check_param(argc,argv,"-log");
   } catch (Xmipp_error XE) {cout << XE; Usage(); exit(1);}
   try {
      Micrograph m;
      m.open_micrograph(fn_micrograph,reverse_endian);
      m.set_window_size(Xdim,Ydim);
      m.read_coordinates(0,fn_pos);
      m.add_label("");
      m.set_transmitance_flag(compute_transmitance);
      m.set_inverse_flag(compute_inverse);
      m.produce_all_images(0,fn_root,startN,fn_orig,alpha);
cout << "compute_transmitance " << compute_transmitance
     << " compute_inverse " << compute_inverse << endl;
   } catch (Xmipp_error XE) {cout << XE;}
}

void Usage() {
   cerr << "Purpose: Cut the images marked with xmipp_mark\n"
        << "Usage: scissor [options]\n"
        << "   -i <input micrograph>      : From which the images will be cutted\n"
        << "  [-orig <original micrograph>: unless this parameter is specified\n"
        << "   -root <root name>          : for the cutted images\n"
        << "   -pos <position file>       : order X,Y\n"
        << "   -Xdim <window X dim>       : in pixels\n"
        << "  [-Ydim <window Y dim>]      : if not given Ydim=Xdim\n"
        << "  [-start <N=1>]              : Number of the first image\n"
        << "  [-reverse_endian]           : of the input micrographs\n"
	<< "  [-alpha <ang>]              : Angle from Y axis to tilt axis\n"
	<< "                                as it comes out from xmipp_mark\n"
	<< "  [-invert]                   : Invert contrast\n"
	<< "  [-log]                      : Compute optical density\n"
	<< "                                from transmitance\n"
   ;
}

/* Colimate menu =========================================================== */
/*Colimate:
   PROGRAM Scissor {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Scissor/Help/scissor.html";
      help="Cut particles from raw micrographs";
      OPEN MENU Scissor;
      COMMAND LINES {
         + usual: xmipp_scissor -i $FILEIN [-orig $FILEORIG]
                   -root $ROOT -pos $POSFILE
                   -Xdim $XDIM [-Ydim $YDIM]
                  [-start $START] [-reverse_endian]
		  [-alpha $ALPHA]
      }
      PARAMETER DEFINITIONS {
         $FILEIN {
            label="Input micrograph";
            help="Raw file";
            type=file existing;
         }
         $FILEORIG {
            label="Original micrograph";
            help="If supplied, particles will be cut from this file
                  although the positions are referrend to the Input
                  micrograph";
            type=file existing;
         }
         $ROOT {
            label="Output file root name";
            help="Cut images will start with this name";
            type=text;
         }
         $POSFILE {
            label="File with positions";
            help="Positions referred to the Input micrograph";
            type=file existing;
         }
         $XDIM {
            label="Output X dim for particles";
            type=natural;
         }
         $YDIM {
            label="Output Y dim for particles";
            type=natural;
            by default=$XDIM;
         }
         $START {
            label="Start numbering images at";
            type=natural;
            by default=1;
         }
         OPT(-reverse_endian) {label="Reverse endian";}
         $ALPHA {
            label="Angle from Y axis to tilt axis";
            type=float;
	    help="As it comes out from xmipp_mark";
         }
      }
   }
   MENU Scissor {
      "I/O parameters"
      $FILEIN
      OPT($FILEORIG)
      $ROOT
      $POSFILE
      "Cut images definition"
      $XDIM
      OPT($YDIM)
      OPT($START)
      OPT(-reverse_endian)
      OPT($ALPHA)
   }
*/

