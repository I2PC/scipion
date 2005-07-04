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
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippImages.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippSelFiles.hh>

/* PROTOTYPES -------------------------------------------------------------- */
void Usage();

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv) {
   FileName        fn_input, fn_output, fn_oext, fn_in, fn_out;
   SelFile         SF;
   ImageXmipp      image;
   VolumeXmipp     volume;
   int             zdim, ydim, xdim;

   // Read arguments --------------------------------------------------------
   try {
     fn_input = get_param(argc,argv,"-i",NULL,1,"Scale: Input file not found");
     fn_out   = get_param(argc,argv,"-o","");
     fn_oext  = get_param(argc,argv,"-oext","");
     if (!Is_ImageXmipp(fn_input) && !Is_VolumeXmipp(fn_input))
       SF.read(fn_input);
     zdim=AtoI(get_param(argc,argv,"-zdim","0"));
     ydim=AtoI(get_param(argc,argv,"-ydim","0"));
     xdim=AtoI(get_param(argc,argv,"-xdim"));
     
     if (ydim==0) ydim=xdim;
     if (zdim==0) zdim=xdim;
   }
   catch (Xmipp_error XE) {cout << XE; Usage(); exit(1);}

   try {
   // Scale a single image -------------------------------------------------
   if (Is_ImageXmipp(fn_input)) {
      image.read(fn_input);
      image().self_scale_to_size_Bspline(3,ydim, xdim);
      if (fn_out=="") image.write(fn_input);
      else            image.write(fn_out);
      
   // Scale a single volume ------------------------------------------------
   } else if (Is_VolumeXmipp(fn_input)) {
      volume.read(fn_input);
      volume().self_scale_to_size_Bspline(3,zdim, ydim, xdim);
      if (fn_out=="") volume.write(fn_input);
      else            volume.write(fn_out);

   // Mask a selection file ------------------------------------------------
   } else {
      SF.read(fn_input);
      
      // Initialise progress bar
      time_config();
      int i=0; init_progress_bar(SF.ImgNo());
      while (!SF.eof()) {
         fn_in=SF.NextImg();
         if (fn_oext=="") fn_out=fn_in;
         else             fn_out=fn_in.without_extension()+"."+fn_oext;
         // Process an image ...............................................
         if (Is_ImageXmipp(fn_in)) {
            image.read(fn_in);
            image().self_scale_to_size_Bspline(3,ydim, xdim);
            image.write(fn_out);
         // Process a volume ...............................................
         } else if (Is_VolumeXmipp(fn_in)) {
            volume.read(fn_in);
            volume().self_scale_to_size_Bspline(3,zdim, ydim, xdim);
            volume.write(fn_out);
         // Not a Spider file ..............................................
         } else
            cout << fn_in << " is not a SPIDER file\n";

         if (i++%25==0) progress_bar(i);
      }
      progress_bar(SF.ImgNo());
   }
   } catch (Xmipp_error XE) {cout << XE;}
   exit(0);
} //main

/* Usage ------------------------------------------------------------------- */
void Usage() {
    cerr << "Purpose:\n";
    cerr << "    Scale images/volumes to a given size\n";
    
    cerr << "Usage: scale <parameters>\n"
         << "   -i <image or volume> [-o <image_out or volume_out]\n"
         << "   -i <selfile> [-oext <output extension>]\n"
         << "   -xdim <new x dimension>\n"
         << "  [-ydim <new y dimension=new x dimension>\n"
         << "  [-zdim <new z dimension=new x dimension>\n";
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Scale {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Scale/Help/scale.html";
      help="Scale volumes and images to a new size";
      OPEN MENU menu_scale;
      COMMAND LINES {
	+ usual: xmipp_scale
               #include "prog_line.mnu"
               -xdim $XDIM [-ydim $YDIM] [-zdim $ZDIM]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        $XDIM {type=natural; label="New X dimension";}
        $YDIM {type=natural; label="New Y dimension"; by default=$XDIM;}
        $ZDIM {type=natural; label="New Z dimension"; by default=$XDIM;}
      }
   }

   MENU menu_scale {
      #include "prog_menu.mnu"
      "Scaling parameters"
      $XDIM
      OPT($YDIM)
      OPT($ZDIM)
   }
*/
