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

#include <data/args.h>
#include <data/volume.h>
#include <data/image.h>
#include <data/selfile.h>

void Usage();

int main(int argc, char **argv) {
   FileName        fn_input, fn_output, fn_oext, fn_in, fn_out;
   SelFile         SF;
   ImageXmipp      image;
	FourierImageXmipp Fourier_image;
   VolumeXmipp     volume;
   bool Force_Reverse_Flag, image_reversed;

   // Read arguments --------------------------------------------------------
   try {
     fn_input = get_param(argc,argv,"-i",NULL,1,"Reverse endian: Input file not found");
     fn_out   = get_param(argc,argv,"-o","");
     fn_oext  = get_param(argc,argv,"-oext","");
     if(check_param(argc,argv,"-force"))
        Force_Reverse_Flag=true;
      else Force_Reverse_Flag=false;
     if (Is_ImageXmipp(fn_input) ||
	      Is_VolumeXmipp(fn_input) ||
			Is_FourierImageXmipp(fn_input) ) {
       SF.insert(fn_input,SelLine::ACTIVE);
       SF.go_beginning();
     } else  {
       SF.read(fn_input);
     }
   }
   catch (Xmipp_error XE) {cout << XE; Usage(); exit(1);}


   try {
   // Mask a single image ---------------------------------------------------
   if(Is_FourierImageXmipp(fn_input)) {
      Fourier_image.read(fn_input);
      if(Force_Reverse_Flag) image_reversed=true;
      else image_reversed=Fourier_image.reversed();
      if (fn_out=="") Fourier_image.write(fn_input,image_reversed);
      else            Fourier_image.write(fn_out,image_reversed);
	
	} else if (Is_ImageXmipp(fn_input)) {
      image.read(fn_input);
      if(Force_Reverse_Flag==true) image_reversed=true;
      else image_reversed=image.reversed();
      if (fn_out=="") image.write(fn_input,image_reversed);
      else            image.write(fn_out,image_reversed);

   // Mask a single volume --------------------------------------------------
   } else if (Is_VolumeXmipp(fn_input)) {
      volume.read(fn_input);
      if(Force_Reverse_Flag==true) image_reversed=true;
      else image_reversed=image.reversed();
      if (fn_out=="") volume.write(fn_input,image_reversed);
      else            volume.write(fn_out,image_reversed);

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
   		if(Is_FourierImageXmipp(fn_in)) {
      		Fourier_image.read(fn_in);
      		if(Force_Reverse_Flag==true) image_reversed=true;
      		else image_reversed=Fourier_image.reversed();
            Fourier_image.write(fn_out,image_reversed);

			} else  if (Is_ImageXmipp(fn_in)) {
            image.read(fn_in);
            if(Force_Reverse_Flag==true) image_reversed=true;
            else image_reversed=image.reversed();
            image.write(fn_out,image_reversed);
         // Process a volume ...............................................
         } else if (Is_VolumeXmipp(fn_in)) {
            volume.read(fn_in);
            if(Force_Reverse_Flag==true) image_reversed=true;
            else image_reversed=image.reversed();
            volume.write(fn_out,image_reversed);
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
    cerr << "    Reverse the little/big endian status of the input files\n"
         << "    Floats are written in the native format for the machine running\n"
         << "    the program\n";

    cerr << "Usage: reverse_endian <parameters>\n"
         << "   -i <image or volume> [-o <image_out or volume_out]\n"
         << "   -i <selfile> [-oext <output extension>]\n"
         << "   [-force change the little/big endian status\n"
	 << "           even if the results is not the native machine format]\n"
	 << endl;
}


/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Reverse_endian {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Reverse_endian/Help/reverse_endian.html";
      help="Reverses endianness of input images and volumes";
      OPEN MENU menu_reverse_endian;
      COMMAND LINES {
	+ usual: xmipp_reverse_endian
                     #include "prog_line.mnu"
      }
      PARAMETER DEFINITIONS {
         #include "prog_vars.mnu"
      }
   }

   MENU menu_reverse_endian {
      #include "prog_menu.mnu"
   }
*/
