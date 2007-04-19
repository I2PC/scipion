/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (1999)
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
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/volume.h>
#include <data/mask.h>

#include <cstdio>

void Usage();

int main(int argc, char **argv) {
   FileName        fn_input, fn_stats;
   SelFile         SF;
   DocFile         DF_stats;
   ImageXmipp      image;
   VolumeXmippT<double> volume;
   VolumeXmippT<int>    volume_int;

   headerXmipp     header;
   Mask_Params     mask_prm(INT_MASK);
   int             stats;            // True if statistics are to be shown
   int             show_old_rot;     // True if old rot is to be shown
   int             short_format;     // True if a short line is to be shown
   int             save_mask;        // True if the masks must be saved
   int             repair;           // True if headers are initialized
   bool            show_angles;      // True if angles are to be shown in stats
   bool            apply_geo;        // True if the header must be taken into account
   matrix2D<int>   mask2D;
   matrix3D<int>   mask3D;
   int             volume_type;
   // Read arguments --------------------------------------------------------
   try {
     if (argc<2) REPORT_ERROR(1,"");
     fn_input=argv[1];
     if (fn_input.get_extension()!="sel") {
       SF.insert(fn_input,SelLine::ACTIVE);
       SF.go_beginning();
     } else  {
       SF.read(fn_input);
     }
     mask_prm.read(argc,argv);
     stats        = check_param(argc,argv,"-stats");
     if (stats)   fn_stats = get_param(argc,argv,"-o","");
     show_old_rot = check_param(argc,argv,"-show_old_rot");
     short_format = check_param(argc,argv,"-short_format");
     save_mask    = check_param(argc,argv,"-save_mask");
     repair       = check_param(argc,argv,"-repair");
     show_angles  = check_param(argc,argv,"-show_angles");
     apply_geo    =!check_param(argc,argv,"-dont_apply_geo");
   }
   catch (Xmipp_error XE) {cout << XE; Usage(); mask_prm.usage(); exit(1);}

   try {
   DF_stats.append_comment((string)"# Statistics of "+fn_input);
   DF_stats.append_comment("# min max avg stddev");

   // Get maximum filename size ---------------------------------------------
   int max_length=0;
   if (stats) max_length=SF.MaxFileNameLength();

   // Process each file -----------------------------------------------------
   #define V VOLMATRIX(volume)
   #define VI VOLMATRIX(volume_int)
   #define I IMGMATRIX(image)
   double min_val, max_val, avg, stddev;
   int min_val_int, max_val_int;
   double mean_min_val=0, mean_max_val=0, mean_avg=0, mean_stddev=0;
   int N=0;

   if (stats && short_format) {
      cout << "Format: Name ZxYxX min max avg stddev ";
      if (show_angles) cout << " <rot tilt psi>";
      if (show_old_rot) cout << " old_rot";
      cout << '>' << endl;
   }

   while (!SF.eof()) {
     FileName file_name = SF.NextImg();

     // For volumes ........................................................
     if( (volume_type=Is_VolumeXmipp(file_name)) ) {
        // Read file
        if (repair || stats) {
           if(volume_type ==headerXmipp::VOL_XMIPP)
             {volume.read(file_name);
              volume().set_Xmipp_origin();
              }
            else if(volume_type ==headerXmipp::VOL_INT)
            {volume_int.read(file_name);
             volume_int().set_Xmipp_origin();
            }
        }

        // Repair?
        if (repair)
           {
           if(volume_type == headerXmipp::VOL_XMIPP)
               {volume.clear_header(); volume.write();}
           else if(volume_type == headerXmipp::VOL_INT)
               {volume_int.clear_header(); volume_int.write();}
           }
        // Statisitics
        if (stats) {
           // Generate mask if necessary
           mask_prm.generate_3Dmask(V);
           const matrix3D<int> &mask3D=mask_prm.get_binary_mask3D();

           if(volume_type == headerXmipp::VOL_XMIPP)
               compute_stats_within_binary_mask(mask3D, V, min_val, max_val,
              avg, stddev);
           else if(volume_type == headerXmipp::VOL_INT)
               compute_stats_within_binary_mask(mask3D, VI, min_val_int,
                max_val_int, avg, stddev);

           // Show information
           cout << AtoA(file_name,max_length+1);
           cout << ItoA(ZSIZE(V),4,' ') << 'x'
                << ItoA(YSIZE(V),4,' ') << 'x'
                << ItoA(XSIZE(V),4,' ') << ' ';
           if (!short_format)
              cout << "min= "    << FtoA(min_val,10) << ' '
                   << "max= "    << FtoA(max_val,10) << ' '
                   << "avg= "    << FtoA(avg    ,10) << ' '
                   << "stddev= " << FtoA(stddev ,10) << ' ';
           else
              cout << FtoA(min_val,10) << ' '
                   << FtoA(max_val,10) << ' '
                   << FtoA(avg    ,10) << ' '
                   << FtoA(stddev ,10) << ' ';
            matrix1D<double> v(4);
            v(0)=min_val;
            v(1)=max_val;
            v(2)=avg;
            v(3)=stddev;
            DF_stats.append_data_line(v);

            // Total statistics
            N++;
            mean_min_val += min_val;
            mean_max_val += max_val;
            mean_avg     += avg;
            mean_stddev  += stddev;
         } else {
            header.read(file_name);
            cout << "FileName     : " << file_name << endl;
            cout << header;
         }

     // For images .........................................................
     } else if(Is_ImageXmipp(file_name)) {
        // Read file
        if (repair || stats) {
           // Read the image applying the header
           image.read(file_name,false,false,apply_geo);
           image().set_Xmipp_origin();
        }

        // Repair?
        if (repair) {image.clear_header(); image.write();}

        // Statisitics
        if (stats) {
           // Generate mask if necessary
           mask_prm.generate_2Dmask(I);
           const matrix2D<int> &mask2D=mask_prm.get_binary_mask2D();

           // Compute statistics
           compute_stats_within_binary_mask(mask2D, I, min_val, max_val,
              avg, stddev);

           // Show information
           cout << AtoA(file_name,max_length+1);
           cout << "    "; // Stands for the ZSIZE in volumes
           cout << ItoA(YSIZE(I),4,' ') << 'x'
                << ItoA(XSIZE(I),4,' ') << ' ';
           if (!short_format) {
              cout << "min= "    << FtoA(min_val,10) << ' '
                   << "max= "    << FtoA(max_val,10) << ' '
                   << "avg= "    << FtoA(avg    ,10) << ' '
                   << "stddev= " << FtoA(stddev ,10) << ' ';
              if (show_angles) {
                 cout << "rot= "    << FtoA(image.rot() ,10) << ' '
                      << "tilt= "   << FtoA(image.tilt(),10) << ' '
                      << "psi= "    << FtoA(image.psi() ,10) << ' ';
                 if(image.Is_flag_set()==1.0f || image.Is_flag_set()==2.0f)
                 cout << "\nrot1= "  << FtoA(image.rot1() ,10) << ' '
                      << "tilt1= "   << FtoA(image.tilt1(),10) << ' '
                      << "psi1= "    << FtoA(image.psi1() ,10) << ' ';
                 if(image.Is_flag_set()==2.0f)
                 cout << "\nrot2= "    << FtoA(image.rot2() ,10) << ' '
                      << "tilt2= "   << FtoA(image.tilt2(),10) << ' '
                      << "psi2= "    << FtoA(image.psi2() ,10) << ' ';
              }

              if (show_old_rot)
                 cout << "old_rot= " << FtoA(image.old_rot() ,10);
           } else {
              cout << FtoA(min_val,10) << ' '
                   << FtoA(max_val,10) << ' '
                   << FtoA(avg    ,10) << ' '
                   << FtoA(stddev ,10) << ' ';
              if (show_angles) {
                 cout << FtoA(image.rot() ,10) << ' '
                      << FtoA(image.tilt(),10) << ' '
                      << FtoA(image.psi() ,10) << ' ';
                 if(image.Is_flag_set()==1.0f || image.Is_flag_set()==2.0f)
                 cout << FtoA(image.rot1() ,10) << ' '
                      << FtoA(image.tilt1(),10) << ' '
                      << FtoA(image.psi1() ,10) << ' ';
                 if(image.Is_flag_set()==2.0f)
                 cout << FtoA(image.rot2() ,10) << ' '
                      << FtoA(image.tilt2(),10) << ' '
                      << FtoA(image.psi2() ,10) << ' ';
              }
              if (show_old_rot)
                 cout << FtoA(image.old_rot() ,10);
           }

           matrix1D<double> v(4);
           v(0)=min_val;
           v(1)=max_val;
           v(2)=avg;
           v(3)=stddev;
           DF_stats.append_data_line(v);

           // Total statistics
           N++;
           mean_min_val += min_val;
           mean_max_val += max_val;
           mean_avg     += avg;
           mean_stddev  += stddev;
         } else {
            cout << "FileName     : " << file_name << endl;
            header.read(file_name);
            cout << header;
            if (show_old_rot)
               cout << "Old rot      : " << header.old_rot() << endl;
            cout << endl;
         }
     // For fourier volumes .................................................
     } else if (Is_FourierVolumeXmipp(file_name)) {
         header.read(file_name);
         cout << "FileName     : " << file_name << endl;
         cout << header;

     // For fourier images .................................................
     } else if (Is_FourierImageXmipp(file_name)) {
         header.read(file_name);
         cout << "FileName     : " << file_name << endl;
         cout << header;

     // Is not an Spider file ..............................................
     } else
       cout << file_name<<" is not a Spider File";

     // Finish information .................................................
     cout << endl;

   } // while

   // Show total statistics ------------------------------------------------
   if (stats) {
      cout << "==================================================\n";
      cout << "Total number of images/volumes: " << N << endl;
      if (N!=0) {
         mean_min_val /= N;
         mean_max_val /= N;
         mean_avg     /= N;
         mean_stddev  /= N;

         cout << AtoA(" ",max_length+13);
         if (!short_format)
            cout << "min= "    << FtoA(mean_min_val,10) << ' '
                 << "max= "    << FtoA(mean_max_val,10) << ' '
                 << "avg= "    << FtoA(mean_avg    ,10) << ' '
                 << "stddev= " << FtoA(mean_stddev ,10) << ' ';
         else
            cout << FtoA(mean_min_val,10) << ' '
                 << FtoA(mean_max_val,10) << ' '
                 << FtoA(mean_avg    ,10) << ' '
                 << FtoA(mean_stddev ,10) << ' ';
         cout << endl;
      }
   }

   // Save masks -----------------------------------------------------------
   if (save_mask) {
      mask_prm.write_2Dmask("mask2D");
      mask_prm.write_3Dmask("mask3D");
   }

   // Save statistics ------------------------------------------------------
   if (fn_stats!="") DF_stats.write(fn_stats);
   } catch (Xmipp_error XE) {cout << XE;}
   exit(0);
} //main

/* Usage ------------------------------------------------------------------- */
void Usage() {
    cerr << "Purpose:\n";
    cerr << "    Displays some of the information stored in the header\n";
    cerr << "    It always calculate max, min and average, even if it is\n";
    cerr << "    stored in the header\n";

    cerr << "Usage: infogeo <Xmipp file or SelFile>" << endl
         << "   [-show_old_rot]   : show old rotational angle\n"
         << "   [-stats]          : show image statistics\n"
         << "      [-o <docfile>] : save the statistics in this docfile\n"
         << "      [-dont_apply_geo]: do not apply geo when the image is read\n"
         << "   [-short_format]   : Don't show labels for statistics\n"
         << "   [-save_mask]      : save 2D and 3D masks with names: \n"
         << "                          mask2D and mask3D respectively\n"
         << "   [-repair]         : reset header\n";
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Infogeo {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Infogeo/Help/infogeo.html";
      help="Display header and statistical information of Xmipp files";
      OPEN MENU Infogeo;
      COMMAND LINES {
         + usual: infogeo $FILE_IN
               #include "binary_mask_line.mnu"
               [-show_old_rot] [-stats] [-short_format] [-save_mask] [-repair]
      }
      PARAMETER DEFINITIONS {
         $FILE_IN {
            label="Input file";
            help="Sel File, Image or Volume";
            type=FILE EXISTING;
         }

         #include "binary_mask_vars.mnu"

         OPT(-show_old_rot) {label="Show Old Xmipp Rotational Angle";}
         OPT(-stats) {label="Show statistics";}
         OPT(-short_format) {label="Show statistics in short format";}
         OPT(-save_mask) {label="Save applied mask"; help="As mask2D or mask3D";}
         OPT(-repair) {label="Clear file SPIDER header";}
      }
   }
   MENU Infogeo {
      "I/O Parameters"
      $FILE_IN
      "Other options"
      OPT(-show_old_rot)
      OPT(-stats)
      OPT(-short_format)
      OPT(-save_mask)
      OPT(-repair)
      #include "binary_mask_menu.mnu"
   }
*/
