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
#include <data/image.h>
#include <data/volume.h>
#include <data/mask.h>
#include <data/histogram.h>

void Usage();

int main(int argc, char **argv)
{
    VolumeXmipp     volume;
    ImageXmipp      image;
    bool            image_mode;
    FileName        fn_in, fn_out;
    Mask_Params     mask_prm(INT_MASK);
    bool            automatic_range;
    double          m, M; // range for histogram
    int             StepsNo;
    histogram1D     hist;
#define         V VOLMATRIX(volume)
#define         I IMGMATRIX(image)

    // Read arguments --------------------------------------------------------
    try
    {
        fn_in      = getParameter(argc, argv, "-i");
        fn_out     = getParameter(argc, argv, "-o", "");

        StepsNo = textToInteger(getParameter(argc, argv, "-steps", "100"));
        int i;
        if ((i = paremeterPosition(argc, argv, "-range")) != -1)
        {
            if (i + 2 >= argc)
                EXIT_ERROR(1, "Histogram: Not enough parameters behind -range\n");
            m = textToFloat(argv[i+1]);
            M = textToFloat(argv[i+2]);
            automatic_range = false;
        }
        else automatic_range = true;

        mask_prm.read(argc, argv);
    }
    catch (Xmipp_error Xe)
    {
        cout << Xe;
        Usage();
        mask_prm.usage();
        exit(1);
    }

    try
    {
        if (Is_ImageXmipp(fn_in))
        {
            image_mode = true;
            image.read(fn_in);
        }
        else if (Is_VolumeXmipp(fn_in))
        {
            image_mode = false;
            volume.read(fn_in);
        }
        else EXIT_ERROR(1, "Histogram: Input file is not an image nor a volume");

        // Compute histogram ----------------------------------------------------
        if (image_mode)
        {
            I.setXmippOrigin();
            mask_prm.generate_2Dmask(I);
            const Matrix2D<int> & mask2D = mask_prm.get_binary_mask2D();
            if (automatic_range)
                compute_hist_within_binary_mask(mask2D, I, hist, StepsNo);
            else
                compute_hist_within_binary_mask(mask2D, I, hist, m, M, StepsNo);
        }
        else
        {
            V.setXmippOrigin();
            mask_prm.generate_3Dmask(V);
            const Matrix3D<int> & mask3D = mask_prm.get_binary_mask3D();
            if (automatic_range)
                compute_hist_within_binary_mask(mask3D, V, hist, StepsNo);
            else
                compute_hist_within_binary_mask(mask3D, V, hist, m, M, StepsNo);
        }

        if (fn_out != "") hist.write(fn_out);
        else            cout << hist;
    }
    catch (Xmipp_error Xe)
    {
        cout << Xe;
    }
    exit(0);
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    cout << "histogram [Parameters]:\n"
    << "   -i <File_in>                : Input Xmipp Volume or Image\n"
    << "  [-o <File_out>]              : Text file with histogram\n"
    << "  [-range <m> <M>]             : range for the histogram\n"
    << "                                 by default, it is automatic\n"
    << "  [-steps <N=100>]             : number of subdivisions\n";
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Histogram {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Histogram/Help/histogram.html";
      help="Compute histograms of images and volumes within a mask";
      OPEN MENU Histogram;
      COMMAND LINES {
         + usual: histogram -i $FILE_IN [-o $FILE_OUT]
             [-range $MIN $MAX]  [-steps $STEPS]
             #include "binary_mask_line.mnu"
      }
      PARAMETER DEFINITIONS {
         $FILE_IN {
            label="Input file";
            help="Image or volume";
            type=FILE EXISTING;
         }
         $FILE_OUT {
            label="Output Text file";
            type=FILE;
            help="Column 1: density value; Column 2: voxel count";
         }
         OPT(-range) {label="Range";}
            $MIN {label="Minimum"; type=FLOAT;}
            $MAX {label="Maximum"; type=FLOAT;}
         $STEPS {label="Steps"; type=NATURAL; by default=100;}
         #include "binary_mask_vars.mnu"
      }
   }
   MENU Histogram {
      "I/O Parameters"
      $FILE_IN
      OPT(-o)
      "Histogram Specification"
      OPT(-range)
      OPT(-steps)
      #include "binary_mask_menu.mnu"
   }
*/
