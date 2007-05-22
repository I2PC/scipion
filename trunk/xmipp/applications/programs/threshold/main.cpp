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

#include <data/volume.h>
#include <data/image.h>
#include <data/args.h>

void Usage();

#define SET_SUBS_VAL(I, subs_val, str_subs_val) \
    if      (str_subs_val=="min") subs_val=I.compute_min(); \
    else if (str_subs_val=="max") subs_val=I.compute_max(); \
    else                          subs_val=AtoF(str_subs_val);

int main(int argc, char **argv)
{
    VolumeXmipp     V, Vdist, Vlabel;
    ImageXmipp      I;
    bool            image_mode;
    FileName        fn_in, fn_out, fn_dist, fn_label;
    double          th_below, th_above, th;
    double          dmin,     dmax;
    int             enable_th_below, enable_th_above;
    int             enable_dmin,     enable_dmax;
    int             binarize;
    bool            enable_substitute;
    double          new_val,     old_val;
    string          str_new_val, str_old_val;
    double          accuracy;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_in      = get_param(argc, argv, "-i");
        fn_out     = get_param(argc, argv, "-o", fn_in.c_str());
        fn_dist    = get_param(argc, argv, "-dist", "");
        fn_label   = get_param(argc, argv, "-label", "");
        if ((enable_th_below = check_param(argc, argv, "-below")))
            th_below = AtoF(get_param(argc, argv, "-below"));
        if ((enable_th_above = check_param(argc, argv, "-above")))
            th_above = AtoF(get_param(argc, argv, "-above"));
        if ((enable_dmin = check_param(argc, argv, "-dmin")))
            dmin = AtoF(get_param(argc, argv, "-dmin"));
        if ((enable_dmax = check_param(argc, argv, "-dmax")))
            dmax = AtoF(get_param(argc, argv, "-dmax"));
        binarize = check_param(argc, argv, "-binarize");
        if (binarize)
        {
            if (enable_th_above) th = th_above;
            else if (enable_th_below) th = th_below;
            else                      th = 0;
        }
        int i;
        if ((i = position_param(argc, argv, "-substitute")) != -1)
        {
            enable_substitute = true;
            if (i + 2 >= argc)
                EXIT_ERROR(1, "Threshold: Not enough parameters behind -substitute\n");
            str_old_val = argv[i+1];
            str_new_val = argv[i+2];
            accuracy = AtoF(get_param(argc, argv, "-accuracy", "0"));
        }
        else enable_substitute = false;
    }
    catch (Xmipp_error Xe)
    {
        cout << Xe;
        Usage();
        exit(1);
    }

    try
    {
        if (Is_ImageXmipp(fn_in))
        {
            image_mode = true;
            I.read(fn_in);
        }
        else if (Is_VolumeXmipp(fn_in))
        {
            image_mode = false;
            V.read(fn_in);
        }
        else EXIT_ERROR(1, "Threshold: Input file is not an image nor a volume");

        // Apply substitution ---------------------------------------------------
        if (enable_substitute)
            if (image_mode)
            {
                SET_SUBS_VAL(I(), new_val, str_new_val);
                SET_SUBS_VAL(I(), old_val, str_old_val);
                I().substitute(old_val, new_val, accuracy);
            }
            else
            {
                SET_SUBS_VAL(V(), new_val, str_new_val);
                SET_SUBS_VAL(V(), old_val, str_old_val);
                V().substitute(old_val, new_val, accuracy);
            }

        // Apply density restrictions -------------------------------------------
        if (image_mode)
        {
            if (enable_th_below) I().threshold("below", th_below, th_below);
            if (enable_th_above) I().threshold("above", th_above, th_above);
        }
        else
        {
            if (enable_th_below) V().threshold("below", th_below, th_below);
            if (enable_th_above) V().threshold("above", th_above, th_above);
        }

        // Apply distance restrictions ------------------------------------------
        if (!image_mode && fn_dist != "")
        {
            Vdist.read(fn_dist);
            if (fn_label == "")
                EXIT_ERROR(1, "Threshold: You must supply a label volume\n");
            Vlabel.read(fn_label);
            FOR_ALL_ELEMENTS_IN_MATRIX3D(Vdist())
            {
                if (Vlabel(k, i, j) > 0) Vdist(k, i, j) *= -1;
                if (enable_dmin && Vdist(k, i, j) < dmin) V(k, i, j) = 0;
                if (enable_dmax && Vdist(k, i, j) > dmax) V(k, i, j) = 0;
            }
        }

        // Binarize -------------------------------------------------------------
        if (binarize)
            if (image_mode) I().binarize(th);
            else            V().binarize(th);

        // Write output volume --------------------------------------------------
        if (image_mode) I.write(fn_out);
        else            V.write(fn_out);
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
    cout << "threshold [Parameters]:\n"
    << "   -i <File_in>                   : Input Xmipp Volume or Image\n"
    << "  [-o <File_out>]                 : If not given, the same as input\n"
    << "  [-dist <Distance volume>        : of the same size as the input\n"
    << "   -label <Label volume>          : of the same size as the input\n"
    << "  [-dmin <dmin>]                  : remove voxels whose distance is smaller\n"
    << "  [-dmax <dmax>]]                 : remove voxels whose distance is greater\n"
    << "  [-below <th>]                   : remove voxels below this threshold\n"
    << "  [-above <th>]                   : remove voxels above this threshold\n"
    << "  [-binarize]]                    : binarize output\n"
    << "  [-substitute <old_val> <new_val>: where a value can be\n"
    << "                                    (<val>|min|max)\n"
    << "    [-accuracy <accuracy=0>]]     : when substituting this value\n"
    << "                                    determines if two values are the same\n";
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Threshold {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Threshold/Help/threshold.html";
      help="Apply a threshold and binarize an image or volume, even with
            distance restrictions";
      OPEN MENU Threshold;
      COMMAND LINES {
         + usual: threshold -i $FILE_IN [-o $FILE_OUT]
             [-below $BELOW_TH] [-above $ABOVE_TH]
             [-binarize]
             [-dist $DISTVOL -label $LABELVOL [-dmin $DMIN] [-dmax $DMAX]]
             [-substitute $OLDVAL $OLDMINMAX $NEWVAL $NEWMINMAX
                [-accuracy $ACCURACY]]
      }
      PARAMETER DEFINITIONS {
         $FILE_IN {
            label="Input file";
            help="Image or volume";
            type=FILE EXISTING;
         }
         $FILE_OUT {
            label="Output file";
            help="If not given, the same as input";
            type=FILE;
         }
         $BELOW_TH {
            label="Below threshold";
            help="All values below this threshold are removed";
            type=FLOAT;
         }
         $ABOVE_TH {
            label="Above threshold";
            help="All values above this threshold are removed";
            type=FLOAT;
         }
         OPT(-binarize) {
            label="Binarize output file";
            help="Non zero values are set to 1";
         }
         OPT(-dist) {
            label="Use distance conditions";
            help="Only for volumes";
         }
            $DISTVOL {
               label="Distance volume";
               help="Minimum distance to a feature/background";
               type=FILE EXISTING;
            }
            $LABELVOL {
               label="Label volume";
               help="Phantom segmentation";
               type=FILE EXISTING;
            }
            $DMIN {
               label="Minimum distance";
               help="Voxels whose distance is smaller than this one are removed";
               type=FLOAT;
            }
            $DMAX {
               label="Maximum distance";
               help="Voxels whose distance is greater than this one are removed";
               type=FLOAT;
            }
         OPT(-substitute) {
            label="Substitute a value by another";
            help="This is the first step before thresholding";
         }
           $OLDVAL {
              label="Old value";
              type=text;
              by default="0";
           }
           $OLDMINMAX {
              label="Value type";
              type=list {
                 "User Defined"
                 "min" {$OLDVAL="min";}
                 "max" {$OLDVAL="max";}
              };
           }
           $NEWVAL {
              label="New value";
              type=text;
              by default="0";
           }
           $NEWMINMAX {
              label="Value type";
              type=list {
                 "User Defined"
                 "min" {$NEWVAL="min";}
                 "max" {$NEWVAL="max";}
              };
           }
           $ACCURACY {
              label="Accuracy";
              help="If the value in the image/volume is closer than ACCURACY
                 to the old value, then it is substituted";
              type=float [0...];
              by default="0";
           }
      }
   }
   MENU Threshold {
      "I/O Parameters"
      $FILE_IN
      OPT($FILE_OUT)
      "Density conditions"
      OPT(-below)
      OPT(-above)
      OPT(-binarize)
      OPT(-substitute)
      "Distance conditions"
      OPT(-dist)
   }
*/
