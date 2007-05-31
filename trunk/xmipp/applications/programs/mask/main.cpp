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

#include <data/mask.h>
#include <data/matrix2d.h>
#include <data/args.h>
#include <data/selfile.h>
#include <data/volume.h>

#include <cstdio>

void Usage();

#define COUNT_ELEMENTS(mask,m,elem_type) \
    { \
        int count; \
        if      (count_above && !count_below) {\
            cout << AtoA(fn_in,max_length) \
            << " number of " << elem_type << " above " << th_above; \
            count=count_with_mask_above(mask,m,th_above); \
        } else if (count_below && !count_above) {\
            cout << AtoA(fn_in,max_length) \
            << " number of " << elem_type << " below " << th_below; \
            count=count_with_mask_below(mask,m,th_below); \
        } else if (count_below && count_above) {\
            cout << AtoA(fn_in,max_length) \
            << " number of " << elem_type << " above " << th_above \
            << " and below " << th_below << " = "; \
            count=count_with_mask_between(mask,m,th_above,th_below); \
        } \
        cout << " = " << count << endl; \
    }

#define SET_SUBS_VAL(I) \
    if      (str_subs_val=="min") subs_val=I.compute_min(); \
    else if (str_subs_val=="max") subs_val=I.compute_max(); \
    else if (str_subs_val=="avg") subs_val=I.compute_avg(); \
    else                          subs_val=AtoF(str_subs_val);

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv)
{
    FileName        oext, fn_out, fn_in, fn_input, fn_mask;
    SelFile         SF_in, SF_out;
    ImageXmipp      image;
    VolumeXmipp     volume;
    Mask_Params     mask_prm;
    Matrix2D<int>   mask2D;
    Matrix3D<int>   mask3D;
    int             save_mask;
    int             create_mask;
    int             count_above;
    double th_above;
    int             count_below;
    double th_below;
    bool            apply_geo;
    double          subs_val;
    string          str_subs_val;
    int             count;
    int             max_length = 0;
#define V VOLMATRIX(volume)
#define I IMGMATRIX(image)

    // Read arguments --------------------------------------------------------
    try
    {
        fn_input     = getParameter(argc, argv, "-i", NULL, 1, "Mask: Input file not found");
        fn_out       = getParameter(argc, argv, "-o", "");
        oext         = getParameter(argc, argv, "-oext", "");
        save_mask    = checkParameter(argc, argv, "-save_mask");
        count_above  = checkParameter(argc, argv, "-count_above");
        apply_geo    = !checkParameter(argc, argv, "-dont_apply_geo");
        if (count_above)
            th_above  = AtoF(getParameter(argc, argv, "-count_above"));
        count_below  = checkParameter(argc, argv, "-count_below");
        if (count_below)
            th_below  = AtoF(getParameter(argc, argv, "-count_below"));
        create_mask  = checkParameter(argc, argv, "-create_mask");
        if (create_mask)
        {
            fn_mask  = getParameter(argc, argv, "-create_mask");
            if (fn_mask[0] == '-')
                REPORT_ERROR(1, "Mask: the output mask filename is missing");
        }
        mask_prm.read(argc, argv);
        str_subs_val = getParameter(argc, argv, "-substitute", "0");

        count = count_below || count_above;
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage();
        mask_prm.usage();
        exit(1);
    }

    try
    {
        // Mask a single image ---------------------------------------------------
        if (Is_ImageXmipp(fn_input))
        {
            image.read(fn_input);
            image().setXmippOrigin();
            if (apply_geo)
            {
                if (mask_prm.x0 + mask_prm.y0 != 0.)
                {
                    REPORT_ERROR(1, "Mask: -center option cannot be combined with apply_geo; use -dont_apply_geo");
                }
                else
                {
                    // Read geometric transformation from the image and store for mask
                    mask_prm.mask_geo = image.get_transformation_matrix();
                }
            }
            mask_prm.generate_2Dmask(I);
            if (!create_mask)
            {
                SET_SUBS_VAL(I);
                mask_prm.apply_mask(I, I, subs_val, apply_geo);
                if (!count)
                    if (fn_out == "") image.write(fn_input);
                    else            image.write(fn_out);
            }
            else mask_prm.write_2Dmask(fn_mask);
            max_length = fn_input.length();
            fn_in = fn_input;
            if (count)
            {
                if (mask_prm.datatype() == INT_MASK)
                {
                    const Matrix2D<int> &mask2D = mask_prm.get_binary_mask2D();
                    COUNT_ELEMENTS(mask2D, I, "pixels");
                }
                else
                    cerr << "Cannot count pixels with a continuous mask\n";
            }

            // Mask a single volume --------------------------------------------------
        }
        else if (Is_VolumeXmipp(fn_input))
        {
            volume.read(fn_input);
            volume().setXmippOrigin();
            mask_prm.generate_3Dmask(V);
            if (!create_mask)
            {
                SET_SUBS_VAL(V);
                mask_prm.apply_mask(V, V, subs_val);
                if (!count)
                    if (fn_out == "") volume.write(fn_input);
                    else            volume.write(fn_out);
            }
            else mask_prm.write_3Dmask(fn_mask);
            max_length = fn_input.length();
            fn_in = fn_input;
            if (count)
            {
                if (mask_prm.datatype() == INT_MASK)
                {
                    const Matrix3D<int> &mask3D = mask_prm.get_binary_mask3D();
                    COUNT_ELEMENTS(mask3D, V, "voxels");
                }
                else
                    cerr << "Cannot count pixels with a continuous mask\n";
            }

        }
        else
        {
            // Mask a selection file ------------------------------------------------
            if (create_mask)
                EXIT_ERROR(1, "Mask: Cannot create a mask for a non-Spider file\n");

            SF_in.read(fn_input);
            SF_out.clear();

            // Initialise progress bar
            time_config();
            int i = 0;
            if (!count) init_progress_bar(SF_in.ImgNo());

            // Get maximum filename size
            if (count) max_length = SF_in.MaxFileNameLength();

            // Process all selfile
            while (!SF_in.eof())
            {
                // In and out filenames ...........................................
                fn_in = SF_in.NextImg();

                if (oext == "") fn_out = fn_in;
                else
                {
                    fn_out = fn_in.without_extension() + "." + oext;
                    SF_out.insert(fn_out);
                }

                // Process an image ...............................................
                if (Is_ImageXmipp(fn_in))
                {
                    image.read(fn_in);
                    image().setXmippOrigin();
                    if (apply_geo)
                    {
                        if (mask_prm.x0 + mask_prm.y0 != 0.)
                        {
                            REPORT_ERROR(1, "Mask: -center option cannot be combined with apply_geo; use -dont_apply_geo");
                        }
                        else
                        {
                            // Read geometric transformation from the image and store for mask
                            mask_prm.mask_geo = image.get_transformation_matrix();
                        }
                    }
                    mask_prm.generate_2Dmask(I);
                    SET_SUBS_VAL(I);
                    mask_prm.apply_mask(I, I, subs_val, apply_geo);
                    if (!count) image.write(fn_out);
                    else
                    {
                        if (mask_prm.datatype() == INT_MASK)
                        {
                            const Matrix2D<int> &mask2D = mask_prm.get_binary_mask2D();
                            COUNT_ELEMENTS(mask2D, I, "pixels");
                        }
                        else
                            cerr << "Cannot count pixels with a continuous mask\n";
                    }

                    // Process a volume ...............................................
                }
                else if (Is_VolumeXmipp(fn_in))
                {
                    volume.read(fn_in);
                    volume().setXmippOrigin();
                    mask_prm.generate_3Dmask(V);
                    SET_SUBS_VAL(V);
                    mask_prm.apply_mask(V, V, subs_val);
                    if (!count) volume.write(fn_out);
                    else
                    {
                        if (mask_prm.datatype() == INT_MASK)
                        {
                            const Matrix3D<int> &mask3D = mask_prm.get_binary_mask3D();
                            COUNT_ELEMENTS(mask3D, V, "voxels");
                        }
                        else
                            cerr << "Cannot count pixels with a continuous mask\n";
                    }

                    // Not a Spider file ..............................................
                }
                else
                    cout << fn_in << " is not a SPIDER file\n";

                if (i++ % 25 == 0 && !count) progress_bar(i);
            }
            if (!count) progress_bar(SF_in.ImgNo());
            if (oext != "")
            {
                fn_out = fn_input.insert_before_extension(oext);
                SF_out.write(fn_out);
            }
        }

        // Save masks -----------------------------------------------------------
        if (save_mask)
        {
            mask_prm.write_2Dmask("mask2D");
            mask_prm.write_3Dmask("mask3D");
        }

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }
    exit(0);
} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    cerr << "Usage: mask <parameters>\n"
    << "   -i <image or volume> [-o <image_out or volume_out]\n"
    << "   -i <selfile> [-oext <output extension>]\n"
    << "   [-save_mask]                       : apply and save mask\n"
    << "   [-dont_apply_geo]                  : dont apply (opposite) geometric transformation as stored in header of 2D-image\n"
    << "   [-create_mask <output mask file>]  : don't apply and save mask\n"
    << "   [-count_above <th>]                : voxels within mask >= th\n"
    << "   [-count_below <th>]                : voxels within mask <= th\n"
    << "   [-substitute <val=0>|min|max|avg]  : value outside the mask\n";
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Mask {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Mask/Help/mask.html";
      help="Apply masks, create mask, create blank images and count voxels";
      OPEN MENU menu_mask;
      COMMAND LINES {
 + apply_single_mask: xmipp_mask
               #include "prog_line.mnu"
               #include "binary_mask_line.mnu"
              [-substitute $VAL $MINMAX]
        + create_mask: xmipp_mask -i $FILEIN -create_mask
               OPT($BMASK_TYPE)
        + create_blank : xmipp_mask -i $FILEIN -create_mask
        + count_voxels: xmipp_mask -i $FILEIN
               [-count_above $ABOVE_TH] [-count_below $BELOW_TH]
        + count_voxels_within_mask: xmipp_mask -i $FILEIN
               OPT($BMASK_TYPE)
               OPT(-count_above) OPT(-count_below)
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        #include "binary_mask_vars.mnu"
        $VAL {
           label="Substitute value";
           help="For values outside the mask";
           type=text;
           by default="0";
        }
           $MINMAX {
              label="Value type";
              type=list {
                 "User Defined"
                 "min" {$VAL="min";}
                 "max" {$VAL="max";}
                 "avg" {$VAL="avg";}
              };
           }
        $ABOVE_TH {label="Count above"; type=float;}
        $BELOW_TH {label="Count below"; type=float;}
      }
   }

   MENU menu_mask {
      #include "prog_menu.mnu"
      "Other parameters"
      OPT(-count_above)
      OPT(-count_below)
      OPT(-substitute)
      #include "binary_mask_menu.mnu"
   }
*/
