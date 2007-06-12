/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (1999)
 *             and modified by Sjors Scheres
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

int main(int argc, char **argv)
{
    FileName        fn_input, fn_stats;
    SelFile         SF;
    DocFile         DF_stats;
    ImageXmipp      image;
    VolumeXmippT<double> volume;
    VolumeXmippT<int>    volume_int;

    headerXmipp     header;
    Mask_Params     mask_prm(INT_MASK);
    int             short_format;     // True if a short line is to be shown
    int             save_mask;        // True if the masks must be saved
    int             repair;           // True if headers are initialized
    bool            show_angles;      // True if angles are to be shown in stats
    bool            apply_geo;        // True if the header must be taken into account
    Matrix2D<int>   mask2D;
    Matrix3D<int>   mask3D;
    int             volume_type;
    // Read arguments --------------------------------------------------------
    try
    {
        fn_input = getParameter(argc, argv, "-i");
        if (Is_VolumeXmipp(fn_input) || Is_ImageXmipp(fn_input))
        {
            SF.insert(fn_input, SelLine::ACTIVE);
        }
        else
            SF.read(fn_input);

        mask_prm.read(argc, argv);
        fn_stats     = getParameter(argc, argv, "-o", "");
        short_format = checkParameter(argc, argv, "-short_format");
        save_mask    = checkParameter(argc, argv, "-save_mask");
        show_angles  = checkParameter(argc, argv, "-show_angles");
        apply_geo    = !checkParameter(argc, argv, "-dont_apply_geo");
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
        DF_stats.append_comment((string)"# Statistics of " + fn_input);
        DF_stats.append_comment("# min max avg stddev");

        // Get maximum filename size ---------------------------------------------
        int max_length = SF.MaxFileNameLength();

        // Process each file -----------------------------------------------------
#define V VOLMATRIX(volume)
#define VI VOLMATRIX(volume_int)
#define I IMGMATRIX(image)
        double min_val, max_val, avg, stddev;
        int min_val_int, max_val_int;
        double mean_min_val = 0, mean_max_val = 0, mean_avg = 0, mean_stddev = 0;
        int N = 0;

        if (short_format)
        {
            cout << "Format: Name ZxYxX min max avg stddev ";
            if (show_angles) cout << " <rot tilt psi>";
            cout << '>' << endl;
        }

        SF.go_beginning();
        while (!SF.eof())
        {
            FileName file_name = SF.NextImg();

            // For volumes ........................................................
            if ((volume_type = Is_VolumeXmipp(file_name)))
            {
                // Read file
                if (volume_type == headerXmipp::VOL_XMIPP)
                {
                    volume.read(file_name);
                    volume().setXmippOrigin();
                }
                else if (volume_type == headerXmipp::VOL_INT)
                {
                    volume_int.read(file_name);
                    volume_int().setXmippOrigin();
                }

                // Generate mask if necessary
                mask_prm.generate_3Dmask(V);
                const Matrix3D<int> &mask3D = mask_prm.get_binary_mask3D();

                if (volume_type == headerXmipp::VOL_XMIPP)
                    computeStats_within_binary_mask(mask3D, V, min_val, max_val,
                                                     avg, stddev);
                else if (volume_type == headerXmipp::VOL_INT)
                    computeStats_within_binary_mask(mask3D, VI, min_val_int,
                                                     max_val_int, avg, stddev);

                // Show information
                cout << stringToString(file_name, max_length + 1);
                cout << integerToString(ZSIZE(V), 4, ' ') << 'x'
                << integerToString(YSIZE(V), 4, ' ') << 'x'
                << integerToString(XSIZE(V), 4, ' ') << ' ';
                if (!short_format)
                    cout << "min= "    << floatToString(min_val, 10) << ' '
                    << "max= "    << floatToString(max_val, 10) << ' '
                    << "avg= "    << floatToString(avg    , 10) << ' '
                    << "stddev= " << floatToString(stddev , 10) << ' ';
                else
                    cout << floatToString(min_val, 10) << ' '
                    << floatToString(max_val, 10) << ' '
                    << floatToString(avg    , 10) << ' '
                    << floatToString(stddev , 10) << ' ';
                Matrix1D<double> v(4);
                v(0) = min_val;
                v(1) = max_val;
                v(2) = avg;
                v(3) = stddev;
                DF_stats.append_data_line(v);

                // Total statistics
                N++;
                mean_min_val += min_val;
                mean_max_val += max_val;
                mean_avg     += avg;
                mean_stddev  += stddev;

                // For images .........................................................
            }
            else if (Is_ImageXmipp(file_name))
            {
                // Read the image applying the header
                image.read(file_name, false, false, apply_geo);
                image().setXmippOrigin();

                // Generate mask if necessary
                mask_prm.generate_2Dmask(I);
                const Matrix2D<int> &mask2D = mask_prm.get_binary_mask2D();

                // Compute statistics
                computeStats_within_binary_mask(mask2D, I, min_val, max_val,
                                                 avg, stddev);

                // Show information
                cout << stringToString(file_name, max_length + 1);
                cout << "    "; // Stands for the ZSIZE in volumes
                cout << integerToString(YSIZE(I), 4, ' ') << 'x'
                << integerToString(XSIZE(I), 4, ' ') << ' ';
                if (!short_format)
                {
                    cout << "min= "    << floatToString(min_val, 10) << ' '
                    << "max= "    << floatToString(max_val, 10) << ' '
                    << "avg= "    << floatToString(avg    , 10) << ' '
                    << "stddev= " << floatToString(stddev , 10) << ' ';
                    if (show_angles)
                    {
                        cout << "rot= "    << floatToString(image.rot() , 10) << ' '
                        << "tilt= "   << floatToString(image.tilt(), 10) << ' '
                        << "psi= "    << floatToString(image.psi() , 10) << ' ';
                        if (image.Is_flag_set() == 1.0f || image.Is_flag_set() == 2.0f)
                            cout << "\nrot1= "  << floatToString(image.rot1() , 10) << ' '
                            << "tilt1= "   << floatToString(image.tilt1(), 10) << ' '
                            << "psi1= "    << floatToString(image.psi1() , 10) << ' ';
                        if (image.Is_flag_set() == 2.0f)
                            cout << "\nrot2= "    << floatToString(image.rot2() , 10) << ' '
                            << "tilt2= "   << floatToString(image.tilt2(), 10) << ' '
                            << "psi2= "    << floatToString(image.psi2() , 10) << ' ';
                    }

                }
                else
                {
                    cout << floatToString(min_val, 10) << ' '
                    << floatToString(max_val, 10) << ' '
                    << floatToString(avg    , 10) << ' '
                    << floatToString(stddev , 10) << ' ';
                    if (show_angles)
                    {
                        cout << floatToString(image.rot() , 10) << ' '
                        << floatToString(image.tilt(), 10) << ' '
                        << floatToString(image.psi() , 10) << ' ';
                        if (image.Is_flag_set() == 1.0f || image.Is_flag_set() == 2.0f)
                            cout << floatToString(image.rot1() , 10) << ' '
                            << floatToString(image.tilt1(), 10) << ' '
                            << floatToString(image.psi1() , 10) << ' ';
                        if (image.Is_flag_set() == 2.0f)
                            cout << floatToString(image.rot2() , 10) << ' '
                            << floatToString(image.tilt2(), 10) << ' '
                            << floatToString(image.psi2() , 10) << ' ';
                    }
                }

                Matrix1D<double> v(4);
                v(0) = min_val;
                v(1) = max_val;
                v(2) = avg;
                v(3) = stddev;
                DF_stats.append_data_line(v);

                // Total statistics
                N++;
                mean_min_val += min_val;
                mean_max_val += max_val;
                mean_avg     += avg;
                mean_stddev  += stddev;

                // Is not an Spider file ..............................................
            }
            else
                cout << file_name << " is not a Spider image nor a volume... ";

            // Finish information .................................................
            cout << endl;

        } // while

        // Show total statistics ------------------------------------------------
        cout << "==================================================\n";
        cout << "Total number of images/volumes: " << N << endl;
        if (N != 0)
        {
            mean_min_val /= N;
            mean_max_val /= N;
            mean_avg     /= N;
            mean_stddev  /= N;

            cout << stringToString(" ", max_length + 13);
            if (!short_format)
                cout << "min= "    << floatToString(mean_min_val, 10) << ' '
                << "max= "    << floatToString(mean_max_val, 10) << ' '
                << "avg= "    << floatToString(mean_avg    , 10) << ' '
                << "stddev= " << floatToString(mean_stddev , 10) << ' ';
            else
                cout << floatToString(mean_min_val, 10) << ' '
                << floatToString(mean_max_val, 10) << ' '
                << floatToString(mean_avg    , 10) << ' '
                << floatToString(mean_stddev , 10) << ' ';
            cout << endl;
        }

        // Save masks -----------------------------------------------------------
        if (save_mask)
        {
            mask_prm.write_2Dmask("mask2D");
            mask_prm.write_3Dmask("mask3D");
        }

        // Save statistics ------------------------------------------------------
        if (fn_stats != "") DF_stats.write(fn_stats);
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
    cerr << "Purpose:\n";
    cerr << "    Displays statistics of images or volumes \n"
    << "      (possibly restricted to the area within a mask)\n";
    cerr << "Usage: statistics " << endl
    << "    -i               : Selfile with images/volumes \n"
    << "                        or individual image or volume \n"
    << "   [-o <docfile>]    : save the statistics in this docfile\n"
    << "   [-dont_apply_geo] : do not apply geo when the image is read\n"
    << "   [-short_format]   : Don't show labels for statistics\n"
    << "   [-show_angles]    : Also show angles in the image header \n"
    << "   [-save_mask]      : save 2D and 3D masks (as \"mask2D\" or \"mask3D\") \n";

}


