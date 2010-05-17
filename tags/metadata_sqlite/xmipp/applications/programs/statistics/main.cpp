/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (1999)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/metadata.h>
#include <data/image.h>
#include <data/mask.h>
#include <cstdio>

void Usage();

int main(int argc, char **argv)
{
    FileName        fn_input, fn_stats;
    MetaData        SF;
    MetaData        DF_stats;
    Image<double>  image;

    Mask_Params     mask_prm(INT_MASK);
    int             short_format;     // True if a short line is to be shown
    int             save_mask;        // True if the masks must be saved
    int             repair;           // True if headers are initialized
    bool            show_angles;      // True if angles are to be shown in stats
    bool            apply_geo;        // True if the header must be taken into account
    MultidimArray<int>   mask;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_input = getParameter(argc, argv, "-i");
        if (fn_input.isMetaData())
        {
            SF.read(fn_input);
        }
        else
        {
            SF.addObject();
            SF.setValue( MDL_IMAGE, fn_input);
            SF.setValue( MDL_ENABLED, 1);
        }

        mask_prm.read(argc, argv);
        fn_stats     = getParameter(argc, argv, "-o", "");
        short_format = checkParameter(argc, argv, "-short_format");
        save_mask    = checkParameter(argc, argv, "-save_mask");
        show_angles  = checkParameter(argc, argv, "-show_angles");
        apply_geo    = !checkParameter(argc, argv, "-dont_apply_geo");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        mask_prm.usage();
        exit(1);
    }

    try
    {
        DF_stats.setComment((std::string)"# Statistics of " + fn_input);
        // Get maximum filename size ---------------------------------------------
        int max_length = SF.MaxStringLength(MDL_IMAGE);

        // Process each file -----------------------------------------------------
        double min_val, max_val, avg, stddev;
        int min_val_int, max_val_int;
        double mean_min_val = 0, mean_max_val = 0, mean_avg = 0, mean_stddev = 0;
        int N = 0;

        if (short_format)
        {
            std::cout << "Format: Name ZxYxX min max avg stddev ";
            if (show_angles)
                std::cout << " <rot tilt psi>";
            std::cout << '>' << std::endl;
        }

        long int ret=SF.firstObject();
        if(ret==MetaData::NO_OBJECTS_STORED)
        {
            std::cerr << "Empty inputFile File\n";
            exit(1);
        }
        FileName file_name;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            SF.getValue(MDL_IMAGE, file_name);
            image.read(file_name);
            image().setXmippOrigin();

            // Generate mask if necessary
            mask_prm.generate_mask(image());
            mask = mask_prm.get_binary_mask();

            computeStats_within_binary_mask(mask, image(), min_val, max_val,
                                            avg, stddev);
            // Show information
            std::cout << stringToString(file_name, max_length + 1);
            if (ZSIZE(image()) > 1)
                std::cout << integerToString(ZSIZE(image()), 4, ' ') << 'x'
                << integerToString(YSIZE(image()), 4, ' ') << 'x'
                << integerToString(XSIZE(image()), 4, ' ') << ' ';
            else
                std::cout << integerToString(YSIZE(image()), 4, ' ') << 'x'
                << integerToString(XSIZE(image()), 4, ' ') << ' ';

            if (!short_format)
            {
                std::cout << "min= "    << floatToString(min_val, 10) << ' '
                << "max= "    << floatToString(max_val, 10) << ' '
                << "avg= "    << floatToString(avg    , 10) << ' '
                << "stddev= " << floatToString(stddev , 10) << ' ';
                if (show_angles)
                {
                    std::cout << "rot= "    << floatToString(image.rot() , 10) << ' '
                    << "tilt= "   << floatToString(image.tilt(), 10) << ' '
                    << "psi= "    << floatToString(image.psi() , 10) << ' ';
                }
            }
            else
            {
                std::cout << floatToString(min_val, 10) << ' '
                << floatToString(max_val, 10) << ' '
                << floatToString(avg    , 10) << ' '
                << floatToString(stddev , 10) << ' ';
            }

            DF_stats.addObject();
            DF_stats.setValue(MDL_MIN,min_val);
            DF_stats.setValue(MDL_MAX,max_val);
            DF_stats.setValue(MDL_AVG,avg);
            DF_stats.setValue(MDL_STDDEV,stddev);

            // Total statistics
            N++;
            mean_min_val += min_val;
            mean_max_val += max_val;
            mean_avg     += avg;
            mean_stddev  += stddev;

            // Finish information .................................................
            std::cout << std::endl;

        } // while

        // Show total statistics ------------------------------------------------
        std::cout << "==================================================\n";
        std::cout << "Total number of images/volumes: " << N << std::endl;
        if (N != 0)
        {
            mean_min_val /= N;
            mean_max_val /= N;
            mean_avg     /= N;
            mean_stddev  /= N;

            std::cout << stringToString(" ", max_length + 13);
            if (!short_format)
                std::cout << "min= "    << floatToString(mean_min_val, 10) << ' '
                << "max= "    << floatToString(mean_max_val, 10) << ' '
                << "avg= "    << floatToString(mean_avg    , 10) << ' '
                << "stddev= " << floatToString(mean_stddev , 10) << ' ';
            else
                std::cout << floatToString(mean_min_val, 10) << ' '
                << floatToString(mean_max_val, 10) << ' '
                << floatToString(mean_avg    , 10) << ' '
                << floatToString(mean_stddev , 10) << ' ';
            std::cout << std::endl;
        }

        // Save masks -----------------------------------------------------------
        if (save_mask)
        {
            mask_prm.write_mask("mask");
        }

        // Save statistics ------------------------------------------------------
        if (fn_stats != "")
            DF_stats.write(fn_stats);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    exit(0);
} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cerr << "Purpose:\n";
    std::cerr << "    Displays statistics of images or volumes \n"
    << "      (possibly restricted to the area within a mask)\n";
    std::cerr << "Usage: statistics " << std::endl
    << "    -i               : Selfile with images/volumes \n"
    << "                        or individual image or volume \n"
    << "   [-o <metadata>]   : save the statistics in this metadata file\n"
    << "   [-dont_apply_geo] : do not apply geo when the image is read\n"
    << "   [-short_format]   : Don't show labels for statistics\n"
    << "   [-show_angles]    : Also show angles in the image header \n"
    << "   [-save_mask]      : save 2D and 3D masks (as \"mask2D\" or \"mask3D\") \n";

}

