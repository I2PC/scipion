/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (1999)
 *             and modified by Sjors Scheres
 *              "     "     "  Joaquin Oton
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

#include <cstdio>
#include "data/args.h"
#include "data/metadata.h"
#include "data/image_generic.h"
#include "data/mask.h"

/* PROGRAM ----------------------------------------------------------------- */
class ProgStatistics: public XmippMetadataProgram
{
protected:
    MetaData        DF_stats;
    ImageGeneric    image;
    Mask            mask_prm;
    MultidimArray<int>   mask;
    int             short_format;     // True if a short line is to be shown
    int             save_mask;        // True if the masks must be saved
    int             repair;           // True if headers are initialized
    bool            show_angles;      // True if angles are to be shown in stats

    double min_val, max_val, avg, stddev;
    double mean_min_val, mean_max_val, mean_avg, mean_stddev;
    int N;
    int max_length;

    void defineParams()
    {
        each_image_produces_an_output = false;
        apply_geo = true;
        allow_time_bar = false;
        addUsageLine("Display statistics of images or volumes. A mask can be applied.");
        XmippMetadataProgram::defineParams();
        addParamsLine("[-o <metadata>]   : Save the statistics in this metadata file.");
        addParamsLine("[-short_format]   : Do not show labels for statistics.");
        addParamsLine("[-show_angles]    : Also show angles in the image header.");
        addParamsLine("[-save_mask]      : Save 2D and 3D masks (as \"mask2D\" or \"mask3D\").");
        mask_prm.defineParams(this,INT_MASK,NULL,"Statistics constrained to the mask area.");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        short_format = checkParam("-short_format");
        save_mask    = checkParam("-save_mask");
        show_angles  = checkParam("-show_angles");
        fn_out = (checkParam("-o"))? getParam("-o"): "";

        mask_prm.allowed_data_types = INT_MASK;
        if (checkParam("--mask"))
            mask_prm.readParams(this);
    }

    void show()
    {
        std::cout << " Statistics of " << fn_in << std::endl;
    }

    void preProcess()
    {
        DF_stats.setComment((std::string)"Statistics of " + fn_in);
        // Get maximum filename size ---------------------------------------------
        max_length = mdIn.MaxStringLength(MDL_IMAGE);

        // Process each file -----------------------------------------------------
        mean_min_val = 0, mean_max_val = 0, mean_avg = 0, mean_stddev = 0;
        N = 0;

        if (short_format)
        {
            std::cout << "Format: Name ZxYxX min max avg stddev ";
            if (show_angles)
                std::cout << " <rot tilt psi>";
            std::cout << '>' << std::endl;
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        image.read(fnImg,true,-1,true);
        image().setXmippOrigin();

        int xDim,yDim,zDim;
        double rot, tilt, psi;
        image.getDimensions(xDim,yDim,zDim);
        image.getEulerAngles(rot,tilt,psi);

        // Generate mask if necessary
        mask_prm.generate_mask(zDim, yDim, xDim);
        mask = mask_prm.get_binary_mask();

        computeStats_within_binary_mask(mask, image(), min_val, max_val,
                                        avg, stddev);
        // Show information
        std::cout << stringToString(fnImg, max_length + 1);
        if (zDim > 1)
            std::cout << integerToString(zDim, 4, ' ') << 'x'
            << integerToString(yDim, 4, ' ') << 'x'
            << integerToString(xDim, 4, ' ') << ' ';
        else
            std::cout << integerToString(yDim, 4, ' ') << 'x'
            << integerToString(xDim, 4, ' ') << ' ';

        if (!short_format)
        {
            std::cout << "min= "    << floatToString(min_val, 10) << ' '
            << "max= "    << floatToString(max_val, 10) << ' '
            << "avg= "    << floatToString(avg    , 10) << ' '
            << "stddev= " << floatToString(stddev , 10) << ' ';
            if (show_angles)
            {
                std::cout << "rot= "    << floatToString(rot , 10) << ' '
                << "tilt= "   << floatToString(tilt, 10) << ' '
                << "psi= "    << floatToString(psi, 10) << ' ';
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
    }

    void postProcess()
    {
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
            mask_prm.write_mask("mask");

        // Save statistics ------------------------------------------------------
        if (fn_out != "")
            DF_stats.write(fn_out);
    }
}
;// end of class ProgStatistics

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgStatistics program;
    program.read(argc, argv);
    program.tryRun();
    return 0;
}
