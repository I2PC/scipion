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
    Mask            mask;
    MultidimArray<int>   maskArray;
    int             short_format;     // True if a short line is to be shown
    int             save_mask;        // True if the masks must be saved
    int             repair;           // True if headers are initialized
    bool            show_angles;      // True if angles are to be shown in stats
    bool            apply_mask;       // True if a mask must be applied

    double min_val, max_val, avg, stddev;
    double mean_min_val, mean_max_val, mean_avg, mean_stddev;
    int N;
    int max_length;

    void defineParams()
    {
        each_image_produces_an_output = false;
        allow_apply_geo = true;
        allow_time_bar = false;
        addUsageLine("Display statistics of images or volumes. A mask can be applied.");
        XmippMetadataProgram::defineParams();
        addParamsLine("[-o <metadata>]   : Save the statistics in this metadata file.");
        addParamsLine("[--short_format]   : Do not show labels for statistics.");
        addParamsLine("[--show_angles]    : Also show angles in the image header.");
        addParamsLine("[--save_mask]      : Save 2D and 3D masks (as \"mask2D\" or \"mask3D\").");
        mask.defineParams(this,INT_MASK,NULL,"Statistics restricted to the mask area.");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        short_format = checkParam("--short_format");
        save_mask    = checkParam("--save_mask");
        show_angles  = checkParam("--show_angles");
        fn_out = (checkParam("-o"))? getParam("-o"): "";

        mask.allowed_data_types = INT_MASK;
        if (apply_mask = checkParam("--mask"))
            mask.readParams(this);
    }

    void show()
    {
        std::cout << " Statistics of " << fn_in << std::endl;
    }

    void preProcess()
    {
        DF_stats.setComment((std::string)"Statistics of " + fn_in);
        // Get maximum filename size ---------------------------------------------
        max_length = mdIn.getMaxStringLength(MDL_IMAGE);

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

        //        // if input is volume do not apply geo
        //        int xDim, yDim, zDim;
        //        size_t nDim;
        //        ImgSize(mdIn, xDim, yDim, zDim, nDim);
        //
        //        if (zDim > 1)
        //            apply_geo = false;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        if (apply_geo)
            image.readApplyGeo(fnImg, mdIn, objId);
        else
            image.read(fnImg, DATA, ALL_IMAGES, true);

        image().setXmippOrigin();

        int xDim,yDim,zDim;
        double rot, tilt, psi;
        image.getDimensions(xDim,yDim,zDim);
        image.getEulerAngles(rot,tilt,psi);

        // Generate mask if necessary
        if (apply_mask)
        {
            mask.generate_mask(zDim, yDim, xDim);
            maskArray = mask.get_binary_mask();

            computeStats_within_binary_mask(maskArray, image(), min_val, max_val,
                                            avg, stddev);
        }
        else
            image().computeStats(avg, stddev, min_val, max_val);


        // Show information
        std::cout << stringToString(fnImg, max_length + 1);
        if (zDim > 1)
            formatString("%4dx%4dx%4d ", zDim, yDim, xDim);
        else
            formatString("%4dx%4d", yDim, xDim);

        if (!short_format)
        {
            std::cout << formatString("min=%10f max=%10f avg=%10f stddev=%10f",
                                      min_val, max_val, avg, stddev);

            if (show_angles)
                std::cout << formatString("rot=%10f tilt=%10f psi=%10f", rot, tilt, psi);
        }
        else
            std::cout << formatString("%10f %10f %10f %10f", min_val, max_val, avg, stddev);

        size_t id;
        id = DF_stats.addObject();
        DF_stats.setValue(MDL_MIN,min_val,id);
        DF_stats.setValue(MDL_MAX,max_val,id);
        DF_stats.setValue(MDL_AVG,avg,id);
        DF_stats.setValue(MDL_STDDEV,stddev,id);

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
                std::cout << formatString("min=%10f max=%10f avg=%10f stddev=%10f",
                                          min_val, max_val, avg, stddev);
            else
                std::cout << formatString("%10f %10f %10f %10f", mean_min_val, mean_max_val, mean_avg, mean_stddev);
            std::cout << std::endl;
        }

        // Save masks -----------------------------------------------------------
        if (save_mask)
            mask.write_mask("mask");

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
    return program.tryRun();
}
