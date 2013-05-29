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
#include <data/xmipp_program.h>
#include <data/xmipp_image_generic.h>
#include <data/metadata.h>
#include <data/mask.h>

/* PROGRAM ----------------------------------------------------------------- */
class ProgStatistics: public XmippMetadataProgram
{
protected:
    MetaData        DF_stats;
    ImageGeneric    image;
    MultidimArray<double>    averageArray;
    MultidimArray<double>    stdArray;
    MultidimArray<double>    dummyArray;

    Mask            mask;
    int             short_format;     // True if a short line is to be shown
    int             save_mask;        // True if the masks must be saved
    int             repair;           // True if headers are initialized
    bool            show_angles;      // True if angles are to be shown in stats
    bool            save_image_stats; // True if average and std images must be computed
    bool            apply_mask;       // True if a mask must be applied

    double min_val, max_val, avg, stddev;
    double mean_min_val, mean_max_val, mean_avg, mean_stddev;
    int max_length;


    FileName maskFileName, statsRoot;

    void defineParams()
    {
        each_image_produces_an_output = false;
        allow_apply_geo = true;
        allow_time_bar = false;
        addUsageLine("Display statistics of 2D/3Dimages. A mask can be applied. Average Images may be computed");
        addUsageLine("All images must have same size");
        XmippMetadataProgram::defineParams();
        addParamsLine("[-o <metadata>]   : Save the statistics in this metadata file.");
        addParamsLine("[--short_format]  : Do not show labels for statistics.");
        addParamsLine("[--show_angles]   : Also show angles in the image header.");
        addParamsLine("[--save_mask        <maskFileName>] : Save 2D and 3D masks.");
        addParamsLine("[--save_image_stats <stats_root=\"\">]: Save average and standard deviation images");
        addUsageLine ("           Mask is ignored  for this operation");
        mask.defineParams(this,INT_MASK,NULL,"Statistics restricted to the mask area.");
        addUsageLine ("NOTE: Geometry will NOT be applied to volumes even if apply_geo flag is on");
        addKeywords("statistics average mean std");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        short_format = checkParam("--short_format");

        if ((save_mask = checkParam("--save_mask")))
            maskFileName = getParam("--save_mask");

        if ((save_image_stats = checkParam("--save_image_stats")))
            statsRoot = getParam("--save_image_stats");

        show_angles  = checkParam("--show_angles");
        fn_out = (checkParam("-o"))? getParam("-o"): "";

        mask.allowed_data_types = INT_MASK;

        if ((apply_mask = checkParam("--mask")))
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
        max_length = getInputMd()->getMaxStringLength(image_label);

        // Process each file -----------------------------------------------------
        mean_min_val = 0, mean_max_val = 0, mean_avg = 0, mean_stddev = 0;

        if (short_format)
        {
            std::cout << "Format: Name ZxYxX min max avg stddev ";
            if (show_angles)
                std::cout << " <rot tilt psi>";
            std::cout << '>' << std::endl;
        }

        // get xdim, ydim,zdim
        //getImageSize(mdIn, xDim, yDim, zDim, nDim, image_label);
        averageArray.resize(ndimOut, zdimOut, ydimOut, xdimOut);
        stdArray.resize(ndimOut,zdimOut,ydimOut,xdimOut);
        averageArray.setXmippOrigin();
        stdArray.setXmippOrigin();

        // Generate mask if necessary
        if (apply_mask)
            mask.generate_mask(zdimOut, ydimOut, xdimOut);

    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
        if (apply_geo)
            image.readApplyGeo(fnImg, rowIn);
        else
            image.read(fnImg, DATA, ALL_IMAGES, true);

        image().setXmippOrigin();

        double rot, tilt, psi;
        if (show_angles)
            image.getEulerAngles(rot,tilt,psi);

        if (apply_mask)
        {
            computeStats_within_binary_mask(mask.get_binary_mask(), image(), min_val, max_val,
                                            avg, stddev);
        }
        else
            image().computeStats(avg, stddev, min_val, max_val);

        if(save_image_stats)
        {
            //copy image from imageGeneric
            image().getImage( dummyArray );
            averageArray     += dummyArray;
            stdArray         += dummyArray * dummyArray;
        }

        // Show information
        std::cout << stringToString(fnImg, max_length + 1);
        if (zdimOut > 1)
            formatString("%4dx%4dx%4d ", zdimOut, ydimOut, xdimOut);
        else
            formatString("%4dx%4d", ydimOut, xdimOut);

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
        DF_stats.setRow(rowIn,id);
        DF_stats.setValue(MDL_MIN,min_val,id);
        DF_stats.setValue(MDL_MAX,max_val,id);
        DF_stats.setValue(MDL_AVG,avg,id);
        DF_stats.setValue(MDL_STDDEV,stddev,id);

        // Total statistics
        mean_min_val += min_val;
        mean_max_val += max_val;
        mean_avg     += avg;
        mean_stddev  += stddev;

        // Finish information .................................................
        std::cout << std::endl;
    }

    void postProcess()
    {

        if (mdInSize > 1)
        {
            // Show total statistics ------------------------------------------------
            std::cout << "==================================================\n";
            std::cout << "Total number of images/volumes: " << mdInSize << std::endl;
            mean_min_val /= mdInSize;
            mean_max_val /= mdInSize;
            mean_avg     /= mdInSize;
            mean_stddev  /= mdInSize;

            std::cout << stringToString(" ", max_length + 13);
            if (!short_format)
                std::cout << formatString("min=%10f max=%10f avg=%10f stddev=%10f",
                                          mean_min_val, mean_max_val, mean_avg, mean_stddev);
            else
                std::cout << formatString("%10f %10f %10f %10f", mean_min_val, mean_max_val, mean_avg, mean_stddev);
            std::cout << std::endl;
            if(save_image_stats)
            {
                averageArray /= mdInSize;
                if (mdInSize > 1)
                {
                    stdArray= stdArray / mdInSize - averageArray * averageArray;
                    stdArray *= mdInSize / (mdInSize - 1);
                    //Do this as an image since it is not define for arrays
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(stdArray)
                    DIRECT_MULTIDIM_ELEM(stdArray,n) = sqrt(fabs(DIRECT_MULTIDIM_ELEM(stdArray,n)));
                }
                else
                    stdArray.initZeros();
            }
        }

        // Save masks -----------------------------------------------------------
        if (save_mask)
            mask.write_mask(maskFileName);
        // Save statistics ------------------------------------------------------
        if (fn_out != "")
            DF_stats.write(fn_out,MD_APPEND);

        //save average and std images
        if(save_image_stats)
        {
            Image<double> dummyImage;
            dummyImage()=averageArray;
            dummyImage.write(statsRoot  + "average.xmp");
            dummyImage()=stdArray;
            dummyImage.write(statsRoot + "stddev.xmp");
        }
    }
};// end of class ProgStatistics

RUN_XMIPP_PROGRAM(ProgStatistics)

