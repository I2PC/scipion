/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2000)
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
#include <data/image.h>
#include <data/metadata.h>
#include <data/image.h>
#include <data/mask.h>
#include <data/histogram.h>

void Usage();

int main(int argc, char **argv)
{
    Image<double>   image;
    MetaData        SF;
    FileName        fn_in, fn_out, fn_sel;
    Mask_Params     mask_prm(INT_MASK);
    bool            automatic_range, do_normalize, is_first = true;
    double          m, M; // range for histogram
    int             StepsNo;
    Histogram1D     hist, histb;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_sel      = getParameter(argc, argv, "-sel","");
        if (fn_sel=="")
            fn_in = getParameter(argc, argv, "-i");
        fn_out     = getParameter(argc, argv, "-o", "");

        StepsNo = textToInteger(getParameter(argc, argv, "-steps", "100"));
        do_normalize = checkParameter(argc, argv, "-norm");
        int i;
        if ((i = paremeterPosition(argc, argv, "-range")) != -1)
        {
            if (i + 2 >= argc)
                REPORT_ERROR(ERR_ARG_MISSING, "Not enough parameters behind -range\n");
            m = textToFloat(argv[i+1]);
            M = textToFloat(argv[i+2]);
            automatic_range = false;
        }
        else
            automatic_range = true;

        if (fn_sel!="" && automatic_range)
            REPORT_ERROR(ERR_ARG_INCORRECT,"Only use selfile in combination with fixed range!");

        mask_prm.read(argc, argv);
    }
    catch (XmippError Xe)
    {
        std::cout << Xe;
        Usage();
        mask_prm.usage();
        exit(1);
    }

    try
    {
        if (fn_sel!="")
            SF.read(fn_sel);
        else
        {
            SF.addObject();
            SF.setValue(MDL_IMAGE,fn_in);
        }
        SF.firstObject();
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            SF.getValue(MDL_IMAGE,fn_in);
            image.read(fn_in);
            image().setXmippOrigin();
            mask_prm.generate_mask(image());

            // Compute histogram ----------------------------------------------------
            if (automatic_range)
                compute_hist_within_binary_mask(mask_prm.get_binary_mask(), image(),
                                                hist, StepsNo);
            else
                compute_hist_within_binary_mask(mask_prm.get_binary_mask(), image(),
                                                hist, m, M, StepsNo);

            if (is_first)
            {
                is_first = false;
                histb = hist;
            }
            else
                histb += hist;
        }

        if (do_normalize)
            histb /= (histb.sum() * histb.step_size);

        if (fn_out != "")
            histb.write(fn_out);
        else
            std::cout << hist;
    }
    catch (XmippError Xe)
    {
        std::cout << Xe;
    }
    exit(0);
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cout << "histogram [Parameters]:\n"
    << "   -i <File_in> / -sel <selfile> : Input Xmipp Volume, Image or selfile\n"
    << "  [-o <File_out>]              : Text file with histogram\n"
    << "  [-range <m> <M>]             : range for the histogram\n"
    << "                                 by default, it is automatic\n"
    << "  [-steps <N=100>]             : number of subdivisions\n"
    << "  [-norm]                      : normalize histogram \n";
}
