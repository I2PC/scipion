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

#include <data/xmipp_image.h>
#include <data/metadata.h>
#include <data/mask.h>
#include <data/histogram.h>

class ProgHistogram: public XmippProgram
{
public:
    Image<double>   image;
    MetaData        mdIn;
    FileName        fn_in, fn_out, fn_sel;
    Mask            mask;
    MultidimArray<int>   maskArray;
    bool            automatic_range, do_normalize, apply_mask;
    double          m, M; // range for histogram
    int             StepsNo;
    Histogram1D     hist, histb;

    void defineParams()
    {
        addUsageLine("Calculate histogram from a list of images/volumes");
        addParamsLine("   -i <input_file>              : input metadata, image, volume or stack");
        addParamsLine("  [-o <text_file>]              : output text file with histogram");
        addParamsLine("  [--range <m> <M>]             : range for the histogram, automatic calculated if not provided");
        addParamsLine("  [--steps <N=100>]             : number of subdivisions");
        addParamsLine("  [--norm]                      : normalize histogram ");
        mask.defineParams(this, INT_MASK, NULL, "Histogram constrained to the mask area.");

        addExampleLine("Create the histogram from an stack setting number of steps to 20", false);
        addExampleLine("xmipp_histogram -i input/images_some.stk -o hist.txt --steps 20");
        addExampleLine("You can use the =hist.txt= file with =gnuplot= to produce graphical histogram", false);
        addExampleLine("After you enter =gnuplot= terminal you can use:", false);
        addExampleLine("gnuplot> plot \"kk.txt\" using 1:2 title \"Histogram\" with steps");

        addSeeAlsoLine("metadata_histogram");
    }

    void readParams()
    {
        fn_in = getParam("-i");
        mdIn.read(fn_in);

        if (mdIn.isEmpty())
            REPORT_ERROR(ERR_PARAM_INCORRECT, "Can't calculate histogram from empty metadata");

        if (checkParam("-o"))
            fn_out = getParam("-o");

        StepsNo = getIntParam("--steps");
        do_normalize = checkParam("--norm");
        automatic_range = true;

        if (checkParam("--range"))
        {
            m = getDoubleParam("--range", 0);
            M = getDoubleParam("--range", 1);
            automatic_range = false;
        }

        mask.allowed_data_types = INT_MASK;
        if ((apply_mask = checkParam("--mask")))
            mask.readParams(this);

        //todo: check this
        //        if (fn_sel!="" && automatic_range)
        //            REPORT_ERROR(ERR_ARG_INCORRECT,"Only use selfile in combination with fixed range!");
        //
        //        mask_prm.read(argc, argv);

    }

    void run()
    {
        MDIterator iter(mdIn);
        image.readApplyGeo(mdIn, iter.objId);
        image().setXmippOrigin();
        double dummy;

        // Generate mask if necessary
        if (apply_mask)
        {
            mask.generate_mask(image());
            maskArray = mask.get_binary_mask();
            if (automatic_range)
                computeStats_within_binary_mask(maskArray, image(), m, M, dummy, dummy);
        }
        else if (automatic_range)
            image().computeDoubleMinMax(m, M);

        if (automatic_range)
        {
            m *= 1.5;
            M *= 1.5;
        }

        if (apply_mask)
          compute_hist_within_binary_mask(maskArray, image(), histb, m, M, StepsNo);
        else
          compute_hist(image(), histb, m, M, StepsNo);

        while (iter.moveNext())
        {
            image.readApplyGeo(mdIn, iter.objId);
            image().setXmippOrigin();
            // Compute histogram ----------------------------------------------------
            if (apply_mask)
              compute_hist_within_binary_mask(maskArray, image(), hist, m, M, StepsNo);
            else
              compute_hist(image(), hist, m, M, StepsNo);
            histb += hist;
        }

        if (do_normalize)
            histb *= 1.0 / (histb.sum() * histb.step_size);

        if (!fn_out.empty())
            histb.write(fn_out);
        else
            std::cout << histb;
    }

}
;//end of class ProgHistogram


void Usage();

int main(int argc, char **argv)
{
    ProgHistogram program;
    program.read(argc, argv);
    return program.tryRun();
}
