/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (1999)
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

#include <data/histogram.h>
#include <data/xmipp_program.h>

class ProgMetadataHistogram: public XmippProgram
{
public:
    MetaData         mdIn;
    FileName         fn_in, fn_out, fn_img;
    MDLabel          col, col2;         // Columns for histogram
    double           m, M, m2, M2;      // range for histogram
    bool             do_hist2d, write_img, automatic_range, automatic_range2;
    int              StepsNo, StepsNo2;
    Histogram1D      hist;
    Histogram2D      hist2;
    double           percentil;

    void defineParams()
    {
        addUsageLine("Calculate histogram from a metadata column(1D) or from a couple of columns(2D)");
        addParamsLine("   -i <input_metadata>          : input metadata");
        addParamsLine("  [-o <text_file=\"/dev/stdout\">] : output text file with histogram, by default standard output");
        addParamsLine("  --col <label>                 : column to create the histogram");
        addParamsLine("  [--range <m> <M>]             : range for the histogram, automatic calculated if not provided");
        addParamsLine("  [--steps <N=100>]             : number of subdivisions");
        addParamsLine("  [--col2 <label=\"\"> ]        : if specified, a 2D histogram is calculated");
        addParamsLine("  [--range2 <m> <M>]            : range for second column in 2D histogram");
        addParamsLine("       requires --col2;");
        addParamsLine("  [--steps2 <N=100>]            : number of subdivisions in second column");
        addParamsLine("       requires --col2;");
        addParamsLine("  [--percentil <p=50.>]         : Only for 1D histograms");
        addParamsLine("  [--write_as_image <image_file>]   : Only for 2D histograms");
        addParamsLine("       requires --col2;");

        addExampleLine("Calculate the histogram of the column 'angleRot' from a metadata:", false);
        addExampleLine("xmipp_metadata_histogram -i images.xmd --col angleRot -o hist.txt");

        addSeeAlsoLine("image_histogram");

    }

    void readColumn(MDLabel &column, bool &automatic, double &m, double &M, int &steps, bool second=false)
    {
        const char* colStr = second ? "--col2" : "--col";
        const char* rangeStr = second ? "--range2" : "--range";
        const char* stepStr = second ? "--steps2" : "--steps";

        column = MDL::str2Label(getParam(colStr));
        steps = getIntParam(stepStr);
        automatic = true;
        if (checkParam(rangeStr))
        {
            m = getDoubleParam(rangeStr, 0);
            M = getDoubleParam(rangeStr, 1);
            automatic = false;
        }

        // Check if valid columns where provided
        if (column == MDL_UNDEFINED)
            REPORT_ERROR(ERR_MD_UNDEFINED, "Metadata Histogram: Column for histogram not valid");
        else if (!MDL::isDouble(column))
            REPORT_ERROR(ERR_MD_BADTYPE, "Metadata Histogram: Column type for histogram should be double");
    }

    void getColumnValues(const MDLabel column, MultidimArray<double> &values, const bool automatic, double &m, double &M)
    {
        std::vector<double> columnValues;
        mdIn.getColumnValues(column, columnValues);
        values = columnValues;
        if (automatic)
           values.computeDoubleMinMax(m, M);
    }

    void readParams()
    {
        fn_in = getParam("-i");
        mdIn.read(fn_in);
        fn_out = getParam("-o");
        percentil = getDoubleParam("--percentil");
        readColumn(col, automatic_range, m, M, StepsNo);

        if ((do_hist2d = checkParam("--col2")))
            readColumn(col2, automatic_range2, m2, M2, StepsNo2, true);

        if ((write_img = checkParam("--write_as_image")))
            fn_img = getParam("--write_as_image");
    }

    void run()
    {
        double avg=0., stddev=0., dummy;
        MultidimArray<double> C;
        getColumnValues(col, C, automatic_range, m, M);

        if (!do_hist2d)
        {
            compute_hist(C, hist, m, M, StepsNo);
            std::cout << formatString("min: %f max: %f steps: %d", m, M, StepsNo) << std::endl;
            C.computeStats(avg, stddev, dummy, dummy);
            std::cout << formatString("mean: %f stddev: %f", avg, stddev) << std::endl;
            std::cout << formatString("percentil (%f)", hist.percentil(percentil)) << std::endl;
            hist.write(fn_out);
        }
        else
        {
            MultidimArray<double> C2;
            getColumnValues(col2, C2, automatic_range2, m2, M2);
            compute_hist(C, C2, hist2, m, M, m2, M2, StepsNo, StepsNo2);
            //stats for column 1
            std::cout << formatString("min1: %f max1: %f steps1: %d", m, M, StepsNo) << std::endl;
            C.computeStats(avg, stddev, dummy, dummy);
            std::cout << formatString("mean: %f stddev: %f", avg, stddev) << std::endl;
            //stats for column 2
            std::cout << formatString("min2: %f max2: %f steps2: %d", m2, M2, StepsNo2) << std::endl;
            C2.computeStats(avg, stddev, dummy, dummy);
            std::cout << formatString("mean: %f stddev: %f", avg, stddev) << std::endl;
            hist2.write(fn_out);

            if (write_img)
            {
              Image<double> img;
              img() = hist2;
              img.write(fn_img);
            }
        }
    }

}
;//end of class ProgHistogram

int main(int argc, char *argv[])
{
  ProgMetadataHistogram program;
  program.read(argc, argv);
  return program.tryRun();
}
