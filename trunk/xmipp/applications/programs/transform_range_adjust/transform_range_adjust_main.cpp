/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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

#include <data/xmipp_program.h>
#include <data/mask.h>

class ProgTransRangeAdjust: public XmippMetadataProgram
{
public:
    // Image
    Image<double> Im;
    /// min_val.
    double min_val;

    /// max_val.
    double max_val;

    /// noise in %.
    double sigma;

    /// Mask
    Mask mask_prm;

protected:

    void defineParams()
    {
        each_image_produces_an_output = true;

        addUsageLine("Adjust the grey level value of the pixels to a specified range.");
        addKeywords("mask, normalization");
        addSeeAlsoLine("transform_normalize");
        XmippMetadataProgram::defineParams();
        addParamsLine("   --range <min_val> <max_val> : Output minimum and maximum values.");
        addParamsLine("   alias -r;");
        addParamsLine("   [--noise <sigma=0>] : Variation of the limit range values.");
        mask_prm.defineParams(this,INT_MASK, NULL, "Pixels in the mask area are guaranteed to be into the specified range.");

        addExampleLine("Adjust an image in the range [-1, 1]:", false);
        addExampleLine("xmipp_transform_range_adjust -i image.xmp -o image_adjust.xmp --range -1 1");
        addExampleLine("Adjust images in a stack applying a circular mask:", false);
        addExampleLine("xmipp_transform_range_adjust -i images.stk -o image_adjust.stk --range -1 1 --mask circular -32");
        addExampleLine("Adjust a selection file in the range [-100, 100] with sigma noise 10:", false);
        addExampleLine("xmipp_transform_range_adjust -i selection.sel --oroot images/im_adjust:spi -o selection_adjust.sel --range -100 100 --noise 10");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();

        min_val = getDoubleParam("--range", 0);
        max_val = getDoubleParam("--range", 1);
        sigma   = getDoubleParam("--noise");
        randomize_random_generator();

        if (checkParam("--mask"))
        {
            mask_prm.allowed_data_types = INT_MASK;
            mask_prm.readParams(this);
        }
    }

    void show()
    {
        XmippMetadataProgram::show();
        if (verbose==0)
            return;
        std::cout << "Min:   " << min_val << std::endl
        << "Max:   " << max_val << std::endl
        << "Noise: " << sigma   << std::endl;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
        double amin = rnd_gaus(0, sigma);
        double amax = rnd_gaus(0, sigma);
        double minval = min_val + amin;
        double maxval = max_val + amax;

        Im.read(fnImg);
        MultidimArray<double> &ImP = Im();
        ImP.setXmippOrigin();

        if (mask_prm.type == NO_MASK)
            ImP.rangeAdjust(minval, maxval);
        else
        {
            mask_prm.generate_mask(ImP);
            ImP.rangeAdjust(minval, maxval,mask_prm.get_binary_mask());
        }

        Im.write(fnImgOut);
    }
}
;

int main(int argc, char *argv[])
{
    ProgTransRangeAdjust program;
    program.read(argc, argv);
    program.tryRun();

    return 0;
}
