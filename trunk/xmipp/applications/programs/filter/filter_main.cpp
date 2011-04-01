/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include "data/program.h"
#include <data/filters.h>


struct FilterParam
{
    Image<char> *mask;
    double     factor;
};

//This define the prototype of filter operations on images
//the result will be left in op
typedef void ImageFilterOperator(Image<double> &op, FilterParam &params);


void forcePos(Image<double> &op, FilterParam &params)
{
    forcePositive(op());
}

void removeBadPixels(Image<double> &op, FilterParam &params)
{
    boundMedianFilter(op(), params.mask->data);
}

void stdDesvFilter(Image<double> &op, FilterParam &params)
{
    pixelDesvFilter(op(), params.factor);
}

/** Apply some filter operation on images, or selfiles */
class ProgFilter2: public XmippMetadataProgram
{
private:
    ///Functions pointers to selected operation
    ImageFilterOperator *filterOperator;

    FileName fn2;
    FileName maskName;
    MetaData md2;
    Image<double> img2;
    bool isValue;
    double value;
    String file_or_value;
    FilterParam params;

protected:
    void defineParams()
    {
        each_image_produces_an_output = true;
        addUsageLine("Apply some filters to images.");
        XmippMetadataProgram::defineParams();
        addParamsLine(" === Median Filter ===");
        addParamsLine("   -forcePositive          :Applied at those negative values. Positive values are untouched.");
        addParamsLine("or -badPixels <mask>       :Applied at those pixels given by mask.");
        addParamsLine("or -desvPixels <factor>    :Applied at those pixels out of the range [mean - factor*std, mean + factor*std].");

    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        filterOperator = NULL;
        isValue = false;
        // Check operation to do
        if (checkParam("-forcePositive"))
            filterOperator = forcePos;
        else if (checkParam("-badPixels"))
        {
            params.mask = new Image<char>;
            maskName = getParam("-badPixels");
            params.mask->read(maskName);
            filterOperator = removeBadPixels;
        }
        else if (checkParam("-desvPixels"))
        {
            params.factor = getDoubleParam("-desvPixels");
            filterOperator = stdDesvFilter;
        }

    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg,mdIn,objId);

        if (filterOperator != NULL)
            filterOperator(img, params);

        img.write(fnImgOut);
    }
}
;//end of class ProgOperate


int main(int argc, char **argv)
{
    ProgFilter2 program;
    program.read(argc, argv);
    program.tryRun();
}




