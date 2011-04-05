/***************************************************************************
 * Authors:     J.M.de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
#include "data/filters.h"
#include "reconstruction/fourier_filter.h"
#include "reconstruction/denoise.h"
#include "reconstruction/mean_shift.h"

/** Apply some filter operation on images, or selfiles */
class ProgFilter: public XmippMetadataProgram
{
private:
    ///Functions pointers to selected operation
    XmippFilter * filter;

protected:
    void defineParams()
    {
        each_image_produces_an_output = true;
        addUsageLine("Apply different type of filters to images or volumes.");
        XmippMetadataProgram::defineParams();
        FourierFilter::defineParams(this);
        DenoiseFilter::defineParams(this);
        BadPixelFilter::defineParams(this);
        MeanShiftFilter::defineParams(this);
        BackgroundFilter::defineParams(this);
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();

        if (checkParam("--fourier"))
          filter = new FourierFilter();
        else if (checkParam("--wavelet"))
          filter = new DenoiseFilter();
        else if (checkParam("--bad_pixels"))
          filter = new BadPixelFilter();
        else if (checkParam("--mean_shift"))
          filter = new MeanShiftFilter();
        else if (checkParam("--background"))
            filter = new BackgroundFilter();
        //Read params
        filter->readParams(this);
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg, mdIn, objId);
        filter->apply(img());
        img.write(fnImgOut);
    }
}
;//end of class ProgOperate


int main(int argc, char **argv)
{
    ProgFilter program;
    program.read(argc, argv);
    return program.tryRun();
}




