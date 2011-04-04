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
        BadPixelFilter::defineParams(this);
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();

        if (checkParam("--fourier"))
        {
          FourierFilter *ff = new FourierFilter();
          ff->readParams(this);
          filter = ff;
        }
        else if (checkParam("--bad_pixels"))
        {
          BadPixelFilter *bpf = new BadPixelFilter();
          bpf->readParams(this);
          filter = bpf;
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg,mdIn,objId);
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




