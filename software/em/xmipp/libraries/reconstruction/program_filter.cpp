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


#include "program_filter.h"
#include "reconstruction/fourier_filter.h"
#include "reconstruction/denoise.h"
#include "reconstruction/mean_shift.h"

ProgFilter::ProgFilter(){}

ProgFilter::~ProgFilter(){}

void ProgFilter::defineParams()
{
    each_image_produces_an_output = true;
    addUsageLine("Apply different type of filters to images or volumes.");
    XmippMetadataProgram::defineParams();
    FourierFilter::defineParams(this);
    WaveletFilter::defineParams(this);
    BadPixelFilter::defineParams(this);
    DiffusionFilter::defineParams(this);
    MeanShiftFilter::defineParams(this);
    BackgroundFilter::defineParams(this);
    MedianFilter::defineParams(this);
    BasisFilter::defineParams(this);
    LogFilter::defineParams(this);
    RetinexFilter::defineParams(this);
    DenoiseTVFilter::defineParams(this);

    //examples
    addExampleLine("Filter a volume using a mask =volumeMask.vol= to remove bad pixels:", false);
    addExampleLine("xmipp_transform_filter -i volume.vol -o volumeFiltered1.vol --bad_pixels mask volumeMask.vol");
    addExampleLine("Remove extreme pixels:", false);
    addExampleLine("xmipp_transform_filter -i volume.vol  -o volumeFiltered2.vol --bad_pixels outliers 0.5");
    addExampleLine("Apply a Fourier filter on a volume:", false);
    addExampleLine("xmipp_transform_filter -i volume.vol -o volumeFiltered.vol --fourier low_pass 0.05");
    addExampleLine("xmipp_transform_filter  -i volume.vol -o volumeFiltered.vol -f band_pass 0.1 0.3");
    addExampleLine("xmipp_transform_filter  -i image.ser  -o imageFiltered.xmp --background plane");
    addExampleLine("xmipp_transform_filter  -i smallStack.stk -o smallFiltered.stk -w DAUB12 difussion");
    addExampleLine("Filter a volume using a wedge mask rotated 10 degress",false);
    addExampleLine("xmipp_transform_filter  --fourier wedge  -60 60 0 0 10 -i ico.spi -o kk0.spi --verbose");
    addExampleLine("Save filtering mask (do not filter)",false);
    addExampleLine("xmipp_transform_filter  --fourier wedge  -60 60 0 0 10 -i ico.spi -o kk0.spi --verbose --save mask.spi");
    addExampleLine("Preprocess image optained in the nikon coolscan",false);
    addExampleLine("xmipp_transform_filter  --log  -i ico.spi -o kk0.spi --fa 4.431 --fb 0.4018 --fc 336.6");
    addExampleLine("xmipp_transform_filter  -i Falcon_2014_07_16-16_32_43_0.mrcs --denoiseTV 50");
}

void ProgFilter::readParams()
{
	readCTF=false;
    XmippMetadataProgram::readParams();

    if (checkParam("--fourier"))
    {
        filter = new FourierFilter();
        String filterType=getParam("--fourier");
        if (filterType=="astigmatism")
        	readCTF=true;
    }
    else if (checkParam("--wavelet"))
        filter = new WaveletFilter();
    else if (checkParam("--bad_pixels"))
        filter = new BadPixelFilter();
    else if (checkParam("--mean_shift"))
        filter = new MeanShiftFilter();
    else if (checkParam("--background"))
        filter = new BackgroundFilter();
    else if (checkParam("--median"))
        filter = new MedianFilter();
    else if (checkParam("--diffusion"))
        filter = new DiffusionFilter();
    else if (checkParam("--basis"))
        filter = new BasisFilter();
    else if (checkParam("--log"))
        filter = new LogFilter();
    else if (checkParam("--retinex"))
        filter = new RetinexFilter();
    else if (checkParam("--denoiseTV"))
        filter = new DenoiseTVFilter();
    else
        REPORT_ERROR(ERR_ARG_MISSING, "You should provide some filter");
    //Read params
    filter->readParams(this);
}

void ProgFilter::preProcess()
{
  if (verbose)
    filter->show();
}

void ProgFilter::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    Image<double> img;
    img.read(fnImg);
    if (readCTF)
    {
    	((FourierFilter *)filter)->ctf.readFromMdRow(rowIn);
    	((FourierFilter *)filter)->ctf.produceSideInfo();
    	((FourierFilter *)filter)->generateMask(img());
    }
    filter->apply(img());
    img.write(fnImgOut);
}
