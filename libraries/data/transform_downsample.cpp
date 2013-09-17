/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "transform_downsample.h"
#include "args.h"
#include "mask.h"
#include "xmipp_fftw.h"
#include "xvsmooth.h"

// Read --------------------------------------------------------------------
void ProgTransformDownsample::readParams()
{

    XmippMetadataProgram::readParams();
    step=getDoubleParam("--step");
    String strMethod=getParam("--method");
    if (strMethod=="fourier")
    {
        nThreads=getIntParam("--method",1);
        method=FOURIER;
    }
    else if (strMethod=="smooth")
        method=SMOOTH;
    else
        method=KER_RECTANGLE;
}

// Usage -------------------------------------------------------------------
void ProgTransformDownsample::defineParams()
{
    each_image_produces_an_output = true;
    addUsageLine("Downsample a micrograph. For volumes use xmipp_transform_geometry.");
    addUsageLine("+There are several downsampling methods. The most general and recommended is Fourier.");
    addUsageLine("+Fourier downsampling puts a window in Fourier space. This is the best downsampling that can be performed.");
    addUsageLine("+Altermatively, smoothing makes color dithering which is pretty good for visualization, ");
    addUsageLine("+but it modifies the particle spectrum. Binning with a rectangle kernel modifies the ");
    addUsageLine("+spectrum of the micrographs and is not recommended. You may see the effects of the different ");
    addUsageLine("+downsampling schemes at [[http://biocomp.cnb.csic.es/~coss/Articulos/Sorzano2009d.pdf][this article]].");
    addUsageLine("+ ");
    addUsageLine("+The downsampling factor (--step) is the factor by which the micrograph will be reduced.");
    addUsageLine("+For instance, a downsampling by 2 will reduce the image size to one half. Using Fourier and smooth ");
    addUsageLine("+you may use non-integer downsampling factors, and the image size will be reduced by 1/factor");
    addSeeAlsoLine("transform_geometry");
    XmippMetadataProgram::defineParams();
    addParamsLine("  --step <factor>    : Downsampling factor. factor=2 reduces the image size to one half.");
    addParamsLine("                     :+Fourier and smooth support non-integer downsampling factors.");
    addParamsLine("                     :+Rectangular binning must use integer factors.");
    addParamsLine(" [--method <mth=fourier>]  : Method for making the downsampling");
    addParamsLine("         where <mth>");
    addParamsLine("               fourier <numThreads=1>: Fourier supports non-integer downsampling factors");
    addParamsLine("                                     :+This is the best choice.");
    addParamsLine("               rectangle: This is simple binning in a square of size factor x factor");
    addParamsLine("                        :+This is not a good choice since it creates aliasing and ");
    addParamsLine("                        :+unequal frequency damping.");
    addParamsLine("               smooth: smooth and colordither");
    addParamsLine("                     :+ Both input and output micrographs must be 8 bits, unsigned char");
    addExampleLine("xmipp_transform_downsample -i micrograph.tif -o downsampledMicrograph.tif --step 2");
}

// Downsample micrograph ---------------------------------------------------
void ProgTransformDownsample::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    // Open input data
    ImageGeneric M_in;
    M_in.readOrReadMapped(fnImg);
    size_t Zdim, Ydim, Xdim;
    M_in.getDimensions(Xdim,Ydim,Zdim);
    if (Zdim!=1)
        REPORT_ERROR(ERR_MULTIDIM_DIM,"This program is not intended for volumes");

    // Open output data, mapped file
    int Xpdim = (int)floor(Xdim/step);
    int Ypdim = (int)floor(Ydim/step);
    ImageGeneric M_out;
    if (method==SMOOTH)
        M_out.setDatatype(DT_UChar);
    else
        M_out.setDatatype(DT_Float);
    M_out.mapFile2Write(Xpdim, Ypdim, 1, fnImgOut,fnImg==fnImgOut);

    // Downsample
    if (method==KER_RECTANGLE)
        downsampleKernel(M_in,step,M_out);
    else if (method==FOURIER)
        downsampleFourier(M_in,step,M_out,nThreads);
    else
        downsampleSmooth(M_in,M_out);

    M_out.write(fnImgOut);
}

/* Downsample -------------------------------------------------------------- */
void downsampleKernel(const ImageGeneric &M, double step, ImageGeneric &Mp)
{
    int istep=(int)step;
    MultidimArray<double> kernel;
    kernel.resizeNoCopy(istep,istep);
    kernel.initConstant(1.0/MULTIDIM_SIZE(kernel));

    size_t Ydim, Xdim, Ypdim, Xpdim;
    M().getDimensions(Xdim, Ydim);
    Mp().getDimensions(Xpdim, Ypdim);

    // Look for input/output ranges
    double a = 1;
    double b = 0;
    double scale = 1;
    size_t ii, jj, i2, j2, i, j, y, x;
    if (Mp.getDatatype() != DT_Float)
    {
        double imin, imax;
        double omin, omax;
        imin=imax=M.getPixel(0,0);
        bool ofirst = true;

        if (M.getDatatype() != DT_Float)
            scale = (pow(2.0, Mp.getDatatypeDepth()) - 1.0) /
                    (pow(2.0, M.getDatatypeDepth()) - 1.0);
        else if (M.getDatatype() == DT_Float)
            scale = 1;
        for (ii = 0, y = 0; y < Ydim && ii < Ypdim; y += istep, ++ii)
            for (jj = 0, x = 0; x < Xdim && jj < Xpdim; x += istep, ++jj)
            {
                double pixval = 0;
                for (i=0, i2=y; i<YSIZE(kernel) && i2<Ydim; ++i, ++i2)
                    for (j=0, j2=x; j<XSIZE(kernel) && j2<Xdim; ++j, ++j2)
                    {
                        double aux=M.getPixel(i2, j2);
                        imin = XMIPP_MIN(imin, aux);
                        imax = XMIPP_MAX(imax, aux);
                        pixval += A2D_ELEM(kernel,i, j) * aux;
                    }
                pixval *= scale;
                if (ofirst)
                {
                    omin = omax = pixval;
                    ofirst = false;
                }
                else
                {
                    omin = XMIPP_MIN(omin, pixval);
                    omax = XMIPP_MAX(omax, pixval);
                }
            }

        // Compute range transformation
        double irange = imax - imin;
        double orange = omax - omin;

        if (M.getDatatype() != DT_Float)
        {
            a = scale * irange / orange;
            b = -omin;
        }
        else if (Mp.getDatatype() != DT_Float)
        {
            a = (pow(2.0, Mp.getDatatypeDepth()) - 1.0) / orange;
            scale = 1;
            b = -omin;
        }
    }

    // Really downsample
    for (ii = 0, y = 0; y < Ydim && ii < Ypdim; y += istep, ++ii)
        for (jj = 0, x = 0; x < Xdim && jj < Xpdim; x += istep, ++jj)
        {
            double pixval = 0;
            for (i=0, i2=y; i<YSIZE(kernel) && i2<Ydim; ++i, ++i2)
                for (j=0, j2=x; j<XSIZE(kernel)&& j2<Xdim; ++j, ++j2)
                    pixval += A2D_ELEM(kernel,i, j) * M.getPixel(i2, j2);
            if (ii < Ypdim && jj < Xpdim)
            {
                if (Mp.datatype != DT_Float)
                    Mp.setPixel(ii, jj, floor(a*(pixval*scale + b)));
                else
                    Mp.setPixel(ii, jj, pixval);
            }
        }
}

void downsampleFourier(const ImageGeneric &M, double step, ImageGeneric &Mp, int nThreads)
{
    size_t Ydim, Xdim, Ypdim, Xpdim;
    M().getDimensions(Xdim, Ydim);
    Mp().getDimensions(Xpdim, Ypdim);

    // Read the micrograph in memory as doubles
    MultidimArray<double> Mmem;
    Mmem.setMmap(true);
    Mmem.resizeNoCopy(Ydim,Xdim);
    MultidimArray<std::complex<double> > MmemFourier;
    M().getImage(Mmem);

    // Perform the Fourier transform
    FourierTransformer transformerM;
    transformerM.setThreadsNumber(nThreads);
    transformerM.FourierTransform(Mmem, MmemFourier, false);

    // Create space for the downsampled image and its Fourier transform
    MultidimArray<double> Mpmem(Ypdim,Xpdim);
    MultidimArray<std::complex<double> > MpmemFourier;
    FourierTransformer transformerMp;
    transformerMp.setThreadsNumber(nThreads);
    transformerMp.setReal(Mpmem);
    transformerMp.getFourierAlias(MpmemFourier);

    int ihalf=YSIZE(MpmemFourier)/2+1;
    for (int i=0; i<ihalf; i++)
        for (size_t j=0; j<XSIZE(MpmemFourier); j++)
            A2D_ELEM(MpmemFourier,i,j)=A2D_ELEM(MmemFourier,i,j);
    for (size_t i=ihalf; i<YSIZE(MpmemFourier); i++)
    {
        size_t ip=YSIZE(MmemFourier)-YSIZE(MpmemFourier)+i;
        for (size_t j=0; j<XSIZE(MpmemFourier); j++)
            A2D_ELEM(MpmemFourier,i,j)=A2D_ELEM(MmemFourier,ip,j);
    }

    // Transform data
    transformerMp.inverseFourierTransform();

    // Find minimun and range in output data
    double omin,omax;
    Mpmem.computeDoubleMinMax(omin,omax);
    double orange = omax - omin;
    double a = (pow(2.0, Mp.getDatatypeDepth()) - 1.0) / orange;
    double b = -omin;
    double scale=1;
    if (M.getDatatype() != DT_Float)
        scale = (pow(2.0, Mp.getDatatypeDepth()) - 1.0) /
                (pow(2.0, M.getDatatypeDepth()) - 1.0);
    else if (M.getDatatype() == DT_Float)
        scale = 1;

    // Copy back data
    if (Mp.datatype!=DT_Float)
        FOR_ALL_ELEMENTS_IN_ARRAY2D(Mpmem)
        Mp.setPixel(i, j, a*(A2D_ELEM(Mpmem,i,j)*scale + b));
    else
        FOR_ALL_ELEMENTS_IN_ARRAY2D(Mpmem)
        Mp.setPixel(i, j, A2D_ELEM(Mpmem,i,j));
}

void downsampleSmooth(const ImageGeneric &M, ImageGeneric &Mp)
{
    if (Mp.datatype!=DT_UChar)
        REPORT_ERROR(ERR_ARG_INCORRECT,formatString("Smooth downsampling is only valid for 8 bit images. \n"
        		"Choose a supporting 8bit file format different from %s",Mp.image->name().c_str()));

    size_t Ydim, Xdim, Ypdim, Xpdim;
    M().getDimensions(Xdim, Ydim);
    Mp().getDimensions(Xpdim, Ypdim);
    byte *inputImage=NULL;
    MultidimArray<unsigned char> Maux;
    if (M.datatype==DT_UChar)
        M().getArrayPointer(inputImage);
    else
    {
        Maux.setMmap(true);
        M().getImage(Maux);
        inputImage=MULTIDIM_ARRAY(Maux);
    }
    unsigned char *outputImage;
    Mp().getArrayPointer(outputImage);
    SmoothResize(inputImage,outputImage, Xdim, Ydim, Xpdim, Ypdim);
}
