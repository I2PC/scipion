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

#include "fourier_filter.h"

#include <data/args.h>
#include <data/image.h>
#include <data/mask.h>
#include <data/fft.h>

/* Clear ------------------------------------------------------------------- */
void ProgFourierFilter::clear()
{
    FilterShape = RAISED_COSINE;
    FilterBand = LOWPASS;
    w2 = w1 = 0;
    raised_w = 0;
    ctf.clear();
    ctf.enable_CTFnoise = false;
    do_correct_phase = false;
    do_generate_3dmask = false;
}

/* Empty constructor ---------------------------------------------------------*/
ProgFourierFilter::ProgFourierFilter()
{
    clear();
}

/* Assignment -------------------------------------------------------------- */
ProgFourierFilter & ProgFourierFilter::operator = (const ProgFourierFilter &F)
{
    if (this != &F)
    {
        clear();

        FilterShape = F.FilterShape;
        FilterBand = F.FilterBand;
        w2 = F.w2;
        w1 = F.w1;
        raised_w = F.raised_w;
        ctf = F.ctf;
    }
    return *this;
}

/* Another function for assignment ----------------------------------------- */
void ProgFourierFilter::assign(const ProgFourierFilter &F)
{
    *this = F;
}


/* Define params ------------------------------------------------------------------- */
void ProgFourierFilter::defineParams()
{
    each_image_produces_an_output = true;
    XmippMetadataProgram::defineParams();
    addParamsLine("   --low_pass  <w1>                   : Cutoff freq (<1/2 or A)");
    addParamsLine("or --high_pass <w1>                   : Cutoff freq (<1/2 or A)");
    addParamsLine("or --band_pass <w1> <w2>              : Cutoff freq (<1/2 or A)");
    addParamsLine("or --stop_band <w1> <w2>              : Cutoff freq (<1/2 or A)");
    addParamsLine("  [--fourier_mask <mask_type=raised_cosine>]: Fourier mask that will be applied.");
    addParamsLine("        where <mask_type> ");
    addParamsLine("         raised_cosine <raisedw=0.02>     : Use raised cosine edges (in dig.freq.)");
    addParamsLine("         wedge <th0> <thF>           : Missing wedge (along y) for data between th0-thF ");
    addParamsLine("         cone <th0>                  : Missing cone for tilt angles up to th0 ");
    addParamsLine("         gaussian                    : sigma=<w1>");
    addParamsLine("         real_gaussian               : convolution with a Gaussian in real-space with sigma=<w1>");
    addParamsLine("         ctf <ctfile>                : Provide a .ctfparam file");
    addParamsLine("         ctfpos <ctfile>             : Provide a .ctfparam file");
    addParamsLine("                                     : The CTF phase will be corrected before applying");
    addParamsLine("         bfactor <B>                 : Exponential filter (positive values for decay) ");
    addParamsLine("           requires --sampling;                                                         ");
    addParamsLine("  [--sampling <sampling_rate>]        : If provided pass frequencies are taken in Ang ");
}

/* Read parameters from command line. -------------------------------------- */
void ProgFourierFilter::readParams()
{
    clear();
    XmippMetadataProgram::readParams();
    // Filter shape .........................................................
    std::string mask_type;
    mask_type = getParam("--fourier_mask");

    if (mask_type == "raised_cosine")
    {
        raised_w = getDoubleParam("--fourier_mask", "raised_cosine");
        FilterShape = RAISED_COSINE;
    }
    else if (mask_type == "wedge")
    {
        w1 = getDoubleParam("--fourier_mask", "wedge", 0);
        w2 = getDoubleParam("--fourier_mask", "wedge", 1);
        FilterShape = WEDGE;
        do_generate_3dmask = true;
    }
    else if (mask_type == "cone")
    {
        w1 = getDoubleParam("--fourier_mask", "cone", 0);
        FilterShape = CONE;
        do_generate_3dmask = true;
    }
    else if (mask_type == "gaussian")
    {
        FilterShape = GAUSSIAN;
        FilterBand = LOWPASS;
    }
    else if (mask_type == "real_gaussian")
    {
        FilterShape = REALGAUSSIAN;
        FilterBand = LOWPASS;
    }
    else if (mask_type == "ctf")
    {
        FilterShape = FilterBand = CTF;
        ctf.enable_CTFnoise = false;
        ctf.read( getParam("--fourier_mask", "ctf") );
        ctf.Produce_Side_Info();
    }
    else if (mask_type == "ctfpos")
    {
        FilterShape = FilterBand = CTFPOS;
        ctf.enable_CTFnoise = false;
        ctf.read( getParam("--fourier_mask", "ctfpos") );
        ctf.Produce_Side_Info();
        do_correct_phase = true;
    }
    else if (mask_type == "bfactor")
    {
        FilterShape = FilterBand = BFACTOR;
        w1 = getDoubleParam("--fourier_mask", "bfactor");
        w2 = getDoubleParam("--sampling");
    }
    else
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE, "This couldn't happen, check argument parser or params definition");

    // Filter band ..........................................................
    if (checkParam("--low_pass"))
    {
        w1 = getDoubleParam("--low_pass");
        FilterBand = LOWPASS;
    }
    else if (checkParam("--high_pass"))
    {
        w1 = getDoubleParam("--high_pass");
        FilterBand = HIGHPASS;
    }
    else if (checkParam("--band_pass"))
    {
        w1 = getDoubleParam("--band_pass", 0);
        w2 = getDoubleParam("--band_pass", 1);
        FilterBand = BANDPASS;
    }
    else if (checkParam("--stop_band"))
    {
        w1 = getDoubleParam("--stop_band", 0);
        w2 = getDoubleParam("--stop_band", 1);
        FilterBand = STOPBAND;
    }

    if (!(FilterBand == BFACTOR) && checkParam("--sampling"))
    {
        double sampling_rate = getDoubleParam("--sampling");
        if (w1 != 0)
            w1 = sampling_rate / w1;
        if (w2 != 0)
            w2 = sampling_rate / w2;
    }
}

/* Show -------------------------------------------------------------------- */
void ProgFourierFilter::show()
{
    if (FilterShape == WEDGE)
    {
        std::cout << "Missing wedge for data between tilting angles of " << w1 << " and " << w2 << " deg\n";
    }
    else if (FilterShape == CONE)
    {
        std::cout << "Missing cone for RCT data with tilting angles up to " << w1 << " deg\n";
    }
    else
    {
        std::cout << "Filter Band: ";
        switch (FilterBand)
        {
        case LOWPASS:
            std::cout << "Lowpass before " << w1 << std::endl;
            break;
        case HIGHPASS:
            std::cout << "Highpass after " << w1 << std::endl;
            break;
        case BANDPASS:
            std::cout << "Bandpass between " << w1 << " and " << w2 << std::endl;
            break;
        case STOPBAND:
            std::cout << "Stopband between " << w1 << " and " << w2 << std::endl;
            break;
        case CTF:
            std::cout << "CTF\n";
            break;
        case CTFPOS:
            std::cout << "CTFPOS\n";
            break;
        case BFACTOR:
            std::cout << "Bfactor "<< w1 <<std::endl;
            break;
        }
        std::cout << "Filter Shape: ";
        switch (FilterShape)
        {
        case RAISED_COSINE:
            std::cout << "Raised cosine with " << raised_w
            << " raised frequencies\n";
            break;
        case GAUSSIAN:
            std::cout << "Gaussian\n";
            break;
        case REALGAUSSIAN:
            std::cout << "Real Gaussian\n";
            break;
        case CTF:
            std::cout << "CTF\n" << ctf;
            break;
        case CTFPOS:
            std::cout << "CTFPOS\n" << ctf;
            break;

        }
    }
}

void ProgFourierFilter::processImage()
{
    static bool first = false;

    img.read(fnImg);
    if (first)
    {
        generateMask(img());
        first = false;
    }
    applyMaskSpace(img());
    img.write(fnImgOut);
}

/* Get mask value ---------------------------------------------------------- */
double ProgFourierFilter::maskValue(const Matrix1D<double> &w)
{
    double absw = w.module();

    // Generate mask
    switch (FilterBand)
    {
    case LOWPASS:
        switch (FilterShape)
        {
        case RAISED_COSINE:
            if (absw<w1)
                return 1;
            else if (absw<w1+raised_w)
                return (1+cos(PI/raised_w*(absw-w1)))/2;
            else
                return 0;
            break;
        case GAUSSIAN:
            return 1/sqrt(2*PI*w1)*exp(-0.5*absw*absw/(w1*w1));
            break;
        case REALGAUSSIAN:
            return exp(-PI*PI*absw*absw*w1*w1);
            break;
        }
        break;
    case HIGHPASS:
        switch (FilterShape)
        {
        case RAISED_COSINE:
            if (absw>w1)
                return 1;
            else if (absw>w1-raised_w)
                return (1+cos(PI/raised_w*(w1-absw)))/2;
            else
                return 0;
            break;
        }
        break;
    case BANDPASS:
        switch (FilterShape)
        {
        case RAISED_COSINE:
            if (absw>=w1 && absw<=w2)
                return 1;
            else if (absw>w1-raised_w && absw<w1)
                return (1+cos(PI/raised_w*(w1-absw)))/2;
            else if (absw<w2+raised_w && absw>w2)
                return (1+cos(PI/raised_w*(w2-absw)))/2;
            else
                return 0;
            break;
        }
        break;
    case STOPBAND:
        switch (FilterShape)
        {
        case RAISED_COSINE:
            if (absw>=w1 && absw<=w2)
                return 0;
            else if (absw>w1-raised_w && absw<w1)
                return 1-(1+cos(PI/raised_w*(w1-absw)))/2;
            else if (absw<w2+raised_w && absw>w2)
                return 1-(1+cos(PI/raised_w*(w2-absw)))/2;
            else
                return 1;
            break;
        }
        break;
    case CTF:
    	ctf.precomputeValues(XX(w)/ctf.Tm,YY(w)/ctf.Tm);
        return ctf.CTF_at();
        break;
    case CTFPOS:
    	ctf.precomputeValues(XX(w)/ctf.Tm,YY(w)/ctf.Tm);
        return ABS(ctf.CTF_at());
        break;
    case BFACTOR:
        double R = absw / w2;
        return exp( - (w1 / 4.)  * R * R);
        break;
    }
}

/* Generate mask ----------------------------------------------------------- */
void ProgFourierFilter::generateMask(MultidimArray<double> &v)
{
    if (do_generate_3dmask)
    {
        transformer.setReal(v);
        MultidimArray< std::complex<double> > Fourier;
        transformer.getFourierAlias(Fourier);
        maskFourier.resize(Fourier);
        maskFourier.setXmippOrigin();

        if (FilterShape==WEDGE || FilterShape==CONE)
        {
            switch (FilterShape)
            {
            case WEDGE:
                {
                    Matrix2D<double> A(3,3);
                    A.initIdentity();
                    BinaryWedgeMask(maskFourier, w1, w2, A);
                    break;
                }
            case CONE:
                BinaryConeMask(maskFourier, 90. - w1);
                break;
            }
        }
    }
}

void ProgFourierFilter::applyMaskSpace(MultidimArray<double> &v)
{
    MultidimArray< std::complex<double> > aux3D;
    transformer.FourierTransform(v, aux3D, false);

    if (do_generate_3dmask)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(aux3D)
        DIRECT_MULTIDIM_ELEM(aux3D,n)*=DIRECT_MULTIDIM_ELEM(maskFourier,n);
    }
    else
    {
        w.resize(3);
        for (int k=0; k<ZSIZE(aux3D); k++)
        {
            FFT_IDX2DIGFREQ(k,ZSIZE(v),ZZ(w));
            for (int i=0; i<YSIZE(aux3D); i++)
            {
                FFT_IDX2DIGFREQ(i,YSIZE(v),YY(w));
                for (int j=0; j<XSIZE(aux3D); j++)
                {
                    FFT_IDX2DIGFREQ(j,XSIZE(v),XX(w));
                    DIRECT_A3D_ELEM(aux3D,k,i,j)*=maskValue(w);
                }
            }
        }
    }
    transformer.inverseFourierTransform();
}

/* Mask power -------------------------------------------------------------- */
double ProgFourierFilter::maskPower()
{
    return maskFourier.sum2()/MULTIDIM_SIZE(maskFourier);
}

// Correct phase -----------------------------------------------------------
void ProgFourierFilter::correctPhase()
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(maskFourier)
    if (DIRECT_MULTIDIM_ELEM(maskFourier,n)< 0)
        DIRECT_MULTIDIM_ELEM(maskFourier,n)*= -1;
}
