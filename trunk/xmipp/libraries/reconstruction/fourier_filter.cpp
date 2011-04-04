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
void FourierFilter::init()
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
FourierFilter::FourierFilter()
{
    init();
}

/* Define params ------------------------------------------------------------------- */
void FourierFilter::defineParams(XmippProgram *program)
{
    program->addParamsLine("   --fourier <cuttoff_type> <mask_type=raised_cosine>    : Filter in Fourier space");
    program->addParamsLine("         where <cuttoff_type>");
    program->addParamsLine("            low_pass  <w1>                   : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            high_pass <w1>                   : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            band_pass <w1> <w2>              : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            stop_band <w1> <w2>              : Cutoff freq (<1/2 or A)");
    program->addParamsLine("         where <mask_type> ");
    program->addParamsLine("            raised_cosine <raisedw=0.02>     : Use raised cosine edges (in dig.freq.)");
    program->addParamsLine("            wedge <th0> <thF>                : Missing wedge (along y) for data between th0-thF ");
    program->addParamsLine("            cone <th0>                       : Missing cone for tilt angles up to th0 ");
    program->addParamsLine("            gaussian                         : sigma=<w1>");
    program->addParamsLine("            real_gaussian                    : convolution with a Gaussian in real-space with sigma=<w1>");
    program->addParamsLine("            ctf <ctfile>                     : Provide a .ctfparam file");
    program->addParamsLine("            ctfpos <ctfile>                  : Provide a .ctfparam file");
    program->addParamsLine("                                             : The CTF phase will be corrected before applying");
    program->addParamsLine("            bfactor <B>                      : Exponential filter (positive values for decay) ");
    program->addParamsLine("               requires --sampling;                                                         ");
    program->addParamsLine("  [--sampling <sampling_rate>]        : If provided pass frequencies are taken in Ang ");
    program->addParamsLine("         alias -s;");
    program->addParamsLine("         requires --fourier;");
}

/* Read parameters from command line. -------------------------------------- */
void FourierFilter::readParams(XmippProgram *program)
{
    init();
    // Filter shape .........................................................
    String cuttoff, mask;
    int maskPos = 2;
    cuttoff = program->getParam("--fourier");

    // Read frequencies cuttoff options
    if (cuttoff == "low_pass")
    {
        w1 = program->getDoubleParam("--fourier", "low_pass");
        FilterBand = LOWPASS;
    }
    else if (cuttoff == "high_pass")
    {
        w1 = program->getDoubleParam("--fourier", "high_pass");
        FilterBand = HIGHPASS;
    }
    else if (cuttoff == "band_pass")
    {
        w1 = program->getDoubleParam("--fourier", "band_pass");
        w2 = program->getDoubleParam("--fourier", "band_pass", 1);
        maskPos = 3;
        FilterBand = BANDPASS;
    }
    else if (cuttoff == "stop_band")
    {
        w1 = program->getDoubleParam("--fourier", "stop_band");
        w2 = program->getDoubleParam("--fourier", "stop_band", 1);
        maskPos = 3;
        FilterBand = STOPBAND;
    }

    if (!(FilterBand == BFACTOR) && program->checkParam("--sampling"))
    {
        double sampling_rate = program->getDoubleParam("--sampling");
        if (w1 != 0)
            w1 = sampling_rate / w1;
        if (w2 != 0)
            w2 = sampling_rate / w2;
    }

    //Read mask options
    mask = program->getParam("--fourier", maskPos);

    if (mask == "raised_cosine")
    {
        raised_w = program->getDoubleParam("--fourier", "raised_cosine");
        FilterShape = RAISED_COSINE;
    }
    else if (mask == "wedge")
    {
        w1 = program->getDoubleParam("--fourier", "wedge", 0);
        w2 = program->getDoubleParam("--fourier", "wedge", 1);
        FilterShape = WEDGE;
        do_generate_3dmask = true;
    }
    else if (mask == "cone")
    {
        w1 = program->getDoubleParam("--fourier", "cone", 0);
        FilterShape = CONE;
        do_generate_3dmask = true;
    }
    else if (mask == "gaussian")
    {
        FilterShape = GAUSSIAN;
        FilterBand = LOWPASS;
    }
    else if (mask == "real_gaussian")
    {
        FilterShape = REALGAUSSIAN;
        FilterBand = LOWPASS;
    }
    else if (mask == "ctf")
    {
        FilterShape = FilterBand = CTF;
        ctf.enable_CTFnoise = false;
        ctf.read( program->getParam("--fourier", "ctf") );
        ctf.Produce_Side_Info();
    }
    else if (mask == "ctfpos")
    {
        FilterShape = FilterBand = CTFPOS;
        ctf.enable_CTFnoise = false;
        ctf.read( program->getParam("--fourier", "ctfpos") );
        ctf.Produce_Side_Info();
        do_correct_phase = true;
    }
    else if (mask == "bfactor")
    {
        FilterShape = FilterBand = BFACTOR;
        w1 = program->getDoubleParam("--fourier", "bfactor");
        w2 = program->getDoubleParam("--sampling");
    }
    else
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE, "This couldn't happen, check argument parser or params definition");
}

/* Show -------------------------------------------------------------------- */
void FourierFilter::show()
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

void FourierFilter::apply(MultidimArray<double> &img)
{
    static bool firstTime = true;

    std::cerr << "DEBUG_JM: firstTime: " << firstTime << std::endl;

    if (firstTime)
    {
        generateMask(img);
        firstTime = false;
    }
    applyMaskSpace(img);
}

/* Get mask value ---------------------------------------------------------- */
double FourierFilter::maskValue(const Matrix1D<double> &w)
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
void FourierFilter::generateMask(MultidimArray<double> &v)
{
    if (do_generate_3dmask)
    {
        transformer.setReal(v);
        MultidimArray< std::complex<double> > Fourier;
        transformer.getFourierAlias(Fourier);
        maskFourier.initZeros(Fourier);
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
    else if (MULTIDIM_SIZE(v)<=1024*1024)
    {
        transformer.setReal(v);
        MultidimArray< std::complex<double> > Fourier;
        transformer.getFourierAlias(Fourier);
        maskFourierd.initZeros(Fourier);
        maskFourierd.setXmippOrigin();

        w.resizeNoCopy(3);
        for (int k=0; k<ZSIZE(Fourier); k++)
        {
            FFT_IDX2DIGFREQ(k,ZSIZE(v),ZZ(w));
            for (int i=0; i<YSIZE(Fourier); i++)
            {
                FFT_IDX2DIGFREQ(i,YSIZE(v),YY(w));
                for (int j=0; j<XSIZE(Fourier); j++)
                {
                    FFT_IDX2DIGFREQ(j,XSIZE(v),XX(w));
                    DIRECT_A3D_ELEM(maskFourierd,k,i,j)=maskValue(w);
                }
            }
        }
    }
}

void FourierFilter::applyMaskSpace(MultidimArray<double> &v)
{
    MultidimArray< std::complex<double> > aux3D;
    transformer.FourierTransform(v, aux3D, false);
    if (XSIZE(maskFourier)!=0)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(aux3D)
        DIRECT_MULTIDIM_ELEM(aux3D,n)*=DIRECT_MULTIDIM_ELEM(maskFourier,n);
    }
    else if (XSIZE(maskFourierd)!=0)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(aux3D)
        DIRECT_MULTIDIM_ELEM(aux3D,n)*=DIRECT_MULTIDIM_ELEM(maskFourierd,n);
    }
    else
    {
        w.resizeNoCopy(3);
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
double FourierFilter::maskPower()
{
    if (XSIZE(maskFourier) != 0)
        return maskFourier.sum2()/MULTIDIM_SIZE(maskFourier);
    else if (XSIZE(maskFourierd) != 0)
        return maskFourier.sum2()/MULTIDIM_SIZE(maskFourierd);
}

// Correct phase -----------------------------------------------------------
void FourierFilter::correctPhase()
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(maskFourier)
    if (DIRECT_MULTIDIM_ELEM(maskFourier,n)< 0)
        DIRECT_MULTIDIM_ELEM(maskFourier,n)*= -1;
}
