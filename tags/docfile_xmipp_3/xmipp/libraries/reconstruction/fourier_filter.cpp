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
#include <data/volume.h>
#include <data/image.h>
#include <data/mask.h>
#include <data/fft.h>

/* Clear ------------------------------------------------------------------- */
void FourierMask::clear()
{
    FilterShape = RAISED_COSINE;
    FilterBand = LOWPASS;
    w2 = w1 = 0;
    raised_w = 0;
    ctf.clear();
    ctf.enable_CTFnoise = false;
    correctPhase = false;
    do_generate_3dmask = false;
}

/* Assignment -------------------------------------------------------------- */
FourierMask & FourierMask::operator = (const FourierMask &F)
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
void FourierMask::assign(const FourierMask &F)
{
    *this = F;
}

/* Read parameters from command line. -------------------------------------- */
void FourierMask::read(int argc, char **argv)
{
    clear();

    // Filter shape .........................................................
    int i = paremeterPosition(argc, argv, "-fourier_mask");
    if (i + 1 >= argc) REPORT_ERROR(3000, "FourierMask: -fourier_mask with no mask_type");
    if (i == -1)
    {
        // The default is to use raised_cosine with width 0.02
        raised_w = 0.02;
        FilterShape = RAISED_COSINE;
    }
    else if (strcmp(argv[i+1], "raised_cosine") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(3000, "FourierMask: Raised cosine needs a number of pixels");
        raised_w = textToFloat(argv[i+2]);
        FilterShape = RAISED_COSINE;
    }
    else if (strcmp(argv[i+1], "wedge") == 0)
    {
        if (i + 3 >= argc)
            REPORT_ERROR(3000, "FourierMask: Wedge needs two angle parameters");
        w1 = textToFloat(argv[i+2]);
        w2 = textToFloat(argv[i+3]);
        FilterShape = WEDGE;
        do_generate_3dmask = true;
    }
    else if (strcmp(argv[i+1], "cone") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(3000, "FourierMask: Cone needs one angle parameter");
        w1 = textToFloat(argv[i+2]);
        FilterShape = CONE;
        do_generate_3dmask = true;
    }
    else if (strcmp(argv[i+1], "gaussian") == 0)
    {
        FilterShape = GAUSSIAN;
        FilterBand = LOWPASS;
    }
    else if (strcmp(argv[i+1], "real_gaussian") == 0)
    {
        FilterShape = REALGAUSSIAN;
        FilterBand = LOWPASS;
    }
    else if (strcmp(argv[i+1], "ctf") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(3000, "FourierMask: CTF needs a CTF file");
        FilterShape = FilterBand = CTF;
        ctf.enable_CTFnoise = false;
        ctf.read(argv[i+2]);
        ctf.Produce_Side_Info();
    }
    else if (strcmp(argv[i+1], "ctfpos") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(3000, "FourierMask: CTF needs a CTF file");
        FilterShape = FilterBand = CTFPOS;
        ctf.enable_CTFnoise = false;
        ctf.read(argv[i+2]);
        ctf.Produce_Side_Info();
        correctPhase=true;
    }
    else if (strcmp(argv[i+1], "bfactor") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(3000, "FourierMask: Bfactor needs a Bfactor in Ang^2");
        FilterShape = FilterBand = BFACTOR;
        w1 = textToFloat(argv[i+2]);
        if (!checkParameter(argc, argv, "-sampling"))
            REPORT_ERROR(3000, "FourierMask: Bfactor needs a -sampling argument");
        w2 = textToFloat(getParameter(argc, argv, "-sampling"));
    }

    else
        REPORT_ERROR(3000, "FourierMask: Unknown filter type");

    // Filter band ..........................................................
    if (checkParameter(argc, argv, "-low_pass"))
    {
        w1 = textToFloat(getParameter(argc, argv, "-low_pass"));
        FilterBand = LOWPASS;
    }
    else if (checkParameter(argc, argv, "-high_pass"))
    {
        w1 = textToFloat(getParameter(argc, argv, "-high_pass"));
        FilterBand = HIGHPASS;
    }
    else if (checkParameter(argc, argv, "-band_pass"))
    {
        if (!getTwoDoubleParams(argc, argv, "-band_pass", w1, w2, 0, 0))
            REPORT_ERROR(3000, "FourierMask: Not enough parameters after -band_pass");
        FilterBand = BANDPASS;
    }
    else if (checkParameter(argc, argv, "-stop_band"))
    {
        if (!getTwoDoubleParams(argc, argv, "-stop_band", w1, w2, 0, 0))
            REPORT_ERROR(3000, "FourierMask: Not enough parameters after -stop_band");
        FilterBand = STOPBAND;
    }

    if (!(FilterBand == BFACTOR) && checkParameter(argc, argv, "-sampling"))
    {
        double sampling_rate = textToFloat(getParameter(argc, argv, "-sampling"));
        if (w1 != 0)       w1 = sampling_rate / w1;
        if (w2 != 0)       w2 = sampling_rate / w2;
    }
}

/* Show -------------------------------------------------------------------- */
void FourierMask::show()
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

/* Usage ------------------------------------------------------------------- */
void FourierMask::usage()
{
    std::cerr << "   -low_pass  <w1>                   : Cutoff freq (<1/2 or A)\n"
              << "   -high_pass <w1>                   : Cutoff freq (<1/2 or A)\n"
              << "   -band_pass <w1> <w2>              : Cutoff freq (<1/2 or A)\n"
              << "   -stop_band <w1> <w2>              : Cutoff freq (<1/2 or A)\n"
              << "   -fourier_mask raised_cosine <raisedw>: Use raised cosine edges (in dig.freq.)\n"
              << "   -fourier_mask wedge <th0> <thF>   : Missing wedge (along y) for data between th0-thF \n"
              << "   -fourier_mask cone <th0>          : Missing cone for tilt angles up to th0 \n"
              << "   -fourier_mask gaussian            : sigma=<w1>\n"
              << "   -fourier_mask real_gaussian       : convolution with a Gaussian in real-space with sigma=<w1>\n"
              << "   -fourier_mask ctf                 : Provide a .ctfparam file\n"
              << "   -fourier_mask ctfpos              : Provide a .ctfparam file\n"
              << "                                       The CTF phase will be corrected before applying\n"
              << "   -fourier_mask bfactor <B>         : Exponential filter (positive values for decay) \n"
              << "  [-sampling <sampling_rate>]        : If provided pass frequencies are taken in Ang \n"
    ;
}

/* Get mask value ---------------------------------------------------------- */
double FourierMask::maskValue(const Matrix1D<double> &w)
{
    double absw=w.module();
    
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
                    else return 0;
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
                    else return 0;
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
                    else return 0;
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
                    else return 1;
                    break;
            }
            break;
        case CTF:
            return ctf.CTF_at(XX(w)/ctf.Tm,YY(w)/ctf.Tm);
            break;
        case CTFPOS:
            return ABS(ctf.CTF_at(XX(w)/ctf.Tm,YY(w)/ctf.Tm));
            break;
        case BFACTOR:
            double R = absw / w2;
            return exp( - (w1 / 4.)  * R * R);
            break;
    }
}

/* Generate mask ----------------------------------------------------------- */
void FourierMask::generate_mask(Matrix2D<double> &v)
{
    transformer.setReal(v);
    Matrix2D< std::complex<double> > Fourier;
    transformer.getFourierAlias(Fourier);
    maskFourier2D.resize(Fourier);
    w.resize(2);
    for (int i=0; i<YSIZE(maskFourier2D); i++)
    {
        FFT_IDX2DIGFREQ(i,YSIZE(v),YY(w));
        for (int j=0; j<XSIZE(maskFourier2D); j++)
        {
            FFT_IDX2DIGFREQ(j,XSIZE(v),XX(w));
            DIRECT_MAT_ELEM(maskFourier2D,i,j)=maskValue(w);
        }
    }
}

/* Generate mask ----------------------------------------------------------- */
void FourierMask::generate_mask(Matrix3D<double> &v)
{
    
    if (do_generate_3dmask)
    {
        transformer.setReal(v);
        Matrix3D< std::complex<double> > Fourier;
        transformer.getFourierAlias(Fourier);
        maskFourier3D.resize(Fourier);
        Matrix3D<double> mask(v);
        mask.setXmippOrigin();
        Matrix2D<double> A(3,3);
        A.initIdentity();
        Matrix3D<int> imask;

        VolumeXmipp Vt;
        switch (FilterShape)
        {
            case WEDGE:
                BinaryWedgeMask(mask, w1, w2, A);
                break;
            case CONE:
                imask.resize(v);
                BinaryConeMask(imask,90. - w1);
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(imask)
                {
                    DIRECT_MULTIDIM_ELEM(mask,n)=(double)(DIRECT_MULTIDIM_ELEM(imask,n));
                }
                Vt=mask;
                Vt.write("cone.vol");
                break;
        }
        
        // Transfer mask to maskFourier3D
        for (int z=0, ii=0, yy, zz; z<ZSIZE(maskFourier3D); z++)
        {
            if ( z > (ZSIZE(maskFourier3D) - 1)/2 ) 
                zz = z-ZSIZE(maskFourier3D);
            else 
                zz = z;
            for (int y=0; y<YSIZE(maskFourier3D); y++)
            {
                if ( y > (YSIZE(maskFourier3D) - 1)/2 ) 
                    yy = y-YSIZE(maskFourier3D);
                else 
                    yy = y;
                for (int xx=0; xx<XSIZE(maskFourier3D); xx++,ii++)
                {
                    DIRECT_MULTIDIM_ELEM(maskFourier3D,ii) = 
                        VOL_ELEM(mask,zz,yy,xx);
                }
            }
        }
    }
 }

/* Apply mask -------------------------------------------------------------- */
void FourierMask::apply_mask_Space(Matrix2D<double> &v)
{
    Matrix2D< std::complex<double> > aux2D;
    transformer.FourierTransform(v, aux2D, false);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(aux2D)
        DIRECT_MULTIDIM_ELEM(aux2D,n)*=DIRECT_MULTIDIM_ELEM(maskFourier2D,n);
    transformer.inverseFourierTransform();
}

void FourierMask::apply_mask_Space(Matrix3D<double> &v)
{
    Matrix3D< std::complex<double> > aux3D;
    transformer.FourierTransform(v, aux3D, false);

    if (do_generate_3dmask)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(aux3D)
            DIRECT_MULTIDIM_ELEM(aux3D,n)*=DIRECT_MULTIDIM_ELEM(maskFourier3D,n);
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
                    DIRECT_VOL_ELEM(aux3D,k,i,j)*=maskValue(w);
                }
            }
        }
    }
    transformer.inverseFourierTransform();
}

/* Mask power -------------------------------------------------------------- */
double FourierMask::mask2D_power()
{
    return maskFourier2D.sum2()/MULTIDIM_SIZE(maskFourier2D);
}

// Correct phase -----------------------------------------------------------
void FourierMask::correct_phase()
{
    FOR_ALL_ELEMENTS_IN_MATRIX2D(maskFourier2D)
        if (maskFourier2D(i, j) < 0) maskFourier2D(i, j) *= -1;
}
