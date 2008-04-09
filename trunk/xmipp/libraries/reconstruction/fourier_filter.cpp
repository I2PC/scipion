/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "fourier_filter.h"

#include <data/args.h>
#include <data/volume.h>
#include <data/image.h>

/* Clear ------------------------------------------------------------------- */
void FourierMask::clear()
{
    FilterShape = RAISED_COSINE;
    FilterBand = LOWPASS;
    w2 = w1 = 0;
    raised_w = 0;
    ctf.clear();
    ctf.enable_CTFnoise = false;
    mask1D.clear();
    mask2D.clear();
    mask3D.clear();
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
        mask1D = F.mask1D;
        mask2D = F.mask2D;
        mask3D = F.mask3D;
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
        FilterBand = LOWPASS;
    }
    else if (strcmp(argv[i+1], "gaussian") == 0)
    {
        FilterShape = GAUSSIAN;
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
    else
    {
        if (i + 1 >= argc)
            REPORT_ERROR(3000, "FourierMask: you haven't supplied a file");
        fn_mask = argv[i+1];
        FilterShape = FilterBand = FROM_FILE;
    }

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
    if (checkParameter(argc, argv, "-sampling"))
    {
        double sampling_rate = textToFloat(getParameter(argc, argv, "-sampling"));
        if (w1 != 0)       w1 = sampling_rate / w1;
        if (w2 != 0)       w2 = sampling_rate / w2;
        /*CO: I think it is more confusing */
        //if (raised_w!=0) raised_w=sampling_rate/raised_w;
    }
}

/* Show -------------------------------------------------------------------- */
void FourierMask::show()
{
    if (FilterShape == WEDGE)
    {
        std::cout << "Missing wedge for data between tilting angles of " << w1 << " and " << w2 << " deg\n";
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
        case FROM_FILE:
            std::cout << "From file " << fn_mask << std::endl;
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
        case CTF:
            std::cout << "CTF\n" << ctf;
            break;
        case FROM_FILE:
            std::cout << "From file " << fn_mask << std::endl;
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
    << "   -fourier_mask <file>              : Provide a Fourier file\n"
    << "   -fourier_mask raised_cosine <raisedw>: Use raised cosine edges (in dig.freq.)\n"
    << "   -fourier_mask wedge <th0> <thF>   : Missing wedge for data between th0-thF \n"
    << "   -fourier_mask gaussian            : sigma=<w1>\n"
    << "   -fourier_mask ctf                 : Provide a .ctfparam file\n"
    << "  [-sampling <sampling_rate>]        : If provided pass frequencies\n"
    << "                                       are taken in Angstroms\n"
    ;
}

// Correct phase -----------------------------------------------------------
void FourierMask::correct_phase()
{
    FOR_ALL_ELEMENTS_IN_MATRIX2D(mask2D)
    if (real(mask2D(i, j)) < 0) mask2D(i, j) = -mask2D(i, j);
}

// Read mask from file -----------------------------------------------------
void FourierMask::read_mask(const FileName &fn)
{
    FilterBand = FilterShape = FROM_FILE;
    fn_mask = fn;
    if (Is_ImageXmipp(fn_mask))
    {
        std::cerr<<"Converting real-valued mask in "<<fn_mask<<" to complex values..."<<std::endl;
        ImageXmipp  I;
        I.read(fn_mask);
        FourierImageXmipp If;
        If().resize(I());
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(If())
        {
            dMij(If(),i,j) = std::complex<double>(dMij(I(),i,j),0.);
            std::cout<< dMij(If(),i,j)<<" "<<dMij(I(),i,j)<<std::endl;
        }
        mask2D = If();
        mask2D.setXmippOrigin();
    }
    else if (Is_VolumeXmipp(fn_mask))
    {
        std::cerr<<"Converting real-valued mask in "<<fn_mask<<" to complex values..."<<std::endl;
        VolumeXmipp  V;
        V.read(fn_mask);
        FourierVolumeXmipp Vf;
        Vf().resize(V());
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Vf())
        {
            dVkij(Vf(),k,i,j) = std::complex<double>(dVkij(V(),k,i,j),0.);
        }
        mask3D = Vf();
        mask3D.setXmippOrigin();
    }
    else if (Is_FourierImageXmipp(fn_mask))
    {
        FourierImageXmipp  I;
        I.read(fn_mask);
        mask2D = I();
        mask2D.setXmippOrigin();
    }
    else if (Is_FourierVolumeXmipp(fn_mask))
    {
        FourierVolumeXmipp V;
        V.read(fn_mask);
        mask3D = V();
        mask3D.setXmippOrigin();
    }
    else
    {
        REPORT_ERROR(1,"Please provide fourier mask as a (Fourier) image/volume\n");
    }
}

/* Save -------------------------------------------------------------------- */
void FourierMask::write_amplitude(const FileName &fn, int dim,
                                  bool do_not_center)
{
    Matrix1D< std::complex<double> > aux1D;
    Matrix2D< std::complex<double> > aux2D;
    Matrix3D< std::complex<double> > aux3D;
    switch (dim)
    {
    case 1:
        aux1D = mask1D;
        break;
    case 2:
        aux2D = mask2D;
        break;
    case 3:
        aux3D = mask3D;
        break;
    }
    if (!do_not_center)
    {
        switch (dim)
        {
        case 1:
            CenterFFT(aux1D, true);
            break;
        case 2:
            CenterFFT(aux2D, true);
            break;
        case 3:
            CenterFFT(aux3D, true);
            break;
        }
    }
    if (dim == 1)
    {
        Matrix1D<double> v;
        FFT_magnitude(aux1D, v);
        FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
        v(i) = log10(1 + v(i) * v(i));
        v.write(fn);
    }
    else if (dim == 2)
    {
        ImageXmipp  I;
        FFT_magnitude(aux2D, I());
        FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
        I(i, j) = log10(1 + I(i, j) * I(i, j));
        I.write(fn);
    }
    else
    {
        VolumeXmipp V;
        FFT_magnitude(aux3D, V());
        FOR_ALL_ELEMENTS_IN_MATRIX3D(V())
        V(k, i, j) = log10(1 + V(k, i, j) * V(k, i, j));
        V.write(fn);
    }
}

void FourierMask::write_mask(const FileName &fn, int dim)
{
    FourierImageXmipp  I;
    FourierVolumeXmipp V;
    switch (dim)
    {
    case 1:
        mask1D.write(fn);
        break;
    case 2:
        I() = mask2D;
        I.write(fn);
        break;
    case 3:
        V() = mask3D;
        V.write(fn);
        break;
    }
}

/* Apply mask -------------------------------------------------------------- */
void FourierMask::apply_mask_Fourier(Matrix1D< std::complex<double> > &v)
{
    v *= mask1D;
}
void FourierMask::apply_mask_Fourier(Matrix2D< std::complex<double> > &v)
{
    v *= mask2D;
}
void FourierMask::apply_mask_Fourier(Matrix3D< std::complex<double> > &v)
{
    v *= mask3D;
}

void FourierMask::apply_mask_Space(Matrix1D<double> &v)
{
    Matrix1D< std::complex<double> > aux1D;
    FourierTransform(v, aux1D);
    aux1D *= mask1D;
    InverseFourierTransform(aux1D, v);
}
void FourierMask::apply_mask_Space(Matrix2D<double> &v)
{
    Matrix2D< std::complex<double> > aux2D;
    FourierTransform(v, aux2D);
    if (XSIZE(mask2D) == 0) generate_mask(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux2D)
    DIRECT_MAT_ELEM(aux2D, i, j) *= DIRECT_MAT_ELEM(mask2D, i, j);
    InverseFourierTransform(aux2D, v);
}
void FourierMask::apply_mask_Space(Matrix3D<double> &v)
{
    Matrix3D< std::complex<double> > aux3D;
    FourierTransform(v, aux3D);
    if (XSIZE(mask3D) == 0) generate_mask(v);
    aux3D *= mask3D;
    InverseFourierTransform(aux3D, v);
}

/* Resize ------------------------------------------------------------------ */
void FourierMask::resize_mask(int Ydim, int Xdim)
{
    mask2D.selfScaleToSize(Ydim, Xdim);
}
void FourierMask::resize_mask(int Zdim, int Ydim, int Xdim)
{
    mask3D.selfScaleToSize(Zdim, Ydim, Xdim);
}

/* Mask power -------------------------------------------------------------- */
double FourierMask::mask2D_power(double wmin, double wmax)
{
    Matrix1D<int>    idx(2);
    Matrix1D<double> freq(2);
    double retval = 0, N = 0;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(mask2D)
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(mask2D, idx, freq);
        double w = freq.module();
        if (w > wmin && w < wmax)
        {
            double mag = abs(mask2D(i, j));
            retval += mag * mag;
            N++;
        }
    }

    if (N != 0) return retval / N;
    else      return 0;
}
