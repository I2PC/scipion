/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2003)
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

#include "resolution_ssnr.h"

#include <data/args.h>
#include <data/projection.h>
#include <data/fftw.h>
#include <data/metadata_extension.h>

// Read parameters from command line ---------------------------------------
void Prog_SSNR_prm::read(int argc, char **argv)
{
    radial_avg = checkParameter(argc, argv, "-radial_avg");
    if (!radial_avg)
    {
        fn_S = getParameter(argc, argv, "-S");
        fn_N = getParameter(argc, argv, "-N");
        fn_Ssel = getParameter(argc, argv, "-selS");
        fn_Nsel = getParameter(argc, argv, "-selN");
        generate_VSSNR = checkParameter(argc, argv, "-generate_VSSNR");
        if (generate_VSSNR)
        {
            fn_VSSNR = getParameter(argc, argv, "-VSSNR");
            fn_out_images = getParameter(argc, argv, "-oimages", "");
        }
    }
    else
        fn_VSSNR = getParameter(argc, argv, "-VSSNR");
    ring_width = textToFloat(getParameter(argc, argv, "-ring", "4"));
    Tm = textToFloat(getParameter(argc, argv, "-sampling_rate", "1"));
    min_power = textToFloat(getParameter(argc, argv, "-min_power", "1e-10"));
    fn_out = getParameter(argc, argv, "-o", "");
}

// Show parameters ---------------------------------------------------------
std::ostream & operator << (std::ostream &out, const Prog_SSNR_prm &prm)
{
    out << "Signal:         " << prm.fn_S       << std::endl
    << "Noise:          " << prm.fn_N       << std::endl
    << "Signal selfile: " << prm.fn_Ssel    << std::endl
    << "Noise  selfile: " << prm.fn_Nsel    << std::endl
    << "Volumetric SSNR:" << prm.fn_VSSNR   << std::endl
    << "Output images:  " << prm.fn_out     << std::endl
    << "Ring width:     " << prm.ring_width << std::endl
    << "Sampling rate:  " << prm.Tm         << std::endl
    << "Generate VSSNR: " << prm.generate_VSSNR << std::endl
    << "Radial average: " << prm.radial_avg << std::endl
    ;
    return out;
}

// Usage -------------------------------------------------------------------
void Prog_SSNR_prm::usage() const
{
    std::cerr << " SSNR Estimation ------------------------------------------------\n"
    << "SSNR\n"
    << "   -S <Volume|Selfile>   : Signal volume or its projections\n"
    << "   -N <Volume|Selfile>   : Noise volume or its projections\n"
    << "   -selS <Selfile>       : Selfile with experimental images\n"
    << "   -selN <Selfile>       : Selfile with noise images\n"
    << "  [-ring <w=4>]          : Ring width for the SSNR computation\n"
    << "  [-sampling_rate <Tm=1>]: Sampling rate A/pixel\n"
    << "  [-o <SSNR file=\"\">]    : Output file\n"
    << " VSSNR Estimation -----------------------------------------------\n"
    << "SSNR\n"
    << "   -generate_VSSNR       : generate VSSNR\n"
    << "   -S <Volume|Selfile>   : Signal volume or its projections\n"
    << "   -N <Volume|Selfile>   : Noise volume or its projections\n"
    << "   -selS <Selfile>       : Selfile with experimental images\n"
    << "   -selN <Selfile>       : Selfile with noise images\n"
    << "   -VSSNR <fn_vol>       : Volume with the Volumetric SSNR\n"
    << "  [-oimages <fn_root>]   : Root name for individual SSNR estimations\n"
    << "  [-ring <w=4>]          : Ring width for the SSNR computation\n"
    << "  [-min_power <th=1e-10>]: Minimum power\n"
    << "  [-sampling_rate <Tm=1>]: Sampling rate A/pixel\n"
    << " Estimation by radial averaging ---------------------------------\n"
    << " SSNR -radial_avg\n"
    << "   -VSSNR <fn_vol>       : Volume with the Volumetric SSNR\n"
    << "  [-ring <w=4>]          : Ring width for the SSNR computation\n"
    << "  [-sampling_rate <Tm=1>]: Sampling rate A/pixel\n"
    << "  [-min_power <th=1e-10>]: Minimum power\n"
    << "  [-o <SSNR file=\"\">]    : Output file\n"
    ;
}

// Produce side Info -------------------------------------------------------
void Prog_SSNR_prm::produce_side_info()
{
    if (!radial_avg)
    {
        if (S.isImage(fn_S))
        {
            S.readApplyGeo(fn_S);
            S().setXmippOrigin();
            N.readApplyGeo(fn_N);
            N().setXmippOrigin();
            if (!S().sameShape(N()))
                REPORT_ERROR(ERR_MULTIDIM_SIZE,
                             "SSNR: Signal and Noise volumes are not of the same size");
        }
        else
        {
            SF_Sth.read(fn_S);
            SF_Nth.read(fn_N);
            if (SF_Sth.size() != SF_Nth.size())
                REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                             "SSNR: the number of projections in both selfiles is different");
        }

        SF_S.read(fn_Ssel);
        SF_N.read(fn_Nsel);
        int sYdim, sXdim;
        ImgSize(SF_S, sXdim, sYdim);
        int nYdim, nXdim;
        ImgSize(SF_N, nXdim, nYdim);
        if (sYdim != nYdim || sYdim != YSIZE(S()) ||
            sXdim != nXdim || sXdim != XSIZE(S()))
            REPORT_ERROR(ERR_MULTIDIM_SIZE,
                         "SSNR: conflict among the projection/projection sizes "
                         "or projection/volume");

        if (SF_S.size() != SF_N.size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                         "SSNR: the number of projections in both selfiles is different");
        if (XSIZE(S()) == 0 && SF_Sth.size() != SF_S.size())
            REPORT_ERROR(ERR_MULTIDIM_SIZE,
                         "SSNR: the number of projections in both selfiles is different");

        if (fn_out_images == "")
            fn_out_images = "individualSSNR";
    }
    else
    {
        VSSNR.read(fn_VSSNR);
        VSSNR().setXmippOrigin();
    }
}

// Estimate SSNR ----------------------------------------------------------
//#define DEBUG
void Prog_SSNR_prm::Estimate_SSNR(int dim, Matrix2D<double> &output)
{
    // These vectors are for 1D
    Matrix1D<double> S_S21D((int)(XSIZE(S()) / 2 - ring_width)),
    S_N21D((int)(XSIZE(S()) / 2 - ring_width)),
    K1D((int)(XSIZE(S()) / 2 - ring_width)),
    S_SSNR1D;
    Matrix1D<double> N_S21D((int)(XSIZE(S()) / 2 - ring_width)),
    N_N21D((int)(XSIZE(S()) / 2 - ring_width)),
    N_SSNR1D;

    // Selfile of the 2D images
    MetaData SF_individual;

    std::cerr << "Computing the SSNR ...\n";
    init_progress_bar(SF_S.size());
    int imgno = 0;
    FOR_ALL_OBJECTS_IN_METADATA(SF_S)
    {
        Image<double> Is, In, Inp;
        FileName fn_img;
        SF_S.getValue(MDL_IMAGE,fn_img);
        Is.readApplyGeo(fn_img,SF_S,objId);
        Is().setXmippOrigin();
        SF_N.getValue(MDL_IMAGE,fn_img);
        In.readApplyGeo(fn_img,SF_N,objId);
        In().setXmippOrigin();
        Inp() = In();

        Projection Iths, Ithn;
        if (XSIZE(S()) != 0)
        {
            project_Volume(S(), Iths, YSIZE(Is()), XSIZE(Is()),
                           Is.rot(), Is.tilt(), Is.psi());
            project_Volume(N(), Ithn, YSIZE(Is()), XSIZE(Is()),
                           Is.rot(), Is.tilt(), Is.psi());
        }
        else
        {
            FileName fn_img;
            SF_Sth.getValue(MDL_IMAGE,fn_img);
            Iths.read(fn_img);
            SF_Nth.getValue(MDL_IMAGE,fn_img);
            Ithn.read(fn_img);
        }

#ifdef DEBUG
        Image<double> save;
        save() = Is();
        save.write("PPPread_signal.xmp");
        save() = In();
        save.write("PPPread_noise.xmp");
        save() = Iths();
        save.write("PPPtheo_signal.xmp");
        save() = Ithn();
        save.write("PPPtheo_noise.xmp");
#endif

        Is() -= Iths();
        In() -= Ithn();

        MultidimArray< std::complex<double> > FFT_Is;
        FourierTransform(Is(), FFT_Is);
        MultidimArray< std::complex<double> > FFT_Iths;
        FourierTransform(Iths(), FFT_Iths);
        MultidimArray< std::complex<double> > FFT_In;
        FourierTransform(In(), FFT_In);
        MultidimArray< std::complex<double> > FFT_Ithn;
        FourierTransform(Ithn(), FFT_Ithn);

#ifdef DEBUG

        Image< std::complex<double> > savec;
        savec() = FFT_Is;
        savec.write("PPPFFTread_signal.xmp");
        savec() = FFT_In;
        savec.write("PPPFFTread_noise.xmp");
        savec() = FFT_Iths;
        savec.write("PPPFFTtheo_signal.xmp");
        savec() = FFT_Ithn;
        savec.write("PPPFFTtheo_noise.xmp");
#endif

        // Compute the amplitudes
        Image<double> S2s;
        FFT_magnitude(FFT_Iths, S2s());
        S2s() *= S2s();
        Image<double> N2s;
        FFT_magnitude(FFT_Is  , N2s());
        N2s() *= N2s();
        Image<double> S2n;
        FFT_magnitude(FFT_Ithn, S2n());
        S2n() *= S2n();
        Image<double> N2n;
        FFT_magnitude(FFT_In  , N2n());
        N2n() *= N2n();

#ifdef DEBUG

        save() = S2s();
        save.write("PPPS2s.xmp");
        save() = N2s();
        save.write("PPPN2s.xmp");
        save() = S2n();
        save.write("PPPS2n.xmp");
        save() = N2n();
        save.write("PPPN2n.xmp");
#endif

        if (dim == 2)
        {
            // Compute the SSNR image
        	Image<double> SSNR2D;
            SSNR2D().initZeros(S2s());
            FOR_ALL_ELEMENTS_IN_ARRAY2D(S2s())
            {
                double ISSNR = 0, alpha = 0, SSNR = 0;
                if (N2s(i, j) > min_power)
                    ISSNR = S2s(i, j) / N2s(i, j);
                if (N2n(i, j) > min_power)
                    alpha = S2n(i, j) / N2n(i, j);
                if (alpha   > min_power)
                    SSNR = XMIPP_MAX(ISSNR / alpha - 1, 0);
                if (SSNR    > min_power)
                    SSNR2D(i, j) = 10 * log10(SSNR + 1);
            }
            CenterFFT(SSNR2D(), true);
#ifdef DEBUG

            save() = SSNR2D();
            save.write("PPPSSNR2D.xmp");
#endif

            // Save image
            FileName fn_img_out = fn_out_images + integerToString(Is.name().getNumber(), 5) + ".xmp";
            SSNR2D.setEulerAngles(Is.rot(),Is.tilt(),Is.psi());
            SSNR2D.write(fn_img_out);
            SF_individual.addObject();
            SF_individual.setValue(MDL_IMAGE,fn_img_out);
        }

        // Average over rings
        Matrix1D<int>    idx(2);
        Matrix1D<double> freq(2);

        STARTINGX(S2s()) = STARTINGY(S2s()) = 0;
        STARTINGX(N2s()) = STARTINGY(N2s()) = 0;
        STARTINGX(S2n()) = STARTINGY(S2n()) = 0;
        STARTINGX(N2n()) = STARTINGY(N2n()) = 0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(S2s())
        {
            XX(idx) = j;
            YY(idx) = i;
            FFT_idx2digfreq(FFT_Is, idx, freq);
            if (XX(freq) < 0)
                continue;

            // Look for the index corresponding to this frequency
            double w = freq.module();
            double widx = w * XSIZE(S());
            if (widx >= VEC_XSIZE(S_S21D))
                continue;
            int l0 = XMIPP_MAX(0, CEIL(widx - ring_width));
            int lF = FLOOR(widx);

            double S_signal = S2s(i, j);
            double S_noise = N2s(i, j);
            double N_signal = S2n(i, j);
            double N_noise = N2n(i, j);
            for (int l = l0; l <= lF; l++)
            {
                S_S21D(l) += S_signal;
                S_N21D(l) += S_noise;
                N_S21D(l) += N_signal;
                N_N21D(l) += N_noise;
                K1D(l)++;
            }
        }

        // Finished with this image
        if (++imgno % 50 == 0)
            progress_bar(imgno);
    }
    progress_bar(SF_S.size());

    // Compute the SSNR
    S_SSNR1D = S_S21D / S_N21D;
    S_S21D /= K1D;
    S_N21D /= K1D;
    N_SSNR1D = N_S21D / N_N21D;
    N_S21D /= K1D;
    N_N21D /= K1D;

    output.resize(VEC_XSIZE(S_SSNR1D), 9);
    int imax = 0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(S_SSNR1D)
    {
        double w;
        FFT_IDX2DIGFREQ(i, XSIZE(S()), w);
        if (w < 0)
        {
            imax = i;
            break;
        }
        output(i, 0) = i;
        output(i, 1) = w * 1 / Tm;
        double SSNR = S_SSNR1D(i) / N_SSNR1D(i);
        if (SSNR > 1)
            output(i, 2) = 10 * log10(SSNR - 1); // Corrected SSNR
        else
            output(i, 2) = -1000;         // In fact it should be -inf
        output(i, 3) = S_SSNR1D(i);
        output(i, 4) = 10 * log10(S_S21D(i) / imgno);
        output(i, 5) = 10 * log10(S_N21D(i) / imgno);
        output(i, 6) = N_SSNR1D(i);
        output(i, 7) = 10 * log10(N_S21D(i) / imgno);
        output(i, 8) = 10 * log10(N_N21D(i) / imgno);
        imax++;
    }
    output.resize(imax, 9);

    // Produce VSSNR ........................................................
    if (dim == 2)
    {
        std::cerr << "Interpolating the VSSNR ...\n";
        SF_individual.write(fn_out_images + ".sel");
        system(((std::string)"xmipp_art -i " + fn_out_images + ".sel -o " + fn_VSSNR +
                " -l 0.1 -R " + integerToString(ROUND(XSIZE(S()) / 3)) + " -ray_length 1 -n 5").c_str());
        system(((std::string)"xmipp_rmsel " + fn_out_images + ".sel ").c_str());
    }
}
#undef DEBUG

// Estimate radial average ---------------------------------------------------
void Prog_SSNR_prm::Radial_average(Matrix2D<double> &output)
{
    // Compute the radial average ...........................................
    Matrix1D<double> VSSNR_avg((int)(XSIZE(VSSNR()) / 2 - ring_width));
    Matrix1D<double> K1D(VSSNR_avg);
    Matrix1D<int>    idx(3);
    Matrix1D<double> freq(3);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(VSSNR())
    {
        VECTOR_R3(idx, j, i, k);
        FFT_idx2digfreq(VSSNR(), idx, freq);
        if (XX(freq) < 0)
            continue;

        // Look for the index corresponding to this frequency
        double w = freq.module();
        double widx = w * XSIZE(VSSNR());
        if (widx >= VEC_XSIZE(VSSNR_avg))
            continue;
        int l0 = XMIPP_MAX(0, CEIL(widx - ring_width));
        int lF = FLOOR(widx);

        double VSSNRkij = pow(10, VSSNR(k, i, j) / 10) - 1;

        for (int l = l0; l <= lF; l++)
        {
            VSSNR_avg(l) += VSSNRkij;
            K1D(l)++;
        }
    }

    FOR_ALL_ELEMENTS_IN_MATRIX1D(VSSNR_avg)
    if (K1D(i) != 0)
        VSSNR_avg(i) /= K1D(i);

    // Produce output .......................................................
    output.resize(VEC_XSIZE(VSSNR_avg), 3);
    for (int i = 0; i < VEC_XSIZE(VSSNR_avg); i++)
    {
        double w;
        FFT_IDX2DIGFREQ(i, XSIZE(VSSNR()), w);
        output(i, 0) = i;
        output(i, 1) = w * 1 / Tm;
        double SSNR = VSSNR_avg(i);
        if (SSNR > 1)
            output(i, 2) = 10 * log10(SSNR - 1); // Corrected SSNR
        else
            output(i, 2) = -1000;         // In fact it should be -inf
    }
}

// Main --------------------------------------------------------------------
void ROUT_SSNR(Prog_SSNR_prm &prm, Matrix2D<double> &output)
{
    std::cout << prm;
    prm.produce_side_info();
    if (!prm.radial_avg)
    {
        if (!prm.generate_VSSNR)
            prm.Estimate_SSNR(1, output);
        else
            prm.Estimate_SSNR(2, output);
        if (prm.fn_out != "")
            output.write(prm.fn_out);
        else
        {
            FileName fn_out = prm.fn_S.insertBeforeExtension("_SSNR");
            fn_out.substituteExtension("vol", "txt");
            output.write(fn_out);
        }
    }
    else
    {
        prm.Radial_average(output);
        if (prm.fn_out != "")
            output.write(prm.fn_out);
        else
            output.write(prm.fn_VSSNR.insertBeforeExtension("_radial_avg"));
    }
}
