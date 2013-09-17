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
#include <data/xmipp_fftw.h>
#include <data/metadata_extension.h>
#include "reconstruct_art.h"


void ProgSSNR::defineParams()
{
    addUsageLine("Evaluate the reconstruction quality using SSNR or its volumetric distribution VSSNR.");
    addUsageLine("+You must produce a reconstruction with gaussian white noise with exactly the same angles and");
    addUsageLine("+reconstruction parameters as the signal reconstruction. Our ART program can produce this noisy");
    addUsageLine("+reconstruction. If you don't use ART as reconstruction algorithm, you will have to produce the");
    addUsageLine("+noisy reconstruction by yourself (see our demo in order to learn how to do this). Finally, the");
    addUsageLine("+VSSNR can be visualized using any volume viewer.");
    addUsageLine("+++ %BR% %BR%");
    addUsageLine("+++* Do not use any mask during the reconstruction since the resolution estimate is biased");
    addUsageLine("+++* Before applying this measure make sure that the noise images and the background in the experimental" \
                 " images have zero average. It is also assumed that the projection of the reconstructed volumes matches the" \
                 " experimental or the noisy projections (i.e., check that both projections are within the same range");
    addUsageLine("+++* The VSSNR is stored in the file as 10*log10(1+VSSNR). Thus, if you want a threshold of VSSNR=1, the" \
                 " threshold in the visualizer must be set to =10*log10(1+1)=3=");
    addUsageLine("+++* If the reconstruction algorithm linearly scales the reconstructed volume so that the reprojected" \
                 " images do not match the gray values of the experimental images, use adjust_volume to correct for the linear" \
                 " transformation before computing the SSNR images");

    addSeeAlsoLine("resolution_fsc");

    addParamsLine("== SSNR or VSSNR Estimation ==");
    addParamsLine("   [--signal <signal_file>]  : Signal volume");
    addParamsLine("  alias -S;");
    addParamsLine("   requires --noise;");
    addParamsLine("   [--noise <noise_file>]    : Noise volume");
    addParamsLine("  alias -N;");
    addParamsLine("   [--sel_signal <signal_mdfile>] : Metadata file with the images used for the signal reconstruction");
    addParamsLine("  alias -selS;");
    addParamsLine("   [--sel_noise  <noise_mdfile>]  : Metadata file with the images used for the noise reconstruction");
    addParamsLine("  alias -selN;");
    addParamsLine("   [-o <SSNR_file=\"\">]     : Output file with the SSNR estimation");
    addParamsLine("                             :+++If the output filename is not given, then the input filename is taken");
    addParamsLine("                             :+++from the -S parameter by inserting _SSNR before the extension.");
    addParamsLine("                             :+++The columns of the output file are as follows: %BR%");
    addParamsLine("                             :+++  =Column 1: Fourier pixel number= %BR%");
    addParamsLine("                             :+++  =Column 2: Corresponding frequency (1/Angstroms)= %BR%");
    addParamsLine("                             :+++  =Column 3: 10*log10(SSNR-1)  (SSNR=ISSNR/alpha)= %BR%");
    addParamsLine("                             :+++  =Column 4: ISSNR    (Input SSNR=Sexp/Nexp)= %BR%");
    addParamsLine("                             :+++  =Column 5: 10*log10(Sexp)= %BR%");
    addParamsLine("                             :+++  =Column 6: 10*log10(Nexp)= %BR%");
    addParamsLine("                             :+++  =Column 7: alpha       (alpha=Snoise/Nnoise)= %BR%");
    addParamsLine("                             :+++  =Column 8: 10*log10(Snoise)= %BR%");
    addParamsLine("                             :+++  =Column 9: 10*log10(Nnoise)= ");

    addParamsLine("   [--ring <w=4>]            : Ring width for the SSNR averaging computation (Measured in Fourier pixels)");
    addParamsLine("   [--sampling_rate <Ts=1>]  : Pixel size (Angstrom)");
    addParamsLine("  alias -s;");
    addParamsLine("   [--min_power <th=1e-10>]  : Minimum power in Fourier space. If at any moment, the SSNR must divide by something");
    addParamsLine("                             : smaller than this value, then the SSNR at that frequency is assumed to be 0. High ");
    addParamsLine("                             : values of this threshold result on speckled images. If this is the case, lower ");
    addParamsLine("                             : this threshold (it must always be positive) ");
    addParamsLine("   [--gen_VSSNR]             : (Allowed global options: --ring, --sampling_rate, --min_power)");
    addParamsLine("                             :+++ %BR%");
    addParamsLine("                             : Generate the individual estimations of the SSNR for each particle and build an ");
    addParamsLine("                             : interpolation volume (VSSNR) that is compatible with the individual SSNRs.");
    addParamsLine("                             :+ In fact, the VSSNR is stored as 10*log10(1+SSNR). Thus, after interpolating the threshold");
    addParamsLine("                             :+ at SSNR = 1 must be shown as the threshold of the output interpolated volume at 3.01 ");
    addParamsLine("                             :+ (10*log10(1+1)). The threshold at SSNR = 4 must be shown as the threshold of the output");
    addParamsLine("                             :+ interpolated volume at 6.99 (10*log10(1+4)). The 1D SSNR is also generated as a side-product");
    addParamsLine("   requires --signal,--noise,--VSSNR;");
    addParamsLine("   [--VSSNR <fn_vol_file>]   : Volume with the Volumetric SSNR");
    addParamsLine("   [--oroot <root=\"\">]     : Root name for individual SSNR estimations");
    addParamsLine("== Estimation by radial averaging of the VSSNR ==");
    addParamsLine("   [--radial_avg]            : Do radial averaging estimation");
    addParamsLine("   : (Allowed global options: --ring, --sampling_rate, --min_power, -o)");
    addParamsLine("   requires --VSSNR;");

    addExampleLine("SSNR Resolution of a reconstructed volume:", false);
    addExampleLine("xmipp_resolution_ssnr -S recFourierSignal.vol -N recFourierNoise.vol -selS projections.xmd  -selN noiseProjections.xmd");
    addExampleLine("VSSNR Resolution of a reconstructed volume:", false);
    addExampleLine("xmipp_resolution_ssnr -S recFourierSignal.vol -N recFourierNoise.vol -selS projections.xmd  -selN noiseProjections.xmd --gen_VSSNR --VSSNR volumeout.vol");
}


void ProgSSNR::readParams()
{
    radial_avg = checkParam("--radial_avg");
    if (!radial_avg)
    {
        fn_S = getParam("-S");
        fn_N = getParam("-N");
        fn_S_sel = getParam("--sel_signal");
        fn_N_sel = getParam("--sel_noise");
        generate_VSSNR = checkParam("--gen_VSSNR");
        if (generate_VSSNR)
        {
            fn_VSSNR = getParam("--VSSNR");
            fn_out_images = getParam("--oroot");
        }
    }
    else
        fn_VSSNR = getParam("--VSSNR");

    ring_width = getDoubleParam("--ring");
    Tm = getDoubleParam("--sampling_rate");
    min_power = getDoubleParam("--min_power");
    fn_out = getParam("-o");
}

void ProgSSNR::show()
{
    std::cout
    << "Signal Volume:     " << fn_S       << std::endl
    << "Noise  Volume:     " << fn_N       << std::endl
    << "Signal selfile:    " << fn_S_sel   << std::endl
    << "Noise  selfile:    " << fn_N_sel   << std::endl
    << "Volumetric SSNR:   " << fn_VSSNR   << std::endl
    << "Output images:     " << fn_out     << std::endl
    << "Ring width:        " << ring_width << std::endl
    << "Sampling rate:     " << Tm         << std::endl
    << "Generate VSSNR:    " << generate_VSSNR << std::endl
    << "Radial average:    " << radial_avg << std::endl;
}

void ProgSSNR::produceSideInfo()
{
    if (!radial_avg)
    {
        S.read(fn_S);
        S().setXmippOrigin();
        N.read(fn_N);
        N().setXmippOrigin();
        if (!S().sameShape(N()))
            REPORT_ERROR(ERR_MULTIDIM_SIZE,
                         "SSNR: Signal and Noise volumes are not of the same size");

        SF_S.read(fn_S_sel);
        SF_N.read(fn_N_sel);

        if (fn_out_images == "")
            fn_out_images = "individualSSNR";
    }
    else
    {
        VSSNR.read(fn_VSSNR);
        VSSNR().setXmippOrigin();
    }
}

//#define DEBUG
void ProgSSNR::run()
{
    show();
    produceSideInfo();

    Matrix2D<double> output;
    if (!radial_avg)
    {
        if (!generate_VSSNR)
            estimateSSNR(1, output);
        else
            estimateSSNR(2, output);
        if (fn_out == "")
        	fn_out=fn_S.insertBeforeExtension("_SSNR").removeLastExtension().addExtension("xmd");
    }
    else
    {
        radialAverage(output);
        if (fn_out == "")
        	fn_out=fn_VSSNR.insertBeforeExtension("_radial_avg").removeLastExtension().addExtension("xmd");
    }
#ifdef DEBUG
    output.write(fn_out);
#endif
    MetaData MD;
    for (size_t i=1; i<MAT_YSIZE(output); ++i)
    {
    	size_t id=MD.addObject();
    	MD.setValue(MDL_RESOLUTION_FREQ,output(i,1),id);
    	if (output(i,2)<-900)
    		MD.setValue(MDL_RESOLUTION_SSNR,0.0,id);
    	else
    		MD.setValue(MDL_RESOLUTION_SSNR,pow(10.0,output(i,2)/10.0),id);
    	MD.setValue(MDL_RESOLUTION_FREQREAL,1.0/output(i,1),id);
    }
    MD.write(fn_out);
}

void ProgSSNR::estimateSSNR(int dim, Matrix2D<double> &output)
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
    int imgno = 1;
    Image<double> Is, In;
    Projection Iths, Ithn;
    MultidimArray< std::complex<double> > FFT_Is, FFT_Iths,  FFT_In, FFT_Ithn;
    MultidimArray<double> S2s, N2s, S2n, N2n;
    FileName fn_img;
    FOR_ALL_OBJECTS_IN_METADATA2(SF_S, SF_N)
    {
        Is.readApplyGeo(SF_S,__iter.objId);
        Is().setXmippOrigin();
        In.readApplyGeo(SF_N,__iter2.objId);
        In().setXmippOrigin();

        if (!XMIPP_EQUAL_REAL(Is.rot(), In.rot()) ||
            !XMIPP_EQUAL_REAL(Is.tilt(), In.tilt()) ||
            !XMIPP_EQUAL_REAL(Is.psi(), In.psi()))
            REPORT_ERROR(ERR_VALUE_INCORRECT, "Noise projection angles are different from signal projections angles.");

        projectVolume(S(), Iths, YSIZE(Is()), XSIZE(Is()),
                      Is.rot(), Is.tilt(), Is.psi());
        projectVolume(N(), Ithn, YSIZE(Is()), XSIZE(Is()),
                      Is.rot(), Is.tilt(), Is.psi());

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

        FourierTransformer FT(FFTW_BACKWARD);

        if (dim == 2)
        {
            FT.completeFourierTransform(Is(), FFT_Is);
            FT.completeFourierTransform(Iths(), FFT_Iths);
            FT.completeFourierTransform(In(), FFT_In);
            FT.completeFourierTransform(Ithn(), FFT_Ithn);
        }
        else
        {
            FT.FourierTransform(Is(), FFT_Is);
            FT.FourierTransform(Iths(), FFT_Iths);
            FT.FourierTransform(In(), FFT_In);
            FT.FourierTransform(Ithn(), FFT_Ithn);
        }

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
        S2s.resizeNoCopy(FFT_Iths);
        N2s.resizeNoCopy(FFT_Iths);
        S2n.resizeNoCopy(FFT_Iths);
        N2n.resizeNoCopy(FFT_Iths);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FFT_Iths)
        {
            DIRECT_MULTIDIM_ELEM(S2s, n) = abs(DIRECT_MULTIDIM_ELEM(FFT_Iths, n));
            DIRECT_MULTIDIM_ELEM(S2s, n) *= DIRECT_MULTIDIM_ELEM(S2s, n);
            DIRECT_MULTIDIM_ELEM(N2s, n) = abs(DIRECT_MULTIDIM_ELEM(FFT_Is, n));
            DIRECT_MULTIDIM_ELEM(N2s, n) *= DIRECT_MULTIDIM_ELEM(N2s, n);
            DIRECT_MULTIDIM_ELEM(S2n, n) = abs(DIRECT_MULTIDIM_ELEM(FFT_Ithn, n));
            DIRECT_MULTIDIM_ELEM(S2n, n) *= DIRECT_MULTIDIM_ELEM(S2n, n);
            DIRECT_MULTIDIM_ELEM(N2n, n) = abs(DIRECT_MULTIDIM_ELEM(FFT_In, n));
            DIRECT_MULTIDIM_ELEM(N2n, n) *= DIRECT_MULTIDIM_ELEM(N2n, n);
        }

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
            SSNR2D().initZeros(S2s);
            const MultidimArray<double> & SSNR2Dmatrix=SSNR2D();
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(S2s)
            {
                double ISSNR = 0, alpha = 0, SSNR = 0;
                double aux = DIRECT_MULTIDIM_ELEM(N2s,n);
                if (aux > min_power)
                    ISSNR = DIRECT_MULTIDIM_ELEM(S2s,n) / aux;
                aux = DIRECT_MULTIDIM_ELEM(N2n,n);
                if (aux > min_power)
                    alpha = DIRECT_MULTIDIM_ELEM(S2n,n) / aux;
                if (alpha > min_power)
                {
                    aux = ISSNR / alpha - 1.0;
                    SSNR = XMIPP_MAX(aux, 0.0);
                }
                if (SSNR    > min_power)
                    DIRECT_MULTIDIM_ELEM(SSNR2Dmatrix,n) = 10.0 * log10(SSNR + 1.0);
            }
            CenterFFT(SSNR2D(), true);
#ifdef DEBUG

            save() = SSNR2Dmatrix;
            save.write("PPPSSNR2D.xmp");
#endif

            // Save image
            FileName fn_img_out;
            fn_img_out.compose(imgno, fn_out_images, "stk");
            SSNR2D.setEulerAngles(Is.rot(),Is.tilt(),Is.psi());
            SSNR2D.write(fn_img_out);
            size_t objId = SF_individual.addObject();
            SF_individual.setValue(MDL_IMAGE,fn_img_out,objId);
        }

        // Average over rings
        Matrix1D<int>    idx(2);
        Matrix1D<double> freq(2);

        STARTINGX(S2s) = STARTINGY(S2s) = 0;
        STARTINGX(N2s) = STARTINGY(N2s) = 0;
        STARTINGX(S2n) = STARTINGY(S2n) = 0;
        STARTINGX(N2n) = STARTINGY(N2n) = 0;
        double Xdim=XSIZE(S());
        double maxwidx=VEC_XSIZE(S_S21D);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(S2s)
        {
            XX(idx) = j;
            YY(idx) = i;
            FFT_idx2digfreq(FFT_Is, idx, freq);
            if (XX(freq) < 0)
                continue;

            // Look for the index corresponding to this frequency
            double w = freq.module();
            double widx = w * Xdim;
            if (widx >= maxwidx)
                continue;
            int l0 = CEIL(widx - ring_width);
            l0=XMIPP_MAX(0, l0);
            int lF = FLOOR(widx);

            double S_signal = A2D_ELEM(S2s, i, j);
            double S_noise = A2D_ELEM(N2s, i, j);
            double N_signal = A2D_ELEM(S2n, i, j);
            double N_noise = A2D_ELEM(N2n, i, j);
            for (int l = l0; l <= lF; l++)
            {
                VEC_ELEM(S_S21D,l) += S_signal;
                VEC_ELEM(S_N21D,l) += S_noise;
                VEC_ELEM(N_S21D,l) += N_signal;
                VEC_ELEM(N_N21D,l) += N_noise;
                VEC_ELEM(K1D,l)++;
            }
        }

        // Finished with this image
        if (++imgno % 50 == 0)
            progress_bar(imgno);
    }
    progress_bar(SF_S.size());

    // Compute the SSNR
    S_SSNR1D = S_S21D / S_N21D;
    S_S21D *= 1.0/K1D;
    S_N21D *= 1.0/K1D;
    N_SSNR1D = N_S21D / N_N21D;
    N_S21D *= 1.0/K1D;
    N_N21D *= 1.0/K1D;

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
        SF_individual.write(fn_out_images + ".xmd");

        ProgReconsART artRecons;
        artRecons.read(formatString("-i %s.xmd -o %s -l 0.1 -R %i --ray_length 1 -n 5",fn_out_images.c_str(),
                                    fn_VSSNR.c_str(), ROUND(XSIZE(S()) / 3)));
        artRecons.run();

        remove(fn_out_images.addExtension("xmd").c_str());
        remove(fn_out_images.addExtension("stk").c_str());
    }
}

// Estimate radial average ---------------------------------------------------
void ProgSSNR::radialAverage(Matrix2D<double> &output)
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
    for (size_t i = 0; i < VEC_XSIZE(VSSNR_avg); i++)
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
