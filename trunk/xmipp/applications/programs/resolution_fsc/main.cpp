/***************************************************************************
 *
 * Authors:     Sjors Scheres (scheres@cnb.uam.es)
 *              Roberto Marabini
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

#include <data/progs.h>
#include <data/args.h>
#include <data/fft.h>
#include <data/selfile.h>

class Resolution_parameters: public Prog_parameters
{
public:
    FileName    fn_ref, fn_root;
    ImageXmipp  refI;
    VolumeXmipp refV;
    float       sam;
    bool        do_phase_residual;

public:
    Resolution_parameters()
    {}
    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        fn_ref = getParameter(argc, argv, "-ref");
        sam = textToFloat(getParameter(argc, argv, "-sam"));
        do_phase_residual = checkParameter(argc, argv, "-phase_residual");
        if (Is_ImageXmipp(fn_ref))
        {
            refI.read(fn_ref, false, false, apply_geo);
            refI().setXmippOrigin();
        }
        else if (Is_VolumeXmipp(fn_ref))
        {
            refV.read(fn_ref);
            refV().setXmippOrigin();
        }
        else exit(0);
    }
    void show()
    {
        Prog_parameters::show();
        std::cout << "Reference file = " << fn_ref << std::endl;
        std::cout << "Sampling rate  = " << sam << std::endl;
    }
    void usage()
    {
        std::cerr << " EITHER:\n";
        std::cerr << "   -ref <input file>        : Filename for reference image/volume \n";
        std::cerr << "   -i <input file>          : either an image/volume or a selection file\n";
        std::cerr << " OR:\n";
        std::cerr << "   -set_of_images <selfile> : SelFile containing a set of 2D-images\n";
        std::cerr << " For both modes:\n";
        std::cerr << "   -sam <sampling rate>     : i.e. pixel size (in Angstrom) \n";
        std::cerr << "  [-phase_residual]         : Calculate phase residual instead of FRC \n";
        std::cerr << "  [-dont_apply_geo]         : for 2D-images: do not apply transformation stored in the header\n";
    }
};

void process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Resolution_parameters *eprm = (Resolution_parameters *) prm;

    Matrix1D<double> freq, frc, dpr, frc_noise;

    fourier_ring_correlation(eprm->refI(), img(), eprm->sam, freq, frc, frc_noise);
    differential_phase_residual(eprm->refI(), img(), eprm->sam, freq, dpr);

    // Write output
    FileName fn_dpr, fn_frc;
    fn_dpr = img.name() + ".dpr";
    fn_frc = img.name() + ".frc";
    std::ofstream out(fn_dpr.c_str(), std::ios::out);
    std::ofstream out2(fn_frc.c_str(), std::ios::out);
    out  << "# Resol. [1/Ang]   DPR [deg]" << std::endl;
    out2 << "# Resol. [1/Ang]      FRC      FRC_random_noise     Resol. [Ang]" << std::endl;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(freq)
    {
        if (i > 0)
        {
            out.width(10);
            out  << VEC_ELEM(freq, i);
            out.width(17);
            out << VEC_ELEM(dpr, i)  << std::endl;
            out2.width(10);
            out2  << VEC_ELEM(freq, i);
            out2.width(17);
            out2  << VEC_ELEM(frc, i);
            out2.width(17);
            out2  << VEC_ELEM(frc_noise, i);
            out2.width(17);
            out2  << 1./VEC_ELEM(freq, i) << std::endl;
        }
    }
    out.close();
    out2.close();
}

void process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Resolution_parameters *eprm = (Resolution_parameters *) prm;

    Matrix1D<double> freq, frc, dpr, frc_noise;
    FileName fn_dpr, fn_frc;

    if (eprm->do_phase_residual)
    {
        differential_phase_residual(eprm->refV(), vol(), eprm->sam, freq, dpr);
        fn_dpr = vol.name() + ".dpr";
        std::ofstream out(fn_dpr.c_str(), std::ios::out);
        out  << "# Resol. [1/Ang]   DPR [deg]" << std::endl;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(freq)
        {
            if (i > 0)
            {
                out.width(10);
                out  << VEC_ELEM(freq, i);
                out.width(17);
                out << VEC_ELEM(dpr, i)  << std::endl;
            }
        }
        out.close();
    }
    else 
    {
        fourier_ring_correlation(eprm->refV(), vol(), eprm->sam, freq, frc, frc_noise);
        fn_frc = vol.name() + ".frc";
        std::ofstream out2(fn_frc.c_str(), std::ios::out);
        out2 << "# Resol. [1/Ang]      FSC      FSC_random_noise     Resol. [Ang]" << std::endl;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(freq)
        {
            if (i > 0)
            {
                out2.width(10);
                out2  << VEC_ELEM(freq, i);
                out2.width(17);
                out2  << VEC_ELEM(frc, i);
                out2.width(17);
                out2  << VEC_ELEM(frc_noise, i);
                out2.width(17);
                out2  << 1./VEC_ELEM(freq, i) << std::endl;
            }
        }
        out2.close();
    }
}


int main(int argc, char **argv)
{
    if (!checkParameter(argc, argv, "-set_of_images"))
    {

        Resolution_parameters prm;
        prm.each_image_produces_an_output = false;
        prm.apply_geo = true;
        SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);

    }
    else
    {

        SelFile     SF, SF1, SF2;
        Image       It, I1, I2, Id;
        FileName    fn_sel;
        float       sam;
        double      dummy;
        bool        apply_geo;
        Matrix1D<double> freq, frc, dpr, frc_noise, ssnr, pixel;

        try
        {
            fn_sel = getParameter(argc, argv, "-set_of_images");
            SF.read(fn_sel);
            sam = textToFloat(getParameter(argc, argv, "-sam"));
            apply_geo = !checkParameter(argc, argv, "-dont_apply_geo");
            SF.split_in_two(SF1, SF2);
            SF.get_statistics(It, Id, dummy, dummy, apply_geo);
            SF1.get_statistics(I1, Id, dummy, dummy, apply_geo);
            SF2.get_statistics(I2, Id, dummy, dummy, apply_geo);
            It().setXmippOrigin();
            I1().setXmippOrigin();
            I2().setXmippOrigin();

            fourier_ring_correlation(I1(), I2(), sam, freq, frc, frc_noise);
            differential_phase_residual(I1(), I2(), sam, freq, dpr);
            //ssnr follows a totally different aproach from FRC and SPR,
            //now we need to pass the whole sel file,
            //The only reason why I send I1 is to keep the
            //feed the template system (ROB March 2005)
            my_ssnr(It(), SF, (double)sam, freq, ssnr, pixel, apply_geo);

            // Write output
            FileName fn_dpr, fn_frc, fn_ssnr;
            fn_dpr = fn_sel + ".dpr";
            fn_frc = fn_sel + ".frc";
            fn_ssnr = fn_sel + ".snr";
            std::ofstream out(fn_dpr.c_str(), std::ios::out);
            std::ofstream out2(fn_frc.c_str(), std::ios::out);
            std::ofstream out3(fn_ssnr.c_str(), std::ios::out);
            out  << "# Resol. [1/Ang]   DPR [deg]" << std::endl;
            out2 << "# Resol. [1/Ang]      FRC      FRC_random_noise" << std::endl;
            out3 << "# Resol. [1/Ang]       SSNR          #Pixels" << std::endl;
            FOR_ALL_ELEMENTS_IN_MATRIX1D(freq)
            {
                out.width(10);
                out  << VEC_ELEM(freq, i);
                out.width(17);
                out << VEC_ELEM(dpr, i)  << std::endl;
                out2.width(10);
                out2  << VEC_ELEM(freq, i);
                out2.width(17);
                out2  << VEC_ELEM(frc, i);
                out2.width(17);
                out2  << VEC_ELEM(frc_noise, i) << std::endl;
                out3.width(10);
                out3  << VEC_ELEM(freq, i);
                out3.width(17);
                out3  << VEC_ELEM(ssnr, i);
                out3.width(17);
                out3  << VEC_ELEM(pixel, i) << std::endl;
            }
            out.close();
            out2.close();
            out3.close();
        }
        catch (Xmipp_error XE)
        {
            std::cout << XE;
            Resolution_parameters prm;
            prm.usage();
        }

    }
}

