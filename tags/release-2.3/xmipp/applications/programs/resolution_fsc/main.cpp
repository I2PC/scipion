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
#include <data/fftw.h>
#include <data/selfile.h>

void writeFiles(const FileName &fnRoot,
    const Matrix1D<double> freq, const Matrix1D<double> frc, 
    const Matrix1D<double> frc_noise, const Matrix1D<double> dpr,
    bool do_dpr,float max_sam)
{
    FileName fn_dpr, fn_frc;
    std::ofstream outDpr;
    fn_dpr = fnRoot + ".dpr";
    if(do_dpr)
    {
        outDpr.open(fn_dpr.c_str(), std::ios::out);
        outDpr  << "# Resol. [1/Ang]   DPR [deg]     Resol. [Ang]" << std::endl;
    }
    fn_frc = fnRoot + ".frc";
    std::cerr << fn_frc;
    std::ofstream outFsc(fn_frc.c_str(), std::ios::out);
    outFsc << "# Resol. [1/Ang]      FRC      FRC_random_noise     Resol. [Ang]" << std::endl;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(freq)
    {
        if (i>0)
        {
            if(max_sam >=0 && ((1./VEC_ELEM(freq, i))<max_sam) )
            {
                if(do_dpr)
                    VEC_ELEM(dpr, i)=0.;
                VEC_ELEM(frc, i)=0.;
                ;
            }
            if(do_dpr)
            {
                outDpr.width(10);
                outDpr  << VEC_ELEM(freq, i);
                outDpr.width(17);
                outDpr << VEC_ELEM(dpr, i);
                outDpr.width(17);
                outDpr << 1./VEC_ELEM(freq, i)  << std::endl;
            }
            outFsc.width(10);
            outFsc  << VEC_ELEM(freq, i);
            outFsc.width(17);
            outFsc  << VEC_ELEM(frc, i);
            outFsc.width(17);
            outFsc  << VEC_ELEM(frc_noise, i);
            outFsc.width(17);
            outFsc  << 1./VEC_ELEM(freq, i) << std::endl;
        }
    }
    if(do_dpr)
        outDpr.close();
    outFsc.close();
}

class Resolution_parameters: public Prog_parameters
{
public:
    FileName    fn_ref, fn_root;
    ImageXmipp  refI;
    VolumeXmipp refV;
    float       sam;
    bool        do_dpr;
    float       max_sam;
public:
    Resolution_parameters()
    {}
    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        fn_ref = getParameter(argc, argv, "-ref");
        sam = textToFloat(getParameter(argc, argv, "-sam"));
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
        std::cout << "Reference file  = " << fn_ref  << std::endl;
        std::cout << "Sampling rate   = " << sam     << std::endl;
        std::cout << "Max Resolution  = " << max_sam << std::endl;
        if(do_dpr)
            std::cout << "-do_dpr is ON  = " << std::endl;
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
        std::cerr << "  [-dont_apply_geo]         : for 2D-images: do not apply transformation stored in the header\n";
        std::cerr << "  [-do_dpr]                 : compute dpr, by default only frc is computed\n";
        std::cerr << "  [-max_sam=-1]             : set fsc to 0 for frequencies above this one (Ang), -1-> all fequencies\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Resolution_parameters *eprm = (Resolution_parameters *) prm;
    Matrix1D<double> freq, frc, dpr, frc_noise;
    frc_dpr(eprm->refI(), img(), eprm->sam, freq, frc, frc_noise, dpr,eprm->do_dpr);
    writeFiles(img.name(), freq, frc, frc_noise, dpr,eprm->do_dpr,eprm->max_sam);
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Resolution_parameters *eprm = (Resolution_parameters *) prm;
    Matrix1D<double> freq, frc, dpr, frc_noise;
    frc_dpr(eprm->refV(), vol(), eprm->sam, freq, frc, frc_noise,dpr,eprm->do_dpr);
    writeFiles(vol.name(), freq, frc, frc_noise, dpr,eprm->do_dpr,eprm->max_sam);
    return true;
}

int main(int argc, char **argv)
{   
    bool do_dpr=checkParameter(argc, argv, "-do_dpr");
    float max_sam = textToFloat(getParameter(argc, argv, "-max_sam","-1"));
    if (!checkParameter(argc, argv, "-set_of_images"))
    {
        Resolution_parameters prm;
        prm.each_image_produces_an_output = false;
        prm.apply_geo = true;
        prm.do_dpr=do_dpr;
        prm.max_sam=max_sam;
        SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
    }
    else
    {
        FileName    fn_sel;
        float       sam;
        bool        apply_geo;
        try
        {
            fn_sel = getParameter(argc, argv, "-set_of_images");
            sam = textToFloat(getParameter(argc, argv, "-sam"));
            apply_geo = !checkParameter(argc, argv, "-dont_apply_geo");
        }
        catch (Xmipp_error XE)
        {
            std::cout << XE;
            Resolution_parameters prm;
            prm.usage();
        }

        try
        {
            SelFile     SF, SF1, SF2;
            Image       I1, I2, Id;
            double      dummy;
            Matrix1D<double> freq, frc, dpr, frc_noise, ssnr, pixel;
            SF.read(fn_sel);
            SF.split_in_two(SF1, SF2);
            SF1.get_statistics(I1, Id, dummy, dummy, apply_geo);
            SF2.get_statistics(I2, Id, dummy, dummy, apply_geo);
            I1().setXmippOrigin();
            I2().setXmippOrigin();
            frc_dpr(I1(), I2(), sam, freq, frc, frc_noise, dpr,do_dpr);
            writeFiles(fn_sel, freq, frc, frc_noise, dpr,do_dpr,max_sam);
        }
        catch (Xmipp_error XE)
        {
            std::cout << XE;
        }
    }
}
