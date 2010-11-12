/***************************************************************************
 *
 * Authors:     Sjors Scheres (scheres@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/progs.h>
#include <data/args.h>
#include <data/fftw.h>
#include <data/metadata_extension.h>


//class ProgResolutionFsc : public XmippProgram
//{
//public:

    void writeFiles(const FileName &fnRoot,
                    const MultidimArray<double> &freq,
                    const MultidimArray<double> &frc,
                    const MultidimArray<double> &frc_noise,
                    const MultidimArray<double> &dpr,
                    const MultidimArray<double> &error_l2,
                    float max_sam, bool do_dpr)
    {
        MetaData MD;

        FileName  fn_frc;
        fn_frc = fnRoot + ".frc";
        FOR_ALL_ELEMENTS_IN_ARRAY1D(freq)
        {
            if (i>0)
            {
                MD.addObject();
                if(max_sam >=0 && ((1./dAi(freq, i))<max_sam) )
                {
                    if(do_dpr)
                        dAi(dpr, i)=0.;
                    dAi(frc, i)=0.;
                }
                MD.setValue(MDL_RESOLUTION_FREQ,dAi(freq, i));
                MD.setValue(MDL_RESOLUTION_FRC,dAi(frc, i));
                if(do_dpr)
                    MD.setValue(MDL_RESOLUTION_DPR,dAi(dpr, i));
                MD.setValue(MDL_RESOLUTION_ERRORL2,dAi(error_l2, i));
                MD.setValue(MDL_RESOLUTION_FRCRANDOMNOISE,dAi(frc_noise, i));
                MD.setValue(MDL_RESOLUTION_FREQREAL,1./dAi(freq, i));
            }
        }
        MD.write(fn_frc);
    }

class Resolution_parameters: public Prog_parameters
    {
    public:
        FileName    fn_ref, fn_root;
        Image<double>  refI;
        float       sam;
        float       max_sam;
        bool        do_dpr;
    public:
        void read(int argc, char **argv)
        {
            Prog_parameters::read(argc, argv);
            fn_ref = getParameter(argc, argv, "-ref");
            sam = textToFloat(getParameter(argc, argv, "-sam"));
            refI.read(fn_ref, true, -1, apply_geo);
            refI().setXmippOrigin();
        }
        void show()
        {
            Prog_parameters::show();
            std::cout << "Reference file  = " << fn_ref  << std::endl;
            std::cout << "Sampling rate   = " << sam     << std::endl;
            std::cout << "Max Resolution  = " << max_sam << std::endl;
            if(do_dpr)
                std::cout << "do_dpr is ON  = " << std::endl;
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

    bool process_img(Image<double> &img, const Prog_parameters *prm)
    {
        Resolution_parameters *eprm = (Resolution_parameters *) prm;
        MultidimArray<double> freq, frc, dpr, frc_noise, error_l2;
        frc_dpr(eprm->refI(), img(), eprm->sam, freq, frc, frc_noise, dpr, error_l2);
        writeFiles(img.name(), freq, frc, frc_noise, dpr, error_l2, eprm->max_sam,eprm->do_dpr);
        return true;
    }

//};

int main(int argc, char **argv)
{
    float max_sam = textToFloat(getParameter(argc, argv, "-max_sam","-1"));
    bool do_dpr=checkParameter(argc, argv, "-do_dpr");
    if (!checkParameter(argc, argv, "-set_of_images"))
    {
        Resolution_parameters prm;
        prm.each_image_produces_an_output = false;
        prm.apply_geo = true;
        prm.do_dpr=do_dpr;
        prm.max_sam=max_sam;
        SF_main(argc, argv, &prm, (void*)&process_img);
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
        catch (XmippError XE)
        {
            std::cout << XE;
            Resolution_parameters prm;
            prm.usage();
        }

        try
        {
            MetaData     SFaux,SF, SF1, SF2;
            std::vector<MetaData> vMD;
            vMD.push_back(SF1);
            vMD.push_back(SF2);
            Image<double>       I1, I2, Id;
            double      dummy;
            MultidimArray<double> freq, frc, dpr, frc_noise, ssnr, error_l2,
            pixel;

            SFaux.read(fn_sel);
            SF.randomize(SFaux);
            SF.split(2,vMD,MDL_IMAGE);
            getStatistics(SF1,I1,Id,dummy,dummy,true);
            getStatistics(SF2,I2,Id,dummy,dummy,true);
            I1().setXmippOrigin();
            I2().setXmippOrigin();
            frc_dpr(I1(), I2(), sam, freq, frc, frc_noise, dpr,error_l2,do_dpr);
            writeFiles(fn_sel, freq, frc, frc_noise, dpr,error_l2,max_sam,do_dpr);
        }
        catch (XmippError XE)
        {
            std::cout << XE;
        }
    }
}
