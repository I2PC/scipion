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


class ProgResolutionFsc : public XmippProgram
{
public:

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

    FileName    fn_ref, fn_root,fn_img;
    float       sam;
    float       max_sam;
    bool        do_dpr, do_set_of_images;

    FileName    fn_sel;
    bool        apply_geo;

    void readParams()
    {
        sam = getDoubleParam("--sam");

        apply_geo = !checkParam("--dont_apply_geo");


        max_sam = getDoubleParam("--max_sam");
        do_dpr = checkParam("--do_dpr");
        do_set_of_images = checkParam("--set_of_images");
        if(do_set_of_images)
        {
            fn_sel = getParam("--set_of_images");
            if (checkParam("-i") || checkParam("--ref"))
                REPORT_ERROR(ERR_ARG_INCORRECT, "--set_of_images should not be provided with -i or --ref");
        }
        else
        {
            fn_ref = getParam("--ref");
            fn_img = getParam("-i");
        }
    }

    void defineParams()
    {

        apply_geo = true;

        addUsageLine("Calculate the resolution of one or more volumes or images with respect to a single reference.");
        addUsageLine("Example of use: Resolution of subset2.vol volume with respect to subset1.vol reference volume using 5.6 pixel size (in Angstrom)");
        addUsageLine("   xmipp_resolution_fsc --ref subset1.vol  -i subset2.vol --sam 5.6 ");
        addUsageLine("Example of use: Resolution of a set of images using 5.6 pixel size (in Angstrom)");
        addUsageLine("   xmipp_resolution_fsc --set_of_images selfile.sel --sam 5.6 ");

        addParamsLine("   -i <input_file>           : either an image/volume or a selection file");
        addParamsLine("   requires --ref;");
        addParamsLine("or --set_of_images <selfile> : SelFile containing a set of 2D-images");
        addParamsLine("   [--ref <input_file>]        : Filename for reference image/volume");
        addParamsLine("   --sam <sampling_rate>     : i.e. pixel size (in Angstrom)");
        addParamsLine("   [--dont_apply_geo]        : for 2D-images: do not apply transformation stored in the header");
        addParamsLine("   [--do_dpr]                : compute dpr, by default only frc is computed");
        addParamsLine("   [--max_sam <max_sr=-1>]   : set fsc to 0 for frequencies above this one (Ang), -1-> all fequencies");
    }


    bool process_img()
    {
        Image<double>  refI,img;
        refI.read(fn_ref, true, -1, apply_geo);
        refI().setXmippOrigin();
        img.read(fn_img);
        img().setXmippOrigin();
        MultidimArray<double> freq, frc, dpr, frc_noise, error_l2;
        frc_dpr(refI(), img(), sam, freq, frc, frc_noise, dpr, error_l2);
        writeFiles(img.name(), freq, frc, frc_noise, dpr, error_l2, max_sam, do_dpr);
        return true;
    }

    bool process_sel()
    {
        MetaData     SFaux,SF, SF1, SF2;
        std::vector<MetaData> vMD;
        vMD.push_back(SF1);
        vMD.push_back(SF2);

        Image<double>       I1, I2, Id;
        double      dummy;
        MultidimArray<double> freq, frc, dpr, frc_noise, ssnr, error_l2, pixel;

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

    void run()
    {

        if (!do_set_of_images)
        {
            process_img();
        }
        else
        {
            process_sel();
        }
    }

};

int main(int argc, char **argv)
{
    try
    {
        ProgResolutionFsc program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
        return 1;
    }
    return 0;

}
