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

#include <data/xmipp_program.h>
#include <data/xmipp_fftw.h>
#include <data/metadata_extension.h>

class ProgResolutionFsc : public XmippProgram
{
public:

    FileName    fn_ref, fn_root,fn_img, fn_out;
    float       sam;
    float       max_sam;
    bool        do_dpr, do_set_of_images, do_o;

    FileName    fn_sel;
    bool        apply_geo;

    void defineParams()
    {
        apply_geo = true;

        addUsageLine("Calculate the resolution of one or more volumes or images with respect to a single reference.");
        addUsageLine("+ ");
        addUsageLine("+Three methods are employed:");
        addUsageLine("+ ");
        addUsageLine("+* Differential Phase Residual (DPR)",true);
        addUsageLine("+In this method the the resolution is defined as the spatial frequency at which the average phase");
        addUsageLine("+discrepancy between the two transforms exceeds 45 degrees.");
        addUsageLine("+* Fourier Ring Correlation (FRC)", true);
        addUsageLine("+In this method the resolution is defined as the spatial frequency at which annular samplings of");
        addUsageLine("+the two Fourier transforms register negligible cross-correlation.");
        addUsageLine("+* Spectral Signal-to-Noise Ratio (SSNR)", true);
        addUsageLine("+This method is based on the measurement of the signal-to-noise ratio as a function of spatial");
        addUsageLine("+frequency. The SSNR is determined by comparing the Fourier transform of individual images with ");
        addUsageLine("+that of the global average image. (this option is only available for 2D images not for volumes,");
        addUsageLine("+ see the program ssnr if you are interested in computing the signal to noise ratio in 3D).");
        addUsageLine("+ ");
        addUsageLine("+To calculate \"the resolution of a reconstruction\", use --set_of_images, where the program divides the");
        addUsageLine("+corresponding set of projection images into two subsets and reconstructs two volumes from these subsets.");
        addUsageLine("+This program may then be used to calculate the DPR and FRC between these two volumes. The resulting plots");
        addUsageLine("+are commonly used to assess the high-resolution limit of the initial reconstruction (see Frank's book for more details).");
        addUsageLine(" ");
        addUsageLine("The program writes out filename.frc files, for each input volume or image, or selfilename.frc, the");
        addUsageLine("set_of_images mode. These ACSII files contain the DPR, FRC and SSNR as a function of resolution (in 1/Angstrom).");
        addUsageLine(" The .frc files also contain a column for the FRC expected for pure noise.");
        addSeeAlsoLine("resolution_ssnr");

        addParamsLine("   -i <input_file>           : either an image/volume or a selection file");
        addParamsLine("   requires --ref;");
        addParamsLine("or --set_of_images <selfile> : selfile containing a set of 2D-images");
        addParamsLine("   [--oroot <root_file=\"\">] : Root of the output metadata. If not set, input file rootname is taken.");
        addParamsLine("   [-o <output_file=\"\">]   : Output file name.");

        addParamsLine("   [--ref <input_file>]      : filename for reference image/volume");
        addParamsLine("   [--sampling_rate <Ts=1>]  : Pixel size (Angstrom)");
        addParamsLine("  alias -s;");
        addParamsLine("   [--dont_apply_geo]        : for 2D-images: do not apply transformation stored in the header");
        addParamsLine("   [--do_dpr]                : compute dpr, by default only frc is computed");
        addParamsLine("   [--max_sam <max_sr=-1>]   : set fsc to 0 for frequencies above this one (Angstrom), -1 -> all fequencies");

        addExampleLine("Resolution of subset2.vol volume with respect to subset1.vol reference volume using 5.6 pixel size (in Angstrom):", false);
        addExampleLine("xmipp_resolution_fsc --ref subset1.vol  -i subset2.vol --sam 5.6 ");
        addExampleLine("Resolution of a set of images using 5.6 pixel size (in Angstrom):", false);
        addExampleLine("xmipp_resolution_fsc --set_of_images selfile.sel --sam 5.6");
    }

    void readParams()
    {
        sam = getDoubleParam("--sampling_rate");

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

        do_o = checkParam("-o");
        if (do_o)
            fn_out = getParam("-o");
        else
            fn_root = getParam("--oroot");
    }

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
        if (do_o)
            fn_frc = fn_out;
        else
            fn_frc = fnRoot + ".frc";
        size_t id;
        FOR_ALL_ELEMENTS_IN_ARRAY1D(freq)
        {
            if (i>0)
            {
                id=MD.addObject();
                if(max_sam >=0 && ((1./dAi(freq, i))<max_sam) )
                {
                    if(do_dpr)
                        dAi(dpr, i)=0.;
                    dAi(frc, i)=0.;
                }
                MD.setValue(MDL_RESOLUTION_FREQ,dAi(freq, i),id);
                MD.setValue(MDL_RESOLUTION_FRC,dAi(frc, i),id);
                if(do_dpr)
                    MD.setValue(MDL_RESOLUTION_DPR,dAi(dpr, i),id);
                MD.setValue(MDL_RESOLUTION_ERRORL2,dAi(error_l2, i),id);
                MD.setValue(MDL_RESOLUTION_FRCRANDOMNOISE,dAi(frc_noise, i),id);
                MD.setValue(MDL_RESOLUTION_FREQREAL,1./dAi(freq, i),id);
            }
        }
        MD.write(fn_frc);
    }

    bool process_img()
    {
        Image<double>  refI,img;
        //if (apply_geo)
        // refI.readApplyGeo(fn_ref);
        //else
        refI.read(fn_ref);
        refI().setXmippOrigin();
        //if (apply_geo)
        // img.readApplyGeo(fn_img);
        //else
        img.read(fn_img);
        img().setXmippOrigin();

        MultidimArray<double> freq, frc, dpr, frc_noise, error_l2;
        frc_dpr(refI(), img(), sam, freq, frc, frc_noise, dpr, error_l2, do_dpr);

        writeFiles((fn_root.empty())?img.name():fn_root, freq, frc, frc_noise, dpr, error_l2, max_sam, do_dpr);
        return true;
    }

    bool process_sel()
    {
        MetaData MD(fn_sel), MDout;
        MultidimArray<double> freq, frc, dpr, frc_noise, ssnr, error_l2;
        getFourierStatistics(MD, sam, MDout, do_dpr, max_sam);
        FileName fnRoot=(fn_root.empty())?fn_sel:fn_root;
        MDout.write(fnRoot+".frc");
        return true;
    }

    void run()
    {
        if (!do_set_of_images)
            process_img();
        else
            process_sel();
    }
};

int main(int argc, char **argv)
{
    ProgResolutionFsc program;
    program.read(argc, argv);
    program.tryRun();
}
