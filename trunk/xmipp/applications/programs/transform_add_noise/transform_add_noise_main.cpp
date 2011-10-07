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

#include <data/xmipp_program.h>

class ProgAddNoise: public XmippMetadataProgram
{
protected:
    double param1, param2;
    double df, limit0, limitF;
    bool   do_limit0, do_limitF;
    std::string noise_type;

    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        //Usage
        addUsageLine("Add random noise to the input images.");
        addUsageLine("Noise can be generated using uniform, gaussian or t-student distributions.");
        //Parameters
        addParamsLine("--type <rand_mode>                : Type of noise to add");
        addParamsLine("       where <rand_mode>");
        addParamsLine("              gaussian <stddev> <avg=0.>     :Gaussian distribution parameters");
        addParamsLine("              student <df> <stddev> <avg=0.> :t-student distribution parameters");
        addParamsLine("              uniform  <min> <max>           :Uniform distribution parameters");
        addParamsLine("  [--limit0 <float> ]               :Crop noise histogram below this value ");
        addParamsLine("  [--limitF <float> ]               :Crop noise histogram above this value ");
        //Examples
        addExampleLine("Add noise to a single image, writing in different image:", false);
        addExampleLine("xmipp_transform_add_noise -i cleanImage.spi --type gaussian 10 5 -o noisyGaussian.spi");
        addExampleLine("Add uniform noise to a volume, overriding input volume:", false);
        addExampleLine("xmipp_transform_add_noise -i g0ta.vol -uniform -0.1 0.1");

    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        do_limit0 = checkParam("--limit0");
        if (do_limit0)
            limit0 =  getDoubleParam("-limit0");
        do_limitF = checkParam("--limitF");
        if (do_limitF)
            limitF =  getDoubleParam("--limitF");

        ///Default value of df in addNoise function
        df = 3.;
        noise_type = getParam("--type");

        if (noise_type == "gaussian")
        {
            param2 = getDoubleParam("--type", 1);
            param1 = getDoubleParam("--type", 2);
        }
        else if (noise_type == "student")
        {
            df = getDoubleParam("--type", 1);
            param1 = getDoubleParam("--type", 2);
            param2 = getDoubleParam("--type", 3);
        }
        else if (noise_type == "uniform")
        {
            param1 = getDoubleParam("--type", 1);
            param2 = getDoubleParam("--type", 2);
        }
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "Unknown noise type");
    }

    void show()
    {
        XmippMetadataProgram::show();
        if (noise_type == "gaussian")
            std::cout << "Noise avg=" << param1 << std::endl
            << "Noise stddev=" << param2 << std::endl;
        else if (noise_type == "student")
            std::cout << "Degrees of freedom= "<< df << std::endl
            << "Noise avg=" << param1 << std::endl
            << "Noise stddev=" << param2 << std::endl;
        else if (noise_type == "uniform")
            std::cout << "Noise min=" << param1 << std::endl
            << "Noise max=" << param2 << std::endl;
        if (do_limit0)
            std::cout << "Crop noise histogram below=" << limit0 << std::endl;
        if (do_limitF)
            std::cout << "Crop noise histogram above=" << limitF << std::endl;
    }


    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
        Image<double> img;
        img.readApplyGeo(fnImg, rowIn);
        img().addNoise(param1, param2, noise_type, df);

        if (do_limit0)
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
        {
            dAkij(img(),k,i,j) = XMIPP_MAX(dAkij(img(),k,i,j), limit0);
        }

        if (do_limitF)
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
        {
            dAkij(img(),k,i,j) = XMIPP_MIN(dAkij(img(),k,i,j), limitF);
        }

        img.write(fnImgOut);

    }

}
;//end of class ProgAddNoise

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgAddNoise program;
    program.read(argc, argv);
    //init random seed
    randomize_random_generator();
    program.tryRun();
    return 0;
}

