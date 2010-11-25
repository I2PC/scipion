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

#include <data/progs.h>
#include <data/args.h>

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
        addUsageLine("Add random noise to the input images");
        addParamsLine("-gaussian <stddev> <avg=0.>        :Gaussian noise parameters");
        addParamsLine("or -student <df> <stddev> <avg=0.> :t-student noise parameters");
        addParamsLine("or -uniform  <min> <max>           :Uniform noise parameters");
        addParamsLine("  [-limit0 <float> ]               :Crop noise histogram below this value ");
        addParamsLine("  [-limitF <float> ]               :Crop noise histogram above this value ");

    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        do_limit0 = checkParam("-limit0");
        if (do_limit0)
            limit0 =  getDoubleParam("-limit0");
        do_limitF = checkParam("-limitF");
        if (do_limitF)
            limitF =  getDoubleParam("-limitF");

        ///Default value of df in addNoise function
        df = 3.;
        if (checkParam("-gaussian"))
        {
            noise_type = "gaussian";
            param1 = getDoubleParam("-gaussian", 0);
            param2 = getDoubleParam("-gaussian", 1);
        }
        else if (checkParam("-student"))
        {
            noise_type = "student";
            df = getDoubleParam("-student", 0);
            param1 = getDoubleParam("-student", 1);
            param2 = getDoubleParam("-student", 2);
        }
        else if (checkParam("-uniform"))
        {
            noise_type = "uniform";
            param1 = getDoubleParam("-uniform", 0);
            param2 = getDoubleParam("-uniform", 1);
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


    void processImage()
    {
        img.read(fnImg);
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
    program.tryRun();
}

