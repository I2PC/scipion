/***************************************************************************
 *
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#include <data/xmipp_strings.h>
#include <data/xmipp_program.h>
#include <data/image_collection.h>

class ProgStackCreate: public XmippMetadataProgram
{
protected:
    bool write_sel, randomize, replace;
    String case_type;
    FileName stackOut;
    ImageCollection *collection;
    int mode;

    void defineParams()
    {
        produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addParamsLine(" [--write_sel <selfile>] : Also produce a new selfile with the images of the stack.");
        addParamsLine(" [--case <case_type> ] : Also produce a new selfile with the images of the stack.");
        addParamsLine(" where <case_type> stack stackc normal");
        addParamsLine(" [--randomize ] : Randomize input metadata.");
        addParamsLine(" [--replace ] : Replace images on stack");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        stackOut = fn_out;

        randomize = checkParam("--randomize");
        replace = checkParam("--replace");

        if (write_sel = checkParam("--write_sel"))
            fn_out = getParam("--write_sel");

        case_type = getParam("--case");

        mode = replace ? WRITE_REPLACE : WRITE_APPEND;

        if (case_type == "stackc")
        {
            collection = new ImageCollection(mdIn, mode);
            std::cout << "using collection......" << std::endl;
        }


    }

    void show()
    {
        XmippMetadataProgram::show();
    }

    void preProcess()
    {
        if (randomize)
        {
            MetaData mdAux(mdIn);
            mdIn.randomize(mdAux);
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, long int objId)
    {
        static int counter = 0;
        std::cerr << "================ IMAGE " << (counter + 1) << "=====" << std::endl;
        if (case_type == "stackc")
        {
            collection->readImage(img, fnImg);

            collection->writeImage(img, stackOut, counter, true);
        }
        else if(case_type == "stack")
        {
            std::cerr << "-------------> Reading image" <<std::endl;
            Image<double> img; img.read(fnImg);;
            //if (replace)
            {

                FileName fn;
                fn.compose("test", counter+1, "xmp");
                std::cerr << "-------------> Saving temp image" << fn <<std::endl;
                img.write(fn);
            }
            std::cerr << "-------------> Writing image to Stack" <<std::endl;
            img.write(stackOut, counter, true, mode);
        }
        else
        {
            Image<double> img; img.read(fnImg);;
            img.write(fnImg+"new");
        }
        //Append image to the newly created stack
        if (write_sel)
        {
            mdOut.addObject();
            fnImgOut.compose(counter, stackOut);
            mdOut.setValue(MDL_IMAGE, fnImgOut);
            mdOut.setValue(MDL_ENABLED, 1);
        }
        ++counter;
    }

}
;//end of class ProgAddNoise

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    try
    {
        ProgStackCreate program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }
}

