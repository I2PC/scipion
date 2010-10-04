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

#include <data/progs.h>
#include <data/args.h>
#include <data/image_collection.h>

class ProgStackCreate: public XmippMetadataProgram
{
protected:
    bool write_sel;
    bool use_collection;
    FileName stackOut;
    ImageCollection *collection;

    void defineParams()
    {
        produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addParamsLine(" [--write_sel <selfile>] : Also produce a new selfile with the images of the stack.");
        addParamsLine(" [--use_collection ] : Also produce a new selfile with the images of the stack.");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        stackOut = fn_out;
        if (write_sel = checkParam("--write_sel"))
            fn_out = getParam("--write_sel");
        if (use_collection = checkParam("--use_collection"))
        {
            collection = new ImageCollection(mdIn);
            std::cout << "using collection......" << std::endl;
        }

    }

    void show()
    {
        XmippMetadataProgram::show();
    }


    void processImage()
    {
        static int counter = 0;
        static int mode = WRITE_APPEND;//(counter == 1) ? WRITE_OVERWRITE : WRITE_APPEND;
        if (use_collection)
        {
            collection->readImage(img, fnImg);
            collection->writeImage(img, stackOut, counter, true, mode);
        }
        else
        {
            img.read(fnImg);
            img.write(stackOut, counter, true, mode);
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

