/***************************************************************************
 *
 * Authors:    Joaquin Oton                (joton@cnb.csic.es)
 *             J.M.de la Rosa Trevin       (jmdelarosa@cnb.csic.es)
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
 * MERCHANTABILITY or FITNESS FO A PARTICULAR PURPOSE.  See the
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

#include <data/args.h>
#include <data/image.h>
#include <data/metadata.h>
#include <data/progs.h>

typedef enum { PRINT, EXTRACT, ASSIGN, RESET } HeaderOperation;

class ProgHeader: public XmippMetadataProgram
{
protected:
    HeaderOperation operation;
    bool round_shifts;
    MDRow row;

    void defineParams()
    {
        produces_an_output = true;
        //each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addUsageLine("Operate with image files headers. By default in Xmipp images files headers are ignored.");
        addUsageLine("This information is read from the images metadata if exist. With this program headers");
        addUsageLine("can be imported/exported, printed or reset.");
        addParamsLine("[   --print <decompose=0>]    : Print the geometrical transformations in image file headers.");
        addParamsLine("                              : if input is stack and decompose=1 print header of each individual image.");
        addParamsLine("       alias -p;");
        addParamsLine("or --extract     : The output is a selfile with geometrical transformations read from image file headers.");
        addParamsLine("       alias -e;");
        addParamsLine("      requires -o;");
        addParamsLine("or --assign     : Write the geometrical transfromations from selfile to the image file headers.");
        addParamsLine("       alias -a;");
        addParamsLine("or --reset      : Reset the geometrical transformations in image file headers.");
        addParamsLine("       alias -r;");
        addParamsLine("   [--round_shifts]    :Round shifts to integers");
    }

    void readParams()
    {
        if (checkParam("--extract"))
        {
            operation = EXTRACT;
            save_metadata_stack = true;
        }
        else if (checkParam("--assign"))
            operation = ASSIGN;
        else if (checkParam("--reset"))
            operation = RESET;
        else
        {
            operation = PRINT;
            allow_time_bar = false;
            decompose_stacks = getIntParam("--print") == 1;
        }
        XmippMetadataProgram::readParams();
        round_shifts = checkParam("--round_shifts");
    }

    void show()
    {
        String msg;
        switch (operation)
        {
        case PRINT:
            msg = "Printing headers...";
            break;
        case EXTRACT:
            msg = "Extracting image(s) geometrical transformations from header to metadata...";
            break;
        case ASSIGN:
            msg = "Assigning image(s) geometrical transformations from metadata to header...";
            break;
        case RESET:
            msg = "Reseting geometrical transformations from headers...";
            break;
        }
        std::cout << msg << std::endl << "Input: " << fn_in << std::endl;

        if (checkParam("-o"))
            std::cout << "Output: " << fn_out << std::endl;

    }

    void roundShifts()
    {
        double aux = 0.;
        if (row.getValue(MDL_SHIFTX, aux))
        {
            aux = (double)ROUND(aux);
            row.setValue(MDL_SHIFTX, aux);
        }
        if (row.getValue(MDL_SHIFTY, aux))
        {
            aux = (double)ROUND(aux);
            row.setValue(MDL_SHIFTY, aux);
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {

        ImageGeneric img;

        switch (operation)
        {
        case PRINT:
            img.read(fnImg, _HEADER_ALL);
            img.print();
            break;
        case EXTRACT:
            img.read(fnImg, _HEADER_ALL);
            row = img.getGeometry();
            if (round_shifts)
                roundShifts();
            mdIn.setRow(row, objId);
            break;
        case ASSIGN:
            mdIn.getRow(row, objId);
            if (round_shifts)
                roundShifts();
            img.readApplyGeo(fnImg, row, HEADER);
            img.setDataMode(_HEADER_ALL);
            img.write(fnImg, -1, fnImg.isInStack(), WRITE_REPLACE);
            break;
        case RESET:
            img.read(fnImg, _HEADER_ALL);
            img.initGeometry();
            img.write(fnImg, -1, fnImg.isInStack(), WRITE_REPLACE);
            break;
        }
    }

    void finishProcessing()
    {
        if (operation == EXTRACT)
        {
            mdOut = mdIn;
            single_image = false;
        }
        XmippMetadataProgram::finishProcessing();
    }
}
;// end of class ProgHeader

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgHeader program;
    program.read(argc, argv);
    program.tryRun();
}


