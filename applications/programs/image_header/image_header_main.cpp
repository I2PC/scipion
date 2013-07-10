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

#include <data/xmipp_image.h>
#include <data/metadata.h>
#include <data/xmipp_program.h>
#include <data/xmipp_hdf5.h>


typedef enum { HEADER_PRINT, HEADER_EXTRACT, HEADER_ASSIGN, HEADER_RESET, HEADER_SAMPLINGRATE, HEADER_TREE } HeaderOperation;

class ProgHeader: public XmippMetadataProgram
{
protected:
    HeaderOperation operation;
    bool round_shifts;
    MDRow row;
    ApplyGeoParams params;
    double sampling;

    void defineParams()
    {
        produces_an_output = true;
        get_image_info = false;

        XmippMetadataProgram::defineParams();
        addUsageLine("Operate with image files headers. By default in Xmipp, geometrical transformations");
        addUsageLine("comming in images files headers are ignored. Instead this information is read from");
        addUsageLine("the images metadata if exist. With this program geometrical transformations can be");
        addUsageLine("extracted to a metadata or assigned to header, also allows print or reset image file headers.");
        addParamsLine("[   --print <decompose=0>]    : Print the geometrical transformations in image file headers.");
        addParamsLine("                              : if input is stack and decompose=1 print header of each individual image.");
        addParamsLine("       alias -p;");
        addParamsLine("or --extract     : The output is a selfile with geometrical transformations read from image file headers.");
        addParamsLine("       alias -e;");
        addParamsLine("      requires -o;");
        addParamsLine("or --assign     : Write the geometrical transformations from selfile to the image file headers.");
        addParamsLine("       alias -a;");
        addParamsLine("or --reset      : Reset the geometrical transformations in image file headers.");
        addParamsLine("       alias -r;");
        addParamsLine("or --tree      : Print the tree scheme from file containers as hdf5 files.");
        addParamsLine("       alias -t;");
        addParamsLine("or --sampling_rate <Ts=-1>  : Change the sampling rate (in Angstrom units) in the image file header.");
        addParamsLine("          : If no value is passed then current value in header is print.");
        addParamsLine("       alias -s;");
        addParamsLine("   [--round_shifts]    :Round shifts to integers");
        addExampleLine("Print the header of the images in metadata: ", false);
        addExampleLine("xmipp_image_header -i images.sel");
        addExampleLine("Extract geometrical transformations from image file headers: ", false);
        addExampleLine("xmipp_image_header -i smallStack.stk --extract -o header.doc");
        addExampleLine("Assign the geometrical transformations from the metadata to header: ", false);
        addExampleLine("xmipp_image_header -i header.doc --assign");
        addSeeAlsoLine("transform_geometry");
        addKeywords("header, geometric, transformation, print");
    }

    void readParams()
    {
        if (checkParam("--extract"))
        {
            operation = HEADER_EXTRACT;
            produces_a_metadata = true;
        }
        else if (checkParam("--assign"))
            operation = HEADER_ASSIGN;
        else if (checkParam("--reset"))
            operation = HEADER_RESET;
        else if (checkParam("--sampling_rate"))
        {
            operation = HEADER_SAMPLINGRATE;
            sampling = getDoubleParam("--sampling_rate");
        }
        else if (checkParam("--tree"))
        {
            operation = HEADER_TREE;
            decompose_stacks = false;
        }
        else
        {
            operation = HEADER_PRINT;
            allow_time_bar = false;
            decompose_stacks = getIntParam("--print") == 1;
        }
        XmippMetadataProgram::readParams();
        round_shifts = checkParam("--round_shifts");
        if (operation != HEADER_EXTRACT && checkParam("-o"))
            REPORT_ERROR(ERR_PARAM_INCORRECT, "Argument -o is not valid for this operation");
        params.datamode = HEADER;
    }

    void show()
    {
        if (verbose == 0)
            return;

        String msg;
        switch (operation)
        {
        case HEADER_PRINT:
            msg = "Printing headers...";
            break;
        case HEADER_EXTRACT:
            msg = "Extracting image(s) geometrical transformations from header to metadata...";
            break;
        case HEADER_ASSIGN:
            msg = "Assigning image(s) geometrical transformations from metadata to header...";
            break;
        case HEADER_RESET:
            msg = "Reseting geometrical transformations from headers...";
            break;
        case HEADER_SAMPLINGRATE:
            if (sampling > 0)
                msg = "Setting sampling rate into headers...";
            else
                msg = "Showing sampling rate from headers...";
            break;
        case HEADER_TREE:
            msg = "Printing tree structure...";
            break;
        }
        std::cout << msg << std::endl << "Input: " << fn_in << std::endl;

        if (checkParam("-o"))
            std::cout << "Output: " << fn_out << std::endl;

    }

    void roundShifts(MDRow &row)
    {
        double aux = 0.;
        if (row.getValue(MDL_SHIFT_X, aux))
        {
            aux = (double)ROUND(aux);
            row.setValue(MDL_SHIFT_X, aux);
        }
        if (row.getValue(MDL_SHIFT_Y, aux))
        {
            aux = (double)ROUND(aux);
            row.setValue(MDL_SHIFT_Y, aux);
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {

        Image<char> img;

        switch (operation)
        {
        case HEADER_PRINT:
            img.read(fnImg, _HEADER_ALL);
            std::cout << img <<std::endl;
            break;
        case HEADER_EXTRACT:
            img.read(fnImg, _HEADER_ALL);
            rowOut = img.getGeometry();
            if (round_shifts)
                roundShifts(rowOut);
            rowOut.setValue(MDL_IMAGE, fnImgOut);
            break;
        case HEADER_ASSIGN:
            rowOut = rowIn;
            if (round_shifts)
                roundShifts(rowOut);
            img.readApplyGeo(fnImg, rowOut, params);
            img.setDataMode(_HEADER_ALL);
            img.write(fnImg, ALL_IMAGES, fnImg.isInStack(), WRITE_REPLACE);
            break;
        case HEADER_RESET:
            img.read(fnImg, _HEADER_ALL);
            img.initGeometry();
            img.write(fnImg, ALL_IMAGES, fnImg.isInStack(), WRITE_REPLACE);
            break;
        case HEADER_SAMPLINGRATE:
            img.read(fnImg, _HEADER_ALL);
            if (sampling < 0)
            {
                img.MDMainHeader.getValue(MDL_SAMPLINGRATE_X, sampling);
                std::cout << sampling << std::endl;
            }
            else
            {
                img.MDMainHeader.setValue(MDL_SAMPLINGRATE_X, sampling);
                img.MDMainHeader.setValue(MDL_SAMPLINGRATE_Y, sampling);
                img.MDMainHeader.setValue(MDL_SAMPLINGRATE_Z, sampling);
                img.write(fnImg, ALL_IMAGES, fnImg.isInStack(), WRITE_REPLACE);
                std::cout << "New sampling rate (Angstrom) = " << sampling << std::endl;
            }
            break;
        case HEADER_TREE:
            {
                XmippH5File H5File;
                FileName filename = fnImg.removeAllPrefixes().removeFileFormat();

                if (H5File.isHdf5(filename.c_str()))
                {
                    H5File.openFile(fnImg.removeAllPrefixes().removeFileFormat(), H5F_ACC_RDONLY);
                    H5File.showTree();
                }
                else
                    REPORT_ERROR(ERR_IMG_UNKNOWN, "Unknown file format to display its data structure.");
            }
            break;
        }
    }

    void finishProcessing()
    {
        if (operation ==HEADER_EXTRACT)
            single_image = false;
        XmippMetadataProgram::finishProcessing();
    }
}
;// end of class ProgHeader

RUN_XMIPP_PROGRAM(ProgHeader)


