/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es  (2007)
 *             Joaquin Oton            joton@cnb.csic.es (2010)
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

#include "metadata_extension.h"
#include "xmipp_image_generic.h"
#include "xmipp_image_convert.h"

ProgConvImg::ProgConvImg()
{
    init();
}

void ProgConvImg::init()
{
    each_image_produces_an_output = true;
    save_metadata_stack = false;
    keep_input_columns = true;
    delete_output_stack = false;

    castMode = CW_CONVERT;
    appendToStack = false;
    // output extension
    oext = "";
    // output type
    type = "auto";
    // Set default write mode
    writeMode = WRITE_OVERWRITE;
    depth = "";
    swap = false;
}


void  ProgConvImg::defineParams()
{
    init();
    CommentList &comments = defaultComments["-i"];
    comments.addComment("++ Supported read formats are:");
    comments.addComment("++ dm3 : Digital Micrograph 3");
    comments.addComment("++ em  : Electron Microscopy");
    comments.addComment("++ jpg : JEPG");
    comments.addComment("++ img : Imagic");
    comments.addComment("++ inf,raw : RAW file with header INF file");
    comments.addComment("++ mrc, map : CCP4");
    comments.addComment("++ pif  : Portable Image Format");
    comments.addComment("++ ser : Tecnai Imaging and Analysis");
    comments.addComment("++ spe : Princeton Instruments CCD camera");
    comments.addComment("++ spi, xmp : Spider");
    comments.addComment("++ tif : TIFF");
    comments.addComment("++ raw#xDim,yDim,[zDim],offset,datatype,[r] : RAW image file without header file");
    comments.addComment("++ where datatype can be: uint8,int8,uint16,int16,uint32,int32,long,float,double,");
    comments.addComment("++                        cint16,cint32,cfloat,cdouble,bool");
    XmippMetadataProgram::defineParams();

    addUsageLine("Convert among stacks, volumes and images, and change the file format.");
    addUsageLine("+Conversion to a lower bit_depth automatically adjusts the gray level range. If it is between same bit ");
    addUsageLine("+depths and different sign, then only a histogram shift is done. If parameter --depth is not passed, then ");
    addUsageLine("+bit_depth is automatically chosen equal to or higher than input bit_depth. For stack output format, ");
    addUsageLine("+a selection file with the images in the stack is optionally created, replicating the labels of the input sel file.");
    addUsageLine("+If output file extension is not set when --oroot is used (neither setting --oext nor :ext), then input format is chosen.");
    addKeywords("conversion, convert, image, stack, volume, format, extension ");
    //Parameters
    addParamsLine("  [--oext <extension=\"\">] :  Output file format extension.");
    addParamsLine("    where <extension>");
    addParamsLine("         img : Imagic (Data types: uint8, int16, float* and cfloat).");
    addParamsLine("         inf : RAW file with header INF file (Data types: (u)int8, (u)int16 and float*).");
    addParamsLine("         raw : RAW file with header INF file (Data types: (u)int8, (u)int16 and float*).");
    addParamsLine("         mrc : CCP4 (Data types: uint8, (u)int16, float* and cfloat).");
    addParamsLine("         spi : Spider (Data types: float* and cfloat).");
    addParamsLine("         xmp : Spider (Data types: float* and cfloat).");
    addParamsLine("         tif : TIFF (Data types: uint8*, uint16, uint32 and float).");
    addParamsLine("         jpg : JPEG (Data types: uint8*).");
    addParamsLine("         custom <ext> : Custom extension name, the real format will be Spider.");
    addParamsLine("  [--type <output_type=auto>] : Force output file type.");
    addParamsLine("          where <output_type>");
    addParamsLine("          auto: Autodetect output type according to output extension and whether --oroot is passed or not.");
    addParamsLine("          img : Image");
    addParamsLine("          vol : Volume");
    addParamsLine("          stk : Stack ");
    addParamsLine("  alias -t;");
    addParamsLine("== Bit options == ");
    addParamsLine("  [--depth+ <bit_depth=default>] : Image bit depth.");
    addParamsLine("          where <bit_depth>");
    addParamsLine("                 default: Default selected value (*)");
    addParamsLine("                 uint8 : Equivalent to uchar");
    addParamsLine("                 int8  : Equivalent to char");
    addParamsLine("                 uint16: Equivalent to ushort");
    addParamsLine("                 int16 : Equivalent to short");
    addParamsLine("                 uint32");
    addParamsLine("                 int32");
    addParamsLine("                 long");
    addParamsLine("                 float");
    addParamsLine("                 double");
    addParamsLine("                 cint16 : Complex int16");
    addParamsLine("                 cint32 : Complex int32");
    addParamsLine("                 cfloat : Complex float");
    addParamsLine("                 cdouble: Complex double");
    addParamsLine("                 bool");
    addParamsLine("  alias -d;");
    addParamsLine("  [--swap]        : Swap the endianess of the image file");
    addParamsLine("  [--range_adjust] : Adjust the histogram to fill the gray level range");
    addParamsLine("  alias -r;");
    addParamsLine("or --dont_convert : Do not apply any conversion to gray levels when writing");
    addParamsLine("                  : in a lower bit depth or changing the sign");
    addParamsLine("== Stack options == ");
    addParamsLine("  [--append]           : Append the input to the output stack instead of overwriting it");
    addParamsLine("  alias -a;");

    //Examples
    addExampleLine("Put a selection file into a stack:",false);
    addExampleLine("xmipp_image_convert -i list.sel -o images.stk");
    addExampleLine("Convert a Spider volume to a MRC stack:",false);
    addExampleLine("xmipp_image_convert -i spider.vol -o stack.mrcs -t stk");
    addExampleLine("Create a stack of volumes with a Spider volume :",false);
    addExampleLine("xmipp_image_convert -i spider.vol -o vol_stack.stk -t vol");
    addExampleLine("Append a volume to a volume stack:",false);
    addExampleLine("xmipp_image_convert -i spider.vol -o vol_stack.stk -a");
    addExampleLine("Substitute a volume in a volume stack:",false);
    addExampleLine("xmipp_image_convert -i spider.vol -o 3@vol_stack.stk");
    addExampleLine("Save images in a stack as independent TIFF files in image directory with \"newimage\" basename in 8bit format:",false);
    addExampleLine("xmipp_image_convert -i stackFile.stk -o tiffImages.sel --oroot images/newimage:tif -d uint8");
    addExampleLine("Convert a selection file of 16bit TIFF images to 8bit and overwrites files and sel file:",false);
    addExampleLine("xmipp_image_convert -i tiff16.sel -d uint8");
    addExampleLine("Append a single image to a stack:",false);
    addExampleLine("xmipp_image_convert -i img.spi -o stackFile.stk --append");
    addExampleLine("Append a selection file to a stack:",false);
    addExampleLine("xmipp_image_convert -i selFile.sel -o stackFile.stk --append");
    addExampleLine("Replace a single image into a stack:",false);
    addExampleLine("xmipp_image_convert -i img.spi -o 3@stackFile.stk");
    addExampleLine("Convert a MRC stack to a MRC volume:",false);
    addExampleLine("xmipp_image_convert -i stack.mrc -o volume.mrc -t vol");
}

void ProgConvImg::readParams()
{
    XmippMetadataProgram::readParams();

    //fn_out = (checkParam("-o"))? getParam("-o") : "";
    castMode = (checkParam("--range_adjust"))? CW_ADJUST: \
               (checkParam("--dont_convert"))? CW_CAST: CW_CONVERT;

    appendToStack = checkParam("--append");

    // output extension
    oext = checkParam("--oext") ? getParam("--oext") : "";
    if ( oext == "custom" )
        oext = getParam("--oext",1);

    // Check output type
    type = getParam("--type");

    if (checkParam("--depth"))
    {
        String depthTemp = (String)getParam("--depth");
        if (depthTemp != "default")
            depth = "%" + depthTemp;
    }
    swap = checkParam("--swap");

}

void ProgConvImg::preProcess()
{
    if (type == "auto")
    {
        if (oroot.empty())
        {
            /* It is stack if extension is stack compatible, or if --oroot is not passed
             * and there are more than one image. Same for volumes.
             * Firstly, we must check extensions, then stack and volume sizes.*/
            if (( mdInSize > 1 || zdimOut > 1 ) && fn_out.hasStackExtension())
                type = "stk";
            else if (( mdInSize > 1 || zdimOut > 1 ) && fn_out.hasVolumeExtension()) // if it is volume
                type = "vol";
            else if (mdInSize > 1 || appendToStack) // If --append we suppose output is stack
                type = "stk";
            else if (zdimOut > 1) // if it is volume
                type = "vol";
            else
                type = "img";
        }
        else
            type = "img";
    }

    // Set write mode
    if (single_image && fn_out.isInStack()) // Replace a single image in a stack
    {
        type = "img";
        writeMode = WRITE_REPLACE;
    }
    else if (type == "stk" && appendToStack)
        writeMode = WRITE_APPEND;
    else
    {
        writeMode = WRITE_OVERWRITE;
        delete_output_stack = true;
    }

    if (delete_output_stack)
    {
        FileName fn_stack_plain = fn_out.removeFileFormat();
        fn_stack_plain.deleteFile();
        delete_output_stack = false;
    }

    create_empty_stackfile = (create_empty_stackfile && !appendToStack);

    convMode = MD2MD;

    if (!single_image && type == "vol")
    {
        convMode = MD2VOL;

        if (zdimOut != 1)
            REPORT_ERROR(ERR_MULTIDIM_DIM,
                         "Only 2D images can be converted into volumes");
        imOut = new ImageGeneric(datatypeOut);
        imOut->mapFile2Write(xdimOut, ydimOut, mdInSize, fn_out);
        k = 0;
    }
    else if (single_image)
    {
        // If --append is set, or fn_out is in a stack, then it is supposed not to convert VOL2MD
        if ( zdimOut > 1 && !(type == "vol" || fn_out.isInStack() || appendToStack))
        {
            convMode = VOL2MD;
            single_image = false;

            MetaData * md = getInputMd();
            md->clear();

            FileName fnTemp;

            /* Fill mdIn to allow XmippMetaDataProgram create the fnImgOut,
            but not used to read input images. Input volume is read here. */
            for (k = 1; k <= zdimOut; k++)
            {
                fnTemp.compose(k, fn_in);
                size_t id = md->addObject();
                md->setValue(MDL_IMAGE, fnTemp, id);
                md->setValue(MDL_ENABLED, 1, id);
            }
            imIn.read(fn_in, DATA, ALL_IMAGES, true);
            imOut = new ImageGeneric(imIn.getDatatype());
            k = 0; // Reset to zero to select the slices when working with volumes
            createEmptyFile(fn_out+depth, xdimOut, ydimOut, 1, zdimOut, true, WRITE_OVERWRITE, swap);
        }
    }
    else if (create_empty_stackfile)
        createEmptyFile(fn_out+depth, xdimOut, ydimOut, zdimOut, mdInSize, true, WRITE_OVERWRITE, swap);

    create_empty_stackfile = false;
}//function preprocess

void ProgConvImg::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    switch(convMode)
    {
    case MD2MD:
        {
            FileName tempName, _fnImgOut;
            if (fnImg == fnImgOut)
            {
                tempName.initUniqueName("tempConvert_XXXXXX");
                _fnImgOut = tempName + ":" + fnImgOut.getExtension();
            }
            else
                _fnImgOut= fnImgOut;

            imIn.read(fnImg, DATA, ALL_IMAGES, true);
            imIn.write(_fnImgOut+depth, ALL_IMAGES, type == "stk", writeMode, castMode, swap);

            if ((fnImg == fnImgOut) && (rename(tempName.c_str(),fnImgOut.c_str())!=0))
                REPORT_ERROR(ERR_IO, formatString("ProgConvImg:: Error renaming the file from %s to %s.",tempName.c_str(), fnImgOut.c_str()));

            break;
        }
    case MD2VOL:
        {
            imIn.read(fnImg,DATA, ALL_IMAGES,true);
            imOut->data->setSlice(k++,imIn.data);
            break;
        }
    case VOL2MD:
        {
            imIn.data->getSlice(k++,imOut->data);
            imOut->write(fnImgOut+depth, ALL_IMAGES, type == "stk", writeMode, castMode, swap);
            break;
        }
    }
}//function processImage

void ProgConvImg::finishProcessing()
{
    switch(convMode)
    {
    case MD2VOL:
        if (swap)
            imOut->image->swapOnWrite();
        imOut->write();
        single_image = true;
        progress_bar(time_bar_size);
        break;
    default:
    	break;
    }

    XmippMetadataProgram::finishProcessing();
}

void ProgConvImg::show()
{
    XmippMetadataProgram::show();
    if (each_image_produces_an_output)
    {
        if (!oext.empty())
            std::cout << "Output Extension: " << oext << std::endl;
    }
}

