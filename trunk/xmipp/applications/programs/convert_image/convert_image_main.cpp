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

#include <data/args.h>
#include <data/image.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/progs.h>

#include <data/image_generic.h>

typedef enum
{
    MD2MD,
    MD2VOL,
    VOL2MD
} ImageConv;

class ProgConvImg: public XmippMetadataProgram
{
private:
    std::string  type, depth;
    Image<char>  imTemp;
    ImageGeneric imIn, *imOut;
    DataType     outDataT;
    MDRow        row;
    ImageConv    convMode;
    bool         adjust;
    int          k;

protected:
    void defineParams()
    {
        each_image_produces_an_output = true;
        save_metadata_stack = true;
        XmippMetadataProgram::defineParams();
        addUsageLine("Convert among stacks, volumes and images, and change the file format.");
        addParamsLine("  [--oext <extension=\"\">] :  Output file format extension.");
        addWhereImageFormat("extension");
        addParamsLine("  [--type <output_type=img>] : Output file type.");
        addParamsLine("          where <output_type>");
        addParamsLine("          img : Image");
        addParamsLine("          vol : Volume");
        addParamsLine("          stk : Stack ");
        addParamsLine("  alias -t;");
        addParamsLine("  [--depth+ <bit_depth=default>] : Image bit depth.");
        addParamsLine("          where <bit_depth>");
        addParamsLine("                 default: Default selected value (*).");
        addParamsLine("                 uint8 : Equivalent to uchar");
        addParamsLine("                 int8  : Equivalent to char");
        addParamsLine("                 uint16: Equivalent to ushort");
        addParamsLine("                 int16 : Equivalent to short");
        addParamsLine("                 uint32");
        addParamsLine("                 int32");
        addParamsLine("                 long");
        addParamsLine("                 float");
        addParamsLine("                 double");
        addParamsLine("                 cint16 : Complex in16");
        addParamsLine("                 cint32 : Complex in32");
        addParamsLine("                 cfloat : Complex float");
        addParamsLine("                 cdouble: Complex double");
        addParamsLine("                 bool");
        addParamsLine("  alias -d;");
        addParamsLine("  [--rangeAdjust] : Adjust the histogram to fill the gray level range.");
        addParamsLine("  alias -r;");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();

        fn_out = (checkParam("-o"))? getParam("-o") : "";
        adjust = checkParam("--rangeAdjust");
        type = getParam("--type");

        oext = checkParam("--oext") ? getParam("--oext") : "";
        if (oext == "custom")
            oext = getParam("--oext",1);

        if (!checkParam("--type"))
        {
            if (fn_out.getExtension() == "vol" || oext == "vol")
                type = "vol";
            else if (fn_out.getExtension() == "stk" || oext == "stk")
                type = "stk";
        }

        if (checkParam("--depth"))
        {
            depth = getParam("--depth");

            outDataT = datatypeString2Int((depth == "default")? "float": depth);

            if (fn_out != "")
                fn_out += "%" + depth;
            else
                oext += "%" + depth;
        }
        else
            outDataT = Float;
    }

    void preProcess()
    {
        FileName fn_stack_plain=fn_out.removeFileFormat();
        if (exists(fn_stack_plain) && type == "stk")
            unlink(fn_stack_plain.c_str());

        convMode = MD2MD;

        if (!single_image && type == "vol")
        {
            convMode = MD2VOL;

            int Xdim, Ydim, Zdim;
            unsigned long Ndim;
            ImgSize(mdIn, Xdim, Ydim, Zdim, Ndim);
            if (Zdim!=1)
                REPORT_ERROR(ERR_MULTIDIM_DIM,
                             "Only 2D images can be converted into volumes");
            imOut = new ImageGeneric(outDataT);
            imOut->newMappedFile(Xdim,Ydim,mdIn.size(),1,fn_out);
            k = 0;
        }
        else if (single_image)
        {
            int Xdim, Ydim, Zdim;
            unsigned long Ndim;
            imTemp.read(fn_in,false);
            imTemp.getDimensions(Xdim,Ydim,Zdim,Ndim);

            //Fill mdIn to allow XmippMetaDataProgram create the fnImgOut
            if (single_image && Zdim > 1 && type != "vol")
            {
                convMode = VOL2MD;
                single_image = false;

                mdIn.clear();

                FileName fnTemp;

                for (k=0;k<Zdim;k++)
                {
                    fnTemp.compose(k, fn_in);

                    size_t id = mdIn.addObject();
                    mdIn.setValue(MDL_IMAGE, fnTemp, id);
                    mdIn.setValue(MDL_ENABLED, 1, id);
                }
                imIn.read(fn_in,true,-1,true);
                imOut = new ImageGeneric(outDataT);
                k = 0;
            }
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        switch(convMode)
        {
        case MD2MD:
            {
                imIn.read(fnImg,true,-1,true);
                imIn.write(fnImgOut,-1, type == "stk",WRITE_APPEND,adjust);
                break;
            }
        case MD2VOL:
            {
                imIn.read(fnImg,true,-1,true);
                imOut->data->setSlice(k++,imIn.data);
                break;
            }
        case VOL2MD:
            {
                imIn.data->getSlice(k++,imOut->data);
                imOut->write(fnImgOut,-1,type == "stk",WRITE_APPEND,adjust);
            }
        }
        mdIn.setValue(MDL_IMAGE,fnImgOut, objId); // to keep info in output metadata
    }

    void finishProcessing()
    {
        switch(convMode)
        {
        case MD2VOL:
            imOut->write();
            break;
        }

        // To keep the mdIn info we overwrite the mdOut done by XmippMetadataProgram
        mdOut = mdIn;
        XmippMetadataProgram::finishProcessing();
    }

    void show()
    {
        XmippMetadataProgram::show();
        if (each_image_produces_an_output)
        {
            if (oext != "")
                std::cout << "Output Extension: " << oext << std::endl;
        }
    }
};

int main(int argc, char *argv[])
{
    ProgConvImg program;
    program.read(argc, argv);
    program.tryRun();

    return 0;
}
