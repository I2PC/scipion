/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2007)
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

typedef enum
{
    IM2IM,
    MD2VOL,
    STK2VOL,
    VOL2STK,
    VOL2MD
} ImageConv;

class ConvImgProg: public XmippProgram
{
private:
    FileName fn_root, fn_oext, fn_in, fn_out;
    std::string type, depth;
    Image<char> in;
    Image<float> out;
    MetaData SF;
    MDRow    row;
    ImageConv convMode;
    bool adjust;
    int k;



protected:
    void defineParams()
    {
        addUsageLine("Convert among stacks, volumes and images, and change the file format.");
        addParamsLine(" -i <metadata>   : Input file: metadata, stack, volume or image.");
        addParamsLine("         :+ Supported read formats are:");
        addParamsLine("         :+ dm3 : Digital Micrograph 3.");
        addParamsLine("         :+ img : Imagic.");
        addParamsLine("         :+ inf,raw : RAW file with header INF file.");
        addParamsLine("         :+ mrc : CCP4.");
        addParamsLine("         :+ spe : Princeton Instruments CCD camera.");
        addParamsLine("         :+ spi, xmp : Spider.");
        addParamsLine("         :+ tif : TIFF.");
        addParamsLine("         :+ ser : tecnai imaging and analysis.");
        addParamsLine("         :+ raw#xDim,yDim,[zDim],offset,datatype,[r] : RAW image file without header file.");
        addParamsLine("         :+ where datatype can be: uint8,int8,uint16,int16,uint32,int32,long,float,double,cint16,cint32,cfloat,cdouble,bool");
        addParamsLine(" alias --input;");

        addParamsLine("  -o <output_file=\"\">  : Output file: metadata, stack, volume or image.");
        addParamsLine("   alias --output;");
        addParamsLine("OR  --oroot <root=\"\">     : Rootname of output individual images (Optional + \":ext\").");
        addParamsLine("  [--oext <extension=spi>] : Output file format extension.");
        addParamsLine("           where <extension>");
        addParamsLine("         img : Imagic (Data types: uint8, int16, float* and cfloat).");
        addParamsLine("         inf raw : RAW file with header INF file (All data types. See -d option).");
        addParamsLine("         mrc : CCP4 (Data types: int8, float* and cfloat).");
        addParamsLine("         spi xmp : Spider (Data types: float* and cfloat).");
        addParamsLine("         tif : TIFF. (Data types: uint8*, uint16, uint32 and float).");



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
        addParamsLine("                 cint16");
        addParamsLine("                 cint32");
        addParamsLine("                 cfloat");
        addParamsLine("                 cdouble");
        addParamsLine("                 bool");
        addParamsLine("  alias -d;");
        addParamsLine("  [--rangeAdjust] : Adjust the histogram to fill the gray level range.");
        addParamsLine("  alias -r;");
    }

    void readParams()
    {
        fn_in = getParam("-i");
        fn_out = getParam("-o");
        fn_root = getParam("--oroot");
        fn_oext = getParam("--oext");

        if (checkParam("--oroot") && !checkParam("--oext"))
            fn_oext = fn_root.getFileFormat();
        fn_root = fn_root.removeFileFormat();

        type = getParam("--type");

        if (!checkParam("--type"))
        {
            if (fn_out.getExtension() == "vol" || fn_oext == "vol")
                type = "vol";
            else if (fn_out.getExtension() == "stk" || fn_oext == "stk")
                type = "stk";
        }

        if (checkParam("--depth"))
        {
            depth = getParam("--depth");

            if (fn_out != "")
                fn_out += "%" + depth;
            else
                fn_oext += "%" + depth;
        }

        adjust = checkParam("--rangeAdjust");

    }

    template <typename T>
    void processImage(Image<T> &ImIn , const FileName &fnIn, const FileName &fnOut ,
                      bool readdata, int select_imgIn,
                      bool apply_geo, bool only_apply_shifts,
                      MDRow * row, bool mapData,
                      int select_imgOut, bool isStack,
                      int mode)
    {
        switch (convMode)
        {
        case IM2IM:
            {
                ImIn.read(fnIn,readdata, select_imgIn,apply_geo,only_apply_shifts,row,true);
                ImIn.write(fnOut,select_imgOut,isStack,mode,adjust);
                break;
            }
        case MD2VOL:
            {
                ImIn.read(fnIn,readdata, select_imgIn,apply_geo,only_apply_shifts,row,true);
                out().setSlice(k++,ImIn());
                break;
            }
        case STK2VOL:
            {
                ImIn.read(fnIn,readdata, select_imgIn,apply_geo,only_apply_shifts,row,true);
                ZSIZE(ImIn())*=NSIZE(ImIn());
                NSIZE(ImIn()) = 1;
                ImIn.write(fnOut,-1,false,WRITE_OVERWRITE,adjust);
                break;
            }
        case VOL2STK:
            {
                ImIn.read(fnIn,readdata, select_imgIn,apply_geo,only_apply_shifts,row,true);
                NSIZE(ImIn())=ZSIZE(ImIn());
                ZSIZE(ImIn())=1;
                int Ndim=NSIZE(ImIn());
                ImIn.MD.resize(Ndim);
                ImIn.write(fnOut,-1,true,WRITE_OVERWRITE,adjust);
                break;
            }
        case VOL2MD:
            {
                MetaData SFout;
                ImIn.read(fnIn,readdata, select_imgIn,apply_geo,only_apply_shifts,row,true);
                for (int k=0; k<ZSIZE(ImIn()); k++)
                {
                    ImIn().getSlice(k,out());
                    FileName fnOut;
                    if (fn_root !="")
                        fnOut.compose(fn_root,k,fn_oext);
                    else
                        fnOut.compose(fn_in.withoutExtension(),k,fn_oext);

                    out.write(fnOut,-1,false,WRITE_OVERWRITE,adjust);
                    SFout.addObject();
                    SFout.setValue(MDL_IMAGE,fnOut);
                    SFout.setValue(MDL_ENABLED,1);
                }
                if (fn_root != "")
                    SFout.write(fn_root+".sel");
                else if (fn_out != "")
                    SFout.write(fn_out);
                else
                    SFout.write(fn_in.withoutExtension() + "_out.sel");
                break;
            }
        }
    }


    void imageConvert( const FileName &fnIn, const FileName &fnOut ,
                       bool readdata=true, int select_imgIn = -1,
                       bool apply_geo = false, bool only_apply_shifts = false,
                       MDRow * row = NULL, bool mapData = false,
                       int select_imgOut=-1, bool isStack=false,
                       int mode=WRITE_OVERWRITE)
    {
        Image<char> Im;
        Im.read(fnIn,false);


        switch (Im.dataType())
        {
        case Unknown_Type:
            REPORT_ERROR(ERR_IMG_UNKNOWN,"");
            break;
        case UChar:
            {
                Image<unsigned char>       IUChar;
                processImage(IUChar,fnIn,fnOut,readdata,select_imgIn,apply_geo,\
                             only_apply_shifts,row,mapData,select_imgOut,isStack,mode);
                break;
            }
        case SChar:
            {
                Image<char>                IChar;
                processImage(IChar,fnIn,fnOut,readdata,select_imgIn,apply_geo,\
                             only_apply_shifts,row,mapData,select_imgOut,isStack,mode);
                break;
            }
        case UShort:
            {
                Image<unsigned short int>  IUShort;
                processImage(IUShort,fnIn,fnOut,readdata,select_imgIn,apply_geo,\
                             only_apply_shifts,row,mapData,select_imgOut,isStack,mode);
                break;
            }
        case Short:
            {
                Image<short int>           IShort;
                processImage(IShort,fnIn,fnOut,readdata,select_imgIn,apply_geo,\
                             only_apply_shifts,row,mapData,select_imgOut,isStack,mode);
                break;
            }
        case Int:
            {
                Image<int>                  IInt;
                processImage(IInt,fnIn,fnOut,readdata,select_imgIn,apply_geo,\
                             only_apply_shifts,row,mapData,select_imgOut,isStack,mode);
                break;
            }
        case Float:
            {
                Image<float>               IFloat;
                processImage(IFloat,fnIn,fnOut,readdata,select_imgIn,apply_geo,\
                             only_apply_shifts,row,mapData,select_imgOut,isStack,mode);
                break;
            }
        default:
            REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Datatype not implemented.");
            break;
        }
    }

public:
    void run()
    {
        convMode = IM2IM;

        // From metadata to ...............................................
        if (fn_in.isMetaData())
        {
            MDRow    row;
            SF.read(fn_in);
            if (type == "stk")
            {
                FileName fn_stack_plain=fn_out.removeFileFormat();
                if (exists(fn_stack_plain))
                    unlink(fn_stack_plain.c_str());
                FOR_ALL_OBJECTS_IN_METADATA(SF)
                {
                    FileName fnImg;
                    SF.getValue(MDL_IMAGE,fnImg);
                    SF.getRow(row);
                    //                    in.read(fnImg,true,-1,false,false,&row);
                    //                    in.write(fn_out,-1,true,WRITE_APPEND);
                    imageConvert(fnImg,fn_out,true,-1,false,false,&row,true,-1,true, WRITE_APPEND);
                }
            }
            else if (type == "vol")
            {
                convMode = MD2VOL;

                int Xdim, Ydim, Zdim;
                unsigned long Ndim;
                ImgSize(SF, Xdim, Ydim, Zdim, Ndim);
                if (Zdim!=1 || Ndim!=1)
                    REPORT_ERROR(ERR_MULTIDIM_DIM,
                                 "Only 2D images can be converted into volumes");
                out().coreAllocate(1,SF.size(),Ydim,Xdim);
                k = 0;
                FOR_ALL_OBJECTS_IN_METADATA(SF)
                {
                    FileName fnImg;
                    SF.getValue(MDL_IMAGE,fnImg);
                    imageConvert(fnImg,fn_out);
                }
                out.write(fn_out);
            }
            else if (type == "img")
            {
                int k=0;
                MetaData SFout;
                FOR_ALL_OBJECTS_IN_METADATA(SF)
                {
                    FileName fnImg;
                    SF.getValue(MDL_IMAGE,fnImg);
                    SF.getRow(row);
                    FileName fnOut;
                    if (fn_root !="")
                        fnOut.compose(fn_root,k++,fn_oext);
                    else
                        fnOut = (fnImg.withoutExtension()).addExtension(fn_oext);

                    //                    in.read(fnImg,true,-1,false,false,&row);
                    //                    in.write(fnOut);
                    imageConvert(fnImg,fnOut,true,-1,false,false,&row);

                    SFout.addObject();
                    SFout.setValue(MDL_IMAGE,fnOut);
                    SFout.setValue(MDL_ENABLED,1);
                }
                if (fn_root != "")
                    SFout.write(fn_root+".sel");
                else if (fn_out != "")
                    SFout.write(fn_out);
                else
                    SFout.write(fn_in.insertBeforeExtension("_out"));
            }
        }
        else // From Image file to.....
        {
            in.read(fn_in,false);
            if (NSIZE(in())>1)   // It's a stack with more than 1 slice

            {
                if (type == "stk")
                {
                    FileName fn_stack_plain=fn_out.removeFileFormat();
                    if (exists(fn_stack_plain))
                        unlink(fn_stack_plain.c_str());
                    int nmax=NSIZE(in());
                    for (int n=0; n<nmax; n++)
                    {
                        //                        in.read(fn_in,true,n);
                        //                        in.write(fn_out,-1,true,WRITE_APPEND);
                        imageConvert(fn_in,fn_out,true,n,false,false,NULL,true,-1,true, WRITE_APPEND);
                    }
                }
                else if (type == "vol")
                {
                    convMode = STK2VOL;
                    imageConvert(fn_in,fn_out);
                }
                else if (type == "img")
                {
                    MetaData SFout;
                    int nmax=NSIZE(in());
                    for (int n=0; n<nmax; n++)
                    {
                        FileName fnOut;
                        if (fn_root !="")
                            fnOut.compose(fn_root,n,fn_oext);
                        else
                            fnOut.compose(fn_in.withoutExtension(),n,fn_oext);

                        //                        in.read(fn_in,true,n);
                        //                        in.write(fnOut);
                        imageConvert(fn_in,fnOut,true,n);

                        SFout.addObject();
                        SFout.setValue(MDL_IMAGE,fnOut);
                        SFout.setValue(MDL_ENABLED,1);
                    }
                    if (fn_root != "")
                        SFout.write(fn_root+".sel");
                    else if (fn_out != "")
                        SFout.write(fn_out);
                    else
                        SFout.write(fn_in.withoutExtension() + "_out.sel");
                }
            }
            else // It's a stack with 1 slice, an image or a volume
            {
                if (type == "stk")
                {
                    if (ZSIZE(in())>1)   // Convert the volume into a stack
                        convMode = VOL2STK;
                    else
                        convMode = IM2IM;
                }
                else if (ZSIZE(in())>1 && type == "img")
                    convMode = VOL2MD;
                else
                    convMode = IM2IM;

                imageConvert(fn_in,fn_out);
            }
        }
    }
};


int main(int argc, char *argv[])
{
    try
    {
        ConvImgProg program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }

    return 0;
}
