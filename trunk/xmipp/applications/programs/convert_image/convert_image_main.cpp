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

class ProgConvImg: public XmippProgram
{
private:
    FileName fn_root, fn_oext, fn_in, fn_out;
    std::string type, bits;
    Image<float> in, out;
    MetaData SF;
    MDRow    row;

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
        addParamsLine("         :+ where datatype can be: uchar, char, ushort, short, uint, int, long, float, double, cshort, cint, cfloat, cdouble, bool");
        addParamsLine(" alias --input;");

        addParamsLine("  [-o <output_file=\"\">]  : Output file: metadata, stack, volume or image.");
        addParamsLine("   alias --output;");
        addParamsLine("  [--oext <extension=spi>] : Output file format extension.");
        addParamsLine("         :+ Supported write formats are:");
        addParamsLine("         :+ img : Imagic");
        addParamsLine("         :+ inf,raw : RAW file with header INF file.");
        addParamsLine("         :+ mrc : CCP4");
        addParamsLine("         :+ spi, xmp : Spider");
        addParamsLine("         :+ tif : TIFF. It supports 8bits, 16bits and float. See -bits option.");
        addParamsLine("  [--oroot <root=\"\">]     : Rootname of output individual images.");

        addParamsLine("  [--type <output_type=img>] : Output file type.");
        addParamsLine("          where <output_type>");
        addParamsLine("          img : Image");
        addParamsLine("          vol : Volume");
        addParamsLine("          stk : Stack ");
        addParamsLine("  alias -t;");

        addParamsLine("  [--bits+ <bit_depth=8>] : Bit depth for TIFF format. Options are: ");
        addParamsLine("                         :   8  : Uint8 (char).");
        addParamsLine("                         :   16 : Uint16 (short).");
        addParamsLine("                         :   32 : Float.");
        addParamsLine("  alias -b;");
    }

    void readParams()
    {
        fn_in = getParam("-i");
        fn_out = getParam("-o");
        fn_oext = getParam("--oext");
        fn_root = getParam("--oroot");

        type = getParam("--type");

        if (!checkParam("--type"))
        {
            if (fn_out.getExtension() == "vol" || fn_oext == "vol")
                type = "vol";
            else if (fn_out.getExtension() == "stk" || fn_oext == "stk")
                type = "stk";
        }

        bits = getParam("--bits");

        if (fn_out.getExtension() == "tif")
          fn_out += "%" + bits;
        else if (fn_oext == "tif")
          fn_oext += "%" + bits;
        else if (checkParam("--bits"))
            REPORT_ERROR(ERR_PARAM_INCORRECT, "-bits option is only valid for TIFF format.");

    }
public:
    void run()
    {
        // From metadata to ...............................................
        if (fn_in.isMetaData())
        {
            MetaData SF;
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
                    in.read(fnImg,true,-1,false,false,&row);
                    in.write(fn_out,-1,true,WRITE_APPEND);
                }
            }
            else if (type == "vol")
            {
                int Xdim, Ydim, Zdim, Ndim;
                ImgSize(SF, Xdim, Ydim, Zdim, Ndim);
                if (Zdim!=1 || Ndim!=1)
                    REPORT_ERROR(ERR_MULTIDIM_DIM,
                                 "Only 2D images can be converted into volumes");
                out().coreAllocate(1,SF.size(),Ydim,Xdim);
                int k=0;
                FOR_ALL_OBJECTS_IN_METADATA(SF)
                {
                    FileName fnImg;
                    SF.getValue(MDL_IMAGE,fnImg);
                    in.read(fnImg);
                    out().setSlice(k++,in());
                }
                out.write(fn_out);
            }
            else if (type == "img")
            {
                fn_root=fn_root.removeFileFormat();
                int k=0;
                MetaData SFout;
                FOR_ALL_OBJECTS_IN_METADATA(SF)
                {
                    FileName fnImg;
                    SF.getValue(MDL_IMAGE,fnImg);
                    SF.getRow(row);
                    in.read(fnImg,true,-1,false,false,&row);
                    FileName fnOut;
                    if (fn_root !="")
                        fnOut.compose(fn_root,k++,fn_oext);
                    else
                        fnOut = (fnImg.withoutExtension()).addExtension(fn_oext);

                    in.write(fnOut);
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
        else // From Stack file to.....
        {
            in.read(fn_in,false);
            if (NSIZE(in())>1)
            {
                // It's a stack with more than 1 slice
                if (type == "stk")
                {
                    FileName fn_stack_plain=fn_out.removeFileFormat();
                    if (exists(fn_stack_plain))
                        unlink(fn_stack_plain.c_str());
                    int nmax=NSIZE(in());
                    for (int n=0; n<nmax; n++)
                    {
                        in.read(fn_in,true,n);
                        in.write(fn_out,-1,true,WRITE_APPEND);
                    }
                }
                else if (type == "vol")
                {
                    in.read(fn_in);
                    ZSIZE(in())*=NSIZE(in());
                    NSIZE(in())=1;
                    in.write(fn_out);
                }
                else if (type == "img")
                {
                    fn_root=fn_root.removeFileFormat();
                    MetaData SFout;
                    int nmax=NSIZE(in());
                    for (int n=0; n<nmax; n++)
                    {
                        in.read(fn_in,true,n);
                        FileName fnOut;
                        if (fn_root !="")
                            fnOut.compose(fn_root,n,fn_oext);
                        else
                            fnOut.compose(fn_in.withoutExtension(),n,fn_oext);
                        in.write(fnOut);
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
            else
            {
                // It's a stack with 1 slice, an image or a volume
                in.read(fn_in);
                if (type == "stk")
                {
                    if (ZSIZE(in())>1)
                    {
                        // Convert the volume into a stack
                        NSIZE(in())=ZSIZE(in());
                        ZSIZE(in())=1;
                        int Ndim=NSIZE(in());
                        in.MD.resize(Ndim);
                    }
                    in.write(fn_out,-1,true);
                }
                else if (ZSIZE(in())>1 && type == "img")
                {
                    fn_root=fn_root.removeFileFormat();
                    MetaData SFout;
                    for (int k=0; k<ZSIZE(in()); k++)
                    {
                        in().getSlice(k,out());
                        FileName fnOut;
                        if (fn_root !="")
                            fnOut.compose(fn_root,k,fn_oext);
                        else
                            fnOut.compose(fn_in.withoutExtension(),k,fn_oext);

                        out.write(fnOut);
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
                else
                    in.write(fn_out);
            }
        }
    }
}
;

int main(int argc, char *argv[])
{
    try
    {
        ProgConvImg program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }

    return 0;
}
