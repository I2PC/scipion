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
#include <data/image.h>

class ProgWindow: public XmippMetadataProgram
{
public:
    typedef enum {CORNERMODE, SIZEMODE, CROPMODE} WindowMode;

    bool physical_coords;
    int sizeX, sizeY, sizeZ;
    int cropX, cropY, cropZ;
    int x0, y0, z0;
    int xF, yF, zF;
    double padValue;
    std::string padType;
    WindowMode mode;
    Image<char>                IChar;
    Image<unsigned char>       IUChar;
    Image<short int>           IShort;
    Image<unsigned short int>  IUShort;
    Image<float>               IFloat;

    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addUsageLine("This program takes a region from the input images");
        addUsageLine("Regions may be larger than the image itself, in that case");
        addUsageLine("padding is applied.");
        addParamsLine("  -corners <...>                        : Windows corners (by default indexes are logical)");
        addParamsLine("                                        : 2D: <x0> <y0> <xF> <yF>");
        addParamsLine("                                        : 3D: <x0> <y0> <z0> <xF> <yF> <zF>");
        addParamsLine("  or -size <sizeX> <sizeY=0> <sizeZ=0>  : Window to a new size");
        addParamsLine("                                        : if only one is given, the other two");
        addParamsLine("                                        : are supposed to be the same");
        addParamsLine("  or -crop <sizeX> <sizeY=0> <sizeZ=0>  : Crop this amount of pixels in each direction");
        addParamsLine("                                        : if only one is given, the other two");
        addParamsLine("                                        : are supposed to be the same");
        addParamsLine("  [-physical]                           : use physical instead of logical coordinates");
        addParamsLine("    requires -corners;");
        addParamsLine("  [-pad <padtype=value>]                : value used for padding");
        addParamsLine("   where <padtype>");
        addParamsLine("         value <v=0>                    : use this value for padding");
        addParamsLine("         corner                         : use the top-left corner for padding");
        addParamsLine("         avg                            : use the image average for padding");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        padType=getParam("-pad");
        padValue=getDoubleParam("-pad","value");
        if (checkParam("-corners"))
        {
            int nparams=getCountParam("-corners");
            if (nparams==4 && nparams==6)
            {
                x0=getIntParam("-corners",0);
                y0=getIntParam("-corners",1);
                if (nparams==4)
                {
                    xF=getIntParam("-corners",2);
                    yF=getIntParam("-corners",3);
                }
                else
                {
                    z0=getIntParam("-corners",2);
                    xF=getIntParam("-corners",3);
                    yF=getIntParam("-corners",4);
                    zF=getIntParam("-corners",5);
                }
            }
            else
                REPORT_ERROR(ERR_ARG_INCORRECT,"Incorrect number of arguments after -corners");
            physical_coords = checkParam("-physical");
            mode=CORNERMODE;
            // TODO Chequear el número de parámetros
        }
        else if (checkParam("-size"))
        {
            sizeX=getIntParam("-size",0);
            sizeY=getIntParam("-size",1);
            sizeZ=getIntParam("-size",2);
            if (sizeY==0)
                sizeY=sizeX;
            if (sizeZ==0)
                sizeZ=sizeX;
            x0 = FIRST_XMIPP_INDEX(sizeX);
            y0 = FIRST_XMIPP_INDEX(sizeY);
            z0 = FIRST_XMIPP_INDEX(sizeZ);
            xF = LAST_XMIPP_INDEX(sizeX);
            yF = LAST_XMIPP_INDEX(sizeY);
            zF = LAST_XMIPP_INDEX(sizeZ);
            mode=SIZEMODE;
            physical_coords = false;
        }
        else if (checkParam("-crop"))
        {
            cropX=getIntParam("-crop",0);
            cropY=getIntParam("-crop",1);
            cropZ=getIntParam("-crop",2);
            if (cropY==0)
                cropY=cropX;
            if (cropZ==0)
                cropZ=cropX;
            mode=CROPMODE;
            physical_coords = false;
        }
    }

    void show()
    {
        XmippMetadataProgram::show();
        switch (mode)
        {
        case SIZEMODE:
            std::cout << "New size: (XxYxZ)=" << sizeX << "x" << sizeY << "x"
            << sizeZ << std::endl;
            break;
        case CROPMODE:
            std::cout << "Crop: (XxYxZ)=" << cropX << "x" << cropY << "x"
            << cropZ << std::endl;
            break;
        case CORNERMODE:
            std::cout << "New window: from (z0,y0,x0)=(" << z0 << ","
            << y0 << "," << x0 << ") to (zF,yF,xF)=(" << zF << "," << yF
            << "," << xF << ")\n"
            << "Physical: " << physical_coords << std::endl;
        }
    }

    template <typename T>
    void processImage(Image<T> &Iin, const FileName &fnImg, const FileName &fnImgOut)
    {
        Iin.readApplyGeo(fnImg);
        double init_value(padValue);
        if (padType=="avg")
            init_value=Iin().computeAvg();
        else if (padType=="corner")
            init_value=DIRECT_MULTIDIM_ELEM(Iin(), 0);
        if (mode==CROPMODE)
        {
            int xl=cropX/2;
            int xr=cropX-xl;
            int yl=cropY/2;
            int yr=cropY-yl;
            int zl=cropZ/2;
            int zr=cropZ-zl;
            if (ZSIZE(Iin())==1)
            {
                zl=zr=0;
            }
            //call to a generic 4D function;
            Iin().window(0,  STARTINGZ(Iin())+zl, STARTINGY(Iin())+yl,  STARTINGX(Iin())+xl,
                         0, FINISHINGZ(Iin())-zr,FINISHINGY(Iin())-yr, FINISHINGX(Iin())-xr);
        }
        else
            if (!physical_coords)
                Iin().window(0, z0, y0, x0, 0, zF, yF,xF, init_value);
            else
                Iin().window(0,STARTINGZ(Iin()) + z0,
                             STARTINGY(Iin()) + y0,
                             STARTINGX(Iin()) + x0,
                             0,STARTINGZ(Iin()) + zF,
                             STARTINGY(Iin()) + yF, STARTINGX(Iin()) + xF,
                             init_value);
        Iin.write(fnImgOut);
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, long int objId)
    {
        Image<double> img;
        img.read(fnImg,false);

        switch (img.dataType())
        {
        case Unknown_Type:
            REPORT_ERROR(ERR_IMG_UNKNOWN,"");
            break;
        case UChar:
            processImage(IUChar,fnImg,fnImgOut);
            break;
        case SChar:
            processImage(IChar,fnImg,fnImgOut);
            break;
        case UShort:
            processImage(IUShort,fnImg,fnImgOut);
            break;
        case Short:
            processImage(IShort,fnImg,fnImgOut);
            break;
        default:
            processImage(IFloat,fnImg,fnImgOut);
            break;
        }
    }
};

int main(int argc, char **argv)
{
    ProgWindow prm;
    prm.read(argc, argv);
    prm.tryRun();
}
