
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Pedro A. de Alarcon (pedro@cnb.csic.es)
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

#include <data/xmipp_program.h>
#include <data/morphology.h>
#include <data/filters.h>

class ProgMorphology: public XmippMetadataProgram
{
public:
#define DILATION     1
#define EROSION      2
#define OPENING      3
#define CLOSING      4
#define SHARPENING   5
#define KEEPBIGGEST  6
#define REMOVESMALL  7

	bool binaryOperation;
    int operation;
    int size;
    int count;
    int neig2D, neig3D;
    double width;
    double strength;
    int smallSize;
public:
    void defineParams()
    {
        each_image_produces_an_output = true;
        addUsageLine("Apply morphological operations to binary or gray images/volumes.");
        addUsageLine("+You may learn about morphological operations in general from");
        addUsageLine("[[http://en.wikipedia.org/wiki/Mathematical_morphology][here]].");
        addUsageLine("");
        XmippMetadataProgram::defineParams();
        addParamsLine("--binaryOperation <op>: Morphological operation on binary images");
        addParamsLine("       where <op>");
        addParamsLine("             dilation   : Dilate white region");
        addParamsLine("             erosion    : Erode white region");
        addParamsLine("             closing    : Dilation+Erosion, removes black spots");
        addParamsLine("             opening    : Erosion+Dilation, removes white spots");
        addParamsLine("             keepBiggest : Keep biggest component");
        addParamsLine("             removeSmall <size=10> : Remove components whose size is smaller than this size");
        addParamsLine("or --grayOperation <op>: Morphological operation on gray images");
        addParamsLine("       where <op>");
        addParamsLine("             sharpening <w=1> <s=0.5>: Sharpening with width (suggested 1 or 2)");
        addParamsLine("                                     : and strength (suggested 0.1-1.0).");
        addParamsLine("                                     : Only valid for volumes, not images.");
        addParamsLine("                                     :++ Implemented according to JGM Schavemaker, MJT Reinders, JJ Gerbrands,");
        addParamsLine("                                     :++ E Backer. Image sharpening by morphological filtering. Pattern Recognition");
        addParamsLine("                                     :++ 33: 997-1012 (2000)");
        addParamsLine("[--neigh2D+ <n=Neigh8>] : Neighbourhood in 2D.");
        addParamsLine("          where <n>");
        addParamsLine("                 Neigh4");
        addParamsLine("                 Neigh8");
        addParamsLine("     requires --binaryOperation;");
        addParamsLine("[--neigh3D+ <n=Neigh18>] : Neighbourhood in 3D.");
        addParamsLine("          where <n>");
        addParamsLine("                 Neigh6");
        addParamsLine("                 Neigh18");
        addParamsLine("                 Neigh26");
        addParamsLine("     requires --binaryOperation;");
        addParamsLine("[--size <s=1>]: Size of the Strutural element.");
        addParamsLine("     requires --binaryOperation;");
        addParamsLine("[--count+ <c=0>]: Minimum required neighbors with distinct value.");
        addParamsLine("     requires --binaryOperation;");
        addExampleLine("xmipp_transform_morphology -i binaryVolume.vol --binaryOperation dilation");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        binaryOperation=checkParam("--binaryOperation");
        if (binaryOperation)
        {
            String strOperation=getParam("--binaryOperation");
            if (strOperation=="dilation")
                operation = DILATION;
            else if (strOperation=="erosion")
                operation = EROSION;
            else if (strOperation=="closing")
                operation = CLOSING;
            else if (strOperation=="opening")
                operation = OPENING;
            else if (strOperation=="keepBiggest")
                operation = KEEPBIGGEST;
            else if (strOperation=="removeSmall")
            {
                operation = REMOVESMALL;
                smallSize = getIntParam("--binaryOperation",1);
            }

            size = getIntParam("--size");
            String neighbourhood=getParam("--neigh2D");
            if (neighbourhood=="Neigh8")
            	neig2D = 8;
            else
            	neig2D = 4;
            neighbourhood=getParam("--neigh3D");
            if (neighbourhood=="Neigh6")
            	neig3D = 6;
            else if (neighbourhood=="Neigh18")
            	neig3D = 18;
            else
            	neig3D = 26;
            count = getIntParam("--count");
        }
        else if (checkParam("--grayOperation"))
        {
            String strOperation=getParam("--grayOperation");
            if (strOperation=="sharpening")
            {
                operation = SHARPENING;
                width=getDoubleParam("--grayOperation",1);
                strength=getDoubleParam("--grayOperation",2);
            }
        }
    }

    void show()
    {
        std::cout << "Performing a ";
        switch (operation)
        {
        case DILATION       :
            std::cout << "Dilation\n";
            break;
        case EROSION       :
            std::cout << "Erosion\n";
            break;
        case OPENING       :
            std::cout << "Opening\n";
            break;
        case CLOSING       :
            std::cout << "Closing\n";
            break;
        case SHARPENING    :
            std::cout << "Sharpening\n"
            << "Width = " << width << std::endl
            << "Strength = " << strength << std::endl;
        case KEEPBIGGEST   :
        	std::cout << "Keeping biggest component\n";
        	break;
        case REMOVESMALL   :
        	std::cout << "Removing small objects\n"
        	<< "Small size < " << smallSize << std::endl;
        	break;
        }
        if (binaryOperation)
            std::cout << "Size=" << size << std::endl
            << "Neighbourhood2D=" << neig2D << std::endl
            << "Neighbourhood3D=" << neig3D << std::endl
            << "Count=" << count << std::endl;
    }
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
        Image<double> img, imgOut;
        img.readApplyGeo(fnImg, rowIn);
        imgOut().resizeNoCopy(img());

        if (binaryOperation)
            std::cout << "Initially the image has " << img().sum()
            << " pixels set to 1\n";
        bool isVolume=ZSIZE(img())>1;
        switch (operation)
        {
        case DILATION:
            if (isVolume)
                dilate3D(img(), imgOut(), neig3D, count, size);
            else
                dilate2D(img(), imgOut(), neig2D, count, size);
            break;
        case EROSION:
            if (isVolume)
                erode3D(img(), imgOut(), neig3D, count, size);
            else
                erode2D(img(), imgOut(), neig2D, count, size);
            break;
        case OPENING:
            if (isVolume)
                opening3D(img(), imgOut(), neig3D, count, size);
            else
                opening2D(img(), imgOut(), neig2D, count, size);
            break;
        case CLOSING:
            if (isVolume)
                closing3D(img(), imgOut(), neig3D, count, size);
            else
                closing2D(img(), imgOut(), neig2D, count, size);
            break;
        case SHARPENING:
            if (isVolume)
                sharpening(img(), width, strength, imgOut());
            else
                REPORT_ERROR(ERR_NOT_IMPLEMENTED,"Sharpening has not been implemented for images");
            break;
        case KEEPBIGGEST:
        	if (isVolume)
        		keepBiggestComponent(img(),0,neig3D);
        	else
        		keepBiggestComponent(img(),0,neig2D);
        	imgOut()=img();
        	break;
        case REMOVESMALL:
        	if (isVolume)
        		removeSmallComponents(img(),smallSize,neig3D);
        	else
        		removeSmallComponents(img(),smallSize,neig2D);
        	imgOut()=img();
        	break;
        }

        if (binaryOperation)
            std::cout << "Finally the image has " << imgOut().sum() << " pixels set to 1\n";
        imgOut.write(fnImgOut);
    }
};

int main(int argc, char **argv)
{
    ProgMorphology prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
