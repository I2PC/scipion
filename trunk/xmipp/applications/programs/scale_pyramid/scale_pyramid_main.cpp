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

#include <data/program.h>
#include <data/args.h>

class ProgPyramid: public XmippMetadataProgram
{
public:
    enum Toperation {Expand, Reduce};
    Toperation operation;
    int levels;
    float scaleFactor;
    Image<double> result;

    void readParams()
    {
        XmippMetadataProgram::readParams();
        levels = getIntParam("-levels");
        if (checkParam("-expand"))
            operation = Expand;
        else if (checkParam("-reduce"))
            operation = Reduce;
    }

    void show()
    {
        if (!verbose)
            return;
        XmippMetadataProgram::show();
        std::cout << "Operation: ";
        switch (operation)
        {
        case Expand:
            std::cout << "Expand\n";
            break;
        case Reduce:
            std::cout << "Reduce\n";
            break;
        }
        std::cout << "Levels: " << levels << std::endl;
    }

    void defineParams()
    {
    	each_image_produces_an_output=true;
    	allow_time_bar=true;
    	XmippMetadataProgram::defineParams();
        addParamsLine("  -expand or -reduce       : Expand or reduce the image");
        addParamsLine(" [-levels <l=1>]           : Expansion/reduction factor");
    }

    void preProcess() {
    	scaleFactor = (float)(pow(2.0, levels));
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, long int objId)
    {
    	Image<double> img;
    	img.readApplyGeo(fnImg,mdIn,objId);
        float Xoff, Yoff, Zoff;
        Xoff=img.Xoff();
        Yoff=img.Yoff();
        Zoff=img.Zoff();
        switch (operation)
        {
        case Expand:
            pyramidExpand(BSPLINE3, result(),img(),levels);
            img.setXoff(Xoff*scaleFactor);
            img.setYoff(Yoff*scaleFactor);
            img.setZoff(Zoff*scaleFactor);
            break;
        case Reduce:
            pyramidReduce(BSPLINE3,result(),img(),levels);
            img.setXoff(Xoff/scaleFactor);
            img.setYoff(Yoff/scaleFactor);
            img.setZoff(Zoff/scaleFactor);
            break;
        }
        result.write(fnImgOut);
    }
};

int main(int argc, char **argv)
{
    ProgPyramid prm;
	prm.read(argc,argv);
    try {
    	prm.run();
    } catch (XmippError XE)
    {
    	std::cerr << XE << std::endl;
    	return 1;
    }
    return 0;
}
