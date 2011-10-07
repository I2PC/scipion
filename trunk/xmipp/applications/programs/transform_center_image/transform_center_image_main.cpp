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

#include <data/xmipp_program.h>
#include <data/filters.h>

class ProgCenterImage: public XmippMetadataProgram
{
public:
    int Niter;
    bool limitShift;

    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        //usage
        addUsageLine("Center a set of images (preferably with a good SNR, e.g., class averages).");
        addUsageLine("After calling the program, all classes will be centered and they will tend ");
        addUsageLine("to show the same orientation. The program centers the images by comparing ");
        addUsageLine("the image with its X, Y, and XY mirrors. The orientation is determined by ");
        addUsageLine("comparing the image with its X mirror. ");
        //examples
        addExampleLine("Center images in smallStack.stk and store results in a different file:", false);
        addExampleLine("xmipp_transform_center_image -i smallStack.stk -o stackCentered.stk ");
        //params
        addParamsLine("  [--iter <n=10>]      : Number of iterations");
        addParamsLine("  [--limit]            : Limit the maximum shift allowed");

    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        Niter = getIntParam("--iter");
        limitShift = checkParam("--limit");
    }

    void show()
    {
        XmippMetadataProgram::show();
        std::cout << "iterations = " << Niter << std::endl;
        std::cout << "limit shift = " << limitShift << std::endl;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
        Image<double> img;
        img.readApplyGeo(fnImg, rowIn);
        img().checkDimensionWithDebug(2,__FILE__,__LINE__);
        centerImage(img(), Niter, limitShift);
        img.write(fnImgOut);
    }

}
;///end of class ProgCenterImage

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgCenterImage program;
    program.read(argc, argv);
    program.tryRun();
}

