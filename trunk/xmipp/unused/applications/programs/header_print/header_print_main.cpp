/***************************************************************************
 *
 * Authors:    Sjors Scheres
 *             Slavica Jonic
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

/* PROGRAM ----------------------------------------------------------------- */

class ProgHeaderPrint: public XmippMetadataProgram
{
protected:
    void defineParams()
    {
        each_image_produces_an_output = false;
        apply_geo = false;
        allow_time_bar = false;
        decompose_stacks = false;
        XmippMetadataProgram::defineParams();
        addUsageLine("Print information from the header of 2D-images.");
        addExampleLine("Print header values in the second image in stack input/smallStack.stk",false);
        addExampleLine("xmipp_header_print -i 2@input/smallStack.stk");
        addExampleLine("Print header values in  image kk.img",false);
        addExampleLine("xmipp_header_print -i kk.img");
        addSeeAlsoLine("header_extract, header_assign");
    }

    void show()
    {
        std::cout << " Printing the header ... " << std::endl;
        std::cout << "***** WARNING: the information stored in header is ignored inside Xmipp, it is only relevant" << std::endl
                  << "***** where importing/exporting data from/to other packages." << std::endl;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        int Xdim, Ydim, Zdim;
        unsigned long Ndim;
        getImageSize(fnImg, Xdim, Ydim, Zdim, Ndim);

        Image<double> img;
        img.read(fnImg, HEADER);
        std::cout << img << std::endl;
    }
}
;// end of class ProgHeaderPrint

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgHeaderPrint program;
    program.read(argc, argv);
    program.tryRun();
 }

