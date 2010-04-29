/***************************************************************************
 *
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include <data/psf_xr.h>
#include <data/args.h>


void usage();

int main(int argc, char **argv)
{
    FileName fnPSF, fnPSFOut, fnImgIn, fnImgOut;
    XmippXRPSF psf;

#define DEBUG

    try
    {
        if (checkParameter(argc, argv, "-psf"))
        {
            fnPSF = getParameter(argc, argv, "-psf");
            psf.read(fnPSF);
        }
        else
            psf.clear();

    }
    catch (Xmipp_error &XE)
    {
        std::cerr << XE << std::endl;
        psf.usage();
        return 1;
    }

    try
    {
        fnImgIn = getParameter(argc, argv, "-i");

        if (checkParameter(argc, argv, "-out"))
            fnImgOut = getParameter(argc, argv, "-out");
        else
            fnImgOut = fnImgIn.without_extension().add_extension("out").add_extension("spi");

        psf.produceSideInfo();

        if (checkParameter(argc, argv, "-v"))
            std::cout << psf << std::endl;

        if (fnImgIn.get_extension()=="vol")
        {
            VolumeXmipp   phantomVol;
            ImageXmipp    imOut;

            phantomVol.read(fnImgIn);

            std::cout << phantomVol().zdim << std::endl;

            project_xr(psf, phantomVol, imOut);

            imOut.write(fnImgOut);


        }
        else if (fnImgIn.get_extension()=="spi")
        {
            ImageXmipp ImXmipp;
            Matrix2D  < double > ImgIn;

            ImXmipp.read(fnImgIn);
            ImXmipp().setXmippOrigin();

            ImgIn.resize(ImXmipp());

            FOR_ALL_ELEMENTS_IN_MATRIX2D(ImgIn)
            ImgIn(i,j) = ImXmipp(i,j);

            psf.generateOTF(ImgIn);

            psf.applyOTF(ImgIn);


            FOR_ALL_ELEMENTS_IN_MATRIX2D(ImgIn)
            ImXmipp(i,j) = abs(ImgIn(i,j));

            ImXmipp.write(fnImgOut);
        }
        else
            usage();


        if (checkParameter(argc, argv, "-psfout"))
        {
            fnPSFOut = getParameter(argc, argv, "-psfout");
            psf.write(fnPSFOut);
        }
    }
    catch (Xmipp_error &XE)
    {
        std::cerr << XE << std::endl;
        usage();
        return 1;
    }
    return 0;
}

void usage()
{
    std::cerr << "Usage: project_xr [options]\n"
              << "   -psf <PSF description file>      : PSF characteristic of the microscope \n"
              << "   -i <Input file>                  : Image or Volume \n"
              << "  [-out <Output file>]              : Resulting Image \n";
}

#undef DEBUG

