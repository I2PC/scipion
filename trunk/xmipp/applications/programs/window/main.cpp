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

class Window_parameters: public Prog_parameters
{
public:
    bool size_mode;
    bool crop_mode;
    bool physical_coords;
    int sizeX, sizeY, sizeZ;
    int cropX, cropY, cropZ;
    int x0, y0, z0;
    int xF, yF, zF;
    double init_value;
    int wrong_parameters;
    int average_pad;
    int corner_pad;

    void read(int argc, char **argv)
    {
        init_value = 0.;
        wrong_parameters = 0;
        average_pad = false;
        corner_pad = false;
        size_mode=false;
        crop_mode=false;

        Prog_parameters::read(argc, argv);

        if (checkParameter(argc, argv, "-pad_value"))
        {
            init_value = textToFloat(getParameter(argc, argv, "-pad_value"));
            wrong_parameters = 1;
        }
        if (checkParameter(argc, argv, "-corner_pad_value"))
        {
            corner_pad = true;
            if (wrong_parameters > 0)  wrong_parameters = -1;
            else  wrong_parameters = 1;
        }
        if (checkParameter(argc, argv, "-average_pad_value"))
        {
            average_pad = true;
            if (wrong_parameters > 0)  wrong_parameters = -1;
            else  wrong_parameters = 1;
        }
        if (wrong_parameters == -1)
            REPORT_ERROR(1, "incompatible options");

        if (checkParameter(argc, argv, "-size"))
        {
            size_mode = true;
            int i = paremeterPosition(argc, argv, "-size");
            if (i + 2 >= argc)
            {
                sizeZ = sizeY = sizeX = textToInteger(argv[i+1]);
            }
            else if (i + 3 >= argc)
            {
                sizeZ = sizeY = textToInteger(argv[i+2]);
                sizeX = textToInteger(argv[i+1]);
            }
            else if (i + 4 >= argc)
            {
                sizeZ = textToInteger(argv[i+3]);
                sizeY = textToInteger(argv[i+2]);
                sizeX = textToInteger(argv[i+1]);
            }
            else REPORT_ERROR(1, "Not enough parameters after -size");

            x0 = FIRST_XMIPP_INDEX(sizeX);
            y0 = FIRST_XMIPP_INDEX(sizeY);
            z0 = FIRST_XMIPP_INDEX(sizeZ);
            xF = LAST_XMIPP_INDEX(sizeX);
            yF = LAST_XMIPP_INDEX(sizeY);
            zF = LAST_XMIPP_INDEX(sizeZ);
            physical_coords = false;
        }
        else if (checkParameter(argc, argv, "-r0"))
        {
            size_mode = false;
            // Get r0
            int i = paremeterPosition(argc, argv, "-r0");
            if (i + 2 >= argc) REPORT_ERROR(1, "Not enough parameters after -r0");
            else
            {
                x0 = textToInteger(argv[i+1]);
                y0 = textToInteger(argv[i+2]);
            }
            if (i + 3 < argc)
                try
                {
                    z0 = textToInteger(argv[i+3]);
                }
                catch (Xmipp_error XE)
                {
                    z0 = 1;
                }

            // Get rF
            i = paremeterPosition(argc, argv, "-rF");
            if (i == -1) REPORT_ERROR(1, "-rF not present");
            if (i + 2 >= argc) REPORT_ERROR(1, "Not enough parameters after -rF");
            else
            {
                xF = textToInteger(argv[i+1]);
                yF = textToInteger(argv[i+2]);
            }
            if (i + 3 < argc)
                try
                {
                    zF = textToInteger(argv[i+3]);
                }
                catch (Xmipp_error XE)
                {
                    zF = 1;
                }

            physical_coords = checkParameter(argc, argv, "-physical");
        }
        else if (checkParameter(argc, argv, "-crop"))
        {
            crop_mode = true;
            int i = paremeterPosition(argc, argv, "-crop");
            if (i + 2 >= argc)
            {
                cropZ = cropY = cropX = textToInteger(argv[i+1]);
            }
            else if (i + 3 >= argc)
            {
                cropZ = cropY = textToInteger(argv[i+2]);
                cropX = textToInteger(argv[i+1]);
            }
            else if (i + 4 >= argc)
            {
                cropZ = textToInteger(argv[i+3]);
                cropY = textToInteger(argv[i+2]);
                cropX = textToInteger(argv[i+1]);
            }
            else REPORT_ERROR(1, "Not enough parameters after -crop");
        }

        else
            REPORT_ERROR(1, "Unknown windowing type");
    }

    void show()
    {
        Prog_parameters::show();
        if (size_mode)
            std::cout << "New size: (XxYxZ)=" << sizeX << "x" << sizeY << "x"
                      << sizeZ << std::endl;
        else if (crop_mode)
            std::cout << "Crop: (XxYxZ)=" << cropX << "x" << cropY << "x"
                      << cropZ << std::endl;
        else
            std::cout << "New window: from (z0,y0,x0)=(" << z0 << ","
                      << y0 << "," << x0 << ") to (zF,yF,xF)=(" << zF << "," << yF
                      << "," << xF << ")\n"
                      << "Physical: " << physical_coords << std::endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-physical]               : Use physical instead of logical\n"
                  << "                             coordinates\n"
                  << "  [-pad_value <val>]        : value used for padding\n"
                  << "  [-corner_pad_value]       : use the value of the upper\n"
                  << "                              left corner for padding\n"
                  << "  [-average_pad_value]      : use the image average for padding\n"
                  << "  [-r0 <x0> <y0> [<z0>]     : Window using window corners\n"
                  << "   -rF <xF> <yF> [<zF>]]    : by default indexes are logical\n"
                  << "  [-size <sizeX> [<sizeY>] [<sizeZ>]]: Window to a new size\n"
                  << "                            : if only one is given, the other two\n"
                  << "                              are supposed to be the same\n"
                  << "  [-crop <sizeX> [<sizeY>] [<sizeZ>]]: Crop this amount of pixels in each direction\n"
                  << "                            : if only one is given, the other two\n"
                  << "                              are supposed to be the same\n"
                  ;
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Window_parameters *eprm = (Window_parameters *) prm;
    if (eprm->average_pad)
        eprm->init_value = (img()).computeAvg();
    else if (eprm->corner_pad)
        eprm->init_value = DIRECT_MAT_ELEM(img(), 0, 0);
    if (eprm->crop_mode)
    {
        int xl=eprm->cropX/2;
        int xr=eprm->cropX-xl;
        int yl=eprm->cropY/2;
        int yr=eprm->cropY-yl;
        img().window(STARTINGY(img())+yl,STARTINGX(img())+xl,
                     FINISHINGY(img())-yr,FINISHINGX(img())-xr);
    }
    else if (!eprm->physical_coords)
        img().window(eprm->y0, eprm->x0, eprm->yF, eprm->xF,
                     eprm->init_value);
    else img().window(STARTINGY(img()) + eprm->y0,
                          STARTINGX(img()) + eprm->x0, STARTINGY(img()) + eprm->yF,
                          STARTINGX(img()) + eprm->xF, eprm->init_value);
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Window_parameters *eprm = (Window_parameters *) prm;
    if (eprm->average_pad)
        eprm->init_value = (vol()).computeAvg();
    else if (eprm->corner_pad)
        eprm->init_value = DIRECT_VOL_ELEM(vol(), 0, 0, 0);
    if (eprm->crop_mode)
    {
        int xl=eprm->cropX/2;
        int xr=eprm->cropX-xl;
        int yl=eprm->cropY/2;
        int yr=eprm->cropY-yl;
        int zl=eprm->cropZ/2;
        int zr=eprm->cropZ-zl;
        vol().window(STARTINGZ(vol())+zl,STARTINGY(vol())+yl,
                     STARTINGX(vol())+xl, FINISHINGZ(vol())-zr,FINISHINGY(vol())-yr,
                     FINISHINGX(vol())-xr);
    }
    else if (!eprm->physical_coords)
        vol().window(eprm->z0, eprm->y0, eprm->x0, eprm->zF, eprm->yF,
                     eprm->xF, eprm->init_value);
    else vol().window(STARTINGZ(vol()) + eprm->z0,
                          STARTINGY(vol()) + eprm->y0,
                          STARTINGX(vol()) + eprm->x0, STARTINGZ(vol()) + eprm->zF,
                          STARTINGY(vol()) + eprm->yF, STARTINGX(vol()) + eprm->xF,
                          eprm->init_value);
    return true;
}

int main(int argc, char **argv)
{
    Window_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
