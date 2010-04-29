/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Pedro A. de Alarcón (pedro@cnb.csic.es)
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
#include <data/morphology.h>

class Morphology_parameters: public Prog_parameters
{
public:
#define DILATION     1
#define EROSION      2
#define OPENING      3
#define CLOSING      4
#define SHARPENING   5
    int operation;

    int size;
    int count;
    int neig;
    double width;
    double strength;
public:
    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        if (checkParameter(argc, argv, "-dil"))      operation = DILATION;
        if (checkParameter(argc, argv, "-ero"))      operation = EROSION;
        if (checkParameter(argc, argv, "-clo"))      operation = CLOSING;
        if (checkParameter(argc, argv, "-ope"))      operation = OPENING;
        if (checkParameter(argc, argv, "-sharp"))
        {
            operation = SHARPENING;
            int i = paremeterPosition(argc, argv, "-sharp");
            if (i+2>=argc)
                REPORT_ERROR(1,"Not enough parameters after -sharp");
            width=textToFloat(argv[i+1]);
            strength=textToFloat(argv[i+2]);
        }

        size = textToInteger(getParameter(argc, argv, "-size", "1"));
        neig = textToInteger(getParameter(argc, argv, "-neig", "-1"));
        count = textToInteger(getParameter(argc, argv, "-count", "0"));
    }

    void show()
    {
        Prog_parameters::show();
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
        }
        if (operation!=SHARPENING)
            std::cout << "Size=" << size << std::endl
                      << "Neighbourhood=" << neig << std::endl
                      << "Count=" << count << std::endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-dil]             : Apply dilation\n"
                  << "  [-ero]             : Apply erosion\n"
                  << "  [-clo]             : Apply closing\n"
                  << "  [-ope]             : Apply opening\n"
                  << "  [-neig <n=8 | 18>] : Neighborhood considered \n"
                  << "                       (2D:4,8 3D:6,18,26)\n"
                  << "  [-size <s=1>]      : Size of the Strutural element\n"
                  << "  [-count <c=0>]     : Minimum required neighbors with \n"
                  << "                       distinct value\n"
                  << "  [-sharp <w> <s>]   : Sharpening with width (suggested 1 or 2)\n"
                  << "                       and strength (suggested 0.1-1.0)\n"
                  ;
    }
};


bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Morphology_parameters *eprm = (Morphology_parameters *) prm;
    if (eprm->neig == -1) eprm->neig = 8;
    ImageXmipp retval;
    retval() = img();

    if (eprm->operation!=SHARPENING)
        std::cout << "Initially the image has " << img().sum()
                  << " pixels set to 1\n";
    switch (eprm->operation)
    {
    case DILATION:
        dilate2D(img(), retval(), eprm->neig, eprm->count, eprm->size);
        break;
    case EROSION:
        erode2D(img(), retval(), eprm->neig, eprm->count, eprm->size);
        break;
    case OPENING:
        opening2D(img(), retval(), eprm->neig, eprm->count, eprm->size);
        break;
    case CLOSING:
        closing2D(img(), retval(), eprm->neig, eprm->count, eprm->size);
        break;
    case SHARPENING:
        REPORT_ERROR(1,"Sharpening has not been implemented for images");
    }

    img() = retval();
    std::cout << "Finally the image has " << img().sum() << " pixels set to 1\n";
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Morphology_parameters *eprm = (Morphology_parameters *) prm;
    if (eprm->neig == -1) eprm->neig = 18;
    VolumeXmipp retval;
    retval() = vol();

    if (eprm->operation!=SHARPENING)
        std::cout << "Initially the volume has " << vol().sum()
                  << " voxels set to 1\n";
    switch (eprm->operation)
    {
    case DILATION:
        dilate3D(vol(), retval(), eprm->neig, eprm->count, eprm->size);
        break;
    case EROSION:
        erode3D(vol(), retval(), eprm->neig, eprm->count, eprm->size);
        break;
    case OPENING:
        opening3D(vol(), retval(), eprm->neig, eprm->count, eprm->size);
        break;
    case CLOSING:
        closing3D(vol(), retval(), eprm->neig, eprm->count, eprm->size);
        break;
    case SHARPENING:
        sharpening(vol(),eprm->width, eprm->strength, retval());
    }

    vol() = retval();
    if (eprm->operation!=SHARPENING)
        std::cout << "Finally the volume has " << vol().sum()
                  << " voxels set to 1\n";
    return true;
}

int main(int argc, char **argv)
{
    Morphology_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
