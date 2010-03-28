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
#include <data/geometry.h>

class Flip_parameters: public Prog_parameters
{
public:
    bool flipX;
    bool flipY;
    bool flipZ;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        flipX = checkParameter(argc, argv, "-flipX");
        flipY = checkParameter(argc, argv, "-flipY");
        flipZ = checkParameter(argc, argv, "-flipZ");
    }

    void show()
    {
        Prog_parameters::show();
        std::cout << "FlipX = " << flipX << std::endl
        << "FlipY = " << flipY << std::endl
        << "FlipZ = " << flipZ << std::endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-flipX]                  : Flip along X\n"
        << "  [-flipY]                  : Flip along Y\n"
        << "  [-flipZ]                  : Flip along Z\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Flip_parameters *eprm = (Flip_parameters *) prm;
    if (eprm->flipX) img().selfReverseX();
    if (eprm->flipY) img().selfReverseY();
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Flip_parameters *eprm = (Flip_parameters *) prm;
    if (eprm->flipX) vol().selfReverseX();
    if (eprm->flipY) vol().selfReverseY();
    if (eprm->flipZ) vol().selfReverseZ();
    return true;
}

int main(int argc, char **argv)
{
    Flip_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

