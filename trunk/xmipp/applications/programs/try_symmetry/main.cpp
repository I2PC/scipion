/***************************************************************************
 *
 * Authors:     Slavica Jonic (slavica.jonic@impmc.jussieu.fr)
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

class Markhan_parameters: public Prog_parameters
{
public:
    int symmetry;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        symmetry = textToInteger(getParameter(argc, argv, "-symorder"));
    }

    void show()
    {
        Prog_parameters::show();
        std::cout << "Symmetry order = " << symmetry << std::endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-symorder <n>]           : Symmetry order\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Markhan_parameters *eprm = (Markhan_parameters *) prm;
    Matrix2D<double> aux = img();
    for (int i = 1; i < eprm->symmetry; i++)
    {
        Matrix2D<double> rotatedImg;
        img().rotateBSpline(3, 360.0 / eprm->symmetry * i, rotatedImg);
        aux += rotatedImg;
    }
    aux /= eprm->symmetry;
    img() = aux;
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    REPORT_ERROR(1, "xmipp_markhan: This program is not intended for volumes. "
                 "Please use xmipp_symmetrize.");
    return true;
}

int main(int argc, char **argv)
{
    Markhan_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
