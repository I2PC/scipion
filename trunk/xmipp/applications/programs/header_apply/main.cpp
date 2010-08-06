/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

class Headerapply_parameters: public Prog_parameters
{
public:
    bool wrap;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        wrap = !checkParameter(argc, argv, "-dont_wrap");
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-dont_wrap]              : By default, the image is wrapped\n";
    }

    void show()
    {
        Prog_parameters::show();
        if (!wrap)
            std::cout << "Do not wrap"<<std::endl;
    }
};

bool process_img(Image<double> &img, const Prog_parameters *prm)
{
    if (ZSIZE(img())!=1 || NSIZE(img())!=1)
        REPORT_ERROR(1,"This program is intended only for images");

    Headerapply_parameters *eprm = (Headerapply_parameters *) prm;
    MultidimArray<double> Maux;
    applyGeometry(BSPLINE3, Maux, img(), img.getTransformationMatrix(), IS_INV, eprm->wrap);
    img()=Maux;

    //Reset in-plane transformations of the header
    img.setShifts(0,0);
    img.setPsi(0.);
    if (img.tilt() == 0)
        img.setRot(0.);
    img.setFlip(0.);

    return true;
}

int main(int argc, char **argv)
{
    Headerapply_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img);
}
