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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/progs.h>
#include <data/args.h>
#include <data/geometry.h>
#include <data/gridding.h>

class Headerapply_parameters: public Prog_parameters
{
public:

    bool gridding;
    bool wrap;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        gridding = checkParameter(argc, argv, "-gridding");
        wrap = !checkParameter(argc, argv, "-dont_wrap");
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-dont_wrap]              : By default, the image is wrapped\n"
	     << "  [-gridding]               : Use reverse gridding for interpolation\n";
    }

    void show()
    {
        Prog_parameters::show();
        if (!wrap)
            std::cout << "Do not wrap"<<std::endl;
        if (gridding)
            std::cout << "Use reverse gridding interpolation"<<std::endl;
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Headerapply_parameters *eprm = (Headerapply_parameters *) prm;
    Matrix2D<double> Maux;
    if (eprm->gridding)
    {
	KaiserBessel kb;
	produceReverseGriddingMatrix2D(img(),Maux,kb);
	applyGeometryReverseGridding(img(), img.get_transformation_matrix(), Maux, kb, IS_INV, eprm->wrap);
    }
    else
    {
	applyGeometryBSpline(Maux, img.get_transformation_matrix(), img(), 3, IS_INV, eprm->wrap);
	img()=Maux;
    }

    //Reset in-plane transformations of the header
    img.set_Xoff(0.);
    img.set_Yoff(0.);
    img.set_psi(0.);
    if (img.tilt() == 0) img.set_rot(0.);
    img.set_flip(0.);

    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    std::cerr << "Error: applygeo does not work with volumes\n";
    exit(0);
    return false;
}

int main(int argc, char **argv)
{
    Headerapply_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
