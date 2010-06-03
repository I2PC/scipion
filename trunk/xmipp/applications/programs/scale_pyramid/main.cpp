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

class Pyramid_parameters: public Prog_parameters
{
public:
    enum Toperation {Expand, Reduce, None};
    Toperation operation;
    int levels;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        levels = textToInteger(getParameter(argc, argv, "-levels", "1"));
        if (checkParameter(argc, argv, "-expand")) operation = Expand;
        else if (checkParameter(argc, argv, "-reduce")) operation = Reduce;
        else                                       operation = None;
    }

    void show()
    {
        if (quiet) return;
        Prog_parameters::show();
        std::cout << "Operation: ";
        switch (operation)
        {
        case Expand:
            std::cout << "Expand\n";
            break;
        case Reduce:
            std::cout << "Reduce\n";
            break;
        case None  :
            std::cout << "None  \n";
            break;
        }
        std::cout << "Levels: " << levels << std::endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  -expand | -reduce         : Expand or reduce the image\n";
        std::cerr << " [-levels=<l=1>]            : Expansion/reduction factor\n";
    }
};

bool process_img(Image<double> &img, const Prog_parameters *prm)
{
    Pyramid_parameters *eprm = (Pyramid_parameters *) prm;
    float Xoff, Yoff, Zoff;
    Xoff=img.Xoff(); Yoff=img.Yoff(); Zoff=img.Zoff();
    MultidimArray<double> result;
    float scale_factor = (float)(pow(2.0, eprm->levels));
    switch (eprm->operation)
    {
    case Pyramid_parameters::Expand:
        pyramidExpand(3,result,img(),eprm->levels);
        img.setXoff(Xoff*scale_factor);
        img.setYoff(Yoff*scale_factor);
        img.setZoff(Zoff*scale_factor);
        break;
    case Pyramid_parameters::Reduce:
        pyramidReduce(3,result,img(),eprm->levels);
        img.setXoff(Xoff*scale_factor);
        img.setYoff(Yoff*scale_factor);
        img.setZoff(Zoff*scale_factor);
        break;
    }
    img() = result;

    return true;
}

int main(int argc, char **argv)
{
    Pyramid_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img);
    return 0;
}
