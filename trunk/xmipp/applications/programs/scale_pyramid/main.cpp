/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
        Prog_parameters::show();
        cout << "Operation: ";
        switch (operation)
        {
        case Expand:
            cout << "Expand\n";
            break;
        case Reduce:
            cout << "Reduce\n";
            break;
        case None  :
            cout << "None  \n";
            break;
        }
        cout << "Levels: " << levels << endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        cerr << "  -expand | -reduce         : Expand or reduce the image\n";
        cerr << " [-levels=<l=1>]            : Expansion/reduction factor\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Pyramid_parameters *eprm = (Pyramid_parameters *) prm;
    float Xoff, Yoff;
    img.get_originOffsets(Xoff, Yoff);
    Matrix2D<double> result;
    float scale_factor = (float)(pow(2.0, eprm->levels));
    switch (eprm->operation)
    {
    case Pyramid_parameters::Expand:
        img().pyramidExpand(result, eprm->levels);
        img.set_originOffsets(Xoff*scale_factor, Yoff*scale_factor);
        break;
    case Pyramid_parameters::Reduce:
        img().pyramidReduce(result, eprm->levels);
        img.set_originOffsets(Xoff / scale_factor, Yoff / scale_factor);
        break;
    }
    img() = result;

    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Pyramid_parameters *eprm = (Pyramid_parameters *) prm;
    Matrix3D<double> result;
    switch (eprm->operation)
    {
    case Pyramid_parameters::Expand:
        vol().pyramidExpand(result, eprm->levels);
        break;
    case Pyramid_parameters::Reduce:
        vol().pyramidReduce(result, eprm->levels);
        break;
    }
    vol() = result;
    return true;
}

int main(int argc, char **argv)
{
    Pyramid_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
