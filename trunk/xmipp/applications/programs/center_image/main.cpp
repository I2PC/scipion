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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/progs.h>
#include <data/args.h>
#include <data/filters.h>

class Center_parameters: public Prog_parameters
{
public:
    int Niter;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        Niter = textToInteger(getParameter(argc, argv, "-iter","10"));
    }

    void show()
    {
        Prog_parameters::show();
        std::cout << "Iterations = " << Niter << std::endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-iter <N=10>]           : Number of iterations\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Center_parameters *eprm = (Center_parameters *) prm;
    centerImage(img(),eprm->Niter);
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    std::cerr << "This function is not implemented for volumes\n";
    return false;
}

int main(int argc, char **argv)
{
    Center_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
