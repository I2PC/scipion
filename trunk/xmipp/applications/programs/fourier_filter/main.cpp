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
#include <reconstruction/fourier_filter.h>

class FourierFilter_parameters: public Prog_parameters
{
public:
    FourierMask fmask;
    bool        first;
    int         dim;
public:
    void read(int argc, char **argv)
    {
        fmask.read(argc, argv);
        Prog_parameters::read(argc, argv);
        first = true;
    }

    void show()
    {
        Prog_parameters::show();
        fmask.show();
    }

    void usage()
    {
        Prog_parameters::usage();
        fmask.usage();
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    FourierFilter_parameters *eprm = (FourierFilter_parameters *) prm;
    if (eprm->first)
    {
        eprm->fmask.generate_mask(img());
        eprm->first = false;
    }
    eprm->fmask.apply_mask_Space(img());
    eprm->dim = 2;
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    FourierFilter_parameters *eprm = (FourierFilter_parameters *) prm;
    if (eprm->first)
    {
        eprm->fmask.generate_mask(vol());
        eprm->first = false;
    }
    eprm->fmask.apply_mask_Space(vol());
    eprm->dim = 3;
    return true;
}

int main(int argc, char **argv)
{
    FourierFilter_parameters prm;
    SF_main(argc, argv, &prm, (void *)&process_img, (void *)&process_vol);
}
