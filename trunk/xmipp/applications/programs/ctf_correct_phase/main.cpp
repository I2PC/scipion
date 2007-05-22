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
#include <reconstruction/ctf_correct_phase.h>

class Prog_CorrectPhase_Params: public Prog_parameters
{
public:
    CorrectPhase_Params cpprm;
public:
    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        cpprm.read(argc, argv);
        cpprm.produce_side_info();
        allow_time_bar = false;
    }

    void show()
    {
        Prog_parameters::show();
        cpprm.show();
    }

    void usage()
    {
        Prog_parameters::usage();
        cpprm.usage();
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Prog_CorrectPhase_Params *eprm = (Prog_CorrectPhase_Params *) prm;
    matrix2D< complex<double> > fft;
    FileName       fn_ctf;
    if (eprm->cpprm.multiple_CTFs)
    {
        fn_ctf = eprm->cpprm.SF_CTF.NextImg();
        eprm->cpprm.ctf.ctf.read(fn_ctf);
        eprm->cpprm.ctf.ctf.Produce_Side_Info();
        cerr << "Correcting " << img.name() << " with " << fn_ctf << endl;
    }
    FourierTransform(img(), fft);
    eprm->cpprm.correct(fft);
    InverseFourierTransform(fft, img());
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    cerr << "This process is not intended for volumes\n";
    return false;
}

int main(int argc, char **argv)
{
    Prog_CorrectPhase_Params prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
