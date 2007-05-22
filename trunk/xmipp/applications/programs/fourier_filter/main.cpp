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
#include <reconstruction/fourier_filter.h>

class FourierFilter_parameters: public Prog_parameters
{
public:
    FourierMask fmask;
    FileName    fn_mask;
    FileName    fn_amplitude;
    bool        do_not_center;

    bool        first;
    int         dim;
public:
    void read(int argc, char **argv)
    {
        fmask.read(argc, argv);
        Prog_parameters::read(argc, argv);
        fn_mask = get_param(argc, argv, "-save_mask", "");
        fn_amplitude = get_param(argc, argv, "-save_amplitude", "");
        do_not_center = check_param(argc, argv, "-do_not_center");

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
        cerr << "  [-save_mask <fn_mask>]                : Save applied filter\n"
        << "  [-save_amplitude <fn_ampl>]           : Save amplitude of filter\n"
        << "  [-do_not_center]                      : For the amplitude file\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    FourierFilter_parameters *eprm = (FourierFilter_parameters *) prm;
    matrix2D< complex<double> > fft;
    FourierTransform(img(), fft);
    if (eprm->first)
    {
        eprm->fmask.generate_mask(fft);
        eprm->first = false;
    }
    eprm->fmask.apply_mask_Fourier(fft);
    InverseFourierTransform(fft, img());
    eprm->dim = 2;
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    FourierFilter_parameters *eprm = (FourierFilter_parameters *) prm;
    matrix3D< complex<double> > fft;
    FourierTransform(vol(), fft);
    if (eprm->first)
    {
        eprm->fmask.generate_mask(fft);
        eprm->first = false;
    }
    eprm->fmask.apply_mask_Fourier(fft);
    InverseFourierTransform(fft, vol());
    eprm->dim = 3;
    return true;
}

int main(int argc, char **argv)
{
    FourierFilter_parameters prm;
    SF_main(argc, argv, &prm, (void *)&process_img, (void *)&process_vol);
    if (prm.fn_mask != "")      prm.fmask.write_mask(prm.fn_mask, 2);
    if (prm.fn_amplitude != "") prm.fmask.write_amplitude(prm.fn_amplitude,
                prm.dim, prm.do_not_center);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM FourierFilter {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/FourierFilter/Help/fourierfilter.html";
      help="Apply a frequency filter to volumes and images";
      OPEN MENU menu_fourierfilter;
      COMMAND LINES {
 + usual: xmipp_fourierfilter
               #include "prog_line.mnu"
               [-save_mask $FN_MASK]
        [-save_amplitude $FN_AMPL]
        [-do_not_center]
               #include "fourier_line.mnu"
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        $FN_MASK {type=file; label="Save applied mask as Xmipp Fourier Image";}
        $FN_AMPL {type=file; label="Save amplitude of applied mask";}
        OPT(-do_not_center) {label="Do not center";}
        #include "fourier_vars.mnu"
      }
   }

   MENU menu_fourierfilter {
      #include "prog_menu.mnu"
      OPT($FN_MASK)
      OPT($FN_AMPL)
      OPT(-do_not_center)
      #include "fourier_menu.mnu"
   }
*/
