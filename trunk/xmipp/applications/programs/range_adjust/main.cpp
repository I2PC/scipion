/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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

#include <data/range_adjust.h>

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Prog_Range_adjust_Parameters *eprm = (Prog_Range_adjust_Parameters *) prm;
    eprm->apply(img());
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Prog_Range_adjust_Parameters *eprm = (Prog_Range_adjust_Parameters *) prm;
    eprm->apply(vol());
    return true;
}

int main(int argc, char **argv)
{
    Prog_Range_adjust_Parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
