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

#include <reconstruction/angular_continuous_assign.h>

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Prog_angular_predict_continuous_prm *eprm =
        (Prog_angular_predict_continuous_prm *) prm;
    double shiftX, shiftY, psi, rot, tilt;
    eprm->get_initial_guess(shiftX, shiftY, rot, tilt, psi);
    double corr = eprm->predict_angles(img, shiftX, shiftY, rot, tilt, psi);
    img.read(img.name());
    img.set_eulerAngles(rot, tilt, psi);
    img.set_originOffsets(shiftX, shiftY);
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    cout << "This program is not intended for volumes\n";
    return false;
}

int main(int argc, char **argv)
{
    Prog_angular_predict_continuous_prm prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
    prm.finish_processing();
}
