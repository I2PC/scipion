/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2010)
 *             Fernando Fuentes
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

#include <reconstruction/enhance_contrast.h>

bool process_vol(Image<double> &vol, const Prog_parameters *prm)
{
    EnhanceContrast_parameters *eprm = (EnhanceContrast_parameters *) prm;
    eprm->enhance(vol());
    return true;
}

int main(int argc, char **argv)
{
    EnhanceContrast_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_vol);
}
