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

#include "range_adjust.h"
#include "args.h"

/* Read parameters --------------------------------------------------------- */
void Prog_Range_adjust_Parameters::read(int argc, char **argv)
{
    Prog_parameters::read(argc,argv);
    min_val = textToFloat(getParameter(argc, argv, "-min"));
    max_val = textToFloat(getParameter(argc, argv, "-max"));
    sigma = textToFloat(getParameter(argc, argv, "-noise", "0"));
    randomize_random_generator();
    mask_prm.read(argc, argv);
}

/* Usage ------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::usage()
{
    Prog_parameters::usage();
    std::cerr << "   -min <min_val>       : Output minimum value\n"
              << "   -max <max_val>       : Output maximum value\n"
              << "  [-noise <sigma %>]    : Variation of the minimum value\n"
              ;
    mask_prm.usage();
}

/* Show -------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::show()
{
    Prog_parameters::show();
    std::cout << "Min:   " << min_val << std::endl
              << "Max:   " << max_val << std::endl
              << "Noise: " << sigma   << std::endl;
}

/* Apply ------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::apply(Matrix2D<double> &I)
{
    double amin = rnd_gaus(0, sigma);
    double amax = rnd_gaus(0, sigma);
    double minval = min_val + amin;
    double maxval = max_val + amax;
    if (mask_prm.type == NO_MASK)
        I.rangeAdjust(minval, maxval);
    else
    {
        mask_prm.generate_2Dmask(I);
        I.rangeAdjust(minval, maxval,mask_prm.get_binary_mask2D());
    }
}

void Prog_Range_adjust_Parameters::apply(Matrix3D<double> &V)
{
    double amin = rnd_gaus(0, sigma);
    double amax = rnd_gaus(0, sigma);
    double minval = min_val + amin;
    double maxval = max_val + amax;
    if (mask_prm.type == NO_MASK)
        V.rangeAdjust(minval, maxval);
    else
    {
        mask_prm.generate_3Dmask(V);
        V.rangeAdjust(minval, maxval,mask_prm.get_binary_mask3D());
    }
}
