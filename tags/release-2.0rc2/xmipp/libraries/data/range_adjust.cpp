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

#include "range_adjust.h"
#include "args.h"

/* Read parameters --------------------------------------------------------- */
void Prog_Range_adjust_Parameters::read(int argc, char **argv)
{
    min_val = textToFloat(getParameter(argc, argv, "-min"));
    max_val = textToFloat(getParameter(argc, argv, "-max"));
    sigma = textToFloat(getParameter(argc, argv, "-noise", "0"));
    randomize_random_generator();
}

/* Usage ------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::usage()
{
    cerr << "   -min <min_val>       : Output minimum value\n"
    << "   -max <max_val>       : Output maximum value\n"
    << "  [-noise <sigma %>]    : Variation of the minimum value\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::show()
{
    cout << "Min:   " << min_val << endl
    << "Max:   " << max_val << endl
    << "Noise: " << sigma   << endl;
}

/* Apply ------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::apply(Matrix2D<double> &I)
{
    double amin = rnd_gaus(0, sigma);
    double amax = rnd_gaus(0, sigma);
    double minval = min_val + amin;
    double maxval = max_val + amax;
    I.range_adjust(minval, maxval);
}

void Prog_Range_adjust_Parameters::apply(Matrix3D<double> &V)
{
    double amin = rnd_gaus(0, sigma);
    double amax = rnd_gaus(0, sigma);
    double minval = min_val + amin;
    double maxval = max_val + amax;
    V.range_adjust(minval, maxval);
}
