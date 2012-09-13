/***************************************************************************
 *
 * Authors:     Debora Gil
                Roberto Marabini
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

#include <data/args.h>
#include <reconstruction/crystal_unbend.h>

void Usage(char *argv[]);

int main(int argc, char **argv)
{

    ImUmbend prm;


    // Get command line parameters ------------------------------------------
    try
    {
        prm.inImfile  = getParameter(argc, argv, "-i");
        prm.outImfile = getParameter(argc, argv, "-o");
        prm.FN_Correlation = getParameter(argc, argv, "-cor");
        prm.cc_peak_factor =  textToFloat(getParameter(argc, argv, "-cc_peak_factor", "0.0"));
        prm.InterpModel = getParameter(argc, argv, "-interp_model", "Linear");

    }
    catch (Xmipp_error XE)
    {
        Usage(argv);
    }

    // Main program ---------------------------------------------------------
    try
    {
        ROUT_Umbend(prm);
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }
    exit(0);
}

void Usage(char *argv[])
{
    cout << "Purpose:\n"
    << "    Computes 2D crystal distorsion AND unbends the crystal\n"
    << "    The input is a cor (MRC) file\n"
    << "Usage:" << argv[0] << " -i InImage -o OutImage -cor filename"
    << " -cc_peak_factor cc_peak_factor -interp_model InterpModel" << endl << endl
    << "\t-i               :  Input Image" << endl
    << "\t-o               :  Output image" << endl
    << "\t-cor             : Correlation file" << endl
    << "\t-cc_peak_factor  : crosscorrelation thershold (0-1)" << endl
    << "\t-interp_model    : interpolation scheme (Linear or Bessel) "
    << "used for computation of crystal distorsion" << endl
    ;
    exit(1);

}
