/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

#include <reconstruction/micrograph_phase_flipping.h>

void Usage(char *argv[]);

int main(int argc, char **argv)
 {

    Prog_micrograph_phase_flipping prm;

    // Get command line parameters ------------------------------------------
    try
    {
        prm.fn_in   = getParameter(argc, argv, "-i");
        prm.fn_out  = getParameter(argc, argv, "-o");
        prm.fnt_ctf = getParameter(argc, argv, "-ctf");
        prm.reversed      = checkParameter(argc, argv, "-reverse_endian");
    }
    catch (Xmipp_error XE)
    {
        Usage(argv);
        exit(1);
    }
    // Main program ---------------------------------------------------------
    try
    {
        prm.show();
        prm.run();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    exit(0);

}

void Usage(char *argv[])
{
    std::cout << "Purpose:\n"
    << "    flip micrograph fase\n"
    << "Usage: "<< argv[0] 
    << "\n"
    << " -i <input_micrograph>            : Either 8, 16 or 32 bits\n"
    << " -o <output_micrograph>           : Spider Format\n"
    << " -ctf <ctf_param_file>            : CTf param file\n"
    << " -reverse_endian                  : assumme no native endianess\n"
    << std::endl;
}
