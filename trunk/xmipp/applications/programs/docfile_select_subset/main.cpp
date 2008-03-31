/***************************************************************************
 *
 * Authors:     Sjors Scheres (scheres@cnb.csic.es)
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
#include <data/docfile.h>

int main(int argc, char **argv)
{
    FileName fn1, fn2, fn_out;
    try
    {
        fn1    = getParameter(argc, argv, "-i");
        fn2    = getParameter(argc, argv, "-sel");
        fn_out = getParameter(argc, argv, "-o");
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE << std::endl;
        std::cerr << "Usage: docfile_select_subset\n"
        << "   -i   <docfile>    : Input docfile\n"
        << "   -sel <selfile>    : Input selfile\n"
        << "   -o   <docfile>    : Output docfile with subset in selfile\n";
        return 1;
    }

    try
    {
        DocFile DFin, DFout;
	SelFile SF;
	DFin.read(fn1);
	SF.read(fn2);
	get_subset_docfile(DFin, SF, DFout);
        DFout.write(fn_out);
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE << std::endl;
        return 2;
    }
    return 0;
}
