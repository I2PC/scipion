/***************************************************************************
 *
 * Authors:     Carlos Oscar S�nchez Sorzano
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
        fn1    = get_param(argc, argv, "-i1");
        fn2    = get_param(argc, argv, "-i2");
        fn_out = get_param(argc, argv, "-o");
    }
    catch (Xmipp_error XE)
    {
        cerr << XE << endl;
        cerr << "Usage: appenddocfile\n"
        << "   -i1 <docfile1>    : Input file 1\n"
        << "   -i2 <docfile2>    : Input file 2\n"
        << "   -o  <docfile1>    : Concatenated file: 1+2\n";
        return 1;
    }

    try
    {
        DocFile DF(fn1);
        DF.append(fn2);
        DF.write(fn_out);
    }
    catch (Xmipp_error XE)
    {
        cerr << XE << endl;
        return 2;
    }
    return 0;
}
