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

#include <data/selfile.h>
#include <data/args.h>

int main(int argc, char **argv)
{
    int mode;
    FileName fn1, fn2, fn_out, fn_meth;
    try
    {
        fn1   = getParameter(argc, argv, "-i1");
        fn2   = getParameter(argc, argv, "-i2");
        fn_out = getParameter(argc, argv, "-o");
	fn_meth = getParameter(argc, argv, "-method","");
	if (fn_meth == "both")
	    mode = 0;
	else if (fn_meth == "file1")
	    mode = 1;
	else if (fn_meth == "file2")
	    mode = 2;
	else 
	    mode = -1;

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        cout << "Usage: compare_selfiles\n"
	     << "   -i1 <selfile1>      : First  selfile to compare\n"
	     << "   -i2 <selfile2>      : Second selfile to compare\n"
	     << "   -o  <selfile_out>   : Output selfile to compare\n"
	     << "  [-method both]       : Output selfile with images that occur in both files\n"
	     << "  [-method file1]      : Output selfile with images that only occur in file 1\n"
	     << "  [-method file2]      : Output selfile with images that only occur in file 2\n"
        ;
    }

    try
    {
        SelFile SF1(fn1), SF2(fn2);
        SelFile SFout = compare(SF1, SF2, mode);
        SFout.write(fn_out);
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }

    return 0;
}
