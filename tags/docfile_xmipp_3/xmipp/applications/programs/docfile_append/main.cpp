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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/docfile.h>

int main(int argc, char **argv)
{
    bool do_sel;
    FileName fn1, fn2, fn_out, remove_multiple;
    std::string s;
    try
    {
	if (checkParameter(argc, argv, "-sel"))
	{
	    fn1 = getParameter(argc, argv, "-sel");
	    do_sel = true;
	}
	else
	{
	    do_sel = false;
	    fn1    = getParameter(argc, argv, "-i1");
	    fn2    = getParameter(argc, argv, "-i2");
	}
        fn_out = getParameter(argc, argv, "-o");
	remove_multiple = getParameter(argc, argv, "-remove_multiple","");

    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE << std::endl;
        std::cerr << "Usage: docfile_append\n"
        << "   -i1 <docfile1>    : Input file 1\n"
        << "   -i2 <docfile2>    : Input file 2\n"
	<< " OR \n"
	<< " -sel <selfile>      : Input selfile of docfiles \n"
	<< "\n" 
        << "   -o  <docfile1>    : Output concatenated file\n"
        << "  [-remove_multiple <string=\"\">] : remove multiple instances of comment lines containing this string\n"; 
	return 1;
    }

    try
    {
	DocFile DF;
	if (do_sel)
	{
	    SelFile SF;
	    SF.read(fn1);
	    SF.go_beginning();
	    DF.read(SF.NextImg());
	    while (!SF.eof())
	    {
                FileName fn_img=SF.NextImg();
                if (fn_img=="") break;
		DF.append(fn_img);
	    }
	}
	else
	{
	    DF.read(fn1);
	    DF.append(fn2);
	}
	if (remove_multiple!="")
	{
	    DF.remove_multiple_strings(remove_multiple);
        }
	DF.write(fn_out);
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE << std::endl;
        return 2;
    }
    return 0;
}
