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

#include <data/selfile.h>
#include <data/args.h>

void Usage();

int main(int argc, char **argv)
{
    FileName fn_in, fn_root;
    int Nsubsets;
    int subsetSize;
    bool bootstrap;

    try
    {
        fn_in = getParameter(argc, argv, "-i");
        subsetSize = textToInteger(getParameter(argc, argv, "-subsetSize","0"));
        Nsubsets = textToInteger(getParameter(argc, argv, "-numberOfSubsets", "1"));
        bootstrap = checkParameter(argc, argv, "-bootstrap");
        fn_root = getParameter(argc, argv, "-o", "");
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE << std::endl;
        Usage();
        exit(1);
    }

    try
    {
        SelFile  SFin, SFout;
        SFin.read(fn_in);
        if (fn_root == "") fn_root = fn_in.without_extension();
        if (subsetSize == 0) subsetSize=SFin.ImgNo();

        for (int i=0; i<Nsubsets; i++)
        {
            SFout=SFin.randomSubset(subsetSize,bootstrap);
            SFout.write(fn_root+"_"+integerToString(i+1,6)+".sel");
        }

    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE << std::endl;
        exit(1);
    }
}

void Usage()
{
    std::cout << "Usage: split_selfile [options]\n"
              << "    -i <selfile>            : Input selfile\n"
              << "  [ -o <rootname=selfile> ] : Rootname for output selfiles\n"
              << "                              output will be: rootname_<n>.sel\n"
              << "  [ -subsetSize <int>]      : Number of images in the subset\n"
              << "                              If not given, the same as the input file\n"
              << "  [ -numberOfSubsets <int> ]: Number of subsets of this kind\n"
              << "  [ -bootstrap ]            : Random subsets with replacement\n"
              ;
}
