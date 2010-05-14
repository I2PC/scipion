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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/metadata.h>
#include <data/args.h>

void Usage();

int main(int argc, char **argv)
{
    FileName fn_in, fn_out, fn_root;
    std::string sortLabel;
    MetaData  SFin, SFout, SFtmp, SFtmp2;
    bool     dont_randomize;
    bool     dont_sort;
    int N;

    try
    {
        fn_in = getParameter(argc, argv, "-i");
        N = textToInteger(getParameter(argc, argv, "-n", "2"));
        fn_root = getParameter(argc, argv, "-o", "");
        dont_randomize = checkParameter(argc, argv, "-dont_randomize");
        dont_sort      = checkParameter(argc, argv, "-dont_sort");
        if(!dont_sort)
            sortLabel = getParameter(argc, argv, "-l","image");
        if (fn_root == "")
            fn_root = fn_in.without_extension();
        SFin.read(fn_in);
    }
    catch (Xmipp_error)
    {
        Usage();
        exit(1);
    }

    try
    {
        if (!dont_randomize)
        {
        	std::cerr << "randomizing" << std::endl;
            SFtmp.randomize(SFin);
            SFtmp.write("test.xmd");
        }
        else
            SFtmp = SFin;

        int Num_images = (int)SFtmp.size();
        int Num_groups = N;
        if (Num_groups > Num_images)
            Num_groups = Num_images;
        int imagesGroup = ceil((double)Num_images / Num_groups);


        for (int i = 0;i < Num_groups;i++)
        {
            SFout.clear();
            SFout.fillWithNextNObjects(SFtmp,i*imagesGroup,imagesGroup);

            if (!dont_sort)
            {
            	std::cerr << "Sorting metada set: "<< i  <<std::endl;
                SFtmp2.sort(SFout,MDL::str2Label(sortLabel));
                SFout  = SFtmp2;
            }
            fn_out = fn_root;
            if (N!=1)
            {
                std::string num = "_" + integerToString(i + 1);
                fn_out +=  num;
            }
            fn_out += ".xmd";
            SFout.write(fn_out);
        }

    }
    catch (Xmipp_error)
    {
        std::cerr << "ERROR, exiting..." << std::endl;
        exit(1);
    }

}

void Usage()
{
    std::cout << "Usage: split_selfile [options]\n"
    << "    -i <selfile>            : Input MetaData File\n"
    << "  [ -n <int=2> ]            : Number of output MetaDatas\n"
    << "  [ -o <rootname=metadata> ]: Rootname for output MetaDatas\n"
    << "                              output will be: rootname_<n>.xpd\n"
    << "  [ -dont_randomize ]       : Do not generate random groups\n"
    << "  [ -dont_sort ]            : Do not sort the output MetaData\n"
    << "  [ -l <image>]     : sort using label sortLabel, default image\n"
    ;
}
