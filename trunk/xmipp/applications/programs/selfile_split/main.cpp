/***************************************************************************
 *
 * Authors:     Sjors Scheres (scheres@cnb.uam.es)
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

void Usage();

int main(int argc, char **argv)
{
    FileName fn_in, fn_out, fn_root;
    SelFile  SFin, SFout, SFtmp, SFtmp2;
    SelLine  line;
    bool     dont_randomize;
    int N;

    try
    {
        fn_in = getParameter(argc, argv, "-i");
        N = textToInteger(getParameter(argc, argv, "-n", "2"));
        fn_root = getParameter(argc, argv, "-o", "");
        dont_randomize = checkParameter(argc, argv, "-dont_randomize");
        if (fn_root == "") fn_root = fn_in.without_extension();
        SFin.read(fn_in);
    }
    catch (Xmipp_error)
    {
        Usage();
        exit(1);
    }

    try
    {
        if (!dont_randomize) SFtmp = SFin.randomize();
        else                 SFtmp = SFin;
        int Num_images = (int)SFtmp.ImgNo();
        int Num_groups = N;
        if (Num_groups > Num_images) Num_groups = Num_images;

        int Nsub_ = (int)Num_images / N;
        int Nres_ = Num_images % N;

        int arr_groups[Num_groups];

        int i, j;

        for (i = 0;i < (Num_groups - Nres_);i++)
        {
            arr_groups[i] = Nsub_;
        }

        for (j = i;j < Num_groups;j++)
        {
            arr_groups[j] = Nsub_ + 1;
        }

        SFtmp.go_beginning();
        for (i = 0;i < Num_groups;i++)
        {
            SFout.clear();
            SFout.reserve(arr_groups[i]);
            for (j = 0;j < arr_groups[i];j++)
            {
                SFout.insert(SFtmp.current());
                SFtmp.NextImg();
            }
            SFtmp2 = SFout.sort_by_filenames();
            SFout = SFtmp2;
            string num = "_" + ItoA(i + 1);
            fn_out = fn_root + num;
            fn_out += ".sel";
            SFout.write(fn_out);
        }

    }
    catch (Xmipp_error)
    {
        cerr << "ERROR, exiting..." << endl;
        exit(1);
    }

}

void Usage()
{
    cout << "Usage: split_selfile [options]\n"
    << "    -i <selfile>            : Input selfile\n"
    << "  [ -n <int=2> ]            : Number of output selfiles\n"
    << "  [ -o <rootname=selfile> ] : Rootname for output selfiles\n"
    << "                              output will be: rootname_<n>.sel\n"
    << "  [ -dont_randomize ]       : Do not generate random groups\n"
    ;
}
