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

#include <data/image.h>
#include <data/volume.h>
#include <data/args.h>
#include <data/filters.h>

void Usage();

/* ------------------------------------------------------------------------- */
/* Main                                                                      */
/* ------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    FileName fn_in, fn_root;
    bool invert;
    double min_size;


    // Get input parameters .................................................
    try
    {
        fn_in   = getParameter(argc, argv, "-i");
        fn_root = getParameter(argc, argv, "-o", "");
        invert  = checkParameter(argc, argv, "-invert");
        min_size  = AtoF(getParameter(argc, argv, "-min_size", "0"));
        if (fn_root == "") fn_root = fn_in.get_root();
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage();
        exit(0);
    }

    // Process ..............................................................
    try
    {
        double number_elements;
        int N = 0;
        FileName fn_out;
        FileName fn_ext = fn_in.get_extension();
        if (Is_ImageXmipp(fn_in))
        {
            ImageXmipp I(fn_in), label;
            int object_no = label_image(I(), label());
            for (int o = 0; o <= object_no; o++)
            {
                I() = label();
                FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
                {
                    I(i, j) = I(i, j) == o;
                    if (invert) I(i, j) = 1 - I(i, j);
                }
                number_elements = I().sum();
                if (number_elements > min_size)
                {
                    fn_out.compose(fn_root, o, fn_ext);
                    I.write(fn_out);
                }

                cout << "Image number " << o << " contains " << number_elements
                << " pixels set to 1\n";
            }
        }
        else if (Is_VolumeXmipp(fn_in))
        {
            VolumeXmipp V(fn_in), label;
            int object_no = label_volume(V(), label());
            for (int o = 0; o <= object_no; o++)
            {
                V() = label();
                FOR_ALL_ELEMENTS_IN_MATRIX3D(V())
                {
                    V(k, i, j) = V(k, i, j) == o;
                    if (invert) V(k, i, j) = 1 - V(k, i, j);
                }
                number_elements = V().sum();
                if (number_elements > min_size)
                {
                    fn_out.compose(fn_root, o, fn_ext);
                    V.write(fn_out);
                }

                cout << "Volume number " << o << " contains " << number_elements
                << " voxels set to 1\n";
            }
        }
        else
        {
            REPORT_ERROR(1, "Separate_objects: Input file is not Spider\n");
        }
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    cerr << "Usage: separate_objects\n"
    << "   -i <fn_in>                      : Input image or volume\n"
    << "  [-o <fn_root>]                   : Root filename for output\n"
    << "                                     By default, the input name\n"
    << "  [-invert]                        : Produce inverse masks\n"
    << "  [-min_size <size=0>]             : Save if size is greater than this\n"
    ;
}
