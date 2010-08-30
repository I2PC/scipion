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

#include <data/image.h>
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
        min_size  = textToFloat(getParameter(argc, argv, "-min_size", "0"));
        if (fn_root == "")
            fn_root = fn_in.get_root();
    }
    catch (XmippError XE)
    {
        std::cout << XE;
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
        Image<double> I, label;
        I.read(fn_in);
        int object_no;
        if (ZSIZE(I())==1)
        	object_no=label_image2D(I(), label());
        else
        	object_no=label_image3D(I(), label());
        for (int o = 0; o <= object_no; o++)
        {
            I() = label();
            MultidimArray<double> &Im=MULTIDIM_ARRAY(I);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Im)
            {
                DIRECT_MULTIDIM_ELEM(Im,n) = (DIRECT_MULTIDIM_ELEM(Im,n) == o);
                if (invert)
                    DIRECT_MULTIDIM_ELEM(Im,n) = 1 - DIRECT_MULTIDIM_ELEM(Im,n);
            }
            number_elements = I().sum();
            if (number_elements > min_size)
            {
                fn_out.compose(fn_root, o, fn_ext);
                I.write(fn_out);
            }

            if (ZSIZE(I())==1)
                std::cout << "Image number " << o << " contains " << number_elements
                << " pixels set to 1\n";
            else
                std::cout << "Volume number " << o << " contains " << number_elements
                << " voxels set to 1\n";
        }
    }
    catch (XmippError XE)
    {
        std::cout << XE;
    }
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cerr << "Usage: separate_objects\n"
    << "   -i <fn_in>                      : Input image or volume\n"
    << "  [-o <fn_root>]                   : Root filename for output\n"
    << "                                     By default, the input name\n"
    << "  [-invert]                        : Produce inverse masks\n"
    << "  [-min_size <size=0>]             : Save if size is greater than this\n"
    ;
}
