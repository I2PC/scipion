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

#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/xmipp_program.h>

class ProgSeparateObjects: public XmippProgram
{
public:
    FileName fn_in, fn_root;
    bool invert;
    double min_size;

    void readParams()
    {
        fn_in   = getParam("-i");
        fn_root = getParam("--oroot");
        invert  = checkParam("--invert");
        min_size  = getDoubleParam("--min_size");
        if (fn_root == "")
            fn_root = fn_in.withoutExtension();
    }

    void defineParams()
    {
        addUsageLine("Separate different disconnected objects in a binary volume. ");
        addUsageLine("+In this way, small objects can be separated from large ones in ");
        addUsageLine("+reconstructions. This is very useful for reconstructions. Objects ");
        addUsageLine("+are written in separate files as binary volumes, too. Output volume ");
        addUsageLine("+number 1 is the background, number 2 corresponds to object 1, number 3 ");
        addUsageLine("+to object 2 ... This program also tells you the number of voxels on ");
        addUsageLine("+each object. ");
        addParamsLine("   -i <fn_in>              : Input image or volume");
        addParamsLine("  [--oroot <fn_root=\"\">] : Root filename for output");
        addParamsLine("                           : By default, the input name");
        addParamsLine("                           : The output masks are <fn_root>_000001.vol, ...");
        addParamsLine("  [--invert]               : Produce inverse masks");
        addParamsLine("  [--min_size <size=0>]    : Save if size is greater than this");
        addSeeAlsoLine("transform_morphology, transform_threshold, transform_mask, volume_segment");
    }

    void run()
    {
        FileName fn_out;
        FileName fn_ext = fn_in.getExtension();
        Image<double> I, label;
        I.read(fn_in);
        int object_no;
        if (ZSIZE(I())==1)
            object_no=labelImage2D(I(), label());
        else
            object_no=labelImage3D(I(), label());
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
            double number_elements = I().sum();
            if (number_elements > min_size)
            {
                fn_out.compose(fn_root, o+1, fn_ext);
                I.write(fn_out);
            }

            if (ZSIZE(I())==1)
                std::cout << "Image number " << o+1 << " contains " << number_elements
                << " pixels set to 1\n";
            else
                std::cout << "Volume number " << o+1 << " contains " << number_elements
                << " voxels set to 1\n";
        }
    }
};

/* ------------------------------------------------------------------------- */
/* Main                                                                      */
/* ------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgSeparateObjects prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
